// worker.js - the crosswordsmith engine Worker: one SWI-Prolog WASM instance,
// off the main thread, speaking the v1 envelope RPC (wasm-sdk-strategy §3-§5).
// The SDK facade (wasm/sdk/crosswordsmith.mjs) is its only intended client.
//
// Message protocol:
//   { type:"ping", id }  -> { type:"pong", id }     idle liveness probe (a
//        running solve does NOT service messages — gate #1 — so a mid-search
//        ping queues until the solve finishes)
//   { type:"warm" }      -> { type:"warmed" }       pre-load SWIPL + the qlf
//        now (the ~2MB wasm fetch + instantiate + qlf load) so a later
//        request on this worker skips that cost; init() is memoized
//   { type:"request", request:{v,id,verb,params,meta?}, stackLimit?, heartbeat? }
//       -> { type:"response", id, envelope }
//
// The envelope is browser_dispatch/3's discriminated JSON (status: success |
// failure | error), JSON.parse'd here. Every request gets EXACTLY ONE
// response; transport-level disasters (init failure, an escaped throw, a
// failed dispatch) are synthesized as {status:"error", error:{type:"internal"}}
// so the facade has a single reply path. Cancellation is worker.terminate()
// from the facade — a yield-free CPU-bound DFS cannot be interrupted
// cooperatively (gate #1: the heartbeat never drains this worker's message
// queue) — and a terminated worker posts nothing, so "cancelled" is a
// facade-side rejection type, never a wire envelope (strategy §4.2/§5.2).
//
// Requires, sitting next to this file: swipl-web.js / swipl-web.wasm /
// swipl-web.data (built with `ninja swipl-web`) and crosswordsmith.qlf (the
// qcompiled app - see wasm/README.md) — resolved through the content-hashed
// names in build-manifest.json when that is present (see loadBuildManifest
// below).

// SWI's WASM URL helpers (e.g. prolog.url_properties, used by consult(URL))
// reach for a browser `window`, which does not exist inside a Worker — without
// this shim they throw "window is not defined" and the consult 404s/fails.
// A Worker's global (self === globalThis) stands in fine: window.location then
// resolves to this worker's URL. Must run BEFORE importScripts.
self.window = self;

// Artifact provenance (supply-chain Batch 2): build-manifest.json, stamped by
// wasm/build/stamp-manifest.sh, maps each artifact to its content-hashed name
// (swipl-web.<sha256:12>.js …). Prefer those when the manifest is present so a
// CDN redeploy can never pair a new js with a long-cached stale wasm/data/qlf;
// fall back to the unhashed names when it is absent (hand-staged dev dir).
// Fetched with SYNCHRONOUS XHR on purpose: importScripts below is synchronous,
// so an async fetch cannot precede it at top level — and sync XHR is legal and
// jank-free off the main thread.
function loadBuildManifest() {
  try {
    const xhr = new XMLHttpRequest();
    xhr.open("GET", new URL("./build-manifest.json", self.location.href).href,
             false /* synchronous */);
    xhr.send();
    if (xhr.status >= 200 && xhr.status < 300) return JSON.parse(xhr.responseText);
  } catch (_) { /* no manifest — unhashed fallback below */ }
  return null;
}
const MANIFEST = loadBuildManifest();
function assetName(name) {
  const entry = MANIFEST && MANIFEST.artifacts && MANIFEST.artifacts[name];
  return (entry && entry.hashed) || name;
}

// Emscripten MODULARIZE output defines the global SWIPL after importScripts.
importScripts("./" + assetName("swipl-web.js"));

let Prolog = null;
let ready = null;

function init() {
  if (ready) return ready;
  ready = SWIPL({
    arguments: ["-q"],
    // Resolve swipl-web.wasm / swipl-web.data relative to this worker,
    // through the manifest's hashed names when stamped.
    locateFile: (file) => new URL("./" + assetName(file), self.location.href).href,
  }).then(async (module) => {
    Prolog = module.prolog;
    // NB: the search's stack cap is NOT set here. Prolog stacks are heap-backed and
    // grow with the wasm heap (ALLOW_MEMORY_GROWTH), so an over-large limit lets a
    // runaway grow the heap until the browser abort()s the tab uncatchably - we want
    // a recoverable resource_error(stack) first. But the search runs in a per-request
    // {engine:true} engine that does NOT inherit this parent engine's stack_limit
    // (it defaults to 1GB - measured), so capping it here would be a no-op for the
    // search. The cap is applied INSIDE the forEach goal instead (see onmessage).
    // Load the app. This .qlf MUST be produced by a same-pointer-size (wasm/node)
    // swipl, not native x86 — see solve_browser.pl / the deployment plan doc.
    // ONE explicit fetch into MEMFS, then a local consult (payload plan Phase
    // 4): Prolog.consult(URL) issues a probe request it then aborts before the
    // real fetch, and neither hits the HTTP cache — the known 2-requests-per-
    // worker / 4-per-cold-load qlf pattern. An ordinary fetch() is issued once
    // and obeys normal HTTP caching, so under immutable headers the second
    // worker's copy is a cache hit, not a server round trip. The URL is
    // resolved ABSOLUTE against this worker (no `window` here to derive a base
    // from), through the manifest's hashed name when stamped.
    const qlfUrl =
      new URL("./" + assetName("crosswordsmith.qlf"), self.location.href).href;
    const resp = await fetch(qlfUrl);
    if (!resp.ok)
      throw new Error("crosswordsmith.qlf fetch failed: HTTP " + resp.status);
    module.FS.writeFile("/crosswordsmith.qlf",
                        new Uint8Array(await resp.arrayBuffer()));
    await Prolog.consult("/crosswordsmith.qlf");
    return Prolog;
  });
  return ready;
}

self.onmessage = async (ev) => {
  const msg = ev.data || {};
  if (msg.type === "ping") { self.postMessage({ type: "pong", id: msg.id }); return; }
  if (msg.type === "warm") { await init(); self.postMessage({ type: "warmed" }); return; }
  if (msg.type !== "request") return;

  const req = msg.request || {};
  const id = req.id ?? null;
  const verb = typeof req.verb === "string" ? req.verb : "";
  // One reply path for the facade: anything that prevented the engine from
  // producing an envelope becomes an internal-error envelope.
  const fail = (message) => self.postMessage({
    type: "response", id,
    envelope: { v: 1, id, verb, status: "error",
                error: { type: "internal", message } },
  });

  try {
    await init();

    // Cap the logical stack UNDER the memory budget so a runaway search throws a
    // recoverable resource_error(stack) — mapped by the Prolog classifier to a
    // {type:"resource_exhausted"} envelope — instead of growing the wasm heap
    // until the browser abort()s the tab (verified by wasm/test/error_probe.mjs).
    // CRUCIAL: the cap MUST be set INSIDE the forEach goal below, not on this
    // parent engine - a {engine:true} throwaway engine does NOT inherit the
    // parent's stack_limit; it defaults to 1GB (measured), above a mobile
    // browser's ~300MB ceiling, so the tab could abort before the limit fires.
    // The limit rides in the goal as a BOUND variable (Limit), NOT string-
    // concat'd, so a caller-supplied value can neither inject Prolog nor arrive
    // as the wrong type. Validate to a POSITIVE INTEGER: a float/string/0/NaN
    // would type-error set_prolog_flag and abort the whole request, so fall
    // back to the 256MB default. Caller-overridable via msg.stackLimit.
    const stackLimit = Number.isInteger(msg.stackLimit) && msg.stackLimit > 0
      ? msg.stackLimit : 256000000;

    // forEach runs the goal asynchronously (yield every `heartbeat` inferences)
    // and engine:true gives each request a throwaway engine so per-run state
    // doesn't accumulate on the reused instance. (browser_dispatch/3 ALSO
    // resets the three instance globals at the top of every dispatch — the
    // throwaway engine does not clear module dynamics.) The heartbeat is not
    // load-bearing (gate #1): 50k is a mild default, caller-overridable.
    //
    // Verb + Payload ride as BOUND variables (never string-concatenated): the
    // JS->Prolog binding delivers the verb string as the atom dispatch/3
    // matches on, and the request JSON as one opaque payload atom that
    // browser_dispatch/3 json_read_dict's so every value lands with the right
    // Prolog type (JSON both ways — see solve_browser.pl).
    const heartbeat = msg.heartbeat || 50000;
    let json = null;
    await Prolog.forEach(
      "set_prolog_flag(stack_limit, Limit), browser_dispatch(Verb, Payload, Json)",
      { Verb: verb, Payload: JSON.stringify(req), Limit: stackLimit },
      (bindings) => { if (json === null) json = bindings.Json; },
      { engine: true, heartbeat }
    );
    if (json === null) {
      // browser_dispatch/3 is total (it backstops even a failing dispatch with
      // an internal-error envelope), so reaching this is an engine bug.
      fail("browser_dispatch/3 failed (engine bug: every dispatch path must tag an outcome or throw)");
      return;
    }
    self.postMessage({ type: "response", id, envelope: JSON.parse(json) });
  } catch (err) {
    fail(String((err && err.message) || err));
  }
};
