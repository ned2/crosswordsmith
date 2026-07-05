// worker.js - SPIKE: run the crosswordsmith solver off the main thread.
//
// The whole SWI-Prolog WASM instance lives inside this Web Worker, so the
// multi-second arrange DFS never touches the page's main thread. The page
// talks to it with two message types:
//
//   { type: "solve", input: {clues:[...], size, mode, bestEffort} }
//        -> posts { type: "result", output }   (the layout dict)
//        or  { type: "error",  message }
//   The page cancels a running search by calling worker.terminate() and
//   spawning a fresh worker (see main.js). There is no reliable way to
//   interrupt a pure-CPU Prolog DFS from outside a yield, so terminate() is
//   the pragmatic "Stop".
//
// Requires, sitting next to this file: swipl-web.js / swipl-web.wasm /
// swipl-web.data (built with `ninja swipl-web`) and crosswordsmith.qlf (the
// qcompiled app - see README). Falls back to consulting solve_browser.pl if
// the qlf is absent.

// SWI's WASM URL helpers (e.g. prolog.url_properties, used by consult(URL))
// reach for a browser `window`, which does not exist inside a Worker — without
// this shim they throw "window is not defined" and the consult 404s/fails.
// A Worker's global (self === globalThis) stands in fine: window.location then
// resolves to this worker's URL. Must run BEFORE importScripts.
self.window = self;

// Emscripten MODULARIZE output defines the global SWIPL after importScripts.
importScripts("./swipl-web.js");

let Prolog = null;
let ready = null;

function init() {
  if (ready) return ready;
  ready = SWIPL({
    arguments: ["-q"],
    // Resolve swipl-web.wasm / swipl-web.data relative to this worker.
    locateFile: (file) => new URL("./" + file, self.location.href).href,
  }).then(async (module) => {
    Prolog = module.prolog;
    // NB: the search's stack cap is NOT set here. Prolog stacks are heap-backed and
    // grow with the wasm heap (ALLOW_MEMORY_GROWTH), so an over-large limit lets a
    // runaway grow the heap until the browser abort()s the tab uncatchably - we want
    // a recoverable resource_error(stack) first. But the search runs in a per-solve
    // {engine:true} engine that does NOT inherit this parent engine's stack_limit
    // (it defaults to 1GB - measured), so capping it here would be a no-op for the
    // search. The cap is applied INSIDE the forEach goal instead (see onmessage).
    // Load the app. This .qlf MUST be produced by a same-pointer-size (wasm/node)
    // swipl, not native x86 — see solve_browser.pl / the plan doc.
    // Use an ABSOLUTE URL: inside a Worker there is no `window`, so SWI cannot
    // derive a base URL for a relative consult path (it logs "window is not
    // defined") and "./crosswordsmith.qlf" would 404. Resolve it ourselves.
    await Prolog.consult(new URL("./crosswordsmith.qlf", self.location.href).href);
    return Prolog;
  });
  return ready;
}

self.onmessage = async (ev) => {
  const msg = ev.data || {};
  // Idle liveness probe (health check for a crashed worker). NB: it is answered
  // only when the worker is IDLE - a running solve does NOT service new messages
  // (gate #1 finding, see below), so a ping sent mid-search is queued until the
  // solve finishes. To cancel a running solve, main.js uses worker.terminate().
  if (msg.type === "ping") { self.postMessage({ type: "pong", id: msg.id }); return; }
  if (msg.type !== "solve") return;
  try {
    await init();

    // Cap the logical stack UNDER the memory budget so a runaway search throws a
    // recoverable resource_error(stack) instead of growing the wasm heap until the
    // browser abort()s the tab (plan §6; verified by wasm/test/error_probe.mjs).
    // CRUCIAL: the cap MUST be set INSIDE the forEach goal below, not on this parent
    // engine - a {engine:true} throwaway engine does NOT inherit the parent's
    // stack_limit; it defaults to 1GB (measured). Setting it on the parent (or at
    // init) leaves the actual search running with a 1GB limit - above a mobile
    // browser's ~300MB ceiling, so the tab could abort before the limit ever fires.
    // Prepending it to the query applies it in the engine that runs the search, and
    // because the engine is thrown away each solve the value can't leak forward.
    // 256MB default; caller-overridable via msg.stackLimit.
    const stackLimit = msg.stackLimit || 256000000;

    // forEach runs the search asynchronously with the yield mechanism (every
    // `heartbeat` inferences), and engine:true gives each solve a throwaway engine
    // so per-run state doesn't accumulate on the reused instance.
    //
    // gate #1 finding (wasm/test/yield_probe.mjs, plan §7): the heartbeat yield
    // cooperates among Prolog engines but does NOT drain this worker's JS message
    // queue - so a running solve can't service pings or a "cancel" message at ANY
    // heartbeat (measured: 0 mid-search pongs at 10k..2e9). The UI stays responsive
    // regardless, because the search runs in this Worker, OFF the page's main
    // thread (measured main-thread drift 0ms). So the heartbeat is not load-bearing
    // here - it costs a few % throughput for no UX gain; 50k is a mild, harmless
    // default. Cancellation is worker.terminate() (prompt hard kill; see main.js),
    // not in-worker messaging. Caller may override heartbeat via msg.heartbeat.
    //
    // JSON in (not a query binding): the JS->Prolog binding delivers JS strings as
    // Prolog atoms, breaking the solver's string "answer" check; solve_browser_str
    // parses the JSON in Prolog so types land right. Symmetric with the JSON-out
    // path (the layout schema's null/bool survive json_write_dict + JSON.parse).
    const heartbeat = msg.heartbeat || 50000;
    let json = null;
    await Prolog.forEach(
      "set_prolog_flag(stack_limit, " + stackLimit + "), solve_browser_str(Payload, Json)",
      { Payload: JSON.stringify(msg.input) },
      (bindings) => { if (json === null) json = bindings.Json; },
      { engine: true, heartbeat }
    );
    if (json === null) {
      self.postMessage({ type: "error", message: "no layout (search failed)" });
    } else {
      self.postMessage({ type: "result", output: JSON.parse(json) });
    }
  } catch (err) {
    self.postMessage({ type: "error", message: String(err && err.message || err) });
  }
};
