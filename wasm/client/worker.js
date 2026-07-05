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
    // Prolog stacks are heap-backed and grow with the wasm heap
    // (ALLOW_MEMORY_GROWTH). Keep the limit UNDER the memory budget (~256MB, not
    // 1GB) so a runaway search throws a recoverable resource_error(stack) rather
    // than growing the heap until the browser abort()s the tab uncatchably.
    Prolog.call("set_prolog_flag(stack_limit, 256 000 000)");
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
  if (msg.type !== "solve") return;
  try {
    await init();
    // forEach runs asynchronously with the yield mechanism, so the worker's
    // own message loop stays live (it can see a future terminate promptly) and
    // the search is cancellable at heartbeat boundaries. engine:true gives each
    // solve a throwaway engine so per-run state does not accumulate on the one
    // reused instance. Tune heartbeat down if a tight DFS starves the loop.
    // Input goes in as a query binding (JS object -> Prolog dict, clean and
    // lossless). Output comes back as a JSON atom we JSON.parse here - faithful
    // for this app because the layout schema contains JSON null (empty cells),
    // which the dict-binding path would turn into the string "null". See
    // solve_browser.pl for the full rationale.
    // Pass the input as a JSON STRING, not a JS object: the JS->Prolog binding
    // turns JS strings into Prolog atoms, which breaks the solver's string
    // "answer" check. solve_browser_str parses the JSON in Prolog so types land
    // correctly (strings/numbers/null/bool). Symmetric with the JSON-out path.
    let json = null;
    await Prolog.forEach(
      "solve_browser_str(Payload, Json)",
      { Payload: JSON.stringify(msg.input) },
      (bindings) => { if (json === null) json = bindings.Json; },
      { engine: true, heartbeat: 50000 }
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
