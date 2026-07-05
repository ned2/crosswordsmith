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
    // Give Prolog room; these are heap-backed and grow with the wasm heap
    // (ALLOW_MEMORY_GROWTH is on by default in the build).
    Prolog.call("set_prolog_flag(stack_limit, 1 000 000 000)");
    // Load the app. Production: one precompiled .qlf fetched from a URL.
    await Prolog.consult("./crosswordsmith.qlf");
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
    let json = null;
    await Prolog.forEach(
      "solve_browser_json(In, Json)",
      { In: msg.input },
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
