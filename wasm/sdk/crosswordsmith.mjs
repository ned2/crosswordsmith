// crosswordsmith.mjs - the client-side SDK facade (wasm-sdk-strategy §5).
//
// createCrosswordsmith() bakes the whole Worker story in: it spawns the engine
// Worker (wasm/client/worker.js) plus one warmed spare, serializes requests
// through an id-correlated single-flight FIFO queue, and dresses cancellation
// as an AbortSignal over worker.terminate() with generation fencing. Consumers
// only ever see typed envelopes:
//
//   import { createCrosswordsmith, randomSeed } from ".../sdk/crosswordsmith.mjs";
//   const cw = await createCrosswordsmith();          // resolves engine-ready
//   const res = await cw.arrange({ clues, size: 15, seed: randomSeed() },
//                                { signal });
//   if (res.status === "success") render(res.result);
//
// Design decisions this file implements (see the plan for the rationale):
// - typed method-per-verb over ONE private RPC; no generic query(goal).
// - the reply is the engine's discriminated envelope (success | failure |
//   error) — all three RESOLVE the promise; rejection is reserved for
//   cancellation (CancelledError) and transport disasters.
// - cancellation = worker.terminate() + promote the warm spare. A yield-free
//   CPU-bound DFS cannot be interrupted cooperatively, and a terminated worker
//   posts nothing — so "cancelled" is a promise-rejection type only, never a
//   wire envelope. Single-flight means at most ONE request is ever posted to
//   the instance, so terminate kills exactly that request; still-queued
//   requests simply continue on the promoted spare. Late replies from a
//   terminated worker are dropped by generation fencing.
// - capabilities() is ENGINE-sourced (a round-trip), not a JS hardcode that
//   would lie after a qlf swap (OQ-4). createCrosswordsmith() performs one
//   capabilities round-trip at init, so it resolves only when the engine is
//   actually ready and cw.version can report engine provenance.

const SDK_VERSION = "0.1.0";

// Default engine-asset location: the in-repo layout, where sdk/ and client/
// are siblings. Routed through a const — NOT written literally inside
// new URL(...) — because bundlers (Vite, webpack 5) statically pattern-match
// `new URL("<literal>", import.meta.url)` and try to resolve the path at
// build time; wasm/client/ is copied out of band and is never in a consumer's
// module graph, so that analysis can only produce a build warning. The const
// keeps resolution a runtime concern, which is what it is.
// NB: if this file is BUNDLED, import.meta.url becomes the bundled chunk's
// URL and this default resolves relative to that — bundling consumers must
// pass assetBaseUrl explicitly.
const DEFAULT_ASSET_DIR = "../client/";

// Draw an app-side seed for "regenerate": an integer in [0, 1e9], matching the
// CLI's --shuffle draw (random_between(0, 1_000_000_000, N)). Keep seeds in
// this range: a float or a value above 2^53 fails the engine's integer guard
// (typed validation error) or risks JSON precision loss.
export function randomSeed() {
  return Math.floor(Math.random() * 1_000_000_001);
}

// The JS-side cancellation rejection. Deliberately NOT an envelope type.
export class CancelledError extends Error {
  constructor(message = "cancelled") {
    super(message);
    this.name = "CancelledError";
  }
}

export async function createCrosswordsmith(options = {}) {
  // Engine assets (worker.js + swipl-web.{js,wasm,data} + crosswordsmith.qlf)
  // live in one directory; default to the in-repo layout. Consumers hosting
  // the assets elsewhere (or bundling this file — see DEFAULT_ASSET_DIR) pass
  // assetBaseUrl.
  const base = options.assetBaseUrl ?? new URL(DEFAULT_ASSET_DIR, import.meta.url);
  const assetBaseUrl = new URL(String(base).endsWith("/") ? base : base + "/",
                               import.meta.url);
  const workerUrl = new URL("worker.js", assetBaseUrl);
  const defaults = {};
  if (Number.isInteger(options.stackLimit) && options.stackLimit > 0) {
    defaults.stackLimit = options.stackLimit;
  }
  if (Number.isInteger(options.heartbeat) && options.heartbeat > 0) {
    defaults.heartbeat = options.heartbeat;
  }

  let seq = 0;          // request id counter (r-1, r-2, ...)
  const queue = [];     // jobs accepted but not yet posted (single-flight FIFO)
  let inflight = null;  // the ONE job currently posted to the active worker
  let active = null;    // { worker, dead } - the serving instance
  let spare = null;     // a warmed standby (already fetched wasm + loaded qlf)
  let disposed = false;

  // --- worker lifecycle ------------------------------------------------------

  function spawnWarm() {
    const record = { worker: new Worker(workerUrl), dead: false };
    record.worker.onmessage = (ev) => onWorkerMessage(record, ev.data || {});
    record.worker.onerror = (ev) => onWorkerError(record, ev);
    record.worker.postMessage({ type: "warm" });
    return record;
  }

  function ensureSpare() {
    if (!disposed && !spare) spare = spawnWarm();
  }

  // Retire a worker for good: no message from it may ever settle anything
  // again (generation fence — a late "response" from a terminated instance
  // must not resolve a cancelled promise).
  function retire(record) {
    if (!record) return;
    record.dead = true;
    record.worker.terminate();
  }

  function promoteSpare() {
    active = spare || spawnWarm();
    spare = null;
    ensureSpare();
  }

  function onWorkerMessage(record, msg) {
    if (record.dead) return;                       // fenced: terminated instance
    if (msg.type !== "response") return;           // warmed/pong: nothing to do
    if (record !== active) return;                 // stale instance
    if (!inflight || msg.id !== inflight.id) return; // stray/duplicate reply
    const job = inflight;
    inflight = null;
    settle(job, "resolve", msg.envelope);
    dispatchNext();
  }

  // A script-level worker error (not a solver outcome — those arrive as error
  // envelopes): fail the posted job and swap in the spare so the SDK recovers.
  function onWorkerError(record, ev) {
    if (record.dead || record !== active) return;
    retire(active);
    promoteSpare();
    if (inflight) {
      const job = inflight;
      inflight = null;
      settle(job, "reject",
             new Error("engine worker error: " + (ev && ev.message || "unknown")));
      dispatchNext();
    }
  }

  // --- the single-flight queue ----------------------------------------------

  function settle(job, how, value) {
    if (job.settled) return;
    job.settled = true;
    if (job.signal && job.onAbort) {
      job.signal.removeEventListener("abort", job.onAbort);
    }
    (how === "resolve" ? job.resolve : job.reject)(value);
  }

  function dispatchNext() {
    if (disposed || inflight) return;
    while (queue.length > 0 && queue[0].settled) queue.shift();
    if (queue.length === 0) return;
    inflight = queue.shift();
    active.worker.postMessage({
      type: "request",
      request: { v: 1, id: inflight.id, verb: inflight.verb, params: inflight.params },
      ...defaults,
    });
  }

  function cancelJob(job) {
    if (job.settled) return;
    if (job === inflight) {
      // The only request posted to the instance: hard-kill it and promote the
      // warmed spare. Queued jobs are untouched — they were never posted.
      inflight = null;
      retire(active);
      promoteSpare();
      settle(job, "reject", new CancelledError());
      dispatchNext();
    } else {
      const at = queue.indexOf(job);
      if (at >= 0) queue.splice(at, 1);
      settle(job, "reject", new CancelledError());
    }
  }

  function request(verb, params, { signal } = {}) {
    if (disposed) {
      return Promise.reject(new Error("createCrosswordsmith instance disposed"));
    }
    if (signal && signal.aborted) {
      return Promise.reject(new CancelledError());
    }
    return new Promise((resolve, reject) => {
      const job = {
        id: "r-" + (++seq), verb, params,
        resolve, reject, signal, onAbort: null, settled: false,
      };
      if (signal) {
        job.onAbort = () => cancelJob(job);
        signal.addEventListener("abort", job.onAbort, { once: true });
      }
      queue.push(job);
      dispatchNext();
    });
  }

  // Unwrap a success envelope or throw (for SDK-internal round-trips).
  async function requireSuccess(verb, params) {
    const env = await request(verb, params);
    if (env.status !== "success") {
      const why = env.status === "error"
        ? env.error.type + ": " + env.error.message
        : "failure";
      throw new Error("crosswordsmith " + verb + " failed (" + why + ")");
    }
    return env.result;
  }

  // --- boot: active + spare, then one engine round-trip ----------------------

  active = spawnWarm();
  ensureSpare();
  const initCaps = await requireSuccess("capabilities", {});

  return {
    // The typed verb surface. The reply is the engine's envelope: check
    // res.status ("success" | "failure" | "error") — a legitimate "no layout"
    // is a failure, not a rejection. { signal } cancels: queued requests are
    // dequeued, an in-flight one hard-kills its worker (terminate + warm-spare
    // promote); either way the promise rejects with CancelledError.
    arrange(input, opts = {}) {
      return request("arrange", input, opts);
    },

    // Validate a layout (an arrange result) against a house-style profile.
    // Always a success envelope when the input is well-formed: a FAIL verdict
    // lives INSIDE the report — it is an answer, not an error.
    lint(input, opts = {}) {
      return request("lint", input, opts);
    },

    // Transform a layout to ipuz (result = the ipuz JSON document) or Exolve
    // (result = {format:"text", body}).
    export(input, opts = {}) {
      return request("export", input, opts);
    },

    // Engine-sourced, fresh each call (OQ-4): the verb list the loaded qlf
    // actually answers for, plus engine provenance.
    capabilities() {
      return requireSuccess("capabilities", {});
    },

    // sdk = this facade; engine = what the init capabilities round-trip
    // reported (the qlf's own answer, not a JS constant).
    version: Object.freeze({
      sdk: SDK_VERSION,
      engine: initCaps.engine,
    }),

    // Terminate both workers and reject anything pending. The instance is
    // dead afterwards; create a new one to continue.
    dispose() {
      if (disposed) return;
      disposed = true;
      const pending = inflight ? [inflight, ...queue] : [...queue];
      inflight = null;
      queue.length = 0;
      for (const job of pending) {
        settle(job, "reject", new CancelledError("disposed"));
      }
      retire(active);
      retire(spare);
      active = spare = null;
    },
  };
}
