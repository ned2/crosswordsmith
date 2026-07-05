// main.js - SPIKE: page controller. Keeps the SWI-Prolog instance in a Worker
// (worker.js) so the UI never blocks, and implements "Stop" as terminate() +
// respawn. A warm standby worker is kept ready - it has already fetched +
// instantiated the wasm and consulted the qlf (via a "warm" message) - so the
// Solve after a Stop starts without paying that ~2MB load cost again.

const $ = (id) => document.getElementById(id);

let active = null;   // the worker currently running a solve (or idle+ready)
let spare = null;    // a warmed standby worker (SWIPL+qlf already loaded)

function spawn() {
  return new Worker("./worker.js");
}

// Spawn a worker AND tell it to pre-load SWIPL+qlf now, so a later solve is fast.
function spawnWarm() {
  const w = spawn();
  w.postMessage({ type: "warm" });
  return w;
}

// Ensure `active` exists (promoting the already-warmed spare if we have one) and
// keep exactly ONE warmed spare on standby. Spawns at most one worker per call and
// never orphans one, so repeated Solve clicks don't leak workers.
function ensureActive() {
  if (!active) { active = spare || spawnWarm(); spare = null; }
  if (!spare) spare = spawnWarm();
  return active;
}

function setBusy(busy) {
  $("solve").disabled = busy;
  $("stop").disabled = !busy;
  $("status").textContent = busy ? "solving…" : "idle";
}

function parseInput() {
  const clues = $("words").value
    .split("\n").map((s) => s.trim()).filter(Boolean)
    .map((answer) => ({ answer }));
  return {
    clues,
    size: parseInt($("size").value, 10),
    mode: "fixed",
    bestEffort: $("best").checked,
  };
}

function renderGrid(layout) {
  const rows = layout.grid.map((row) =>
    row.map((cell) =>
      // Empty cells are JSON null; a filled cell is an object with .letter.
      cell === null
        ? "·"                       // middle dot for a block/empty
        : cell.letter
    ).join(" ")
  );
  $("grid").textContent = rows.join("\n");
}

$("solve").onclick = () => {
  const worker = ensureActive();
  const input = parseInput();
  setBusy(true);
  $("grid").textContent = "";
  worker.onmessage = (ev) => {
    const msg = ev.data;
    if (msg.type === "result") {
      renderGrid(msg.output);
      $("status").textContent =
        `placed ${msg.output.words.length} words on ${msg.output.gridLength}×${msg.output.gridLength}`;
    } else if (msg.type === "error") {
      $("status").textContent = "error: " + msg.message;
    } else {
      return;   // ignore non-terminal replies (e.g. a warmed/pong from this worker)
    }
    $("solve").disabled = false;
    $("stop").disabled = true;
  };
  worker.postMessage({ type: "solve", input });
};

$("stop").onclick = () => {
  if (active) active.terminate();     // hard-cancel a CPU-bound DFS
  active = null;
  ensureActive();                      // promote the warmed spare + re-arm one
  setBusy(false);
  $("status").textContent = "stopped";
};

// Warm the first instance (and its spare) up front so the first Solve is fast.
ensureActive();
