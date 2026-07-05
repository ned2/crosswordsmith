// main.js - SPIKE: page controller. Keeps the SWI-Prolog instance in a Worker
// (worker.js) so the UI never blocks, and implements "Stop" as terminate() +
// respawn. A warm spare worker is kept ready so the next Solve starts instantly
// instead of paying the wasm fetch/instantiate + qlf-load cost again.

const $ = (id) => document.getElementById(id);

let active = null;   // the worker currently running a solve (or idle+ready)
let spare = null;    // a pre-warmed standby worker

function spawn() {
  return new Worker("./worker.js");
}

// Ensure `active` exists; promote the spare if we have one, then re-arm a spare.
function ensureActive() {
  if (!active) active = spare || spawn();
  spare = spawn();               // warm the next one in the background
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
    } else {
      $("status").textContent = "error: " + msg.message;
    }
    $("solve").disabled = false;
    $("stop").disabled = true;
  };
  worker.postMessage({ type: "solve", input });
};

$("stop").onclick = () => {
  if (active) active.terminate();     // hard-cancel a CPU-bound DFS
  active = null;
  ensureActive();                      // promote/re-arm so the next Solve is warm
  setBusy(false);
  $("status").textContent = "stopped";
};

// Warm the first instance up front so the first Solve is responsive.
ensureActive();
