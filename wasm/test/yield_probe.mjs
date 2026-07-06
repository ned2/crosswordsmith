// yield_probe.mjs - gate #1: how the browser build behaves under a genuinely
// large search (~38M inferences), across four questions the 5-word toy never
// exercised:
//
//   T1 idle ping      - does the worker's ping/pong path work at all (sanity)?
//   T2 heavy search   - does a big search COMPLETE in-browser, and does the PAGE
//                       (main thread) stay responsive while it runs? measured as
//                       setInterval drift on the main thread + mid-search pongs
//                       from the worker, swept over the heartbeat.
//   T3 terminate      - is worker.terminate() a prompt cancel mid-search?
//
// The heavy word set is fixtures/ladder_15x15_36w.pl (36 words, 15x15,
// satisfiable; ~38.3M inferences, parity-certified). Prereq: server for the
// wasm/ tree on $URL's origin + Playwright in wasm/test.
// Run: node wasm/test/yield_probe.mjs  (exit 0 = all three gates pass; each of
// T1/T2/T3 asserts and a failure makes the process exit 1).

import { chromium } from 'playwright';
import { HEAVY as WORDS } from './heavy_words.mjs';

const URL = process.env.URL || 'http://127.0.0.1:8080/client/harness.html?noauto=1';
const SEARCH_INF = 38275505;
// Responsiveness gate: gate #1 measured 0ms main-thread drift (the Worker isolates
// the UI). Anything near the search's multi-second wall would mean the main thread
// was blocked; 200ms is a generous margin that still catches that regression.
const DRIFT_MAX_MS = 200;
// terminate() is a synchronous thread kill; gate #1 measured it returns in ~0ms.
const TERMINATE_MAX_MS = 500;
const WARM = ['CAT','CAR','ARC','RAT','TAR'];

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();
const pageLogs = [];
page.on('pageerror', e => pageLogs.push(`[pageerror] ${e.message}`));
await page.goto(URL, { waitUntil: 'load' });

const fails = [];   // each entry: a human-readable reason this run failed
const check = (cond, reason) => { if (!cond) fails.push(reason); };

// Shared page-side helper source: make a worker and warm it (load wasm+qlf, so
// later timing excludes load). Returned as a string spliced into each evaluate.
const HELPERS = `
  let reqSeq = 0;
  function payload(words, opts){ return { clues: words.map(a=>({answer:a})), size: opts.size, mode:'fixed', bestEffort: !!opts.best }; }
  function requestMsg(params, heartbeat){
    return { type:'request', request: { v:1, id:'y-'+(++reqSeq), verb:'arrange', params }, heartbeat };
  }
  async function makeWarmWorker(heartbeat){
    const w = new Worker('./worker.js');
    await new Promise((res, rej) => {
      w.onmessage = (ev)=>{ const m=ev.data; if(m.type!=='response') return;
        const e=m.envelope; if(e.status==='success') res(); else rej(new Error(JSON.stringify(e))); };
      w.postMessage(requestMsg(payload(${JSON.stringify(WARM)}, {size:5, best:true}), heartbeat));
    });
    return w;
  }
`;

// T1 - idle ping sanity: 5 pings 40ms apart to an IDLE (post-warm) worker.
const t1 = await page.evaluate(async ({ HELPERS }) => {
  eval(HELPERS);
  const w = await makeWarmWorker(50000);
  let pongs = 0;
  w.onmessage = (ev) => { if (ev.data.type === 'pong') pongs++; };
  for (let i = 0; i < 5; i++) { w.postMessage({ type: 'ping', id: i }); await new Promise(r => setTimeout(r, 40)); }
  await new Promise(r => setTimeout(r, 100));
  w.terminate();
  return { sent: 5, pongs };
}, { HELPERS });
const t1ok = t1.pongs === t1.sent;
check(t1ok, `T1: idle ping/pong broken (sent ${t1.sent}, got ${t1.pongs})`);
console.log(`T1 idle ping/pong: sent ${t1.sent}, got ${t1.pongs}  -> ${t1ok ? 'PASS' : 'FAIL (path broken)'}`);

// T2 - heavy search: completion + main-thread responsiveness + mid-search worker
// servicing, swept over heartbeat (last value >> total inferences = no-yield control).
console.log('\nT2 heavy search (36 words / 15x15, ~38.3M inf):');
console.log('  heartbeat      solveMs  placed  mainDriftMax  midSearchPongs   note');
for (const hb of [10000, 50000, 500000, 2000000000]) {
  const r = await page.evaluate(async ({ HELPERS, WORDS, hb }) => {
    eval(HELPERS);
    const w = await makeWarmWorker(hb);
    let placed = null, err = null, solveMs = 0, midPongs = 0;
    const start = performance.now();
    // main-thread heartbeat: fire every 50ms, record drift + ping the worker
    const ticks = [];
    let id = 0;
    const timer = setInterval(() => {
      ticks.push(performance.now());
      const k = ++id; w.__sent = w.__sent || {}; w.__sent[k] = performance.now();
      w.postMessage({ type: 'ping', id: k });
    }, 50);
    await new Promise((res) => {
      w.onmessage = (ev) => {
        const m = ev.data;
        if (m.type === 'pong') { midPongs++; }
        else if (m.type === 'response') {
          solveMs = performance.now() - start;
          const e = m.envelope;
          if (e.status === 'success') placed = e.result.words.length;
          else err = e.status === 'error' ? e.error.message : JSON.stringify(e.detail);
          res();
        }
      };
      w.postMessage(requestMsg(payload(WORDS, { size: 15, best: false }), hb));
    });
    clearInterval(timer);
    w.terminate();
    // main-thread drift: worst deviation of consecutive ticks from the 50ms target
    let driftMax = 0;
    for (let i = 1; i < ticks.length; i++) driftMax = Math.max(driftMax, Math.abs((ticks[i] - ticks[i - 1]) - 50));
    return { solveMs: Math.round(solveMs), placed, err, driftMax: Math.round(driftMax), midPongs, ticks: ticks.length };
  }, { HELPERS, WORDS, hb });
  check(r.err === null && r.placed === 36,
    `T2: heavy search did not complete at heartbeat ${hb} (placed=${r.placed}, err=${r.err})`);
  check(r.driftMax < DRIFT_MAX_MS,
    `T2: main-thread drift ${r.driftMax}ms >= ${DRIFT_MAX_MS}ms at heartbeat ${hb} (UI blocked)`);
  const note = r.err ? `ERROR ${r.err}` : (hb >= SEARCH_INF ? 'no-yield control' : '');
  console.log(`  ${String(hb).padStart(11)}  ${String(r.solveMs).padStart(7)}  ${String(r.placed).padStart(6)}  ${String(r.driftMax + 'ms').padStart(12)}  ${String(r.midPongs).padStart(14)}   ${note}`);
}

// T3 - terminate() as cancel: start the heavy search, terminate at 800ms, then
// watch for 3s. `resultPhase` records WHEN (relative to terminate) a result first
// arrived: 'after' => terminate did not cancel; 'before' => the search finished
// before we could cancel (inconclusive - too fast to test); null => cancelled
// promptly (the pass case). This distinction matters because a faster machine or a
// perf win could finish the search under 800ms and make a naive "any result" flag
// misreport a working cancel as a failure.
const t3 = await page.evaluate(async ({ HELPERS, WORDS }) => {
  eval(HELPERS);
  const w = await makeWarmWorker(50000);
  let resultPhase = null, terminated = false, killMs = 0;
  const start = performance.now();
  w.onmessage = (ev) => {
    if (ev.data.type === 'response' && resultPhase === null)
      resultPhase = terminated ? 'after' : 'before';
  };
  w.postMessage(requestMsg(payload(WORDS, { size: 15, best: false }), 50000));
  await new Promise(r => setTimeout(r, 800));
  const t0 = performance.now(); w.terminate(); terminated = true; killMs = Math.round(performance.now() - t0);
  await new Promise(r => setTimeout(r, 3000));   // the search would need ~4s more to finish
  return { terminateReturnedMs: killMs, resultPhase, watchedMs: Math.round(performance.now() - start) };
}, { HELPERS, WORDS });
const t3cancelled = t3.resultPhase === null;
check(t3cancelled, t3.resultPhase === 'after'
  ? 'T3: terminate() did NOT cancel - a result arrived after terminate'
  : `T3: search finished before the 800ms terminate (inconclusive - resultPhase=${t3.resultPhase}); lower heartbeat/raise word count`);
check(t3.terminateReturnedMs < TERMINATE_MAX_MS,
  `T3: terminate() took ${t3.terminateReturnedMs}ms (>= ${TERMINATE_MAX_MS}ms)`);
console.log(`\nT3 terminate() cancel @800ms: terminate() returned in ${t3.terminateReturnedMs}ms, ` +
  `resultPhase=${t3.resultPhase} (watched ${t3.watchedMs}ms)  -> ` +
  `${t3cancelled ? 'PASS (cancelled promptly)' : 'FAIL'}`);

if (pageLogs.length) { console.log('\nPAGE:'); console.log(pageLogs.join('\n')); }
await browser.close();

if (fails.length) {
  console.log(`\ngate #1 FAILED (${fails.length}):`);
  for (const f of fails) console.log(`  - ${f}`);
  process.exit(1);
}
console.log('\ngate #1 OK (T1 ping, T2 completion+responsiveness, T3 prompt cancel)');
process.exit(0);
