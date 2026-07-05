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
// satisfiable; ~38.3M inferences, parity-certified). Prereq: server for
// wasm/client on $URL (default http://127.0.0.1:8080/) + Playwright in wasm/test.
// Run: node wasm/test/yield_probe.mjs

import { chromium } from 'playwright';

const URL = process.env.URL || 'http://127.0.0.1:8080/';
const SEARCH_INF = 38275505;
const WORDS = ['DFAD','FCC','FCB','BEED','CBED','AFCC','DED','DEF','BFF','BCD','DBED','CDF',
  'FFCC','CCD','EAB','FCF','EAA','EFAF','ABFA','BBEF','FFE','EFEE','ABCB','EFD','DACA','FAFB',
  'ACA','CBF','DEAA','AFBF','AEFD','EADF','EDDE','CEF','CADF','FDD'];
const WARM = ['CAT','CAR','ARC','RAT','TAR'];

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();
const pageLogs = [];
page.on('pageerror', e => pageLogs.push(`[pageerror] ${e.message}`));
await page.goto(URL, { waitUntil: 'load' });

// Shared page-side helper source: make a worker and warm it (load wasm+qlf, so
// later timing excludes load). Returned as a string spliced into each evaluate.
const HELPERS = `
  function payload(words, opts){ return { clues: words.map(a=>({answer:a})), size: opts.size, mode:'fixed', bestEffort: opts.best }; }
  async function makeWarmWorker(heartbeat){
    const w = new Worker('./worker.js');
    await new Promise((res, rej) => {
      w.onmessage = (ev)=>{ const m=ev.data; if(m.type==='result') res(); else if(m.type==='error') rej(new Error(m.message)); };
      w.postMessage({ type:'solve', input: payload(${JSON.stringify(WARM)}, {size:5, best:true}), heartbeat });
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
console.log(`T1 idle ping/pong: sent ${t1.sent}, got ${t1.pongs}  -> ${t1.pongs === t1.sent ? 'path works' : 'PATH BROKEN'}`);

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
        else if (m.type === 'result') { solveMs = performance.now() - start; placed = (m.output.words || []).length; res(); }
        else if (m.type === 'error') { solveMs = performance.now() - start; err = m.message; res(); }
      };
      w.postMessage({ type: 'solve', input: payload(WORDS, { size: 15, best: false }), heartbeat: hb });
    });
    clearInterval(timer);
    w.terminate();
    // main-thread drift: worst deviation of consecutive ticks from the 50ms target
    let driftMax = 0;
    for (let i = 1; i < ticks.length; i++) driftMax = Math.max(driftMax, Math.abs((ticks[i] - ticks[i - 1]) - 50));
    return { solveMs: Math.round(solveMs), placed, err, driftMax: Math.round(driftMax), midPongs, ticks: ticks.length };
  }, { HELPERS, WORDS, hb });
  const note = r.err ? `ERROR ${r.err}` : (hb >= SEARCH_INF ? 'no-yield control' : '');
  console.log(`  ${String(hb).padStart(11)}  ${String(r.solveMs).padStart(7)}  ${String(r.placed).padStart(6)}  ${String(r.driftMax + 'ms').padStart(12)}  ${String(r.midPongs).padStart(14)}   ${note}`);
}

// T3 - terminate() as cancel: start the heavy search, terminate at 800ms, then
// watch for 3s. A prompt cancel means NO result arrives after terminate.
const t3 = await page.evaluate(async ({ HELPERS, WORDS }) => {
  eval(HELPERS);
  const w = await makeWarmWorker(50000);
  let resultAfterTerminate = false, killMs = 0;
  const start = performance.now();
  w.onmessage = (ev) => { if (ev.data.type === 'result' || ev.data.type === 'error') resultAfterTerminate = true; };
  w.postMessage({ type: 'solve', input: payload(WORDS, { size: 15, best: false }), heartbeat: 50000 });
  await new Promise(r => setTimeout(r, 800));
  const t0 = performance.now(); w.terminate(); killMs = Math.round(performance.now() - t0);
  await new Promise(r => setTimeout(r, 3000));   // the search would need ~4s more to finish
  return { terminateReturnedMs: killMs, resultAfterTerminate, watchedMs: Math.round(performance.now() - start) };
}, { HELPERS, WORDS });
console.log(`\nT3 terminate() cancel @800ms: terminate() returned in ${t3.terminateReturnedMs}ms, ` +
  `result after terminate = ${t3.resultAfterTerminate} (watched ${t3.watchedMs}ms)  -> ` +
  `${t3.resultAfterTerminate ? 'NOT cancelled' : 'cancelled promptly'}`);

if (pageLogs.length) { console.log('\nPAGE:'); console.log(pageLogs.join('\n')); }
await browser.close();
