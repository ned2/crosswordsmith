// error_probe.mjs - gate the browser build's UNHAPPY paths (happy path is
// headless.mjs). Each case runs in a fresh worker, then re-solves a KNOWN-GOOD
// request on the SAME worker (liveness) - proving the failure was recoverable, not
// a dead tab / wedged instance:
//
//   throw : missing `size` -> solve_browser_json throws browser_missing_size ->
//           worker posts {type:"error"} (the forEach-rejects path).
//   fail  : infeasible STRICT input -> arrange_solve fails (no exception) -> worker
//           posts the generic "no layout (search failed)".
//   heavy-default : the ~38M search under the normal 256MB cap COMPLETES (36 words)
//           - the cap doesn't harm real work.
//   heavy-capped  : the same search under a tiny 300KB stack cap -> a RECOVERABLE
//           resource_error(stack), NOT a wasm-heap-growth tab abort (plan §6). This
//           only holds because the cap is applied INSIDE the forEach engine; a
//           parent-set limit is ignored by the {engine:true} search (defaults 1GB).
//
// Prereq: server for wasm/client on $URL (default http://127.0.0.1:8080/) +
// Playwright in wasm/test. Run: node wasm/test/error_probe.mjs (exit 0 = all pass).

import { chromium } from 'playwright';

const URL = process.env.URL || 'http://127.0.0.1:8080/';
const TOY = ['CAT', 'CAR', 'ARC', 'RAT', 'TAR'];
const HEAVY = ['DFAD','FCC','FCB','BEED','CBED','AFCC','DED','DEF','BFF','BCD','DBED','CDF',
  'FFCC','CCD','EAB','FCF','EAA','EFAF','ABFA','BBEF','FFE','EFEE','ABCB','EFD','DACA','FAFB',
  'ACA','CBF','DEAA','AFBF','AEFD','EADF','EDDE','CEF','CADF','FDD'];
const clues = (ws) => ws.map(a => ({ answer: a }));

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();
page.on('pageerror', e => console.log(`[pageerror] ${e.message}`));
await page.goto(URL, { waitUntil: 'load' });

// Run one solve message, capture its reply, then send a KNOWN-GOOD solve on the
// SAME worker and report whether it still succeeds (liveness).
async function run(msg) {
  return await page.evaluate(async ({ msg, TOY }) => {
    const good = { clues: TOY.map(a => ({ answer: a })), size: 5, mode: 'fixed', bestEffort: true };
    const w = new Worker('./worker.js');
    const once = (send) => new Promise((resolve) => {
      w.onmessage = (ev) => {
        const m = ev.data;
        if (m.type === 'result') resolve({ type: 'result', words: (m.output.words || []).length });
        else if (m.type === 'error') resolve({ type: 'error', message: String(m.message).replace(/\s+/g, ' ').slice(0, 160) });
      };
      setTimeout(() => resolve({ type: 'timeout' }), 60000);
      w.postMessage(send);
    });
    const reply = await once(msg);
    const alive = await once({ type: 'solve', input: good });   // same worker, known-good
    w.terminate();
    return { reply, alive };
  }, { msg, TOY });
}

const recoverable = (a) => a.type === 'result' && a.words === 5;

const cases = [
  { name: 'throw:missing-size',
    msg: { type: 'solve', input: { clues: clues(TOY), mode: 'fixed' } },
    ok: ({ reply, alive }) => reply.type === 'error' && recoverable(alive) },

  { name: 'fail:infeasible-strict',
    msg: { type: 'solve', input: { clues: clues(['AAAA', 'BBBB', 'CCCC']), size: 4, mode: 'fixed', bestEffort: false } },
    ok: ({ reply, alive }) => reply.type === 'error' && /no layout/i.test(reply.message) && recoverable(alive) },

  { name: 'heavy-default-completes',
    msg: { type: 'solve', input: { clues: clues(HEAVY), size: 15, mode: 'fixed', bestEffort: false } },
    ok: ({ reply, alive }) => reply.type === 'result' && reply.words === 36 && recoverable(alive) },

  { name: 'heavy-capped-recoverable',   // 300KB « the search's ~0.5-1MB peak
    msg: { type: 'solve', input: { clues: clues(HEAVY), size: 15, mode: 'fixed', bestEffort: false }, stackLimit: 300000 },
    ok: ({ reply, alive }) => reply.type === 'error' && /stack limit|exceeded|resource_error/i.test(reply.message) && recoverable(alive) },
];

let fail = 0;
console.log(`error paths (fresh worker each; liveness = same worker solves a good request after):\n`);
for (const c of cases) {
  const r = await run(c.msg);
  const pass = c.ok(r);
  if (!pass) fail++;
  const replyStr = r.reply.type === 'error' ? `error("${r.reply.message}")` : r.reply.type === 'result' ? `result(${r.reply.words})` : r.reply.type;
  const aliveStr = r.alive.type === 'result' ? `alive(${r.alive.words})` : `NOT-ALIVE(${r.alive.type})`;
  console.log(`  ${pass ? 'PASS' : 'FAIL'}  ${c.name}`);
  console.log(`         reply=${replyStr}  after=${aliveStr}`);
}
console.log(`\n${fail === 0 ? 'error-paths OK' : `error-paths FAILED (${fail})`}`);
await browser.close();
process.exit(fail === 0 ? 0 : 1);
