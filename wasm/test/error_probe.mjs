// error_probe.mjs - gate the browser build's UNHAPPY paths at the WORKER level
// (the happy path + facade behaviour is headless.mjs). Each case runs in a
// fresh worker, then re-runs a KNOWN-GOOD request on the SAME worker
// (liveness) - proving the failure was recoverable, not a dead tab / wedged
// instance. Asserts the v1 envelope taxonomy end-to-end under the real wasm
// build:
//
//   validation:bad-seed    seed:-1 -> {status:"error", error.type:"validation"}
//                          (a failing resolver must THROW-and-map, never the
//                          old bare-fail "no layout" masquerade)
//   failure:infeasible     strict no-crossing input -> {status:"failure",
//                          detail.reason:"unplaceable_words"} - a legitimate
//                          answer, not an error
//   heavy-default          the ~38M search under the normal 256MB cap
//                          COMPLETES (36 words) - the cap doesn't harm real work
//   heavy-capped           the same search under a tiny 300KB stack cap ->
//                          {error.type:"resource_exhausted"} (the classifier's
//                          resource_error(stack) mapping), NOT a wasm-heap-
//                          growth tab abort. Only holds because the cap is
//                          applied INSIDE the forEach engine; a parent-set
//                          limit is ignored by the {engine:true} search.
//
// Prereq: server for the wasm/ tree on $URL's origin + Playwright in
// wasm/test. Run: node wasm/test/error_probe.mjs (exit 0 = all pass).

import { chromium } from 'playwright';

const URL = process.env.URL || 'http://127.0.0.1:8080/client/harness.html?noauto=1';
const TOY = ['CAT', 'CAR', 'ARC', 'RAT', 'TAR'];
const HEAVY = ['DFAD','FCC','FCB','BEED','CBED','AFCC','DED','DEF','BFF','BCD','DBED','CDF',
  'FFCC','CCD','EAB','FCF','EAA','EFAF','ABFA','BBEF','FFE','EFEE','ABCB','EFD','DACA','FAFB',
  'ACA','CBF','DEAA','AFBF','AEFD','EADF','EDDE','CEF','CADF','FDD'];
const clues = (ws) => ws.map(a => ({ answer: a }));

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();
page.on('pageerror', e => console.log(`[pageerror] ${e.message}`));
await page.goto(URL, { waitUntil: 'load' });

// Send one request message, capture its envelope, then send a KNOWN-GOOD
// request on the SAME worker and report whether it still succeeds (liveness).
async function run(msg) {
  return await page.evaluate(async ({ msg, TOY }) => {
    const good = {
      type: 'request',
      request: { v: 1, id: 'good', verb: 'arrange',
                 params: { clues: TOY.map(a => ({ answer: a })), size: 5, bestEffort: true } },
    };
    const w = new Worker('./worker.js');
    const once = (send) => new Promise((resolve) => {
      w.onmessage = (ev) => {
        const m = ev.data;
        if (m.type === 'response') resolve(m.envelope);
      };
      setTimeout(() => resolve({ status: 'timeout' }), 60000);
      w.postMessage(send);
    });
    const reply = await once(msg);
    const alive = await once(good);
    w.terminate();
    return { reply, alive };
  }, { msg, TOY });
}

const recoverable = (a) => a.status === 'success' && a.result.words.length === 5;
const req = (params, extra = {}) =>
  ({ type: 'request', request: { v: 1, id: 'probe', verb: 'arrange', params }, ...extra });

const cases = [
  { name: 'validation:bad-seed',
    msg: req({ clues: clues(TOY), size: 5, seed: -1 }),
    ok: ({ reply, alive }) => reply.status === 'error'
      && reply.error.type === 'validation'
      && !/no layout/i.test(reply.error.message)
      && recoverable(alive) },

  { name: 'failure:infeasible-strict',
    msg: req({ clues: clues(['AAAA', 'BBBB', 'CCCC']), size: 4 }),
    ok: ({ reply, alive }) => reply.status === 'failure'
      && reply.detail.reason === 'unplaceable_words'
      && recoverable(alive) },

  { name: 'heavy-default-completes',
    msg: req({ clues: clues(HEAVY), size: 15 }),
    ok: ({ reply, alive }) => reply.status === 'success'
      && reply.result.words.length === 36 && recoverable(alive) },

  { name: 'heavy-capped-recoverable',   // 300KB « the search's ~0.5-1MB peak
    msg: req({ clues: clues(HEAVY), size: 15 }, { stackLimit: 300000 }),
    ok: ({ reply, alive }) => reply.status === 'error'
      && reply.error.type === 'resource_exhausted' && recoverable(alive) },
];

let fail = 0;
console.log(`error paths (fresh worker each; liveness = same worker serves a good request after):\n`);
for (const c of cases) {
  const r = await run(c.msg);
  const pass = c.ok(r);
  if (!pass) fail++;
  const e = r.reply;
  const replyStr = e.status === 'error' ? `error(${e.error.type}: "${String(e.error.message).slice(0, 80)}")`
    : e.status === 'failure' ? `failure(${e.detail.reason})`
    : e.status === 'success' ? `success(${e.result.words.length})`
    : e.status;
  const aliveStr = r.alive.status === 'success' ? `alive(${r.alive.result.words.length})` : `NOT-ALIVE(${r.alive.status})`;
  console.log(`  ${pass ? 'PASS' : 'FAIL'}  ${c.name}`);
  console.log(`         reply=${replyStr}  after=${aliveStr}`);
}
console.log(`\n${fail === 0 ? 'error-paths OK' : `error-paths FAILED (${fail})`}`);
await browser.close();
process.exit(fail === 0 ? 0 : 1);
