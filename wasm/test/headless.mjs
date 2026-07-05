// headless.mjs — end-to-end browser test for the crosswordsmith SDK facade.
//
// Drives the real harness page (wasm/client/harness.html — the first
// createCrosswordsmith() consumer) in headless Chrome (system Chrome via
// Playwright's `channel: 'chrome'` — no browser download):
//
//   smoke        the page's own auto-solve renders a placed 5×5 grid
//   envelope     id/verb/v echo + a typed validation error (never the
//                generic "no layout" masquerade)
//   determinism  seed:42 twice through the same facade → identical results,
//                seed provenance in diagnostics
//   concurrency  two overlapping arrange() calls each resolve to THEIR OWN
//                result (id-correlated single-flight queue)
//   cancel       abort an in-flight heavy solve → prompt CancelledError +
//                the promoted warm spare serves the next request; abort a
//                QUEUED request → only it rejects, the in-flight one lands
//   capabilities engine-sourced verb list
//
// Prereq: a server for the wasm/ tree on $URL's origin (sdk/ and client/ are
// siblings, so the web root must be wasm/):
//   ( cd wasm && python3 -m http.server 8080 )
// and Playwright installed here (see wasm/test/README.md):
//   ( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm install )
// Run:  node wasm/test/headless.mjs

import { chromium } from 'playwright';

const url = process.env.URL || 'http://127.0.0.1:8080/client/harness.html';
// fixtures/ladder_15x15_36w.pl's word set: ~38.3M inferences, several seconds
// under wasm — long enough to cancel mid-flight.
const HEAVY = ['DFAD','FCC','FCB','BEED','CBED','AFCC','DED','DEF','BFF','BCD','DBED','CDF',
  'FFCC','CCD','EAB','FCF','EAA','EFAF','ABFA','BBEF','FFE','EFEE','ABCB','EFD','DACA','FAFB',
  'ACA','CBF','DEAA','AFBF','AEFD','EADF','EDDE','CEF','CADF','FDD'];

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();

const logs = [];
page.on('console', m => logs.push(`[console.${m.type()}] ${m.text()}`));
page.on('pageerror', e => logs.push(`[pageerror] ${e.message}`));
page.on('requestfailed', r => logs.push(`[reqfail] ${r.url()} :: ${r.failure()?.errorText}`));
page.on('response', r => { if (r.status() >= 400) logs.push(`[HTTP ${r.status()}] ${r.url()}`); });

const fails = [];
const check = (cond, name, detail = '') => {
  console.log(`  ${cond ? 'PASS' : 'FAIL'}  ${name}${detail ? `  ${detail}` : ''}`);
  if (!cond) fails.push(name);
};

console.log(`navigating to ${url} …`);
await page.goto(url, { waitUntil: 'load' });

// --- smoke: the page's own auto-solve (first load pays wasm + qlf) -----------
await page.waitForFunction(() => {
  const s = document.getElementById('status')?.textContent || '';
  return /placed|error|failure/.test(s);
}, { timeout: 120000 });
const status = await page.textContent('#status');
const grid = await page.textContent('#grid');
console.log('\nsmoke (harness auto-solve):');
check(/placed 5 words on 5×5/.test(status), 'smoke solve', JSON.stringify(status));
check(grid.trim().length > 0, 'grid rendered');

// --- facade scenarios through window.harness.cw ------------------------------
const r = await page.evaluate(async ({ HEAVY }) => {
  const toy = (extra = {}) => ({
    clues: [{ answer: 'CAT' }, { answer: 'CAR' }, { answer: 'ARC' },
            { answer: 'RAT' }, { answer: 'TAR' }],
    size: 5, ...extra,
  });
  const heavy = { clues: HEAVY.map(a => ({ answer: a })), size: 15 };
  const cw = await window.harness.boot();
  const out = {};

  // envelope echo + typed validation
  const ok = await cw.arrange(toy());
  out.echo = { status: ok.status, verb: ok.verb, v: ok.v, idIsString: typeof ok.id === 'string' };
  const bad = await cw.arrange(toy({ seed: -1 }));
  out.validation = {
    status: bad.status,
    type: bad.status === 'error' ? bad.error.type : null,
    masquerade: bad.status === 'error' && /no layout/i.test(bad.error.message),
  };

  // determinism + seed provenance through the same facade instance
  const s1 = await cw.arrange(toy({ seed: 42 }));
  const s2 = await cw.arrange(toy({ seed: 42 }));
  out.determinism = {
    bothSuccess: s1.status === 'success' && s2.status === 'success',
    identical: JSON.stringify(s1.result) === JSON.stringify(s2.result),
    provenance: s1.status === 'success' && s1.result.diagnostics.arrange.seed === 42,
  };

  // concurrency: overlapping calls, distinguishable results, routed by id
  const [c1, c2] = await Promise.all([
    cw.arrange(toy()),                    // 5×5
    cw.arrange(toy({ size: 7 })),         // 7×7
  ]);
  out.concurrency = {
    bothSuccess: c1.status === 'success' && c2.status === 'success',
    ownResults: c1.result?.gridLength === 5 && c2.result?.gridLength === 7,
    distinctIds: c1.id !== c2.id,
  };

  // cancel a QUEUED request: heavy in-flight, toy queued+aborted; the heavy
  // one must still land (only the aborted promise dies)
  {
    const ctrl = new AbortController();
    const heavyP = cw.arrange(heavy);
    const queuedP = cw.arrange(toy(), { signal: ctrl.signal });
    setTimeout(() => ctrl.abort(), 50);
    const queued = await queuedP.then(
      () => ({ rejected: false }),
      (e) => ({ rejected: true, name: e.name }));
    const heavyRes = await heavyP;
    out.cancelQueued = {
      queuedCancelled: queued.rejected && queued.name === 'CancelledError',
      heavyLanded: heavyRes.status === 'success' && heavyRes.result.words.length === 36,
    };
  }

  // cancel IN-FLIGHT: heavy solve aborted mid-search → prompt rejection +
  // the promoted spare serves the follow-up
  {
    const ctrl = new AbortController();
    const heavyP = cw.arrange(heavy, { signal: ctrl.signal });
    await new Promise(res => setTimeout(res, 300));   // let it get going
    const t0 = performance.now();
    ctrl.abort();
    const inflight = await heavyP.then(
      () => ({ rejected: false }),
      (e) => ({ rejected: true, name: e.name, ms: performance.now() - t0 }));
    const after = await cw.arrange(toy());
    out.cancelInflight = {
      cancelled: inflight.rejected && inflight.name === 'CancelledError',
      promptMs: Math.round(inflight.ms ?? -1),
      recovered: after.status === 'success' && after.result.words.length === 5,
    };
  }

  // capabilities: engine-sourced
  const caps = await cw.capabilities();
  out.capabilities = {
    hasArrange: caps.verbs.includes('arrange'),
    swipl: caps.engine.swipl,
  };

  // randomSeed hygiene: integers in [0, 1e9]
  out.seedHygiene = Array.from({ length: 1000 }, window.harness.randomSeed)
    .every(n => Number.isInteger(n) && n >= 0 && n <= 1_000_000_000);

  return out;
}, { HEAVY });

console.log('\nenvelope:');
check(r.echo.status === 'success' && r.echo.verb === 'arrange' && r.echo.v === 1 && r.echo.idIsString,
      'id/verb/v echo', JSON.stringify(r.echo));
check(r.validation.status === 'error' && r.validation.type === 'validation' && !r.validation.masquerade,
      'typed validation (no masquerade)', JSON.stringify(r.validation));

console.log('\ndeterminism:');
check(r.determinism.bothSuccess && r.determinism.identical, 'seed:42 ×2 identical');
check(r.determinism.provenance, 'seed provenance in diagnostics');

console.log('\nconcurrency:');
check(r.concurrency.bothSuccess && r.concurrency.ownResults && r.concurrency.distinctIds,
      'overlapping calls get their own results', JSON.stringify(r.concurrency));

console.log('\ncancel:');
check(r.cancelQueued.queuedCancelled && r.cancelQueued.heavyLanded,
      'queued abort: only the aborted promise dies', JSON.stringify(r.cancelQueued));
check(r.cancelInflight.cancelled && r.cancelInflight.recovered,
      'in-flight abort: CancelledError + spare recovers',
      `(rejected in ${r.cancelInflight.promptMs}ms)`);
check(r.cancelInflight.promptMs >= 0 && r.cancelInflight.promptMs < 1500,
      'in-flight abort is prompt', `${r.cancelInflight.promptMs}ms`);

console.log('\ncapabilities:');
check(r.capabilities.hasArrange, 'engine-sourced verb list', JSON.stringify(r.capabilities));

console.log('\nseed hygiene:');
check(r.seedHygiene, 'randomSeed() ∈ [0, 1e9] integers');

if (logs.length) { console.log('\n===== PAGE LOGS ====='); console.log(logs.join('\n')); }
await browser.close();

if (fails.length) {
  console.log(`\nheadless FAILED (${fails.length}): ${fails.join(', ')}`);
  process.exit(1);
}
console.log('\nheadless OK');
process.exit(0);
