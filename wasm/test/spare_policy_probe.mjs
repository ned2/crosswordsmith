// spare_policy_probe.mjs — gate the SDK's spare-worker POLICIES (payload plan
// Phase 5). The spare exists so cancellation (worker.terminate()) recovers by
// promoting a warmed standby; the sparePolicy option trades that recovery
// latency against idle workers. For EACH policy (eager | idle | lazy) this
// asserts, on the real wasm build in headless Chrome:
//
//   worker census   eager boots 2 workers; idle boots 1 and warms the spare
//                   shortly after the first response; lazy stays at 1 —
//                   observed via Playwright's page.workers() (a fresh page
//                   per policy keeps the census honest)
//   prompt cancel   aborting an in-flight heavy solve rejects with
//                   CancelledError promptly (well under the solve's runtime)
//   queued work     a request QUEUED behind the cancelled one still resolves
//                   success — under lazy this is the "a queued request
//                   requires a replacement" spawn path
//   lazy drain      lazy only: a cancel with nothing queued leaves ZERO
//                   workers (nothing respawns), and the next request
//                   cold-boots one and succeeds
//   dispose         terminates every worker (census → 0)
//   validation      an unknown policy throws TypeError at create()
//
// It also MEASURES abort→next-result recovery latency per policy
// (informational, printed not asserted — the documented numbers in
// wasm/README.md come from here).
//
// Prereq: server for the wasm/ tree on $URL's origin + Playwright in
// wasm/test. Run: node wasm/test/spare_policy_probe.mjs (exit 0 = all pass).

import { chromium } from 'playwright';
import { HEAVY } from './heavy_words.mjs';

const URL = process.env.URL || 'http://127.0.0.1:8080/client/harness.html?noauto=1';

const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });

const fails = [];
const check = (cond, name, detail = '') => {
  console.log(`  ${cond ? 'PASS' : 'FAIL'}  ${name}${detail ? `  ${detail}` : ''}`);
  if (!cond) fails.push(name);
};

// Poll the page's dedicated-worker census until it hits `want` (worker
// creation and terminate cleanup are both async) — returns the final count.
async function workersSettle(page, want, ms = 5000) {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    if (page.workers().length === want) return want;
    await new Promise(r => setTimeout(r, 100));
  }
  return page.workers().length;
}

for (const policy of ['eager', 'idle', 'lazy']) {
  console.log(`\npolicy: ${policy}`);
  const page = await browser.newPage();
  page.on('pageerror', e => console.log(`[pageerror] ${e.message}`));
  await page.goto(URL, { waitUntil: 'load' });

  // --- boot + worker census ---------------------------------------------------
  await page.evaluate(async ({ policy }) => {
    window.cw = await window.harness.createCrosswordsmith({ sparePolicy: policy });
  }, { policy });
  // eager: active+spare from create(); idle: the create() capabilities
  // round-trip is the first response, so the spare warms right after, via
  // requestIdleCallback (bounded 2s); lazy: active only, forever.
  const bootWant = policy === 'lazy' ? 1 : 2;
  const boot = await workersSettle(page, bootWant);
  check(boot === bootWant, `${policy}: boot census`, `${boot} worker(s), want ${bootWant}`);

  // --- prompt cancel + queued work proceeds + recovery latency ----------------
  const r = await page.evaluate(async ({ HEAVY }) => {
    const toy = {
      clues: [{ answer: 'CAT' }, { answer: 'CAR' }, { answer: 'ARC' },
              { answer: 'RAT' }, { answer: 'TAR' }],
      size: 5, bestEffort: true,
    };
    const heavy = { clues: HEAVY.map(a => ({ answer: a })), size: 15 };
    const cw = window.cw;
    const out = {};

    // in-flight heavy + a QUEUED toy behind it, then abort the heavy
    const ac = new AbortController();
    const heavyP = cw.arrange(heavy, { signal: ac.signal });
    const queuedP = cw.arrange(toy);                 // queued (single-flight)
    await new Promise(r => setTimeout(r, 300));      // let the heavy engage
    const t0 = performance.now();
    ac.abort();
    out.cancel = await heavyP.then(
      () => ({ how: 'resolved' }),
      (e) => ({ how: 'rejected', name: e.name, ms: performance.now() - t0 }));
    const queued = await queuedP;                    // must proceed regardless
    out.queued = { status: queued.status };
    out.recoveryMs = performance.now() - t0;         // abort → queued result
    return out;
  }, { HEAVY });
  check(r.cancel.how === 'rejected' && r.cancel.name === 'CancelledError',
        `${policy}: prompt CancelledError`, `${r.cancel.ms?.toFixed(1)}ms`);
  check(r.cancel.ms < 1000, `${policy}: cancel latency < 1s`, `${r.cancel.ms?.toFixed(1)}ms`);
  check(r.queued.status === 'success', `${policy}: queued request proceeds`,
        `status=${r.queued.status}`);
  console.log(`  info  ${policy}: abort -> queued-result recovery ${r.recoveryMs.toFixed(1)}ms`);

  // --- lazy drain: cancel with an empty queue leaves zero workers -------------
  if (policy === 'lazy') {
    const drained = await page.evaluate(async ({ HEAVY }) => {
      const heavy = { clues: HEAVY.map(a => ({ answer: a })), size: 15 };
      const ac = new AbortController();
      const p = window.cw.arrange(heavy, { signal: ac.signal });
      await new Promise(r => setTimeout(r, 300));
      ac.abort();
      return await p.then(() => 'resolved', (e) => e.name);
    }, { HEAVY });
    check(drained === 'CancelledError', 'lazy: drain cancel rejected');
    const zero = await workersSettle(page, 0);
    check(zero === 0, 'lazy: no respawn after cancel with empty queue',
          `${zero} worker(s)`);
    const cold = await page.evaluate(async () => {
      const t0 = performance.now();
      const env = await window.cw.arrange({
        clues: [{ answer: 'CAT' }, { answer: 'CAR' }, { answer: 'ARC' },
                { answer: 'RAT' }, { answer: 'TAR' }],
        size: 5, bestEffort: true,
      });
      return { status: env.status, ms: performance.now() - t0 };
    });
    check(cold.status === 'success', 'lazy: next request cold-boots a worker',
          `${cold.ms.toFixed(1)}ms (full engine boot)`);
  }

  // --- dispose terminates everything -------------------------------------------
  await page.evaluate(() => window.cw.dispose());
  const gone = await workersSettle(page, 0);
  check(gone === 0, `${policy}: dispose census`, `${gone} worker(s)`);

  await page.close();
}

// --- option validation ---------------------------------------------------------
const vPage = await browser.newPage();
await vPage.goto(URL, { waitUntil: 'load' });
const bad = await vPage.evaluate(async () => {
  try { await window.harness.createCrosswordsmith({ sparePolicy: 'zealous' }); return 'created'; }
  catch (e) { return e.constructor.name; }
});
console.log('\noption validation:');
check(bad === 'TypeError', 'unknown sparePolicy throws TypeError', bad);
await vPage.close();

await browser.close();
if (fails.length) {
  console.error(`\nFAILED: ${fails.join(', ')}`);
  process.exit(1);
}
console.log('\nall spare-policy checks passed');
