// headless.mjs — end-to-end browser test for the crosswordsmith WASM client.
//
// Drives the real page in headless Chrome (system Chrome via Playwright's
// `channel: 'chrome'` — no browser download): load, click Solve, wait for the
// grid. Exits 0 iff the status line reports a placed layout. This is the
// regression guard for the browser path — re-run it after any swipl bump or
// change under wasm/client/.
//
// Prereq: a server for wasm/client/ on $URL (default http://127.0.0.1:8080/):
//   ( cd wasm/client && python3 -m http.server 8080 )
// and Playwright installed here (see wasm/test/README.md):
//   ( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm install )
// Run:  node wasm/test/headless.mjs

import { chromium } from 'playwright';

const url = process.env.URL || 'http://127.0.0.1:8080/';
const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();

const logs = [];
page.on('console', m => logs.push(`[console.${m.type()}] ${m.text()}`));
page.on('pageerror', e => logs.push(`[pageerror] ${e.message}`));
page.on('requestfailed', r => logs.push(`[reqfail] ${r.url()} :: ${r.failure()?.errorText}`));
page.on('response', r => { if (r.status() >= 400) logs.push(`[HTTP ${r.status()}] ${r.url()}`); });

console.log(`navigating to ${url} …`);
await page.goto(url, { waitUntil: 'load' });
await page.waitForSelector('#solve');

console.log('clicking Solve (first click pays wasm instantiate + qlf load)…');
await page.click('#solve');

let terminal = 'TIMEOUT';
try {
  await page.waitForFunction(() => {
    const s = document.getElementById('status')?.textContent || '';
    return /placed|error|stopped/.test(s);
  }, { timeout: 90000 });
  terminal = 'REACHED';
} catch {
  logs.push('[timeout] #status never reached a terminal state within 90s');
}

const status = await page.textContent('#status').catch(() => '(no #status)');
const grid = await page.textContent('#grid').catch(() => '(no #grid)');

console.log('\n===== RESULT =====');
console.log('terminal:', terminal);
console.log('STATUS :', JSON.stringify(status));
console.log('GRID   :\n' + grid);
console.log('\n===== PAGE LOGS =====');
console.log(logs.length ? logs.join('\n') : '(none)');

await browser.close();
process.exit(/placed/.test(status) ? 0 : 1);
