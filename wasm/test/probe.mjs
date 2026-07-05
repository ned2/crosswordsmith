// probe.mjs — diagnostic harness for the low-level swipl-web load path.
//
// Drives wasm/client/probe.html (a main-thread page, NOT the Worker) to check
// library resolution and json autoload in isolation. Use this first when the
// browser path breaks — it separates "did the wasm image load + json work" from
// "is the Worker plumbing right". Same Playwright/system-Chrome setup as
// headless.mjs.
//
// Prereq: server for wasm/client/ on $URL (default http://127.0.0.1:8080/).
// Run:  node wasm/test/probe.mjs

import { chromium } from 'playwright';

const base = process.env.URL || 'http://127.0.0.1:8080/';
const browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });
const page = await browser.newPage();
const logs = [];
page.on('console', m => logs.push(`[${m.type()}] ${m.text()}`));
page.on('pageerror', e => logs.push(`[pageerror] ${e.message}`));
page.on('response', r => { if (r.status() >= 400) logs.push(`[HTTP ${r.status()}] ${r.url()}`); });
await page.goto(new URL('probe.html', base).href, { waitUntil: 'load' });
await page.waitForFunction(() => window.__done === true, { timeout: 60000 }).catch(() => logs.push('[timeout]'));
console.log('===== PROBE OUTPUT =====');
console.log(await page.textContent('#out'));
console.log('\n===== LOGS =====');
console.log(logs.join('\n') || '(none)');
await browser.close();
