// payload_bench.mjs - the payload/request/startup benchmark for the browser
// bundle (docs/plans/wasm-payload-performance.md §4, Phase 0).
//
// Measures, against the CURRENTLY STAGED wasm/client/ artifact set:
//
//   sizes      raw / gzip-9 / brotli-11 bytes of every shipped artifact
//              (worker.js, build-manifest.json, hashed js/wasm/data/qlf),
//              plus the engine-trio and shipped-set aggregates
//   requests   actual per-URL request counts observed through a local
//              cache-aware counting server, for the scenario matrix
//              {sdk-eager, raw-single-worker} x {no-cache-headers, immutable},
//              and the modeled Brotli cold-transfer for each scenario
//              (sum over observed requests of the artifact's brotli-11 size)
//   timing     fresh-Chrome-process-per-sample (cold HTTP + compiled-wasm
//              cache) medians for: engine readiness (worker construction to
//              the capabilities envelope) and the first small arrange after
//              readiness - for both the SDK facade (eager spare) and one raw
//              worker
//   memory     page-side performance.memory after the first arrange,
//              best-effort and INFORMATIONAL ONLY (Chrome's page JS-heap
//              numbers are quantized and exclude worker heaps)
//   provenance buildId, swipl commit, package list, toolchain pins (from
//              build-manifest.json) + machine/chrome/node versions
//
// Wall-clock numbers are machine-specific and reporting-only - variance is
// not yet characterised, so nothing here gates CI time (plan §4). What IS
// deterministic and gateable now: artifact names, raw sizes, and sha256s.
// `--check` verifies exactly those against the committed baseline.
//
// Prereq: staged + stamped wasm/client/ (run wasm/test/run_all.sh or
// wasm/build/build-wasm.sh first) and Playwright in wasm/test (npm ci).
//
// Run (from anywhere; paths are script-relative):
//   node wasm/test/payload_bench.mjs                 # measure + human report
//   node wasm/test/payload_bench.mjs --samples 20    # more timing samples
//   node wasm/test/payload_bench.mjs --json          # full JSON on stdout
//   node wasm/test/payload_bench.mjs --record        # + write payload-baseline.json
//   node wasm/test/payload_bench.mjs --sizes-only    # no browser (sizes table only)
//   node wasm/test/payload_bench.mjs --check         # gate sizes vs the baseline
//                                                    # (deterministic, no browser)

import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import zlib from 'node:zlib';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const WASM_ROOT = path.resolve(HERE, '..');            // web root: sdk/ + client/
const CLIENT = path.join(WASM_ROOT, 'client');
const BASELINE_PATH = path.join(HERE, 'payload-baseline.json');

// --- args ---------------------------------------------------------------------
const argv = process.argv.slice(2);
const has = (f) => argv.includes(f);
const SAMPLES = (() => {
  const at = argv.indexOf('--samples');
  return at >= 0 ? Math.max(1, parseInt(argv[at + 1], 10) || 0) : 10;
})();
const MODE_CHECK = has('--check');
const MODE_RECORD = has('--record');
const MODE_JSON = has('--json');
const SIZES_ONLY = has('--sizes-only') || MODE_CHECK;

const fail = (msg) => { console.error(`payload_bench FAILED: ${msg}`); process.exit(1); };

// --- the shipped set (manifest-driven) ------------------------------------------
// Everything a cold browser deploy serves for the engine, in the plan's table
// order. Hashed names come from the manifest so the bench always measures the
// stamped set the worker actually loads.
const manifestPath = path.join(CLIENT, 'build-manifest.json');
if (!existsSync(manifestPath)) {
  fail('wasm/client/build-manifest.json missing - stage + stamp first (wasm/test/run_all.sh)');
}
const manifest = JSON.parse(readFileSync(manifestPath, 'utf8'));

const shipped = [
  { name: 'worker.js', file: 'worker.js' },
  { name: 'build-manifest.json', file: 'build-manifest.json' },
  ...Object.entries(manifest.artifacts).map(([name, a]) => ({
    name, file: a.hashed, sha256: a.sha256, bytes: a.bytes,
  })),
];

for (const s of shipped) {
  const p = path.join(CLIENT, s.file);
  if (!existsSync(p)) fail(`${s.file} missing from wasm/client/ - restage (wasm/test/run_all.sh)`);
  if (s.sha256) {
    const digest = createHash('sha256').update(readFileSync(p)).digest('hex');
    if (digest !== s.sha256) {
      fail(`${s.file} does not match its manifest sha256 - the staged set is inconsistent; restamp`);
    }
  }
}

// --- sizes ----------------------------------------------------------------------
// gzip-9 / brotli-11 are diagnostic comparators (plan §2): a production host
// will not emit byte-identical streams, but the levels are fixed so the
// numbers are comparable run-to-run.
console.error('measuring artifact sizes (brotli-11 on the wasm takes a few seconds)...');
const sizes = {};
for (const s of shipped) {
  const buf = readFileSync(path.join(CLIENT, s.file));
  sizes[s.name] = {
    file: s.file,
    raw: buf.length,
    gzip: zlib.gzipSync(buf, { level: 9 }).length,
    brotli: zlib.brotliCompressSync(buf, {
      params: { [zlib.constants.BROTLI_PARAM_QUALITY]: 11 },
    }).length,
  };
}
const sum = (names, k) => names.reduce((t, n) => t + sizes[n][k], 0);
const TRIO = ['swipl-web.js', 'swipl-web.wasm', 'swipl-web.data'];
const ALL = shipped.map((s) => s.name);
const aggregates = {
  engineTrio: { raw: sum(TRIO, 'raw'), gzip: sum(TRIO, 'gzip'), brotli: sum(TRIO, 'brotli') },
  shippedSet: { raw: sum(ALL, 'raw'), gzip: sum(ALL, 'gzip'), brotli: sum(ALL, 'brotli') },
};

// --- --check: the deterministic gate ---------------------------------------------
// Artifact membership + raw bytes + manifest-sha consistency (verified above)
// are reproducible for identical build inputs; wall-clock and request counts
// are not, so only sizes gate (plan §4).
if (MODE_CHECK) {
  if (!existsSync(BASELINE_PATH)) fail(`no baseline at ${BASELINE_PATH} - record one first (--record)`);
  const base = JSON.parse(readFileSync(BASELINE_PATH, 'utf8'));
  const problems = [];
  for (const [name, b] of Object.entries(base.sizes)) {
    if (!sizes[name]) { problems.push(`baseline artifact missing from staged set: ${name}`); continue; }
    if (sizes[name].raw !== b.raw) {
      problems.push(`${name}: raw ${sizes[name].raw} bytes != baseline ${b.raw}`);
    }
  }
  for (const name of Object.keys(sizes)) {
    if (!base.sizes[name]) problems.push(`staged artifact not in baseline: ${name}`);
  }
  if (problems.length) {
    console.error('payload size drift vs the committed baseline (intentional? re---record and review the diff):');
    for (const p of problems) console.error(`  ${p}`);
    process.exit(1);
  }
  console.log(`payload sizes match the baseline (${Object.keys(base.sizes).length} artifacts, buildId ${base.build.buildId.slice(0, 12)}...)`);
  process.exit(0);
}

// --- counting static server -------------------------------------------------------
const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.json': 'application/json',
  '.wasm': 'application/wasm',
  '.data': 'application/octet-stream',
  '.qlf': 'application/octet-stream',
  '.pl': 'text/plain; charset=utf-8',
};
const HASHED_RE = /\.[0-9a-f]{12}\.(js|wasm|data|qlf)$/;

function startServer(cacheMode /* 'none' | 'immutable' */) {
  const requests = [];   // { path, status } in arrival order
  const server = http.createServer((req, res) => {
    const urlPath = decodeURIComponent(new URL(req.url, 'http://x').pathname);
    const filePath = path.join(WASM_ROOT, path.normalize(urlPath).replace(/^([/\\])+/, ''));
    if (!filePath.startsWith(WASM_ROOT) || !existsSync(filePath)) {
      requests.push({ path: urlPath, status: 404 });
      res.writeHead(404); res.end('not found'); return;
    }
    requests.push({ path: urlPath, status: 200 });
    const headers = { 'Content-Type': MIME[path.extname(filePath)] || 'application/octet-stream' };
    if (cacheMode === 'immutable') {
      // The plan's "correctly configured" model: hashed artifacts are immutable
      // for a year; worker.js + the manifest revalidate every deploy. No
      // validators (ETag/Last-Modified) are sent, so a cold profile never
      // sees a 304 - every observed request is a full transfer.
      headers['Cache-Control'] = HASHED_RE.test(filePath)
        ? 'public, max-age=31536000, immutable'
        : 'no-cache';
    }
    res.writeHead(200, headers);
    res.end(readFileSync(filePath));
  });
  return new Promise((resolve) => {
    server.listen(0, '127.0.0.1', () => {
      resolve({ server, port: server.address().port, requests });
    });
  });
}

// --- browser drivers ---------------------------------------------------------------
// Playwright lives in wasm/test/node_modules (committed lockfile) - imported
// lazily so --sizes-only/--check need no browser toolchain. System Chrome,
// one FRESH browser process per sample: a new launch gets a fresh temp
// profile, so the HTTP cache and the compiled-wasm code cache are cold -
// without this, sample 2+ silently measures the repeat-visit path (plan §4).
let chromium = null;
const launchBrowser = () => chromium.launch({ channel: 'chrome', headless: true, args: ['--no-sandbox'] });

const TOY = ['CAT', 'CAR', 'ARC', 'RAT', 'TAR'];

// One cold sample: navigate the harness host page, boot the engine, time
// readiness + the first small arrange. `sdk` drives createCrosswordsmith()
// (active + eager warmed spare; readiness = the facade's init capabilities
// round-trip resolving). `raw` constructs ONE bare worker and speaks the
// envelope RPC directly (readiness = the capabilities envelope arriving).
async function coldSample(page, mode) {
  return await page.evaluate(async ({ mode, TOY }) => {
    const toy = { clues: TOY.map((a) => ({ answer: a })), size: 5 };
    const out = {};
    if (mode === 'sdk') {
      const t0 = performance.now();
      const cw = await window.harness.boot();
      out.readyMs = performance.now() - t0;
      const t1 = performance.now();
      const res = await cw.arrange(toy);
      out.arrangeMs = performance.now() - t1;
      out.ok = res.status === 'success' && res.result.words.length === 5;
    } else {
      const w = new Worker('./worker.js');
      const once = (msg) => new Promise((resolve) => {
        w.onmessage = (ev) => { if (ev.data.type === 'response') resolve(ev.data.envelope); };
        w.postMessage(msg);
      });
      const req = (id, verb, params) =>
        ({ type: 'request', request: { v: 1, id, verb, params } });
      const t0 = performance.now();
      const caps = await once(req('r-1', 'capabilities', {}));
      out.readyMs = performance.now() - t0;
      const t1 = performance.now();
      const res = await once(req('r-2', 'arrange', toy));
      out.arrangeMs = performance.now() - t1;
      out.ok = caps.status === 'success'
        && res.status === 'success' && res.result.words.length === 5;
      w.terminate();
    }
    // Best-effort, informational only: Chrome's page JS heap is quantized and
    // does NOT include the workers' wasm heaps.
    out.pageJsHeap = performance.memory
      ? { used: performance.memory.usedJSHeapSize, total: performance.memory.totalJSHeapSize }
      : null;
    return out;
  }, { mode, TOY });
}

async function runScenario(mode, cacheMode, settleMs = 1500) {
  const { server, port, requests } = await startServer(cacheMode);
  const browser = await launchBrowser();
  try {
    const page = await browser.newPage();
    await page.goto(`http://127.0.0.1:${port}/client/harness.html?noauto=1`, { waitUntil: 'load' });
    const sample = await coldSample(page, mode);
    if (!sample.ok) fail(`scenario ${mode}/${cacheMode}: engine did not produce a valid toy solve`);
    // Let stragglers land (the eager spare's fetches, the known duplicate qlf
    // probe) before freezing the count.
    await new Promise((r) => setTimeout(r, settleMs));
    const counts = {};
    for (const r of requests) counts[r.path] = (counts[r.path] || 0) + 1;
    // Modeled cold transfer: every observed request for a shipped artifact at
    // its brotli-11 size (the plan's §2 model - aborted duplicates count full,
    // matching how the baseline figures were derived).
    let modeled = 0;
    for (const r of requests) {
      const name = shipped.find((s) => r.path === `/client/${s.file}`)?.name;
      if (name) modeled += sizes[name].brotli;
    }
    return { counts, modeledBrotli: modeled, sample };
  } finally {
    await browser.close();
    server.close();
  }
}

const stats = (xs) => {
  const s = [...xs].sort((a, b) => a - b);
  const q = (p) => s[Math.min(s.length - 1, Math.ceil(p * s.length) - 1)];
  return {
    n: s.length,
    median: s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2,
    mean: s.reduce((a, b) => a + b, 0) / s.length,
    min: s[0], max: s[s.length - 1], p95: q(0.95),
  };
};
const round1 = (o) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, Math.round(v * 10) / 10]));

// --- measure -------------------------------------------------------------------------
let scenarios = null;
let timing = null;
let memory = null;
let chromeVersion = null;

if (!SIZES_ONLY) {
  ({ chromium } = await import('playwright'));

  // request-count matrix (one cold run per cell)
  scenarios = {};
  for (const mode of ['sdk-eager', 'raw-single']) {
    for (const cacheMode of ['no-cache-headers', 'immutable']) {
      console.error(`requests: ${mode} / ${cacheMode} ...`);
      const m = mode === 'sdk-eager' ? 'sdk' : 'raw';
      const { counts, modeledBrotli } = await runScenario(m, cacheMode === 'immutable' ? 'immutable' : 'none');
      scenarios[`${mode}/${cacheMode}`] = { requests: counts, modeledColdBrotli: modeledBrotli };
    }
  }

  // timing: SAMPLES fresh-Chrome samples per mode, alternating so slow drift
  // (thermal, background load) hits both modes evenly. Immutable headers -
  // the realistic production config; localhost makes cache mode a non-factor
  // for wall-clock anyway.
  const { server, port } = await startServer('immutable');
  const samples = { sdk: [], raw: [] };
  const heaps = [];
  try {
    for (let i = 0; i < SAMPLES; i++) {
      for (const mode of ['sdk', 'raw']) {
        const browser = await launchBrowser();
        if (!chromeVersion) chromeVersion = browser.version();
        try {
          const page = await browser.newPage();
          await page.goto(`http://127.0.0.1:${port}/client/harness.html?noauto=1`, { waitUntil: 'load' });
          const s = await coldSample(page, mode);
          if (!s.ok) fail(`timing sample ${mode} #${i + 1}: engine did not produce a valid toy solve`);
          samples[mode].push(s);
          if (mode === 'sdk' && s.pageJsHeap) heaps.push(s.pageJsHeap.used);
        } finally {
          await browser.close();
        }
        console.error(`timing: sample ${i + 1}/${SAMPLES} ${mode} ready=${samples[mode].at(-1).readyMs.toFixed(1)}ms arrange=${samples[mode].at(-1).arrangeMs.toFixed(1)}ms`);
      }
    }
  } finally {
    server.close();
  }
  timing = {
    samplesPerMode: SAMPLES,
    coldProfilePerSample: true,
    cacheHeaders: 'immutable',
    sdkEager: {
      readinessMs: round1(stats(samples.sdk.map((s) => s.readyMs))),
      firstArrangeMs: round1(stats(samples.sdk.map((s) => s.arrangeMs))),
    },
    rawSingleWorker: {
      readinessMs: round1(stats(samples.raw.map((s) => s.readyMs))),
      firstArrangeMs: round1(stats(samples.raw.map((s) => s.arrangeMs))),
    },
  };
  memory = heaps.length
    ? { note: 'page JS heap after first arrange (sdk mode); quantized, EXCLUDES worker wasm heaps - informational only',
        usedJsHeapBytes: round1(stats(heaps)) }
    : { note: 'performance.memory unavailable', usedJsHeapBytes: null };
}

// --- report ----------------------------------------------------------------------------
const report = {
  v: 1,
  recorded: new Date().toISOString(),
  machine: {
    platform: `${os.platform()} ${os.release()}`,
    cpu: os.cpus()[0]?.model || 'unknown',
    node: process.version,
    chrome: chromeVersion,
  },
  build: {
    buildId: manifest.buildId,
    swiplCommit: manifest.swipl.commit,
    packages: Object.keys(manifest.swipl.submodules).sort(),
    toolchain: manifest.toolchain,
  },
  sizes,
  aggregates,
  scenarios,
  timing,
  memory,
};

if (MODE_RECORD) {
  writeFileSync(BASELINE_PATH, JSON.stringify(report, null, 2) + '\n');
  console.error(`baseline written: ${BASELINE_PATH}`);
}

if (MODE_JSON) {
  console.log(JSON.stringify(report, null, 2));
} else {
  const kb = (b) => `${(b / 1024).toFixed(1)} KB`;
  console.log(`\nbuild ${manifest.buildId.slice(0, 12)}... (swipl ${manifest.swipl.commit.slice(0, 9)}, emsdk ${manifest.toolchain.emsdk})`);
  console.log('\nartifact sizes (bytes; gzip-9 / brotli-11 diagnostic comparators):');
  const w = Math.max(...ALL.map((n) => n.length));
  console.log(`  ${'artifact'.padEnd(w)}  ${'raw'.padStart(10)}  ${'gzip'.padStart(10)}  ${'brotli'.padStart(10)}`);
  for (const n of ALL) {
    const s = sizes[n];
    console.log(`  ${n.padEnd(w)}  ${String(s.raw).padStart(10)}  ${String(s.gzip).padStart(10)}  ${String(s.brotli).padStart(10)}`);
  }
  for (const [label, a] of Object.entries(aggregates)) {
    console.log(`  ${label.padEnd(w)}  ${String(a.raw).padStart(10)}  ${String(a.gzip).padStart(10)}  ${String(a.brotli).padStart(10)}`);
  }
  if (scenarios) {
    console.log('\ncold request scenarios (observed counts; modeled transfer = requests x brotli-11):');
    for (const [name, sc] of Object.entries(scenarios)) {
      console.log(`  ${name}: ${kb(sc.modeledColdBrotli)} modeled`);
      for (const [p, c] of Object.entries(sc.requests).sort()) {
        if (p.startsWith('/client/') && !p.endsWith('.html')) console.log(`      ${c}x ${p}`);
      }
    }
  }
  if (timing) {
    console.log(`\ntiming (fresh Chrome process per sample, n=${SAMPLES} per mode, immutable headers):`);
    for (const [label, t] of [['sdk-eager', timing.sdkEager], ['raw-single', timing.rawSingleWorker]]) {
      console.log(`  ${label}: readiness median ${t.readinessMs.median}ms (min ${t.readinessMs.min}, max ${t.readinessMs.max}); first arrange median ${t.firstArrangeMs.median}ms`);
    }
  }
  if (memory?.usedJsHeapBytes) {
    console.log(`\npage JS heap after first arrange (informational): median ${kb(memory.usedJsHeapBytes.median)}`);
  }
}
