# Benchmark and probe tooling inventory (2026-07)

Status: PHASE 0 COMPLETE (2026-07-17). This inventory was taken from
`cleanup/benchmark-probe-tooling-review` at
`1119899358a196a86d792d057df33050c4dd61d1` before correctness changes. It is a
classification record, not authorization to perform the listed removals.

## 1. Scope and method

The inventory starts from `git ls-files benchmarks` and gives one row to every
tracked `.pl`, `.py`, or `.sh` file. Physical line counts include comments and
blank lines and are exclusive to the named file. The 12 `cliff_*` Prolog files
contain only comments and one ground `clues/1` fact; they are classified below
as fixture data rather than executable source, but remain in the tracked-file
inventory so no `.pl` file is silently omitted.

Closeout arithmetic correction (2026-07-19): the original subtotal omitted the
202-line `run_arrange.pl` from source lines and added those lines to fixture data;
it also overstated the source and total file counts by one. The corrected source
inventory is 73 executable/source files and 11,225 exclusive lines. The 12
fixture-data files add 774 lines, for 85 tracked `.pl`/`.py`/`.sh` files and
11,999 lines total. The historical-source review ceiling reconciles as:

| Area | Lines |
|---|---:|
| `probe_arrange/*.{pl,py}` | 4,011 |
| `probe_f1/*.pl` | 572 |
| `probe_fh2/*.pl` | 539 |
| Fill MAC/B0/frontier campaign tools | 1,320 |
| `probe_backtrack.pl` | 183 |
| Top-level A-G2 launchers | 112 |
| **Historical review ceiling** | **6,737** |

Reachability abbreviations in the tables are `M` (Make target), `T` (test),
`D` (living documentation or dated report), and `S` (source caller). Evidence
keys are defined here:

| Key | Evidence location |
|---|---|
| `AB`, `AI` | `benchmarks/baseline.json`, `benchmarks/arrange_identity.sha256` |
| `GB`, `GI` | `benchmarks/greedy_baseline.json`, `benchmarks/greedy_identity.sha256` |
| `FB`, `FI` | `benchmarks/fill_baseline.json`, `benchmarks/fill_identity.sha256` |
| `P0` | `benchmarks/results/2026-07-16-p0-arrange-probes.md` |
| `D0` | `benchmarks/results/2026-07-17-p-d0-support-delta.md` |
| `C0` | `benchmarks/results/2026-07-17-p-c0-duplicate-failed-work.md` |
| `R0` | `benchmarks/results/2026-07-17-p-r0-fixed-instance-pilot.md` and its JSONL |
| `G1` | `benchmarks/results/2026-07-17-a-g1-legality-before-score-premise.md` |
| `G2` | `benchmarks/results/2026-07-17-a-g2-transpose-premise.md` and `2026-07-17-a-g2-transpose-product.md` |
| `D1` | `benchmarks/results/2026-07-17-a-d1-direct-buckets.md` |
| `D2` | `benchmarks/results/2026-07-17-a-d2-newest-source-delta.md` |
| `CO` | `benchmarks/results/2026-07-17-arrange-campaign-closeout.md` |
| `F1` | `benchmarks/results/2026-07-05-p-f1-attribution.md` |
| `FH2G`, `FH2B` | `benchmarks/results/2026-07-05-f-h2-gate-probe.md`, `2026-07-05-f-h2-bitset-count.md` |
| `FQ`, `B0` | `benchmarks/fill_quality/README.md`, `benchmarks/results/2026-07-16-fill-b0-mac-instrumentation.md` |

Each row records the requested owner, primary ownership, reachability, product
model, evidence, drift risk, shared dependencies, replacement coverage,
disposition, and exclusive line count. `REMOVE` and `EXTRACT THEN REMOVE` remain
blocked until their evidence and replacement gates are satisfied. `HOLD
AC-FILL-12` is stronger: no frontier source may be removed until a permanent
gate covers `blocked_13a` at score thresholds 30 and 1, fails on missing or
budget-exhausted fills, and enforces the recorded quality floors.

## 2. Permanent and top-level source

| Path | Lines | Owner / primary ownership | Reachability | Product model | Evidence | Drift risk / shared dependencies | Replacement coverage | Disposition |
|---|---:|---|---|---|---|---|---|---|
| `benchmarks/analyze_layout.py` | 108 | Layout analyzer / layout analysis | D: `fixtures/README.md`, cryptic layout spec | Current canonical JSON | Cryptic layout spec results | Medium; schema and stale `crossword.pl` example | No independent replacement | RETAIN; repair usage |
| `benchmarks/bench_core.pl` | 142 | Measurement core / shared benchmark infrastructure | S: arrange, greedy, fill, matrix runners; D: root README | Current shared tooling | AB, GB, FB | High; SWI statistics/process and GNU `time`; C1/C2/C6 | All permanent benches consume it | RETAIN; correctness fixes |
| `benchmarks/check_arrange_identity.sh` | 101 | Strict identity / strict arrange | M: `bench-arrange-verify`; S: workloads and identity runner | Current product | AI | Medium; shell lifecycle and manifest mechanics | Strict goldens are narrower | RETAIN |
| `benchmarks/check_baseline.pl` | 518 | Strict ratchet/storage/history / strict arrange | M: strict check/promote/record/log/history; T: `tests/benchmarks.plt` | Current product | AB | High; runner subprocess, JSONL and duplicated storage | Focused promotion tests; no ratchet replacement | RETAIN; C2/C7 |
| `benchmarks/check_fill_baseline.pl` | 487 | Fill ratchet/storage/history / fill | M: fill check/promote/record/log/history | Current product | FB | High; unsafe direct record plus duplicated process/storage | No existing focused fill persistence test | RETAIN; C2/C3/C7 |
| `benchmarks/check_fill_identity.sh` | 103 | Raw CLI identity / fill | M: `bench-fill-verify` | Current product | FI | High; partial record can replace manifest | Fill goldens and artifact identity are narrower/different | RETAIN; C5 |
| `benchmarks/check_fill_identity_artifact.sh` | 104 | Artifact CLI identity / fill | M: `bench-fill-verify` | Current product artifact path | FI, F-H2 reports | Medium; shell lifecycle and artifact flags | Fill artifact unit tests plus shared FI manifest | RETAIN |
| `benchmarks/check_greedy_baseline.pl` | 241 | Greedy ratchet/storage/history / greedy arrange | M: greedy check/promote/record/log/history; T: greedy self-test | Current product | GB | High; subprocess and duplicated storage | Self-test covers record/read-back | RETAIN; C2/C7 |
| `benchmarks/check_greedy_identity.sh` | 29 | Greedy semantic identity / greedy arrange | M: greedy identity/verify | Current product | GI | High; direct record and producer dual-pipe capture | Complete identity has no replacement | RETAIN; C2/C5 review |
| `benchmarks/fill_subjects.pl` | 109 | Fill samplers / fill | S: `run_fill.pl`; historical F1/FH2 callers | Current exported benchmark seams | FB | Medium; fresh-slot/destructive-state invariant | Fill ratchet and product tests | RETAIN |
| `benchmarks/fill_workloads.pl` | 88 | Fill workload manifest / fill | S: runner and two identity scripts | Current product | FB, FI | Medium; stale comments and artifact alignment | Baseline/manifest mirror but do not select | RETAIN; repair comments |
| `benchmarks/fixtures.pl` | 26 | Strategy manifest / strategy research | M: matrix; S: matrix/start sensitivity; D: README | Current research path | Dated matrix/start results | Medium; fixture manifest/docs synchronization | Matrix and start sweep cover different axes | RETAIN |
| `benchmarks/gen_fill_dict.py` | 89 | Frozen dictionary generator / reproducible inputs | D: fixture dictionary README | Current input generator | Recorded seeds/counts/hashes in fixture README | Medium; Python RNG and external ENABLE bytes | Built-in check and committed hashes | RETAIN |
| `benchmarks/gen_mesh_fixture.py` | 168 | Synthetic fixture generator / reproducible inputs | S: historical fixture checker; D: fixture README | Current copied legality model | Generated fixture headers | High; Python legality clone of product | Committed fixtures and strict identities | RETAIN |
| `benchmarks/gen_real_fixture.py` | 194 | Real-word fixture generator / reproducible inputs | D: root/fixture README and generated headers | Current copied legality model | Generated fixture headers | High; copied mesh/legality logic and ENABLE | Committed fixtures and strict identities | RETAIN |
| `benchmarks/greedy_subjects.pl` | 494 | Greedy samplers/identity / greedy arrange | S: greedy runner; T: greedy benchmark; historical G2 callers | Current path plus A-G1 replay | GB, GI, G1 | High; private hooks, A-G1 replay, sequential dual pipes | Ratchet, identity, current phase tests | RETAIN; remove A-G1-only pass later |
| `benchmarks/greedy_workloads.pl` | 53 | Greedy manifest / greedy arrange | S: greedy runner and historical G2 tools | Current product | GB, GI | Medium; baseline/identity synchronization | GB/GI do not replace selection owner | RETAIN |
| `benchmarks/probe_backtrack.pl` | 183 | Historical P1 museum file / historical P1 | D: `docs/experiments.md`, cleanup plan | Removed private hooks; explicitly non-runnable | Experiments log and `docs/experiments/p1-backtrack-instrumentation.patch` | Terminal; removed `probe_*` hooks | Durable report and patch | REMOVE after reference confirmation |
| `benchmarks/run_arrange.pl` | 202 | Strict product runner / strict arrange | M: `bench`; S: strict checker | Current product | AB | High; C1/C4 and duplicated setup | Ratchet executes it, no focused selection tests | RETAIN; C1/C4 |
| `benchmarks/run_arrange_g2_probe.pl` | 20 | A-G2 premise launcher / historical A-G2 | D: G2 premise only | Retired four-direct-corner replay | G2; source commit `5500ebd` | High; campaign `g2_transpose.pl` | Focused product G2 tests and GI/GB | REMOVE after extraction |
| `benchmarks/run_arrange_g2_wall.pl` | 92 | A-G2 wall launcher / historical A-G2 | D: G2 product report only | Retired paired worktree comparison | G2; source commit `7a14082` | Terminal; two historical worktrees | G2 report and current GB | REMOVE |
| `benchmarks/run_arrange_greedy.pl` | 192 | Greedy product/identity runner / greedy arrange | M: greedy bench/check/identity | Current product plus A-G1 pass | GB, GI | High; C4 and mandatory campaign replay | Permanent ratchet/identity | RETAIN; C4 and later V5 cleanup |
| `benchmarks/run_arrange_identity.pl` | 110 | Strict identity payload / strict arrange | S: strict identity shell | Current private diagnostic path | AI | Medium; seed/check-target state and private builder | Full strict identity shell | RETAIN |
| `benchmarks/run_fill.pl` | 216 | Fill product runner / fill | M: fill bench; S: fill checker | Current product | FB | High; C1/C4 and stale gate text | Fill ratchet | RETAIN; C1/C4 |
| `benchmarks/run_matrix.pl` | 133 | Strategy matrix runner / strategy research | M: matrix; D: root README | Current research path | Dated matrix CSVs | Medium; C1 and duplicated fixture setup | No replacement | RETAIN; C1 |
| `benchmarks/start_sensitivity.pl` | 85 | Start-position sweep / strategy research | D: experiments; direct command | Current research path | `2026-06-20-start-sensitivity.csv` | Medium; duplicated fixture loading | Matrix does not sweep every start | RETAIN |
| `benchmarks/subjects.pl` | 82 | Strict command/search samplers / strict arrange | S: strict runner | Current sanctioned benchmark seams | AB | Medium; C1 and exact `/6` seam | Strict ratchet/identity | RETAIN |
| `benchmarks/workloads.pl` | 115 | Strict workload manifest / strict arrange | S: strict runner/check/identity; D: README | Current product | AB, AI | Medium; stale “15x15 ladder” text | AB/AI do not replace selector | RETAIN; repair comments |

Permanent/top-level source subtotal: 28 files, 4,484 lines. The three historical
top-level files account for 295 lines; the other 25 files account for 4,189.

## 3. Arrange campaign source

All files in this section are primarily owned by the completed arrange campaign.
No source is retained merely as a museum copy. Removal remains ordered and
conditional on the assertion migrations in section 7 and report provenance.

| Path | Lines | Owner / primary ownership | Reachability | Product model | Evidence | Drift risk / shared dependencies | Replacement coverage | Disposition |
|---|---:|---|---|---|---|---|---|---|
| `benchmarks/probe_arrange/ad1_buckets.pl` | 212 | A-D1 differential / A-D1 | T: probe suite; S: A-D1/A-D2 runners | Retired assoc vs full-direct recursion | D1 | Very high; private recount and copied recursion | Core bucket tests and assoc full-tree oracle | EXTRACT/CONFIRM THEN REMOVE |
| `benchmarks/probe_arrange/ad2_delta.pl` | 240 | A-D2 differential / A-D2 | T: probe suite; S: D2/closeout tools | Copy of shipping-at-adoption A-D2 recursion | D2, CO | Very high; D0 scanner, A-D1 reference, private state | Focused core A-D2 tests, strict ratchet/identity | EXTRACT/CONFIRM THEN REMOVE |
| `benchmarks/probe_arrange/analyze_pr0.py` | 146 | P-R0 analysis / P-R0 | D: campaign README/report; S: schema | Historical fixed 131-row analysis | R0 | High; fixed base, row count, seeds and schema | Final report; raw JSONL pending data decision | REMOVE after provenance/data decision |
| `benchmarks/probe_arrange/check_fixtures.py` | 67 | Frozen campaign corpus checker / P-R0 | M: probe check; D: campaign README | Historical fixture regeneration | P0, R0 | Medium; permanent generator output can drift | Parameters and hashes in reports | REMOVE |
| `benchmarks/probe_arrange/check_seeds.py` | 61 | Frozen seed checker / P-R0 | M: probe check; D: campaign README; S: seed generator | Historical SplitMix manifest check | P0, R0 | Medium; private product SplitMix hook | Product seeded permutation tests | REMOVE |
| `benchmarks/probe_arrange/closeout_direct.pl` | 140 | Final direct attribution / closeout | S: closeout runner/test | Copy of shipping A-D2 direct recursion | CO | Very high; exact private recursion/event positions | Core A-D2 tests and strict sentinels | REMOVE |
| `benchmarks/probe_arrange/d0_support.pl` | 631 | P-D0 observer/scanner / P-D0 | T: probe suite; S: D0 and D2 tools | Retired assoc strict replay | D0 | Extreme; copied DFS, watches, private helpers | Core A-D2 tests and assoc full-tree oracle | EXTRACT/CONFIRM THEN REMOVE |
| `benchmarks/probe_arrange/g2_transpose.pl` | 396 | A-G2 four-corner oracle / A-G2 | T: greedy benchmark; S: G2 launcher | Retired four-direct-corner greedy sweep | G2 | Extreme; private scoring/order and old pool | Migrate focused product tests; GI/GB | EXTRACT THEN REMOVE |
| `benchmarks/probe_arrange/make_seed_manifest.py` | 14 | One-time seed generator / P-R0 | S: imports seed checker; no live command | Historical regeneration script | P0, R0 | Low; closed campaign partition policy | Manifest and reports | REMOVE |
| `benchmarks/probe_arrange/measure_ad2.pl` | 46 | A-D2 control launcher / A-D2 | D: D2 report | Historical full-vs-delta measurement | D2 | High; probe root, D0, A-D1/A-D2 copies | D2 report and permanent tests | REMOVE |
| `benchmarks/probe_arrange/measure_d0.pl` | 109 | P-D0 report generator / P-D0 | D: campaign README/report | Historical premise measurement | D0 | Very high; probe root and D0 replay | D0 tables/report | REMOVE |
| `benchmarks/probe_arrange/measure_overhead.pl` | 71 | Phase-0 calibration / Phase-0 campaign | D: campaign README/report | Historical authority/twin calibration | P0 | High; twin and host-specific wall | Recorded calibration | REMOVE |
| `benchmarks/probe_arrange/measure_process.py` | 34 | Phase-0 process envelope / Phase-0 campaign | D: campaign README; S: `run.pl` | Historical schema wrapper | P0 | Medium; GNU `time`, closed row schema | Permanent process samplers | REMOVE |
| `benchmarks/probe_arrange/pc0_duplicate_work.pl` | 94 | P-C0 cache measurement / P-C0 | D: C0 report | Rejected cache observer | C0 | Very high; probe root and closed counters | C0 report; no product cache invariant | REMOVE |
| `benchmarks/probe_arrange/probe_arrange.pl` | 546 | Authority/twin framework / Phase-0 campaign | T: probe suite; S: most campaign runners | Mixed current wrappers and retired assoc DFS | P0, C0 | Extreme; copied engine/private seams/dynamic tables | Product tests, strict sentinels, assoc oracle | EXTRACT THEN REMOVE last |
| `benchmarks/probe_arrange/profile.pl` | 14 | Twin profiler / Phase-0 campaign | D: campaign README | Historical twin, not current product | P0 | High; profiler attribution follows retired twin | Re-profile current product if needed | REMOVE |
| `benchmarks/probe_arrange/run.pl` | 137 | Phase-0/P-R0 row runner / Phase-0 campaign | D: campaign README; S: process/R0/closeout/schema tools | Historical positional schema over authority/twin | P0, R0, CO | High; private authority and obsolete schema | Permanent runners and dated reports | REMOVE after dependants |
| `benchmarks/probe_arrange/run_ad1_wall.pl` | 66 | A-D1 paired wall / A-D1 | D: D1 report | Historical assoc/full-direct comparison | D1 | Very high; host-specific and A-D1 copy | D1 report/current strict ratchet | REMOVE |
| `benchmarks/probe_arrange/run_ad2_wall.pl` | 68 | A-D2 paired wall / A-D2 | D: D2 report | Historical full/delta comparison | D2 | Very high; A-D1/A-D2/D0 copies | D2 report/current strict ratchet | REMOVE |
| `benchmarks/probe_arrange/run_closeout_authority.py` | 88 | 500M authority envelope / closeout | D: CO report; S: `run.pl` | Current-at-closeout authority through old schema | CO | High; hard-coded cases/budget | Strict workloads/ratchet and recorded table | REMOVE |
| `benchmarks/probe_arrange/run_closeout_direct.pl` | 73 | Direct closeout orchestrator / closeout | D: CO report; S: three copied paths | Product/copy/differential comparison | CO | Extreme; three drifting implementations | CO counters/digests and permanent gates | REMOVE |
| `benchmarks/probe_arrange/run_pr0.py` | 249 | P-R0 campaign launcher / P-R0 | D: campaign README/report | Refuses runs outside `probe/a-r0@1bccf47` | R0 | Terminal; setup base predates runner source | Report and raw JSONL | REMOVE after reconstruction recipe |
| `benchmarks/probe_arrange/schema.py` | 124 | Campaign JSONL schema / P-R0 | M: schema test; S: R0 runner/analyzer/tests | Historical authority/instrumented schema | P0, R0 | Medium; no retained producer after cleanup | Reports describe schema; raw data immutable | REMOVE with consumers |
| `benchmarks/probe_arrange/test_closeout_direct.pl` | 64 | Closeout copy verification / closeout | D: CO report; standalone source | Campaign copies against product | CO | Extreme; duplicate private recursion | Core A-D2 tests and strict sentinels | REMOVE |
| `benchmarks/probe_arrange/test_d0_support.pl` | 77 | D0 observer tests / P-D0 | M: probe check | Retired observer mechanics | D0 | Very high; probe root and D0 copy | Core A-D2 tests/assoc oracle | REMOVE after coverage confirmation |
| `benchmarks/probe_arrange/test_schema.py` | 176 | Campaign schema tests / P-R0 | M: probe check; S: `run.pl` | Obsolete schema and runner emissions | P0, R0 | High; no retained schema consumer | No replacement needed after removal | REMOVE with schema |
| `benchmarks/probe_arrange/verify_controls.pl` | 68 | Phase-0/P-C0 replay verifier / Phase-0 campaign | D: campaign README/reports | Authority vs none/lean/full twins | P0, C0 | Extreme; copied recursion/instrumentation | Product tests, assoc oracle, strict identity | REMOVE |

Arrange campaign source subtotal: 27 files, 4,011 lines.

## 4. Arrange campaign fixture data

These files are not loaded as Prolog programs. `read_file_to_terms/3` reads their
single ground `clues/1` term. They are counted here because they are tracked
`.pl` files. The common owner is the P-R0/cliff corpus, reachability is the
historical fixture checker and campaign reports, product model is reproducible
input data, drift risk is generator/hash drift, shared dependencies are
`gen_mesh_fixture.py` and the campaign seed manifest, and replacement coverage
is the recorded generation parameters and hashes. Their data disposition is a
separate Phase 3 decision: retain only for independent reanalysis, otherwise
remove and rely on reports plus Git.

| Path | Exclusive lines | Primary ownership | Disposition |
|---|---:|---|---|
| `benchmarks/probe_arrange/fixtures/cliff_09x09_18w_seed11.pl` | 24 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_09x09_18w_seed12.pl` | 24 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_09x09_18w_seed13.pl` | 24 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_15x15_44w_seed11.pl` | 50 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_15x15_44w_seed12.pl` | 50 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_15x15_44w_seed13.pl` | 50 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_84w_seed11.pl` | 90 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_84w_seed12.pl` | 90 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_84w_seed13.pl` | 90 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_88w_seed11.pl` | 94 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_88w_seed12.pl` | 94 | P-R0 corpus | HOLD DATA DECISION |
| `benchmarks/probe_arrange/fixtures/cliff_21x21_88w_seed13.pl` | 94 | P-R0 corpus | HOLD DATA DECISION |

Fixture-data subtotal: 12 files, 774 lines.

## 5. Fill campaign source

| Path | Lines | Owner / primary ownership | Reachability | Product model | Evidence | Drift risk / shared dependencies | Replacement coverage | Disposition |
|---|---:|---|---|---|---|---|---|---|
| `benchmarks/probe_f1/dict_stats.pl` | 79 | Dictionary shape probe / P-F1 | D: F1/experiments only | Current calls for retired sizing question | F1 | Medium; loader/index/fixtures | Fill loader tests and load ratchet | HOLD evidence repair, then REMOVE |
| `benchmarks/probe_f1/load_stages.pl` | 89 | Loader attribution / P-F1 | D: F1/experiments only | Copied 2026-07-05 pre-current loader | F1 | High; private loader pipeline | Load ratchet and loader tests | HOLD evidence repair, then REMOVE |
| `benchmarks/probe_f1/oneoff_baseline_repro.pl` | 44 | `call_time/2` offset diagnosis / P-F1 | D: F1 conclusion only | One-off old-baseline reproduction | F1 | Medium; SWI-version-locked counts | Permanent fill ratchet | REMOVE |
| `benchmarks/probe_f1/profile_search.pl` | 48 | Search profiler cross-check / P-F1 | D: F1 only | Current entry over retired MRV interpretation | F1 | High; profiler tree and old search model | Fill tests/identity/ratchet | HOLD evidence repair, then REMOVE |
| `benchmarks/probe_f1/run_search.pl` | 114 | Attribution runner / P-F1 | D: F1; S: search instrument | Current MAC vs copied pre-MAC search | F1 | Very high; non-equivalent search trees | Fill identity/ratchet | HOLD reconstruction recipe, then REMOVE |
| `benchmarks/probe_f1/search_instr.pl` | 198 | Copied search instrument / P-F1 | S: F1 runner; D: F1 | Retired full-recount MRV engine | F1 | Very high; copied engine and unscoped globals | Product tests and F1 findings | HOLD with runner, then REMOVE |
| `benchmarks/probe_fh2/build_probe.pl` | 34 | Artifact build/size probe / F-H2 | D: FH2B only | Stale schema-v2 current API caller | FH2B | High; hard-coded missing `/tmp/claude-1000/f-h2` | Artifact tests/identity | HOLD recipe, then REMOVE |
| `benchmarks/probe_fh2/index_stats.pl` | 58 | Index sizing / F-H2 | D: FH2G only | Historical proposed mask index | FH2G | Medium; dictionary and SWI dependence | Current mask equality tests | HOLD recipe, then REMOVE |
| `benchmarks/probe_fh2/kernel_bench.pl` | 112 | Native/WASM kernel gate / F-H2 | D: FH2G only | Historical ordset-vs-bignum microbench | FH2G | Medium; backend/GMP/host | Mask tests and artifact identity | HOLD decision/reference repair, then REMOVE |
| `benchmarks/probe_fh2/phase_a.pl` | 262 | Post-F-H1 attribution / F-H2 | D: FH2B/experiments only | Copied pre-MAC path vs current search | FH2B | Very high; private hooks and missing stale artifacts | Current MAC tests/identity/ratchet | HOLD recipe, then REMOVE |
| `benchmarks/probe_fh2/phase_b.pl` | 73 | Controlled kernel ablation / F-H2 | D: FH2B only | Labels no longer match current MAC path | FH2B | Very high; refused schema-v2 artifacts and `/tmp` | Artifact mask tests/identity | HOLD recipe, then REMOVE |

`probe_f1` subtotal: 6 files, 572 lines. `probe_fh2` subtotal: 5 files,
539 lines. The F1 report must use the available source commit `22e8c21` with
historical base `a178923`; the currently named `0329b2f` object is absent from
this clone. The F-H2 reports must state their phase-specific patch recipes before
removal.

## 6. Fill-quality source

| Path | Lines | Owner / primary ownership | Reachability | Product model | Evidence | Drift risk / shared dependencies | Replacement coverage | Disposition |
|---|---:|---|---|---|---|---|---|---|
| `benchmarks/fill_quality/gen_grids.py` | 91 | Grid/mask generator / permanent quality comparison | S: `run.sh`, `matrix.sh`; D: FQ | Current input generator with duplicated `amer11` | FQ | High; embedded mask and stale engine prose | Stock-grid tests do not prove generator parity | RETAIN; use shipped grid as source |
| `benchmarks/fill_quality/score_fill.py` | 165 | Independent scorer / permanent quality comparison | S: `run.sh`; D: FQ/plans | Current independent post-hoc oracle | FQ | Medium; canonical schema and STW folding assumptions | Engine report is intentionally not independent | RETAIN; add focused tests later |
| `benchmarks/fill_quality/run.sh` | 43 | External quality comparison / permanent quality comparison | D: FQ/plans | Current CLI on easy grids; missing fills are reportable | FQ | High; external STW/ingrid and shell normalization | No complete AC-FILL-12 gate | RETAIN; harden/pair with gate |
| `benchmarks/fill_quality/matrix.sh` | 64 | FS-3(b)/FS-4 frontier / fill-quality campaign | D: FQ/design decision | Current CLI but non-gating historical frontier | FQ | High; external tools and collapsed failures | No AC-FILL-12 replacement yet | HOLD AC-FILL-12, then REMOVE |
| `benchmarks/fill_quality/probe_mac.pl` | 705 | MAC/dom-wdeg adoption prototype / fill-quality campaign | D: FQ and implementation plan | Historical 705-line engine twin | FQ | Very high; private hooks/global state/copied engine | Current MAC tests/identity/ratchet, but quality gate gap | HOLD AC-FILL-12 and evidence, then REMOVE |
| `benchmarks/fill_quality/probe_mac_b0.pl` | 521 | B0 instrumentation twin / fill-quality campaign | S: B0 launcher; D: B0/FQ | Historical exact-replay product copy | B0 | Very high; many private MAC hooks/globals | B0 counters are durable; gate gap remains | HOLD AC-FILL-12 and provenance, then REMOVE |
| `benchmarks/fill_quality/run_b0_instrument.sh` | 30 | B0 launcher / fill-quality campaign | D: B0; S: B0 twin | Historical campaign launcher | B0 | High; machine-local STW path | B0 report records commands/results | HOLD AC-FILL-12, then REMOVE |

Fill-quality subtotal: 7 files, 1,619 lines. The retained permanent portion is
3 files and 299 lines. The 4 historical frontier/adoption files are 1,320 lines
and cannot be removed yet.

## 7. Campaign-test assertion dispositions

### 7.1 `tests/probe_arrange.plt`

Every assertion in this file imports campaign modules. `REPLACE` means the named
permanent assertion must exist and be reviewed before this assertion disappears;
`RETIRE` means the assertion protects only a campaign mechanism or removed
schema and its dated evidence is the durable disposition.

| Assertion | Lines | Disposition | Named retained or replacement coverage |
|---|---:|---|---|
| `fixture_exact_count` | 9-11 | RETIRE campaign loader wrapper | `core:load_clues_from_prolog_fixture` |
| `fixture_wrong_count_throws` | 13-15 | RETIRE campaign exact-count wrapper | C4 permanent explicit-selection validation owns workload selection |
| `authority_tiny_budget_is_censored_not_proven` | 17-20 | REPLACE | `arrange:strict_budget_exhausted_not_proven`; `strict_budget_not_proven_is_not_infeasible` |
| `authority_matches_direct_product_tiny_and_500m` | 22-30 | RETIRE redundant self-differential | Same budget tests plus strict identity/ratchet |
| `twin_replays_both_nontranspose_corners` | 32-41 | RETIRE campaign replay | `core:full_tree_solution_count_matches_assoc_reference`; strict identity |
| `twin_replays_seeded_layout` | 43-51 | RETIRE | `arrange:seed_reproducible`; `seed_places_all_words` |
| `full_observer_preserves_exact_decision_order` | 53-61 | RETIRE mechanism counter contract | P0 report; no product observer invariant |
| `ad1_assoc_and_direct_count_selection_decision_traces_match` | 63-80 | REPLACE broad historical differential | Focused `direct_buckets` tests and permanent assoc full-tree oracle |
| `ad2_full_and_delta_count_selection_decision_traces_match` | 82-99 | REPLACE broad historical differential | Focused core A-D2 tests and permanent assoc full-tree oracle |
| `ad2_seeded_full_and_delta_traces_match_both_corners` | 101-112 | RETIRE after focused coverage review | Core A-D2 tests, seeded product tests, strict heavy identity |
| `canonical_state_key_is_absolute_and_orientation_sensitive` | 114-126 | RETIRE rejected P-C0 cache | C0 report; no product cache exists |
| `fingerprint_bucket_never_establishes_equality` | 128-139 | RETIRE rejected P-C0 cache | C0 report; no product cache exists |
| `semantic_cutoff_never_marks_dead` | 141-147 | RETIRE campaign cutoff/cache | C0 report; no corresponding product mechanism |
| `twin_operation_preserves_shared_seed_stream` | 149-156 | REPLACE | `arrange:arrange_budget_shared_across_corners` and seeded product tests |
| `seed_cleanup_on_failure_exception_and_interruption` | 158-165 | RETIRE campaign wrapper | `arrange:seed_cleared_restores_deterministic`; retained owners need local cleanup tests |
| `decision_limit_is_censored` | 167-171 | RETIRE semantic campaign cap | P0 schema/report only |
| `trace_row_keeps_configured_cutoff_separate_from_termination` | 173-184 | RETIRE campaign schema | P0 report; no retained schema |
| `trace_row_interruption_preserves_configured_cutoff` | 186-197 | RETIRE campaign schema | P0 report; product failure contracts remain elsewhere |

The compact independent assoc/full-tree oracle
`tests/core.plt:full_tree_solution_count_matches_assoc_reference` is explicitly
retained. The focused A-D2 coverage for bucket classification, proof
multiplicity, residue restoration, and stale non-sharing state is also retained.
The one older core assertion that still calls `refresh_direct_counts/8` must be
rewritten against current delta behavior or shown redundant with the stronger
stale-residue assertion before that probe-only product recount is removed.

### 7.2 Campaign-dependent `tests/greedy_benchmark.plt`

| Assertion | Lines | Disposition | Named retained or replacement coverage |
|---|---:|---|---|
| `record_readback_and_retention` | 10-14 | RETAIN | Permanent greedy baseline persistence |
| `replay_equivalence_easy` | 16-19 | RETIRE with A-G1 replay | GI and greedy inference ratchet |
| `replay_equivalence_heavy` | 23-25 | RETIRE with A-G1 replay | GI and greedy inference ratchet |
| `semantic_counter_partitions` | 27-44 | RETIRE hard-coded A-G1 counters | Sweep inference ratchet and GI |
| `candidate_phases_match_product` | 46-63 | RETAIN while permanent phase twins exist | Current raw-pool/postprocess differential |
| `best_effort_phases_match_product` | 65-85 | RETAIN while permanent phase twins exist | Current raw-pool/postprocess differential |
| `direct_attempt_slots_and_strict_omission` | 87-109 | RETIRE with four-direct replay | GI raw-pool/selected identity plus focused G2 tests |
| `g2_transpose_is_involution_with_fresh_clue_vars` | 111-119 | REPLACE probe helper | Move/extend product transpose geometry/fresh-variable test in `tests/arrange.plt` |
| `g2_product_transpose_geometry_and_fresh_clue_vars` | 121-132 | RETAIN and move | Product transpose geometry, involution, fresh clue variables |
| `g2_product_preserves_four_corner_block_order` | 134-138 | RETAIN and move | Product visible block order |
| `g2_product_omits_symmetric_setup_failures` | 140-147 | RETAIN and move with local fixture | Product setup-failure symmetry |
| `g2_product_preserves_dropped_terms_order_and_fresh_copies` | 149-161 | RETAIN and move | Product dropped-term order/copy freshness |
| `g2_additional_samples_cover_required_behaviors` | 163-179 | RETIRE broad campaign differential | The three focused product tests above |

The permanent A-G2 set must therefore cover transpose geometry and involution,
fresh clue variables, block order, setup-failure symmetry, and dropped-term
ordering/fresh copying before `g2_transpose.pl` is removed.

## 8. Frozen sentinels

### 8.1 Environment and tests

| Sentinel | Phase-0 value |
|---|---|
| Branch | `cleanup/benchmark-probe-tooling-review` |
| Commit | `1119899358a196a86d792d057df33050c4dd61d1` |
| SWI-Prolog | `10.1.10` (`x86_64-linux`) |
| Native suite | 453 passed, 0 failed, 0 timed out, 0 blocked (414 named tests) |
| Goldens | All passed |
| CLI/stderr checks | All passed |

### 8.2 Strict identity digests

`benchmarks/check_arrange_identity.sh --heavy` passed all 15 rows:

| Row | SHA-256 |
|---|---|
| `ladder_09x09_08w.pl` | `a5a6a03997794053ecbad4ea994f6132e3f84e754f9f456e1a5e5299663343d1` |
| `ladder_15x15_12w.pl` | `b6e8d412a6016c5f384aaf31fc066fc2c3e2a38db20e8ea5bcdc38986c2ecdd4` |
| `ladder_21x21_25w.pl` | `5c56768a610b7df00e5f436a0cf2feb205c19382032723476697af88097e6193` |
| `ladder_15x15_28w.pl` | `437410f510fdca91ee4a13c3574ba2cc0ddeb536a7509044e913c0a4b08726ae` |
| `ladder_15x15_32w.pl` | `005590d960ed304e14c3855344a51d8af4be0189e90bcf52c64ae26e79f23b60` |
| `real_13x13_12w.pl` | `e00925d58820f7db6886d76616ed523bb7e7115f1032f85083066ed94dabf4d4` |
| `real_15x15_18w.pl` | `faa15b1b95a0665756669941cc8ad6bd7cf11f8b60580e1665af670c5ef974ea` |
| `ladder_09x09_16w.pl` | `5cecff18a4b5becc0f8d52cec1b39b0b8b8c5b9a7d2466cd1849d6873cacb532` |
| `ladder_21x21_80w.pl` | `5c3f74a40fdc3aa6908e5489b59a0645d6d8b5462b2265f4a52642f37a211291` |
| `ladder_15x15_34w.pl` | `6485a6d8f5b1f16a946f9d2c642b371b494c4ebcb21d5f66e6205087d5433c75` |
| `ladder_15x15_36w.pl` | `9f0291ed76a97c934147d49d487d2de87a3d8f67a39b2dfeac942b30d25851f2` |
| `ladder_09x09_17w.pl` | `03812aeec0b3c1e237c6df391ab7dddd1eb2c1ab6e9d021d39dcc23feb517902` |
| `ladder_15x15_40w.pl` | `8db2391aeec891654f86cbc2c4d9022652cc256f18b8f17ccb5fee3401ee038a` |
| `ladder_21x21_82w.pl` | `a5f631627002f4fe77e525d220e568e30aaf4f9a9c3ed4ae8f2d219f438388a2` |
| `benchmark_16_dense_words.pl` | `90289af7db529b0132bae8bb910a18e90daa6f7b1c19de7cb97508883a56c81b` |

The last value is the strict latency-probe identity digest cited by the
campaign reports.

### 8.3 Greedy and fill identity digests

Greedy complete semantic document:
`ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.

Both raw and artifact fill identity checks passed all rows against the same
domain manifest:

| Row | SHA-256 |
|---|---|
| `g09_full` | `abb5479cd08e24e6996e161fc383107d8a9398e205ed5cf243a3dcb8394f7f51` |
| `g11_full` | `3212674691ad9f1ce5346f7505bd38947c63c5bb96dfafe0c5e2d834c3e832e9` |
| `g11_full_seed` | `62788e2eb2bda9c27bc514bdb65fca4f5c5838a24bcb91f6fd132c9eb8aea040` |
| `g13_full` | `18a81f18b177644d97f332079e938943fbad91a1dd0641becba10e0b7eee3e40` |
| `g15_full` | `fb04fcfe5391f8a7b1a84efc24f09538341053036c4e14571d645b71d27b5bfe` |
| `g17_50k` | `0049488b4ac75511232659cf25b308caaf29b9ac1466bee98f26d4573cb4bedb` |
| `g17_full` | `ce22ca04ecd0dd8df14530d61f88d1ccb2c734501c1bed0854b885ffbb88d5db` |
| `g21_full` | `f0fe6113b9743ab002264c8a6f80923790621c97e617ff1dc11c839a0651fd71` |
| `sq04_50k` | `b938a4a51bd358b19ef80b0640641c5956b333ed966493d6c22cf7ccec5a4da5` |
| `sq04_full` | `d9446f195a3f114b20308d00384b6a7b815904ff82f0def576ac1f3d46b4e24a` |
| `sq05_full` | `fac4c9314fe4d6ca4aa5f36dbcccd13510bcfade8feb01a31a55545b3784701a` |

### 8.4 Strict core plus heavy ratchet

`make bench-check BENCH_ARGS=--heavy` passed. Every gated row was an exact
same-SWI match; the latency probe is informational.

| Row | Baseline and measured search inferences | Status |
|---|---:|---|
| `ladder_09x09_08w.pl` | 22,297 | exact |
| `ladder_15x15_12w.pl` | 74,355 | exact |
| `ladder_21x21_25w.pl` | 216,653 | exact |
| `ladder_15x15_28w.pl` | 240,483 | exact |
| `ladder_15x15_32w.pl` | 627,744 | exact |
| `real_13x13_12w.pl` | 2,091,663 | exact |
| `real_15x15_18w.pl` | 172,537 | exact |
| `ladder_09x09_16w.pl` | 431,007 | exact |
| `ladder_21x21_80w.pl` | 3,107,968 | exact |
| `ladder_15x15_34w.pl` | 7,928,254 | exact |
| `ladder_15x15_36w.pl` | 14,074,891 | exact |
| `ladder_09x09_17w.pl` | 19,700,703 | exact |
| `ladder_15x15_40w.pl` | 4,480,185 | exact |
| `ladder_21x21_82w.pl` | 3,770,401 | exact |
| `benchmark_16_dense_words.pl` | 500,008,116 | exact, latency-only |

### 8.5 Greedy core plus heavy ratchet

`make bench-greedy-check BENCH_ARGS=--heavy` passed. Each tuple is
`construction / sweep / postprocess` inferences, baseline and measured exactly.

| Row | Exact inference tuple |
|---|---|
| `bundled_17_candidates` | 12,811 / 127,545 / 6,910 |
| `bundled_11_best_effort` | 2,098 / 10,037 / 29 |
| `benchmark_08_candidates` | 43,676 / 164,166 / 5,974 |
| `real_13x13_12w_best_effort` | 36,087 / 358,826 / 90 |
| `real_15x15_18w_best_effort` | 98,467 / 627,096 / 130 |
| `ladder_15x15_32w_best_effort` | 368,867 / 1,241,744 / 233 |
| `ladder_21x21_80w_best_effort` | 239,344 / 4,910,962 / 593 |

### 8.6 Fill core plus heavy ratchet

`make bench-fill-check BENCH_ARGS=--heavy` passed. Each tuple is
`search_inf / load_inf`, baseline and measured exactly.

| Row | Exact gated inference tuple |
|---|---|
| `sq04_full` | 3,302,946 / 12,395,209 |
| `g11_full` | 3,458,153 / 12,395,209 |
| `g11_full_seed` | 3,414,404 / 12,395,209 |
| `sq05_full` | 3,379,784 / 12,395,209 |
| `g17_full` | 3,766,202 / 12,395,209 |
| `g21_full` | 7,537,238 / 12,395,209 |
| `g13_full` | 3,608,308 / 12,395,209 |
| `sq04_50k` | 1,153,139 / 3,614,719 |
| `g15_full` | 3,701,888 / 12,395,209 |
| `g17_50k` | 9,756,833 / 3,614,719 |
| `g09_full` | 3,691,031 / 12,395,209 |

No `--record`, `--promote`, or `--log` mode was used. At the end of Phase 0,
Git reported only this inventory and the plan status edit; no baseline, history,
identity manifest, golden, or dated result report was modified.

## 9. Phase-0 exit gate

- Every tracked benchmark source or Prolog fixture-data file has an owner,
  primary ownership, reachability, model, evidence, drift/dependency assessment,
  replacement coverage, disposition, and exclusive line count.
- Every assertion in `tests/probe_arrange.plt` and every campaign-dependent
  assertion in `tests/greedy_benchmark.plt` has an explicit retained, replaced,
  or intentionally retired disposition.
- Native tests, strict identity, greedy identity, fill raw/artifact identity,
  and all three core plus heavy ratchets are green.
- The assoc/full-tree oracle and the required A-D2/A-G2 focused invariants are
  explicitly protected.
- The AC-FILL-12 gap blocks fill frontier deletion.
- Phase 1 correctness work may begin; Phase 2 deletion may not begin until all
  Phase 1 gates are green.
