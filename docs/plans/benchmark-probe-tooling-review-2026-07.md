# Candidate plan: benchmark and probe tooling review

Status: IN PROGRESS (2026-07-18). Phases 0 and 1 are complete; Phase 2 is in
progress. This is a review and cleanup plan, not an authorization to delete
every item listed below. Each candidate must pass its own reachability, evidence,
and replacement-coverage gate before removal.

## 1. Goal

Reduce the maintenance cost of `benchmarks/` without weakening the permanent
performance ratchets, semantic identity gates, standing quality checks, or
research tools that still answer live questions.

The work has three ordered concerns:

1. Remove historical apparatus that models retired product paths, is locked to
   a completed campaign setup, or exists only to reproduce evidence already
   preserved in a dated report and Git history.
2. Refactor duplicated mechanics that remain after deletion into small shared
   utilities, while keeping arrange, greedy, and fill policy separate.
3. Correct safety and clarity problems in the retained tooling using the local,
   version-matched SWI-Prolog manual as the authority.

The order is load-bearing: do not build abstractions for campaign code that is
about to be deleted.

## 2. Scope and non-goals

In scope:

- all source under `benchmarks/`, including `probe_arrange/`, `probe_f1/`,
  `probe_fh2/`, and `fill_quality/`;
- benchmark and probe tests under `tests/`;
- benchmark Make targets and living README/spec references;
- product code retained solely as a historical probe reference;
- process capture, CLI selection, baseline persistence, history, identity
  manifests, fixture loading, and report-envelope mechanics.

Out of scope:

- changing solver behavior or benchmark workload composition;
- re-recording baselines, histories, or identity manifests to make cleanup pass;
- deleting dated Markdown result reports;
- merging strict arrange, greedy, fill, strategy-matrix, and fill-quality row
  semantics into one generic schema;
- adding product exports solely to make a test or historical probe convenient;
- optimizing hot solver paths as part of a tooling cleanup.

## 3. Verified starting point

### 3.1 Permanent regression infrastructure

Retain these roles even if their implementation is refactored:

| Role | Current files |
|---|---|
| Shared measurements | `bench_core.pl` |
| Strict arrange product bench | `run_arrange.pl`, `subjects.pl`, `workloads.pl` |
| Strict ratchet and identity | `check_baseline.pl`, `run_arrange_identity.pl`, `check_arrange_identity.sh`, `baseline.json`, `history.jsonl`, `arrange_identity.sha256` |
| Greedy product bench | `run_arrange_greedy.pl`, `greedy_subjects.pl`, `greedy_workloads.pl` |
| Greedy ratchet and identity | `check_greedy_baseline.pl`, `check_greedy_identity.sh`, `greedy_baseline.json`, `greedy_history.jsonl`, `greedy_identity.sha256` |
| Fill product bench | `run_fill.pl`, `fill_subjects.pl`, `fill_workloads.pl` |
| Fill ratchet and identity | `check_fill_baseline.pl`, `check_fill_identity*.sh`, `fill_baseline.json`, `fill_history.jsonl`, `fill_identity.sha256` |
| Strategy research | `run_matrix.pl`, `fixtures.pl`, `start_sensitivity.pl` |
| Reproducible fixture generation | `gen_fill_dict.py`, `gen_mesh_fixture.py`, `gen_real_fixture.py` |
| Layout analysis | `analyze_layout.py` |
| Independent fill-quality comparison | `fill_quality/run.sh`, `score_fill.py`, `gen_grids.py` |

The strategy matrix and start-sensitivity tool are live research consumers, not
vestiges. `docs/benchmark-rework-plan.md` and
`docs/plans/legacy-surface-dissolution.md` explicitly preserve them.

### 3.2 Historical-source review pool

The source review pool is up to 6,737 lines before accounting for code that may
need to remain as a small permanent oracle:

| Area | Source lines | Initial assessment |
|---|---:|---|
| `probe_arrange/*.{pl,py}` | 4,011 | Completed arrange campaign machinery, with some current test dependencies |
| `probe_f1/*.pl` | 572 | Pre-current-fill attribution and copied-search probes |
| `probe_fh2/*.pl` | 539 | F-H2 campaign probes, including machine-local paths |
| Fill MAC/B0/frontier campaign tools | 1,320 | Completed adoption and rejection evidence |
| `probe_backtrack.pl` | 183 | Explicitly non-runnable museum file |
| A-G2 top-level launchers | 112 | Completed premise/product comparison runners |

This is an audit ceiling, not a promised line reduction. The review must report
actual deleted, retained, and refactored lines at closeout.

### 3.3 Evidence policy

The existing benchmark architecture already states the right policy in
`docs/benchmark-rework-plan.md` section 2: reproducibility comes from Git plus a
versioned result log, not live museum code.

Apply that policy as follows:

- Keep dated Markdown reports under `benchmarks/results/`.
- Before deleting a file named by a report, make the report or the living
  historical-reconstruction index identify the commit containing runnable source.
  Keep corrections out of immutable dated evidence.
- For a runner locked to a base commit that predates the runner itself, record
  the exact reconstruction recipe: checkout the measurement base, apply the
  runner source commit as a patch without advancing `HEAD`, then invoke the
  command. Naming only the source-containing commit is insufficient.
- Keep accepted baseline/history/identity artifacts.
- Keep raw result data only when it supports independent reanalysis not captured
  by the report. Raw data does not by itself justify retaining a live runner.
- Do not create a second in-tree source archive. Git history is the source
  archive.

## 4. Classification method

Create one inventory row per tracked source file with these fields:

| Field | Question |
|---|---|
| Path and owner | Which supported command, test, or report owns it? |
| Reachability | Is it called by Make, README instructions, tests, another source file, or only a dated report? |
| Product model | Does it exercise the current product path, copy a retired path, or call removed/private hooks? |
| Evidence | Is unique evidence available only by rerunning it, or already captured durably? |
| Drift risk | Which product internals or setup assumptions can change underneath it? |
| Replacement coverage | Which permanent test, identity gate, or ratchet protects the same invariant? |
| Disposition | Retain, refactor, remove, or hold pending a named extraction |

Disposition rules:

- REMOVE when there is no live command or invariant, the setup is obsolete, and
  dated evidence plus Git preserve the result.
- EXTRACT THEN REMOVE when tests still import a campaign twin or the twin is the
  only reason retired product code remains.
- RETAIN when it is a supported ratchet, identity gate, standing quality check,
  live research tool, or reproducible-input generator.
- REFACTOR only after the deletion pass, and only when a shared helper reduces
  net code and semantic drift across at least two retained consumers.

No-file-left-unclassified is the review-completeness gate.

## 5. Vestigial-code candidate register

### V1. Explicitly dead artifacts

Initial disposition: REMOVE.

Status: COMPLETE (2026-07-18). All four gates were confirmed. P1 provenance was
repaired in the append-only experiment ledger and the historical reconstruction
index; the G2 test oracle remains for Phase 3, independent of its two deleted
launchers. Deleted source: 339 lines.

| Candidate | Evidence | Gate before deletion |
|---|---|---|
| `benchmarks/probe_backtrack.pl` | Lines 2-13 declare it historical and intentionally non-runnable; it calls removed `probe_*` hooks at lines 59-64 | Confirm `docs/experiments.md` and the preserved patch identify the result and implementation |
| `probe_f1/oneoff_baseline_repro.pl` | One-time 44-line explanation of a `call_time/2` inference offset | Confirm the conclusion remains in `2026-07-05-p-f1-attribution.md` |
| `run_arrange_g2_probe.pl`, `run_arrange_g2_wall.pl` | No supported Make entry; only dated A-G2 reports call them | Separate any still-used test oracle from the top-level launchers |

### V2. Retired fill-path probes

Initial disposition: EXTRACT EVIDENCE REFERENCES, THEN REMOVE.

`probe_f1/` reconstructs the 2026-07-05 fill loader and copied pre-MAC search.
The current engine is MAC plus dom/wdeg, so these probes no longer attribute the
shipping path. Their findings are in
`benchmarks/results/2026-07-05-p-f1-attribution.md` and `docs/experiments.md`.

`probe_fh2/` reconstructs the pre-current count/index path and includes hard-coded
`/tmp/claude-1000/f-h2/` inputs in `build_probe.pl`, `phase_a.pl`, and
`phase_b.pl`. Its gate and attribution results are already in the two dated F-H2
reports. `kernel_bench.pl` should remain only if the review identifies a current
decision or standing gate that consumes it.

Deletion gate:

- Replace living-doc claims that the source is "merged to document the method"
  with result-report and historical-commit references.
- Verify no current test or Make target imports these directories.
- Do not rewrite these probes to the current engine merely to keep them alive.

### V3. Completed fill-quality campaign rigs

Initial disposition: SPLIT PERMANENT CHECK FROM HISTORICAL CAMPAIGNS.

Retain the independent scorer and comparison driver. They are not, by
themselves, the complete AC-FILL-12 standing gate: `run.sh` covers only the easy
quality grids and treats failed fills as reportable output, while the contract in
`docs/design-spec.md` requires the `blocked_13a` reference rows at score
thresholds 30 and 1 to complete under the default budget and retain quality.

Before deleting any historical frontier machinery, add one explicit permanent
quality-gate command that:

- runs both AC-FILL-12 reference rows with an STW-class dictionary;
- fails on a missing, budget-exhausted, or malformed fill;
- enforces recorded completion and quality floors;
- retains independent post-hoc scoring and sidecar cross-checking;
- clearly reports its external STW dependency and is run whenever that dependency
  is available.

Review for removal after that gap is closed:

- `fill_quality/probe_mac.pl`, a 705-line file that labels itself
  "NOT SHIPPABLE ENGINEERING" and predates the accepted product engine;
- `fill_quality/probe_mac_b0.pl` and `run_b0_instrument.sh`, the completed B0
  instrumentation twin and launcher;
- `fill_quality/matrix.sh`, the completed FS-3(b)/FS-4 frontier tool unless a
  current product decision still uses it.

Before deletion, split the 533-line `fill_quality/README.md` into a short living
operational guide and dated evidence references. Update plans that link directly
to `probe_mac.pl` to identify its historical commit instead.

Preserve scorer independence: do not replace `score_fill.py` with engine-side
score reporting. Make `grids/amer11.json` the single source of truth rather than
manually synchronizing an embedded copy in `gen_grids.py`.

### V4. Completed arrange campaign apparatus

Initial disposition: EXTRACT PERMANENT INVARIANTS, THEN REMOVE IN DEPENDENCY
ORDER.

The `probe_arrange/README.md` explicitly calls the directory campaign-only. The
campaign is closed, and its mechanism evidence is preserved in dated reports.
The directory nevertheless remains coupled to the normal test suite through
`tests/probe_arrange.plt` and `tests/greedy_benchmark.plt`.

Candidate groups:

| Group | Principal files | Why it is historical now |
|---|---|---|
| Phase-0 authority/twin framework | `probe_arrange.pl`, `run.pl`, `verify_controls.pl`, `measure_overhead.pl`, `measure_process.py`, `profile.pl` | It embeds a second strict DFS plus campaign counters; the accepted results and final authority envelope are recorded |
| P-D0 | `d0_support.pl`, `measure_d0.pl`, `test_d0_support.pl` | The support-delta premise led to accepted A-D1/A-D2 and is no longer the product algorithm |
| P-C0 | `pc0_duplicate_work.pl` and P-C0 tables/counters in `probe_arrange.pl` | Both cache avenues are measured closed |
| P-R0 | `run_pr0.py`, `analyze_pr0.py`, `schema.py`, `test_schema.py`, seed/corpus checks | `run_pr0.py` refuses to run outside historical branch `probe/a-r0` at commit `1bccf47`; Track R is closed |
| A-G2 | `g2_transpose.pl` and top-level launchers | It replays the retired four-direct-corner path; product now derives transpose partners |
| A-D1 | `ad1_buckets.pl`, `run_ad1_wall.pl` | It retains assoc and full direct-recount references superseded by A-D2 |
| A-D2 | `ad2_delta.pl`, `measure_ad2.pl`, `run_ad2_wall.pl` | It copies the now-shipping A-D2 recursion and imports D0 for an independent scanner |
| Closeout | `closeout_direct.pl`, `run_closeout_direct.pl`, `test_closeout_direct.pl`, `run_closeout_authority.py` | Final attribution is complete and durably reported |

These groups express conceptual ownership, not disjoint line counts.
`probe_arrange.pl` is the Phase-0 root and also contains P-C0 tables/counters;
`g2_transpose.pl` is included in the 4,011-line `probe_arrange/` total while its
two top-level launchers are the separate 112-line inventory row. The final
inventory must record primary ownership, shared dependencies, and exclusive line
counts so deletion totals are not double-counted.

Migration requirements:

- Move genuine product invariants into the owning `tests/core.plt` or
  `tests/arrange.plt`; tests must not keep a full old search engine alive.
- Produce an assertion-by-assertion disposition table for
  `tests/probe_arrange.plt` and the campaign-dependent portions of
  `tests/greedy_benchmark.plt`. Name the retained or replacement assertion before
  deleting each test.
- Preserve strict full-ladder identity and inference ratchets as the end-to-end
  equivalence guards.
- Preserve the compact independent full-tree assoc oracle in `tests/core.plt`;
  current-product identity alone is not a substitute for every independent
  semantic differential.
- Preserve focused A-D2 tests for bucket classification, proof multiplicity,
  residue restoration, and stale non-sharing state.
- Preserve focused A-G2 product tests for transpose geometry, fresh clue
  variables, block order, setup-failure symmetry, and dropped-term copying.
- Remove campaign-only tests rather than relabeling mechanism counters as product
  contracts.
- Remove the groups from leaves to roots: launchers and closeout, P-C0/P-R0,
  A-G2 after its product tests are probe-independent, A-D2, P-D0, A-D1, then the
  generic Phase-0 twin.

After the A-D1 differential disappears, remove
`core.pl:778-800` (`refresh_direct_counts/8` and `direct_recount/8`) if the
reachability audit confirms the probes are its only callers. Its source comment
already states that the old product path remains solely as the A-D1 reference.

Do not remove the 12 cliff fixtures or the P-R0 raw JSONL merely because their
runners are removed. Decide them as data: retain only if they support meaningful
independent reanalysis beyond the final reports; otherwise rely on Git and the
recorded generator hashes.

### V5. Completed A-G1 counters in the permanent greedy runner

Initial disposition: REFACTOR.

`greedy_subjects.pl:152-329` and `run_arrange_greedy.pl:87-99` run the accepted
A-G1 semantic-counter replay on every permanent greedy benchmark. The ratchet
gates construction, sweep, and postprocess inferences; the campaign counters are
informational. Retain raw-pool and selected-output identity, but remove the
A-G1-only counter pass and hard-coded mechanism counts after confirming they are
not needed to distinguish product output.

### V6. Campaign-only commands and stale instructions

Remove Make targets at `Makefile:229-243` when their probe dependencies are
gone. Replace `probe_arrange/README.md` with no active source-side archive; dated
reports and living benchmark docs should identify permanent commands.

Correct stale living text found during the audit:

- Root README says there are two benchmarks and omits permanent greedy and fill
  ratchets.
- Strict Make comments call a 9/13/15/21 workload set a "15x15 ladder".
- Fill comments still describe `load_inf` as informational although it is gated.
- Fill workload and fill-quality comments still describe the retired scoreless
  or fixed-budget MRV engine as current.
- `fixtures/README.md` uses obsolete benchmark Make variables.
- `analyze_layout.py` examples reference retired `crossword.pl` usage.

## 6. Retained-code correctness work

These fixes precede broad refactoring because they protect the cleanup process.

### C1. Declare closure arguments as meta-predicates

Priority: P0.

`bench_core:measure/3` calls its first argument as `call(Sampler, Sample)`, and
`inproc_sampler/2` passes its first argument to `call_time/2`, but neither has a
`meta_predicate` declaration. Calls currently depend on `user`-module script
loading or explicit qualification and can resolve in `bench_core` when called
from a normal module.

Add and test the appropriate closure declarations, including representative
module-local callers. Audit retained wrappers such as seed or limit wrappers for
the same issue. The authority is `docs/reference/swi-manual/metapred.md`, which
requires dynamically called arguments to retain caller qualification.

### C2. Make process capture exception-safe

Priority: P0.

The baseline checkers and several runners repeat `process_create/3`, pipe reads,
closes, and `process_wait/2` without one cleanup boundary. If a pipe read or close
throws before the wait, a requested PID can remain unreaped. JSON parsing in the
three current baseline checkers occurs after their normal close/wait sequence and
is not itself that leak window. Greedy CLI identity also drains stdout before
stderr, which can deadlock if the child fills stderr.

Build one retained process-capture utility using `setup_call_cleanup/3` that
always closes streams, waits for requested PIDs, returns status, and avoids
sequential two-pipe deadlocks. The local authorities are:

- `metacall.md:89-153`: cleanup runs on failure, exception, commit, and success;
- `packages/clib.md:84-147`: requesting `process(PID)` requires
  `process_wait/2`.

### C3. Make fill baseline recording match its contract

Priority: P0.

`check_fill_baseline.pl:43-49` promises read-back verification, but
`do_record/3` at lines 239-255 writes directly to the live baseline. Port the
temporary sibling, parse-back, completeness verification, and rename mechanics
already used by strict and greedy ratchets. Add fill promotion tests matching
`tests/benchmarks.plt` before sharing the implementation.

### C4. Reject false-green command selections

Priority: P0.

Arrange, fill, and greedy runners can accept an unmatched `--fixture` and emit a
successful empty result. Ratchets then describe rows as unmeasured rather than
failing. They also ignore positional arguments and validate output format only
after running workloads.

Validate positional arguments, format, count ranges, mutually exclusive checker
modes, and nonempty explicit selections before measurement. Follow
`optparse.md:179-199`: unknown flags are errors and positional arguments are a
separate result that callers must handle.

### C5. Make identity recording complete and atomic

Priority: P0.

`check_fill_identity.sh:51-84` records and exits successfully even if one or more
CLI runs set `status=1`; it can replace the manifest with a partial file. Bring
fill to the completeness behavior already present in
`check_arrange_identity.sh:75-82`. Review all identity scripts for unknown-option
handling and atomic replacement while preserving each domain's digest semantics.

### C6. Validate sample schemas

Priority: P1.

`bench_core:summarize_samples/2` chooses numeric keys from the first sample and
silently allows later samples to omit them. Require a stable numeric metric key
set for the current permanent samplers, or explicitly report per-metric sample
counts if a retained consumer truly needs heterogeneity. Do not add this check to
the timed hot path; validate during summarization.

### C7. Add an exact inference comparison mode

Priority: P0 for pure tooling refactors.

The normal ratchets intentionally allow tolerance, accept improvements, and
downgrade regressions when SWI versions differ. They therefore cannot prove the
hard no-change claim required for a refactor that should not affect measured
work.

Add an exact mode to the retained ratchet machinery that requires:

- the same SWI-Prolog version as the reference;
- every selected core and heavy row to be present;
- every gated inference metric to equal the reference exactly;
- any increase or decrease to fail rather than classify as a regression or win.

Use this mode for all pure tooling and shared-utility migrations. Continue using
the normal tolerance policy for actual solver-performance work.

## 7. Shared-utility candidates

Refactor only the permanent survivors of section 5.

### S1. Process lifecycle

One Prolog module should own checked subprocess capture, status propagation, and
stream/PID cleanup. It should support only the retained capture shapes that have
at least two callers. In particular, greedy identity must retain separate stdout
and stderr bytes without sequentially draining two bounded pipes.

Atomic temporary-file replacement is a separate storage concern under S2, not
part of the process module. Bash identity scripts may share a shell helper if
that is smaller than moving their lifecycle into Prolog.

### S2. Ratchet storage and history

The three baseline checkers duplicate JSON loading, row retention, new-rung
insertion, temporary writes, read-back, provenance, JSONL append, and history
rendering. Extract storage mechanics behind a small module. Each domain retains
its own:

- workload key and required metadata;
- gated and informational metric definitions;
- regression classification and SWI-version policy;
- report wording and domain-specific row conversion.

Do not create a generic benchmark configuration DSL. Prefer a few explicit
callback predicates or data arguments with PlDoc contracts.

### S3. Runner setup and selection

Arrange, fill, and greedy runners duplicate repository-root discovery, option
parsing, tier/filter selection, warmup/iteration overrides, SWI metadata, format
validation, and JSON envelope construction. Share these mechanics while retaining
domain-specific samplers and output columns.

Fixture loading can be shared between arrange runners, identity, matrix, and
start-sensitivity. Do not fold fill grid/dictionary setup into the arrange
fixture helper.

### S4. Identity manifest lifecycle

Share complete-record, stale/missing-ID checks, digest comparison, and atomic
replacement only if the resulting helper is smaller and clearer than the three
scripts. Keep separate payload producers and formats:

- strict uses per-rung diagnostics-bearing rows;
- fill uses CLI and artifact bytes;
- greedy hashes one complete semantic document.

### S5. Report envelopes

Share tool/version metadata and format validation. Keep strict search/rest,
fill load/grid/search, and greedy phase columns separate. A universal result row
would erase useful semantics and is explicitly rejected.

## 8. General Prolog cleanup register

Apply these only to retained code, unless a correctness fix is needed before a
safe deletion.

| ID | Cleanup | Manual basis |
|---|---|---|
| P1 | Put mutable seed, dynamic table, and `nb_*` state behind `setup_call_cleanup/3`; remove globals with `nb_delete/1` | `metacall.md`, `gvar.md` |
| P2 | Do not prepend growing traces with `nb_setarg/3`, which duplicates non-atomic values; use reverse accumulation or scoped `nb_setval/2` where appropriate | `manipterm.md:133-161`, `gvar.md:25` |
| P3 | Aggregate profiler results per predicate rather than taking the first call-tree node | `profile.md`, especially `profile_procedure_data/2` |
| P4 | Serialize JSONL directly; do not pass serialized JSON through `normalize_space/2`, which is not JSON-aware | `packages/json.md:200-237` |
| P5 | Report malformed history lines with line numbers instead of silently dropping ledger entries | `packages/json.md` |
| P6 | Give exported benchmark predicates proper `%!` PlDoc headers, closure modes, preconditions, and verified `det`/`semidet` tags | `preddesc.md:34-43` |
| P7 | Add explicit imports and remove unused imports so `autoload(false)` claims are true | `usingmodules.md`, `consulting.md` |
| P8 | Replace one-shot asserted root/dir predicates with explicit arguments or immutable setup where this reduces state | `dynamic.md` |
| P9 | Describe inference counts as SWI-version-locked regression signals, not absolute machine-independent quantities; remove hard-coded `call_time` corrections | `statistics.md`, `builtin-statistics.md` |
| P10 | Use ordered-set/assoc predicates only where their sorted/ground preconditions are established; do not churn hot paths for style | `ordsets.md`, `assoc.md` |

The review must not replace direct compounds, bitmasks, or measured hot-path
operations with generic maps/options merely for uniformity. Any cleanup touching
a gated sampler requires exact inference comparison.

## 9. Execution phases

### Phase 0: freeze sentinels and finish the inventory

- Record the branch/commit, SWI version, tracked-file inventory, source-line
  totals, and full caller/document reachability table.
- Record primary ownership and exclusive line counts for overlapping campaign
  groups.
- Record native test count and strict, greedy, and fill identity digests.
- Run all three ratchets, including heavy rows, without recording.
- Classify every tracked source file under `benchmarks/`.
- Build the assertion-by-assertion disposition table for tests that import
  campaign modules.

Exit gate: every source file has one owner and disposition; all sentinels are
green before cleanup.

### Phase 1: correctness before consolidation

- Implement C1-C6 with focused tests.
- Implement C7 and prove it rejects both a one-inference increase and decrease.
- Correct only documentation that materially misstates current safety or gate
  behavior.
- Do not introduce the broad shared ratchet module yet.

Exit gate: failure-path tests prove no empty-selection pass, partial identity
record, corrupt fill baseline, unresolved module-local closure, or unreaped
normal child process.

### Phase 2: delete zero-value and retired fill apparatus

- Remove V1 after report-reference checks.
- Retire `probe_f1/` and `probe_fh2/` after replacing living-source references
  with historical commit/report references.
- Split permanent fill-quality operations from completed campaign evidence and
  remove obsolete rigs only after the AC-FILL-12 gate gap in V3 is closed.

Exit gate: permanent fill bench, identity, ratchet, independent scorer, and
AC-FILL-12 quality gate remain discoverable; no living doc instructs users to run
removed code.

### Phase 3: dissolve arrange campaign dependencies

- Migrate permanent A-D2 and A-G2 invariants to owning product tests.
- Remove A-G1-only counters from routine greedy measurements while retaining
  semantic identity.
- Remove arrange probe groups in the dependency order in V4.
- Remove campaign-only Make targets and test-suite imports.
- Remove the old A-D1 product recount only after proving it has no remaining
  caller.

Exit gate: normal tests no longer import a campaign search twin; strict and
greedy identities and all gated inference counts are unchanged.

### Phase 4: refactor permanent shared mechanics

- Extract process lifecycle first.
- Extract ratchet storage/history second.
- Extract runner selection/metadata and arrange fixture loading third.
- Consolidate identity mechanics only if the result is a net simplification.

Each extraction gets equivalence tests before callers move. Delete old code in
the same commit that moves the last caller; do not leave compatibility wrappers
without a concrete consumer.

Exit gate: each shared module has an explicit narrow export list, `%!` PlDoc,
verified determinism, and at least two retained callers.

### Phase 5: retained-code cleanup and docs

- Apply the section 8 cleanup register.
- Reduce fill-quality README to operating instructions plus result links.
- Document strict, greedy, fill, strategy, and quality commands in the root
  README without campaign archaeology.
- Report before/after source lines, supported entry points, duplicate
  implementations, and test coverage.

## 10. Commit strategy

Keep reviewable concerns separate:

1. Safety fixes and focused regression tests.
2. Explicit dead-file deletion and reference repair.
3. Retired fill probe deletion and fill-quality split.
4. Arrange invariant migration.
5. Arrange campaign-probe and probe-only product-path deletion.
6. Process utility extraction.
7. Ratchet storage/history extraction.
8. Runner/identity cleanup where justified.
9. General cleanup, living docs, and closeout measurements.

Do not combine baseline changes, solver optimization, or output-contract changes
with these commits.

## 11. Verification

Run after every phase as applicable, and run the complete battery at closeout:

```sh
make unit
make test
make bench-arrange-verify BENCH_ARGS=--heavy
make bench-fill-verify BENCH_ARGS=--heavy
make bench-greedy-verify BENCH_ARGS=--heavy
make bench-matrix
make test-wasm
```

Additional focused gates:

- module-local tests for `measure/3` and `inproc_sampler/2` closure resolution;
- fill baseline temporary-write, read-back, new-rung, and heavy-row-retention
  tests;
- identity record failure tests proving manifests remain unchanged;
- invalid CLI mode, positional argument, format, and empty-filter tests;
- process read/close-exception tests proving closure and wait cleanup;
- a dual-pipe child that writes more than pipe capacity to stderr before closing
  stdout, proving retained capture cannot deadlock;
- exact same-SWI inference comparison over every core and heavy gated row;
- temporary-ledger tests for append, render, malformed-line diagnostics,
  whitespace inside JSON strings, and failure atomicity, plus all three history
  render commands;
- `autoload(false)` load and undefined-predicate checks over retained benchmark
  modules;
- Python unit tests for retained schemas/scorers;
- `swipl -q benchmarks/start_sensitivity.pl`, the retained research consumer
  with no Make target;
- the permanent AC-FILL-12 quality gate whenever its external STW-class
  dictionary is available; the ingrid head-to-head remains optional and is
  reported honestly when `ingrid_core` is unavailable.

Hard no-change gates:

- no baseline or history rewrite;
- strict, greedy, and fill identity digests unchanged unless a separately
  reviewed manifest-format migration is the explicit commit purpose;
- every gated inference count exactly unchanged under C7, not merely accepted by
  the tolerance ratchet, for pure tooling refactors;
- native and WASM product outputs unchanged;
- no dated result report deleted.

## 12. Completion criteria

The review is complete when:

- every tracked benchmark/probe source file has a documented disposition;
- no live source models a retired product path solely as museum code;
- no product predicate remains solely for a deleted differential probe;
- normal tests do not depend on full campaign-local solver twins;
- retained subprocess and baseline writes are exception-safe and tested;
- duplicated permanent mechanics have either one shared owner or an explicit
  documented reason to remain separate;
- living documentation names supported commands and historical reports name
  historical commits;
- closeout reports line-count change, entry-point change, test-count change, and
  all no-change sentinels.

The desired result is not the smallest possible `benchmarks/` directory. It is a
small, trustworthy permanent measurement system plus durable historical evidence,
with no active-looking campaign machinery left to drift.
