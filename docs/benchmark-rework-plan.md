# Benchmark rework — plan & spec

**Status:** DONE (2026-07-06). Phases 1–3 built (Phase 4 was built *as amended*
during the arrange campaign — see "Phase 4 — as built"); the Phase-2 fixture
realism correction and the §11 decisions were closed out last, on branch
`bench/infra-finish` (see [`plans/bench-infra-finish.md`](./plans/bench-infra-finish.md)).
Adversarially stress-tested against the codebase; findings folded in (search for
"stress-test" — the amendments are inline). Verdict: mechanism sound, Phase 1
safe to start; the corrections are concentrated in Phase-2 workload realism (§6),
the memory metric (§7), and a dependency the first draft missed (§3, M4).
**Charter:** give `arrange` a principled performance-measurement setup that
answers *two distinct questions without letting them tangle*:

1. **Product:** how fast is the `arrange` *command* a user runs (end-to-end
   latency + peak memory), and how much of that is the search itself?
2. **Research:** which *search strategy* is best, and did a change to the core
   solver regress it? (existing concern — preserved, not extended)

**Metric of record:** `inferences` (deterministic, machine-independent) for all
*search* comparisons; `wall` + `maxRSS` for *command* latency, valid only within
a fixed `(machine, commit)` pair.

---

## 1. Why rework (the accretion this prevents)

`benchmarks/run_benchmarks.pl` began as "measure the search's efficiency," then
accreted a `--strategy` knob, then a second harness (`run_matrix.pl`) to sweep
strategies × fixtures. The two harnesses now re-implement the same timing loop
with *subtly different definitions* — different CSV columns, and different
medians (`run_matrix` takes the upper-of-two central values; `run_benchmarks`
averages them; noted in `docs/prolog-audit-findings.md` F025).

Root cause: **two questions (product efficiency vs. algorithm research) shared
one tool with no structural seam**, and the *measurement loop was never
factored out*, so a new subject meant a new harness. The rework fixes the
mechanism, not just the symptom.

## 2. Principles

1. **Factor the measurement core; subjects are plug-ins.** One engine does
   warmup → iterate → summarize over an opaque *sampler*. Adding something to
   benchmark = registering a sampler, never editing the loop, never a new
   harness.
2. **Product-bench and research-bench are separated *by construction*.** The
   product bench's subjects are "the `arrange` command" and "the `arrange`
   search (mrv_inc)". It has **no concept of strategy**, so it *cannot* accrete
   one.
3. **Reproducibility = git + a versioned result log, not live museum code.**
   Record each batch with its commit SHA / machine / config; reproduce by
   checkout + rerun. A strategy stays in live code only if it's an *active
   comparison anchor* (`baseline`) or *load-bearing for a test* (it is — see
   §8), not merely "so we could re-run it."
4. **One result schema, one metric definition, one charter per bench.** Every
   result is the same record; every bench file states its question / metric of
   record / valid-comparison rule in a header.

## 3. Target architecture

```
benchmarks/
  bench_core.pl        NEW  measurement engine: measure/3 + samplers + summarize
  subjects.pl          NEW  goal/sampler adapters (arrange_command, arrange_search,
                            strategy_search)
  workloads.pl         NEW  PRODUCT-bench manifest: rows of
                            (fixture, size, mode, iters, warmup, expected)
  run_arrange.pl       NEW  PRODUCT bench — command latency + search, no research knobs
  run_matrix.pl        KEEP (retargeted onto bench_core) — RESEARCH strategy matrix
  fixtures.pl          KEEP as-is — RESEARCH manifest for run_matrix AND
                            start_sensitivity.pl (do NOT delete/absorb it — see below)
  start_sensitivity.pl KEEP as-is — second consumer of fixtures.pl (bench_fixture/5)
  run_benchmarks.pl    RETIRE — folded into the two above
```

**Do not "absorb" `fixtures.pl`.** The stress-test found a *second* consumer,
`benchmarks/start_sensitivity.pl` (`consult(fixtures.pl)` + iterates
`bench_fixture/5`), beyond `run_matrix`. The clean split is: `fixtures.pl` stays
the RESEARCH manifest (both research consumers keep reading it, unchanged), and
`workloads.pl` is a *new, product-only* manifest. No migration of the research
consumers is needed, and nothing breaks.

| File | Question it answers | Metric of record | Subjects |
|---|---|---|---|
| `run_arrange.pl` | "Is the shipping command fast?" | inferences (search); wall+RSS (command) | `arrange_command`, `arrange_search` |
| `run_matrix.pl` | "Which strategy is best / did the core regress?" | inferences | `strategy_search` (all 4 strategies) |

## 4. `bench_core.pl` — the shared engine

The engine measures an opaque **sampler**: a closure that, invoked once,
produces one sample as a dict of `metric → number`. It knows nothing about
crosswords, strategies, or the CLI.

```prolog
% measure(+Sampler, +Opts, -Summary)
%   Sampler : called as call(Sampler, Sample); Sample is a dict of metric->number,
%             e.g. _{wall:0.03, cpu:0.03, inferences:319183}
%   Opts    : _{warmup:W, iterations:N}
%   Summary : _{iterations:N, warmup:W, stats:_{Metric:_{min,median,mean}, ...}}
%
% Runs W warmup samples (discarded), then N measured samples; for EACH metric
% key present in the samples, computes {min, median, mean}. Heterogeneous metric
% sets are fine (command samples carry wall+rss; search samples wall+cpu+inferences).
measure(Sampler, Opts, Summary).
```

**Failure/timeout is NOT `measure/3`'s job** (stress-test m1 — the Phase-1
must-fix). A sampler can fail (an infeasible cell where `once(Goal)` fails) or a
research cell can time out. `measure/3` measures a sampler that is *already known
to succeed*; it does not decide solved/no/timeout. `run_matrix` must keep its
existing `solve_status/6` gate (`call_with_time_limit(60, …) → solved|no|timeout`,
`run_matrix.pl:96-101`) OUTSIDE `measure/3`, calling `measure/3` only for
`solved` cells. Contract: if the sampler ever fails/throws inside `measure/3`,
that is a bug in the caller's gating, and `measure/3` should error loudly rather
than record a bogus zero-work sample. Nail this before writing the loop.

**Two built-in samplers** (subjects wrap these):

```prolog
% In-process Prolog goal — the search layer. Uses call_time/2 (returns
% time{cpu,inferences,wall}); inferences is the deterministic metric of record.
inproc_sampler(+Goal, -Sample) :-
    call_time(once(Goal), T),
    Sample = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

% External process — the command layer. Runs the CLI under
%   /usr/bin/time -f "%e %M"   (NOT -v: stress-test m4)
% which prints wall-seconds and max-RSS-KiB as two plain numbers — no m:ss / h:mm:ss
% parsing. Peak memory is only available at process level; see §7 for what it does
% and does NOT tell you (it is dominated by swipl's ~10.7 MB base — it is a process
% footprint, not the search's working set). GNU-time-only; fine for this Linux bench.
process_sampler(+Argv, -Sample) :-
    Sample = _{wall:Seconds, rss:KiBMax}.
```

**One median definition** (average the two central values for even N — the
`run_benchmarks.pl` semantics). This unifies the two harnesses' disagreement;
immaterial for `inferences` (deterministic ⇒ min == median == mean), a
sub-noise change for wall on the matrix.

## 5. `subjects.pl` — goal adapters

```prolog
% PRODUCT — end-to-end command (Layer A). --out /dev/null so serialization +
% write syscalls are measured but no disk is touched.
arrange_command_sampler(Fixture, Size, Mode, Sample) :-
    size_flag(Mode, Size, Flag),        % --size N | --max-size N
    process_sampler(['arrange','--input',Fixture,Flag,'--out','/dev/null'], Sample).

% PRODUCT — the search itself (Layer B). The REAL 4-corner path
% (construct_corners + rescore + select + number), module-qualified — verified to
% reproduce arrange_solve's search exactly (identical placement + deterministic
% inferences), minus emit and minus check_unique_answers/1 (negligible cost; only
% a duplicate-answer fixture would diverge — such fixtures are excluded). Size
% framing (fixed/max) is irrelevant here: both frame emit only; the search always
% runs on the N×N canvas (verified). Outcome is BOUND and returned so the
% workload's Expected can be asserted (m2) — arrange_best_layout ALWAYS succeeds,
% even on infeasible inputs (Outcome=infeasible, reward=-1), so discarding it
% would silently measure the wrong thing.
arrange_search_sampler(Words, GridLen, Sample) :-
    inproc_sampler(crosswordsmith_arrange:arrange_best_layout(Words, GridLen, _, _, Outcome), S0),
    Sample = S0.put(outcome, Outcome).

% RESEARCH — single-corner core search for one strategy (the existing matrix cell).
strategy_search_sampler(Strategy, Words, Grid, Start, Sample) :-
    inproc_sampler(find_crossword(Strategy, Grid, Words, Start, _, _), Sample).
```

**Warmup ≥ 1 is load-bearing, not cosmetic** (stress-test m5). The first
in-process call is a measured outlier (verified: iter 1 = 681,016 inferences;
iters 2+ = 673,335, identical thereafter). The `min == median == mean` invariant
the plan leans on for `inferences` holds *only* with `warmup ≥ 1`. Any product
workload row therefore needs `warmup ≥ 1`; the research `warmup=0` rows
(`benchmark_16_dense`, `benchmark_70_mesh`) only get away with it because they
run `iters=1`.

**Search-vs-rest breakdown** (the decomposition you asked for): per (fixture,
size), report `command.wall` (A), `search.wall` (B), and `rest ≈ A − B`
(startup + clue-load + emit + JSON-write + I/O). Flagged in output as
*indicative, not exact* — A is a process measurement and B is in-process, so the
subtraction is a decomposition estimate, not an identity. **`rest` can go
negative** for search-dominated cells where A ≈ B (stress-test m3 measured
`benchmark_70_mesh`: A ≈ 1.08 s, B ≈ 1.13 s ⇒ rest ≈ −0.05 s — inter-run noise).
For the trial, clamp `rest` at 0 and annotate `A≈B (search-bound)` when
`|A−B|` is within noise. A cleaner future option (Phase 4-ish) is to get the
search share from *inside* the command process — a bench-only path where the CLI
emits its own `call_time` of the search to stderr — rather than subtracting two
independent measurements; deferred, since it means a bench affordance on the
product CLI.

## 6. `workloads.pl` — a NEW product-only manifest

A new, separate manifest (NOT an absorption of `fixtures.pl` — see §3). One
predicate:

```prolog
% PRODUCT: (Fixture, Size, Mode, Iterations, Warmup, Expected)
%   Mode ∈ {size, max_size};  Warmup ≥ 1 (load-bearing, §5)
%   Expected ∈ {placed, infeasible} — ASSERTED against the sampler's bound Outcome
%     (m2); a mislabeled row is a hard error, not a silently-wrong measurement.
arrange_workload('fixtures/benchmark_14_words.pl', 17, size, 20, 3, placed).
...
```

`bench_fixture/5` stays in `fixtures.pl`, unchanged — it is the *research*
manifest and has its own consumers (§3).

*[As built the manifest is `arrange_workload/9`:
`(Fixture, Size, Mode, Iterations, Warmup, Expected, Tier, Gate, Budget)` —
`Tier` (core|heavy) selects the default vs `--heavy` set, `Gate` (inf|latency)
encodes the budget-saturation rule below, and `Budget` lets the search layer
run hard rungs to true completion above the shipped 500 M. See the
`workloads.pl` header for the authoritative column semantics.]*

**The existing `benchmark_*` fixtures are mostly UNSUITABLE as product
workloads** (stress-test M1/M2 — the biggest correction to this plan):

- `benchmark_20_words` places only at `--size 37`; `benchmark_26_words` only at
  `--size 49` (infeasible at every smaller size tried). They are size-brittle
  single-corner skeletons that interlock only on their sparse hand-tuned grid —
  not grids any user would request for 20/26 words. **Drop them from product
  scope.**
- `benchmark_16_dense` **saturates the 500M inference budget** (`arrange.pl:190`):
  three corners exhaust it, ~32 s/run, and its inference count is pinned to the
  budget *constant* (≈ 500,000,197) — degenerate as an algorithm signal. Include
  it only as **latency-only, `iters=1, warmup=0`, and NOT regression-gated on
  inferences**. (At the illustrative `iters=20` it would be ~12 min for one cell.)
- `benchmark_14_words @17` is a genuine, fast, place-all product case (verified:
  places 14, reward 588) — keep this class.

So the product manifest is not "tweak the old sizes"; it needs **new, realistic
blocked-13/15 fixtures** authored in Phase 2. Rule to encode: **any cell whose
measured `inferences ≈ arrange_budget` is budget-saturated → `iters=1`,
latency-only, excluded from the inference regression gate.** Today the *only*
clean search-dominated (non-saturated) stressor is `benchmark_70_mesh @21`
(~1.1 s, ~19.7 M inferences).

## 7. Result schema, metrics & comparison

One record per (subject, workload):

```prolog
_{ subject:  arrange_command,                         % | arrange_search | strategy_search
   workload: _{fixture:'benchmark_14_words.pl', size:17, mode:size},
   metrics:  _{wall:_{min,median,mean}, rss:_{min,median,mean}},   % keys depend on subject
   iterations: 20, warmup: 3,
   env: _{machine:'…', swiProlog:'10.1.10', commit:'86328f5', epochSeconds:…} }
```

| Metric | Source | Portable? | Use |
|---|---|---|---|
| `inferences` | `call_time` (search) | **yes** — deterministic (INV-2) | search regression gate; cross-machine |
| `wall`, `cpu` | `call_time` / `time -f` | no — machine-dependent | latency reporting; within `(machine,commit)` |
| `rss` | `/usr/bin/time -f` max RSS (command) | no | command footprint; NOT a search-memory metric |

**`rss` does NOT answer "how much memory does the search use"** (stress-test M3).
It is dominated by swipl's constant base (~10.7 MB); the search-attributable
delta measured across fixtures was only 0.9–3.0 MB (8–28% of the total). So `rss`
is "peak footprint of a swipl process running the command" — useful as a coarse
command-footprint / CI-trend number, useless for isolating the search's working
set. SWI exposes **no per-goal peak**; the only in-process handle is
`statistics(globalused/localused/trailused)` sampled at solution time, which is a
snapshot after GC, not a true peak. Bottom line for the original question:
*inferences* is the honest proxy for search cost; there is no clean search-memory
number, and the plan does not pretend otherwise.

**Comparison rule (stated in each bench header):** compare `inferences` freely
across runs/machines; compare `wall`/`cpu`/`rss` only within one machine+commit,
or as a CI trend. This is what stops "keep the old variant around to compare" —
the comparison is the recorded number, not resurrected code.

## 8. What does NOT change

- **The 4 strategies stay in `core.pl`.** Two are genuinely load-bearing:
  `baseline` is a *correctness oracle* in `tests/crossword.plt` (enumeration-count
  agreement, lines ~208–218) and `mrv_inc` is the production search
  (`arrange.pl:249`). `mrv`/`mrv_capped` rest on the weaker "comparison anchor"
  rationale (`strategies([baseline,mrv,mrv_capped,mrv_inc])`, `core.pl:91`) — kept
  for now, but they are the candidates to retire-to-the-log if the strategy
  question is ever declared closed. Deleting strategies is off the table for
  *this* rework regardless.
- **`experiments.md`** remains the research log `run_matrix` feeds.
- **The test suite, goldens, determinism fuzz** — untouched; benchmarks are not
  on the test path. The only observable change is `run_matrix`'s CSV: identical
  `inferences`, a sub-noise `wall`-median shift from the unified median def.

## 9. Phased migration (each phase independently verifiable)

**Phase 1 — extract `bench_core`; retarget `run_matrix`.** (Unblocked — touches
neither `fixtures.pl` nor product workloads.)
- Write `bench_core.pl` (`measure/3`, both samplers, `summarize`).
- **Must-fix (m1):** define `measure/3`'s success-only contract and keep
  `run_matrix`'s `solve_status/6` timeout/`no` gate OUTSIDE it (§4). This is the
  one thing the retarget must not drop.
- Replace `run_matrix.pl`'s inline measured-loop with `bench_core`; keep its CLI,
  CSV, and `solve_status` gating.
- *Verify:* capture `run_matrix` CSV before; after, `inferences` columns
  byte-identical; `wall` medians differ only by the intended median-def change;
  `solved|no|timeout` gating preserved. `./run_tests.sh` green (unaffected).

**Phase 2 — product bench.** *[DONE]*
- Add `subjects.pl`, `workloads.pl` (new product manifest), `run_arrange.pl`
  (Layers A+B + breakdown; text/csv/json out). Wire `Expected` vs the sampler's
  bound `Outcome` (m2). Warmup ≥ 1 on all product rows (m5). *[built during the
  arrange campaign, with the mesh cost ladder as the manifest]*
- **Author realistic fixtures** — do NOT rely on the `benchmark_*` set for
  product (M1/M2): keep `benchmark_14 @17`, mark `benchmark_16_dense` latency-only
  `iters=1`, drop `benchmark_20/26`, and add a couple of real blocked-13/15
  puzzles. `benchmark_70_mesh @21` is the one usable search-dominated stressor.
  *[CLOSED 2026-07-06, as amended by the as-built ladder: the realism gap is
  filled by `fixtures/real_13x13_12w.pl` / `real_15x15_18w.pl` — real ENABLE
  dictionary words planted on a legal witness (`benchmarks/gen_real_fixture.py`)
  at the blocked daily sizes, as core rungs. `benchmark_16_dense` runs as a
  heavy `gate=latency` rung at the SHIPPED 500 M budget (verified: `placed`,
  count pinned at 500,000,190 ≈ the budget constant). `benchmark_20/26` stay
  research-only (`fixtures.pl`). `benchmark_14`/`benchmark_70_mesh` were
  superseded by the ladder, which covers their cost points with planted-witness
  instances; they too remain research-matrix fixtures. Real-word instances near
  the density cliff measured wildly instance-noisy (15×15: 12w=1.7M, 14w=85K,
  16w=853M, 18w=248K warm inf) — anchors are deliberately well-behaved
  instances, per the same rule as the ladder's envelope guards.]*
- *Verify:* `arrange_search` inferences reproduce (min==median==mean, warmup≥1);
  `command.wall > search.wall` on non-search-bound cells; `rest` plausible
  (clamped/annotated when A≈B, m3); budget-saturated cells flagged latency-only.
  *[holds for the final manifest: both real rungs min==median, reproduced across
  independent runs; `gate=latency` implemented in workloads/run_arrange/
  check_baseline — a latency rung can neither fail nor win the inference gate]*

**Phase 3 — retire `run_benchmarks.pl`.** *[DONE]*
- Delete the file; repoint `make bench` → `run_arrange.pl`; keep
  `make bench-matrix` → `run_matrix.pl`. *[done during the arrange campaign]*
- **Rewrite** the README benchmark section (currently ~L473–507, `BENCH_GRID` /
  `--grid` / strategy examples) — a ~35-line rewrite to the product-bench flags,
  not a one-line repoint. Leave audit/history docs as point-in-time records.
  *[done during the campaign; residual drift (heavy tier described as "~26 s
  budget-saturating probes", a dead `--fixture bundled` example, no ratchet
  paragraph) fixed 2026-07-06]*

**Phase 4 — regression baseline (DEFERRED until Phase 2 numbers are seen).**
- Commit `benchmarks/baseline.json` (per subject/workload: `inferences` exact +
  `wall`/`rss` machine-tagged) and a `make bench-check` that diffs with
  tolerances (inferences: tight; wall/rss: loose, same-machine). Not before we
  know what's stable.

**Phase 4 — as built (2026-07): hill-climbing ratchet + cost ladder + history.**
The exact regression *gate* became a *hill-climbing instrument* (the goal shifted
to tracking arrange perf toward a WASM runtime, where the inference count is the
portable signal — identical native vs WASM, only wall-per-inference changes).
- **Workloads → a cost ladder** (`benchmarks/workloads.pl`): planted-witness mesh
  fixtures (`gen_mesh_fixture.py`) whose warm search cost climbs ~0.1M→182M across
  three grid sizes — 9×9, 15×15 (deep density ladder), 21×21 — since size changes
  branching geometry / 4-corner start / N² memory, so a win at one size can regress
  at another. Difficulty is driven by word density + small alphabet K, NOT size.
  Each size has a difficulty *cliff* (fast-resolve or thrash); each "hard" rung is
  the densest that still completes deterministically. `Budget` raised to 2e9 so hard
  rungs run to true completion (a ratchetable count) instead of saturating.
- **`baseline.json` is now a ratchet** (`mode:"ratchet"`, `metric:"search_inf"`,
  `regression_tolerance_pct`): a DROP is a win (`make bench-record` to accept), a
  RISE past tolerance fails; SWI-version mismatch downgrades a regression to WARN.
  `core` rungs run by default; `--heavy` adds the hard tail. `check_baseline.pl`
  hosts check / `--record` / `--log` / `--history`.
- **`benchmarks/history.jsonl`** — append-only, git-tracked ledger (one JSON line
  per run, stamped with git commit + timestamp). The baseline is a *moving*
  reference so it can't show a trajectory; the ledger keeps per-rung `search_inf`
  comparable over time. `make bench-record`/`bench-log` append; `make bench-history`
  renders the per-rung trend (latest, step Δ, cumulative Δ).

## 10. Makefile / invocation (post-rework)

```make
bench:         ## product: command latency + search, over the arrange workloads
	swipl -q benchmarks/run_arrange.pl -- \
	    $(if $(BENCH_FIXTURE),--fixture $(BENCH_FIXTURE),) \
	    $(if $(BENCH_SIZE),--size $(BENCH_SIZE),) \
	    --mode $(BENCH_MODE) --iterations $(BENCH_ITERATIONS) \
	    --warmup $(BENCH_WARMUP) --format $(BENCH_FORMAT)

bench-matrix:  ## research: strategy × fixture inference matrix (CSV)
	swipl -q benchmarks/run_matrix.pl -- $(BENCH_STRATEGIES)
```

## 11. Open decisions — ALL RESOLVED (closed 2026-07-06)

1. **Command I/O:** RESOLVED — `--out /dev/null`, as recommended (the as-built
   `subjects.pl` uses it). *Verified safe:* `with_output/2` is
   capture-then-plain-write (`core.pl:118-123`), no temp-file/rename, so writing
   to `/dev/null` works and leaves the device intact (the suspected blocker was a
   non-issue).
2. **Framings to bench in Layer A:** RESOLVED — `--size` only as the rung mode.
   Measured 2026-07-06 on `ladder_15x15_12w` (7 runs each, same host, under
   concurrent load — indicative): `--size` 91–94 ms, `--max-size` 91–96 ms;
   the crop/emit delta is sub-noise, and the search count is identical by
   construction (both modes frame emit only — `subjects.pl size_flag/4`).
   `mode: max_size` stays supported in the manifest schema for ad-hoc use.
3. **Canonical product fixture/size set:** RESOLVED — the mesh cost ladder
   (9×9/15×15/21×21, hill-climb instrument) + the real-word realism anchors
   `real_13x13_12w` / `real_15x15_18w` (core) + the `benchmark_16_dense @17`
   budget-saturation latency probe (heavy, `gate=latency`). `benchmark_20/26`
   are OUT of product scope (research matrix only).
4. **`arrange_best_layout/5` access:** RESOLVED — module-qualified call, as
   recommended, onto the budget-explicit **/6** (the bench must raise the budget
   above the shipped 500 M so hard rungs complete). This is a sanctioned
   white-box reach: purity-audit **C25** annotation at `benchmarks/subjects.pl`,
   with the export deferred in STATUS.md's de-accretion roadmap (trigger: the
   next change to `arrange.pl`'s export surface). Every recorded `search_inf`
   is defined against `/6`, so the access point must not silently move.
   *(2026-07-06: the trigger fired — the thin-form work touched the export
   surface, so `/6` is now exported and `subjects.pl` plain-imports it. The
   access point is unchanged, ratchet-verified +0.00% on every gated rung;
   the white-box annotation is retired. C25 fully closed.)*

## 12. Non-goals

Config-driven bench DSL; pluggable reporters beyond text/csv/json; a committed
regression baseline before Phase 2; deleting strategies; benchmarking `fill` /
`lint` / `export` (this rework is `arrange`-scoped).
