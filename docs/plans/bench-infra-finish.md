# Benchmark-infra finish-up

**Branch:** `bench/infra-finish` · **Source of truth:** [`docs/benchmark-rework-plan.md`](../benchmark-rework-plan.md)
**Goal:** close out the rework plan's unfinished phases — the product-fixture
realism correction (Phase 2 §6/§9) and the cleanup/decision residue (Phase 3,
§11) — against the *as-built* state left by the arrange campaign (ratchet +
mesh cost ladder, plan "Phase 4 — as built").

## 1. Audit: plan vs tree (what is already done)

Verified on `main` (78c086a) before this branch:

| Plan item | State |
|---|---|
| Phase 1 (bench_core, run_matrix retarget) | **done** — `benchmarks/bench_core.pl` exists; `run_matrix.pl` consumes it |
| Phase 2 harness (`subjects.pl`, `workloads.pl`, `run_arrange.pl`) | **done** — plus the Phase-4-as-built ratchet (`check_baseline.pl`, `baseline.json` mode:ratchet, `history.jsonl`) |
| Phase 2 **fixtures** ("author realistic blocked-13/15; drop 20/26; 16_dense latency-only") | **NOT done** — `workloads.pl` is 100% planted-mesh `ladder_*` rungs (`gen_mesh_fixture.py`, alphabet ≤ 8, lengths 3–4). No real-word rung exists; `benchmark_16_dense` is absent rather than latency-only |
| Phase 3 retire `run_benchmarks.pl` | **done** — file gone; `make bench` → `run_arrange.pl`; `make bench-matrix` kept |
| Phase 3 README rewrite | **mostly done** — product/matrix split documented; BUT the heavy-tier prose ("budget-saturating latency probes (~26 s each)") describes rungs that no longer exist (the as-built heavy tier is completion-run ladder rungs, 0.6–8.7 s, ratcheted). Same drift in `run_arrange.pl --heavy` help text |
| §11.1 `--out /dev/null` | **resolved** (as-built uses it; plan already records the verification) |
| §11.2 `--size` vs `--max-size` in Layer A | **open** — all rungs use `mode=size`; no recorded decision |
| §11.3 canonical product fixture/size set | **open** — the ladder exists but no decision text; realism gap above |
| §11.4 `arrange_best_layout/5` vs `/6` access | **effectively resolved** — as-built uses module-qualified **/6** with a sanctioned white-box annotation (purity-audit C25, commit 085a44c); STATUS.md de-accretion roadmap defers the export until `arrange.pl`'s export surface next changes. **No code change here** — this plan only records the resolution |
| Plan doc header | says "Status: proposed … No code written yet" — stale |

`benchmark_20/26` were never added to the product manifest (they live only in
`benchmarks/fixtures.pl`, the research manifest, which stays as-is per plan §3)
— so "drop from product scope" is a documentation act, not a deletion.

## 2. Deliverables

### D1 — real-word product fixtures (the Phase-2 correction)

The stress-test's point (plan §6 M1/M2): the product bench must include inputs
a real user would feed `arrange` — real English words at the blocked-standard
daily sizes 13/15 — not only tiny-alphabet synthetic meshes. `arrange` itself
produces a sparse interlocking mesh (not a blocked grid), so "real blocked
13/15" means: **real-dictionary word sets with a blocked-puzzle length profile
(4–9 letters + a long seed anchor), at `--size 13` and `--size 15`, that
arrange places in full, deterministically, without budget saturation.**

Mechanism: `benchmarks/gen_real_fixture.py` — same planted-witness algorithm as
`gen_mesh_fixture.py` (seed word at (0,0) across; attach perpendicular crossing
words only by legal solver moves ⇒ satisfiable + reachable by construction;
shuffled output order), except candidate words are drawn from the frozen
public-domain ENABLE list (`fixtures/dict/enable1.txt`, see
`fixtures/dict/README.md`) via a `(length, position, letter)` index instead of
random small-alphabet strings. Deterministic for a given (dict bytes, G, N,
Lmin, Lmax, seed).

Rungs: one 13×13 and one 15×15 fixture (`fixtures/real_13x13_NNw.pl`,
`fixtures/real_15x15_NNw.pl`), word count chosen **empirically** — the densest
generated instance whose warm `arrange_best_layout/6` completes with a
deterministic count at core-tier cost (target ≲ 0.3 s / low-millions of
inferences), probed across several seeds/densities. Both join `workloads.pl`
as `core` rungs (`Expected = placed`) and `baseline.json`.

### D2 — `benchmark_16_dense` as a latency-only rung; 20/26 dropped

Plan §6 rule to encode: *a cell whose measured inferences ≈ budget is
budget-saturated → iters=1, latency-only, excluded from the inference gate.*

- Extend `arrange_workload` /8 → **/9** with a `Gate` column: `inf` (default,
  ratchet-gated) | `latency` (search_inf reported info-only, never gates).
  Touches: `workloads.pl` (all rows), `run_arrange.pl` (pass `gate` through
  into the result row, like `tier`), `check_baseline.pl` (a `gate:latency` row
  can neither fail nor win the inference ratchet; `--record` copies the flag
  into new baseline specs). Baseline entries without `gate` default to `inf`.
- Add `fixtures/benchmark_16_dense_words.pl` @17 as a `heavy` +
  `gate=latency` rung at the **shipped 500 M budget** (the point is the
  user-felt worst-case cliff the README already describes: four corners share
  one budget; a deep corner burns it even after a layout was found). iters=1,
  warmup=0. Expected outcome verified empirically before committing (§6 says
  it places; the count is pinned to the budget constant — degenerate as a
  signal, which is exactly why it is latency-only).
- `benchmark_20/26`: record (here + rework-plan §6) that they are **out of
  product scope** — size-brittle skeletons (place only @37/@49); they remain
  research-matrix fixtures in `benchmarks/fixtures.pl` untouched.

### D3 — §11 decisions recorded; doc drift fixed

- **§11.2:** probe `--size` vs `--max-size` command wall on one fixture (the
  modes differ only in emit framing; the search count is identical by
  construction — `subjects.pl size_flag/4`). Expected sub-noise ⇒ decision:
  `--size` stays the canonical rung mode; `max_size` remains supported for
  ad-hoc runs. Probe numbers recorded in the rework plan (indicative only —
  this machine is under concurrent load).
- **§11.3:** canonical product set = mesh cost ladder (hill-climb instrument)
  + `real_13x13`/`real_15x15` (realism anchors) + `benchmark_16_dense`
  (latency-cliff probe); `benchmark_20/26` dropped.
- **§11.4:** record the as-built resolution (module-qualified `/6`, C25
  sanction, export deferred to the next export-surface change per STATUS.md).
- Rework plan header updated (no longer "proposed / no code written"); Phase
  2/3 marked done-as-built with pointers here.
- README heavy-tier prose + `run_arrange.pl --heavy` help fixed to describe
  the actual heavy tier (completion-run ladder rungs + the one opt-in
  budget-saturating latency probe).

## 3. Verification

1. `swipl -q benchmarks/check_baseline.pl` (= `make bench-check`) **before**
   any change: all core rungs `ok` (search_inf is load-independent; wall/rss
   informational).
2. After: same check green — existing rungs' `search_inf` byte-identical
   (adding rungs must not move them); new rungs appear (added to baseline via
   a `--record` run or hand-authored spec matching `new_rung_spec/2`'s schema;
   diff reviewed so ONLY new entries + info-only fields change).
3. `make bench BENCH_ARGS="--fixture real"` exercises the new rungs; a single
   `--heavy` pass (or `--fixture 16_dense`) exercises the latency rung once
   (~1 min — budget-saturating by design).
4. Determinism: each new rung's `search_inf` min == median across iterations
   (warm), and a repeated run reproduces the count exactly.
5. `make unit` is **not** part of this surface (benchmark files are not on the
   test path; no `prolog/` or `tests/` file changes) — run only if any shared
   .pl file is touched.

## 4. Non-goals

Changing `arrange.pl` (incl. its export list — C25 stays deferred); touching
`benchmarks/fixtures.pl` / `run_matrix.pl` / the fill bench; re-recording any
existing rung's baseline number; wasm/ and the fragment work (owned by
sibling branches).
