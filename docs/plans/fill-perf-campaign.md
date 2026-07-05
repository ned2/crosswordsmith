# Plan: the fill performance campaign

Status: DRAFT (pending critical review). Author: arrange-campaign orchestrator,
2026-07-05. Companion docs: `docs/research/arrange-perf-campaign-2026-07.md`
(the playbook this reuses), `docs/experiments.md` (ledger; fill entries will
append there), `docs/research/swi-vm-wasm-performance.md` (WASM cost model).

## Goal

Make `fill` (Flavour-B dictionary fill of a blocked stock grid, `fill.pl`)
fast enough for interactive use at realistic scale — 15x15/21x21 grids against
10k-100k-word dictionaries — both as a native CLI and in the browser via WASM,
using the same measure-first, byte-identity-gated, ratcheted methodology that
took arrange down 71-84%.

## Why there is headroom (assessment, 2026-07-05)

Fill shares only the periphery with arrange (pw/8 records, numbering, emit,
fragment seeds, stockgrid masks). Its search core is separate and pre-dates
the campaign's lessons, except one it already embodies: the shared-cell-
variable grid (crossings by unification, undo by trail) is the E-H2-native
design. The identified gaps:

1. **No incremental counting.** `select_mrv/6` recomputes candidate counts
   for EVERY unfilled slot at EVERY node. Binding a word only constrains
   slots crossing the filled slot; all other counts are provably unchanged —
   and unlike arrange's over-estimate carry-forward, carried counts here
   would be EXACT, so the tree and output stay identical by construction.
2. **Ordset/list machinery at dictionary scale.** Candidate counting does
   `ord_intersection` over integer ordsets that scale with dictionary size;
   the literature-standard alternative (this IS the well-studied
   dictionary-fill problem) is bitset intersection: one unbounded integer
   per (len,pos,char) index set, `/\` to intersect, `popcount` for the MRV
   count. Candidate materialization walks `nth0` over list buckets
   (O(bucket) per candidate); buckets should be compound terms read by
   `arg/3`.
3. **Dictionary load is an every-invocation startup tax** (read, normalize,
   sort, full positional index over all lengths — including lengths the
   grid has no slot for). WASM roughly triples it.

## Constraints (inherited from the arrange campaign charter)

- Deterministic output is a product requirement. All perf work must be
  order-preserving: same candidate enumeration order, same MRV tie-breaks,
  byte-identical output. Output-changing ideas escalate to the orchestrator.
- The ratchet is strict: relative tolerance (0.5%), zero regressions, no
  input-size-dependent switches or magic constants.
- Inference counts are the tracked metric (machine-independent, WASM-proxy).

## Phase 0 — benchmarking build-out (blocking; nothing exists today)

Fill's only fixtures are a 3x3 toy grid and a 26-word list; the bench
machinery (`workloads.pl`, `run_arrange.pl`, `check_baseline.pl`) is
arrange-only. Deliverables:

- **Stock-grid fixtures**: legal blocked masks (pass `stockgrid_validate`)
  at 9x9, 15x15, 21x21 with realistic block density/symmetry; a small
  deterministic generator or hand-authored masks committed under fixtures/.
- **Benchmark dictionaries**: a seeded deterministic generator (analogous to
  `gen_mesh_fixture.py`) producing wordlists at ~5k and ~50k words with
  roughly natural letter/length distributions. Rationale: license-clean,
  reproducible, tunable difficulty. (Open question for review: is synthetic
  letter distribution realistic enough for fill difficulty, which is
  dictionary-structure-dependent? Fallback: vendor a public-domain list,
  e.g. ENABLE, as one rung.)
- **`benchmarks/run_fill.pl`** mirroring run_arrange: command-layer wall/rss
  AND in-process `fill_attempt` inferences under `call_time`, with a THIRD
  measured layer fill has that arrange lacks: `load_dict` cost, reported
  separately from search (so a search win is never masked by load noise and
  vice versa).
- **Workload manifest + ratchet**: a fill ladder (grid size x dictionary
  size x seed density) spanning ~0.1M to ~100M search inferences, chosen by
  measured cost like the arrange ladder; `fill_baseline.json`; extend
  `check_baseline.pl` (parameterize the runner+baseline paths) rather than
  fork it. Same warmup discipline (arrange found a fixed ~7.6k cold/warm
  JIT wobble; verify fill's equivalent).
- **A scale golden**: at least one 15x15 fill golden locking byte-identity
  (the 3x3 golden exercises nothing).
- Acceptance: counts reproduce at +0.00% across repeated runs; budget-
  saturating rungs excluded from the ratchet (the P2 lesson: near-cliff
  marginal instances flip on any change).

## Phase 1 — instrumentation probe P-F1 (before any optimization)

The arrange campaign's highest-ROI step. Measure on the new ladder:
- Inference attribution: dictionary load vs search; within search, the
  select_mrv recount share vs candidate materialization vs unification.
- Backtrack profile: fill may backtrack far more than arrange (dictionary-
  scale branching factors); if so, tree-shape techniques are NOT dead here
  the way P1 killed them for arrange — value ordering is off-limits
  (determinism), but the cost balance decides F-H1 vs F-H2 priority.
- Fan-out stats: slots per grid, candidates per slot over time, index-set
  sizes at each dictionary scale.

## Phase 2 — the experiments (serial, ratcheted, byte-identity-gated)

- **F-H1 — incremental candidate counts** (the inc_count port): carry
  per-slot counts in threaded state; recount only slots sharing a cell with
  the just-filled slot. Soundness note recorded up front: `Used`-word
  exclusion happens at try time, NOT count time, in the current code — the
  port must preserve that exactly (counts that ignore Used), or the tree
  changes. Expected: large win on grids with many slots.
- **F-H2 — bitset dictionary index**: ordsets -> unbounded-integer bitmasks;
  count = popcount(mask), candidates = ascending bit iteration (same order
  as the ascending ordset -> order-preserving). Verify popcount/1 is an
  arithmetic function in pinned SWI 10.1.10 and measure GMP bigint AND/
  popcount cost at 50k-word scale before building the full change.
- **F-H3 — buckets as compound terms** (arg/3 candidate lookup); may fold
  into F-H2.
- **F-H4 — Used as a set** (order-preserving; same accept/reject).
- Sequence F-H1 vs F-H2 by P-F1's attribution; measure separately even if
  both look certain (the E-H9 attribution discipline).

## Phase 3 — startup + WASM posture

- Slot-length-filtered index build (only lengths the mask has); lazy
  per-length building if attribution supports it.
- Dictionary-load strategy for the browser (precompiled artifact aligned
  with the .qlf plan in the WASM research doc). Measured by the load layer
  added in Phase 0.

## Explicit non-goals

- Envelope/completability work (restarts etc.) — latency campaign only.
- Any change to fill's output or its infeasible/not_proven semantics.
- Real-dictionary curation beyond what benchmarking needs.

## Risks / open questions (for the critical review)

1. Synthetic-dictionary realism (letter distribution drives crossing
   compatibility, which drives both difficulty and index-set shapes).
2. Fill's difficulty distribution is unmapped — the ladder may need
   seed-pinning density as a difficulty knob, not just size; heavy tails
   may make deterministic rungs hard to pick (arrange P2's lesson).
3. Is per-node recounting actually dominant at dictionary scale, or does
   ord_intersection so dominate that F-H2 must go first? (P-F1 answers.)
4. Does fill_search's `\+ memberchk(Word, Used)` interact with any planned
   change in a tree-visible way?
5. How much of end-to-end fill latency is load_dict at realistic scale? If
   >50%, Phase 3 outranks Phase 2 for perceived performance.
6. bigint bitset cost model under WASM (GMP is optional there; LibBF
   fallback has worse bignum performance per the WASM research doc) — a
   native-only win that regresses WASM would be a bad trade; needs a
   WASM-side sanity check before F-H2 ships.

## Execution model

Same as the arrange campaign: orchestrator adjudicates; one Opus agent per
phase/experiment in an isolated worktree (verify base commit first — three
stale-worktree incidents last campaign); serial composition against
re-recorded baselines; every accept/reject gets an experiments.md entry;
research findings land as durable docs.
