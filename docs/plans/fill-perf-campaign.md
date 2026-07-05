# Plan: the fill performance campaign

Status: EXECUTED 2026-07-05 — campaign complete. Four experiments accepted
(F-L1, F-H1, F-L2, F-H2 masks-opt-in), all output-byte-identical; g09_full
end-to-end 5.18s → 0.87s. Full record + close-out (residue, corrections):
`docs/experiments.md` "Fill campaign" section. Previous status: REVIEWED —
GO-WITH-CHANGES (red-team review 2026-07-05; all amendments folded in below,
each marked [RT]). Author: arrange-campaign orchestrator. Companion docs:
`docs/research/arrange-perf-campaign-2026-07.md` (the playbook this reuses),
`docs/experiments.md` (ledger; fill entries append there),
`docs/research/swi-vm-wasm-performance.md` (WASM cost model).

## Goal

Make `fill` (Flavour-B dictionary fill of a blocked stock grid, `fill.pl`)
fast enough for interactive use at realistic scale — 15x15/21x21 grids against
10k-100k-word dictionaries — native CLI and browser/WASM, using the
measure-first, byte-identity-gated, ratcheted methodology from the arrange
campaign.

## Strategic thesis [RT — reframed after empirical refutation]

Fill is NOT arrange. An empirical probe (red-team review) showed fill is
**search-tree-bound and heavy-tailed**: a trivial open 5x5 against a 500-word
dictionary saturates a 50M-inference budget; every cell is doubly constrained
and candidates come from a huge pool, so backtracking is deep and
content-dependent (the opposite of arrange's shallow-retreat, counting-
dominated profile). Consequences, stated up front rather than as risks:

- **Per-node wins buy latency on completing instances, never completion.**
  Halving node cost on a heavy-tail instance just hits the budget cap at 2x
  the nodes (arrange's P2/P3 decoupling, more severe here — and WASM makes
  backtracking disproportionately expensive on top).
- **The determinism charter freezes the tree.** Fill's candidate count is
  "words matching this slot's bound cells" — deliberately ignoring crossing
  consistency. True forward-checking/AC pruning would change counts, change
  MRV selection, and change output: off-limits. So no shippable lever in
  this campaign reduces the node count of a heavy-tail instance; the
  envelope lever for fill (like arrange) is a future randomized-restart
  product feature, out of scope here.
- Therefore this campaign's honest deliverable is: **fast completion of the
  completable regime** (and cheap proof of infeasibility), measured on a
  ladder of instances that complete deterministically.

## Why there is headroom

Fill shares only the periphery with arrange (pw/8 records, numbering, emit,
seeds, stockgrid masks); its search core pre-dates the campaign's lessons,
except the shared-cell-variable grid, which is already the unification-native
design. The gaps, ranked by measured evidence [RT]:

1. **Candidate materialization is O(bucket x candidates) per node** —
   `candidates/4` walks a linear `nth0` per candidate over list buckets
   (fill.pl:161-164). Measured: 30 candidates over a 5,000-word bucket =
   1.42 ms/op; the 200-candidate/50k case timed out at micro-bench reps.
   Likely the dominant search cost at dictionary scale.
2. **MRV counting does `ord_intersection` at dictionary scale** every node
   for every unfilled slot. Bigint AND + popcount measured 250-600x faster
   for the count at 5k-50k-word scale.
3. **No incremental counting** — `select_mrv/6` recounts EVERY unfilled slot
   at EVERY node; only slots sharing a newly-bound cell can change.
4. **Everything is recomputed for the winning slot too** — the full
   candidate list is materialized before `member/2` + used-word rejection
   discard most of it (the E-H9 "don't materialize what you'll discard"
   lesson).
5. **Dictionary load is an every-invocation startup tax** (full positional
   index over all lengths, even lengths the grid lacks). WASM ~triples it.

## Constraints (inherited)

- Deterministic, byte-identical output; order-preserving changes only;
  output-changing ideas escalate to the orchestrator.
- Strict ratchet (relative 0.5%, zero regressions, no input-size switches).
- Inference counts are the tracked metric; budget-saturating rungs are
  excluded from the ratchet (they flip on any change).

## Phase 0 — benchmarking build-out (blocking)

Fill has no bench: a 3x3 toy grid, a 26-word list. Deliverables:

- **Stock-grid fixtures**: legal masks (pass `stockgrid_validate`) at 9x9,
  15x15, 21x21, realistic block density/symmetry.
- **Dictionaries** [RT — inverted from draft]: vendor a REAL public-domain
  wordlist (ENABLE/SCOWL family) as the PRIMARY rung dictionary — fill
  difficulty is crossing-completion structure, which random synthetic words
  do not model (they inter-complete far too rarely; this is what saturated
  the 5x5 probe). A seeded synthetic generator is a SIZE-SCALING knob only.
  Byte-freeze whatever is committed (load-layer counts depend on exact
  bytes). ASCII-only for now: `normalize_word`'s `char_type(alpha)` +
  `string_upper` are locale/Unicode-fragile (fill.pl:129-131) — note at the
  fixture level, fix out of scope.
- **Fixture traps** [RT]: `\+ memberchk(Word, Used)` forbids the same word
  across AND down (fill.pl:196) — planted fixtures must have asymmetric
  solutions and respect the no-duplicate rule.
- **`benchmarks/run_fill.pl`**: FOUR attribution buckets [RT]: command
  (process), dict-load (`inproc` on `load_dict/3`), search (`inproc` on the
  budget-explicit `fill_attempt/8` with pre-loaded dict), and grid/slot
  derivation (`fill_grid/4`) either as its own bucket or explicitly inside
  "rest" — never silently absorbed. The seams already exist (fill.pl:243,
  :114, :53).
- **Budget story** [RT]: bench drives `fill_attempt/8` (budget-explicit);
  the /7 arity hardcodes `fill_budget(800_000_000)` (fill.pl:46) with NO
  command-layer override — add a `--budget` CLI hook (or reconcile the
  hardcode with arrange's shipped 500M) as a small pre-campaign fix.
- **Ratchet** [RT — fork, don't parameterize]: `check_fill_baseline.pl` +
  `fill_baseline.json` forked from the arrange ratchet. The coupling is the
  schema, not paths: fill wants TWO measured layers (search_inf gated;
  load_inf reported, gate-decision explicit) and different rung semantics.
  Keep the two ratchets independent.
- **Ladder**: grid x dictionary x seed-density rungs spanning ~0.1M-100M
  search inferences, chosen by measured cost. Expect this to be Phase 0's
  hardest problem [RT]: the proven-fast/saturates cliff is steep and
  dictionary-content-dependent; every rung must complete deterministically
  with headroom (the arrange P2 lesson, harder here).
- **A scale golden** (15x15 fill) locking byte-identity.
- Acceptance: +0.00% reproduction across repeated runs; warmup/JIT wobble
  characterized as in arrange.

## Phase 1 — instrumentation probe P-F1 (gate: no experiment is built first)

- Attribution: load vs slot-derivation vs search; within search,
  materialization (`candidates/4`) vs counting (`ord_intersection`) vs
  unification vs used-list scans. [RT: the review's micro-benches predict
  materialization leads; P-F1 confirms on real rungs.]
- Backtrack profile per rung (depth, churn) — quantifies the heavy-tail
  thesis on the actual ladder.
- Fan-out stats: slots, candidates/slot over time, index-set sizes per
  dictionary scale.

## Phase 2 — experiments (serial, ratcheted, byte-identity-gated; order set by P-F1)

- **F-H3 — buckets as compound terms + O(1) candidate access** [RT —
  promoted to lead candidate]: replace list buckets + per-candidate `nth0`
  with `arg/3` on compound-term buckets. Expected: the largest win at
  dictionary scale.
- **F-H5 — lazy candidate generation for the selected slot** [RT — new]:
  generate candidates on backtracking from the index instead of
  materializing the full list `member/2` will mostly discard. Must preserve
  enumeration order exactly.
- **F-H1 — incremental candidate counts**: carry per-slot counts in
  threaded state; recount exactly the slots sharing a NEWLY-BOUND cell with
  the just-placed word — INCLUDING slots thereby completed by crossings,
  which stay in the counted set (a 0-count completed slot must remain
  selectable so MRV fails the branch; dropping it diverges) [RT — sharpened
  wording]. Counts never read `Used` (try-time only, fill.pl:196), so
  carried counts are exact and the tree is preserved by construction.
- **F-H2 — bitset index for COUNTING ONLY** [RT — scope cut]: bigint per
  (len,pos,char) index set; count = `popcount(A /\ B)` (measured 250-600x
  on the count). Candidate ENUMERATION stays on ordsets/buckets — bit
  iteration measured ~7x SLOWER than an ordset walk (each `lsb` rescans the
  bignum, each clear allocates). Keep both representations or convert only
  the winning mask. GATE [RT]: measure bignum AND+popcount under
  WASM/LibBF (GMP is optional there) BEFORE building — this is the one
  change most likely to be a native win and a WASM loss. Distinction from
  arrange's rejected E-H8 letter bitmask, recorded so the ledger doesn't
  look self-contradicting: E-H8 gated a cheap O(len) check with ~2-14% hit
  rates; F-H2 replaces the EXPENSIVE O(dictionary) intersection itself,
  with a measured 250-600x kernel speedup — different economics entirely.
- **F-H4 — Used as a set** (O(placed x wordlen) memberchk per try today;
  compounds on deep trees).

## Phase 3 — startup + WASM posture

- Slot-length-filtered index build (only lengths the mask has); lazy
  per-length build if attribution supports it.
- Dictionary-load strategy for the browser (precompiled artifact aligned
  with the .qlf plan). Measured by the load bucket from Phase 0.
- If P-F1 shows load_dict > ~50% of end-to-end latency at realistic scale,
  Phase 3 outranks Phase 2 for perceived performance.

## Explicit non-goals

- Envelope/completability work (restarts, AC-strength pruning — both
  output-changing or out of scope; see thesis).
- Any change to fill's output or infeasible/not_proven semantics.
- Unicode/locale hardening of dictionary normalization (noted, deferred).

## Execution model

As the arrange campaign: orchestrator adjudicates; one Opus agent per
phase/experiment in an isolated worktree (verify base commit first — three
stale-worktree incidents last campaign); serial composition against
re-recorded baselines; every accept/reject gets an experiments.md entry;
research findings land as durable docs. Two do-not-skip gates from the
review: P-F1 before any experiment; WASM bignum measurement before F-H2.
