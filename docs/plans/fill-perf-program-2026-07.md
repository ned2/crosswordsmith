# Fill performance program (2026-07): DP-9 follow-ups + Balafoutis candidates

A program of work for exploring fill performance uplifts, folding together
(a) the two engineering rows surfaced by the DP-9 grounding runs and (b) the
six search-power candidates from the Balafoutis thesis review. Status:
PROPOSED — nothing here is spec'd or locked; Track D is where any contract
change would be decided.

Evidence base:

- `benchmarks/fill_quality/README.md`, "CWL on the §8.4c engine (2026-07-16)
  — DP-9 grounding" — the two stack envelopes, the quality table (§8.4c+CWL
  beats ingrid_core on all five completable masks), and the list-specificity
  of the DP-8 reference-row completion (STW yes, CWL no at any loadable
  floor).
- `docs/design-spec.md` §8.5 — the capacity-robustness row and the seeded
  O(n²) load-walk row, each with numbers and fix-shape constraints; §10 DP-9
  for why neither was built in that pass.
- `docs/plans/fill-mac-dwd-implementation.md` — the §8.4c build, the C3
  findings (lazy candidates, diversified restarts, pinned load shuffle,
  and where the seeded walk's O(n²) form comes from), and the measured
  baselines.
- `docs/research/balafoutis-thesis-review-2026-07.md` — candidates F1–F6
  with thesis-section grounding.

## Standing constraints (read before touching anything)

1. **Byte-identity is contract (AC-FILL-3/13).** Any landed change to fill's
   search or load must keep `make test` goldens,
   `benchmarks/fill_identity.sha256` (11 rungs), and the inference ratchet
   (`benchmarks/check_fill_baseline.pl`) diff-clean. "Faster but different
   fills" is a **decision pass**, not a patch (design-spec §10, OD-8).
2. **§8.4c / DP-8 is LOCKED.** The MAC + dom/wdeg + restarts core is not
   relitigated here; every Track B/C item is a refinement measured *against*
   it, adopted (if ever) as one amendment pass.
3. **DP-6 pins stand.** Budget buys latency, never completion (measured ×20);
   `--min-score` is never presented as a completion lever;
   `blocked_13b`/`blocked_15a` are report-don't-chase — Track C's spike is
   the one sanctioned, explicitly-framed exception and its go/no-go must
   weigh that pin.
4. **Probes live outside the default path.** The DP-7 pattern
   (`probe_mac_dwd/6`: a throwaway variant rig in benchmark-land calling
   engine internals, default path untouched) is the required shape for all
   Track B/C measurement. Threading strategy flags through the shipped
   `mac_*` hot path before a GO decision risks ratchet churn for zero user
   value.

## The framing that shapes the program

Two axes classify every item:

- **Identity-churning vs identity-clean.** All six thesis candidates change
  learned weights or branching order, hence which fill is found, hence every
  golden/identity/ratchet artifact — so none can land incrementally. The
  economical shape is the one DP-7→DP-8 already proved: measure everything
  behind probes, then batch the winners into **one** engine amendment
  (Track D) so the churn (engine version bump, golden regen, identity +
  ratchet re-baseline, wasm parity) is priced **once**. The two §8.5 rows
  (Track A) are the opposite: user-facing paper cuts fixable with zero
  identity churn — land them first and independently.
- **Latency vs envelope.** Both prior campaigns converged on the same
  decoupling: node-cost and even ×20 budget never flip completion. F1–F5
  are latency/search-efficiency bets; **F6 is the only envelope bet** in the
  thesis basket. The program keeps them in separate tracks with separate
  success metrics so a latency win is never oversold as search power.

## Track A — contract-safe engineering (land first, no decision-pass batching)

### A1. Capacity robustness at the two stack envelopes (§8.5 row) — effort S/M

Facts: (a) *load* — the full 566k cleaned CWL blows the default 1GB stack in
`build_index`; (b) *search* — the ≥30 band (437,400 words) loads but blows
the global stack inside `mac_support`/`mac_revise` bignum mask arithmetic
~17s into full-13 search (`blocked_13a`). Both die as raw overflows, not the
§5.1 failure contract (INV-3 tension). Reproducer:
`scripts/fetch-cwl.sh --min-score 30` + `fill --grid grids/blocked_13a.json`.

Shape: a micro decision pass (OD-8) choosing between **clean capacity
error** (catch/report at the measured envelopes; document the envelope in
README + `dicts/README.md`) and **streamed/segmented index build** (raises
the load envelope for real). Recommendation: take the clean error first —
small, identity-clean, converts both crashes into ordinary §5.1 outcomes —
and leave segmentation as a stretch item that only activates on a concrete
user need (per the best-effort STOP discipline). Note the search-side
envelope is mask-*width* (bits per length-bucket), so a segmented-mask
representation is the only real fix there; that is engine-internal surgery
and, if ever done, belongs in Track D's batch (it must re-verify
`make test-wasm` — LibBF bignum parity is part of the battery).

Side benefit: A1 unblocks Track B/C probing at wide bands (today a probe at
`--min-score 30` on full-13 grids just crashes).

### A2. Seeded load-walk O(n²) (§8.5 row) — effort M

Facts: `--seed`/`--shuffle` on `dicts/cwl50.dict` costs ~2m45s **at load**
regardless of grid (85,800-entry score-50 band ≈ 7.4G walk steps; a trivial
3×3 pays the same 2m43s). The default path's pinned shuffle is already an
O(n log n) keysort (`prng_shuffle/2`); the seeded path kept
`seeded_permutation/2`'s selection-walk form because the draw *sequence* is
DP-6 contract.

Shape: the §8.5 row's micro decision pass, with the analysis question posed
sharply: **is the existing draw sequence replayable at O(n log n)?** If each
walk step draws an index into the *remaining* items and selects the k-th
survivor, an order-statistic structure (e.g. Fenwick-indexed selection)
consumes the identical PRNG stream and emits the identical permutation —
keep-sequence-optimise, zero golden churn, the recommended arm. Only if the
walk's draw consumption turns out not to be replayable does the re-contract
arm (one-time seeded-fill churn, `fill_3_rng_seed7` regenerated, AC-FILL-10
wording) get priced. Acceptance either way: seeded goldens byte-identical
(keep-sequence) or regenerated under an explicit decision record
(re-contract); `--seed` on cwl50 drops from minutes to seconds.

## Track B — Balafoutis probe campaign (measurement only; no default-path change)

All items run on the DP-7-style probe rig (extend the `probe_mac_dwd`
lineage in benchmark-land). Standard measurement matrix for every arm, in
descending order of authority:

| Signal | Source | What it arbitrates |
|---|---|---|
| Reference row completion + quality (AC-FILL-12) | `blocked_13a` × STW `@30`/`@1` (baseline: 2m20s/7 attempts, mean 45.0 / 19s, mean 38.7) | search power + §8.4a quality contract |
| Ladder inferences/nodes/attempts | the 11 identity rungs' workloads (read-only — probes don't touch the sha256 file) | latency; regression watch on easy rungs |
| CWL clean-floor row | `blocked_13a` × `cwl50` @50 (baseline: not-proven, 3m57s; defeats ingrid too) | open stretch signal — any completion here is news |
| Quality gate | `benchmarks/fill_quality/run.sh` five-mask table | no quality regression vs the DP-9 table |

### B0. Instrumentation probe first — effort S, do before any variant

The arrange campaign's single best move (P1) was measuring before building;
its E-H6 rejection shows how a theoretically-sound weight signal dies on an
unmeasured structural fact. Fill needs the equivalent numbers, all cheap
counters on the probe rig:

- **Worklist length distribution** in `mac_propagate/7` at pop time. If
  worklists are almost always length ≤ 2, F1 (revision ordering) is dead on
  arrival — this is the kill-test.
- **Revision outcome split** (redundant / fruitful / DWO) — the thesis
  measured 96.5% / 3.3% / 0.2% on its benchmarks; ours will differ and the
  fruitful fraction sizes F4's addressable surface.
- **Weight concentration**: do conflict weights concentrate on few edges
  (learnable structure) or smear (F2/F3 have nothing to sharpen)?
- **Where wall-time goes at wide bands** (mask width vs node count) — sizes
  whether Track B latency wins matter at CWL scale or the cost is all
  bignum arithmetic (which only A1-adjacent representation work touches).

### B1. F1 — weight-ordered revision queue (thesis Ch 5) — effort S, first variant

Pop the queued slot with smallest dom/wdeg (`v_dom/wdeg`; arms for `v_wdeg`
and descending-weight edge order within a revision). Gated by B0's worklist
kill-test. Thesis evidence: up to ~5× cpu swings and node-count cuts from
revision order alone under a conflict-driven VOH.

### B2. F2 + F4 — weight-credit variants (thesis §3.4.2/3.4.4/3.4.6) — effort S/M

One multi-arm experiment: H1/H2/H3 (deletion-responsible credit at DWO),
`alldel` (credit on every deletion), fully-assigned (credit fruitful
revisions of the DWO episode) — each with and without B1's ordering, since
the thesis's §5.5 variance study says credit scheme and revision order
interact (alldel/fully-assigned are order-robust; plain dom/wdeg is not).
Thesis evidence: plain dom/wdeg (ours) won zero instances outright in its
tables.

### B3. F3 + F5 — aging and probing-init arms — effort S

F3: intra-attempt aging (÷2 every ~20 backtracks) vs our between-attempt
×0.99. F4-style caveat applies: high per-class variance in the thesis (5×
win to outright worst), a pure bench gamble. F5: pinned probing runs to warm
weights before the greedy attempt 1 (which currently searches with uniform
weights, i.e. degraded dom/wdeg). Both are small arms appended to B2's rig.

Stop conditions for Track B: an arm that loses or ties on the reference row
AND the ladder is closed with its numbers recorded (experiments-ledger
discipline); no arm graduates on a single-rung win. If B0's kill-tests fire
(short worklists, smeared weights), record and close the corresponding arms
without building them — that IS the completed outcome.

## Track C — the envelope bet

### C1. F6 — set branching on crossing-letter projections (thesis Ch 7) — effort M/L

The only candidate aimed at completion rather than latency: candidates
agreeing on all checked cells are propagation-interchangeable; branch on the
checked-cell projection and refute whole classes. Classes are large exactly
on UK blocked grids (~half the cells unchecked) — the shape of every
currently-uncompletable row (`blocked_13b`/`blocked_15a` at all thresholds;
`blocked_13a` × CWL @50).

Discipline: this deliberately probes grids under the DP-6 report-don't-chase
pin, so it runs as a **bounded throwaway spike** (the DP-7 shape): fixed
effort box, probe rig only, and it graduates to a Track D proposal *only* on
a completion at spike node rates — the same falsification bar
`probe_mac_dwd` cleared before DP-8. Design wrinkles to resolve in the
spike: the all-different filter and final word scoring distinguish class
members (resolve at solution assembly); class enumeration must stay lazy
(the C3 lesson — eager materialization is how the C2 core died on open
grids). If the spike fails, the recorded outcome is "F6 measured, does not
close the blocked rows" and the pins stand unchanged.

## Track D — the adoption decision pass (DP-10 candidate)

One batched go/no-go over every Track B/C arm that beat baseline, priced as
a single §8.4c amendment: engine version bump, golden regeneration,
`fill_identity.sha256` re-baseline, ratchet re-baseline, `make test-wasm`
parity, seeded-goldens audit. Inputs it must carry (the DP-8 template):
completion evidence per arm, quality table vs the DP-9 baseline, determinism
posture (all arms must be RNG-free or pinned-stream — nothing in F1–F6
requires new randomness), and honesty pins for anything an arm does *not*
buy. NO-GO is a fine outcome: the evidence attaches to the §8.5 backlog the
way the DP-8 probe evidence did, and the default engine stays byte-stable.

## Sequencing

1. **A1** (capacity error) — smallest item, converts crashes to contract
   outcomes, unblocks wide-band probing.
2. **A2** (seeded walk) — micro decision pass, then the keep-sequence
   rewrite if replayable.
3. **B0** (instrumentation) — cheap, and its kill-tests may delete half of
   Track B before anything is built.
4. **B1 → B2 → B3** as surviving arms, on the shared probe rig, each with
   recorded numbers.
5. **C1** in parallel with B2/B3 once B0's rig exists (it reuses the
   counters).
6. **D** only when the arms are exhausted or a clear winning basket exists.

Tracks A and B/C are independent; A can land while B is still measuring.

## Out of scope (locked or measured-dead)

- Relitigating §8.4c/DP-8, the §8.4a score-first quality contract, or the
  fillability value order (measured quality inversion, DP-8 honesty pin).
- Default-path/default-dict changes (DP-9 decision: defaults untouched).
- Budget-default changes (DP-6: ×20 budget bought no completion).
- Thesis ideas already skipped with reasons in the research note (adaptive
  2-way branching, impact heuristics, freeze/undo multi-DWO, generic
  clustered set branching).
- Arrange: this program is fill-only; the arrange verdict and the
  arrange-as-CSP conceptual exercise live in the research note.
