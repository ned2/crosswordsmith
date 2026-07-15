# crosswordsmith — Design Specification

**Status: normative.** This is the single source of truth for *what* crosswordsmith is, what it is not, and how you know a piece is done. It exists to stop the program growing by accretion: nothing gets built that is not specified here, and this document only grows by a **deliberate decision pass**, never organically during implementation.

## How to read this document

- Every component is tagged **LOCKED**, **PARTIAL**, or **DEFERRED**:
  - **LOCKED** — decided; has acceptance criteria; ready to implement against.
  - **PARTIAL** — direction decided, some sub-decisions still open (flagged inline + in the register).
  - **DEFERRED** — design *intent* recorded and bounded, but **not buildable** until a future decision pass resolves its open decisions (see §10). Writing code for a DEFERRED component is itself scope creep.
- Acceptance criteria are written as observable, testable statements with stable IDs (`AC-<area>-<n>`). They are the test contract: a component is "done" when its ACs pass.
- Keywords **MUST / SHOULD / MAY** carry their usual normative force.

### Relationship to the other docs

| Doc | Role | Layer |
|---|---|---|
| [`cryptic-setter-research.md`](./cryptic-setter-research.md) | **Why** — discovery, evidence, rationale (exploratory) | Cited *from* here |
| **`design-spec.md`** (this) | **What** + acceptance criteria (normative) | Source of truth |
| [`arrange-implementation-plan.md`](./arrange-implementation-plan.md) | **How / phasing** for the `arrange` engine (derived) | Implements §7 |

When this spec and the research doc disagree, **this spec wins** (the research doc is a snapshot of exploration; decisions made since supersede it). Implementation plans derive from this spec, never the reverse.

### Change discipline

1. A change to a **LOCKED** contract or acceptance criterion requires editing this doc *first*, in its own commit, with the rationale.
2. New scope (a new verb, flag, feature, or output key) MUST appear here before it is implemented. "We can just add a flag" is the exact failure mode this doc exists to prevent.
3. Promoting a **DEFERRED** component to **LOCKED** requires a recorded decision pass resolving every open decision in §10 for that component.
4. Removing scope is also a change: move it to §3 (Non-goals) with a one-line reason rather than deleting silently.

---

## 1. Vision & scope

crosswordsmith is a **free, cross-platform, scriptable, open-format, CLI-first crossword layout-and-fill engine** that emits deterministic, diffable output. It feeds clue-writing tools; it is not one. (Research doc, Executive Summary §8; Design principle #1.)

It is delivered as **one CLI with two product surfaces** ("flavours"), sharing a common substrate but **not** a common placement engine:

- **Flavour A — `arrange`.** "Arrange *these specific* words into a nice crossword-shaped layout." Closed word set; aesthetic interlock; does **not** aspire to authentic cryptic legality (no symmetry guarantee, no ≥50%-checking guarantee). This is the program's genuine niche — no Exet/Crossword-Compiler equivalent exists.
- **Flavour B — cryptic setting workflow tools.** Validation (`lint`), standard-format interchange (`export`), stock-grid libraries/profiles, and — the one heavy piece — a grid-first, open-dictionary auto-`fill` engine. Flavour B is *grid-first* (legality is designed into the black-square pattern up front), the opposite inversion from Flavour A's closed set.

**Why two engines and not one with flags:** the legality cores are structurally incompatible. Flavour A is closed-set-emergent (the grid emerges from how the supplied words interlock); Flavour B `fill` is grid-first-open-dictionary (a pre-validated grid, words drawn from a ~270k lexicon). Forcing one engine to serve both is the trap. (See research doc, "The closed-set vs. authentic-layout tension".) **Share the substrate, not the solver.**

---

## 2. Build sequence & status at a glance

Sequence: **A → B-`lint`/`export` (cheap shared wins) → B-`fill` (blocked) → B-`fill` (barred, furthest out, may never ship).** This table is a reading aid; the per-component tags in §§5–8 are authoritative.

| Component | Section | Tag | State |
|---|---|---|---|
| CLI contract | §5 | LOCKED | spec'd; built at the Phase-7 cutover |
| Shared substrate | §6 | LOCKED | mostly exists in `crossword.pl` |
| `arrange` (Flavour A) | §7 | LOCKED | + build plan ([`arrange-implementation-plan.md`](./arrange-implementation-plan.md)) |
| `lint` | §8.1 | LOCKED | toc/blocked-uk/american + barred-ximenean (OD-7 resolved: Ximenean band) |
| `export` (ipuz/Exolve) | §8.2 | LOCKED | transformations of canonical JSON |
| Stock-grid library / profiles | §8.3 | LOCKED | mask schema (OD-5) + grid set (OD-6) resolved; ships 3 lint-validated grids |
| `fill` (Flavour B engine) | §8.4 | LOCKED | OD-1…4 resolved (DP-1/DP-2); blocked-only, in-memory index, stock-grid profiles |
| Backlog features | §8.5 | unspec'd | each needs its own decision pass |

---

## 3. Non-goals / out of scope

These are **deliberately excluded**. Adding any of them requires a §10 decision pass and a scope change here, not a quiet PR.

| Excluded | Reason |
|---|---|
| **Auto clue-writing / surface-reading judgement** | Even Exet and Wordplay Wizard stop at presenting candidates. Be an assistant, not a generator. (Principle #3.) |
| **Hand-rolled binary `.puz` writer** (CRC-16 in Prolog) | Brittle, and `.puz` cannot represent bars. Reach `.puz`/`.jpz`/PDF via `kotwords` from the native `.ipuz` emitter. |
| **Solver-facing consumption** — HTML5 applet, hosting, embedding, solve analytics/heatmaps, blogs | All require a hosted backend; orthogonal to a CLI. |
| **Bundling paywalled dictionaries** (Chambers/Collins/Oxford) | License. Ship only freely-redistributable data — MIT/CC0/CC-BY word lists, or UKACD18 (J Ross Beresford's *redistributable freeware* lexicon; **not** BSD-3 — if bundled, ship its copyright notice + license text verbatim). (Principle #9.) |
| **Rectangular / non-square grids** (8×10 New Yorker, etc.) | Cut from the engine; square-only `--size N`. Revisit only if a Flavour-B engine ever needs it. |
| **`free`/auto-size mode** for `arrange` | Replaced by `--max-size N` (ceiling + crop). No size-search outer loop. |
| **Priority / per-word importance scores** in `arrange` | De-scoped; anchoring is handled by the fragment grid instead. |
| **Barred-grid *generation*** | Different cell/edge model entirely — a separate engine, deferred indefinitely. The barred checking-math *lint* is a cheap exception (§8.3). |
| **Real-time co-editing / community voting platform** | Large product surfaces orthogonal to a CLI. Git-diffable source is the collaboration story instead. |
| **Selection-from-an-open-pool in `arrange`** | Flavour A is closed-set by definition; open-dictionary fill is exclusively Flavour B `fill`. |

---

## 4. Architecture — shared substrate vs. divergent tops

The discipline: **one substrate, two solver tops.**

| Layer | Shared | Flavour A (`arrange`) | Flavour B |
|---|---|---|---|
| Input parsing / `meta` passthrough | ✅ | | |
| Blocked-grid cell/coordinate model | ✅ | | |
| Clue numbering + enumeration | ✅ | | |
| Metric predicates (`word_meets_half/2`, `max_unch_run`, `checked_cells`, connectivity) | ✅ | used as **optimizer signals** | used as **validators** |
| Emit layer (JSON canonical; later ipuz/Exolve) | ✅ | | |
| Fragment-grid seeding primitive | ✅ (one primitive) | anchors | ninas / seed-then-fill |
| **Placement engine** | ✗ | unified B&B optimizer (new) | grid-first + dictionary fill (new, separate) |
| Legality core (`adj_is_free`, `no_word_merge`, `check_prev/next_cell`) | ✗ | free-canvas (blocked) | template-slot / barred |

**Module layout (target):** an SWI-Prolog library layout under `prolog/crosswordsmith/`, one module per file (flat-namespace module names prefixed `crosswordsmith_*`), with explicit export lists. Migration is staged in [`source-structure-migration-plan.md`](source-structure-migration-plan.md).

- `prolog/crosswordsmith/core.pl` (`crosswordsmith_core`) — shared substrate: input, grid model, legality core for the free-canvas case, clue numbering, emit.
- `prolog/crosswordsmith/metrics.pl` (`crosswordsmith_metrics`) — the shared metric predicates (`crossing_count`/`placed_bbox`/`word_meets_half`, checked-cell/unch-run helpers); today's `quality.pl` renamed, minus the greedy constructor. Metrics remain a **separate module** — they are *not* lifted into core (supersedes earlier language here): `lint` depends only on the JSON contract plus metrics, never on the solver substrate, and that boundary is deliberate and preserved.
- `prolog/crosswordsmith/arrange.pl` (`crosswordsmith_arrange`) — Flavour-A engine, **including the greedy constructor** (moved from `quality.pl`): its only consumer is arrange, and its project dependencies are core's search primitives, not the metrics.
- Flavour-B modules — `prolog/crosswordsmith/{lint,export,stockgrid,fill}.pl` (`crosswordsmith_lint`/`_export`/`_stockgrid`/`_fill`).
- Root `load.pl` — search-path alias + single owner of load order; the CLI driver, tests, and benchmarks load the project through it.
- `quality.pl` — retired by the rename to `metrics.pl` (the stale `--quality`-era name goes away; the metric layer itself stays).

**Cross-cutting invariants (apply to all components):**
- **INV-1 — Metadata-agnostic engine.** Clue text, annotations, cluing scores, review flags ride in an opaque `meta` object joined to the answer at emit time. They MUST NOT be threaded through placement or clue-numbering. The JSON contract changes only by *additive* keys. (Principle #2.)
- **INV-2 — Deterministic by default, diffable output.** Output is stable sorted-key JSON; under the default flags, identical input ⇒ byte-identical output. The *only* source of randomization is the opt-in `--seed N` / `--shuffle` perturbation of the strict `arrange` search (§7.6): it never touches the default flow, and a fixed `--seed N` is itself reproducible (identical input + `--seed N` ⇒ byte-identical). `--shuffle` deliberately draws a fresh seed each run and reports it under `--verbose` so a liked layout can be reproduced with `--seed`. (Principle #8.)
- **INV-3 — Report, never silently degrade.** Dropped words, layout compromises, lint failures, and infeasible seeds are always reported, never swallowed. (Principle #5.) Reporting channel: per-layout quality caveats ride the emitted payload's `diagnostics` property (json-output-spec §6.4); failures and request-level shortfalls report on stderr (§5.1).
- **INV-4 — License-clean by default.** Only permissively-licensed bundled data. (Principle #9.)

---

## 5. CLI contract — explicit verbs, no catch-all default  **[LOCKED]**

Every capability is a **verb**. A bare invocation prints usage and exits non-zero — it never performs an action.

```
crosswordsmith arrange [--strict | --best-effort] [--size N | --max-size N]
                       [--fragment grid.json] [--candidates N] [--enumerate]
                       [--check-target N] [--seed N | --shuffle]
                       [--verbose] --input words.json [--out file.json]
crosswordsmith lint    --profile <name> [--allow-asymmetry] layout.json     # Flavour B
crosswordsmith export  --to ipuz|exolve  layout.json [--out file]           # Flavour B / shared
crosswordsmith fill    --grid template.json --seeds seeds.json [--verbose]  # Flavour B (DEFERRED)
crosswordsmith                                                              # → usage, exit ≠ 0
```

### 5.1 Conventions
- `--input <file>` is the word/clue source (JSON or `.pl` fixture), required by `arrange`.
- `--out <file>` writes to a file instead of stdout; composes with any verb. Absent ⇒ stdout.
- `--flag=value` and `--flag value` are both accepted.
- Unknown flags/verbs are an error with a usage hint, not silently ignored.
- **Quiet on clean success; `--verbose` opts into summaries.** A successful run
  writes nothing to stderr by default; `--verbose` (on `arrange`/`fill`) turns
  on the routine success summary (grid, placed, reward / filled slots) so the
  engine composes silently in pipes. **Quality caveats ride the payload, not
  stderr**: arrange output is best-effort by nature (most real clue sets fall
  short of perfect health), so per-layout compromises — dropped words
  (AC-ARR-2), cap-inert degeneration (§7.2) — are reported in the emitted
  payload's `diagnostics` property ([`json-output-spec.md`](./json-output-spec.md) §6.4) where consumers
  can act on them, satisfying INV-3 without polluting the terminal (`--verbose`
  echoes them on stderr for interactive runs). Request-level shortfalls with no
  payload home (fewer-than-K candidates, AC-ARR-7) and every failure report
  remain **unconditional** on stderr per INV-3. Stdout data is identical either
  way.

### 5.2 Exit codes (LOCKED)
| Code | Meaning |
|---|---|
| `0` | Success: a layout/report was produced (for `lint`, PASS or WARN only), or an explicit `--help`/`-h`/`help` request (usage to stdout — GNU/POSIX convention). |
| non-zero | Failure: unplaceable input under `--strict`, unsatisfiable fragment, lint FAIL, malformed input, or a bare/unknown-verb invocation (usage to stderr). |

A non-zero exit with `--out` MUST NOT write a partial output file.

### 5.3 Migration from the current CLI (LOCKED, one deliberate pass)
The current interface — `./crossword.pl --input F [--strategy S] [--shuffle] [--all] [--quality ...] [--out O] <grid_length> [<start_loc>]` — is replaced by subcommands. This is a **breaking change**, done as one pass:
- old `solve` (fixed-grid, place-all) ≡ `arrange --strict --size N` (exact N×N is the default framing)
- old `--quality` (engine-size, drop) ≡ `arrange --best-effort --max-size N`
- old `--all` (enumerate) ≡ `arrange --enumerate`
- The required `GridLen` positional becomes `--size N` (resolving the current wart where `GridLen` is required but silently ignored under `--quality`).
- old `--shuffle` / `solve_shuffled` (a whole randomized *strategy*) is **removed**; a later, narrower `--seed N` / `--shuffle` opt-in perturbs only the strict search's branch order, leaving the deterministic default flow untouched (INV-2, §7.6).
- Old-style invocations are not specially handled: they fall through to the standard unknown-verb usage error (the dedicated "did you mean `arrange`?" migration hint has been retired now that there are no pre-migration consumers).
- README, `run_tests.sh`, and golden fixtures are updated in the same pass.

**AC-CLI-1** Bare `crosswordsmith` prints usage and exits non-zero.
**AC-CLI-2** Each documented verb/flag combination behaves as specified or errors with a usage hint; no undocumented flag exists.
**AC-CLI-3** No old-style invocation silently succeeds: a numeric-first or flags-first invocation errors with the standard usage hint and a non-zero exit; no silent behavioural change. (The §5.3 table documents the one-time mapping to `arrange` equivalents.)

---

## 6. Shared substrate  **[LOCKED]**

### 6.1 Input & metadata passthrough
- Accepts a word/clue set as JSON or `.pl` fixture; each entry is `{answer, meta?}` where `meta` is opaque (clue text, url, scores, flags).
- Answers MUST be unique within a set (existing `check_unique_answers/1`); duplicates are an error.
- `meta` is joined to the answer **at emit time only**; it never influences placement (INV-1).

**AC-IN-1** A set with duplicate answers is rejected with a clear error.
**AC-IN-2** Arbitrary `meta` keys round-trip from input to emitted output unchanged.

### 6.2 Grid model
- Square N×N canvas; assoc `cell(1..N²) → empty | Letter`; row-major numbering; geometry via `cell_coord/3`, `next_cell`/`prev_cell`/`calc_num`/`calc_start`/`fits_on_grid`/`is_start_cell`/`is_end_cell`.
- No `(W,H)` plumbing (square-only, §3).

### 6.3 Clue numbering & enumeration
- Standard crossword numbering via `assign_clue_numbers/2`.
- **Enumeration** (e.g. `(4,3)`, `(4-2)`) is derived from preserved spaces/hyphens in the answer and carried on the emitted word. Required by every clue-bearing export format.

**AC-NUM-1** Numbering matches conventional crossword rules for any legal layout.
**AC-ENUM-1** A multi-word/hyphenated answer emits the correct enumeration string.

### 6.4 Metric predicates (measurement, shared)
`crossing_count/6`, `placed_bbox/4`, `word_meets_half/2` (`checked ≥ ceil(L/2)`), `checked_cells`, `max_unch_run`, connectivity. Flavour A consumes these as **optimizer signals**; Flavour B `lint` consumes them as **validators**. They are pure measurement — no side effects, no placement.

### 6.5 Emit / output contract (canonical JSON)
- Canonical output is the existing stable sorted-key JSON: `{gridLength, grid: rows, words: [{number, direction, answer, cells, meta}...]}` (the authoritative per-key contract is [`json-output-spec.md`](./json-output-spec.md) §6.2; `dir`/`num`/`len`/`start` are *internal* placed-word fields and are not emitted).
- Output is deterministic (INV-2) and the **canonical fragment schema** (§6.6) — i.e. emit is round-trippable back in as input.
- Export formats (ipuz/Exolve) are *transformations of* this canonical form (§8.2), not separate emitters.

**AC-EMIT-1** Identical input + flags ⇒ byte-identical JSON across runs and machines.
**AC-EMIT-2** Emitted JSON re-ingested as a `--fragment` reproduces the identical layout (round-trip).

### 6.6 Fragment-grid primitive
A **fragment grid** is a partial layout the engine solves *from*. It is the single mechanism behind three features: **anchors** (Flavour A — pre-placed words), **ninas** (Flavour B — fixed letters along a path), and **interactive iteration** (place some, let the engine finish).

- **Schema = the emit format made partial**, with the convention **presence = fixed**: any word/cell present is pinned; everything absent is the engine's to place. There is no separate "anchor" field — *being in the fragment is the anchor*.
- An optional **thin, hand-authorable convenience form** (e.g. `{answer, row, col, dir}` tuples + optional fixed cells) **desugars into** the canonical form. *(Implemented for `arrange --fragment` and `fill --seeds`, words-only: a top-level JSON list of `{answer, row, col, dir}` entries, 0-based start cell. It carries no `gridLength`; the consumer frames the desugar — arrange by `--size`/`--max-size` (default 15), fill by its own grid's side, so thin seeds land on exactly the grid they pin.)*
- **Words and letters both supported.** A placed word = a run of fixed cells of known identity; a lone fixed cell = a nina constraint. Words-only is a valid v1; the schema *permits* letter-level so Flavour-B ninas drop in later. Fixed letters not belonging to a list word are **free nina constraints, not entries** — they do not count toward `--strict`'s "place all words".
- **Reconciliation by answer string** against `--input`; the engine places the unmatched remainder. A fragment word absent from `--input` is an error.
- **Validate up front** against the legality core (overlaps, adjacency/merge, unsatisfiable fixed letters) and report conflicts *before* searching.
- **Size frame:** the fragment carries grid dimensions, which set `N`; `--size`/`--max-size` are then redundant (if given and they disagree with the fragment, error). Which flag is used still frames the result — `--size` exact, `--max-size` cropped.

**AC-FRAG-1** A fragment word not present in `--input` is rejected with a clear error.
**AC-FRAG-2** A self-conflicting fragment (overlap/merge/unsatisfiable letter) is reported before any search runs.
**AC-FRAG-3** Pinned words/cells appear at exactly their fragment positions in every produced layout.
**AC-FRAG-4** Thin convenience form and canonical form for the same fragment produce identical results.

*(Resolved by DP-1 (§10): the thin convenience form was **deferred** — v1 shipped canonical-only (the emit format made partial). The deferral has since been lifted in two steps: first for `arrange --fragment` (satisfying AC-FRAG-4: plunit + a golden check that byte-compares the two forms against the same golden), then extended to `fill --seeds` (same byte-identity shape, desugared on the fill grid's own side). Duplicate `--input` answers are **rejected** (answers are unique, `check_unique_answers/1`), so no fragment duplicate-answer disambiguation is needed.)*

---

## 7. Flavour A — `arrange`  **[LOCKED]**

A single verb that unifies the old `pack`/`solve` engines into one **deterministic MRV-first layout engine** that maximizes a capped interlock objective over an N×N grid (with admissible-bound pruning + within-budget improvement where it pays). Rationale: research doc "Two flavours, one CLI"; build plan + the empirical design-exploration reframe: [`arrange-implementation-plan.md`](./arrange-implementation-plan.md). *(Earlier called a "branch-and-bound optimizer"; probes showed the bound prunes little and the cap is often inert at our densities, so that label oversold the search — see §7.3 and the plan's "Empirical reframe.")*

### 7.1 Contract
| Dimension | Decision |
|---|---|
| **Scenario** | General-purpose (any word list; website-TOC is one consumer). |
| **Drop contract** | Configurable, **default `--strict`**. Strict fails only on a *genuinely unplaceable* word (no legal intersection anywhere). `--best-effort` drops unplaceable words and reports them. |
| **Sizing** | Two mutually-exclusive square framings: `--size N` = an exact N×N grid (**default**, N=**15**, the dominant blocked-cryptic size); `--max-size N` = build up to N×N, then crop to the tight enclosing square (≤ N). |
| **Aesthetic goal** | High interlock + even/balanced. Explicitly *not* raw density (old `quality.pl`) and *not* target-shape. |
| **Weak slots** | Soft-penalty only: always place any legally-placeable word; stubs are penalised, never refused. No per-word hard floor. |
| **Objective** | Capped per-word reward + small interlock tail (§7.2). |
| **Outputs** | Stable single best by default; opt-in `--candidates N` returns *meaningfully distinct* layouts (§7.4). |
| **Word distinction** | Anchors via the fragment grid (§6.6); **priority scores de-scoped** (§3). |
| **Enumerate** | `--enumerate` retained — exhaustive count/emit of *all* feasible full placements (old `--all`). |

### 7.2 Objective (LOCKED form; constants empirically calibrated)
`maximize  Σ_w [ WCap · min(checked(w), target(w)) + WTail · checked(w) ]`
- `target(w) = ceil(L/2)` by default (reuses `word_meets_half/2`).
- Default integer weights `WCap:WTail = 5:1` ⇒ `ε = 0.2`. The **cap creates balance** (stops rewarding already-checked words, redirecting effort to laggards, anti-stub in a length-aware way); the **tail keeps an interlock gradient above the cap** and breaks ties toward crossier layouts.
- Integer arithmetic for deterministic tiebreaks (mirrors `placement_key`'s idiom).
- **Reachability caveat (load-bearing — empirically confirmed).** The cap only balances when `target` is attainable at the densities Flavour A actually produces. Probes confirm this bites on realistic inputs: toc_demo reaches `ceil(L/2)` on only **1/16** words (and `--min-half` makes it infeasible outright), while the dense quality_22 mesh reaches it on **19/22** — so the cap is *inert on sparse link-sets, active on dense meshes*. Where it is inert, the objective **silently degenerates to plain total-crossings**. ⇒ `--check-target` MUST be a **real tunable**, and `target` and `ε` MUST be **calibrated against the fixtures *before weights are locked*** (calibration fixtures: toc_demo = sparse/inert, quality_22 = dense/active), lowering `target` below `ceil(L/2)` where half is unreachable, and **reporting when the objective has degenerated**. This is a calibration/tunability obligation, not a change to the objective *form*. *(Implemented as `--check-target N` (§5): it lowers the per-word target to `min(ceil(L/2), N)` — N only ever lowers, never raises — and every arrange payload reports `capInert: true` in its `diagnostics` property (json-output-spec §6.4) when no placed word reaches its target, i.e. the objective has degenerated to total-crossings; the `--verbose` summary echoes the caveat on stderr per §5.1.)*
- **`maximin`/leximin are NOT the search driver** (non-decomposable, useless early bound, defeats pruning). Reserve them only as a tiebreak among final top candidates, or as an optional hard per-word floor (a cheap forward-checkable *constraint*, not the objective).

### 7.3 Engine properties
- MRV-first DFS over `crossword.pl`'s skeleton with incumbent tracking + admissible-bound pruning; branch step reuses `select_inc`/`find_intersecting_word`/`assign_word` unchanged. Whether the search layer *beyond the first incumbent* earns its keep is gated by a measurement (plan Phase 1.5) — probes suggest the bound prunes little at our sizes.
- The objective is **additively decomposable** ⇒ incremental per-placement update + cheap admissible bound for pruning.
- **Anytime with a node/inference-count budget (not wall-clock):** optimization is strictly harder than first-solution search, so the contract is **"best within budget,"** not "proven optimum." MRV-first finds a strong incumbent early; budget exhaustion returns best-so-far. The budget MUST be a node/inference count for determinism (INV-2); a wall-clock cap may exist only as a safety valve that marks output `truncated`. **Budget-exhausted ≠ unplaceable** (see AC-ARR-1).

### 7.4 Sizing & candidates
- **`--size N`** (default framing; bare `arrange` ⇒ `--size 15`): `N` is the *canvas*; emit a full N×N grid, uncovered cells as blocks.
- **`--max-size N`**: `N` is a *ceiling*; emit the tight crop ≤ N×N (auto-shrinks for small inputs). Implemented as the tight enclosing **square**, side `max(H,W)`, anchored at the content's top-left — the output schema is single-`gridLength` square (§6.5), so the crop is squared rather than rectangular.
- The two are **mutually exclusive** — one selects exact framing, the other cropped — and both are orthogonal to `--strict`/`--best-effort` (either can fail-strict or drop-best-effort when words don't fit N×N).
- **`--candidates N`**: default is a single deterministic best. `N` opts into alternatives selected for **diversity**. Distinctness is sourced from **constructor breadth (multiple deterministic seeds) + greedy diversity**, *not* from a single deterministic B&B's leaves (which are near-duplicates — the global optimum ± one nudged word). Then greedily pick the best, then each next-best whose **placement-distance ≥ τ** from all already-picked (distance = fraction of words placed at a different cell/orientation; τ default ~0.3, tunable/calibrated). This is why candidates ride the greedy path (plan Phase 6).

### 7.5 Acceptance criteria
**AC-ARR-1** `arrange --strict` resolves to exactly one of three outcomes and never drops silently: **(a)** all words placed in a legal layout; **(b)** search *completed* without placing a word ⇒ exit non-zero **naming the unplaceable word(s)**; **(c)** budget exhausted before a complete placement ⇒ exit non-zero reporting **"not proven within budget"** (NOT "unplaceable").
**AC-ARR-2** `arrange --best-effort` places a maximal subset, emits the layout, and reports the dropped set in the payload's `diagnostics.arrange.dropped` (json-output-spec §6.4; echoed on stderr under `--verbose`); never fails when ≥1 word is placeable.
**AC-ARR-3** Every emitted layout re-validates as legal against the shared legality core (free property test).
**AC-ARR-4** For the benchmark fixtures, the optimizer's reward ≥ the old `solve`'s arbitrary-fill reward (no regression in interlock quality).
**AC-ARR-5** `--max-size N` output is the tight enclosing square (side `max(H,W)` ≤ N), anchored at the content's top-left; `--size N` output is exactly N×N with blocks for empty cells. The two flags are mutually exclusive; a bare `arrange` defaults to `--size 15`.
**AC-ARR-6** Output is deterministic by default (INV-2): with no `--seed`/`--shuffle`, repeated runs are byte-identical. The opt-in perturbation (§7.6) is the sole randomization: `--seed N` is reproducible; `--shuffle` varies per run and prints its seed under `--verbose`.
**AC-ARR-7** `--candidates K` returns up to K layouts, each pairwise ≥ τ apart; fewer than K only if fewer ≥τ-distinct layouts exist (and that is reported).
**AC-ARR-8** `--enumerate` counts/emits all feasible full placements, matching the old `--all` count on shared fixtures.
**AC-ARR-9** The carried incremental reward equals a from-scratch `layout_reward` recompute on the final placement (delta correctness).
**AC-ARR-10** A run that hits the node/inference budget returns the best layout found so far and **flags** that the optimum was not proven.

*(Calibration **locked by DP-1 (§10)** at the v1 defaults: `WCap:WTail = 5:1` (ε = 0.2), `target = ceil(L/2)`, `τ = 0.30`. The Phase-1.5 gate showed the cap is largely inert on realistic inputs (objective ≈ total-crossings), so fine-tuning has low leverage; `--check-target` remains the tunable escape hatch where `ceil(L/2)` is unreachable.)*

### 7.6 Search perturbation — `--seed` / `--shuffle` (opt-in, off the deterministic path)

The default `arrange` search is fully deterministic (INV-2, AC-ARR-6). Two opt-in flags perturb **only** the strict single-layout search, for pseudo-random variety without abandoning reproducibility:

- **`--seed N`** (N ≥ 0) seeds the RNG to `N` and shuffles branch order at two seams — the seed-word choice, and within each equal-MRV-count bucket the candidate order. Same input + same `N` ⇒ byte-identical output. Distinct seeds give distinct layouts, all still placing every word: the fail-first MRV heuristic and reward selection are preserved, only *ties* are reordered.
- **`--shuffle`** draws a fresh seed each run (a different layout every time) and, under `--verbose`, prints it so a liked layout reproduces with `--seed N`. Mutually exclusive with `--seed`.

The seed never touches the deterministic flow: with neither flag, no RNG is seeded or consulted (INV-2 holds by construction). The perturbation rides the strict path only — it is rejected in combination with `--best-effort`, `--candidates`, and `--enumerate` (none reach the seam, so silently doing nothing would mislead). When set, the seed is recorded at `diagnostics.arrange.seed` (json-output-spec §6.4).

---

## 8. Flavour B — cryptic setting tools

### 8.1 `lint` — grid validator / profiles  **[LOCKED]**
Consumes a canonical layout (what `arrange` emits) and reports **PASS / WARN / FAIL per rule, per word**, plus a summary verdict. Reuses the shared metric predicates; needs no engine change.

Rules (per the research doc + symmetry/size addenda):
| Rule | Default severity | Notes |
|---|---|---|
| Min length ≥ 3 | profile-dependent | FAIL for cryptic profiles; WARN for the TOC default. Always *reported*. |
| Checked fraction ≥ `ceil(L/2)` per word | FAIL (cryptic) | Reuses `word_meets_half/2`. |
| Longest unchecked run | FAIL if ≥ 3 (triple-unch ban) | `max_unch_run`. |
| Times positional double-unch (no double-unch at word ends) | WARN/FAIL by profile | Optional floor predicate. |
| Odd-vs-even checking | WARN | Advisory. |
| Connectivity | PASS by construction | Already enforced. |
| **Symmetry deficit** | **profile-set, never unconditional FAIL** | 180° rotational symmetry is a **strong default with theme-gated exceptions**, NOT a universal rule. Report a deficit; severity by profile; honour a per-puzzle `--allow-asymmetry` override (mirrors Exet's toggle). |

Profiles (`--profile`): `blocked-uk` (enforce symmetry, ≥half checking, triple-unch ban), `barred-ximenean` (per-length unch band; relaxed symmetry), `american` (every cell checked), `toc` (advisory-only — the relaxed Flavour-A default). Inherently-judgemental checks (surface, loose synonyms) are **advisory**, never hard FAIL (Principle #7: Ximenean-configurable, never Ximenean-only).

**AC-LINT-1** Given a layout + profile, `lint` emits a per-rule per-word PASS/WARN/FAIL report and a summary verdict.
**AC-LINT-2** Exit code is 0 for PASS/WARN, non-zero only when a FAIL-severity rule trips under the chosen profile.
**AC-LINT-3** Symmetry never hard-FAILs under `--allow-asymmetry`, and never FAILs at all under a profile that doesn't enforce it.
**AC-LINT-4** Each profile applies exactly its documented rule/severity set; advisory checks never escalate to FAIL.

*(OD-7 resolved (DP-3): the barred-ximenean profile ships the **primary-sourced** Ximenean per-length unchecked-letter band — none in a 3, one in 4–5, two in 6–7, three in 8 (Azed self-limits to 2), ⌊L/3⌋ in 9+ — from Ximenes' 1966 *Ximenes on the Art of the Crossword* + the Azed slip conventions (`lint.pl` `barred_max_unch/2`). Symmetry is **relaxed to advisory** (WARN); exact per-publication barred symmetry codes are not modelled — a documented v1 simplification, not a barred-grid bar/edge model (§3).)*

### 8.2 `export` — standard-format interchange  **[LOCKED]**
Transformations of the canonical JSON (§6.5), not new emitters.
- **ipuz** (`--to ipuz`): invert `null`→`#`, split the merged cell object into parallel `puzzle`/`solution` arrays, group words by direction into a `clues` dict, add `version`/`kind`/`dimensions`. Emit **ipuz v2** (`"version":"http://ipuz.org/v2"`, `"kind":["http://ipuz.org/crossword#1"]`) plus a default `title` (optional in ipuz, so kotwords is unaffected; but Exet imports ipuz→Exolve and crashes Save on a null title, so a non-null one keeps it Exet-round-trippable — mirrors the Exolve fix, V1). ipuz carries no symmetry field — symmetry stays a crosswordsmith-side lint concept.
- **Exolve** (`--to exolve`): plain-text, git-diffable; the only common format carrying bars, ninas, definition-span marking, annotations. Line-per-element granularity (round-trips to Exet). Emits a default `exolve-title` — Exet's Save path does `title.replace(...)` and crashes on a null/absent title, so a non-null one is required to round-trip (found + fixed during the 2026-07-01 AC-EXP-2 verification; see [`exet-verification.md`](./exet-verification.md)).
- **`.puz`/`.jpz`/PDF**: reached via off-the-shelf `kotwords` from the ipuz output — **not** emitted natively (§3).

**AC-EXP-1** `export --to ipuz` produces spec-valid ipuz v2 that a third-party consumer (e.g. via kotwords) ingests without error.
**AC-EXP-2** `export --to exolve` produces text that round-trips through Exet (load → save → equivalent grid+entries). *(Verification: a byte-exact golden (`tests/golden/export_bundled_17.exolve`) pins the emitter, and Exet ingestion is a manual check — there is no in-repo Exolve parser, so third-party round-trip cannot be automated in plunit. Run the step-by-step round-trip via [`exet-verification.md`](./exet-verification.md); see revamp-audit R13.)*
**AC-EXP-3** Export preserves enumerations and `meta`-borne clue text where the target format has a field for them; no data invented.

### 8.3 Stock-grid library + house-style profiles  **[LOCKED]**
A **bundled, curated** set of pre-validated legal grid templates (black-square patterns / slot lists), **not** a generator — real publications curate only dozens (Times 64, Guardian ~72), so this is a small static asset. Realistic presets for the *current block-only* world: **15×15** (primary), **13×13** (quick), **23×23** (jumbo) — all blocked, all odd. The `~12×12 barred` cluster and rectangular/thematic shapes are out of scope (need the barred engine / are different products).
- Templates double as `lint` profile inputs and (later) `fill` inputs.
- Validated once at design time with the shared metric predicates repurposed as **template validators**.

**Template schema (OD-5, resolved by DP-1): a black-square mask is the single source of truth** — required fields `{name, size, mask:["#.....", ...]}` (one string per row; `#` = block, any other char = light). An optional `symmetry` annotation MAY appear (the shipped grids carry `"rot180"`) but is **not** trusted or required — the validator re-derives symmetry from the mask, so a declared value is descriptive only (revamp-audit R10). Slots are *derived* on load by the existing run-scanning / clue-numbering, not stored — no redundant slot list. Each shipped grid is **validated at design time by `lint --profile blocked-uk`** (the metric predicates as template validators).

**OD-6 (resolved): the v1 library ships three original, `lint`-validated, license-clean grids** — `blocked_13a`, `blocked_13b` (13×13), `blocked_15a` (15×15), all 180°-rotational, all PASS under `lint --profile blocked-uk`. A small starter set, as intended (real publications curate only dozens); it grows by adding more validated masks under `grids/`. Loaded + validated by `stockgrid.pl` (the mask → lights derivation + the blocked-uk validator); the set's legality is a CI regression (`tests/stockgrid.plt`).

### 8.4 `fill` — grid-first, open-dictionary auto-fill  **[LOCKED]**
The path to *authentic* cryptic layout: take a pre-validated legal blocked grid (a §8.3 stock template or a user grid), fill its slots from an open dictionary subject to crossing constraints, with the user's words pinned as **seeds** (via the fragment primitive). It is a **separate engine** from `arrange` — it shares the substrate/metric/emit layers, not the legality core (black squares delimit runs explicitly).

**Architecture (decided — DP-1/DP-2 §10):**
1. Grid template = a §8.3 black-square mask; slots derived by run-scanning (reuses `stockgrid.pl`).
2. Dictionary = an **in-memory pattern index** (by length + a (length, position, letter) index), normalised to A–Z. Default lexicon **UKACD18** (redistributable freeware — J Ross Beresford; ship its copyright + license text verbatim when bundled; **not** BSD-3, contrary to earlier drafts — see the license-pass note in [`revamp-audit-findings.md`](./revamp-audit-findings.md)); `--dict FILE` overrides; a small fixture wordlist ships for tests until the full data is bundled.
3. **MRV/backtracking** search selecting *dictionary* words per slot (not a fixed list); deterministic (dictionary order + lowest-start tiebreak; INV-2), bounded by a node/inference budget.
4. **Seeds** are HARD PINS via the §6.6 fragment-grid primitive (OD-3): pin-and-fill-around.
5. **No-fill contract** (OD-4, INV-3): report the unfillable slot(s) and exit non-zero — `arrange` strict's outcomes (filled / infeasible-naming-slots / not-proven-within-budget).
6. Output is the canonical layout (so it composes with `lint`/`export`).

**Bifurcation:** `fill`-blocked is built; `fill`-barred (different cell/edge model) is out of scope (§3) and may never ship.

**AC-FILL-1** Given a legal grid + dictionary, `fill` either emits a canonical layout whose every slot is a dictionary word with consistent crossings, or exits non-zero naming the unfillable slot(s) / reporting "not proven within budget".
**AC-FILL-2** Seed words pinned via `--seeds` appear at exactly their slots in any produced fill (reuses AC-FRAG-3).
**AC-FILL-3** Output is deterministic (INV-2): identical grid + dictionary + seeds ⇒ byte-identical layout.
**AC-FILL-4** A produced fill re-validates as a legal layout against the shared metric predicates (the grid's `lint` verdict is unchanged by the fill — letters do not affect structure).

### 8.4a `fill --min-score` — score-aware candidate selection  **[LOCKED — DP-4]**
Scored fill closes the #1 *measured* competitive gap (research §D): scoreless MRV places non-word junk (`AAAAA`/`AAAAQ`, [`benchmarks/fill_quality/`](../benchmarks/fill_quality/README.md)), while every serious filler (Crossword Compiler / CrossFire / ingrid_core) drives search by per-word score. This is an **extension of §8.4, not a new engine** — the same MRV/backtracking substrate with one added candidate filter and one added sort key. It stays inside the identity: score is intrinsic to the *dictionary* (word→score), **never** clue metadata (INV-1 untouched — the solver reads no `meta`), and never appears in the output layout.

**Scored-dictionary format (DP-4):** `fill` accepts a `word;score` list (one entry per line, `score` an integer 0–100 — the Spread-the-Wordlist / Broda convention). Plain wordlists (no `;score`) read as before and are treated as a **uniform score** (every word equal), so an unscored dict is exactly backward-compatible. Malformed lines are reported, not silently dropped (INV-3). **License (INV-4):** the bundled default stays UKACD18 (permissive freeware, unscored) — the core distribution ships **no** NC/SA data. Scored fill is opt-in via `--dict FILE` pointing at a user-supplied scored list; Spread the Wordlist (CC BY-NC-SA 4.0) is the recommended license-clean scored source, and if ever bundled it ships as a clearly-flagged **optional** asset with its attribution + ShareAlike notice verbatim (never the default). A small **scored** fixture wordlist ships for tests.

**`--min-score N` (DP-4):** a **hard candidate prune** — every candidate word scoring `< N` is removed from each slot's domain *before* search (identical to the benchmark's `score≥50`-dict column, which is the regression target). Default `N = 0` (no prune): scored fill never silently changes feasibility — a grid fillable today stays fillable. A prune that empties a slot's domain is an ordinary no-fill outcome (AC-FILL-1: name the unfillable slot, exit non-zero, INV-3).

**Ordering (DP-4, determinism-preserving):** within each slot's surviving domain, candidates are ordered **score-descending, then the existing §8.4 tiebreak** (dictionary order, lowest-start). Score is a *total* order key layered above the current one, so: (a) the first legal fill found prefers high-score words → better quality even at `--min-score 0`; (b) for an unscored (uniform-score) dict the sort collapses to exactly the §8.4 order, so AC-FILL-3 still holds byte-for-byte. INV-2 preserved.

**Fill-quality report (DP-4, INV-3):** `fill` reports the produced fill's per-slot score distribution — `n` slots, `mean`, `min`, and a below-threshold count — mirroring `lint`'s report model and the benchmark's `score_fill.py`. The report is a **sidecar** (human summary on stderr; `--report-json` for machines); the canonical layout on stdout is unchanged, so export/round-trip and AC-FILL-4 are untouched (letters and scores never affect structure).

**AC-FILL-5** `fill --min-score N` excludes every candidate scoring `< N` from each slot's domain; a produced fill contains no word scoring `< N`, matching the benchmark `>=N`-dict result ([`benchmarks/fill_quality/`](../benchmarks/fill_quality/README.md)).
**AC-FILL-6** Determinism (INV-2) holds with scoring: score-descending is a total order over the §8.4 tiebreak, and for an unscored/uniform-score dict every produced fill is byte-identical to the pre-scoring engine (AC-FILL-3 unchanged).
**AC-FILL-7** `fill` emits a fill-quality report (n / mean / min / below-threshold count) per INV-3 without altering the canonical layout on stdout.
**AC-FILL-8** A `word;score` dict is ingested with scores 0–100; unscored lines default to the uniform score; malformed lines are reported (INV-3); the bundled default dictionary stays permissively licensed (INV-4).

*Implementation-time refinements (surfaced by [`research/wordlist-scoring-2026.md`](research/wordlist-scoring-2026.md), to resolve when built — a small DP-5 if either changes contract):* (i) **Scale units** — real scored lists are **not** uniformly 0–100 (STW/Broda 0–100, XWord Info 5–60, Crossword Compiler ≤50); v1 treats `--min-score N` as being in **the dict's own native units** (documented, deterministic, no auto-normalisation), not a fixed 0–100 assumption. (ii) **Score-0 semantics** — some lists (STW) reserve **score 0 for a deliberate blocklist** (harmful/X-rated), so the default prune should exclude score-0 entries (`score ≥ 1`) rather than literally "include everything," keeping the spirit of "default never *removes usable* words." Both are documented + evidence-cited; neither reopens the LOCKED core (DP-4).

### 8.5 Backlog — tagged, not yet specified
These are recognised future capabilities. Each is **out of scope until specified here** (a §10 decision pass per feature). Listed so they are tracked, not so they are built. **Evidence & prioritization** for these against the 2026 competitive field live in [`research/setter-tool-landscape-2026.md`](research/setter-tool-landscape-2026.md); the ★-tagged rows below are the three the 2026 research flags as the current highest-leverage backlog and which need a decision pass + spec + implementation plan next.

| Feature | Flavour | Note |
|---|---|---|
| ★ **Scored fill** (`fill --min-score` + fill-quality report) | B | **2026 research: the #1 competitive gap.** Every serious filler (Crossword Compiler / CrossFire / ingrid_core) drives fill by per-word score; crosswordsmith's scoreless MRV places non-word junk (measured — `AAAAA`/`AAAAQ`). Closable **license-clean** via Spread the Wordlist (CC BY-NC-SA 4.0). A `score≥50` dict prefilter already recovers ingrid-parity quality in the prototype. **Spec'd — DP-4 → §8.4a; ready to build (not yet implemented).** Evidence: research §D + [`../benchmarks/fill_quality/`](../benchmarks/fill_quality/README.md). |
| Per-entry cluing-potential annotation (`meta.cluing`) | B | Computed at emit-time answer→meta join; engine stays metadata-agnostic. 2026 research: MyCrossword's anagram-%/device-balance heuristic is a ready blueprint. |
| Nina seeding | Shared | **Absorbed into the fragment-grid primitive** (§6.6) — not a separate feature. (Qxw "free lights" confirm the capability is served nowhere open.) |
| `min_len:3` hard floor option | B | Soft default for TOC; always *report* sub-3. |
| ★ Clue stockpile keyed by answer (`meta`) | B | Pure data plumbing on the answer key. **2026 research: unserved across the entire surveyed field (both passes) — a genuine differentiator, and low-effort.** Needs decision pass + spec + plan. |
| Prolog-native pattern/anagram query endpoint | B | Maps onto unification/backtracking. 2026 research: Nutrimatic's `A/C/V/#/_` + `<…>` + `&` syntax is a ready model. |
| Wordplay-decomposition engine (`meta.clueCandidates`) | B | Generate + verify via backtracking; prune hard. |
| Batch/deck interface + shell exit codes | Shared | Pipeline embedding. |
| Scored `arrange` drop-order/tiebreak | A (legacy) | Secondary to the scored-fill row above; scores were de-scoped for `arrange` (§3), so this is a tiebreak/drop-order signal only (must preserve determinism + the closed-set model). Currently dormant. |

> The third ★ priority — **`stats`/`inspect` + `diff` verbs on the `xword` companion** — is tracked in [`xword-status.md`](xword-status.md) (xword lives on its own spec/status, not §8.x); the full xword breadth shortlist is [`plans/xword-breadth-expansion.md`](plans/xword-breadth-expansion.md).

---

## 9. Cross-cutting requirements

Re-stated as testable obligations (the invariants from §4 plus determinism/license):
**AC-X-1** No component threads `meta` through placement or numbering (INV-1) — verifiable by a placement run being identical with `meta` stripped.
**AC-X-2** All output is stable sorted-key JSON; no nondeterministic ordering (INV-2).
**AC-X-3** Every drop/compromise/infeasibility is reported in output, never swallowed (INV-3).
**AC-X-4** No bundled data carries a non-permissive license (INV-4).

---

## 10. Open-decisions register

The **only** sanctioned places where scope is still undecided. A component cannot be promoted from DEFERRED/PARTIAL to LOCKED until its rows here are resolved in a recorded decision pass.

| # | Component | Open decision | Status |
|---|---|---|---|
| OD-1 | `fill` (§8.4) | Blocked-only v1, or design barred-compatibility in from the start? | **resolved (DP-1): blocked-only v1** — reuse the existing blocked cell model; barred is a separate engine, deferred indefinitely (§3). |
| OD-2 | `fill` | Dictionary integration shape: in-memory index, external index, query protocol? Which lexicon as default (UKACD18 confirmed license-clean)? | **resolved (DP-2): in-memory pattern index** (by length + a (length,position,letter) index), normalised to A–Z; **default lexicon UKACD18** (redistributable freeware — **not** BSD-3; ship its license text verbatim when bundled) via `--dict FILE` (bundled when the data is on hand; a small fixture wordlist ships for tests). |
| OD-3 | `fill` | How seeds/fragment-grid semantics carry into open-dictionary fill (pin-and-fill-around vs. seed-as-hint). | **resolved (DP-1): pin-and-fill-around** — seeds are HARD PINS via the §6.6 fragment-grid primitive, identical to `arrange`. |
| OD-4 | `fill` | Which house-style profiles ship in v1, and the failure contract when no fill exists. | **resolved (DP-2): the §8.3 stock-grid templates ARE the v1 profiles** — `fill` fills any pre-validated blocked grid (stock or user-supplied) from the dictionary, seeds pinned. **No-fill contract:** report the unfillable slot(s) and exit non-zero (INV-3), mirroring `arrange` strict's outcomes (filled / infeasible-naming-slots / not-proven-within-budget). |
| OD-5 | Stock-grid library (§8.3) | Template schema (black-square mask vs. explicit slot list vs. both). | **resolved (DP-1): black-square mask** is the single source of truth; slots are derived on load (no redundant slot list). |
| OD-6 | Stock-grid library | Which specific grids seed the bundled library, and provenance/license of each. | **resolved (§8.3 build): ships `blocked_13a`, `blocked_13b`, `blocked_15a`** — original, 180°-symmetric, all PASS under `lint --profile blocked-uk`. Set grows by adding validated masks. |
| OD-7 | `lint` barred profile (§8.1) | Exact barred-Ximenean per-length unch table; per-publication barred symmetry codes — primary-source before building. | **resolved (DP-3): the Ximenean band is primary-sourced + built** (Ximenes 1966 + Azed); symmetry relaxed to advisory (per-publication codes not modelled — documented simplification). |
| OD-8 | Backlog (§8.5) | Each backlog feature needs its own decision pass + spec section before implementation. | open |
| OD-10 | `fill` scored extension (§8.4a) | Scored-dict format + default lexicon; `--min-score` as hard prune vs soft tiebreak; determinism under scoring; fill-quality report shape; license of any bundled scored list. | **resolved (DP-4): see §8.4a** — `word;score` (0–100) dict; `--min-score N` hard prune (default 0, so feasibility never silently changes) + score-descending ordering layered on the §8.4 tiebreak (INV-2 held; unscored dict byte-identical to the pre-scoring engine); fill-quality sidecar report (INV-3); bundled default stays permissive UKACD18, STW (CC BY-NC-SA) opt-in via `--dict`, never bundled as default (INV-4). |
| OD-9 | `arrange` (impl, not product) | Empirical calibration of `ε`/`target`/`τ`; thin fragment-form syntax; duplicate-answer disambiguation. | **resolved (DP-1):** calibration locked at `WCap:WTail=5:1` (ε=0.2) / `target=ceil(L/2)` / `τ=0.30`, with `--check-target` the tunable escape hatch; thin fragment-form **deferred** (canonical-only) *(since implemented for `arrange --fragment` and `fill --seeds` — AC-FRAG-4, see §6.6)*; duplicate input answers **rejected** (answers are unique). |

### Decision passes

- **DP-4 (2026-07-15).** Resolved OD-10, promoting the ★ scored-fill backlog item (§8.5) through its required decision pass into the LOCKED **§8.4a** extension. Rationale: the 2026 landscape research (§D) + the [`benchmarks/fill_quality/`](../benchmarks/fill_quality/README.md) prototype *measured* scoreless MRV placing non-word junk while a `score≥50` prefilter recovered ingrid_core-parity quality — making scored fill the #1 competitive gap and its fix a small, no-regret extension of the existing MRV engine (one filter + one sort key, not a new engine). Score is treated as an intrinsic dictionary property (word→score), keeping the solver metadata-agnostic (INV-1) and the output layout score-free. `--min-score` defaults to 0 so feasibility is never silently changed; score-descending ordering is layered *above* the §8.4 tiebreak so determinism (INV-2) holds and an unscored dict stays byte-identical to the pre-scoring engine. The bundled default lexicon stays the permissively-licensed UKACD18 (INV-4) — Spread the Wordlist (CC BY-NC-SA) is an opt-in `--dict`, never bundled as default — and the `benchmarks/fill_quality/` harness is the standing regression (native `--min-score 50` must match its `>=50`-dict column). **Not yet implemented** — spec'd and ready to build. *Post-DP-4 evidence:* [`research/wordlist-scoring-2026.md`](research/wordlist-scoring-2026.md) (2026-07-15) confirmed the "50 = clean floor" convention and the licensing verdicts (no cleanly-bundleable scored list yet — the default-stays-unscored choice is validated), and surfaced two implementation-time refinements now recorded in §8.4a (scale-units + score-0 blocklist).
- **DP-1 (2026-06-30).** Resolved OD-1, OD-3, OD-5, OD-9 (see rows above). Rationale: these are no-regret — OD-9 was already settled by the shipped `arrange` (Phases 1–7); OD-1/OD-3 reuse existing primitives (the blocked cell model, the fragment-grid pins) so `fill`, *when built*, inherits `arrange`'s seed semantics; OD-5 picks the smallest asset (a mask; slots are derivable). **`fill` (§8.4) stays DEFERRED** — promotion to LOCKED still needs OD-2 (dictionary shape/lexicon) and OD-4 (profiles + no-fill contract). **Stock-grid (§8.3)** moves to buildable: OD-5 fixed here, OD-6 closes as its grids are authored.
- **DP-3 (2026-06-30).** Resolved OD-7. The Ximenean per-length unchecked-letter band was **primary-sourced** (Ximenes' 1966 *Ximenes on the Art of the Crossword* + the Azed slip conventions, two corroborating sources) and built as the `barred-ximenean` `lint` profile (`checked_band` rule over `barred_max_unch/2`). Barred symmetry is **relaxed to advisory** (WARN); exact per-publication symmetry codes are *not* modelled and barred linting applies the checking math to the canonical layout's words — a dedicated bar/edge model stays out of scope (§3). Documented simplification, recorded so the deviation is explicit.
- **DP-2 (2026-06-30).** Resolved OD-2, OD-4 (see rows above), completing the OD-1…4 set, so **`fill` (§8.4) is promoted DEFERRED → LOCKED** and built. Rationale: a single-puzzle CLI fill wants the simplest self-contained machinery — an in-memory pattern index over a bundled, license-clean lexicon (UKACD18, `--dict`-overridable) — and it should reuse what already exists: the stock-grid templates as the grid/profile, the fragment-grid primitive as hard pins (OD-3), and `arrange` strict's report-and-fail outcome contract (INV-3). The full UKACD18 is added when the data is on hand; a small fixture wordlist ships so the engine is tested end-to-end now.

---

## 11. Glossary

- **Flavour A / Flavour B** — the two product surfaces under one CLI (§1).
- **Closed-set / closed-set inversion** — placing only the user-supplied words, no filler dictionary (Flavour A).
- **Grid-first** — designing a legal black-square pattern up front, then filling it from a dictionary (Flavour B `fill`).
- **Checking / checked cell** — a cell shared by an across and a down word (occupancy ≥ 2).
- **Unch / unchecked run** — consecutive cells belonging to only one word; runs ≥ 3 are the triple-unch FAIL.
- **Fragment grid** — a partial layout the engine solves from; presence = fixed (§6.6).
- **Capped reward / cap target** — `min(checked, ceil(L/2))`; the cap is what produces balance (§7.2).
- **Anytime / budget** — return best-so-far within a node/inference-count budget (not wall-clock, for determinism) rather than proving the optimum (§7.3).
- **Blocked vs. barred** — black-square grid vs. bar-delimited grid; barred is a separate engine (out of scope for now).
- **Stock grid** — a pre-validated, bundled legal grid template (§8.3).
- **LOCKED / PARTIAL / DEFERRED** — implementation-readiness tags (see "How to read this document").
