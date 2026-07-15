# Plan: fill-scoring uplift — ranked proposals

**Status: PROPOSAL (analysis only, 2026-07-15).** Nothing in this document is
implemented, and nothing here changes the spec by itself — where a proposal
needs a decision pass it says so explicitly (spec
[§10 change discipline](../design-spec.md)). This is the ranked, evidence-cited
menu for uplifting crosswordsmith along the **fill scoring** axis, built on the
adversarially-verified 2026 research; it does not re-derive any claim from the
web.

**Evidence base (source of truth for every claim below):**

- [`research/wordlist-scoring-2026.md`](../research/wordlist-scoring-2026.md) —
  band definitions, scoring signals, per-list license verdicts, §5 spec impact.
- [`research/setter-tool-landscape-2026.md`](../research/setter-tool-landscape-2026.md)
  §D (the measured scored-fill gap) + §F (ranked roadmap). §D's license claims
  were patched by the scoring doc.
- [`design-spec.md`](../design-spec.md) — §8.4 (fill, LOCKED), §8.4a (scored
  fill, LOCKED via DP-4, incl. two unresolved implementation-time refinements),
  §8.5 backlog, §10 register, INV-1..4, §7 (capped reward, anytime/budget model).
- [`benchmarks/fill_quality/`](../../benchmarks/fill_quality/README.md) — the
  standing fill-quality regression (measured 2026-07, STW snapshot 2026-07-01).
- [`README.md`](../../README.md) "Implementation notes / limitations" — the
  "fill is unscored, and its search has a budget ceiling on hard grids" bullet.
- [`fill-perf-campaign.md`](./fill-perf-campaign.md) — executed 2026-07-05;
  records the determinism charter's constraints on fill-search changes.
- [`STATUS.md`](../STATUS.md) — research-derived backlog priorities table.

---

## 1. The two axes — quality vs. search power (read this first)

The benchmark separated two **orthogonal** competitive gaps, and every proposal
below is tagged with which one it addresses:

- **QUALITY** — scoreless MRV places *legal junk*: on grids it completes, mean
  STW score 27–42 with entries at 0 (`AAAAA`, `AAAAQ`, `ABAAB`, `AEAEA`), where
  ingrid_core at `--min-score 50` holds mean/min 50
  ([benchmark results table](../../benchmarks/fill_quality/README.md)). A crude
  `score>=50` dict prefilter **fully recovers ingrid-parity quality** on every
  completable benchmark grid — quality is the cheap axis.
- **SEARCH POWER** — crosswordsmith **cannot complete `blocked_13a`** (standard
  13×13, full-length slots) within its 800M-inference budget *with any dict*,
  while ingrid completes it in ~10s at `--min-score ≤30`
  ([benchmark, completion table](../../benchmarks/fill_quality/README.md)).
  This is the harder engineering gap, and scoring does not touch it.

**The interaction cuts the wrong way: a `--min-score` prune makes a hard grid
*less* fillable, not more.** Measured: ingrid completes `blocked_13a` at
min-score 20/30 in ~10s but **times out at 50**. Pruning shrinks slot domains,
which trims per-node branching but removes completions. So the quality uplift
(FS-1/FS-2) *raises* the value of the search-power work (FS-4): once users can
ask for clean fills, hard grids will hit "no fill within budget" *more* often
at high thresholds, and that must be reported as the honest outcome it is
(INV-3), never advertised as a speed feature.

| Proposal | Quality | Search power |
|---|---|---|
| FS-1 core `--min-score` build | ✅ primary | ❌ (ordering may shift trajectories; never lifts the ceiling) |
| FS-2 refinements (i)+(ii) via DP-5 | ✅ (correct defaults) | ⚠️ guards against accidental feasibility loss |
| FS-3 benchmark harness growth | 📏 measures both | 📏 measures both |
| FS-4 fill search-power pass | ❌ | ✅ primary |
| FS-5 xword fill-quality report | ✅ (visibility) | ❌ |
| FS-6 dict tooling + license errands | ✅ (enabler) | ❌ |
| FS-7 scored `arrange` signal | ✅ (marginal) | ❌ |

## 2. Guardrails held throughout

Every proposal was checked against, and none breaks: **INV-1** (score is an
intrinsic dictionary property, word→score; never `meta`, never in the output
layout), **INV-2** (total-order sort keys only; no wall-clock, no unseeded RNG),
**INV-3** (every prune/drop/no-fill reported), **INV-4** (no NC/SA or
personal-use data bundled — the research confirmed **no cleanly-bundleable
scored list yet exists**), the **closed-set inversion** (scoring = prune +
ordering in `fill`; tiebreak/drop-order only in `arrange`; never word
invention), and the CLI-first / not-a-clue-generator identity. Directions that
would break one are listed in §4 as explicit cuts, not proposals.

---

## 3. Ranked proposals

### FS-1 — Build §8.4a exactly as locked: `word;score` ingestion, `--min-score`, score-descending ordering, fill-quality report

**What / surface (`fill`).** Implement the LOCKED §8.4a/DP-4 contract, verbatim:
`word;score` dict ingestion (plain lists read as a uniform score, so unscored
dicts stay byte-identical — AC-FILL-6); `--min-score N` as a **hard candidate
prune** before search; candidate ordering **score-descending, then the existing
§8.4 tiebreak** (a total order — INV-2); and the **fill-quality sidecar report**
(n / mean / min / below-threshold; human summary on stderr, `--report-json` for
machines) with the canonical stdout layout untouched (AC-FILL-7). Malformed
lines reported, never dropped (INV-3, AC-FILL-8). A small **original,
license-clean scored fixture** ships for tests (INV-4 — no STW bytes in the
repo, ever).

**Why now.** The #1 measured competitive gap (landscape §D, verified 3-0):
every serious filler scores (Crossword Compiler PGF, CrossFire, ingrid_core),
and the benchmark proved scoreless MRV fills junk while a `>=50` prefilter
recovers full parity. Both §8.4a mechanisms are established, documented filler
behaviour — Crossword Compiler ships exactly this pair (relative-rank ordering
*plus* a "Minimum word score" hard floor;
[scoring doc §1](../research/wordlist-scoring-2026.md)). The decision pass is
already done (DP-4 → §8.4a LOCKED); STATUS lists it as "spec'd; not started —
needs implementation plan + build". This document is the proposal ranking; the
implementation plan is the next artifact.

**Tags.** Effort **medium** · leverage **very high** · identity-fit **direct**
(landscape §F ranks it #1).

**Spec impact.** None — fits §8.4a as-is (ACs FILL-5..8). The only contract
question it *touches* is the default-prune value, which FS-2 resolves; build
FS-1 and FS-2 in the same milestone so the default never ships twice.

**Measurement (benchmarks/fill_quality/).** The harness's README already names
the acceptance test: add the `cs_minscore_<g>.json` variant (native
`--min-score 50`) and it **must match the `>=50`-dict column** — mean/min 50,
zero below-clean, on open4/open5/mini7/mini9 (AC-FILL-5). The `--report-json`
numbers must agree with `score_fill.py`'s post-hoc stats on the same fill.
Unscored-dict byte-identity is checked by the existing goldens + `make fuzz`.

**Risks / traps.**
- *Determinism:* equal scores are common (STW has 120k words at exactly 50) —
  the score key must sit strictly *above* the full existing tiebreak so ties
  collapse to today's deterministic order.
- *Scale trap (defused by docs, not code):* `--min-score 50` against an
  XWord Info-scale (5–60) list is near-maximum strictness — see FS-2(i).
- *Uniform-score value:* the spec leaves the unscored-line uniform score
  unspecified; it must land **≥ 1** or FS-2(ii)'s default prune would empty
  unscored dicts. Pin it in the DP-5 note (FS-2).
- *Degenerate flag use:* `--min-score N` above a uniform-score dict's score
  empties every domain — an ordinary all-slots-unfillable outcome per
  AC-FILL-1, but worth a targeted hint in the error text (INV-3 spirit).
- *License:* recommend STW via `--dict` in docs (point-at is fine; CC BY-NC-SA
  forbids bundling only). The test fixture must be authored, not derived from
  any scored list.

---

### FS-2 — Resolve refinements (i) scale-units and (ii) score-0 default, as a recorded micro-pass (DP-5)

**What / surface (`fill` spec + docs).** Adopt the scoring research's two
recommendations and record them: **(i)** `--min-score N` is in **the
dictionary's own native units** — documented, deterministic, no
auto-normalisation and no `--score-scale` flag in v1
([scoring §5.2, option (a)](../research/wordlist-scoring-2026.md));
**(ii)** the default prune becomes **`score ≥ 1`** — excluding only the score-0
blocklist floor — instead of the literal "include everything" default 0
([scoring §5.3, option (b)](../research/wordlist-scoring-2026.md)). Alongside
(ii), pin the two adjacent details: unscored lines' uniform score is ≥ 1
("unrated" ≠ "blocklisted"), and docs *recommend* `--min-score 50` as the
practical clean floor (distinct from the code default).

**Why now.** Both are flagged in §8.4a itself as "to resolve when built — a
small DP-5 if either changes contract". (ii) **does** change LOCKED contract
text (DP-4 says "Default `N = 0` (no prune)"), so it is not an implementation
detail — it needs the micro-pass. The evidence is direct: scales are
heterogeneous (STW/Broda 0–100, XWord Info 5–60, CC ≤50, christophsjones
~1–50 — scoring §1), and STW reserves 0 for a deliberate harmful/X-rated
blocklist, so a default that includes score-0 admits intentionally-excluded
entries. The reframed invariant is the research's own wording: the default
"never removes *usable* words" (rather than "never changes feasibility").
The `>=50` recommendation is documented convention, not folklore — three
independent primary sources (STW FAQ, XWord Info FAQ, ingrid's default).

**Tags.** Effort **low** (spec edit + a default constant + docs + 2–3 tests) ·
leverage **medium-high** (it's the difference between correct and subtly-wrong
defaults) · identity-fit **direct**.

**Spec impact.** A recorded **DP-5** entry in §10 + one-line amendments to
§8.4a (the default-prune sentence and a scale-units sentence) and AC-FILL-8's
"scores 0–100" phrasing (→ "integer scores in the dict's native units").
OD-10 stays resolved; DP-4's core (prune + ordering + report + license
posture) is untouched.

**Measurement.** Extend the harness's scored fixture with a score-0 entry and
assert it never appears in a default-flags fill; add a feasibility check that
the ≥1 default changes no benchmark grid's completability (on STW-derived
dicts, only the blocklist floor is removed). For (i): a tiny 5–60-scale fixture
asserting no normalisation occurs (byte-identical fills whether the same ranks
are expressed as 5–60 or 0–100? — no: assert the *threshold semantics*, i.e.
`--min-score 50` against the 5–60 fixture prunes to the ≥50 tail exactly).

**Risks / traps.** (ii) technically makes a grid fillable *only* by
blocklisted words fail by default — that is the intended behaviour, but it
must be reported per INV-3 and called out in the DP-5 rationale. Deferring
`--score-scale` normalisation is a best-effort STOP: revisit only if real
multi-list workflows demand it (scoring §5.2), don't build it speculatively.

---

### FS-3 — Grow `benchmarks/fill_quality/` into the standing regression it claims to be

**What / surface (benchmark).** Three additions, no methodology change:
**(a)** the native `--min-score 50` column (`cs_minscore_<g>.json`) wired into
`run.sh` + `score_fill.py` — the FS-1 acceptance gate; **(b)** a broader mask
spread — the harness's own caveats section calls the current set "four easy
grids + one hard grid" and asks for standard 11×11/15×15 American masks (the
bundled `blocked_15a` and `blocked_13b` are free candidates) with a
**completion-rate column** so the search-power axis is tracked per-grid,
per-threshold (the quality/feasibility frontier: completion × min-score);
**(c)** a **CI-runnable subset**: today the whole harness is gated on two
external deps (ingrid_core via cargo, STW via download), so nothing about fill
quality runs in `make test` — a bundled-fixture mode (the FS-1 original scored
fixture + the small grids) can assert AC-FILL-5/-7 in plunit/golden form,
keeping the STW/ingrid head-to-head as the on-demand deep benchmark.

**Why now.** Every other proposal in this document says "measured against
benchmarks/fill_quality/" — that only means something if the harness (1) has
the native column, (2) covers grids beyond the easy four, and (3) can't drift
silently because no CI runs it. The README's "next step" already commits to
(a); (b) is its own recorded caveat; (c) follows the repo's discipline that
regressions live in `make test` (AGENTS.md Verification).

**Tags.** Effort **low→medium** · leverage **high** (it multiplies every other
row; a quality regression the harness can't see doesn't exist per INV-3's
spirit) · identity-fit **direct** (deterministic, measured-first — the house
methodology from the arrange/fill perf campaigns).

**Spec impact.** None — benchmarks aren't spec surface. The CI-subset fixtures
ride FS-1's test assets. If (b) adds masks to `grids/`, they follow §8.3's
existing rule (lint-validated, license-clean, original).

**Measurement.** It *is* the measurement. Concrete deliverable: the results
table gains rows per new mask and a completion matrix
(grid × {no-prune, ≥1, ≥30, ≥50}); numbers must be reproducible against a
pinned STW snapshot date (already the convention: "snapshot 2026-07-01").

**Risks / traps.** STW drifts (snapshot-date the results, as now). Don't let
the CI subset quietly *replace* the head-to-head — ingrid remains the external
reference point (landscape §A calls it the sharpest benchmark target). Keep
`junk=0 by construction` documented: the cs dict is STW-derived, so junk shows
up as score-0/low-mean, not as a junk count.

---

### FS-4 — Search power: a decision pass for output-changing fill search (the harder, separate gap)

**What / surface (`fill` engine + spec).** A **decision pass (DP-6-class), not
an implementation**: choose which output-changing search lever(s) to spec so
`fill` can complete standard full-slot grids. The menu, from weakest to
strongest medicine:

1. **`--budget N` CLI knob** — the perf campaign flagged the hardcoded
   `fill_budget(800_000_000)` (verified still hardcoded, `fill.pl:90`) with no
   CLI override as a pre-campaign fix; it never landed as a flag. Cheapest
   honest lever: lets a user *buy* completion with inferences,
   deterministically. Not a fix (the tail is exponential), but the right
   escape hatch. Needs a §5 CLI-contract edit (AC-CLI-2: no undocumented
   flags).
2. **Deterministic seeded restarts** (`--seed N`, arrange-§7.6 precedent) —
   the fill perf campaign explicitly named the "randomized-restart product
   feature" as the envelope lever for heavy-tail instances and put it out of
   its own scope. Opt-in, reproducible per seed, INV-2-clean by the same
   construction as arrange's.
3. **Crossing-aware pruning / forward-checking** — the real ingrid-class fix.
   The campaign recorded that fill's candidate counting *deliberately* ignores
   crossing consistency because true forward-checking "would change counts,
   change MRV selection, and change output: off-limits" under the perf
   ratchet's byte-identity charter. Off-limits *for a perf campaign* — but a
   product decision pass can change the contract: regenerate goldens once,
   version the engine, keep INV-2 (still deterministic — same input ⇒ same
   bytes; just *different* bytes than v1).

**Why now.** The benchmark's second finding: crosswordsmith cannot complete
`blocked_13a` within budget **with any dict**, scored or not — a completion
gap the README now names as a limitation, orthogonal to scoring and "the
harder engineering gap". Without it, FS-1 is a quality feature for easy/medium
grids only, and the FS-1×hard-grid interaction (prune reduces feasibility)
makes high thresholds fail more often. The perf campaign already proved the
complementary negative: per-node speedups "buy latency on completing
instances, never completion" — so no further constant-factor work will close
this; only tree-shape changes can.

**Tags.** Effort **high** (option 3; options 1–2 are low/medium) · leverage
**high** (completion on standard grids is table stakes — ingrid does 13×13 in
~10s) · identity-fit **yes with care** (deterministic CSP improvement is
squarely in-identity; the *contract churn* is the cost).

**Spec impact.** **Needs a new decision pass — say so and stop there.**
*(Since resolved by **DP-6**, 2026-07-15: options 1–2 adopted as design-spec
§8.4b [LOCKED, build pending] — measured as control/variety levers, not
completion fixes (a ×2/×4/×10/×20 budget ladder and an 8-salt ordering probe
both fail `blocked_13a` at min ≤ 30, the FS-3(b) frontier's sole vs-ingrid gap
row); option 3 elevated to the named necessary path and deferred to its own
decision pass with exactly the contract-churn plan below. Probe evidence:
[`../../benchmarks/fill_quality/`](../../benchmarks/fill_quality/README.md).)*
Option 1 edits the LOCKED §5 CLI contract; option 2 mirrors §7.6 into §8.4;
option 3 revises §8.4's "deterministic (dictionary order + lowest-start
tiebreak)" search description and AC-FILL-3's byte-stability across the
version bump, plus golden regeneration. None of this can ride FS-1's build
(scope discipline: §8.4a is one filter + one sort key, explicitly "not a new
engine").

**Measurement.** The harness's completion table is the gate: `blocked_13a`
(and FS-3's added 15×15 masks) must flip from "budget-exhausted ~20s" to
completed, at stated `--min-score` levels, within stated inference budgets.
Track against ingrid's ~10s/min-30 reference. The perf ratchet
(`fill_baseline.json`) must be re-baselined, not silently broken, if option 3
changes counts.

**Risks / traps.** *Determinism:* wall-clock anything stays banned (budget is
inference-count; restarts schedule on inference milestones, not time).
*Golden churn:* option 3 invalidates every fill golden + `fill_identity.sha256`
— one deliberate pass, documented, INV-2 preserved going forward. *Scope
creep:* don't let "search power" swallow FS-1; they are separate axes and
separate passes by design. *Honesty:* even ingrid times out on `blocked_13a`
at min-50 — some grid×threshold combinations are genuinely infeasible-in-
practice; report, don't chase (INV-3).

---

### FS-5 — `xword`: post-hoc fill-quality scoring of any imported puzzle (rider on the planned `stats` verb)

**What / surface (`xword`).** Productize the benchmark's post-hoc trick
(`score_fill.py`): given any puzzle `xword` can read (native/ipuz/exolve) plus
a **user-supplied** `word;score` list, report n / mean / min / below-threshold
per entry and in aggregate — the same shape as FS-1's engine report, but for
*imported* puzzles the engine never filled. Rides the already-prioritized
`stats`/`inspect` verb (STATUS priority #2;
[xword-breadth-expansion](./xword-breadth-expansion.md) item 1) as one more
deterministic, read-only stat block, rather than being a new verb.

**Why now.** The mechanism is already proven in-repo — `score_fill.py` scores
fills post hoc "even though crosswordsmith emits no scores" — and CrossFire's
Grid/Final quality metrics (landscape §D) show constructors expect a
whole-grid quality readout, not only a fill-time one. It extends the quality
axis to the entire interchange hub (score a downloaded `.puz`-via-convert, an
Exet export, a hand-made grid), which no convert-only library offers — exactly
the "compete on the hub, not conversion breadth" strategy (landscape §B/§F).

**Tags.** Effort **low** (port `score_fill.py`'s logic onto xword's existing
parsers) · leverage **medium** · identity-fit **direct** (read-only,
deterministic, hub-amplifying).

**Spec impact.** xword-side decision pass — it lands inside the `stats` verb's
(still-needed) spec in `xword-status.md`/breadth plan, not in design-spec §8.
No engine change; engine↔xword parity is unaffected (report only, no layout
mutation).

**Measurement.** Parity gate against the harness: `xword stats --dict <stw>`
over the benchmark's `cs_*.json` fills must reproduce `score_fill.py`'s
numbers exactly. Longer term the harness can delegate its scoring step to this
verb (one scorer, two consumers).

**Risks / traps.** The scored list is always user-supplied (`--dict`-style
point-at); bundle nothing (INV-4 — STW is NC/SA, and there is still no
bundleable scored list). Normalisation must match the engine's (A–Z,
strip spaces/hyphens) or the two reports will disagree on multi-word entries
— pin with shared fixtures.

---

### FS-6 — Scored-dict inspection tooling + the two license errands

**What / surface (wordlist tooling; plus pure research).** Two cheap enablers:

- **(a) `dict-stats` style inspection** (home: most naturally a second rider
  on xword's `stats` machinery, or a small script under `benchmarks/`/`tools/`
  — the decision pass picks): given a `word;score` file, report its scale
  (observed min/max/percentiles), score-0 count, unscored-line count, and
  malformed lines. This *operationalizes* FS-2(i)'s "native units, no magic"
  stance: crosswordsmith refuses to auto-normalise, so give users a one-command
  way to *see* their list's scale before choosing `--min-score`. It also picks
  the explicit, documented aggregation rule if lists are ever merged — the
  research found every upstream's merge formula undocumented, so "imitate the
  ecosystem" is not an option (scoring §2).
- **(b) The two license errands** (research, hours, no code): verify the
  **MIT Collaborative Word List's scored-status + license** — the *only*
  candidate for a cleanly-bundleable scored list, named the "highest-value
  follow-up" (scoring §6) — and re-fetch **Peter Broda's terms** (unfetchable
  this run, expired TLS cert; currently UNCONFIRMED — do not rely on it until
  read). Keep a watch on christophsjones (currently unlicensed →
  all-rights-reserved, NOT permissive — §D stands corrected).

  > **(b) DONE (2026-07-15)** — recorded in the scoring doc's
  > [follow-up section](../research/wordlist-scoring-2026.md#follow-up-verification-2026-07-15--the-two-license-errands):
  > **CWL confirmed MIT + scored 0–100** (567,657 entries; the first
  > bundleable scored option — stale since 2023-02, band meanings
  > undocumented; bundling remains its own decision pass), **Broda hardened
  > to NO** (page carries no terms text at all), christophsjones re-checked
  > still unlicensed. Bonus corroboration for FS-2(ii): CWL also uses a
  > score-0/1 junk-blocklist floor (`ABUSE;0`, `ACKROYD;0`).

**Why now.** (a): the scale-heterogeneity finding makes user error the top
usability risk of FS-1 (a 5–60-scale list + `--min-score 50` = near-max
strictness, silently). (b): if the MIT list is confirmed MIT *and* scored, the
"bundled default stays unscored" posture gets its first real alternative —
that's a strategic unlock worth a few hours now, and it can run in parallel
with everything else.

**Tags.** Effort **low** · leverage **medium** (a), **medium-high option
value** (b) · identity-fit **direct** (deterministic reporting; license-clean
diligence is the INV-4 discipline itself).

**Spec impact.** (a) needs a small decision pass wherever it lands (xword plan
or a tools note) — it is new scope, however small. (b) is not scope: it
resolves open questions already carried in scoring §6. **If (b) confirms the
MIT list, actually bundling it is a separate, new decision pass** — §8.4a's
"bundled default stays UKACD18" is LOCKED and changing the default lexicon
reopens the DP-2/DP-4 lexicon posture; this document deliberately does not
propose that ahead of the facts.

**Measurement.** (a): byte-deterministic reports; validated against known
lists (its STW report must show the documented shape — 0–100 with a score-0
blocklist mass; its christophsjones report ~1–50). (b): a dated update to
`wordlist-scoring-2026.md` §4/§6 — the same doc-discipline as every verdict
there.

**Risks / traps.** (a) must embed no list data (fixtures stay original). (b):
treat "no LICENSE file" and "cert expired" as hard NOs until positively
resolved — the research's exact posture; do not soften verdicts to unblock
bundling.

---

### FS-7 — Scored `arrange`: tiebreak + best-effort drop-order signal (keep dormant; spec only when demanded)

**What / surface (`arrange`).** An optional scored-dict lookup (`--dict FILE`,
same `word;score` reader as FS-1) used **only** as (a) a drop-order signal in
`--best-effort` — among candidate drops that are otherwise equal under the
existing lexicographic most-placed → reward rule, drop the lowest-scored word
first — and (b) a final tiebreak between equal-reward layouts. Words absent
from the dict get the uniform score; score **never** gates placement — every
legally-placeable user word is still placed (closed-set inversion intact), and
score never invents words.

**Why now — honestly, weakly.** It's the recorded natural follow-on
(landscape §F item 2, "effort low"; §8.5 carries the backlog row, explicitly
"dormant"). But the same research flags the open question of whether scored
`arrange` is even coherent given the closed-set model (landscape, open
question 2), and the mechanism only ever acts on *ties* — in a verb whose
objective (§7.2 capped reward) already went through a calibration pass that
found its own cap largely inert on realistic inputs. Expected behavioural
surface: small.

**Tags.** Effort **low→medium** (new `--dict` plumbing into a verb that has
none) · leverage **low→medium** · identity-fit **yes, narrowly** (tiebreak/
drop-order only — the §D boundary).

**Spec impact.** **Needs its own OD-8 decision pass** (the §8.5 row says so),
and that pass must explicitly distinguish this from §3's excluded "priority /
per-word importance scores in `arrange`": that exclusion is *user-asserted
importance as placement input*; this is *dictionary-intrinsic quality as a
tie/drop key*. Adjacent enough that skipping the distinction would read as
silently reopening a §3 non-goal. INV-1 requires the score to come from a
dictionary lookup keyed on the answer, never from `meta` — even though
`arrange` users own their `meta`, a `meta.score` route is the trap to refuse.

**Measurement.** `benchmarks/fill_quality/` doesn't cover `arrange`; this
would be plunit-first (drop-order: a fixture where two words are droppable and
the lower-scored one goes; tiebreak: two equal-reward layouts split by score)
plus, if pursued, a small post-hoc scoring column for best-effort arrange
fixtures. Determinism via existing goldens/fuzz (unscored path byte-identical,
mirroring AC-FILL-6's construction).

**Risks / traps.** Determinism is easy (total order over existing keys) but
the payoff may not exist — recommendation: **leave dormant until FS-1 ships
and real usage asks for it.** Building it now is speculative scope.

---

## 4. Explicitly cut (guardrail violations or negative-value)

- **Bundling any currently-known scored list.** STW / Abate: CC BY-NC-SA
  (point-at only). XWord Info: paid + personal-use-only — cannot even be
  referenced as a `--dict` asset to ship. christophsjones: no LICENSE →
  all-rights-reserved. Broda: no published terms (verified 2026-07-15 — hard
  NO). *Exception opened by FS-6(b)'s completed verification:* the
  **Collaborative Word List is MIT + scored**, so bundling *it* is no longer
  an INV-4 violation — but it stays cut from this plan because it needs its
  own decision pass against §8.4a's LOCKED default-lexicon posture (and a
  staleness/band-semantics assessment).
- **Auto-detecting and normalising score scales.** Rejected in the research
  itself ("non-deterministic-feeling, hides data", scoring §5.2 option (c)).
  Even the opt-in `--score-scale MIN:MAX` (option (b)) is deferred until
  multi-list workflows demand it — FS-6(a)'s inspection report is the v1
  answer.
- **Score in the output layout, or score via `meta`.** INV-1: score is a
  dictionary property; the emitted layout stays score-free (the report is a
  sidecar). Threading a score through placement from `meta` is the exact
  failure mode INV-1 exists to prevent.
- **Whole-grid score *optimization*** (CrossFire "Final Score"-style objective
  search: fill, then hill-climb the mean). A new engine and a new anytime-
  optimization contract, in a search space the perf campaign showed is already
  budget-bound on hard instances; the §8.4a prune+order+report gets the
  measured parity result without it. Not proposed; if ever wanted, it is a
  full decision pass on a par with FS-4 option 3.
- **Scoring as word *selection* in `arrange`** (drop a placeable word because
  it scores low, or pull replacement words from a dictionary). Breaks the
  closed-set inversion (§1/§3); `arrange` must place all user words.

## 5. Recommended sequencing

1. **DP-5 micro-pass (FS-2) first, then build FS-1 + the FS-3 harness growth
   as one milestone.** Rationale: FS-2 is a spec-text edit that determines
   FS-1's default constant — resolving it first means the default ships once.
   FS-1 is the highest-leverage item on the board and everything else is
   measured against it; FS-3(a)+(c) are its acceptance gate and belong in the
   same commit series (the harness README already specifies the gate).
2. **Kick off FS-6(b) — the two license errands — immediately and in
   parallel.** Hours of research, no code dependency, and the MIT-list answer
   changes the long-term lexicon strategy.
3. **FS-4 decision pass next** (not necessarily the build): with FS-1 landed,
   the harness's completion×min-score matrix (FS-3(b)) will show exactly how
   often quality thresholds hit the budget ceiling — that evidence sizes the
   DP-6 options (budget knob and seeded restarts are cheap and likely worth
   spec'ing regardless; forward-checking is the big swing to decide
   deliberately).
4. **FS-5 rides the xword `stats` pass** whenever that (already-★'d) item gets
   its decision pass — it should be in that spec from day one, since it's the
   cheapest way to spread the quality axis across the whole hub. FS-6(a) lands
   with it if `stats` is its home.
5. **FS-7 stays dormant** until FS-1 usage produces a concrete ask.

The quality axis is steps 1–2 (+4 for visibility); the search-power axis is
step 3 alone. Do not merge them: §8.4a is locked as "not a new engine", and
the benchmark's central lesson is that they fail independently and must be
measured independently.

## 6. Reconciling the open refinements and license gaps

- **Refinement (i) — scale units:** adopt research option (a): `--min-score`
  is in the dict's native units; documented, no normalisation. Recorded via
  DP-5 (FS-2) as a one-line §8.4a amendment + an AC-FILL-8 wording fix
  ("0–100" → "native units"). FS-6(a) gives users the scale-visibility that
  makes this honest in practice.
- **Refinement (ii) — score-0:** adopt research option (b): default prune
  `score ≥ 1` (exclude only the blocklist floor), uniform score for unscored
  entries pinned ≥ 1, docs recommending `--min-score 50` as the practical
  clean floor. This *changes* DP-4's "default N = 0" sentence → it is the part
  of DP-5 that genuinely needs the recorded pass, with the reframed invariant:
  the default "never removes *usable* words".
- **MIT Collaborative Word List:** ~~scored-status and license both
  UNCONFIRMED~~ **VERIFIED 2026-07-15 (FS-6(b) executed): MIT + scored
  (0–100, 567,657 entries)** — the only bundleable scored option now exists.
  **Bundling it is still a new decision pass** against §8.4a's locked lexicon
  posture, weighing staleness (content frozen 2023-02-12) and undocumented
  band semantics — not an automatic consequence.
- **Peter Broda's terms:** ~~UNCONFIRMED (expired cert blocked the fetch)~~
  **VERIFIED 2026-07-15: hard NO** — the page publishes no terms at all;
  Broda-derived data is treated as unusable, same as christophsjones (no
  LICENSE → all-rights-reserved; re-checked 2026-07-15, unchanged).
- **With both gaps now closed,** the shipped posture stands and is *validated*
  by the research, not merely tolerated: bundled default = unscored UKACD18;
  scored fill = `--dict` point-at opt-in; recommended list = STW; recommended
  threshold = 50 (documented convention: STW FAQ, XWord Info FAQ, ingrid's
  default).
