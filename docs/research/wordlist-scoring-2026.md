# Word-list scoring — how "quality score" is defined, sourced & licensed (2026)

**Status:** research deliverable, 2026-07-15. Deep-research pass (5 angles, 18
primary sources fetched, 25 claims adversarially verified 3-0/2-1, 2 refuted).
**Follow-up license errands completed later the same day** — see the
[follow-up section](#follow-up-verification-2026-07-15--the-two-license-errands):
the Collaborative Word List is **CONFIRMED MIT + scored** (the first cleanly-
bundleable scored option), Peter Broda's list hardens to **NO** (no published
terms at all).
Feeds the just-spec'd **scored fill** (`fill --min-score`, design-spec
[§8.4a / DP-4](../design-spec.md)) and refines
[`setter-tool-landscape-2026.md` §D](setter-tool-landscape-2026.md#d-american-scored-fill--the-load-bearing-gap-and-how-to-close-it).
Where a claim is convention/folklore rather than documented, it says so.

## TL;DR — the four things that change our design

1. **"50 = usable/clean baseline" is a *real documented convention*, not folklore.**
   Three independent primary sources converge on 50 as the trust floor (STW FAQ,
   XWord Info FAQ, ingrid_core's `--min-score` default of 50). Crucially it means
   **"clean / acceptable," explicitly NOT "lively / exciting"** — STW: *"we trust
   answers scored at 50 to be clean (not necessarily exciting, possibly needing
   fair crosses)."* → **Recommending `--min-score 50` as the practical setting is
   well-grounded, and matches our benchmark's `>=50` column.**

2. **There is NO universal 0–100 scale. Our spec's "integer 0–100" is only one
   list's convention.** Documented ranges: STW/Broda **0–100**, XWord Info **5–60**,
   Crossword Compiler default list tops out at **50**, christophsjones **~1–50**.
   The hypothesis that Crossword Compiler originated a "50 = neutral midpoint"
   0–100 convention was **refuted** (in CC's default list 50 is the *ceiling*). →
   **`--min-score N` cannot assume a fixed scale; ingesting a foreign list needs
   per-list scale awareness. This is a genuine spec gap (see §5).**

3. **Score 0 is not always "unrated" — on STW it's a deliberate blocklist
   (harmful/X-rated).** → **Our spec default `--min-score 0` ("include everything")
   would admit STW's intentionally-excluded entries. Needs a decision (see §5).**

4. **The license reality is starker than §D assumed.** No cleanly-*bundleable*
   scored list was confirmed: STW and Matt Abate's list are CC BY-NC-SA (point-at
   only); XWord Info is paid + personal-use-only (unusable); and
   `christophsjones/crossword-wordlist` **has no LICENSE file at all** (⚠️ correction
   to §D, which called it permissive) → all-rights-reserved, not bundleable. →
   **This *validates* the spec choice: bundled default stays the permissive,
   unscored UKACD18; scored lists are `--dict` point-at only.**
   *Update (same-day follow-up): the MIT-list gap is now closed — the Crossword
   Nexus **Collaborative Word List** is verified **MIT and scored (0–100)**, so
   a cleanly-bundleable scored option now exists (content stale since 2023-02;
   bundling it would still need its own decision pass — the shipped default is
   unchanged). See the follow-up section.*

## 1. Band-definition documentation (Workstream 1)

The single best-documented band table is XWord Info's. Note the **heterogeneous
scales** — read the "Scale" column before comparing numbers across rows.

| List / tool | Scale | Documented band meanings | 50 = ? |
|---|---|---|---|
| **XWord Info Word List** | **5–60** | 60 = assets/well-above-average · 50 = fine/acceptable · 25 = 3–4-letter hard-to-defend (partials, odd abbrevs, esoteric names) · 20 = 5-letter "gluey" · 15 = long random Roman numerals · 5 = editor-identified puzzle-killers ("avoid at all costs") | acceptable baseline (nominal default; curator downgrades from 50, rarely up) |
| **Spread the Wordlist (STW)** | **0–100** | 50 = "clean" (trusted, *not necessarily exciting*) · 40 ≈ "a mixed bag" · ≤30 ≈ "throwing darts blindfolded" · **0 = blocklisted** (harmful/X-rated) | the "clean" trust floor; STW recommends "setting a minimum score of 50" |
| **Crossword Compiler (default list)** | 0–100 *allowed*, default list uses **~10–50** | 50 = fairly common root words (top main tier) · 25 = less common words / derivatives (‑ing/‑ed) · 10 = vulgar words | the **ceiling** of the default list, not a midpoint |
| **christophsjones/crossword-wordlist** | **~1–50** | 50 = "common word or phrase you wouldn't hesitate to put in the grid" · 25 = "an acceptable word" · 2 = "lowest score in the list" · 1 = reserved semi-word placeholder | "wouldn't hesitate" tier |
| **Collaborative Word List (Crossword Nexus)** *(added 2026-07-15 follow-up)* | **0–100** (observed: full range in live use; heavy modes at 20/30/45/50/75/90) | Band meanings **undocumented** — README/blog define only the `word;score` format. Observed floor semantics: 0/1 hold misspellings + flagged-sensitive entries (`ACKROYD;0`, `ABUSE;0`, `AACK;1`) — a junk/blocklist floor in practice | 85,900 entries sit at exactly 50; 253k of 568k are ≥50 — 50-convention-compatible in practice |
| **ingrid_core (filler, not a list)** | consumes the list's scale | `--min-score <N>` prunes candidates below N; **default 50** | de-facto usable floor |

**Crossword Compiler's model matters for us:** its help states *"The absolute
magnitude of the word scores is not very important. Mostly what matters is the
relative scores, with a word higher than another being used much more frequently."*
It also exposes a *"Minimum word score"* hard floor (and per-length minima) — a
direct analogue of `--min-score`. So both mechanisms in §8.4a (score-descending
ordering **and** a hard prune) are established, documented filler behaviours. The
semantic subtlety: CC treats score as a **relative rank**; our spec treats
`--min-score` as an **absolute threshold**. Both are legitimate; they are not the
same thing.

## 2. Signals beyond curator-typed scores (Workstream 2)

Modern scores are **largely computed, not hand-typed per word** (🟢 verified 3-0):

- **STW:** *"cleanliness-based and largely automated based on frequency and venue,"*
  with 0 reserved as a blocklist floor.
- **XWord Info:** of 253,122 entries, **118,157 come from published NYT puzzles**,
  the rest manual additions — i.e. published-usage-seeded.
- **christophsjones:** *"A score is computed for each word based on the number of
  occurrences"* across NYT, WSJ, WaPo, UKACD, Peter Broda's list, Peter Norvig's
  frequency counts, smaller papers, plus ~5k handpicked words.
- **Matt Abate's list (ML):** *"An SVM is trained to separate liked words from
  disliked words, and final word scores are generated from this SVM"* — ~42k
  manually-labelled (binary approve/dislike) words over ~413k, using text
  embeddings.

The full spectrum, then: pure frequency counts → venue weighting → hand-curated
floors/blocklists → ML auto-scoring.

**Documented vs. NOT documented.** The *inputs* above are documented. The
**aggregation formula** when a word appears in several corpora (max vs. average vs.
weighted) and any **length-normalization** were **not found in any primary source**
— treat as unknown, not assumed. (This is why a naive list-merge in crosswordsmith
should pick an explicit, documented rule rather than imitate an unstated one.)

## 3. Constructor / setter spoken heuristics (Workstream 3) — **thin, declared**

This workstream produced **no surviving verified quotable heuristics from
constructor interviews / podcasts / blogs** (Nediger, Husic, Agard, Steinberg,
Cruciverb threads). The only durable "heuristics" are the ones *baked into the list
documentation* already captured above:

- **Cleanliness ≠ liveliness.** The load-bearing distinction the whole ecosystem
  encodes: a 50 is *clean/defensible*, not *sparkly*. High scores (60+, or 70–100
  on 0–100 lists) are where liveliness lives. A scored filler optimising for a
  minimum floor buys *cleanliness*; sparkle needs a *higher* target, not just a
  floor.
- **The workflow is downgrade-biased.** XWord Info's curator *"erred on the side of
  keeping as many at 50 as possible, only downgrading entries he believes will not
  be received well."* Scoring in practice is mostly *demotion of junk*, not
  promotion of gems.
- **Junk is scored, not deleted.** Bad entries stay in the list at low scores (or 0)
  so the filler *can* reach them as a last resort — the list is exhaustive; the
  *score* does the gating. This is exactly the `--min-score` model.

> **Gap, explicitly declared:** if we want quotable working-constructor heuristics
> beyond list docs, that is a separate targeted pass (podcast show-notes / blogs) —
> this run did not surface verifiable ones. Not pursued further here (best-effort
> STOP).

## 4. License & redistribution verdicts (Workstream 4)

Verdict legend: **BUNDLE** = may ship inside a permissive distribution · **POINT-AT**
= usable only via `--dict` against a user-supplied copy (cannot embed) · **NO** =
not usable for distribution at all.

| List | Scored? | License / terms | Verdict for crosswordsmith |
|---|---|---|---|
| **Spread the Wordlist** | ✅ 0–100 | **CC BY-NC-SA 4.0** (attribution, non-commercial, share-alike). Puzzles made with it *may* be sold; the list itself may not. | **POINT-AT.** NC + ShareAlike would contaminate an otherwise-permissive bundle. (ingrid_core *does* embed it — but ingrid_core is itself NC-compatible, not a license workaround.) |
| **Matt Abate's ML list** | ✅ (SVM-derived) | **CC BY-NC-SA 4.0** | **POINT-AT.** Same NC/ShareAlike constraint as STW. |
| **XWord Info Word List** | ✅ 5–60 | **Paid** ($50/yr Angel + $200 one-time for Jeff Chen's personal list) **AND personal-use-only**: *"for your personal use only… not allowed to give it away or sell it… or create derivative lists except for your own use"* under a "DO NOT DISTRIBUTE" heading. | **NO.** Cannot bundle, cannot even ship as a `--dict` reference. Only an end-user who bought their own copy could load it locally. |
| **christophsjones/crossword-wordlist** | ✅ ~1–50 | ⚠️ **No LICENSE file** (checked 2026-07-15) → all-rights-reserved by default. | **NO (as of check).** *Correction to §D, which called it permissive.* Reclassify to BUNDLE only if the author adds a permissive license. |
| **Peter Broda's wordlist** | ✅ (a scored grid-text download is offered) | **No published terms at all** *(follow-up, 2026-07-15)*: the page was finally fetched (TLS cert still expired; fetched with verification disabled — integrity caveat noted) and carries **zero license/usage/redistribution text** — only download links (Crossword Compiler `.txt` / CrossFire `.dict`, scored grid text; other variants "currently unavailable" after a 2024-01-07 server-data-loss notice). | **NO** (hardened from UNCONFIRMED). No stated terms ⇒ all-rights-reserved by default: not bundleable, and not recommendable even as a point-at source while terms are unstated. Revisit only if terms are ever published. |
| **Collaborative Word List** (Crossword Nexus — the "MIT list") | ✅ **CONFIRMED scored** *(follow-up, 2026-07-15)*: `word;score`, integer **0–100** — 567,657 entries, 0 malformed on parse; 253,136 ≥ 50; 324 at score 0 + 166 at 1 (misspellings + flagged-sensitive entries) | ✅ **CONFIRMED MIT** — in-repo LICENSE file, "Copyright (c) 2021 Crossword-Nexus"; README: "free for everyone"; launch blog: "free for all to use in any way they please". | **BUNDLE — the first confirmed cleanly-bundleable scored list.** Two caveats: band semantics are undocumented (scores observed, meanings inferred), and the **content is stale** — `xwordlist.dict` last changed **2023-02-12** (last repo push 2024-01-07; not archived). Actually bundling it is a **separate decision pass**: §8.4a's default-lexicon posture is LOCKED. |
| **UKACD18** (current default) | ❌ **unscored** | Redistributable freeware (per DP-2). | **BUNDLE (unscored).** It is the permissive default dictionary, *not* a scoring source. |
| **Ginsberg clue database** | n/a (clue DB, not a scored wordlist) | Not verified this run. | Out of scope for scored fill. |

**Bottom line on a bundleable scored default:** none confirmed. STW/Abate are
point-at-only; XWord Info is off-limits; christophsjones is unlicensed; the MIT
Collaborative Word List's scored-status is unconfirmed. **Until the MIT list is
verified as scored, crosswordsmith has no cleanly-bundleable scored dictionary** —
which is exactly why §8.4a keeps the bundled default *unscored* (UKACD18) and makes
scoring a `--dict` opt-in. The research validates that choice.

> **Superseded same day (follow-up errands, 2026-07-15):** the verification
> happened — the **Collaborative Word List is MIT + scored**, so a
> cleanly-bundleable scored option now **exists**. The shipped posture
> (unscored UKACD18 default; scored via `--dict`) remains correct *until a
> decision pass weighs bundling CWL*, and that pass must weigh the freshness
> caveat: CWL's content has been frozen since 2023-02-12, where STW is
> actively curated. Broda meanwhile hardens to NO (no published terms), and
> christophsjones re-checks as still unlicensed.

## 5. What this means for `fill --min-score` (§8.4a / DP-4)

Three findings touch the spec directly. Two are **new design questions** the
current §8.4a does not yet resolve.

1. **✅ Validated: default lexicon & licensing.** §8.4a's "bundled default stays
   permissive UKACD18 (unscored); STW opt-in via `--dict`, never bundled" is the
   correct call — there is no confirmed cleanly-bundleable scored list. *Action:*
   drop §D's "christophsjones … permissive" and "Collaborative Word List … " as
   settled facts; treat the MIT list's scored-status as the one gap worth closing
   if we want a bundleable scored default. *(Closed by the same-day follow-up:
   CWL confirmed MIT + scored — see §4 and the follow-up section.)*

2. **⚠️ New: `--min-score N` assumes a fixed 0–100 scale, but scales are
   heterogeneous (5–60, 1–50, 0–100).** A user pointing `--dict` at an XWord Info
   export and passing `--min-score 50` gets *near-maximum* strictness (60 is the
   cap), not a mid-band floor. *Options:* (a) document that `--min-score` is in the
   *dict's own units* and leave scale-awareness to the user (simplest, honest,
   deterministic); (b) add an optional `--score-scale MIN:MAX` normalisation; (c)
   auto-detect range and normalise (rejected — non-deterministic-feeling, hides
   data). **Recommendation: (a) for v1** — `--min-score` is in the dictionary's
   native units, documented, no magic. Revisit (b) only if multi-list workflows
   demand it. *This warrants a one-line amendment to §8.4a (or a DP-5 if we adopt
   normalisation).*

3. **⚠️ New: default `--min-score 0` + STW's score-0 blocklist.** On STW, 0 is
   reserved for deliberately-excluded (harmful/X-rated) entries. A default of 0 that
   *includes* score-0 words would admit them. *Options:* (a) keep default 0 but
   document that 0-scored words are included and that STW users should pass
   `--min-score 1`; (b) make the default prune `score > 0` (exclude only the
   blocklist floor) — a safer default that still "changes feasibility only by
   removing explicitly-blocklisted words"; (c) distinguish "unrated" from
   "blocklisted" in ingestion. **Recommendation: (b)** — default to excluding
   score-0 (blocklist) entries while leaving everything ≥1 in; it matches the
   ecosystem's intent and barely touches feasibility. *This too is a §8.4a
   amendment / small DP.*

**Unchanged & reinforced:** the "50 = clean floor" convention means our benchmark's
`>=50` column and a recommended `--min-score 50` are the right regression target and
the right documented default *suggestion* (distinct from the code default in #3).

## 6. Open questions / gaps carried forward

- ~~**MIT Collaborative Word List:** is it scored, and confirmed MIT?~~
  **RESOLVED (2026-07-15 follow-up): yes on both** — in-repo MIT LICENSE +
  `word;score` 0–100 across 567,657 entries. The bundleable-scored-default path
  is now open, gated on its own decision pass (and on the staleness caveat:
  content frozen 2023-02-12). *New residual question:* CWL band **meanings**
  are undocumented — scores observed, semantics inferred.
- ~~**Peter Broda's list terms:** unfetchable this run (expired cert); confirm
  before relying.~~ **RESOLVED (2026-07-15 follow-up): hard NO** — the page
  (fetched with TLS verification disabled; cert still expired) publishes no
  license/usage terms whatsoever.
- **christophsjones:** watch for a license being added; currently unlicensed
  *(re-checked 2026-07-15 follow-up: still only `README.md` +
  `crossword_wordlist.txt`, no LICENSE)*.
- **Score aggregation across corpora** (max/avg/weighted) and length-normalization:
  undocumented in every primary source — pick an explicit rule if crosswordsmith
  ever merges lists.
- **Constructor spoken heuristics:** this pass found none verifiable beyond list
  docs; a targeted blog/podcast pass would be needed to fill Workstream 3.

## Follow-up verification (2026-07-15) — the two license errands

Same-day follow-up closing §6's two license gaps (dispatched from
[`plans/fill-scoring-uplift.md`](../plans/fill-scoring-uplift.md) FS-6(b));
single-agent verification against primary sources (not the main pass's 3-vote
protocol — flagged for honesty; every claim below is from a fetched primary
artifact, not a summary).

**1. Collaborative Word List (Crossword Nexus) — CONFIRMED MIT + scored.**

- **License:** the repo (`Crossword-Nexus/collaborative-word-list`) carries an
  in-repo `LICENSE` file — standard **MIT**, "Copyright (c) 2021
  Crossword-Nexus" (raw file fetched). Corroborated by the README ("it's free
  for everyone") and the 2021-06-08 launch post on the Crossword Nexus blog
  ("This list is free for all to use in any way they please").
- **Scored-status:** the full `xwordlist.dict` was downloaded and parsed:
  **567,657 entries, all `word;score`, 0 malformed**, integer scores
  **0–100** (observed min 0, max 100). Distribution: 253,136 entries ≥ 50;
  heavy modes at 20 (97k), 30 (62k), 45 (59k), 50 (85.9k), 90 (46.6k).
  **Score 0 (324 entries) / 1 (166 entries) function as a junk/blocklist
  floor** — observed contents are misspellings (`ACKROYD;0`,
  `ACHEIVEMENTTEST;1`) and flagged-sensitive words (`ABUSE;0`, `ADDICT;0`) —
  which independently corroborates the §5.3 "default prune `score ≥ 1`"
  recommendation for a third list beyond STW.
- **Band semantics: undocumented.** Neither README, launch post, nor repo
  files define what the numbers mean; the contribution pipeline (GitHub
  Actions) validates format/duplicates/size (>425,000 entries), not meaning.
- **Freshness caveat (load-bearing for any bundling decision):**
  `xwordlist.dict` last changed **2023-02-12**; last repo push 2024-01-07 (a
  trivial merge); repo not archived. The list is maintained-in-name,
  **stale-in-practice** — vs. STW's active curation.
- **Implication:** the "no cleanly-bundleable scored list exists" premise is
  now false. §8.4a's shipped posture (unscored UKACD18 default) stays valid —
  **bundling CWL is a new decision pass** (default-lexicon posture is LOCKED),
  which must weigh MIT-cleanliness + 0–100-scale fit against staleness and
  undocumented bands.

**2. Peter Broda's wordlist — hardened UNCONFIRMED → NO.**

- The terms page's TLS cert is *still* expired (same failure as the main
  pass); the page was fetched read-only with certificate verification
  disabled — an integrity caveat, acceptable for reading public prose, noted.
- The page contains **no license, usage, or redistribution text at all**:
  only download links (Crossword Compiler `.txt` / CrossFire `.dict` "grid
  text SCORED"; full-text and unscored variants "currently unavailable") and
  a 2024-01-07 notice of server data loss.
- No stated terms ⇒ all-rights-reserved by default. **Verdict: NO** — not
  bundleable, and not recommendable as a documented point-at source while its
  terms are unstated.

**3. christophsjones/crossword-wordlist — re-checked, unchanged.** Top-level
files are still `README.md` + `crossword_wordlist.txt` only; **no LICENSE** ⇒
all-rights-reserved stands.

## Sources (all primary unless noted)

- Spread the Wordlist FAQ — https://www.spreadthewordlist.com/faq · /faq-old
- XWord Info Word List + FAQ — https://www.xwordinfo.com/WordList · /WordListFAQ
- Crossword Compiler help (scoring / grid filling) —
  https://www.crossword-compiler.com/en/help/html/aboutwordlistscoring.htm ·
  /gridfilling.htm · /changingwordscores.htm
- ingrid_core — https://github.com/rf-/ingrid_core (README + `--min-score` default 50)
- christophsjones/crossword-wordlist — https://github.com/christophsjones/crossword-wordlist (no LICENSE file, 2026-07-15)
- Matt Abate wordlist (ML/SVM scoring, CC BY-NC-SA) — https://github.com/mattabate/wordlist
- Peter Broda wordlist — http://www.peterbroda.me/crosswords/wordlist/ (fetched 2026-07-15 follow-up with TLS verification disabled — cert still expired; page carries no terms text)
- Collaborative Word List — https://github.com/Crossword-Nexus/collaborative-word-list (LICENSE = MIT © 2021 Crossword-Nexus; README; `xwordlist.dict` raw download parsed 2026-07-15; freshness via GitHub API) · launch post: https://crosswordnexus.com/blog/2021/06/08/the-collaborative-word-list/
- T Campbell, "Shopping for Wordlists" (blog, context) — https://tcampbell.substack.com/p/shopping-for-wordlists
