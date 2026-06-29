# Evolving crosswordsmith for Cryptic Crossword Setters

*A decision-useful report distilling setter methodology, the tooling landscape, and concrete feature guidance for a small Prolog layout tool.*

---

## Executive Summary

1. **Clue-writing is the heart of setting, and it is barely automated anywhere.** Every leading tool (Exet, Crossword Compiler, the late Sympathy) automates grid design and fill but leaves clue *composition* to humans — they only provide *research* (wordplay breakdowns, indicator/abbreviation lists, prior-clue lookup). This is the single biggest gap between crosswordsmith (pure layout) and what setters need, but it is also a *bounded, assistive* surface — not an AI-clue-writer.

2. **Cryptic grids obey hard, codifiable rules that crosswordsmith currently does not enforce or even validate.** 15×15, 180° rotational symmetry, ≥50% checking in every word, no triple unches, double-unches not at word ends (the Times rule), min length 3, full interconnectivity. crosswordsmith already enforces *some* of these (min-half checking with correct `ceil` rounding, connectivity, no-merge), but explicitly disclaims symmetry and min-length. A **grid linter/validator** is the highest-leverage, lowest-risk, most in-scope first feature. *(One caveat established by a follow-up discovery pass: 180° rotational symmetry is a **strong default with theme-gated exceptions**, not an exceptionless universal rule — it should be a per-publication-profile setting, not a hard FAIL. See "Addendum — Symmetry: how universal is the constraint?" below.)*

3. **Fill quality in the real ecosystem comes from *scored word lists*, not just the search algorithm.** Crossword Compiler, Crossfire, Exet and Crosshare all rank candidate fill by a per-word quality/popularity score with a minimum-score cutoff. crosswordsmith ranks only by placed-letter footprint.

4. **The closed-set model inverts the ecosystem.** Because crosswordsmith places a *fixed user-supplied word list* (a website's links, a setter's seed words) with no filler dictionary, the standard "prefer high-scoring words from an over-complete dictionary" pattern does **not** apply. Scoring is only actionable for *drop-order* and *placement tiebreaks* and *human-facing annotation* — a genuine design constraint, and also a differentiator worth keeping.

5. **Interoperability is a cheap, high-value win.** crosswordsmith emits JSON today; **.ipuz** is open JSON and a near-mechanical conversion, and **Exolve** is plain-text, git-diffable, and the *only* widely-used format carrying cryptic-specific data (bars, ninas, definition-span marking, annotations). Do **not** hand-roll binary **.puz** in Prolog; reach it (and .jpz/PDF) by shelling out to `kotwords` from a native .ipuz emitter.

6. **Ninas are a pure layout feature** (fixed letters along an ordered cell path), squarely in scope — but the closed-set model makes guaranteed placement hard, so the realistic deliverable is *seed + best-effort + feasibility report*, not a dictionary-backed fill-around-the-Nina.

7. **Nobody offers git-native, diffable puzzle source with branch/PR review.** crosswordsmith is a CLI emitting deterministic, sorted-key JSON — it is unusually well placed to own this collaboration niche.

8. **Stay a layout-and-fill engine that *feeds into* clue tools; don't clone Crossword Compiler.** The defensible position is: free, cross-platform, scriptable, open-format, JSON-/text-native. Lean into that rather than chasing GUI parity, dictionaries, or a solving applet.

---

## How cryptic setters actually work

The setter pipeline is a distinct, ordered workflow — but with one crucial twist: **clue-writing is not strictly post-fill**. Setters write a handful of clues *during* the fill and treat "cluing potential" as a fill-time selection criterion.

| Stage | What happens | Notes for crosswordsmith |
|---|---|---|
| **0. Theme / Nina (optional)** | Choose a hidden message or theme; theme/Nina words become *seeds*. | Pure layout constraint — in scope. |
| **1. Grid choice / design** | Pick a 15×15 (or stock) grid; most newspapers *force* a pre-designed grid library (legacy of hot-metal typesetting). The Independent lets setters design their own. | Grid legality is codifiable; stock-grid libraries are a real asset. |
| **2. Seed placement** | Place theme/Nina words and any answers with a clue idea already in mind. | "Seed-then-fill" is universal. |
| **3. Grid fill** | Fill **longest-to-shortest**, selecting words for **cluing potential** (charade-breakability, container/contents, vowel ratio, common letters, anagram fodder), avoiding obscure/compound words. Run candidates through the Internet Anagram Server. | crosswordsmith fills for *density*, not cluability — a key divergence. |
| **4. Concurrent clue writing** | Lock in ~6 of 28 clues during fill (Anax); finish the grid; then clue the rest systematically, spreading device types and capping any one (anagrams <~25%, hidden 1–2/puzzle). | Out of crosswordsmith's current scope, but the device-distribution *count* is layout-adjacent. |
| **5. Clue editing** | Refine surfaces; check fairness against a Ximenean or Libertarian house style. | Linting is partly mechanical, partly judgement. |
| **6. Test-solving** | Share a private **preview link** with a test solver to catch "ambiguities, unfairness, and outright errors." Capture per-clue triage (liked / iffy / flawed + reason). | Feedback is currently social/prose; a structured record is an open niche. |
| **7. Submission** | Email a `.puz`/`.ccw`/`.cfp` to an editor (NYT wants a specific PDF). De-facto audition against house style; first puzzles rarely accepted. | Export formats matter; submission infra does not. |

Key methodologies named in the research:
- **Seed-then-fill** and **longest-to-shortest** ordering (reduces dead-end backtracking; long slots are harder to clue, so they need the most options).
- **Cluing-potential fill selection** (Stella Zawistowski's explicit checklist: does it break into charades? nest a word? anagram to something fun? high vowel ratio + common letters? avoid compounds).
- **Nina embedding** (hide a message along perimeter / column / diagonal / unchecked cells; abandon it if it forces unacceptable obscurity; puzzle must be solvable *without* noticing it).
- **Test-solve-and-iterate**, structured peer critique (Rookie Corner's liked/disliked/flawed triage + numeric "commentometer"; CCCWC's ranked vote-and-comment; Azed's Prize/VHC/HC tiers).

---

## Landscape of existing tools

| Tool | Platform | Cost | Grid fill | Clue tools | Cryptic-specific support | Export formats |
|---|---|---|---|---|---|---|
| **Crossword Compiler (Pro)** | Windows (Mac via VM) | Subscription ~A$45–450/yr; one-off bundles to ~A$1,455 | AutoFill + Pro Grid Filler (scored 400k list, per-length min scores, multi-list, diagnostics) | Anagram gen, indicator pop-ups, clue DB (SQLite), AI simple-clue gen, WordWeb/UKACD/Chambers add-ons | Blocked **and barred** grids; American/Cryptic lattice modes; rich symmetry; grid library manager | .puz, .ipuz, .jpz, CCXML, PDF, RTF, images, HTML5 applet |
| **Sympathy + TEA + Wordplay Wizard** *(WITHDRAWN ~2019)* | Windows | n/a (off market) | Interactive + up to 10 auto fills; 8 dicts in priority order | **Wordplay Wizard**: per-answer wordplay-type breakdown + indicators; clue DB | Blocked/barred/diagramless; gimmicks (letters latent, misprints, Playfair); TEA pattern language | .ipuz, RTF, HTML, images |
| **Exet** | Browser (local-only) | **Free, MIT** | Beam search (K=64), viability scoring, popularity threshold, stem-dedup | **Research tab** (anagram/charade/deletion/homophone/Spoonerism gens; hidden via nutrimatic; prior clues via georgeho; indicator/abbrev lists) — *no auto-clue* | British/US/**barred**/3-D grids; first-class **ninas**; definition-span marking; UKACD18-derived lexicon (270k, Wikipedia popularity) | **Exolve** (lossless), .puz, .ipuz, PDF, SVG |
| **Qxw** | Linux/Win/Mac | **Free, GPL** | Backtracking + **MRV + forward-checking** (closest to crosswordsmith); 9 dicts, PCRE filters | none (layout/fill only) | **Barred + blocked + mixed**; ~11 answer treatments (letters latent, ciphers, misprints); C plug-ins; **free lights** (Nina primitive); configurable checking math; **batch CLI `qxw -b`** | text/grid; no JSON |
| **Crossfire** | Win/Mac/Linux (Java) | $50 | Scored fill; Word Score + Grid Score; min-score cutoff; fail diagnostics | Matt Ginsberg clue DB | American-first; no native cryptic mode | .puz, .jpz, .ipuz |
| **Crosshare** | Web (open source, AGPL) | Free | Instant incremental autofill | none (puzzle-level "cryptic" tag only) | **First-class barred + mixed bar/block**; rebus, Schrödinger; solver heatmaps/analytics | .puz import/export; embed |
| **MyCrossword** | Web (UK) | Free | Fill assist + pattern matcher + anagram builder | **~25 cryptic device tags**, indicator list, abbreviations, DEF underline, explanation field, device-frequency + anagram-% score | Cryptic-first; preview-link **test-solve loop**; annotated solutions | Guardian-format (solving component MIT) |
| **ingrid / ingrid_core** | Desktop / Rust CLI (MIT core) | Free | Backtracking CSP (Spread-the-Wordlist, min-score, max-shared-substring) | none | American blocked only | .jpz; pushes to Crosshare |
| **Phil** | Browser (Apache-2.0) | Free | SAT-based suggestions | none | American only; unmaintained (2017) | .xw (JSON), .puz, PDF |

**The leaders, briefly:**

- **Crossword Compiler** is the commercial incumbent setters compare against — Windows-only, subscription-priced, with its best filling engine (scored 400k list, per-length thresholds, exhaustive recursion, fail diagnostics) *paywalled*. Its quality story rests entirely on the scored word list.
- **Exet** is the closest open-source analogue to where crosswordsmith could go, and the **reference capability map** for a setter tool. It proves what's wanted (clue-aware research, definition marking, ninas, scored fill) and where the line is (it never auto-writes a clue; all data stays local; no collaboration).
- **Qxw** is the tool *architecturally* closest to crosswordsmith — backtracking + most-constrained-first + forward-checking, with a non-interactive **deck/batch CLI** and meaningful exit codes. It validates crosswordsmith's engine and shows the empty niche: Qxw is GPL/scriptable but emits text, not JSON; Sympathy/TEA were GUI-only and are now off the market. **A CLI-first, JSON-native Prolog engine occupies genuinely empty ground.**
- **Sympathy/TEA/Wordplay Wizard** (Ross Beresford) were the UK cryptic specialists and are now **withdrawn from sale**, leaving a gap in dedicated cryptic-aware desktop tooling. Notably, UKACD18 — the standard British cryptic base list — was compiled by the same author and is freely redistributable under BSD-style terms.

---

## What features matter

### Grid construction
- **Symmetry** (180° rotational standard; 90°, axial mirror, and diagonal also seen) — a *defining* aesthetic of publishable blocked grids, but a **strong per-profile default, not a hard universal constraint**: near-mandatory for blocked dailies, relaxed for barred thematics, and overridable for theme/Nina puzzles even in mainstream outlets. Best modelled as a publication-profile setting with a per-puzzle override (see the symmetry addendum).
- **Checking rules**: ≥50% of every word checked; no triple unches; double-unches not at word starts/ends (Times); avoid the 5-letter unch-check-unch-check-unch pattern; odd-vs-even checking style.
- **Profiles**: BLOCKED-UK, BARRED/Ximenean (≥75% checking, per-length unch table, no blocks), and AMERICAN (every cell checked) are genuinely distinct rule sets.
- **Min length ≥3**, full interconnectivity (no islands), ~1:3 block ratio, odd dimensions (**15×15 dominates** blocked dailies; the full observed range is ~8×10–24×24 — see "Addendum — Grid sizes in the wild").
- **Stock grid libraries** (named publication templates) and **per-publication house-style profiles**.

### Fill / word lists
- **Scored word lists** in the standard `word;score` line format; minimum-score cutoff; per-length thresholds.
- **Popularity/commonness scoring** (Wikipedia-derived, as Exet) and **stem-duplication avoidance**.
- **Cluing-potential scoring** (charade-breakability, container potential, vowel ratio, common-letter weight, compound penalty) — *no existing tool emits this as a number*.
- **Multiple dictionaries with per-light selection** and theme/regex filters.
- **Fill diagnostics**: which slot/area is over-constrained, why a word couldn't be placed.

### Clue-writing (assistive, not generative)
- **Wordplay decomposition** of an answer (charade, container, deletion, reversal, hidden, alternation, anagram, homophone/Spoonerism).
- **Indicator + abbreviation datasets** (device-specific; directionality matters for reversals).
- **Prior-clue precedent lookup** (georgeho corpus).
- **Definition-span marking** (`~{...}~`) + wordplay annotation field.

### Test-solving / QA
- **Preview/test-solve loop** with per-clue triage (liked / iffy / flawed + reason).
- **Fairness linting** (Ximenean): definition present + POS match; definition-by-example flagged; recognised abbreviations; no indirect anagrams; indicator-present for required devices; link-word directionality.
- **Device-overuse analytics** (anagram <~25%, hidden 1–2, no repeated wordplay element).

### Collaboration
- **Clue stockpile** keyed by answer (Anax's ~71-clue Excel; CC's SQLite clue DB), with rough/ready status and freshness.
- **Git-native diffable puzzle source** + branch/PR review (an open niche).

### Publishing / export
- **.ipuz** (open JSON), **Exolve** (cryptic-aware text), **.jpz** (cryptic-capable XML), **.puz** (binary, widely accepted, *cannot represent bars*), PDF.
- Hosting handoff: Exost (ipuz/Exolve/puz), Crosshare (puz/jpz), Amuse Labs PuzzleMe (jpz/xml).

---

## Where crosswordsmith fits today + gaps

**What it is.** Two engines: (A) a fixed-grid **MRV-backtracking solver** that places *all* words, and (B) a `--quality` **greedy density-construction** engine that picks its own grid size, packs words densely toward a "British lattice" aesthetic, and *drops* words it can't place. It consumes a plain word list (originally a website's links → a crossword-shaped table of contents) and emits JSON.

**What it already does right (confirmed against `quality.pl` / `crossword.pl`):**
- `word_meets_half/2` computes checking as `CC ≥ ceil(L/2)` — the *strict / Times* rounding, which correctly flags the odd-length unch-check-unch-check-unch trap. This is more rigorous than a naive `floor(L/2)`.
- `max_unch_run` exists as a tunable drop-to-satisfy floor.
- Connectivity (no islands) is enforced by construction — every placed word must cross an existing one.
- `no_word_merge` and adjacency invariants prevent flush parallel words and word-inside-word.
- Output JSON is deterministic (sorted keys via `json_write_dict`), so golden-file diffs are stable — a property that makes git-friendly review viable.
- Architecturally it is **Qxw's family** (MRV + backtracking), and CLI-first + JSON-native is an empty niche.

**Honest gaps vs setter needs (some are deliberate non-goals per `cryptic-layout-spec.md`):**
- **No symmetry** (spec D4) — the biggest conceptual gap *for blocked-daily output specifically*. Free-form packing essentially never yields the 180° rotational symmetry blocked dailies treat as near-mandatory. *But* discovery (symmetry addendum) reframes this as a **strong per-profile default with theme-gated exceptions**, not a universal hard rule — so the realistic deliverable is a *profile-selectable* symmetry validator (warn/fail per profile, with a per-puzzle `allow-asymmetry` override à la Exet), not enforced symmetry everywhere. The spec itself reports current real-input layouts sit at ~0.09–0.14 checked fraction, *far* from the ~0.5 cryptic ideal.
- **Min length 3 is a soft penalty, not a floor** (spec D5) — right for the link-TOC use case, divergent from a hard cryptic convention.
- **No scored fill** — ranks only by placed-letter footprint (`placement_key = Crossings*10000 - Growth`); no popularity, cluability, or stem-dedup.
- **No Times positional double-unch rule**, and `max_unch_run` defaults to *off* rather than 2 (triple-unch ban).
- **No bars** — the cell model is block-only (`null` = empty), so barred grids and the British barred aesthetic can't be represented.
- **No clue layer at all** — no clue text, enumeration, annotation, stockpile, or device tagging.
- **No standard-format export** — JSON only; not consumable by any solver app or editor.
- **The closed-set inversion**: with no filler dictionary, crosswordsmith *cannot* invent words to satisfy a forced Nina or a min-score cutoff the way Qxw/Exet do. This bounds what scoring and Nina features can promise. *(This is deeper than scoring/Ninas: the closed-set, place-everything, emergent-grid model is in **structural tension** with authentic cryptic layout — symmetry and ~50% checking essentially cannot be hit while placing an arbitrary fixed word set. See "Addendum — The closed-set vs. authentic-layout tension" below.)*

**Where it genuinely differentiates:** word-driven, *size-choosing*, drop-if-unplaceable construction has **no Exet/CCW equivalent** (they fill a predetermined grid and drop nothing). It's a "pack *these specific* words" tool, not a "fill *this* grid" tool — worth keeping.

---

## Recommended features for crosswordsmith

Sized for a small Prolog tool. Each line is grounded in the research.

### Must-have (highest leverage, lowest risk, clearly in scope)

1. **Cryptic grid linter / validator subcommand.** Consume the layout crosswordsmith already emits; report PASS/WARN/FAIL per rule, per word: min-length, checked fraction vs `ceil(L/2)`, longest unch run (≥3 = FAIL), Times positional double-unch, odd-vs-even checking, connectivity (already PASS), and a *symmetry deficit* report. *Reuses existing metric code (`checked_cells/2`, `word_max_unch_run`); needs no solver change; the benchmark analyzer is already ~80% of this.*
2. **Native .ipuz emitter.** Localized change in `emit_json`: invert `null`→`#`, split the merged cell object into parallel `puzzle`/`solution` arrays, group words by direction into a `clues` dict, add `version`/`kind`/`dimensions`. *Open CC-licensed spec, no new dependencies, instantly makes output consumable by Exost and convertible (via kotwords) to .puz/.jpz/PDF.*
3. **Compute and carry enumeration per word.** The answer already preserves spaces, so `(4,3)` / `(4-2)` enumerations are derivable for free — *required by every clue-bearing export format.*

### High-value (in scope, modest effort)

4. **Default `max_unch_run` to 2 (triple-unch ban) + add an optional Times positional double-unch floor.** *Matches a near-universal hard convention; one-line default change plus one floor predicate.*
5. **Scored-fill support for the `--quality` engine.** Ingest `word;score` lists; use score as the **drop-order tiebreak** (drop least-valuable first) and a **placement_key tiebreak** (equal-density placements keep high-score words better-checked). *Honours the closed-set inversion — scoring drives drop/order, not selection-from-a-pool.*
6. **Per-entry cluing-potential annotation (`meta.cluing`).** Zawistowski's checklist made numeric: charade-breakability, container potential, vowel-to-consonant ratio (toward ~0.4), common-letter weight, compound penalty. Computed at the emit-time answer→meta join so the **solver stays metadata-agnostic**. *No existing tool emits this number — a genuine differentiator. Vowel-ratio/common-letter sub-scores ship first with zero data dependencies; charade/container need a small bundled word+abbreviation list.*
7. **Exolve emitter.** Plain-text, git-diffable, and the only format carrying bars, ninas, and definition-span marking. *Serves the cryptic/setter audience and round-trips to Exet.*
8. **Nina seeding (as a layout constraint).** Accept a Nina = `{target, ordered cell path}` (perimeter / column / diagonal / explicit cells), pre-bind those cells before MRV in the **solver** engine, and emit a **feasibility report** (placed / placed-with-compromise / infeasible, with conflicting crossings). *Pure layout; prefer unchecked cells; the solver path is the clean target since the closed-set greedy engine can't guarantee placement.*

### Nice-to-have (later, scope-dependent)

9. **`min_len:3` hard floor option** for setters (keep the soft default for the link-TOC use case; always *report* sub-3 entries).
10. **House-style / grid profiles** (BLOCKED-UK, BARRED-Ximenean per-length unch table, AMERICAN every-cell-checked) selectable for the linter.
11. **Clue stockpile keyed by answer** in `meta` (`{clueText, status: rough|ready, date}`) that auto-surfaces when a word is placed — modelled on Anax + CC's clue DB. *Pure data plumbing on the answer key; stays on the layout side of the line.*
12. **Prolog-native pattern/anagram query endpoint** over a loaded wordlist (Nutrimatic-style `A/C/V/#/_` classes, `<...>` anagram, `&` intersection), ranked by score. *Maps naturally onto unification/backtracking; reuses existing string handling.*
13. **Wordplay-decomposition engine** (charade via append/3, hidden/deletion via sub_atom/5, reversal, anagram via msort) emitting ranked candidates per answer as `meta.clueCandidates`. *Backtracking gives exhaustive enumeration for free; the same code verifies and generates. Prune hard (cap fodder length, require real abbreviation/synonym leaves, top-N only).*
14. **Batch/deck-style declarative interface + shell exit codes** (0 = filled, non-zero = unfillable/error), à la `qxw -b`, for pipeline embedding.

### Out of scope (don't build)

- **Auto clue-writing / surface-reading judgement** — even Exet and Wordplay Wizard stop at presenting candidates. Be an assistant, not a generator.
- **Hand-rolled binary .puz writer** (CRC-16 in Prolog) — and it can't represent bars anyway. Reach .puz/.jpz/PDF via `kotwords`.
- **Solver-facing consumption**: interactive HTML5 applet, hosting, embedding, solve analytics/heatmaps, custom-URL blogs — all require a hosted backend.
- **Bundling paywalled dictionaries** (Chambers/Collins/Oxford). Ship/recommend only UKACD18 (BSD-style), the MIT Collaborative Word List, or CC-BY-NC-SA Spread-the-Wordlist; never the personal-use-only XWord Info list.
- **Full barred-grid *generation*** (different cell model entirely) — though the barred *checking-math lint* is a cheap, free add.
- **Real-time co-editing / a community voting platform** — large product surfaces orthogonal to a CLI tool.

---

## Design principles

1. **Stay a layout-and-fill engine that feeds clue tools; don't become a suite.** The defensible identity is free, cross-platform, scriptable, CLI-first, open-format. Compete on automation and interoperability, not GUI parity or dictionaries.

2. **Keep the solver metadata-agnostic.** Clue text, annotations, cluing scores, and review flags ride in `meta` (an opaque passthrough joined to the answer at emit time) — never thread them through placement or clue-numbering. This preserves the single-shot model and changes the JSON contract only by additive keys.

3. **Separate layout from clueing explicitly.** Layout features (grids, checking, ninas, enumerations, stockpile *storage*) are in scope; clue *composition* and subjective fairness are not. crosswordsmith is the answer-keyed container; humans (or Exet) do the craft.

4. **Honour cryptic grid legality as first-class, parameterised constraints — not heuristics.** Model symmetry, the unchecked-letter lattice, and the per-length unch rules as validators (and, where feasible, generators), with selectable publication profiles. A vague "cryptic aesthetic" must become measurable acceptance criteria.

5. **Respect the closed-set inversion.** With no filler dictionary, scoring drives *drop-order and tiebreaks*, not selection-from-a-pool; a Nina is best-effort + reported, not guaranteed. Always report compromises rather than silently dropping words or producing illegal fill.

6. **Interoperate via standard formats.** Native **.ipuz** (open JSON, near-mechanical) and **Exolve** (text, cryptic-aware, git-diffable); reach .puz/.jpz/PDF via off-the-shelf converters. This unlocks every major open and hosted consumer with at most one converter in between.

7. **Be Ximenean-configurable, never Ximenean-only.** Any fairness lint must support a Libertarian (Guardian/Araucaria) profile too, plus house-style toggles, or it will wrongly reject valid clues. Mark inherently-judgemental checks (surface, loose synonyms, positional indicators) as *advisory* to avoid false-positive over-rejection.

8. **Diffable, deterministic output is a feature.** crosswordsmith already emits stable sorted-key JSON; lean into "puzzle source as code" so multi-setter review can happen as a PR diff. Follow Exolve's line-per-element granularity for any canonical source format.

9. **License-clean by default.** Ship only permissively-licensed data (UKACD18 BSD-style, MIT/CC-BY-NC-SA word lists, Wikipedia CC-BY-SA abbreviations, CMUdict, ODbL georgeho corpus with attribution/share-alike); never embed paywalled IP.

---

## Open questions & confidence

Overall the research is **high-confidence**: across eleven verification objects, nearly every load-bearing, numeric, or attribution claim was confirmed verbatim against primary sources (Exet/Exolve/Qxw READMEs, Crossword Compiler help, crypticcrosswords.net, georgeho, ipuz spec). The corrections and caveats worth flagging:

- **Exet autofill internals, lexicon counts, CCW pricing & scoring scales, the Qxw deck/batch interface, the grid-rule conventions, and the Nina mechanics** were all confirmed near-verbatim. Treat these as reliable.
- **Disputed/corrected (hedge accordingly):**
  - **MyCrossword device tags**: the headline "~27" is **25** in the actual user guide. Minor.
  - **CCW Dual-symmetry as "90°"**: the *slot-restriction* facts are verbatim, but "90°" is an interpretive gloss the help page doesn't state.
  - **Exet lexicon file** is `importance-and-words.txt`, not `.tsv`; the UKACD18 license is a **custom permissive "Advanced Cryptics Dictionary" license** (not strictly BSD-2), and the "completely unrestricted" phrase actually describes **CMUdict**, not UKACD18 — but the practical "freely bundleable with notice" conclusion holds.
  - **ipuz version string**: the cited spec shows `http://ipuz.org/v1`, not v2; target the version your consumers expect.
  - **Datamuse/OneLook API**: now requires an **API key after Jan 1 2027** and caps at 100k requests/day — not the open "no key, experimental rate-limiting" the older note implied. Relevant if you ever integrate concept search.
  - **Exost "launched late 2025"** and **"Azed ended after No. 2,776"** are imprecise (the Azed *competition* ended with 2775; the puzzle continued). Non-load-bearing.
  - **Clue Clinic Nina/Ghost-Theme quotes** were mis-worded in the source synthesis, though the underlying *concept* (solvable-without-noticing) is sound.
  - **ingrid JPZ export** could not be confirmed from the cited FAQ — treat as unverified.
- **Thinner areas (medium confidence):** the exact **stock-grid libraries** each UK newspaper forces (sourced indirectly); the precise **barred Ximenean per-length unch table** (corroborated but partly from secondary summaries of Ximenes); **Qat's exact operator glyphs** (quinapalus.com blocks automated fetch); and **publisher submission formats** (.ccw/.cfp/NYT PDF, from a single secondary guide). Verify these against primary sources before building features that depend on their specifics.
- **Genuinely unexplored:** whether crosswordsmith's MRV/greedy engine can be *retargeted* to barred layouts (a separate engine), and a head-to-head **fill-quality benchmark** of crosswordsmith vs Exet's beam search / Qxw / scored-greedy on real cryptic grids — the load-bearing question for whether the fill engine is competitive.

---

## Sources

**Setter methodology & conventions**
- Stella Zawistowski, *Constructing a Cryptic* (grid design + fill): https://inteltainment.org/2022/04/07/constructing-a-cryptic-grid-design/ and https://inteltainment.org/2022/04/15/constructing-a-cryptic-3-grid-filling/
- *Actually Setting* (Anax): http://crypticcrosswords.net/crosswords/actually-setting/
- Michael Callaghan, *A brief guide to the construction of cryptic crossword clues* (PDF): http://crypticcrosswords.net/wp-content/uploads/2018/06/Cryptic-Crossword-Clues_v1-2.1.pdf
- Cryptic crossword — Wikipedia: https://en.wikipedia.org/wiki/Cryptic_crossword
- Crossword (American vs British) — Wikipedia: https://en.wikipedia.org/wiki/Crossword
- Crossword Unclued (grid symmetry / checking / barred): https://www.crosswordunclued.com/2009/09/crossword-grid-checking.html
- Ninas — Alberich: https://www.alberich-crosswords.com/articles/ninas
- Ximenean clueing revisited — Alberich: https://www.alberich-crosswords.com/articles/ximenean-clueing-revisited
- Ximenes vs Araucaria — David Astle: https://davidastle.com/da-blog/ximenes-vs-araucaria
- Times Crossword House Style: https://times-xwd-times.livejournal.com/174088.html

**Tools**
- Exet: https://github.com/viresh-ratnakar/exet and https://viresh-ratnakar.github.io/about-exet.html
- Exolve (format spec): https://github.com/viresh-ratnakar/exolve
- Crossword Compiler: https://www.crossword-compiler.com/ (Pro Grid Filler: /ProFill.html; word scoring: /en/help/html/aboutwordlistscoring.htm)
- Qxw guide: https://quinapalus.com/qxw-guide-20200708.pdf and https://www.quinapalus.com/qxw.html
- Crossword Man / Sympathy / TEA (withdrawn): https://www.rossberesford.com/crosswordman
- Crosshare: https://crosshare.org/construct
- MyCrossword: https://www.mycrossword.co.uk/crosswords/user-guide
- ingrid_core (Rust CSP filler): https://github.com/rf-/ingrid_core
- Nutrimatic: https://nutrimatic.org/2024/usage.html | Qat: https://www.quinapalus.com/qat.html | OneLook/Datamuse: https://www.datamuse.com/api/

**Data, formats & corpora**
- ipuz spec: https://www.puzzazz.com/ipuz/v1 and http://www.ipuz.org/
- libipuz extensions (barred kind): https://www.libipuz.org/misc/ipuz-extensions.html
- .puz format: http://justsolve.archiveteam.org/wiki/PUZ_(crossword_puzzles)
- kotwords (converter): https://github.com/jpd236/kotwords | puzpy: https://github.com/alexdej/puzpy
- cryptics.georgeho.org (clue corpus + datasheet): https://cryptics.georgeho.org/ and /datasheet
- Scored word lists: https://www.xwordinfo.com/WordList | https://github.com/Crossword-Nexus/collaborative-word-list | https://www.spreadthewordlist.com/ | https://github.com/christophsjones/crossword-wordlist
- Crossword abbreviations — Wikipedia: https://en.wikipedia.org/wiki/Crossword_abbreviations
- Daily Cryptic indicator directory: https://dailycryptic.co/tools/indicators

**Community / QA**
- Big Dave's Rookie Corner: http://crypticcrosswords.net/puzzles/rookie-corner/
- CCCWC: https://www.andlit.org.uk/cccwc/cwc_about.php | Azed slip archive: https://www.andlit.org.uk/azed/about.php
- Fifteensquared: https://www.fifteensquared.net/
- The Clue Clinic: https://clueclinic.com/index.php/resources/

---

## Addendum — Open Questions Resolved (follow-up verification pass)

A second pass took the 13 open / low-confidence items flagged above and resolved each against **primary sources**, with an independent skeptical cross-check that re-fetched citations (10 items carried a full second-agent cross-check; the remaining 3 — ipuz version, ingrid export, publisher formats — were resolved directly in a follow-up because their original agents malfunctioned). **Net result: 5 items confirmed as-stated, 8 corrected, 0 left unverified** — though a few corrected items retain a thin edge (noted below). The single most consequential correction for the roadmap is the **UKACD18 license**: it is plain **3‑clause BSD** (freely bundleable with attribution), not a custom or BSD‑2 license — so shipping it with crosswordsmith is cleanly permitted. Two other corrections change concrete build details: emit **ipuz v2** (not v1), and treat **barred-grid support as a separate engine**, not a patch.

### Corrections

| Item | Original claim | Corrected / confirmed fact | Verdict | Conf. |
|---|---|---|---|---|
| MyCrossword device tags | "~27" tags | **Exactly 25** named device tags in the user guide (enumerated: &Lit, Alternation, Anagram, Charade, Container, Cryptic definition, DBE, Deletion, Hidden, Homophone, Lift & separate, Link, Middles, Movement, Multiple definitions, Partial anagram, Reversal, Reverse engineer, Sides, Spoonerism, Sub anagram, Substitution, Tails, Theme, Tops) | corrected | high |
| Crossword Compiler symmetry | "Dual-symmetry = 90°" | CCW's own terms are **Normal (S), Dual-Symmetry (square grids only), left/right, top/bottom, Dual-mirror, maximum**. CCW docs **never** say "90°" (grep of the manual: one unrelated hit). "90°/quarter-turn" is a correct *external* gloss, not CCW wording. | corrected | high |
| **UKACD18 license / Exet lexicon** | bundled list `importance-and-words.txt`; "BSD-2 / Advanced Cryptics Dictionary license"; "completely unrestricted" | Exet's shipped lexicon is **`lufz-en-lexicon.js`** (270,372 entries; `.tsv` is a build intermediate). UKACD18 (©2009 J Ross Beresford) is **3-clause BSD** — retain notice + non-endorsement clause. **"Completely unrestricted" describes CMUdict, not UKACD18.** ✅ **Bundleable with attribution.** | corrected | high |
| **ipuz version string** | v1 vs v2 ambiguity | Current spec is **`http://ipuz.org/v2`** ("replaces v1 and v1.1… recommended that all implementors support v2"). Standard crossword kind = **`["http://ipuz.org/crossword#1"]`**. Emit v2. | corrected | high |
| Datamuse / OneLook API | key required after 2027-01-01, 100k/day cap | **Confirmed verbatim**: free + key-free until 2027-01-01; thereafter API key required, 100,000 req/day/key. No published paid tier. | confirmed | high |
| Exost launch / Azed numbering | "Exost late 2025"; "Azed ended after No. 2,776" | Exost launched **mid-Dec 2025** (repo commits Dec 18, archive Dec 19). The Azed **competition** ended after **No. 2,775** (4 Jan 2026); **No. 2,776 (1 Feb 2026) was the first plain, no-competition puzzle — the puzzle continues** (Crowther; competition to return under "Gemelo"). | corrected | high |
| Nina "solvable without noticing" | quoted as if verbatim | Concept is sound but the phrase is a **paraphrase**. Verbatim sources split: Clue Clinic — a Nina "does not need to be identified in order to solve the puzzle"; Alberich (personal view) — "possible to solve the puzzle without even noticing a Nina." | corrected | high |
| **ingrid JPZ export** | "JPZ export unconfirmed" | **`ingrid_core` (the Rust CLI — the comparator to crosswordsmith's engine) exports no puzzle file at all** — it prints a plain-text filled grid. The Ingrid GUI app's sharing path is **Crosshare / squares.io**; app-side JPZ export remains unconfirmed from primary docs. | corrected | high (core); low (app JPZ) |
| UK stock-grid libraries | "most UK papers force a grid library; Independent excepted" | **First-hand confirmed for the Times** (64 editor-managed stock grids, all setters get the set) and the **Independent** (setters design their own) via setter Anax + the Times insider blog. Industry norm well-attested. **Guardian/Telegraph/FT specifics remain anecdotal.** | confirmed | high (Times/Indy); low (others) |
| Barred / Ximenean unch table | corroborated only via secondary summaries | **Confirmed against primary Ximenes text + Listener 2011 setter notes**: 3 letters → 0 unch; 4–5 → 1; 6–7 → 1–2; 8–9 → 2–3; 10–12 → 3–4; ~75 %+ checking overall. Treat per-length figures as a **band**, not a hard cap. | confirmed | high |
| Qat operators | unconfirmed (site blocked bots) | **Full operator set recovered** from quinapalus.com/qat.html (the 403 was a bot filter; curl with a browser UA → 200). Wildcards `.` `*` `[…]` `@` `#`; anagram `/`; reversal `~`; boolean `! & |`; misprint backquote; length `n:`/`n-:`/`-n:`/`n-m:`; equation solver (A–Z vars, `=`, `!=`, `;`); qategories `{def:…}`. **No `&lit` primitive exists.** | confirmed | high |
| Publisher submission formats | ".ccw/.cfp/NYT-PDF" from one secondary guide | **No universal standard.** US indies prefer **`.puz`** (AVCX: "preferably in .puz format"); AVCX Cryptic has a separate spec. UK papers (Times/FT/Independent) request **Crossword Compiler files**; submission is editor-mediated, not publicly format-specced. **The New Yorker discontinued cryptics in March 2024.** | corrected | medium |
| **Barred-grid engine feasibility** *(code analysis)* | (previously unexplored) | **Separate engine / new cell model — not a localized change.** The grid is block-only (`cell → empty\|letter`, `init_grid`/`new_tile`, crossword.pl:951); there is **no bar/edge datum** between filled cells. `adj_is_free/4` (682–710) actively **bans** the flush-parallel adjacency that barred grids require; `no_word_merge`/`check_prev_cell`/`check_next_cell`/`find_intersecting_word` all assume runs terminate only at an empty cell or edge. Supporting bars means a per-cell-edge bar map + rewriting the legality core. | confirmed | high |

### Still thin / unverified

- **Guardian, Telegraph, FT grid policies** — the "forced stock-grid library" norm is first-hand only for the Times (and the Independent's custom-grid exception). The Guardian's reuse of ~72 grids is *consistent with* a mandate but never stated as one; Telegraph/FT have no first-hand source. *To close: a setter statement or editor guidance per paper.*
- **Ingrid (GUI app) JPZ export** — a secondary snippet suggested "export as JPZ from the File menu," but the primary FAQ only documents Crosshare/squares.io upload. *Immaterial to crosswordsmith* (the relevant comparator is `ingrid_core`, which exports no file). *To close: the app's own docs/changelog.*
- **UK newspaper submission file specs** — beyond "they use Crossword Compiler files," no paper publishes an exact format/mechanism; submission is editor-mediated. *To close: direct contributor guidance, which papers generally don't post publicly.*
- Minor: the **6–7 letter unch row** differs between Ximenes ("1 or 2") and the Listener summary ("2") — use a band; and UKACD18's exact distributed artifact on quinapalus.com couldn't be fetched (403), though its BSD-3 license text is verbatim-confirmed from two independent repos.

### Implications for the roadmap

- ✅ **Bundle UKACD18.** Cleared: it is BSD-3, redistributable with the copyright notice + non-endorsement clause. This unblocks the "license-clean scored/fill wordlist" direction without paywalled IP. (Pair with CMUdict — 2-clause BSD — for pronunciations if homophone tooling is ever wanted.)
- ✅ **ipuz emitter targets v2.** Emit `"version":"http://ipuz.org/v2"` and `"kind":["http://ipuz.org/crossword#1"]`. (v1 is still parsed by many consumers, but v2 is the recommended current target.)
- ✅ **Barred grids = separate engine, schedule accordingly.** Do **not** plan barred support as an option on the existing engine; it needs a new cell/edge model and a rewrite of the placement-legality core (`adj_is_free`, `no_word_merge`, `check_prev/next_cell`, `find_intersecting_word`). The "barred-checking-math lint" remains a cheap, in-scope add even without a barred *generator*.
- **Grid-linter unch rules are now firmly sourced.** Use the Ximenes/Listener per-length band for a *barred/Ximenean* profile, and the "≥ half checked, ≤2 consecutive unch, double-unch not at word ends" set for the *blocked-UK* profile — both now primary-sourced.
- **Export realism:** `.puz` is the lingua franca for US indie submission and `.ipuz`/Exolve for open interchange; UK paper submission stays Crossword Compiler-centric and editor-mediated. Reaffirms: emit `.ipuz`/Exolve natively, reach `.puz` via a converter, and don't over-invest in publisher-specific submission packaging.
- **Pure-fact corrections** (MyCrossword 25 tags; CCW "Dual-Symmetry" not "90°"; Nina quote attribution; Exost/Azed dates; New Yorker cryptic end) tighten the main report but don't move features.

### Sources (follow-up pass)

- MyCrossword User Guide — https://www.mycrossword.co.uk/crosswords/user-guide
- Crossword Compiler Help — Grid symmetry — https://www.crossword-compiler.com/en/help/html/symmetry.htm ; PDF manual — https://cdn.crossword-compiler.com/cc/help/Crossword%20Compiler%20Help.pdf
- Exet README (lexicon, UKACD18 BSD-3, Exost) — https://raw.githubusercontent.com/viresh-ratnakar/exet/master/README.md ; lufz README — https://raw.githubusercontent.com/viresh-ratnakar/lufz/master/README.md ; UKACD18 license — https://github.com/pbevin/lexicon/blob/master/LICENSE.ukacd18
- ipuz spec (v2 + crossword kind) — https://libipuz.org/ipuz-spec.html ; http://www.ipuz.org/
- Datamuse API — https://www.datamuse.com/api/
- Azed No. 2776 (Fifteensquared) — https://fifteensquared.net/ ; Clue Clinic — https://clueclinic.com/index.php/the-setting-room/ ; Alberich, Ninas — https://www.alberich-crosswords.com/articles/ninas
- ingrid_core — https://github.com/rf-/ingrid_core ; Ingrid FAQ — https://ingrid.cx/faq/
- Anax, "Actually Setting" — http://crypticcrosswords.net/crosswords/actually-setting/ ; Times grids — https://timesforthetimes.co.uk/more-than-you-ever-want-to-know-about-grids
- Ximenes Ch. 11 — https://xotaotc.nfshost.com/chapter-11-composing-a-ximenes-bars/ ; Listener Notes for Setters (2011) — https://www.listenercrossword.com/PDF/notes2011web.pdf
- Qat syntax — https://www.quinapalus.com/qat.html
- AVCX submissions — https://avxwords.com/submit-a-puzzle/
- crosswordsmith code — `crossword.pl` (init_grid 951, adj_is_free 682, no_word_merge 609, check_prev/next_cell 632); `quality.pl`

---

## Addendum — Symmetry: how universal is the constraint? (discovery pass)

A targeted discovery pass took the single question the main report leaned on hardest — *is 180° rotational symmetry actually a hard universal rule?* — and pressure-tested it against primary and practitioner sources (17 sources fetched, 69 candidate claims, 25 adversarially verified with a 2-of-3-refute kill threshold, 21 confirmed / 4 killed). **Net finding: the main report's "hard, codifiable rule" framing overstates it. 180° rotational symmetry is a *strong default with theme-gated exceptions*, and its strictness varies sharply by grid type and publication. For the engine, it should be a per-publication-profile setting with a per-puzzle override — not a hard constraint.**

### Verdict

| Modelling option | Holds? |
|---|---|
| HARD universal constraint (fail anything asymmetric) | ❌ No — overstates real practice and tooling |
| STRONG default with exceptions | ✅ Yes — accurate for blocked dailies |
| **Per-publication-profile setting (default 180° for blocked dailies; selectable 90°/axial/diagonal/none; per-puzzle override for theme/Nina/barred)** | ✅ **Recommended** |

### The gradient (strictest → loosest)

- **Blocked daily cryptics (Times, Guardian, Telegraph, FT, Independent).** 180° rotational is the standard and *editorially near-mandatory in practice*. Strongest primary statement, from working setter **Phi (Paul Henderson** — Independent/FT, "Kcit" Toughie at the Telegraph, a Listener setter): *"I think it would be hard to get a truly asymmetric grid into an outlet for blocked daily puzzles."* **Confidence: high** (but see limitations — this is a setter's blog, not a newspaper house-style page).
- **Barred thematics (Azed, Listener, Mephisto, …).** The constraint visibly relaxes — same source: *"They're commoner in barred thematics."* Asymmetric grids are a recognized category, though symmetry remains the strong default even here. **Confidence: high for the divergence; low for any per-publication specifics.**
- **Mainstream blocked dailies permit theme-gated exceptions** (documented, not folklore): **NYT** — *"yet rarer are asymmetrical puzzles, usually when an unusual theme requires breaking the symmetry rule"*; **LA Times** — rotational preferred, left/right "acceptable if the theme requires it," and "no symmetry or unusual symmetry will only be considered if the pattern is theme-related"; **Daily Beast** — published its *first-ever* asymmetric grid because theme-entry lengths "required bending the symmetry rules." **USA Today** has been *normalizing* asymmetry since ~2021 (even mainstream daily practice can shift). **Confidence: high** (caveat: these are US blocked grids; NYT/LAT are non-cryptic).
- **Other symmetries are real and supported alternatives** — 90° quarter-turn, left-right/up-down axial mirror, diagonal — appearing in themed/variety and diagramless puzzles, and offered as first-class options in construction tools. **Confidence: high.**

### What the tools and specs encode (strongest signal for an engine builder)

Every tool models symmetry as **default-on-but-overridable — never an immutable invariant**:

| Tool / spec | Symmetry model |
|---|---|
| **Qxw** | Defaults to 180° rotational; offers 90°/mirror/axial "combined at will"; explicit **`Symmetry-None`** menu option. *"does not particularly distinguish between barred and blocked grids."* |
| **Exet** (Exolve author's tool) | Auto-applies symmetric edits by default; **`Allow asymmetry`** per-puzzle checkbox that **resets to enforcing symmetry for each new crossword** — the cleanest model of exactly the recommended strong-default + per-puzzle-override behaviour. |
| **Crossword Compiler** | Six symmetry generators: Normal (rotational), Dual-Symmetry, Left/right, Top/bottom, Dual-mirror, Maximum. *(Its labelled default being the rotational one was **not** confirmed — refuted 1-2; don't assume it.)* |
| **ipuz spec (v2)** | **No symmetry field at all** — a fully asymmetric grid is valid ipuz. If crosswordsmith emits ipuz, the format will neither require nor preserve a symmetry declaration. |
| **libipuz** | Symmetry is on-demand enforce/detect/check API helpers, not a stored invariant or validity gate. |

### The deeper "why"

Symmetry is an **aesthetic/craft + editorial-gatekeeping convention, not a technical necessity of the fill**. Phi: *"It's perfectly possible — though I'd say harder — to meet your length and checking criteria in an asymmetric grid."* Crossword Unclued: symmetry *"doesn't affect the solving in any way, just adds to visual appeal."* The convention traces to the Amateur Crossword Puzzle League of America (1924) — partly a hot-metal-typesetting-era legacy. **Implication: crosswordsmith's fill engine never needs symmetry to produce valid fill; symmetry is a *validation/profile* concern, not a solver concern** — consistent with design principle #2 (keep the solver metadata-agnostic).

### Limitations of this pass (honest)

- **No British newspaper *house-style page* was obtained.** The strongest British primary source is one setter's blog (Phi). The hard documentary exceptions (NYT/LAT/Daily Beast) are **US blocked grids, and NYT/LAT are non-cryptic** — solid for the "why/enforcement" question, weaker as British-cryptic house style.
- **Per-barred-publication detail is thin.** No confirmed source pins down *which* symmetry Azed vs Listener vs Mephisto each use; a claim attempting to characterize barred mixed (vertical-rotational / horizontal-mirror) symmetry was **refuted 0-3**.
- The absolutist phrasing *"all mainstream grids are 180°"* was **killed 0-3** in verification; the confirmed framing is "tradition / typical / nearly all" — i.e. a strong default.
- Tooling behaviour reflects 2020–2025 docs and can change.

### Implication for the roadmap

- **Downgrade symmetry from "biggest hard gap" to a *profile-able default*.** This folds into the existing **grid-linter** (must-have #1) and **house-style / grid profiles** (nice-to-have #10) rather than demanding a symmetry-aware *generator*.
- **Linter behaviour:** report a **symmetry deficit**, with severity set by profile — `BLOCKED-UK` → enforced 180° (warn/fail), `BARRED` → relaxed, plus a per-puzzle **`allow-asymmetry`** override (mirroring Exet's toggle) for theme/Nina work. Don't hard-FAIL asymmetry unconditionally.
- **Export:** emitting ipuz neither requires nor carries a symmetry flag, so symmetry stays purely a crosswordsmith-side lint/profile concept.

### Open questions (carried forward)

1. Do the major British blocked dailies have *written* house-style statements on 180° symmetry, and do any permit theme/Nina exceptions the way NYT/LAT do — or is it an unwritten near-absolute rule enforced editorially?
2. What symmetry does each specific *barred* publication actually use/require (Azed, Listener, Mephisto, Spectator, Inquisitor, EV), and how is it notated to solvers (e.g. the Listener's per-puzzle symmetry codes)?
3. In US *cryptic* venues specifically (AVCX/Out of Left Field, the old New Yorker cryptics, The Nation, Harper's), what symmetry conventions apply — British barred-tolerance or US blocked-grid rotational?
4. Base-rate question: how *frequently* are theme/Nina-driven asymmetric grids actually accepted and published in mainstream blocked outlets — theoretically permitted but vanishingly rare, or genuinely seen?

### Sources (symmetry discovery pass)

- Phi (Paul Henderson), *Symmetry* (setter's own site, primary) — http://phionline.net.nz/56-2/symmetry/
- Crossword Unclued, *Crossword Grid Symmetry* and *Barred-Grid Crosswords* — https://www.crosswordunclued.com/2009/09/crossword-grid-symmetry.html ; https://www.crosswordunclued.com/2009/03/barred-grid-crosswords.html
- Wikipedia — *Crossword*, *Cryptic crossword*, *The New York Times crossword* — https://en.wikipedia.org/wiki/Crossword ; https://en.wikipedia.org/wiki/Cryptic_crossword ; https://en.wikipedia.org/wiki/The_New_York_Times_crossword
- LA Times grid guidelines (via CommuniCrossings) — https://communicrossings.com/constructing-crosswords-grid
- Daily Beast, *For the Good of the Puzzle, Sometimes Rules Must Be Broken* (Matt Gaffney) — https://www.thedailybeast.com/for-the-good-of-the-puzzle-sometimes-rules-must-be-broken/
- Grid Alchemy, *How USA Today is changing crossword symmetry* — https://grid-alchemy.com/2021/11/20/how-usa-today-is-changing-crossword-symmetry/
- The Straight Dope, *Why are crossword puzzles symmetrical?* (1924 ACPLA origin) — https://www.straightdope.com/21341628/why-are-crossword-puzzles-symmetrical
- Crossword Compiler Help — Grid symmetry — https://www.crossword-compiler.com/en/help/html/symmetry.htm
- Exet README (Allow-asymmetry toggle) — https://github.com/viresh-ratnakar/exet
- Qxw guide (Symmetry-None; combinable symmetries) — https://quinapalus.com/qxw-guide-20200708.pdf
- ipuz spec (no symmetry field) — https://libipuz.org/ipuz-spec.html ; libipuz Crossword class — https://libipuz.org/libipuz-1.0/class.Crossword.html

---

## Addendum — The closed-set vs. authentic-layout tension (structural analysis)

The symmetry discovery above exposes a more fundamental point that the main report circles but never states bluntly: **crosswordsmith's founding framing — "given a fixed set of words, lay them out on a crossword grid" — is in structural (near-pathological) tension with the goal of producing an *authentic* cryptic layout.** This section makes the argument explicit, grounds it in the actual engine code, and lays out the resulting architectural fork. It is analysis, not new external research.

### Why the ecosystem doesn't have this problem — and crosswordsmith does

Every other tool (Exet, Crossword Compiler, Qxw) and every human setter **decouples** two problems:

1. **Grid legality** — symmetry, ≥50% checking, min-length, no-triple-unch — is a property of the **black-square pattern alone**. It is satisfied at *design time, before any word exists*, by choosing a legal template (or a stock grid). It is entirely word-independent.
2. **Fill** — finding words for the slots — is a search over an **open dictionary** (~270k entries), so a letter-consistent fill almost always exists.

crosswordsmith makes **two inversions at once**, and crucially *both* engines do so — confirmed against the code:

- The grid is **emergent, not designed.** `crossword.pl`'s `init_grid(GridLen, …)` builds a **blank square canvas of a fixed *size*** (a single integer dimension), *not* a pre-laid black-square template; words are placed by interlock and the "black squares" are simply wherever no word landed. `quality.pl` goes further and **chooses the canvas size itself** (`grid_candidates/2`), greedily packing for density. *(Note: the main report calls `crossword.pl` a "fixed-grid" solver — more precisely it is fixed-**size**, blank-canvas; neither engine consumes a legal grid template.)*
- The dictionary is **closed** — both engines only ever place the supplied `Words`; there is no filler lexicon and no candidate-from-dictionary lookup anywhere in the code.

The consequence: the grid's legality is now a **downstream consequence of how your specific words happen to interlock at shared letters.** The constraints the ecosystem gets *for free at design time* become emergent properties of an arbitrary letter-overlap graph — which is why they cannot be dialled in.

### The constraints, ranked by how badly they fight a closed set

| Constraint | Severity vs. closed-set place-all | Why |
|---|---|---|
| **180° rotational symmetry** | **Pathological** (≈ infeasible) | See length-pairing argument below — needs an exact, dictionary-free set-cover. |
| **~50% checking** | Severe; empirically ~4–5× short | Spec reports real inputs at **0.09–0.14 checked fraction** vs. the ~0.5 ideal. Dense interlock needs a rich pool a small fixed set can't supply. |
| **No triple-unch / double-unch position** | Fails *downstream* of checking | These only bite once checking is decent; sparse checking ⇒ long unchecked runs by default. |
| **Min length ≥3, connectivity, no-merge** | Manageable | Input-dependent (short labels); connectivity already enforced by construction. |

**The symmetry wall, precisely (the crux).** Under 180° rotation, every Across slot maps to *another Across slot of equal length* (its rotated image), except a possible central self-paired slot — likewise for Down. So a fully-symmetric, fully-filled blocked grid requires your **multiset of word lengths to pair up** to match some legal template's slot inventory, *and* every crossing to be letter-consistent, *with no filler dictionary to patch a missing length*. That is an exact set-cover with crossing constraints over a tiny, uncurated dictionary — almost never satisfiable for an arbitrary input (a site's nav labels, a setter's seed list). "You'd have to be lucky" is exactly right: the probability is ≈0, and search cannot rescue an infeasible instance. The measured 0.09–0.14 checked fraction is the empirical fingerprint of the same under-supply.

**Blunt conclusion:** you cannot have *"place exactly these words"* **and** *"authentic cryptic layout (symmetry + dense checking)"* in one engine. The current `--quality` "British lattice aesthetic" is the *relaxed* version of authenticity, and 0.09–0.14 is what relaxation buys.

### The resulting fork

- **A — Two products.** Keep crosswordsmith as the "pack *these specific* words densely" tool (its genuine niche — no Exet/CCW equivalent, per the main report) and stop benchmarking it against cryptic authenticity. A crossword-*shaped* table-of-contents is a different artifact from a publishable grid; lean in rather than apologise.
- **B — Go grid-first + dictionary; demote the word set to *seeds*.** The *only* path to authentic layout, and how every other tool works. See the architecture sketch below.
- **C — Honest hybrid (≈ today's `--quality`).** Closed set as seeds, dropping allowed, *relaxed* profile with symmetry/checking as **advisory** lints. Label it "crossword-shaped," not "publishable cryptic." Cheapest; matches current behaviour.

What to *avoid*: any roadmap item that implicitly assumes both — e.g. chasing ~0.5 checking, or *generating* symmetric grids, out of the closed-set packer. That is the one direction the math forbids.

### Option B — what grid-first would actually take (architecture sketch)

Grounded in the code, Option B is a **new engine**, not a modification of the existing legality core — because nothing today consumes a legal grid or a filler lexicon. It would need:

1. **A legal-grid template representation + library.** A black-square pattern (or explicit slot list) for pre-validated symmetric grids — the main report's "stock grid libraries." Validate symmetry/checking/min-len/unch *once, at design time*, reusing the existing metric predicates (`word_meets_half/2`, `max_unch_run`) as **template validators** rather than per-word floors. *(How big must this library be? Real publications curate only **dozens** — the Times 64, the Guardian ~72 — so this is a small, static, **bundle-don't-enumerate** asset, not a generation problem. The full legal grid space is explosively large (≥10⁸ at 15×15) and not enumerable as a list, but you never need it. See "Addendum — Stock-grid libraries & the size of the legal-grid space.")*
2. **Slot enumeration from a fixed pattern.** Derive Across/Down runs from the template's black squares. (Today, slots *emerge* from placement; there is no enumerate-slots-from-a-pattern predicate.)
3. **A filler dictionary + pattern index.** Load a large lexicon, indexed by length and by partial pattern (the letters forced by crossings). **The report already cleared UKACD18 (BSD-3) for bundling** — so this ingredient is license-ready.
4. **A dictionary fill search.** MRV/backtracking selecting dictionary words per slot subject to crossing constraints, with the user's words **pinned as seeds**. The MRV machinery in `crossword.pl` is conceptually reusable, but it must search the *dictionary* per slot, not a fixed `Words` list.
5. **Drop the free-canvas legality core.** `adj_is_free/4` (bans flush-parallel adjacency), `no_word_merge`, and `check_prev/next_cell` exist *because* the canvas is blank and runs must self-delimit. In a template grid the black squares delimit runs explicitly, so that whole core is replaced by "fill the enumerated slots."

**Reused:** JSON/ipuz emit, clue numbering, and the checking/unch metric predicates (repurposed as design-time validators). **New:** template format + library, slot enumeration, filler dictionary + index, and a dictionary-driven fill search. Both *new data* ingredients (grid templates, UKACD18) are already identified in the main report as available and license-clean — so Option B is gated on engine work, not on data acquisition.

**Recommendation:** the main report's "differentiator" framing already leans toward **A/C** — keep the closed-set packer honest and relaxed, and treat symmetry/checking as advisory, profile-able lints (consistent with the symmetry addendum and design principle #5). Pursue **B only if authentic, publishable cryptic output becomes a first-class goal** — and price it as a second engine sharing only the emit/metric layers, not a patch to the existing one.

---

## Addendum — Stock-grid libraries & the size of the legal-grid space (discovery pass)

This pass answers two quantitative questions raised while scoping a "select from pre-validated stock grids" feature: **(A)** how many grids do real libraries actually offer, and **(B)** once the cryptic constraints are applied, is the space of *legal* grids practically enumerable or combinatorially explosive? Method: throttled multi-source web search + 3-vote adversarial verification (14 sources, 56 candidate claims, 14 verified → 13 confirmed / 1 killed). **Headline: real libraries hold *dozens* of grids; the legal space is *explosively large* and not enumerable as a list — and those two facts sit ~6 orders of magnitude apart, which is the whole point.**

### Part A — How many grids do real libraries offer? Dozens.

| Source | Count | Provenance / confidence |
|---|---|---|
| **The Times** | **64** stock grids "in current use" | Insider data — setter Anax (Dean Mayer) via Sunday Times puzzles editor Peter Biddlecombe, *dated April 2008*. Built up over decades: editor Akenhead's original 25 → culled/extended by Greer & Laws (only 13 of Akenhead's 25 survive). Internal breakdown given (49 odd-row/odd-col, 8 odd/even, 7 double-unch; 17 quarter-turn-symmetric; 25 grids with 30 answers) — too specific to be folklore. **High** (but an 18-yr-old snapshot — see caveats). |
| **The Guardian** | **~72** distinct grids across ~1,938 puzzles | Independent enumeration by Steven Tattersall (puzzles 24,000–25,938); most-used grid appeared 108×. **High** (one analyst's measurement, not an official statement). |
| **The Independent** | setters design their own | The documented exception — accommodated by "most publications." **High.** |
| **Crossword Compiler** | fixed "Grid Library" / "Grid Library Manager" **+ algorithmic Grid Pattern Generator** | Vendor docs: the generator produces new patterns on demand from parameters (symmetry type, min word length, min/max block %, max consecutive edge blocks, min words per length), supporting American (fully checked) vs cryptic/quick (alternate-unchecked) styles. So software treats the space as *parameterized/generated*, not a fixed enumerated list. Exact size of the *fixed* library: not found (open). **High** (primary vendor docs). |

**Verdict (A):** a setter at a given publication picks from roughly **50–80** pre-validated grids. The order of magnitude is firmly **dozens** — not hundreds or thousands. A "bundle a stock-grid library" feature is therefore modest and realistic, not a moonshot.

### Part B — How big is the *legal* grid space? Explosively large; not enumerable as a list.

The only published enumeration of *British/cryptic* legal grids is **George Ho, "Counting Cryptics"**, which enforced exactly this report's Part-B constraint set (180° rotational symmetry, odd length, full white-cell connectivity, ≥50% checking rounded up, no >2 consecutive unches, no double-unch at word start/end):

| Grid size | Cryptic legal grids (Ho) | American legal grids (contrast) |
|---|---|---|
| 5×5 | 17 | 12 |
| 7×7 | 346 | 312 |
| 9×9 | **9,381** | 31,187 |
| 11×11 | *intractable — not computed* | 17,438,702 (≈1.7×10⁷) |
| 13×13 | — | 40,575,832,476 (≈4.06×10¹⁰) |
| 15×15 | — *(no published count exists)* | **404,139,015,237,875 (≈4.04×10¹⁴)** |

American figures: Michael Kleber / OEIS **A323838** (the cited 15×15 value is the "no black edge rows/columns" convention; Kleber's all-grids total is 409,764,131,469,787). Independently verified against Kleber's own post and a DP implementation that matches A323838 through 15×15.

Salient points:
- **Ho's enumeration died at 11×11** ("computationally intractable to enumerate all possible grids… we can only count them"). So there is **no published 15×15 cryptic count** — only the American ruleset reaches 15×15.
- **Extrapolation (this report's estimate, *not* sourced):** Ho's per-step growth ratio climbs (17→346 is ×20.4; 346→9,381 is ×27.1). Continuing at ~30–36× per size step puts 15×15 cryptic *very roughly* at **10⁸–10⁹** — about 5–6 orders of magnitude *below* the American 4×10¹⁴ (the stricter checking/unch rules prune hard), yet still **astronomically beyond any browsable/curatable list.** Treat 10⁸–10⁹ as an order-of-magnitude guess, not a computed figure.
- **Crossover quirk:** at *small* sizes cryptic has *more* legal grids than American (17>12, 346>312); American only overtakes at 9×9. The stricter cryptic rules dominate only as the grid grows.
- **⚠️ Bogus number to avoid:** a figure of **1.7×10¹⁹** for A323838 at 15×15 circulates and was **refuted 3-0** in verification (an off-by-one in grid size, conflating term offsets). The correct American 15×15 count is ≈4.04×10¹⁴. Re-verify A323838 term values directly against OEIS before quoting any of them.

**Verdict (B):** the legal 15×15 cryptic grid space is **explosively large and not practically enumerable as a list** (intractable past 9×9 in the one published attempt). *Counting* it may be tractable with a good dynamic program; *listing/browsing* it is not.

### The key distinction (and the direct answer to "enumerable vs. explosive")

These are **two spaces separated by ~6+ orders of magnitude:**

- **(i) Legal space** — ≥10⁸ at 15×15; not enumerable as a list.
- **(ii) Conventionally-accepted subset** — **dozens** (Times 64, Guardian ~72); tiny and fully curatable.

Publications reuse a fixed handful **not because the legal space is small** but because **(a)** most legal grids are not *aesthetically/conventionally* acceptable, and **(b)** historically new grids were expensive (hot-metal typesetting legacy — Wikipedia, verbatim: "In the past this was because hot metal typesetting meant that new grids were expensive"). The gap between "what's legal" and "what's used" is the entire reason stock libraries exist.

### Implications for crosswordsmith

- **The stock-grid feature is cheap and realistic.** You **bundle ~50–100 curated templates** — exactly what the Times and Guardian themselves do — rather than enumerating or generating the legal space. A small, static, license-checkable asset. **Bundle, don't enumerate.**
- **Or follow Crossword Compiler:** ship a modest fixed library *plus* a constraint-parameterized generator (symmetry, min-len, block %) for variety. Both are tractable; neither needs the legal space enumerated.
- **It sharpens Option B (grid-first + dictionary, from the previous addendum):** template *supply* was never the bottleneck — dozens suffice and are trivially curatable. The real Option-B work is the **closed-set/dictionary fill** of whatever template you pick. Grid templates ≈ free; the fill ≈ the hard part.
- **For the linter/validator (must-have #1):** Ho's exact constraint list and the Times house-style rules (≥50% checking; never >2 unches in succession; double-unches never the first or last two letters) are a ready-made, primary-sourced acceptance spec for a `BLOCKED-UK` profile.

### Caveats

- **Time-sensitivity:** the Times "64" is an April-2008 snapshot; the set has evolved and editorship changed (Mike Laws died 2011), so the *current* count may differ.
- **Single-source enumeration:** the cryptic counts (17 / 346 / 9,381) come from one personal blog (George Ho), not peer-reviewed, using a deliberately *simplified* rule set (omits clue-count and black-island-size limits) — so they are specific to his ruleset. Rated **medium**.
- **No 15×15 cryptic count exists** anywhere; the 10⁸–10⁹ figure here is this report's extrapolation only.
- **Guardian ~72** is one analyst's measurement over a fixed puzzle range, not an official figure.

### Open questions

1. The Times' *current* (2026) stock-grid count, given conventions/editorship changed since 2008's 64.
2. Counts for the other named publications — Telegraph, FT, Daily Mail (and the Independent's own-grid practice) — no published figures confirmed.
3. Exactly how many patterns ship in Crossword Compiler's *fixed* Grid Library (UK cryptic vs American), and whether Sympathy/Qxw/Exet/Crosshare/ingrid ship fixed libraries and of what size.
4. The true legal 15×15 *cryptic* count with an optimized algorithm and a complete (non-simplified) British ruleset — merely large (10⁸–10⁹) or genuinely American-scale (10¹⁴)?

### Sources (grid-library & grid-space pass)

- Times "64 grids" — Peter Biddlecombe / Anax, *More than you ever want to know about grids* — https://timesforthetimes.co.uk/more-than-you-ever-want-to-know-about-grids ; https://times-xwd-times.livejournal.com/227851.html
- Times Crossword House Style (checking / unch rules) — https://times-xwd-times.livejournal.com/174088.html
- Akenhead grid history — https://crosswordsakenhead.com/crossword-history
- Guardian ~72 grids (Steven Tattersall enumeration) — https://www.clarets.org/steve/crosswords/2013_05_guardian_grids_since_2007.html
- *Cryptic crossword* (stock grids; hot-metal rationale) — https://en.wikipedia.org/wiki/Cryptic_crossword
- George Ho, *Counting Cryptics* (cryptic enumeration 5×5 / 7×7 / 9×9) — https://www.georgeho.org/counting-cryptics/
- OEIS **A323838** (American legal-grid enumeration, Michael Kleber) — https://oeis.org/A323838 ; Kleber post (404T figure) — https://x.com/Log3overLog2/status/1093152564911181824 ; DP matching A323838 — https://github.com/akshayravikumar/crosswords
- Crossword Compiler — Grid Pattern Generator — https://www.crossword-compiler.com/en/help/html/GridPatternGenerator.htm ; newspaper grids — https://www.crossword-compiler.com/en/newspaper.html
- Anax, *Actually Setting* (stock grids) — http://crypticcrosswords.net/crosswords/actually-setting/ ; Crossword Unclued, grid checking — https://www.crosswordunclued.com/2009/09/crossword-grid-checking.html

---

## Addendum — Grid sizes in the wild (discovery pass)

A focused single-agent pass (sequential web search, kept deliberately light for memory) mapping the grid **dimensions** actually used across published cryptics — not just "15×15 is standard" but the full range, the publication/puzzle-type at each size, and blocked-vs-barred. **Headline: the distribution is strongly *bimodal* — 15×15 blocked (mainstream dailies) and ~12×12 barred (advanced/variety) are the twin standards, with thin tails running from 8×10 up to 24×24.** This bears directly on grid-size defaults, stock-grid presets, and what is out of scope for the current block-only engine.

### The two standards

| Cluster | Size | Block/bar | Who |
|---|---|---|---|
| **① Mainstream daily** | **15×15** | blocked | Times, Guardian, Telegraph (+ Toughie), FT, Independent, Observer *Everyman*; US blocked cryptics (The Nation, Out of Left Field). Guardian *Quiptic* too (easier, but full-size). |
| **② Advanced / "variety"** | **~12×12** (flexes to 13×13) | barred | Mephisto, Azed, Listener, Inquisitor, Enigmatic Variations, Spectator, Magpie; US "variety" cryptics (Hex / Atlantic + WSJ, Harper's, AVCX). |

Wikipedia, verbatim: standard cryptics are "generally 15×15, with half-turn rotational symmetry"; variety/advanced cryptics "typically use a barred grid with no black squares and a slightly smaller size; 12×12 is typical." Barred grids accept *even* dimensions (no odd-forcing block-symmetry constraint), hence the 12×12 / 13×13 flex.

### The tails / outliers

| Size | Where | Block/bar | Confidence |
|---|---|---|---|
| **13×13** | Times **Quick Cryptic** (the one daily that shrinks the grid) | blocked | high |
| **23×23** | Times Saturday **Jumbo** — largest well-attested cryptic | blocked | high |
| **24×24** | Globe and Mail Canada Day puzzle (discontinued after 2015) | barred | high |
| **8×10** | **New Yorker** cryptic (Tina Brown era 1997–99 + 2021–Mar 2024 revival) — the small end | barred | high |
| **variable** (rectangular, circular, irregular — e.g. a "double-L") | Listener / thematic puzzles where the grid *is* the theme | barred | high (variability attested by a working setter) |
| **3-D** (dimensions not cited) | Sirius 3D Calendar Puzzles | — | format real, dimensions unverified |

**Verified observed range: ~8×10 to 24×24** — but it is two tight clusters with thin tails, not a smooth spread.

### Two corrections to common assumptions
- **The Guardian Quiptic is 15×15, not a small grid** — "quick" means *easier*, not *smaller*. Among the regulars, only the **Times Quick Cryptic** (13×13) actually shrinks the grid.
- **US cryptics split in two:** *blocked British-style 15×15* (The Nation; Out of Left Field) vs *barred "variety"* (Cox & Rathvon "Hex" in the Atlantic/WSJ; Harper's; AVCX; the New Yorker's 8×10). Crucially, American cryptics adopt the **British** checking convention (≈half-checked), not the fully-checked American style — so "American cryptic" does **not** mean an American-style grid.

### Confidence flags
- **Solid (specific sources, often verbatim):** 15×15 standard blocked; 12×12 typical barred; Times Quick 13×13; Times Jumbo 23×23; New Yorker 8×10 barred; Mephisto 12×12; Globe and Mail 24×24 barred; Independent 15×15.
- **Weaker / inferred:** Azed's exact square size (confirmed barred, number not pinned); per-title sizes for Listener / Inquisitor / EV / Spectator / Magpie (confirmed barred + ~12×12 family, but from the general barred standard, not title-specific cites — and they vary by theme); FT/Telegraph 15×15 (solid by convention, not per-paper-cited); **Daily Mail 15×15 ≈ folklore** (no dedicated source); exact rows×cols for Hex / Harper's / AVCX variety grids (vary puzzle-to-puzzle); Guardian *Genius* 15×15 (inferred from Guardian blocked house style).
- **Not found:** any verified **27×27** cryptic — largest confirmed is the 23×23 Jumbo; treat 27×27 as unconfirmed. Exact 3-D dimensions.

### Implications for crosswordsmith
- **Odd square sizes dominate the blocked world** (15×15, 13×13, 23×23) — an odd-dimension default is correct and consistent with the engine's existing assumptions.
- **Realistic stock-grid / profile presets for the *current* block-only engine:** **15×15** (primary), **13×13** (quick), **23×23** (jumbo) — all blocked, all odd. These are the curatable, dozens-per-size targets from the stock-grid addendum.
- **Out of scope for the current engine:** the entire **~12×12 / 13×13 *barred* cluster** (Mephisto / Azed / Listener / etc.) — it needs the separate barred engine (closed-set tension addendum), not merely a size option. Likewise the **8×10 New Yorker** and **variable / thematic shapes** are "different product" territory.
- The bimodality reframes "support cryptic grids" as **two distinct targets**, not one: a blocked-15×15 path (in reach for the current engine) and a barred-~12×12 path (a separate build).

### Open questions
1. Azed's exact grid dimensions, and confirmed per-title sizes for Listener / Inquisitor / Enigmatic Variations / Spectator / Magpie (vs. the general 12×12 barred default).
2. Exact rows×cols (if fixed at all) for the US variety cryptics (Hex Atlantic/WSJ, Harper's, AVCX).
3. Whether any cryptic uses a verified size larger than 23×23 (the 27×27 figure was not found).
4. The dimensions/convention of 3-D cryptic grids (Sirius series).

### Sources (grid-sizes pass)
- *Cryptic crossword* — Wikipedia (15×15 standard; 12×12 barred variety; Quick Cryptic 13×13; Times Jumbo 23×23; Globe and Mail 24×24 barred; New Yorker 8×10 barred; American cryptics use British-style grids) — https://en.wikipedia.org/wiki/Cryptic_crossword
- Phi Crosswords (Paul Henderson), *Types of crossword* (15×15 daily profile; 12×12 barred; thematic shapes vary) — https://phionline.net.nz/types-of-crossword/
- bestforpuzzles — *The Times* (23×23 Jumbo) — https://bestforpuzzles.com/people/the-times.html ; *Advanced Cryptics* — https://bestforpuzzles.com/people/advanced-cryptic-crosswords.html ; *The Magpie* — https://bestforpuzzles.com/people/the-magpie.html
- Times for the Times — *Jumbo Cryptic* category — https://timesforthetimes.co.uk/category/jumbo-cryptic
- Fifteensquared — Everyman / Azed / Inquisitor categories — https://fifteensquared.net/
- Crossword Unclued — grid symmetry/checking; 3-D crosswords — https://www.crosswordunclued.com/2009/09/crossword-grid-symmetry.html ; https://www.crosswordunclued.com/2009/01/3-dimensional-crosswords.html
- Wikipedia — *Azed*; *Emily Cox and Henry Rathvon* — https://en.wikipedia.org/wiki/Azed ; https://en.wikipedia.org/wiki/Emily_Cox_and_Henry_Rathvon
- AVCX Cryptic Solving Guide (PDF; barred variety) — https://avxwords.com/wp-content/uploads/2020/06/AVxword_CrypticSolvingGuide.pdf
- The Nation crossword — https://www.thenation.com/content/crossword-puzzle/ ; Out of Left Field — https://www.leftfieldcryptics.com/

---

## Addendum — Two flavours, one CLI: the Flavour-A engine spec (design decisions)

Unlike the discovery addenda above, this section is a **design-decision record**, not external research. It operationalises the **Option A / Option B fork** from the closed-set-tension addendum into a concrete product shape: two "flavours" under one CLI, an explicit-subcommand surface, and an opinionated engine contract + objective for **Flavour A**. The decisions were reached by working through the design questions deliberately; where a choice closes off an alternative, that's noted.

### The two flavours

- **Flavour A — "arrange my words into a nice crossword-shaped layout."** Closed-set; takes a word list, returns an aesthetically-interlocked layout. Does **not** aspire to authentic cryptic legality (no symmetry guarantee, no ≥50%-checking guarantee). This is Option A from the closed-set addendum.
- **Flavour B — cryptic setting workflow tools.** Linter/validator, standard-format export, stock-grid library/profiles, and (the one heavy piece) a grid-first + open-dictionary auto-fill engine. Itself bifurcates into **blocked-cryptic** (reachable from the existing blocked cell model) and **barred-cryptic** (a separate engine — different cell/edge model). This is Option B.

**One CLI, two product *surfaces* — not two repos.** A shared core; divergent tops. Sequence: **A → B-blocked → B-barred**, B-barred furthest out.

### Shared substrate vs. divergent tops

| Layer | Shared | Flavour A | Flavour B |
|---|---|---|---|
| Input parsing / `meta` passthrough | ✅ | | |
| Blocked-grid cell/coordinate model | ✅ | | |
| Clue numbering | ✅ | | |
| Metric predicates (`word_meets_half/2`, `max_unch_run`, `checked_cells`, connectivity) | ✅ (as *measurement*) | used as optimizer signals | used as *validators* |
| Emit layer (JSON; later ipuz/Exolve) | ✅ | | |
| Fragment-grid seeding (see below) | ✅ (one primitive) | anchors | ninas / seed-then-fill |
| **Placement engine** | | the unified optimizer (new) | grid-first + dictionary fill (new, separate) |
| Legality core (`adj_is_free`, `no_word_merge`, `check_prev/next_cell`) | | free-canvas (blocked) | template-slot / barred (different) |

The discipline: **share the substrate, not the solver.** Forcing one engine to serve both via flags is the trap (the legality cores are incompatible — closed-set-emergent vs. grid-first-open-dictionary).

### CLI surface — explicit verbs, no catch-all default

Every capability is a verb; a bare invocation prints usage (never an action). Working shape:

```
crosswordsmith arrange [--strict | --best-effort] [--size N] [--size-mode fixed|max]
                       [--fragment grid.json] [--candidates N] --input words.json
crosswordsmith lint    --profile blocked-uk  layout.json          # Flavour B
crosswordsmith fill    --grid template.json --seeds seeds.json     # Flavour B (future)
crosswordsmith export  --to ipuz|exolve  layout.json              # Flavour B / shared
crosswordsmith                                                     # → usage/help
```

**`arrange` is the whole of Flavour A — a single verb, not the old `pack`+`solve` pair.** The engine-spec decisions below unify those two path-dependent engines into one; the old contracts become flag combinations:
- old `solve` (fixed-grid, place-all) ≡ `arrange --strict --size-mode fixed`
- old `pack` (engine-size, drop) ≡ `arrange --best-effort --size-mode max`

This also untangles the current wart where `GridLen` is a required positional but silently ignored in `--quality` mode. (`solve`'s only genuinely distinct capability is exhaustive enumeration/counting — **kept** as `arrange --enumerate`: emit/count *all* feasible full placements via `count_solutions`/`all_crossword`, not a separate engine.) **Migration:** moving from today's flag+positional CLI to subcommands is a breaking change — do it as one deliberate pass (rename, helpful "did you mean `arrange`?" error for old-style calls, update README / `run_tests.sh` / golden fixtures).

### Flavour-A engine: contract + objective

The decisions, and what each closes off:

| Dimension | Decision | Notes / what it rules out |
|---|---|---|
| **Scenario** | General-purpose (any word list; the website-TOC is one consumer) | Not TOC-only. |
| **Drop contract** | **Configurable, default strict** | Strict = fail only on a *genuinely unplaceable* word (no legal intersection anywhere). `--best-effort` drops + reports instead. |
| **Sizing** | Square `--size N` (default **15**) + `--size-mode` `fixed`/`max` (default `max`) | No free/auto-size, no rectangular grids — both deliberately cut. See "Sizing" below. |
| **Aesthetic goal** | **High interlock + even/balanced** | Explicitly *not* raw density/compactness (today's `quality.pl` objective) and *not* target-shape. |
| **Weak slots** | **Soft-penalty only** | Always include any legally-placeable word; stubs are penalised, never refused. Rules out a per-word *hard floor* that would drop/fail on weak (but legal) attachment. |
| **Interlock objective** | **Capped per-word reward** `Σ min(checked, ceil(L/2))` **+ small interlock tail `ε`** | Tunable `--check-target` / `ε`; balance comes from the cap; see detail + reachability caveat below. |
| **Outputs** | **Stable best + opt-in `--candidates N`** | Deterministic single by default (TOC / golden-diffable); `N` returns *meaningfully distinct* layouts (needs a diversity guard, not raw top-N). Rules out seeded-shuffle as the variety mechanism. |
| **Word distinction** | **Anchors yes (via fragment grid); priority scores NO** | Scores explicitly de-scoped. |

**Why this is a new engine.** Neither current engine does "place every word **and** optimise interlock": `quality.pl` is one-pass greedy (no search → can't guarantee place-all while optimising), `crossword.pl` is complete search but has *no objective* (pure feasibility). Flavour A is **`crossword.pl`'s backtracking search turned into branch-and-bound** — optimise an interlock objective over feasible full placements — reusing the shared legality core, emit, and metrics. This is the principled resolution of the "`pack`/`solve` feel path-dependent" observation: there is **one** Flavour-A engine; strict/drop and `--size`/`--size-mode` are *flags*, not separate engines.

**The objective, and why this exact form (performance + balance both point here):**
- The objective is a **capped per-word reward + small interlock tail**: `maximize Σ_w [ min(checked(w), target) + ε·checked(w) ]`, with **`target = ceil(L/2)` checked cells** by default (reuses `word_meets_half/2`; the reward is 0-marginal once a word meets-half). The **per-word cap is what creates *balance***: it stops rewarding already-checked words, redirecting the optimizer's effort to the laggards, and is anti-stub in a length-aware way. The tail `ε` (default ~0.2) keeps an interlock gradient *above* the cap and breaks ties toward crossier layouts. Both `target` and `ε` are tunable (`--check-target`). *(This is algebraically the "additive shortfall penalty" — `min(checked,target) = target − max(0,target−checked)` — just expressed as a reward, which makes the balance mechanism legible.)*
- This form is **additively decomposable** (each crossing changes only the two touched words' terms), so it updates incrementally per placement and yields tight, cheap admissible bounds for branch-and-bound pruning — the *same* cheap class as the raw total-crossings `placement_key` already uses.
- **Reachability caveat (load-bearing):** the cap only *balances* when `target` is attainable in the operating regime. Set so high that almost no word reaches it (closed-set fills can be sparse — the doc notes ~0.09–0.14 checked fraction), the cap never activates and the objective silently **degenerates to plain total-crossings (no balancing)**. So `target`/`ε` are **empirically tuned against the benchmark fixtures**, lowering `target` below `ceil(L/2)` if half proves unreachable at the densities Flavour A actually produces.
- **Avoid bare `maximin`/leximin as the search driver.** `min` is non-decomposable and its bound is useless early (dominated by unplaced words ≈ 0%), so it defeats pruning and degrades branch-and-bound toward exhaustive search; it also has plateau/"hostage" degeneracies. Reserve maximin/leximin only as a **last-resort tiebreak among the few top candidates** (where the search is already done) or as an **optional hard per-word floor** (which, as a *constraint* "min ≥ T" ≡ "every word ≥ T", is cheaply forward-checkable — the cheap form). The cost only bites when `min` is the objective being maximised.
- Net: the performant choice and the stated goals (soft-penalty, balanced, stub-averse) coincide — no quality-vs-tractability trade here.

### Sizing

Deliberately minimal — `free` (auto-size) and rectangular grids were both cut as unnecessary:
- **Square only, single `--size N`, default `15`** (the dominant blocked-cryptic size per the grid-sizes addendum). Keeps the existing `init_grid(GridLen)` square model — **no `(W,H)` plumbing.**
- **No `free`/auto-size mode.** The engine is *given* `N`, so there is no size-search outer loop (`quality.pl`'s `grid_candidates` retires) — simpler, faster, more deterministic. Auto-down-sizing is recovered by `max` mode (below).
- **`--size-mode {fixed|max}`, default `max`** — the *same* within-N placement search; the modes differ only in **output framing**:
  - **`max`** (default): `N` is a *ceiling*; emit the tight cropped bounding box ≤ N×N. Never sprawls; auto-shrinks for small inputs (the TOC / general case, and what `free` used to give — but bounded).
  - **`fixed`**: `N` is the *canvas*; emit a full N×N grid, uncovered cells as blocks ("give me a 15×15 puzzle").
- **Orthogonal to `--strict`/`--best-effort`:** either mode can fail (strict) or drop (best-effort) when the words don't fit within N×N.
- Lineage: `fixed` is the (square-only) generalisation of the old `solve`; `max` recovers what `free` gave, with a ceiling. Rectangular/non-square grids (e.g. the 8×10 New Yorker) are explicitly out of scope for this engine — revisit only if Flavour B ever needs them.

### Candidate selection (`--candidates N`)

Default output is a single deterministic best layout (golden-diffable — the TOC case). `--candidates N` opts into alternatives, but the optimizer's raw top-N are near-duplicates (the global optimum ± one nudged word), so they are selected for **diversity**:
- During search, retain a bounded buffer of the best ~K *complete* layouts (cheap — branch-and-bound encounters many).
- Select N greedily: take the best, then repeatedly add the next-best layout whose **distance from every already-picked layout ≥ τ**; stop at N.
- **Distance** = fraction of words placed differently (cell or orientation) — a Jaccard-style distance over the `word → (cell, dir)` placement set.
- **τ** tunable (default ~0.3: "≥30% of placements differ"), empirically calibrated alongside the objective knobs.
- Surfaces *diverse-and-good* layouts (a post-filter on the retained frontier), not "the N most diverse in existence" — the right tradeoff for "give me a few to choose from."

### The fragment-grid primitive (the anchor mechanism)

Anchors are expressed **not** as per-word annotations but as a **fragment grid**: the user supplies a *partial grid* (some words/cells pre-placed) and the engine solves from that starting state. Decided because it is cleaner *and* more powerful:

- **One primitive subsumes three features:** anchors (pre-placed words), **ninas** (pre-fixed letters along a path — Flavour B), and **interactive iteration** (place some, let the engine finish). Replaces "anchors + (later) nina-seeding" with a single *solve-from-a-partial-grid* mechanism.
- **It generalises machinery that already exists:** both engines already start from a seed — `quality.pl`'s `seed_word/7` (one word at a fixed cell/dir) and `crossword.pl`'s `StartLoc` first word. A fragment grid is just "user supplies the seed(s) and their positions, possibly many." So it's *less* invasive than threading anchor metadata through word selection: set `Grid = fragment`, `Placed = [fragment words]`, then run the same search.
- **Input/output symmetry / round-tripping:** crosswordsmith already *emits* a grid; ingesting one as a fragment gives **emit → user edits (pin a word, fix a cell) → feed back → re-solve**, converging input and output on one schema. Serves the git-diffable "puzzle source as code" ethos, and composes with `--candidates N` (tweak a candidate, re-seed, re-solve).
- **Separates *where* (grid) from *what* (word list).**

Schema & reconciliation decisions:
0. **Schema = the emit format, made partial.** The canonical fragment schema **is the JSON `arrange` emits** (so emit → edit → re-ingest round-trips), with the convention **presence = fixed**: any word/cell present is pinned, everything absent is for the engine to place. There is no separate "anchor" field — *being in the fragment is the anchor*. An optional **thin, hand-authorable convenience form** (e.g. `{answer, row, col, dir}` tuples + optional fixed cells) **desugars into** the canonical form — the same "core mechanism + sugar on top" pattern as the anchor shorthand.
1. **Words vs. letters:** support **both** (a placed word is a run of fixed cells whose identity is known; a lone fixed cell is a nina constraint). Words-only is a fine v1; the schema *allows* letter-level so Flavour-B ninas drop in later.
2. **Reconciliation by answer string:** fragment words are matched to `--input` by their **answer string** (the emit metadata join is already answer-keyed); the engine places the unmatched remainder. A fragment word absent from `--input` is an error, not a silent placement. Fixed letters not belonging to a list word are free nina constraints — **not entries**, so they don't count toward the strict "place all words" contract.
3. **Validation up front:** validate the fragment against the legality core and report conflicts (overlaps, adjacency/merge violations, unsatisfiable fixed letters) *before* searching — reusing existing validators.
4. **Size frame:** the fragment carries grid dimensions, which **set `N`** (so `--size` is redundant when a fragment is given; if both are passed and disagree, error). `--size-mode {fixed|max}` still governs framing of the *result*.

- **Tradeoff (accepted):** for the trivial "feature one link" case a whole fragment is heavier than a `--anchor HOME:center` hint — but a minimal fragment is one word at one cell, and lightweight anchor *sugar* can desugar into a fragment later. Make the fragment grid the **core mechanism**; treat per-word anchor shorthand as optional convenience.

*Deferred to implementation (not product decisions): the exact thin-form syntax, and how to disambiguate a fragment match when `--input` contains duplicate answer strings.*

### Re-tag of the recommended-features list, by flavour

| # | Feature (from "Recommended features") | Flavour |
|---|---|---|
| 1 | Cryptic grid linter / validator | **B** (uses shared metrics) |
| 2 | Native .ipuz emitter | **Shared** (export) |
| 3 | Enumeration per word | **Shared** |
| 4 | Default `max_unch_run`=2 + Times double-unch | **B** (cryptic lint) |
| 5 | ~~Scored-fill / drop-order by score~~ | **De-scoped** — priority scores ruled out for Flavour A |
| 6 | Cluing-potential annotation (`meta.cluing`) | **B** (computed at emit; metadata) |
| 7 | Exolve emitter | **Shared** / B-leaning (export) |
| 8 | Nina seeding | **Absorbed into the fragment-grid primitive** (Shared mechanism, B use) |
| 9 | `min_len:3` hard floor | **B** (cryptic) |
| 10 | House-style / grid profiles | **B** (lint) |
| 11 | Clue stockpile | **B** |
| 12 | Pattern/anagram query | **B** |
| 13 | Wordplay decomposition | **B** |
| 14 | Batch/deck interface + exit codes | **Shared** (CLI plumbing) |
| — | Unified optimising engine (branch-and-bound + interlock objective) | **A** (new core) |
| — | Fragment-grid seeding | **Shared** primitive (A anchors / B ninas) |
| — | `--candidates N` (diverse top-N; ≥ τ placement-distance) | **A** |
| — | `--strict`/`--best-effort`, `--size N` + `--size-mode {fixed\|max}`, `--enumerate` | **A** |
| — | Verb = **`arrange`** (single Flavour-A command) | **A** |

### Decisions resolved (follow-up design pass)

All open decisions from the first draft of this spec are now settled:
- **Verb:** Flavour A is **`arrange`** — a single command; old `pack`/`solve` become flag combinations.
- **`--enumerate`:** kept — exhaustive count/emit of *all* feasible full placements (the old `solve`'s one distinct capability).
- **Interlock objective:** capped per-word reward `Σ min(checked, ceil(L/2))` + small interlock tail `ε` (both tunable; balance comes from the cap — reachability caveat above).
- **Sizing:** square `--size N` (default 15) + `--size-mode {fixed|max}` (default `max`); no free mode, no rectangular grids.
- **Fragment-grid schema:** the emit format made partial (presence = fixed), round-trippable; optional thin convenience form desugars into it; reconciliation by answer string.
- **`--candidates N`:** diverse top-N by greedy selection (placement-distance ≥ τ, τ tunable).

*Deferred to implementation (not product decisions): the exact thin fragment-form syntax; duplicate-answer disambiguation; and empirical calibration of `ε`, `target`, and the diversity threshold `τ` against the benchmark fixtures.*
