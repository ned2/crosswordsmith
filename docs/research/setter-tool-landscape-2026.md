# Crossword-setter tool landscape (2026) — where crosswordsmith competes

**Status:** research snapshot, 2026-07-15. **Method:** fan-out web search →
primary-source fetch → 3-vote adversarial verification (24/25 claims confirmed
3-0 against primary sources; 1 refuted) → synthesis. **Scope:** both of
crosswordsmith's surfaces — the SWI-Prolog CLI engine (`arrange`/`lint`/`fill`/
`export` + the WASM browser SDK) and the `xword` Python conversion/rendering
companion.

**Relationship to [`docs/cryptic-setter-research.md`](../cryptic-setter-research.md):**
that doc is thorough but **stale** — it predates the current module layout
(`quality.pl`/`crossword.pl`) and lists as "must-have/high-value" a pile of
things now *shipped*: the `lint` subcommand + profiles, native ipuz + Exolve
export, per-word enumeration, and the `fill` engine. It is also cryptic-only and
says nothing about the `xword` tool or the browser/WASM surface. **This doc
re-baselines the competitive picture against what crosswordsmith actually ships
today.** Where a prior-doc claim has changed, it is flagged inline. Treat the
cryptic doc as still-valid for setter *methodology* and cryptic conventions; use
*this* doc for the competitive/gap picture.

> **Honesty note on coverage.** The first pass's verified corpus was heavy on
> **fill / word-list / format-conversion tooling** and lighter on cryptic-specific
> and rendering competitors. A **[second and third verification pass](#second-and-third-verification-passes-2026-07-15)**
> (see the section at the foot of this doc) has since re-verified **all** the
> carried-over tools for 2026 against primary sources — Qxw, Nutrimatic,
> MyCrossword, Phil, Black Ink, Amuse PuzzleMe (second pass) and Qat,
> Sympathy/TEA/Wordplay Wizard, xword-dl, Crossword Nexus (third pass). No scoped
> tool remains unverified.

---

## Executive summary

crosswordsmith occupies a **genuinely distinctive position** — CLI-first,
JSON-native, deterministic, closed-set "place *these* words" `arrange`, plus a
client-side WASM constructor SDK — but it **trails the field on one load-bearing
axis: scored American fill.**

Three headline conclusions, each verified 3-0 against primary sources:

1. **Scored fill is table stakes for construction; crosswordsmith's scoreless
   MRV fill is its biggest competitive gap.** Every serious constructor —
   Crossword Compiler's Pro Grid Filler (400k+ scored words), CrossFire
   (`word;score` 1–100 + Grid/Final quality metrics), ingrid_core (embedded
   Spread the Wordlist, `--min-score` default 50), Exet and Crosshare — drives
   fill with per-word scores and quality heuristics. crosswordsmith places words
   *legally* but has no notion of fill *quality*.

2. **The fix is already license-clean and free.** Spread the Wordlist (CC
   BY-NC-SA 4.0, explicitly permits selling puzzles made with it) is the exact
   `word;score` data model the ecosystem standardised on — ingrid_core embeds it,
   Crosshare builds on it. crosswordsmith can close the gap with **no** paywalled
   (Chambers/Collins/Oxford) or personal-use-only (XWord Info) dependency.

3. **The scriptable client-side WASM *constructor* SDK is a near-empty niche.**
   The only comparable competitor is `wasm_crossword_generator` (Rust→WASM,
   generation-only). Exet and Crosshare construct in-browser but are **end-user
   apps, not embeddable SDKs**. This is defensible open ground worth naming as a
   first-class positioning.

Conversely, **format conversion is *not* a differentiator** — puzpy, kotwords,
xd-crossword-tools and Exolve's own scripts already cover the interchange graph
with published round-trip fidelity. `xword`'s edge is narrow: the deterministic
native-JSON hub plus multi-target rendering breadth (terminal/HTML/SVG/PNG).

---

## A. The 2026 landscape, re-scored against today's crosswordsmith

**Legend:** ✅ crosswordsmith leads / is distinctive · 🟰 matches · ❌ trails.
Confidence: 🟢 verified 3-0 this pass · 🟡 prior-doc, not re-verified 2026.

| Tool | Platform / cost | Fill model | Closest to crosswordsmith on… | vs crosswordsmith |
|---|---|---|---|---|
| **Crossword Compiler + Pro Grid Filler** 🟢 | Windows-only (v11, no Mac); subscription + paid PGF add-on | **Scored** (400k+ American list, per-length min-score thresholds, non-stop optimal search) | `fill` | ❌ trails on scored fill; ✅ leads on CLI/deterministic/open-format/cross-platform. *Prior-doc "dominant for most NYT puzzles" is **REFUTED** — 0-3.* |
| **CrossFire** 🟢 | Win/Mac/Linux (Java); **$50** one-time, lifetime updates, 1-hr no-save trial | **Scored** (`word;score` 1–100; Grid Score = neighbourhood quality, Final Score = whole-grid; bundled ~183k dict, scores 5–50); *why-unfillable* diagnostics; rebus | `fill` | ❌ trails on scored fill + fail diagnostics; ✅ leads on CLI/batch/deterministic/free |
| **ingrid_core** 🟢 | Rust **CLI + library, MIT**; `Ingrid` desktop app is proprietary freeware | **Scored** backtracking CSP (Balafoutis-derived); embeds Spread the Wordlist; `--min-score` default 50; ASCII `#`/`.` grids | **`fill` — closest architectural analog & sharpest benchmark target** | ❌ trails on scoring; 🟰 matches on grid-first CLI + license-clean; ✅ leads on JSON-native output + closed-set `arrange` + WASM |
| **Exet** 🟢 | Browser, **free MIT**, local-only | **Beam search** (light-specific regexp, preferred/excluded lists ≤50k), popularity threshold | **Overall surface — closest single competitor**; its built-in analysis overlaps `lint` | ❌ trails on cryptic wordplay-variety analysis; 🟰 matches on lint/symmetry/unch/connectivity; ✅ leads on CLI/scriptable/deterministic/closed-set |
| **Crosshare** 🟢 | Web, free, open-source | **Live** incremental autofill (STW + Peter Broda + expanded names); per-entry wordlist panel; Enter-to-regenerate | browser construction | ❌ trails on interactive UX + hosting; ✅ leads on CLI-first automation, diffable JSON, embeddable SDK |
| **Qxw** 🟡 | Linux/Win/Mac, free GPL | Backtracking + MRV + forward-checking; free lights; ~11 answer treatments; batch `qxw -b` | **Engine architecture** (MRV + backtracking) | 🟰 matches on MRV/batch; ❌ trails on bars + hex/circular topologies + answer treatments; ✅ leads on JSON + WASM. *Carried from prior doc; not re-verified 2026.* |
| **MyCrossword / Phil / Black Ink / Amuse PuzzleMe / Sympathy-TEA** 🟡 | (various) | — | — | Not re-verified this pass; see prior doc. Sympathy/TEA remained withdrawn as of the prior doc. |

**Net position.** crosswordsmith **leads** on being CLI-first, JSON-native,
deterministic/diffable, license-clean, and on the closed-set `arrange` model
(place *these specific* words, size-choosing, drop-if-unplaceable) that **no
grid-first competitor offers**. It **matches** the field on layout legality
(`lint`) and standard-format export. It **trails** on exactly one thing that
matters to constructors day-to-day: **scored fill quality.**

---

## B. The `xword` companion's competitors — conversion is table stakes

This is the category the prior research skipped entirely. The verdict: the
format-interchange graph is **already crowded and solved**, so byte-parity /
round-trip fidelity is **table stakes, not a differentiator** (🟢 verified).

| Tool | Covers | Note |
|---|---|---|
| **puzpy** (`alexdej/puzpy`) | `.puz` binary + `.txt` read/write | Claims **100% round-trip fidelity** on 9700+ real NYT/WaPo/Onion/WSJ puzzles (self-reported). Round-trip byte-fidelity is a *published baseline*. |
| **kotwords** (`jpd236`) | Reads PUZ/JPZ + many publisher/online formats; writes `.puz`/`.jpz`/`.ipuz`/**PDF** | Already the tool crosswordsmith's own docs cite for `.puz`/`.jpz`/PDF handoff. (ipuz-as-*read*-input is 2-1, not fully confirmed.) |
| **xd-crossword-tools** (`puzzmo-com`) | `puz ↔ xd ↔ json` | Covers the `.xd` text format + JSON. |
| **Exolve** (`viresh-ratnakar`) | Its **own** bidirectional `exolve ↔ puz` and `exolve ↔ ipuz` scripts + text importer | The exact exolve/ipuz interchange `xword` offers is already covered by Exolve upstream. |

**Where `xword` actually sits.** Its residual edge is narrow but real:

- **The deterministic native-JSON hub.** `xword` treats crosswordsmith's
  canonical JSON as a first-class interchange centre with `native ⇄ ipuz ⇄
  exolve` and format auto-detection — nobody else centres on that hub because
  nobody else emits it. This is a *consequence* of the engine, not an
  independent selling point.
- **Multi-target rendering breadth.** `view` (terminal) + `render`
  (terminal/HTML/SVG/**raster PNG**) in one tool is broader than the
  convert-only libraries. Rendering (esp. SVG/PNG/terminal) is the less-crowded
  half of `xword`'s surface.

**Format gaps that matter (in priority order):** `.puz` (the universally-accepted
binary; can't represent bars) — but reach it via kotwords, do **not** hand-roll a
CRC-16 writer (prior-doc guidance still holds). `.jpz`/PDF likewise via kotwords.
`.xd` (git-diffable text, Puzzmo ecosystem) is the one interchange target `xword`
doesn't touch and could add cheaply if `.xd` adoption matters to you.

**Recommendation for `xword`:** do not invest in conversion *breadth* to compete
— it's a losing race against mature libraries. Invest in the two things that are
actually differentiated: (1) rendering quality/targets, and (2) being the
faithful, deterministic bridge *out of* the native format. Keep engine↔xword
parity best-effort (as the spec already states); don't chase byte-parity with
`json_write_dict` vs `json.dumps`.

---

## C. The client-side WASM *constructor* SDK — a near-empty niche

crosswordsmith is now a scriptable client-side construction engine (SWI→WASM in a
Web Worker, typed JS SDK for `arrange`/`lint`/`export`). **This is crosswordsmith's
clearest open ground** (🟢 verified).

- The **only** comparable non-Prolog client-side WASM *construction* engine found
  is **`wasm_crossword_generator`** (Rust→WASM, generated TS types, a
  `CrosswordClient` class with `generate_crossword_puzzle` + a `SolutionConf
  {height,width,max_words,words,requirements}` config). It is **generation-oriented
  and less feature-complete** than a layout+lint+fill+export engine.
- **Exet** and **Crosshare** both construct in-browser but are **end-user apps,
  not embeddable scriptable SDKs** — you cannot `npm install` them into another
  web app and call `arrange()`.

So a typed JS SDK exposing `arrange`/`lint`/`export`(/`fill`) for **embedding into
other web apps** is a defensible niche with almost no direct competition.

**Open risk (worth naming, not resolved):** near-solitude in a niche can mean
"empty" *or* "no demand." The market visibly wants *apps* (Exet, Crosshare), and
whether there's real pull for an embeddable *constructor* SDK is unproven.
Adoptability levers to investigate before over-investing: **bundle size** (the
WASM payload plan in [`docs/plans/wasm-payload-performance.md`](../plans/wasm-payload-performance.md)
is directly on-point), a **clean typed API shape**, and **framework integration**
(React/Svelte examples).

---

## D. American scored-fill — the load-bearing gap, and how to close it

**The finding (🟢 verified 3-0).** crosswordsmith's `fill` is scoreless MRV: it
finds a *legal* fill but has no fill-*quality* signal, so it cannot compete on the
"clean fill" quality American constructors expect. Every competitor scores:

- **Crossword Compiler PGF:** 400k+ scored list, per-length minimum-score thresholds.
- **CrossFire:** `word;score` 1–100 (0 excludes); **Grid Score** (per-word
  neighbourhood quality, >0.9 = good) + **Final Score** (whole-grid).
- **ingrid_core:** filters candidates by `--min-score` (default 50) over the
  embedded scored list.

**The word-list is a solved, license-clean dependency.** **Spread the Wordlist**
(spreadthewordlist.com) is **CC BY-NC-SA 4.0**, explicitly permitting free-product
use with attribution *and selling puzzles made with it* — compatible with
crosswordsmith's license-clean identity. Snapshot (2026-07-01, will drift):
**314,276 entries; 120,178 scored 50** (the "clean" benchmark). Alternatives, all
permissive: the **Collaborative Word List** (MIT, crowdsourced, largest surveyed);
`christophsjones/crossword-wordlist` (~170k, 1–50 band). **Never** bundle
Chambers/Collins/Oxford or the personal-use-only XWord Info list.

**How scoring fits the closed-set model (design tension — flagged as open).** The
subtlety: scored fillers assume an **open dictionary they can reject words from**.
crosswordsmith's `arrange` is the opposite — a **closed set it must place all of**.
So scoring belongs cleanly in **`fill`** (open-dictionary, grid-first — the exact
ingrid_core shape) and only as a **tiebreak/drop-order signal** in `arrange`
(equal-density placements prefer better-checked high-score words; best-effort
drops lowest-value first). This matches the prior doc's "respect the closed-set
inversion" principle — scoring drives *drop-order and tiebreaks in `arrange`*, and
*candidate selection + `--min-score` in `fill`*, never magic word invention.

Concretely, in scope and license-clean:
1. **`fill`: ingest `word;score` lists; add `--min-score N`.** Rank/prune fill
   candidates by score. This is the direct ingrid_core-parity feature.
2. **`fill`/`lint`: report a fill-quality summary** (mean/min word score, count of
   sub-threshold entries) — matches CrossFire's Grid/Final score idea, and slots
   naturally into `lint`'s existing per-word report.
3. **`arrange`: score as tiebreak only** (see above), preserving determinism.

**Head-to-head benchmark — measured (prototype).** Run and quantified against
**ingrid_core** (CLI, grid-first, MIT, same STW list); full harness + results in
[`benchmarks/fill_quality/`](../../benchmarks/fill_quality/README.md). Every
placed entry scored *post hoc* against STW, so the gap is visible before
crosswordsmith implements any scoring. Headline numbers (score 50 = "clean"):

- **On grids both tools complete** (small, short-slot): crosswordsmith's scoreless
  MRV fill scores **mean 27–42, with entries at 0** — it grabs the
  alphabetically-first fitting strings (`AAAAA`, `AAAAQ`, `ABAAB`, `AEAEA`).
  ingrid at `--min-score 50` scores mean/min 50, zero below-clean. **A crude
  `score>=50` *dict prefilter* on crosswordsmith fully recovers ingrid-parity
  quality** — the cheapest form of scored fill already closes the quality gap.
- **On a hard grid** (`blocked_13a`, full 13-length slots): crosswordsmith
  **cannot complete** (budget-exhausted ~20s, any dict); ingrid completes in ~10s
  at `--min-score ≤30` but times out at 50. So there is a **second, orthogonal
  gap — search power on hard grids — separate from scoring**, and the harder one.

The original benchmark design (unchanged, for the fuller sweep still worth doing):

- **Inputs:** a fixed set of legal blocked grids (crosswordsmith's own stock
  grids + a few standard 15×15 American masks), same **STW dictionary** for both.
- **Procedure:** fill each grid with (a) crosswordsmith `fill` scoreless, (b)
  crosswordsmith `fill` *with* `--min-score` once implemented, (c) `ingrid_core`.
- **Metrics:** completion rate (fillable?), **mean/min word score of the fill**
  (the quality axis), count of sub-50 entries, and wall/inference time. Track
  score metrics under crosswordsmith's existing inference-based, machine-independent
  benchmark discipline.
- **The load-bearing question** (open until measured): does scoreless MRV produce
  *materially* worse fill, or just occasionally worse corners? That number decides
  how urgent scored fill really is.

---

## E. Cryptic-specific residual gaps (post-`lint`/`export`)

Re-checking the prior doc's cryptic wishlist against what now ships. **Now
covered** (were gaps, now done): the grid linter + house-style profiles, native
ipuz + Exolve export, per-word enumeration. **Still genuinely unserved:**

| Gap | Status | Notes |
|---|---|---|
| **Barred-grid representation** | Still absent | The cell model is block-only (`null`=empty). `barred-ximenean` lint applies the *checking math* to a blocked-shaped layout; a true bar/edge cell model is still out of scope. Qxw remains the free tool that does bars 🟡. |
| **Ninas** (constrained cell-path targets) | Still absent | Prior doc's design (pre-bind cells before MRV in the solver, emit a feasibility report) still stands as the clean approach. |
| **Cluing-potential scoring** (`meta.cluing`) | Still absent | *No existing tool emits this as a number* — still a genuine differentiator if built. Vowel-ratio/common-letter sub-scores are data-free; charade/container need a small bundled list. |
| **Clue stockpile keyed by answer** | Still absent | Pure `meta` data plumbing; stays on the layout side of the line. |
| **Nutrimatic-style pattern/anagram query** | Still absent 🟡 | Maps naturally onto Prolog unification/backtracking; Nutrimatic itself not re-verified 2026. |
| **Wordplay decomposition** | Still absent (correctly out of scope as *generation*) | Exet ships cryptic **wordplay-variety analysis** (🟢) — an *assistive* analog crosswordsmith lacks. Note: even Exet stops short of auto-writing clues, consistent with crosswordsmith's "feed clue tools, don't be one" identity. |

None of these are *newly urgent* — they were gaps before and remain so — but
cluing-potential scoring and the pattern/anagram query are the two that best fit
Prolog's strengths and the "assistive, not generative" line.

---

## F. Ranked gap list & roadmap, by surface

Tagged **leverage** (how much it closes a real competitive gap) × **effort** ×
**in-scope** (per crosswordsmith's design identity — layout/validate/fill/export
that *feeds* clue tools; CLI-first; deterministic; license-clean; metadata-agnostic).

### Prolog engine

1. **Scored `fill` (ingest `word;score`, `--min-score`, quality summary).**
   Leverage **very high** · effort **medium** · in-scope **yes**. *The one gap
   the whole field has and crosswordsmith doesn't. License-clean via STW. Start
   here.*
2. **Score as `arrange` tiebreak / drop-order.** Leverage medium · effort low ·
   in-scope yes (preserves determinism + closed-set model). *Natural follow-on to #1.*
3. **Cluing-potential annotation (`meta.cluing`).** Leverage medium (unique
   differentiator) · effort low→medium · in-scope yes (emit-time, solver stays
   agnostic). *Vowel/common-letter sub-scores ship data-free first.* MyCrossword's
   anagram-% / device-balance heuristic (second pass) is a ready blueprint.
4. **Clue stockpile keyed by answer** (`meta`: `{clueText, status, date}`, surfaced
   when a word is placed). Leverage medium (**unserved across both research
   passes** — a genuine differentiator) · effort **low** (pure `meta` plumbing) ·
   in-scope yes.
5. **Nina seeding + feasibility report** (solver path). Leverage medium (cryptic
   audience) · effort medium · in-scope yes. *Second pass: Qxw "free lights"
   confirm this is a real, still-served-nowhere-open capability.*
6. **Prolog-native pattern/anagram query endpoint.** Leverage low→medium · effort
   medium · in-scope yes. *Plays to unification; assistive not generative. Second
   pass: Nutrimatic's `A/C/V/#/_` + `<...>` + `&` syntax is a ready model.*
7. **Barred-grid cell model.** Leverage medium (cryptic) · effort **high**
   (different cell model) · in-scope debatable. *Only Qxw serves it (no JSON/ipuz
   export); build native, or a lossy Qxw-XML import as an on-ramp. Defer.*

### `xword` tool
(Full shortlist + backlog in [`docs/plans/xword-breadth-expansion.md`](../plans/xword-breadth-expansion.md).)

8. **Rendering quality/targets** (the differentiated half). Leverage medium ·
   effort low→medium · in-scope yes. *Compete on render, not convert.*
9. **`stats`/`inspect` + `diff` verbs** over the native-JSON hub. Leverage medium ·
   effort low · in-scope yes. *Read-only, deterministic, hub-amplifying.*
10. **`.puz` read (puzpy extra) + `.xd` interchange.** Leverage low→medium · effort
    low→medium · in-scope yes. *The two format gaps that matter.*
    — **Do NOT** chase conversion breadth or byte-parity (table stakes; mature
    libraries own it).

### Interop / export

11. **Fill-quality report surfaced through `lint`.** Leverage medium · effort low ·
    in-scope yes. *Pairs with #1; reuses the per-word report machinery.*
    — `.puz`/`.jpz`/PDF stay a **kotwords handoff**, not native (prior-doc guidance
    unchanged; `.puz` can't represent bars anyway).

### Browser SDK (the open niche)

12. **Name & harden the embeddable WASM constructor SDK as first-class
    positioning** (bundle size, typed API, framework examples). Leverage
    **high/strategic** (nearly uncontested niche) · effort medium (payload work
    already planned) · in-scope yes. *Bring `fill` to the SDK once #1 lands, so the
    browser surface is arrange+lint+export+fill.*

**One-line roadmap:** *scored `fill` first (closes the only real gap, license-clean
via STW), then surface fill-quality through `lint`/`arrange`, then invest the
browser SDK into the niche only it occupies — and stop competing on format
conversion.*

---

## Confidence, caveats, open questions

**Confidence.** 24 of 25 verified claims confirmed 3-0 against **primary** sources
(vendor pages, tool READMEs, official docs). One prior-doc-era claim **refuted**
0-3: *Crossword Compiler is dominant for most NYT puzzles* — not supported.

**Vendor/self-report caveats:** CC's "400,000+", CrossFire's feature copy, and
puzpy's "100% round-trip on 9700+" are self-reported and undated. STW's
120,178/314,276 is a 2026-07-01 snapshot and drifts. kotwords reading ipuz as
*input* is 2-1 (soft), not fully confirmed.

**Not re-verified this pass** (carried from prior doc at *its* confidence): Qxw,
MyCrossword, Phil, Black Ink, Amuse PuzzleMe, Sympathy/TEA, xword-dl, Crossword
Nexus, Nutrimatic. The competitive matrix is only as complete as the verified
corpus — heavy on fill/word-list/conversion, lighter on cryptic-specific and
rendering competitors. A follow-up pass on those tools would firm up sections A
and E.

**Open questions (load-bearing, unresolved):**
1. ~~The measured fill-quality/speed delta between scoreless MRV and ingrid_core~~
   **ANSWERED (prototype, [`benchmarks/fill_quality/`](../../benchmarks/fill_quality/README.md)):**
   *materially worse* — scoreless fill drops to non-word junk (mean 27–42, entries
   at 0) where scored fill holds mean/min 50; a `score>=50` dict prefilter recovers
   parity. Plus a *second* gap: crosswordsmith can't complete standard 13×13
   full-slot grids ingrid solves in ~10s (search power, not scoring). Remaining
   open: the delta across a proper spread of standard 11×11/15×15 masks.
2. Is scored *`arrange`* even coherent given the closed-set model, or does scoring
   belong *only* in `fill`? (Section D takes the position that `arrange` gets
   tiebreak-only; worth validating in practice.)
3. Real adoption demand for an embeddable client-side constructor SDK — genuine
   empty niche, or empty because the market wants apps not SDKs? What bundle
   size / framework story makes it adoptable?
4. Do any of the unverified cryptic tools (Qxw et al.) already close ninas /
   barred grids / wordplay analysis better than assumed?

---

## Second and third verification passes (2026-07-15)

Two focused follow-ups re-verified every tool the first pass carried over
unverified, asking one question per tool: **is it maintained in 2026, and does it
close any of crosswordsmith's residual gaps** — (1) barred-grid model, (2) ninas,
(3) cluing-potential scoring, (4) clue stockpile keyed by answer, (5)
pattern/anagram query, (6) wordplay decomposition, (7) scored fill? All claims
below verified 3-0 against primary sources unless noted. The **second pass**
covered Qxw/Nutrimatic/MyCrossword/Phil/Black Ink/PuzzleMe; the **third pass**
closed out the last four (Qat, Sympathy/TEA, xword-dl, Crossword Nexus).

| Tool | 2026 status | Closes which gap? |
|---|---|---|
| **Qxw** | GPL v2; last release **2020-07-08**, but site curated to Dec 2025 (semi-maintained). Rect/hex/circular grids, blocks+bars+mixed, auto-fill, answer treatments, batch `-b`. Exports EPS/SVG/HTML/PNG + an *experimental, lossy* Crossword-Compiler XML — **no JSON/ipuz**. | **Strongest evidence for (1) barred grids + (2) ninas/"free lights".** Does *not* close the JSON/ipuz conversion niche — that stays crosswordsmith's. |
| **Nutrimatic** | Open-source (`egnor/nutrimatic`, now PuzzleTechHub), self-hostable, maintained through Jun 2025 (sporadic). Syntax: `A/C/V/#/_/-` classes, `<...>` anagram, `&` intersection. | **Fully demonstrates (5) pattern/anagram query** — validates a Prolog-native equivalent (maps onto unification/backtracking). |
| **MyCrossword** | **Live & maintained** (open-source `t-blackwell/mycrossword` v2.3.0 Apr 2026; puzzles dated Jul 2026), free. Device tagging across ~28 clue types, indicator/abbrev references (Chambers filter), **anagram-% score** (advises ≤25%). | **Clearest exemplar of (3) cluing-potential scoring** — the *only documented, concrete scoring formula* found in the wild. Does **not** do (6) wordplay decomposition (that framing was refuted 1-2). |
| **Phil** | **Unmaintained since 2017** (2025 commit is just a redirect to a successor). SAT fill (Glucose 3.0), American blocked only. | None. |
| **Black Ink** | Current (Black Ink 2, $29.95, Mac App Store) but a **solver, not a constructor**; `.puz` only. | None. |
| **Amuse PuzzleMe** | Actively operated hosted construction/publishing platform; embeddable (JS/WordPress/iframe), ~€99/mo premium tier, plugin v1.4.0 May 2026. | None (it's a platform, not a gap-closer). |
| **Qat** (Mark Owen, quinapalus) | **Live/maintained** (page updated Dec 2025). *Server-hosted* web query engine (not self-hostable — the downloadable GPL sibling is Qxw the *constructor*); full operator set `. [] @ # * / & \| ! ~`, defaults to UKACD. Its CGI endpoint is one of Exet's autofill backends. | **Demonstrates (5) pattern/anagram query** — a third independent witness alongside Nutrimatic + umiaq. |
| **Sympathy / TEA / Wordplay Wizard** (Ross Beresford) | **Still withdrawn from sale** (official site confirms 2026); no successor/re-release. The associated **UKACD** word list stays **freely redistributable freeware** (notice-must-be-included; silent on modification, so *not* OSI-open). | None directly — but their withdrawal leaves the **UK-cryptic-construction niche empty** (context for crosswordsmith's cryptic direction). |
| **xword-dl** (thisisparker) | **Actively maintained** (v2025.10.14; commits into 2026). A **`.puz`-only downloader** scraping ~40 outlets. | None (a *downloader* — orthogonal). **Adjacency:** it emits `.puz`, so an `xword` `.puz` *read* would let `xword-dl → xword` pipeline (reinforces the xword breadth item). |
| **Crossword Nexus** (open-source, MIT; GitHub user, actively updated 2026) | Ships **umiaq** (MIT Python "Qat clone" pattern engine), **pypuz** (MIT Python — reads/writes `.puz`/`.ipuz`/`.jpz`/`.cfp`/AmuseLabs JSON, *not* `.xd`), **jscrossword** (JS, PUZ/JPZ/iPUZ/CFP + `puz2pdf`), **html5-crossword-solver** (browser player). | **Reusable prior art / partial competitor to `xword`.** umiaq → (5) pattern query; pypuz/jscrossword overlap `xword`'s convert surface (shared `.ipuz`) and are candidate reusable components for `.puz`/`.jpz` read. Closes none of (1)–(4),(6),(7). |

**Roadmap implications (net of both passes):**

- **Gaps (4) clue stockpile keyed by answer and (7)-as-*cryptic*-fill surfaced in
  NO re-verified tool.** Scored fill *is* ubiquitous among American *fillers*
  (pass 1: CC/CrossFire/ingrid) but absent from the cryptic tools; a **clue
  stockpile keyed by answer is unserved across both passes** — a genuine
  differentiator crosswordsmith could own (and it's pure `meta` plumbing, squarely
  in-scope).
- **Qxw validates barred grids + ninas but exports no JSON/ipuz.** Open strategic
  fork: implement the barred-grid model *natively*, or add Qxw-XML *import* to
  `xword`'s converter (cheaper, but Qxw's XML is lossy — drops topologies/answer
  treatments/free-lights, so it can't faithfully carry the very features that make
  Qxw interesting). Leaning native for a real barred model; XML import only as a
  low-fidelity on-ramp.
- **MyCrossword's anagram-% / device-balance heuristic is a ready blueprint for
  gap (3)** — the concrete first cut of cluing-potential scoring, if pursued.
- **Gap (5) pattern/anagram query is now triple-validated** (Nutrimatic, Qat,
  umiaq) with *converging* operator syntax — a well-understood capability, and
  **umiaq (MIT Python) is a reusable reference** for a Prolog-native equivalent.
- **The `xword` `.puz`-read item gains two reusable, license-clean options**
  (both MIT Python): **puzpy** (the plan's pick) and **pypuz** (also does
  `.jpz`/`.cfp`). And **xword-dl** (`.puz`-only downloader) makes `.puz` read the
  natural on-ramp for a `xword-dl → xword convert/render` pipeline.
- **UKACD stays freely redistributable freeware** — reconfirms crosswordsmith's
  license-clean lexicon story, but note it's silent on modification (not OSI-open),
  so ship its notice verbatim (as the repo already does for UKACD18).

**Caveats:** `quinapalus.com` returns HTTP 403 to default fetchers — Qxw/Qat facts
were confirmed via curl with a browser user-agent (+ manpages / the guide PDF for
Qxw); a few exact quoted strings not byte-verified. "Actively maintained" is
generous for Qxw (2020 release), Nutrimatic (Jun 2025), and Qat (Dec 2025 page
edit). Crossword Nexus is a GitHub *user* not an org; jscrossword's SPDX is
NOASSERTION though its README states MIT (treat as README-stated). All four scoped
tools are now verified; nothing scoped remains outstanding.

## Sources (primary unless noted)

**Construction tools / fill:**
- Crossword Compiler — Pro Grid Filler: https://www.crossword-compiler.com/ProFill.html · main: https://www.crossword-compiler.com/
- CrossFire (Beekeeper Labs): https://beekeeperlabs.com/crossfire/ · docs: https://beekeeperlabs.com/crossfire/docs/index.html
- Exet: https://github.com/viresh-ratnakar/exet
- Crosshare construct: https://crosshare.org/construct
- ingrid_core (Rust CLI/lib): https://github.com/rf-/ingrid_core · Ingrid app FAQ: https://ingrid.cx/faq/
- Tool overview [blog]: https://communicrossings.com/constructing-crosswords-tools

**Second-pass verification (2026-07-15):**
- Qxw: https://www.quinapalus.com/qxw.html · guide: https://www.quinapalus.com/qxw-guide-20200708.pdf
- Nutrimatic: https://nutrimatic.org/2024/usage.html · source: https://github.com/egnor/nutrimatic
- MyCrossword: https://www.mycrossword.co.uk/crosswords/user-guide · https://mycrossword.co.uk/about
- Phil: https://github.com/keiranking/Phil
- Black Ink: https://redsweater.com/blackink/
- Amuse PuzzleMe: https://amuselabs.com/features/

**Third-pass verification (2026-07-15):**
- Qat: https://www.quinapalus.com/qat.html
- Sympathy/TEA/Wordplay Wizard (withdrawn): https://www.rossberesford.com/crosswordman · UKACD licence: https://www.wirdz.com/ukacdasc.html
- xword-dl: https://github.com/thisisparker/xword-dl
- Crossword Nexus — umiaq: https://github.com/crosswordnexus/umiaq · pypuz: https://github.com/crosswordnexus/pypuz · jscrossword: https://github.com/crosswordnexus/jscrossword

**Word lists (license/scoring):**
- Spread the Wordlist (CC BY-NC-SA 4.0): https://www.spreadthewordlist.com/
- crossword-wordlist (christophsjones): https://github.com/christophsjones/crossword-wordlist
- Wordlist shopping [blog]: https://tcampbell.substack.com/p/shopping-for-wordlists

**Formats / conversion / rendering (`xword` competitors):**
- puzpy: https://github.com/alexdej/puzpy
- kotwords: https://github.com/jpd236/kotwords
- xd-crossword-tools: https://github.com/puzzmo-com/xd-crossword-tools
- Exolve (own puz/ipuz scripts): https://github.com/viresh-ratnakar/exolve

**Browser / WASM constructor SDK:**
- wasm_crossword_generator: https://github.com/krhoda/wasm_crossword_generator

**Cryptic / query [mixed confidence]:**
- Nutrimatic usage [unverified 2026]: https://nutrimatic.org/2024/usage.html
- Puzzle-tool index [secondary]: https://signals.mysteryleague.com/tools/
