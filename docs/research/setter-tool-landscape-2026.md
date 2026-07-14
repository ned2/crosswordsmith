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

> **Honesty note on coverage.** The verified claim corpus is heavy on
> **fill / word-list / format-conversion tooling** and lighter on cryptic-specific
> and rendering competitors. Tools that were in scope but returned no
> *independently verified* claims this pass — **Qxw, MyCrossword, Phil, Black Ink,
> Amuse Labs PuzzleMe, Sympathy/TEA, xword-dl, Crossword Nexus, Nutrimatic** — are
> carried forward from the prior doc at its confidence, not re-confirmed here, and
> are marked `[prior-doc, unverified 2026]` below. Do not treat their 2026 status
> as researched.

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

**Head-to-head benchmark design (proposed — not yet run).** The sharpest,
fairest target is **ingrid_core** (CLI, grid-first, open-source, MIT, uses the
*same* STW list):

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
   agnostic). *Vowel/common-letter sub-scores ship data-free first.*
4. **Nina seeding + feasibility report** (solver path). Leverage medium (cryptic
   audience) · effort medium · in-scope yes.
5. **Prolog-native pattern/anagram query endpoint.** Leverage low→medium · effort
   medium · in-scope yes. *Plays to unification; assistive not generative.*
6. **Barred-grid cell model.** Leverage medium (cryptic) · effort **high**
   (different cell model) · in-scope debatable. *Prior doc rightly flags this as
   near-a-separate-engine; defer.*

### `xword` tool

7. **Rendering quality/targets** (the differentiated half). Leverage medium ·
   effort low→medium · in-scope yes. *Compete on render, not convert.*
8. **`.xd` interchange** (only missing text format). Leverage low · effort low ·
   in-scope yes. *Optional; only if Puzzmo/`.xd` matters to you.*
   — **Do NOT** chase conversion breadth or byte-parity (table stakes; mature
   libraries own it).

### Interop / export

9. **Fill-quality report surfaced through `lint`.** Leverage medium · effort low ·
   in-scope yes. *Pairs with #1; reuses the per-word report machinery.*
   — `.puz`/`.jpz`/PDF stay a **kotwords handoff**, not native (prior-doc guidance
   unchanged; `.puz` can't represent bars anyway).

### Browser SDK (the open niche)

10. **Name & harden the embeddable WASM constructor SDK as first-class
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
1. **The measured** fill-quality/speed delta between scoreless MRV and
   ingrid_core's scored CSP on identical grids + STW — materially worse, or just
   occasionally worse corners? (Decides urgency of D#1.)
2. Is scored *`arrange`* even coherent given the closed-set model, or does scoring
   belong *only* in `fill`? (Section D takes the position that `arrange` gets
   tiebreak-only; worth validating in practice.)
3. Real adoption demand for an embeddable client-side constructor SDK — genuine
   empty niche, or empty because the market wants apps not SDKs? What bundle
   size / framework story makes it adoptable?
4. Do any of the unverified cryptic tools (Qxw et al.) already close ninas /
   barred grids / wordplay analysis better than assumed?

---

## Sources (primary unless noted)

**Construction tools / fill:**
- Crossword Compiler — Pro Grid Filler: https://www.crossword-compiler.com/ProFill.html · main: https://www.crossword-compiler.com/
- CrossFire (Beekeeper Labs): https://beekeeperlabs.com/crossfire/ · docs: https://beekeeperlabs.com/crossfire/docs/index.html
- Exet: https://github.com/viresh-ratnakar/exet
- Crosshare construct: https://crosshare.org/construct
- ingrid_core (Rust CLI/lib): https://github.com/rf-/ingrid_core · Ingrid app FAQ: https://ingrid.cx/faq/
- Qxw guide [prior-doc, unverified 2026]: https://www.quinapalus.com/qxw-guide-20200708.pdf
- Tool overview [blog]: https://communicrossings.com/constructing-crosswords-tools

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
