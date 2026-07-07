# Spec: `xword` — a terminal viewer & format multitool for crossword layouts

Status: **building** — Phases 0–3 done, see [`xword-status.md`](xword-status.md).
This file is the single source of truth; it absorbed the original standalone
decision record (2026-07-03).

`xword` is a standalone **Python** tool, a downstream consumer of the
crosswordsmith engine — **not** part of the engine's four-verb design-spec. The
Prolog core stays pure and dependency-light; `xword` is a separate toolchain.

## Background & rationale — why this tool, this shape

Origin (2026-07-03): we wanted to *see* the boards `crosswordsmith arrange` (and
the other verbs) produce, and weighed keeping rendering in Prolog against
building a separate tool on the pipe. The analysis that fixed this design:

- **Rendering is not the hard part.** All three formats carry the grid as a
  near-render-ready structure — native `grid[][]` (`null` = block, else
  `{letter, number, across, down}`) *is* the board; Exolve is a letter/`.` text
  grid; ipuz is `puzzle[][]` + `solution[][]`. A box-drawn grid is ~50–80 lines
  against any of them. The real decision axes are **static-vs-interactive**,
  **distribution**, and **where this sits relative to the pure-Prolog core**.
- **"Render the board" hides two tools** — a static pipe filter (`… | xword`)
  and an interactive TUI, which *takes over* the terminal and cannot sit in a
  pipe. They pull opposite ways, so they are **sequenced**, not conflated: the
  Rich-based static filter first (Phases 1–3), the Textual TUI later (Phase 4),
  both over the same `Board`.
- **Leaving Prolog is right — not because Prolog can't render** (it can; it's
  `format/2`) but for (a) **architectural cleanliness** — the engine "feeds
  tools; it is not one," so a viewer is a downstream consumer and presentation
  stays out of the deterministic core; (b) a **format-agnostic** viewer
  shouldn't be re-implemented in Prolog; (c) the Rich/Textual ecosystem, which
  Prolog lacks.
- **Python over Go/Rust** matches the maintainer's preference and maps onto the
  two tiers (Rich → Textual — one rendering lineage, one model reused).
  Distribution — Python's weak spot vs a single binary — is handled by shipping
  a package run via `uv`/`uvx`/`pipx` (§12); a Go rewrite is revisited only if
  we ever distribute to non-Python users, not before.

## 1. Problem statement

The engine emits deterministic layout data (native canonical JSON) and can
transform it to interchange formats (ipuz, Exolve via its `export` verb), but
there is **no way to *see* a layout** short of reading JSON, and no way to move
*between* the interchange formats or into visual artefacts. A setter wants to
preview "how this looks in print," convert a layout to whatever a downstream
tool eats, and produce a shareable image — none of which belongs inside a
deterministic Prolog engine whose job is to *feed* tools, not *be* one.

`xword` is that tool: pipe a layout in, and view it, convert it, or render it.

## 2. Goals

- **G1.** **View** any supported layout in the terminal, emulating the printed
  puzzle — grid + Across/Down clue lists — with a solved and a blank mode.
- **G2.** **Convert** between supported data formats in any direction
  (*any → any*), preserving semantics and staying git-diffable.
- **G3.** **Render** a layout to a visual artefact (SVG master; HTML; PNG/PDF
  derived from the SVG behind an optional dependency).
- **G4.** Be a **good Unix citizen**: read stdin, write stdout, compose in
  pipes, honour TTY detection and `NO_COLOR`, `--out` to a file.
- **G5.** **Grow cheaply**: one normalised model so a new input format is one
  parser, a new data format one serializer, a new visual target one renderer.
- **G6.** Open a clean path to an **interactive Textual TUI** as a later
  renderer over the *same* model — no second parse/layout path.

## 3. Non-goals

- Not a clue-writing tool, and not a solver or setter (consistent with the
  engine's boundary).
- Does **not** modify the Prolog engine or its output formats.
- Adds **no** dependency to the Prolog core.
- v1 does not read/write binary/print formats directly (`.puz`, `.jpz`, PDF as
  *input*); those stay a later parser or an offload to kotwords (§9).

## 4. Decisions taken

- **D1. Three verbs: `view`, `convert`, `render`.** Deliberately avoids the
  engine's `export` (which means *data* export). `view` presents to the
  terminal; `render` presents to a visual file; `convert` is data interchange.
  `view` and `render` are the **same operation, different backend** (§5.2).
- **D2. One `Board` model is the spine.** Every command is
  `parse → Board → {serialize | present}`. Data plane = `convert`; presentation
  plane = `view`/`render`/(future Textual). See §5.
- **D3. `view` defaults to the *solved* grid**, with `--blank` for the empty
  reader's view. (A setter most often wants to see the fill; the blank preview
  is a toggle, not the default.)
- **D4. SVG is the master visual format.** HTML is a separate semantic page;
  **PNG and PDF are conversions *of the SVG***, gated behind an optional
  `xword[raster]` extra — one geometry/layout implementation feeds terminal,
  HTML, SVG, and thence PNG/PDF.
- **D5. v1 input/output formats: native canonical JSON, ipuz v2, Exolve.**
  `.puz` and friends are deferred (§9).
- **D6. Deterministic output everywhere.** No timestamps, no random ids, stable
  key/element ordering — so `convert`/`render` output is git-diffable and
  golden-testable, mirroring the engine's contract.
- **D7. Strict on *structure*, drop-and-warn on *metadata*.** Source properties
  split in two:
  - **Puzzle structure** (grid shape, blocks, bars, rebus, cell styling,
    numbering, letters, clues, enumeration) — **strict**: if the target cannot
    represent it, the conversion **fails**. Dropping it would yield a
    *different puzzle*.
  - **Metadata annotations** (native's per-word `link`, arbitrary `meta` keys) —
    **dropped** when the target has no home for them, with a **warning on
    stderr** (never stdout). Rationale: they don't change the puzzle, and
    blocking on them would make `convert` frustrating for a user who already
    knows the target can't hold them and accepts the loss.

  So there are no `--force`/`--best-effort` fidelity flags in v1, only a
  `-q/--quiet` to silence the drop warnings. The two former "deferred escape
  hatches" are now decided (D9): **uplifting native's model** is adopted on
  lossless-only terms (title/author landed P1+P2; styling + enumeration
  separators follow), while **best-effort *structural* conversion**
  (crop/flatten with warnings) is **rejected** — it would produce a different
  or invalid puzzle.
- **D8. CLI built with [`cyclopts`](https://github.com/BrianPugh/cyclopts).**
  Python-native subcommand dispatch for `view`/`convert`/`render` with typed
  params, enum `--from`/`--to`, and per-verb help. Its idioms and the resulting
  app skeleton (package layout, `xword` console entry point, stdin-as-default,
  exit codes) are pinned down by spike **S1** before Phase 1 (§13,
  [`xword-spikes.md`](xword-spikes.md)).
- **D9. Native-schema uplift — additive-optional anchors + riders,
  lossless-only.** (Ratified 2026-07-07,
  [`plans/native-schema-uplift-plan.md`](plans/native-schema-uplift-plan.md).)
  The native schema is extended where the engine can *compute-or-carry*, every
  conversion stays lossless, and the lossy best-effort escape hatch (§10 (a)) is
  **rejected**, not deferred. **Landed:** optional top-level `title`/`author`
  (P1 engine + P2 xword). **Coming under the same rule:** closed-subset cell
  styling (P3) and the `'`/`.` enumeration separators (P4). Every addition is
  **additive-optional** — an absent field means today's semantics, old layouts
  stay valid, and nothing new is ever required; since neither the native payload
  nor json-output-spec carries a `version` field, additive-optional *is* the
  compatibility story. **Permanent fail-strict** (no uplift, no lossy fallback):
  rectangular grids and rebus cells — native's model is square-N×N single-letter,
  and flattening either yields a *different* or *invalid* puzzle (ipuz is the
  rectangular/rebus carrier). **Bars** stay **gated** on
  [`cryptic-layout-spec.md`](cryptic-layout-spec.md) landing.

## 5. Architecture

```
                          ┌── serialize → native / ipuz / exolve   →  convert   (data plane)
  any input ── parse ───► │ Board │
  native/ipuz/exolve      └── present ──► terminal (Rich)          →  view      (presentation
                                     ├──► svg / html / png / pdf   →  render     plane)
                                     └──► interactive (Textual)    →  (future)
```

### 5.1 The `Board` model

One normalised structure both planes consume. It is the *union superset* of
what the supported formats carry, tolerant of missing pieces (a stock-grid
template has blocks and numbers but no letters or clues).

| Field | Meaning |
|---|---|
| `height`, `width` | grid dimensions. Native is square (`gridLength`); ipuz/Exolve may be rectangular, so the model carries both dims. |
| `grid[r][c]` | a **block** or a **cell** `{letter?, number?, across_num?, down_num?}`. Row-major, `[row,col]`, 0-indexed — matching the native `grid`. |
| `words[]` | `{number, direction, cells[], answer?, enumeration, clue?, meta?}`. |
| `meta` | puzzle-level: `title?`, `author?` — **all three formats carry these** (native as optional top-level anchors since the D9 uplift; ipuz/Exolve natively). |

**Numbering is derivable, not required from input.** Standard crossword
numbering (scan row-major; a cell is numbered iff an across or down word starts
there) is computed by `xword` when a source omits it — so a bare block template
still renders correct labels.

**Enumeration** (`(5,5)`, `(4-5)`) is normalised from the source: derived from
the answer's spaces/hyphens for native (as the engine's `export` does), taken
from the field for ipuz, parsed from the clue for Exolve. "Invent nothing": a
word with no clue carries an empty clue, never a fabricated one.

The **Exolve parse-back rule** (S5) is normative, because Exolve embeds the
enumeration in the clue string (`… (5,5)`): the enumeration is the **last**
parenthesised group anchored at end-of-string whose interior is only the enum
alphabet — digits and the separators `, - . ' ` and space
(`\s*\(([0-9]+(?:[,\-. ']*[0-9]+)*)\)\s*$`); everything before it (rstripped) is
the clue. This keeps clue-internal parens intact (`Send (a letter) to HQ (4)` →
clue `Send (a letter) to HQ`, enum `(4)`). It is irreducibly ambiguous when the
clue text itself ends in a bare `(n)` (`Beethoven's Fifth (5)`, `See 12 (4)`);
the **disambiguator is the grid** — accept the enum only if its integer runs sum
to the word's cell count, else treat the parens as clue text. This reduces but
does not fully eliminate the ambiguity (`(10)` on a real 10-cell entry still
collides); on mismatch, derive the enum from the cell count.

**Cell-encoding gotchas for the parsers** (S4/S5): ipuz `block`/`empty` are
*declared* top-level fields with defaults `"#"`/`"0"` and must be **read, not
hardcoded** — the engine emits an **integer** `0` for `empty` where real-world
files use the string `"0"`; and a styled/labelled ipuz cell is an **object**
(`{cell, style}`), so a block test must special-case dicts. In Exolve, an
**unfilled light cell is `0`** (matters for `--blank`/stock-grid templates) and
**decorators are cell suffixes** (bars `|`/`_`, circle `@`, prefill `!`, rebus
via `exolve-option: rebus-cells`), so grid rows must be **tokenised**
(main-char + trailing decorator run), never column-indexed.

### 5.2 Planes

- **Data plane** — `parse(fmt) → Board` and `Board → serialize(fmt)`. `convert`
  is their composition.
- **Presentation plane** — `Board → present(renderer)`, where a renderer is
  terminal (Rich), svg, html, or (later) Textual. `view` = present→terminal;
  `render` = present→file. The **grid geometry and clue-list layout are
  computed once** and shared by all renderers; a renderer only decides how a
  cell/number/block/clue is *drawn* in its medium. The terminal renderer is a
  **pure `BoardGeom → Rich renderable`** function built from `Text`/`Segment`
  (not Rich `Table` — see S2); the future Textual widget's `render()` returns
  the same renderable over the identical `BoardGeom`, so the seam (S7) is only
  interactivity (cursor/highlight = restyling cells in the same geometry) —
  validated low-risk, no separate spike needed.

## 6. Commands

Common CLI conventions (§8) apply to all three.

### 6.1 `view` — see it in the terminal

    crosswordsmith arrange --size 17 --input clues.json | xword view
    xword view --from ipuz --blank puzzle.ipuz

- Reads a layout (stdin or a positional file), detects the format (§7) or takes
  `--from`, and renders **grid + Across/Down clue lists** laid out like the
  printed puzzle.
- **`--blank`** renders the empty reader's view (no letters, numbers kept);
  default shows the **solved** fill (D3).
- Colour is optional, on only for a TTY and when `NO_COLOR` is unset —
  `enabled = stdout.isatty() and not os.environ.get("NO_COLOR")`, passed to
  `Console(no_color=not enabled, force_terminal=enabled)` (S2).
- **Cell geometry (resolved S2):** 3-char-wide × 1-row-tall cell interior with
  single-line box borders (`─ │ ┼ …`); the 4×2 char pitch cancels the terminal's
  ~2:1 aspect so the grid reads square (a 17-grid is 69w×36h, fits 80 cols).
  Blocks are solid `█`×3 (**not** reverse-video — an SGR attribute vanishes in
  `NO_COLOR`/piped/golden output, making a block indistinguishable from empty).
  Corner numbers are **superscript digits inline before the letter** (`¹O `);
  colour = `bold cyan` numbers / default letters / `grey42` borders.
- Terminal-width aware: clue lists wrap; the grid **never wraps**. It is rendered
  at fixed geometry and printed with Rich reflow disabled (`console.print(grid,
  soft_wrap=True)`); if the grid is wider than `console.width` it overflows /
  scrolls (optionally a stderr note), never garbled by wrapping.

### 6.2 `convert` — data-format interchange

    xword convert --to exolve < layout.json > puzzle.exolve
    xword convert --from ipuz --to native puzzle.ipuz --out layout.json

- `--to native|ipuz|exolve` **required**; `--from` optional (else detect).
- *Any → any* among the three (the reverse/cross directions the engine's
  one-way `export` cannot do are the point).
- **Strict on structure, drop-and-warn on metadata (D7, §10)** — a conversion
  *fails* if the target cannot hold a **structural** property; native's
  **metadata annotations** (`link`, arbitrary `meta`) are *dropped* with a
  **stderr** warning when the target has no home. `-q/--quiet` silences the
  warning; stdout is never touched by it.
- **Invent nothing** — enumeration/clue are *derived* where absent, never
  fabricated; nothing puzzle-level (e.g. a placeholder title) is manufactured.
  **One sanctioned exception (Q5, §14):** a title-less board targeting
  **Exolve** gains the engine's default `exolve-title: Untitled` at the
  `convert` boundary (warned on stderr, `-q`-silenceable) — an ecosystem
  requirement (Exet's Save crashes on a null title), not parity-chasing. The
  serializers themselves stay invent-nothing; ipuz output invents nothing.
- **`native → native` is payload-lossless.** The native target holds
  everything metadata-class (per-word `link`/`meta` verbatim, top-level
  `diagnostics` passed through), so D7's drop-when-homeless rule never
  triggers and no warning is emitted. `diagnostics` is carried on the `Board`
  and re-emitted only when present (a payload without it stays without it).
- **No `exolve-id` is emitted** (invent-nothing; resolved Phase 2). Exolve
  documents the directive as optional: when absent it "create[s] one from a
  signature of the grid and the clues", so consumers require nothing — the
  only cost is that solving state keyed on the auto-id resets if the puzzle
  is edited, which an invented constant id would not fix anyway. Since the
  uplift's P1 (D9) the **engine agrees** — it dropped its constant
  `exolve-id: crosswordsmith-export`, so neither side emits one.
- **Engine cross-check** — dropping `link` now matches the engine's own
  `export`. Since the native-schema uplift's P1 (D9) the engine invents a
  default title on **Exolve only** (`exolve-title: Untitled`, confirmed S5) —
  the invented ipuz `title:"Untitled"` was retired — and `xword` matches that
  Exolve default at the convert boundary, so the title line is byte-identical
  on both sides while ipuz carries no title on either. With `exolve-id` gone
  from both (above), the cross-check is **structural** agreement (§11) — the
  guaranteed contract. **Byte-level parity is best-effort, not claimed either
  way** (§14): `xword`'s Exolve header now leads with title/setter like the
  engine's, so in practice the two currently match byte-for-byte; ipuz differs
  only in JSON serializer formatting (SWI `json_write_dict` vs Python
  `json.dumps`). Neither is pinned as a byte-identity contract.

### 6.3 `render` — visual artefact

    xword render --to svg  < layout.json > puzzle.svg
    xword render --to html --blank layout.json --out puzzle.html
    xword render --to png  layout.json --out puzzle.png     # needs xword[raster]

- `--to svg|html|png|pdf` **required**. `--blank` applies here too (D3).
- **SVG** = self-contained master: title + grid + two-column clue lists,
  deterministic element ids (`cell-r{r}-c{c}`, `word-{n}-{dir}`); glyphs are
  `<text>` (Q3 resolved — cairosvg embeds outlines in the PDF regardless).
- **HTML** = semantic page: grid as inline SVG + clue lists as HTML (`<ol>` with
  per-clue `value`/id), styled by **one built-in inline `<style>`** — the page
  is self-contained and deterministic (Q2 resolved); a `--css` hook is deferred.
- **PNG/PDF** = rasterised/converted *from the SVG* (`xword[raster]` extra,
  `cairosvg`); no second layout path. Tested by dimensions/magic, not
  byte-golden (host-font `<text>`, §11).
- **Binary-to-TTY guard**: `png`/`pdf` refuse to write to a terminal — require
  `--out` or a redirect.

## 7. Format detection

An explicit ladder; `--from` always overrides.

1. **Binary magic** — reserved for future `.puz` (`ACROSS&DOWN`); absent in v1.
2. **JSON parses?** → check **ipuz first** (S4): the signal is a `version` or
   `kind` **value that contains the `ipuz.org` host** (`kind` is a list of such
   URLs) — not the mere presence of a `version` key. Else **native**
   (`grid` + `words` + `gridLength`).
3. **Otherwise text** → **Exolve** markers (`exolve-begin` / `exolve-width`).
4. No match → a clear error naming the tried formats and suggesting `--from`.

## 8. CLI conventions (all verbs)

- **stdin → stdout** by default; a positional file overrides stdin; `--out
  FILE` overrides stdout (written only on success — a failed run leaves no
  partial file, matching the engine).
- **`--from FMT`** overrides detection; **`--to FMT`** selects the target where
  applicable.
- **TTY / `NO_COLOR`** honoured for `view`; colour off when piped.
- **Diagnostics on stderr, data on stdout.** Metadata-drop warnings (§10) and
  all errors go to stderr, so a piped/redirected stdout is always clean data.
  **`-q/--quiet`** silences the drop warnings.
- **`--help`/`-h`** per verb; a bare `xword` prints usage and exits non-zero
  (mirrors the engine's CLI shape). **Impl note (S1):** this is *not* cyclopts'
  default — a bare app prints help and exits 0. Register an `@app.default`
  handler that prints help to a stderr `Console` and raises `SystemExit`.
  cyclopts already renders parse/validation errors to **stderr** with exit 1
  automatically, so only *domain* errors need a manual hook in a thin `main()`.
- Errors are clear and actionable (bad format, non-square grid on `→native`,
  missing raster extra), non-zero exit. The **missing-raster** case has two
  distinct forms the error must tell apart (S6): the `[raster]` extra not
  installed (`ImportError: cairosvg`) vs the extra installed but the **system
  cairo library missing** (`cairocffi` `OSError` on `libcairo.so.2`) — the
  latter points the user at `apt install libcairo2` / distro equivalent.

## 9. Input-format scope & the `.puz` boundary

v1 parses/serializes **native, ipuz v2, Exolve** (§5.1). `.puz`/`.jpz`/PDF are
**deferred**: the engine already reaches them *outbound* via kotwords, so
`xword` gains most by first covering the three text/JSON formats it can handle
purely. When added, `.puz` is one more parser/serializer (via `puzpy` or a
kotwords offload) behind the same `Board`.

## 10. Conversion fidelity

Loss is not random per-format quirks — it follows from the formats having
different **expressive power**, with native the narrowest box:

| Property | Class | native | ipuz v2 | Exolve |
|---|---|---|---|---|
| Rectangular grid | structure | ✗ (square N×N) | ✓ | ✓ |
| Multi-letter (rebus) cells | structure | ✗ | ✓ | ✓ |
| Circled / shaded cells, **bars** | structure | ✗ | ✓ | ✓ |
| Puzzle title / author | structure | ✓ (D9) | ✓ | ✓ |
| Blocks, numbering, clues, enumeration, single letters | structure | ✓ | ✓ | ✓ |
| Per-word `link` + arbitrary `meta` | **metadata** | ✓ (only native) | ✗ | ✗ |
| Top-level `diagnostics` (engine quality caveats, json-output-spec §6.4) | **metadata** | ✓ (only native) | ✗ | ✗ |

**The rule (D7)** splits on the *Class* column:

- **structure** the target shows as `✗` → the conversion **fails**. No silent
  crop or flatten; the error names the property and suggests a target that
  supports it. (Dropping structure would produce a *different puzzle*.)
- **metadata** the target shows as `✗` → **dropped**, with a **stderr** warning
  naming what was dropped (`-q/--quiet` silences it). It doesn't change the
  puzzle, and a user targeting ipuz/Exolve is assumed to accept that those
  formats can't carry native's annotations.

Consequences of the table:

- **The one metadata case is native → {ipuz, Exolve}** — native alone carries
  `link`/`meta` and `diagnostics`. So converting **current crosswordsmith
  output (every word has a `link`) to ipuz/Exolve now succeeds**, dropping the
  links (and any `diagnostics`) with a warning. **native → native drops
  nothing** — the target holds it all (§6.2), so the identity conversion is
  payload-lossless and silent.
- **Structural loss is the blocking direction, mostly `→ native`** (native can't
  hold rectangular grids, rebus, styling, bars) — the native-schema uplift (D9)
  splits that list: rectangular and rebus are **permanent** fail-strict, bars
  wait on cryptic-layout-spec, cell styling gets a lossless closed-subset in P3,
  and puzzle **title/author** already landed (P1+P2, 2026-07-07) — optional
  top-level anchors that convert into native losslessly.
- **Enumeration survives `→ native` by answer reconstruction.** Native encodes
  the enumeration in the answer's display form, so an answer-less source word
  (ipuz/Exolve) has its answer rebuilt from the grid letters split per the
  enumeration (`,` → space, `-` → hyphen) — `(5,5)` over `OMEGAPOINT` regains
  `"OMEGA POINT"`. On enum/grid mismatch the grid wins (the enum re-derives as
  `(n)`, same principle as the Exolve disambiguator, §5.1); separators native
  answers cannot encode (`.`, `'`, space) are a **structural failure**.
- **Round-trips are *puzzle*-lossless, not *payload*-lossless.**
  native→ipuz→native preserves the whole puzzle but drops `link`/`meta` at the
  ipuz hop; native→Exolve→native likewise. Payload-lossless round-trips need a
  format that carries the metadata (today only native→native qualifies).
- **The three former v1 serializer coverage gaps are now closed** — they each
  emit the target's native representation instead of fail-striking as a `✗`
  would (the D7 rule that a serializer gap behaves like a `✗` cell still governs
  anything genuinely unmappable):
  - **prefilled cells → ipuz** — a given cell serializes to a puzzle-cell
    `value` and round-trips back (the ipuz parser reads that `value` as a
    prefilled letter).
  - **rebus → Exolve** — a multi-char cell flips the whole grid into
    `exolve-option: rebus-cells` mode (space-separated tokens) and round-trips
    back through the rebus-mode tokeniser; a grid with no rebus cell stays in
    the compact non-spaced form byte-for-byte (D6).
  - **shaded/styled cells → Exolve** (**serialize-only**) — a background-shade
    `color` style emits a top-level `exolve-colour` directive. Exolve's colour
    model is strictly coarser than the ipuz StyleSpec, so `exolve→ipuz` cannot
    reconstruct the exact StyleSpec — that direction is **puzzle-lossless only**.
    A shade is *structure* not metadata, so it never warns: the serializer maps
    the `color` shade or fails-strict on an unmappable style key (text/border
    colour, `imagebg`, `mark`, boolean `highlight`, … — a narrow, not a removal).

**Escape hatches — resolved (D9, 2026-07-07).** The two former deferred hatches
are decided: (a) **best-effort *structural* conversion** (crop/flatten with
warnings) is **rejected** — dropping structure yields a different or invalid
puzzle, which the D7 strict rule exists to prevent; (b) **native-model uplift**
is **adopted on lossless-only terms** — the schema grows where the engine can
compute-or-carry, additive-optionally. Landed: `title`/`author` (P1 engine + P2
xword — now ✓ for native above). Coming under the same rule: closed-subset cell
styling (P3), the `'`/`.` enumeration separators (P4). Staying **permanent**
fail-strict: rectangular grids and rebus (uplift would break native's
square-N×N single-letter model — ipuz is their carrier). Bars stay **gated** on
[`cryptic-layout-spec.md`](cryptic-layout-spec.md). Engine byte-parity: P1
retired the invented ipuz title and the constant `exolve-id`; the Exolve header
was then aligned to lead with title/setter (Stretch), so Exolve currently
matches the engine byte-for-byte and ipuz differs only in JSON serializer
formatting. **Byte-level parity is best-effort and not claimed either way** —
the guaranteed contract is *structural* agreement (§11, §14).

## 11. Testing

- **Parser round-trips** — `native → Board → native` is identity where
  lossless; same for ipuz, Exolve.
- **`convert` fidelity** — assert *structural failures* (rectangular/rebus/bars
  → native error, each naming the property) **and** *metadata drops*
  (`link`-bearing native → ipuz/Exolve succeeds, drops `link`, warns on stderr;
  `-q` silences the warning but stdout is byte-identical either way). Since the
  native-schema uplift (D9) **title/author no longer fail-strict** — assert they
  round-trip *into* native instead (a title-less source that gains the Q5 Exolve
  default carries `title: "Untitled"` back through native: gains the ecosystem
  default, loses nothing). Assert puzzle-lossless round-trips over an
  intersection fixture. **JSON output rule (S3, D6):** serialize with
  `json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"` —
  deterministic, idempotent (`dump==redump`), and **key-order-identical to the
  engine** (whose output is already fully key-sorted at every depth), so
  `convert`'s native/ipuz output matches the engine's key order. Full ipuz
  byte-identity is **not** a goal: the engine's SWI `json_write_dict` formatting
  differs from `json.dumps`, and matching it would abandon this clean-JSON rule
  — byte-level parity is best-effort only (§10, §14).
- **Engine cross-check (structural + anchor-field parity)** — the engine's
  `export` drops `link` and injects the Exolve `Untitled` default (P1), so v1
  asserts that `xword`'s ipuz/Exolve *structurally* agrees with
  `crosswordsmith export` (grid, numbering, clues, enumeration) — not full
  byte-identity. After P1 retired the invented ipuz title and the constant
  `exolve-id`, **title/author agree on both sides**: ipuz carries neither for a
  title-less layout, Exolve carries the byte-identical `  exolve-title: Untitled`
  line (both sides), and neither emits `exolve-id`. **Byte-level parity is
  best-effort, not asserted as a contract** (§14): the test checks structural
  agreement plus a few order-agnostic line checks, deliberately *not* a
  whole-file byte-equality (which would couple the suite to the engine's exact
  output). In practice Exolve currently matches byte-for-byte (header aligned to
  lead with title/setter); ipuz differs only in JSON serializer formatting.
- **`view`** — render to a string and compare golden snapshots (solved and
  `--blank`). **Exact capture recipe (S3):** `Console(width=<pinned>,
  no_color=True, force_terminal=False, record=True, file=io.StringIO(),
  legacy_windows=False)` → `print(renderable)` → `export_text(styles=False)`.
  The **pinned `width` is required** — an unpinned width leaks `COLUMNS` and
  de-stabilises the golden (clue lists wrap to console width); use one
  project-wide constant (S2/S3 used 200) so grids never wrap. Capturing to
  `io.StringIO()` (not real stdout) makes the golden invariant under
  TTY/pipe/`NO_COLOR`/`TERM`/`LANG`. An ANSI-preserving variant uses
  `force_terminal=True, color_system="truecolor"` + `export_text(styles=True)`.
- **`render`** — golden SVG/HTML (deterministic, so byte-comparable). **SVG-id
  discipline (S3, D6):** element ids are pure functions of grid coordinates
  (`cell-r{r}-c{c}`, `word-{n}-{dir}`); `uuid`/`datetime`/counters are banned
  from the serializer, and a one-line lint (grep output for `uuid`/`urn:`/4-digit
  years) guards the rule. PNG/PDF
  tested lightly (produced, correct dimensions/magic), not byte-golden —
  *because* the master SVG's `<text>` glyphs render against the **host font**
  (S6), so raster bytes are deterministic per-machine but not guaranteed
  cross-machine. A future cross-machine byte-golden would need text-as-`<path>`
  or a bundled/pinned font (§14).
- Fixtures reuse the engine's (`bundled_17`, a stock-grid mask, an ipuz and an
  Exolve sample). `pytest`.

## 12. Dependencies, layout & distribution

- **Core deps:** `rich` (rendering) + **`cyclopts`** (CLI, D8) + stdlib
  (`json`). **`xword[raster]`** extra = **`cairosvg`** (SVG→PNG *and* SVG→PDF),
  fixed by spike **S6**. Chosen over `resvg-py` (emits PNG only; PDF is required
  by D4/§6.3) and a headless browser (weight + non-deterministic PDF metadata).
  cairosvg output is byte-deterministic (no dates/ids to scrub) and its PDF is
  true vector. Its **one non-pip prerequisite is the system cairo library**
  (`cairocffi` dlopens `libcairo.so.2` at import) — `apt install libcairo2` or
  the distro equivalent; present on most dev machines, an extra line on minimal
  CI/containers. **Textual** joins later for the TUI.
- **Repo layout:** a self-contained Python package under **`xword/`** at the
  repo root (own `pyproject.toml`, `tests/`), isolating the Python toolchain
  from the Prolog tree. Use **src-layout** (`xword/src/xword/…`, hatchling with
  `[tool.hatch.build.targets.wheel] packages=["src/xword"]`); the console entry
  point is `xword = "xword.cli:main"` — a `main()` wrapper (**not** `:app`), so
  domain errors map to a non-zero exit (§8). The `[raster]` extra is
  `[project.optional-dependencies].raster`, exercised via `uv run --extra raster`
  / `pipx install 'xword[raster]'` (S1).
- **Distribution:** installed/run via **`uv` / `uvx` / `pipx`** so it stays a
  drop-on-`PATH` command with deps materialised on first run (`uv sync` for
  editable dev, `uv run xword …`, `uvx --from . xword …` for ephemeral,
  `uv tool install .` / `pipx install .` for a PATH command — all confirmed in
  S1). (The decision-record's single-file PEP 723 idea suits an initial spike
  only; a three-verb tool with parsers/renderers ships as a package.)

## 13. Phasing

| Phase | Deliverable |
|---|---|
| **0 (spikes)** | Ephemeral technical spikes to resolve tool choices + surface gotchas; **findings** (not code) fold back into concrete stories + spec corrections. **Gates Phase 1.** Seed set S1–S6 + backlog in [`xword-spikes.md`](xword-spikes.md); agents may add gating spikes. |
| **1 (MVP)** | `Board` + parsers (native/ipuz/Exolve) + `view` (solved default, `--blank`), grid + clue lists, `--from`, detection, Unix conventions. |
| **2** | `convert` any→any (D7): structural failures block; metadata drops warn on stderr (`-q`); structural engine cross-check. |
| **3** | `render` to SVG + HTML; PNG/PDF behind `xword[raster]`. |
| **4** | Interactive **Textual** renderer over the same `Board`. |
| **later** | **Native-model uplift** (D9) — additive-optional riders (cell styling P3, enumeration separators P4) after title/author (P1+P2). Best-effort *structural* conversion is **rejected** (D9). Engine byte-parity is best-effort (§14), not a phase deliverable. |

A `STATUS`-style tracker — [`xword-status.md`](xword-status.md) — is the system
of record for build progress (Phase 0 done; Phase 1 staged). This spec's phasing
is the plan; the tracker says where we are.

## 14. Open questions

- ~~**Colour scheme & cell geometry** for the terminal grid~~ — **resolved by
  S2** (see §6.1): 3w×1h cell interior, single-line box borders (pitch 4×2 ⇒
  square), `█` blocks, superscript inline corner numbers, `bold cyan`/default/
  `grey42` colour scheme gated on TTY-and-not-`NO_COLOR`.
- ~~**HTML styling surface** — a single built-in stylesheet vs a `--css`
  hook~~ — **resolved (Phase 3)**: one built-in inline `<style>`, so the HTML
  page is self-contained and deterministic (D6). A `--css` hook stays a cheap
  future add, not v1.
- ~~**SVG glyphs for the raster/PDF path** (S6) — `<text>` vs `<path>`~~ —
  **resolved (Phase 3)**: `<text>`. cairosvg already embeds glyph outlines into
  the PDF (S6), so the PDF is print-quality regardless; `<path>` glyphs (and the
  cross-machine raster byte-goldens they would enable) stay deferred (§11). The
  raster path is therefore tested by dimensions/magic, not byte-identity.
- **Rectangular native** — confirmed **permanent** hard-error (structural,
  D7/D9); the native-schema uplift (2026-07-07) **rejected** the
  square-padding/uplift option outright — ipuz is the rectangular carrier.
  Decided, not deferred.
- ~~**Engine byte-parity** — reachable once we decide whether `xword` should
  match the engine's *invented* default title on **both** ipuz (`"Untitled"`)
  and Exolve (`exolve-title: Untitled`) (currently: no, invent nothing)~~ —
  **resolved (2026-07-07)**: match on **Exolve only**. A title-less board
  targeting Exolve gains the engine's default at the `convert` boundary
  (warned, `-q`-silenceable; the serializers stay invent-nothing so library
  round-trips are untouched). Rationale: `docs/exet-verification.md` records
  Exet's *Save* crashing on a null title — the same evidence that made the
  engine inject its default — so this is a format-ecosystem requirement, not
  parity-chasing; it also closes the "Exet re-save of `xword`'s title-less
  Exolve" flag once considered for Phase 2. ipuz stays invent-nothing (title
  optional there, no breakage evidence). **Update (native-schema uplift P1+P2,
  2026-07-07, D9):** the engine itself dropped the invented ipuz `"Untitled"`
  and the constant `exolve-id`, and native now carries title/author, so
  byte-parity is no longer gated on either. **Byte-parity endgame (Stretch,
  2026-07-07) — closed as best-effort:** `xword`'s Exolve header was reordered
  to lead with title/setter (like `export.pl`), so Exolve output currently
  matches the engine byte-for-byte. For ipuz, a full-payload probe brought the
  two to within **7 of 697 lines** using engine options `[step(2), tab(80),
  width(1)]` + Python `separators=(",", ":")`; the residual is SWI
  `json_write_dict`'s *type-dependent* colon-spacing (`": "` before object/array
  values, `":"` before scalars), which no `json.dumps` flag reproduces —
  closing it would need a bespoke `xword` encoder that abandons its own
  clean-JSON rule (§11, D6) and diverges ipuz from native. **Decision:
  byte-level parity is best-effort, and no compatibility claim is made either
  way; the guaranteed contract is *structural* parity (§11).** Not chased
  further absent a consumer that needs identical bytes.
- **Textual scope** (candidate cycling, `--watch`, clue panes) — deferred to
  Phase 4 design.
