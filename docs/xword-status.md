# `xword` — Build status

Progress tracker for the work specified in [`xword-spec.md`](./xword-spec.md).
This is the **system of record** for what is done / in progress / blocked in the
`xword` tool. The spec says *what*; [`xword-spikes.md`](./xword-spikes.md) recorded
the Phase-0 *findings* (tool choices + gotchas); this file says *where we are*.

`xword` is a **standalone Python toolchain** — a downstream consumer of the
crosswordsmith engine, not part of its four-verb design-spec. It has its own
`docs/STATUS.md`-style tracker (this file) so the Prolog engine's status stays
separate (spec §12).

**Update discipline:** move a row's status in the **same commit/PR** as the work.
A component is `done` only when its referenced spec section is satisfied **and**
its tests pass. Don't add rows for unspecified scope — spec it first, then track
it here. `xword-spec.md` has no formal `AC-*` codes (unlike the engine spec); the
"Verified by" column names the concrete test/golden that pins each row.

**Status legend:**

| Mark | Meaning |
|---|---|
| `not started` | Specified, nothing built yet |
| `in progress` | Actively being built |
| `done` | Built; spec section satisfied + tests pass |
| `blocked` ⊘ | Waiting on an open decision (see Open questions) |
| `deferred` | Spec'd as later-phase — not in the current phase's scope |

---

## Phase 0 — spikes (spec §13) — **done (2026-07-03)**

All six seed spikes complete and gating cleared; findings + working snippets in
[`xword-spikes.md`](./xword-spikes.md), spec corrections folded into
`xword-spec.md`. **No new gating spikes were filed**, so Phase 1 is unblocked.

| Spike | Outcome (one line) | Gates |
|---|---|---|
| S1 CLI skeleton | cyclopts 4.20, **src-layout**, entry `xword.cli:main`; bare-`xword` non-zero needs `@app.default` | P1 · done |
| S2 grid render | manual `Text`/`Segment` (not `Table`); 3w×1h cells, `█` blocks, superscript numbers; resolves §14 geometry | P1 · done |
| S3 determinism | pinned-width `Console`+`StringIO`+`export_text` recipe; `json.dumps(sort_keys=True…)` = engine key order | P1 · done |
| S4 ipuz | §10 table holds; full native↔ipuz mapping; structural parity with engine `export --to ipuz` | P1 · done |
| S5 exolve | §10 holds; normative enumeration parse-back rule (end-anchored regex + cell-count check); decorator suffixes | P1 · done |
| S6 raster | **`cairosvg`** chosen (PNG+PDF, deterministic); one prereq = system cairo | P3 · done |

---

## Phase 1 — MVP: `Board` + parsers + `view` (spec §13) — **done (2026-07-03)**

**Deliverable:** the `Board` model + parsers (native / ipuz / Exolve) + the `view`
command (solved default, `--blank`), grid + Across/Down clue lists, `--from`,
format detection (§7), Unix conventions (§8). NOT `convert`/`render`.

| # | Component | Spec | Status | Verified by |
|---|---|---|---|---|
| 1 | **Package skeleton** — `xword/` src-layout, `pyproject.toml`, `xword.cli:main` entry point, `[raster]` extra stub, `uv`/`pytest` dev-run (S1) | §12, D8 | done | `uv run xword --help`; `uv run pytest` green (57 tests) |
| 2 | **CLI shell (cyclopts)** — `view`/`convert`/`render` verbs (convert/render stubbed), `--from`/`--to` enums, `--blank`/`-q`, per-verb `--help`, bare-`xword`→usage+non-zero, stdin default / positional override / `--out` on success (S1) | §6, §8, D1 | done | `tests/test_cli.py` (subprocess exit codes + stream split) |
| 3 | **`Board` model** — union-superset dataclasses; derived row-major numbering; enumeration normalisation; "invent nothing" (S2/S4/S5) | §5.1 | done | `tests/test_board.py` (numbering + enumeration units) |
| 4 | **native parser** — `native JSON → Board` + `Board → native` | §5.1, json-output-spec §6 | done | `tests/test_formats.py::TestNative` (round-trip identity on `bundled_17`) |
| 5 | **ipuz parser** — `ipuz → Board` + `Board → ipuz`; read declared `empty`/`block`; object-form styled cells; prefilled `value` (given letter) read+emit; 3 clue-item shapes (S4) | §5.1, §10 | done | `tests/test_formats.py::TestIpuz` (round-trip is structurally identical to engine `export`) |
| 6 | **Exolve parser** — `Exolve → Board` + `Board → Exolve`; tokenised grid rows + decorator suffixes; `exolve-option: rebus-cells` mode (space-separated multi-char tokens) read+emit; `exolve-colour` shade emit (serialize-only); enumeration parse-back rule (S5) | §5.1, §10 | done | `tests/test_formats.py::TestExolve` + enum units in `test_board.py` |
| 7 | **Format detection ladder** — JSON→ipuz-first (`ipuz.org` in `version`/`kind`)→native; else Exolve markers; `--from` overrides; clear no-match error (S4/S5) | §7 | done | `tests/test_detect.py` (all three fixtures + garbage) |
| 8 | **`view` renderer** — `BoardGeom → Rich renderable` (manual `Text`); grid + clue lists; solved default + `--blank`; `soft_wrap` never-wraps; colour gating `isatty() and not NO_COLOR` (S2) | §5.2, §6.1 | done | `tests/test_view.py` goldens (`tests/golden/view_bundled_17_{solved,blank}.txt`, S3 recipe) + never-wraps test |
| 9 | **Unix conventions** — diagnostics→stderr / data→stdout; `NO_COLOR`/TTY; `--out` atomic-on-success; deterministic output (D6) | §8, D6 | done | `tests/test_cli.py` (pipe/no-ANSI, `--out` on success only, stderr-only errors) |

Fixtures regenerated from the engine (`crosswordsmith arrange --input
fixtures/bundled_17_clues.pl --size 17`, then `crosswordsmith
export --to …`) into `xword/tests/fixtures/`, plus the recreated S4/S5 keepsakes
(`sample.ipuz.json` — block/circle/rebus/title/both-clue-shapes;
`sample_decorated.exolve` — block/unfilled/circle/bars/enum-collision clue).
`pytest` under `xword/tests`.

---

## Phase 2 — `convert` (spec §13) — **done (2026-07-03)**

**Deliverable:** `xword convert --to native|ipuz|exolve` any→any (D7): structural
failures block with the property + a capable target named; metadata drops warn
on stderr (`-q` silences, stdout untouched); structural engine cross-check.

| # | Component | Spec | Status | Verified by |
|---|---|---|---|---|
| 1 | **`convert` verb wiring** — parse → fidelity → serialize; `--to` required, `--from` optional (detect), stdin/stdout + `--out` on success, `-q` | §6.2, §8 | done | `tests/test_convert.py::TestCli` |
| 2 | **D7 strict structure** — `→ native` fails on rectangular / circled / barred / prefilled / styled / unfilled / rebus / unencodable enum, each naming a capable target; the three former v1 serializer gaps are now closed — rebus→exolve (whole-grid `rebus-cells` mode, round-trip), prefill→ipuz (puzzle `value`, round-trip), styling→exolve (background `color` as `exolve-colour`, serialize-only). **title/author dropped off the fail list** — the native-schema uplift (D9, 2026-07-07) made them optional top-level anchors native parses+serializes | §10, D7, D9 | done | `tests/test_convert.py::TestStructuralFailures` + `TestPrefilledToIpuz` + `TestShadedToExolve` + `TestRebusToExolve` |
| 3 | **D7 metadata drop-and-warn** — per-word `link`/`meta` keys + top-level `diagnostics` dropped on `→ ipuz/exolve` with one stderr warning each; `-q` silences, stdout byte-identical; native target drops nothing | §6.2, §10 | done | `tests/test_convert.py::TestMetadataDrops` + `TestCli` |
| 4 | **`diagnostics` passthrough** — carried on `Board`, native→native re-emits verbatim (payload-lossless identity, json-output-spec §6.4) | §6.2, §10 | done | fixture `arrange_bundled_17_fixed.json` (engine golden verbatim) in `TestMetadataDrops`/`TestCli` |
| 5 | **Puzzle-lossless round-trips** — native→ipuz→native and native→exolve→native = source minus `link`; display answers reconstructed from enumeration (`(5,5)` → `"OMEGA POINT"`), grid wins on enum mismatch; exolve→ipuz→exolve keeps bars/circles/unfilled (barred StyleSpec synthesized) + prefill (`!`↔`value`) + rebus (multi-char↔`rebus-cells`). Since D9 **title/author survive the round-trip** too: a title-less exolve source that gained the Q5 `Untitled` default carries it back through native rather than fail-stricting the way home | §10, D9 | done | `tests/test_convert.py::TestRoundTrips` + `test_native_exolve_native` + `test_titled_exolve_native` + `TestPrefilledToIpuz` + `TestRebusToExolve` |
| 6 | **Engine cross-check (structural + anchor parity)** — convert native→ipuz/exolve agrees with `crosswordsmith export` on grid/numbering/clues/enumeration; the exolve title line is byte-identical (both `exolve-title: Untitled`). After P1 (D9) retired the invented ipuz title **and** the constant `exolve-id`, title/author agree on both sides and neither emits `exolve-id`; residual gaps are cosmetic — ipuz JSON whitespace + the exolve title line's header position (engine after `exolve-begin`, xword after width/height) — for the byte-parity Stretch | §11 | done | `tests/test_convert.py::TestEngineCrossCheck` |

---

## Phase 3 — `render` (spec §13) — **done (2026-07-03)**

**Deliverable:** `xword render --to svg|html|png|pdf` over the shared `BoardGeom`
(spec §5.2) — an **SVG** master (title + grid + two-column clue lists), a
semantic **HTML** page (inline SVG grid + HTML clue lists, one built-in
stylesheet), and **PNG/PDF** rasterised *from that SVG* behind the
`xword[raster]` extra (`cairosvg`, S6). Deterministic (D6): coordinate-derived
ids, no uuid/timestamps. Binary-to-TTY guard on png/pdf.

| # | Component | Spec | Status | Verified by |
|---|---|---|---|---|
| 1 | **`render` verb wiring** — `--to svg\|html\|png\|pdf` required, `--from` optional (detect), `--blank`, stdin/stdout + `--out` on success; text vs binary write | §6.3, §8 | done | `tests/test_cli.py` (svg/html→stdout, `--out` on success, png `--out` binary) |
| 2 | **SVG renderer** — `master_svg` (standalone: title + grid + clue lists) + `grid_svg` (grid-only, embedded by HTML); `<text>` glyphs (Q3); circle/bar/rebus; coord-derived ids | §5.2, §6.3, D6 | done | `tests/test_render.py` goldens (`render_bundled_17_{solved,blank}.svg`) + id / decoration / determinism units |
| 3 | **HTML renderer** — self-contained page: one built-in inline `<style>` (Q2), inline `grid_svg`, semantic `<ol>` clue lists with per-clue `value` + `word-{n}-{dir}` id | §6.3, §14 Q2 | done | `tests/test_render.py` goldens (`render_bundled_17_{solved,blank}.html`) |
| 4 | **Raster backend** — `to_png`/`to_pdf` convert the master SVG via **lazy** cairosvg; the two failure modes reported distinctly (extra missing vs system cairo missing) | §6.3, §8, §12 | done | `tests/test_render.py` (PNG magic + IHDR width; PDF magic; monkeypatched `ImportError`/`OSError`) |
| 5 | **Binary-to-TTY guard** — png/pdf refuse a terminal stdout (checked *before* any raster work); require `--out`/redirect | §6.3 | done | `tests/test_cli.py::test_render_png_refuses_terminal` (real pty) |
| 6 | **Determinism (D6)** — coord-derived ids (`cell-r{r}-c{c}`, `word-{n}-{dir}`), no uuid/datetime; SVG/HTML byte-golden; raster dims/magic only (host-font `<text>`, §11) | D6, §11 | done | `tests/test_render.py` (goldens + the `uuid`/`urn:`/year lint, modulo the fixed SVG namespace URI) |

Geometry stayed the shared superset (§5.2): `GeomCell` gained `bar_right`/
`bar_below` (populated from `Cell`) so the SVG master draws bars — the terminal
renderer ignores them, unchanged. cairosvg 2.9.0 + system `libcairo.so.2` were
present, so the full PNG/PDF path is exercised (raster tests `importorskip`
cairosvg, so the base suite still runs without the extra).

---

## Later phases (spec §13) — deferred

| Phase | Deliverable | Status |
|---|---|---|
| 4 | Interactive **Textual** renderer over the same `Board` (S2 confirmed the seam) | deferred |
| — | **Native-model uplift** (D9) — additive-optional, lossless-only: `title`/`author` **landed** (P1 engine + P2 xword, 2026-07-07); closed-subset cell styling (P3) + `'`/`.` enumeration separators (P4) are the next riders, each with its own go/no-go | in progress (title/author done) |
| — | Best-effort **structural** conversion (crop/flatten + warnings) | **rejected (D9)** — yields a different/invalid puzzle |
| Stretch | Byte-parity endgame — close the last two cosmetic engine↔xword gaps (ipuz JSON whitespace; exolve title-line header position) | deferred (own go/no-go) |

---

## Open questions (spec §14) — these gate future rows

| # | Question | Status |
|---|---|---|
| Q1 | Colour scheme & cell geometry for the terminal grid | **resolved by S2** (§6.1) |
| Q2 | HTML styling surface — built-in stylesheet vs `--css` hook | **resolved (Phase 3)**: one built-in inline `<style>` (self-contained + deterministic, D6); a `--css` hook stays a cheap future add |
| Q3 | SVG glyphs for raster/PDF: `<text>` vs `<path>` | **resolved (Phase 3)**: `<text>` — cairosvg embeds glyph outlines in the PDF regardless (S6); `<path>` / cross-machine raster byte-goldens stay deferred (§11/§14) |
| Q4 | Rectangular native — confirmed hard-error; square-padding uplift is deferred best-effort | **resolved (permanent hard-error)**: the native-schema uplift (D9, 2026-07-07) **rejected** square-padding/uplift outright — ipuz is the rectangular carrier |
| Q5 | Engine byte-parity — whether `xword` matches the engine's *invented* default title (ipuz `"Untitled"` + `exolve-title: Untitled`) | **resolved (2026-07-07)**: match on **Exolve only** — title-less boards gain the engine's default at the `convert` boundary (warned, `-q`-silenceable; serializers stay invent-nothing) because Exet's Save crashes on a null title (§6.2/§14). ipuz stays bare. **Follow-on (native-schema uplift P1+P2, D9):** the engine dropped the invented ipuz title **and** `exolve-id`, so byte-parity is no longer gated on either — only the two cosmetic Stretch gaps remain (ipuz JSON whitespace; exolve title-line header position) |
| Q6 | Textual scope (candidate cycling, `--watch`, clue panes) | open (Phase 4 design) |

---

## At a glance

- **Phase 0 (spikes): done (2026-07-03)** — all S1–S6 findings folded into
  `xword-spikes.md`, spec corrections applied, no new gating spikes.
- **Phase 1 (MVP): done (2026-07-03)** — all 9 components built under `xword/`
  (src-layout package, `uv sync` / `uv run xword` / `uv run pytest`); 57 tests
  green. Notable outcome: `xword`'s ipuz serialization is *structurally
  identical* to the engine's `export --to ipuz` on `bundled_17` (the §11
  cross-check holds already at the parser level).
- **Phase 2 (`convert`): done (2026-07-03)** — any→any among the three formats
  with the D7 policy; 83 tests green. Two flagged judgment calls resolved and
  spec'd (§6.2): native→native **preserves** `diagnostics` (payload-lossless
  identity), and no `exolve-id` is emitted (Exolve auto-derives one from a
  grid/clue signature, so nothing downstream requires it).
- **Phase 3 (`render`): done (2026-07-03)** — `render --to svg|html|png|pdf`
  over the shared `BoardGeom`; **101 tests green** (base suite `importorskip`s
  cairosvg so it runs without the extra; the full raster path was exercised
  here). SVG/HTML are deterministic byte-goldens; PNG/PDF are tested by
  dims/magic (host-font `<text>`, §11). The two Phase-3 open questions are
  resolved and spec'd: **Q2** → one built-in inline stylesheet (self-contained,
  D6); **Q3** → `<text>` glyphs (cairosvg embeds outlines in the PDF anyway).
- **Native-schema uplift (D9): title/author landed (P1+P2, 2026-07-07).** Both
  the engine and `xword` now carry optional top-level `title`/`author`
  losslessly across native/ipuz/exolve; the invented ipuz title and the constant
  `exolve-id` are retired. `xword` suite green (114 tests). Next riders under the
  same lossless-only rule: closed-subset cell styling (P3) and `'`/`.`
  enumeration separators (P4), each with its own go/no-go.
- **Phase 4 (Textual TUI) is next**; Phase 4 + later: deferred, not yet started.
- **Nothing blocked; the native-model uplift is the one thing in progress**
  (title/author done, styling/enum riders pending). §14 open questions resolved
  so far: Q1 (grid geometry), Q2 (HTML styling), Q3 (SVG glyphs), Q4 (rectangular
  native — permanent hard-error, uplift rejected), Q5 (engine title: Exolve-only
  at the convert boundary, plus P1/P2 retiring the ipuz title & `exolve-id`).
  Q6 gates Phase 4.
