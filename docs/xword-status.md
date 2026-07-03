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
| 5 | **ipuz parser** — `ipuz → Board` + `Board → ipuz`; read declared `empty`/`block`; object-form styled cells; 3 clue-item shapes (S4) | §5.1, §10 | done | `tests/test_formats.py::TestIpuz` (round-trip is structurally identical to engine `export`) |
| 6 | **Exolve parser** — `Exolve → Board` + `Board → Exolve`; tokenised grid rows + decorator suffixes; enumeration parse-back rule (S5) | §5.1, §10 | done | `tests/test_formats.py::TestExolve` + enum units in `test_board.py` |
| 7 | **Format detection ladder** — JSON→ipuz-first (`ipuz.org` in `version`/`kind`)→native; else Exolve markers; `--from` overrides; clear no-match error (S4/S5) | §7 | done | `tests/test_detect.py` (all three fixtures + garbage) |
| 8 | **`view` renderer** — `BoardGeom → Rich renderable` (manual `Text`); grid + clue lists; solved default + `--blank`; `soft_wrap` never-wraps; colour gating `isatty() and not NO_COLOR` (S2) | §5.2, §6.1 | done | `tests/test_view.py` goldens (`tests/golden/view_bundled_17_{solved,blank}.txt`, S3 recipe) + never-wraps test |
| 9 | **Unix conventions** — diagnostics→stderr / data→stdout; `NO_COLOR`/TTY; `--out` atomic-on-success; deterministic output (D6) | §8, D6 | done | `tests/test_cli.py` (pipe/no-ANSI, `--out` on success only, stderr-only errors) |

Fixtures regenerated from the engine (`crosswordsmith arrange --input
fixtures/bundled_17_clues.pl --size 17 --size-mode fixed`, then `crosswordsmith
export --to …`) into `xword/tests/fixtures/`, plus the recreated S4/S5 keepsakes
(`sample.ipuz.json` — block/circle/rebus/title/both-clue-shapes;
`sample_decorated.exolve` — block/unfilled/circle/bars/enum-collision clue).
`pytest` under `xword/tests`.

---

## Later phases (spec §13) — deferred

| Phase | Deliverable | Status |
|---|---|---|
| 2 | `convert` any→any (D7): structural failures block; metadata drops warn on stderr (`-q`); structural engine cross-check | deferred |
| 3 | `render` to SVG + HTML; PNG/PDF behind `xword[raster]` (`cairosvg`, S6) | deferred |
| 4 | Interactive **Textual** renderer over the same `Board` (S2 confirmed the seam) | deferred |
| later | Best-effort **structural** conversion (crop/flatten + warnings) and/or native-model uplift; restores engine byte-parity | deferred |

---

## Open questions (spec §14) — these gate future rows

| # | Question | Status |
|---|---|---|
| Q1 | Colour scheme & cell geometry for the terminal grid | **resolved by S2** (§6.1) |
| Q2 | HTML styling surface — built-in stylesheet vs `--css` hook | open (Phase 3) |
| Q3 | SVG glyphs for raster/PDF: `<text>` vs `<path>` | open (Phase 3 sub-decision, S6) |
| Q4 | Rectangular native — confirmed hard-error; square-padding uplift is deferred best-effort | resolved (hard-error); uplift deferred |
| Q5 | Engine byte-parity — whether `xword` matches the engine's *invented* default title (ipuz `"Untitled"` + `exolve-title: Untitled`) | open (currently: no, invent nothing) |
| Q6 | Textual scope (candidate cycling, `--watch`, clue panes) | open (Phase 4 design) |

---

## At a glance

- **Phase 0 (spikes): done (2026-07-03)** — all S1–S6 findings folded into
  `xword-spikes.md`, spec corrections applied, no new gating spikes.
- **Phase 1 (MVP): done (2026-07-03)** — all 9 components built under `xword/`
  (src-layout package, `uv sync` / `uv run xword` / `uv run pytest`); 57 tests
  green. Notable outcome: `xword`'s ipuz serialization is *structurally
  identical* to the engine's `export --to ipuz` on `bundled_17` (the §11
  cross-check holds already at the parser level); `convert`/`render` are
  registered as stubs that exit 1 with a "not implemented (Phase 2/3)" error.
- **Phase 2 (`convert`) is next** — the parsers/serializers it needs already
  exist and are round-trip-tested; what remains is the verb wiring, the D7
  strict-structure/drop-metadata policy (+ `-q`), and the engine cross-check
  tests.
- **Phases 2–4 + later: deferred** — not yet started.
- **Nothing in progress; nothing blocked.** Q1 (grid geometry) is the only §14
  open question resolved so far; the rest gate Phases 3–4.
