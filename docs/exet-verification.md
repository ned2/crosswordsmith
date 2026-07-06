# AC-EXP-2 — Exolve ↔ Exet round-trip verification

**Status: manual procedure.** `AC-EXP-2` ("`export --to exolve` produces text
that round-trips through Exet: load → save → equivalent grid+entries") cannot be
checked in-repo: [Exet](https://viresh-ratnakar.github.io/exet.html) is a
third-party browser app and there is no Exolve *parser* in this codebase (an
in-program emit-then-reparse would only prove self-consistency, not third-party
ingestion). The automated coverage we *do* have is a byte-exact golden
(`tests/golden/export_bundled_17.exolve`) that pins the emitter against
regression. This checklist is the manual complement; run it when the Exolve
emitter changes, and record the outcome in the log at the foot of this file.

Exet and Exolve are by Viresh Ratnakar; Exet's **Open** menu imports an Exolve /
`.puz` / `.ipuz` file, and its **Save** menu downloads the puzzle back as Exolve
— that is the round-trip under test.

## The contract — what "equivalent" means

The round-trip must preserve the **puzzle**, not the bytes. After Exet re-exports,
these must match the crosswordsmith output:

1. **Grid dimensions** (`exolve-width` / `exolve-height`).
2. **Block pattern** — every blocked cell (`.`) in the same position.
3. **Solution letters** — every light cell carries the same letter (crosswordsmith
   emits a fully-filled solution grid).
4. **Clue numbering** — Exet re-derives numbers from the grid; they must equal
   crosswordsmith's across/down numbers.
5. **Entry set** — the same across entries and the same down entries (same cells ⇒
   same answers).
6. **Enumerations** — each clue's `(5,5)` / `(4-5)` enumeration is preserved.
7. **Clue text** — preserved verbatim for every entry that had a clue; entries with
   no clue stay blank (Exet/crosswordsmith must not invent text — AC-EXP-3).

## Automated load-half check (exolve engine) — reproducible

`tests/exolve_ingest_check.sh` covers the **load** half of the round-trip
automatically. It exports a layout to Exolve, then has the **real Exolve engine**
(`exolve-m.js` — the parser Exet is built on, by the same author) ingest it in a
headless browser and confirms it reconstructs the source layout's *exact* grid,
block pattern, solution letters, entry set, enumerations, and clue text:

```sh
tests/exolve_ingest_check.sh [layout.json]   # needs Chrome/Chromium + curl + python3
```

This is the substantive risk in AC-EXP-2 — *does a third-party tool read our
Exolve correctly?* — and it is verified. It does **not** cover the **save-back**
half (Exet UI Open → Save → download), which stays the manual procedure below;
that half is low-risk (Exet re-serialises its own exolve-parsed model to standard
Exolve), so the manual run is a confirmation, not the primary evidence. Optional /
not part of `make test`.

## Prerequisites

- A current build of `crosswordsmith` in this repo.
- A browser; open <https://viresh-ratnakar.github.io/exet.html> (all data stays
  local — Exet sends nothing to a server). Note the Exet version (shown in-app) for
  the log.

## Procedure

### 1 — Produce the Exolve artifact

Any canonical layout works. The pinned golden's source layout is convenient:

```sh
./crosswordsmith export --to exolve \
    tests/golden/arrange_bundled_17_fixed.json --out /tmp/cs.exolve
```

> **Symmetry note.** `arrange` (Flavour A) layouts are often **asymmetric**
> (no symmetry guarantee). Exet can open them, but if you want a grid in Exet's
> house style (180°-symmetric blocked), export a **`fill`** output on a stock grid
> instead — e.g. `./crosswordsmith fill --grid grids/blocked_15a.json --dict
> UKACD18.txt --out /tmp/fill.json` then `export --to exolve /tmp/fill.json --out
> /tmp/cs.exolve`. An asymmetry warning from Exet is an Exet construction-mode
> nicety, **not** a round-trip failure.

Keep `/tmp/cs.exolve` as the reference.

### 2 — Load into Exet

In Exet: **Open** → open `/tmp/cs.exolve`. The grid, solution letters, clue
numbers, clues, and enumerations should appear. (crosswordsmith emits a bare
`exolve-begin … exolve-end` block; if your Exet build expects the block pasted
rather than uploaded, copy the file contents into Exet's Exolve input.)

### 3 — Save back from Exet

In Exet: **Save** → download as **Exolve** → `/tmp/exet.exolve`. If Exet wraps the
Exolve in an HTML file, extract just the `exolve-begin … exolve-end` block.

### 4 — Compare

```sh
diff /tmp/cs.exolve /tmp/exet.exolve   # orientation only; expect benign diffs
```

`diff` will show Exet's added metadata — that's expected (see below). Walk the
**equivalence checklist** instead of requiring byte-identity:

- [ ] width/height unchanged
- [ ] block (`.`) pattern identical
- [ ] every light cell's solution letter identical
- [ ] across/down clue numbers identical
- [ ] same across entries; same down entries
- [ ] every enumeration identical (`(5,5)`, `(7,6)`, …)
- [ ] every clue's text preserved; no clue text invented for blank entries

If all boxes hold, **AC-EXP-2 PASSES** for this Exet version.

## Acceptable differences (NOT failures)

Exet legitimately adds/normalizes things that don't affect puzzle equivalence:

- New metadata directives: `exolve-title`, `exolve-setter`/`-author`,
  `exolve-copyright`, `exolve-maker` ("created with Exet"), `exolve-version`, a
  changed `exolve-id`, `exolve-option …` lines.
- Whitespace / indentation / section-ordering normalization.
- HTML wrapping around the Exolve block.
- Annotations Exet may attach (prefill, colours, ninas) if you touched the grid in
  Exet — re-export without editing to keep the comparison clean.

## If it fails — triage

| Symptom | Likely cause / check |
|---|---|
| Light cells lost their letters | crosswordsmith emits a full solution grid; confirm the artifact has letters (not just `.`/`0`), and that Exet imported the solution, not a blank grid. |
| Enumerations dropped or wrong | crosswordsmith carries the *original* answer (with spaces/hyphens) and derives the enumeration explicitly (`answer_enumeration/2`); compare against Exet's. A multi-word answer with no separator can't enumerate — check the source answer. |
| Clue numbers misaligned | a grid-interpretation difference; recheck the `.`-as-block mapping (`exolve_cell_char/2`) and that the block pattern survived import. |
| Exet refuses / warns on the grid | usually asymmetry (Flavour A) — re-run with a symmetric `fill`-on-stock-grid artifact (see step 1 note), or enable Exet's asymmetry allowance. |

A genuine failure (a checklist box that can't be satisfied) is an **emitter** bug —
fix `layout_to_exolve/2` (or the enumeration deriver), update the golden, and note
it here.

## Verification log

Record each manual run so AC-EXP-2's status is auditable.

| Date | Engine / tool | crosswordsmith commit | Artifact | Result | Notes |
|---|---|---|---|---|---|
| 2026-07-01 | exolve engine (`exolve-m.js`, headless Chrome) | `aaf1e88` | `arrange_bundled_17_fixed` (17×17) + `arrange_toc_demo_max` (22×22) | **Load half: PASS** (automated) | The real Exolve engine ingests the export and reconstructs grid + block pattern + solution letters + entry set + enumerations + clue text *exactly*, via `tests/exolve_ingest_check.sh`. **Save-back half (Exet UI Open→Save) NOT run** — it needs a human at a browser; low-risk, still pending. |
| 2026-07-01 | Exet v1.06 (browser UI) | `aaf1e88` + export title fix | `arrange_bundled_17_fixed` → Exet → `.ipuz` + `.puz` | **Save-back: PASS** | Full Exet Open→Save→export; both re-exports match the source exactly (dimensions, solution grid, 4 across + 2 down entries, all enumerations, all clue text). **Note:** Exet's Save crashed on the null title (`fileTitle`/`updateSavePanel`) until a title was set — export now emits a default `exolve-title` (see revamp-audit V1). AC-EXP-2 fully verified (load + save-back). |

**Addendum (2026-07-07 — P1 native-schema uplift):** the default-title behavior
recorded above has narrowed. The default `exolve-title: Untitled` is now emitted
**only when the layout carries no title** (a title-less layout still gets it, so
the Exet Save defence this log established is intact); a layout that *does* carry
a top-level `title` now exports that title verbatim, and a top-level `author`
becomes `exolve-setter`. The **ipuz** side no longer invents a title at all
(previously `"title": "Untitled"`) — ipuz makes title optional and Exet's crash
was specific to the Exolve conversion, so the defence stays Exolve-only. Net:
the load + save-back PASS above is unaffected for the title-less
`arrange_bundled_17_fixed` artifact (byte-identical `exolve-title` line); the new
titled/authored path (`fixtures/titled_layout.json`) is golden-pinned but not yet
manually re-verified through the Exet UI — a low-risk follow-up if a titled
artifact is ever round-tripped. Ref: `docs/plans/native-schema-uplift-plan.md` P1.
