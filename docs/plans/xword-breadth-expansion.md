# xword breadth expansion — proposal

**Status:** proposal, 2026-07-15. Companion to
[`docs/research/setter-tool-landscape-2026.md`](../research/setter-tool-landscape-2026.md)
§B/§F. Grounded in the current `xword/src/xword/` source and
[`docs/xword-spec.md`](../xword-spec.md). Nothing here is built yet.

## The strategic frame (from the landscape research)

The 2026 research verdict for the `xword` companion is explicit: **format
conversion is table stakes, not a differentiator** — puzpy, kotwords,
xd-crossword-tools, and Exolve's own scripts already own the interchange graph
with published round-trip fidelity. So breadth investment should go where `xword`
is actually differentiated:

1. **Rendering quality/targets** (the less-crowded half of its surface), and
2. **Being the faithful, deterministic multitool over the native-JSON hub** —
   read-only analysis/inspection verbs that no convert-only library offers.

Do **not** spend effort chasing convert *coverage* or engine↔xword byte-parity
(spec §14 already closes byte-parity as best-effort). Keep the `.puz`/`.jpz`/PDF
*binary write* boundary delegated (kotwords for the engine; puzpy for reads) —
never hand-roll a CRC-16 `.puz` writer.

## What xword already does (so we don't re-propose it)

Three verbs (`cli.py`): `view` (terminal), `convert` (native/ipuz/Exolve any→any,
D7 fidelity policy), `render` (svg/html/png/pdf off one master SVG). Reads
native/ipuz/Exolve with auto-detection (`detect.py`; binary-magic branch reserved
but inert). Already covered: clue-list-alongside-grid, solution-vs-blank
(`--blank`), stdin→stdout piping, `NO_COLOR`/TTY guards, `--out`-on-success,
metadata round-trip fidelity. **Do not** re-file these as gaps.

## Shortlist — "double down on breadth" (top 5)

Tagged effort × leverage × identity-fit.

1. **`stats` / `inspect` verb** — parse any supported format → deterministic
   descriptive metrics (dimensions, block/word counts, checked/unchecked cells,
   length histogram, enumeration distribution, symmetry, connectivity, solved vs
   blank). Human table default; `--json` for machines.
   *Effort low · leverage high · fit strong.* Reuses `Board`/`derive_words`;
   read-only; amplifies the native-hub differentiator; extends inspection to
   ipuz/Exolve the **engine** can't ingest. Keep it *descriptive* metrics only —
   no quality scoring (that's the engine's lane; avoid a second source of truth
   vs engine `lint`).

2. **Print/PDF rendering breadth** — page-oriented PDF (Letter/A4, margins),
   `--grid-only` / `--clues-only` splits, combined blank+solution sheet, optional
   `--css` hook (spec Q2). *Effort medium · leverage medium-high · fit strong.*
   The research-endorsed differentiated half ("compete on render, not convert");
   cairosvg output stays byte-deterministic.

3. **`.puz` *read* via an optional `puzpy` extra** — `.puz` → `Board` (puzpy is
   pure-Python **MIT**), as an `xword[puz]` extra mirroring the `[raster]` pattern,
   plus wiring the reserved binary-magic detect branch. *Effort medium · leverage
   medium · fit ok.* Highest-priority real-world format gap (xword can't read a
   `.puz` at all today). Spec §9 already sanctions "one more parser via puzpy."
   The kotwords boundary was about the engine's *outbound binary write*, not
   inbound reads — MIT puzpy respects license-cleanliness and hand-rolls no CRC.
   `.puz` can't hold bars, so barred boards fail-strict on `.puz` *write* (D7).

4. **`diff` verb** — parse two layouts (any formats) → structural diff (grid
   shape, blocks, per-cell letters, numbering, clues, enumerations); human default,
   `--json`, exit non-zero on difference. *Effort medium · leverage medium · fit
   strong.* No surveyed tool offers cross-format semantic diff; serves the
   git-diffable/deterministic identity ("did this convert change the puzzle?").

5. **`.xd` read/write** — the one text-interchange format xword doesn't touch
   (git-diffable, Puzzmo ecosystem; xd-crossword-tools' format). *Effort low ·
   leverage low-med · fit strong.* Research rates leverage low ("only if `.xd`
   matters to you") but it's the cheapest genuine *format*-breadth add and fits the
   diffable ethos. STOP condition: if `.xd`'s quirks fight the `Board` model,
   document the gap and stop.

## Backlog (lower priority / higher strain)

- **Themeable `--css` HTML hook** (spec Q2) — subsumed into #2 if pursued.
- **Alternate terminal styles** — plain-ASCII fallback for non-Unicode terminals;
  low effort/low leverage (fixed box geometry is a deliberate spec choice).
- **Batch/glob** (`convert --to ipuz *.json --out-dir DIR`) — mild strain: current
  model is strictly one-in-one-out; needs an output-naming scheme.
- **`.jpz` read** — XML; kotwords already reads it; weaker than `.puz`.
- **`.puz`/`.jpz`/PDF *write* via kotwords offload** — shell-out handoff; lower
  priority than puzpy `.puz` read.
- **Format `validate` verb** — *strains identity* (edges toward a setter tool /
  overlaps engine `lint`); if built, scope strictly to format well-formedness and
  fold into #1.

## Do NOT pursue

- Convert *coverage* or ipuz/Exolve byte-parity with the engine — table stakes.
- A hand-rolled `.puz` CRC-16 writer — use puzpy/kotwords.
- Any clue generation, solving, or fill-quality scoring in xword — the engine's
  lane and the hard identity line ("feeds clue tools; is not one").

## Recommended first cut

`stats`/`inspect` + `diff` (both cheap, deterministic, read-only, hub-amplifying)
and **print/PDF rendering breadth** — then `.puz` read via a puzpy extra and
`.xd` I/O for the two format gaps that actually matter. Explicitly *not*
competing on convert coverage or byte-parity.
