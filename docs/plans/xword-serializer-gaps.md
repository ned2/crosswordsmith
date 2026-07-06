# Plan — Close the three v1 serializer coverage gaps (`xword`)

Status: proposed (plan only — no code/spec changes made by this document).
Scope: `xword/` (the standalone Python tool). Owner-facing spec references are to
`docs/xword-spec.md` and the tracker `docs/xword-status.md`.

## 1. The item

`docs/xword-spec.md` §10 (`docs/xword-spec.md:364-369`) records three v1 serializer
coverage gaps that currently **FAIL-STRICT** under policy **D7** (they hard-error
naming a capable target instead of serializing), flagged as "ordinary later work,
not a deferred decision":

- **rebus → Exolve** — Exolve can't currently be emitted for a multi-char cell.
- **prefilled cells → ipuz** — ipuz can't currently be emitted for a given/prefilled cell.
- **shaded / arbitrary-styled cells → Exolve** — Exolve can't currently be emitted for an ipuz cell style.

Goal: make each of these actually serialize into its **native representation** in
the target format, remove/narrow the corresponding fail-strict guard, and add
coverage. Native stays out of scope — native is a genuine `✗` in the §10 table
(`docs/xword-spec.md:324-328`), i.e. a real structural limit, not a serializer gap.

## 2. The load-bearing discovery — the parsers don't read these features either

The naive framing is "the serializer is the inverse of how the parser reads the
feature." That only half-holds: for a genuine **round-trip** (serialize → *parse
back* → assert fidelity), the *target-format parser* must read what we newly emit,
and today it does not:

- **Exolve parser** tokenises one char per main + decorator suffixes
  (`_tokenize_row`, `xword/src/xword/formats/exolve.py:106-127`); it **cannot read
  rebus**, **drops every `exolve-option`** (`exolve.py:43-59`, see the comment at
  line 58), and has **no `exolve-colour`/`exolve-color` handling** at all.
- **ipuz parser** reads a `value` only from the **solution** array
  (`xword/src/xword/formats/ipuz.py:98`) and **never sets `Cell.prefilled`**
  (`ipuz.py:79-86`); it **cannot read a prefilled `value` off a puzzle cell**.

Consequence: "close the gap" is a **scoping fork** per gap (§4). The serializer
change alone closes the fail-strict gap that the tracker measures; a symmetric
reader is an *additional* piece required only for a true parse-back round-trip.

## 3. Where everything lives (evidence)

### Serializers
| Serializer | Function | File:line |
|---|---|---|
| Exolve | `serialize` | `xword/src/xword/formats/exolve.py:162-207` |
| ipuz | `serialize` | `xword/src/xword/formats/ipuz.py:170-220` |
| native | `serialize` | `xword/src/xword/formats/native.py:100-191` |

Dispatch registry: `xword/src/xword/formats/__init__.py:10-11`
(`parse_board`/`serialize_board`).

### The three fail-strict guards (this task removes/narrows)
| Gap | Guard | File:line |
|---|---|---|
| rebus → Exolve | `if cell.letter is not None and len(cell.letter) > 1:` raise | `exolve.py:177-181` |
| prefilled → ipuz | `if cell.prefilled:` raise | `ipuz.py:182-186` |
| shaded/styled → Exolve | `if cell.style and not _style_expressible(cell.style):` raise (helper `_style_expressible`) | `exolve.py:182-186` (helper `exolve.py:148-156`) |

Native's parallel guards (`native.py:129-133` prefilled, `134-138` style,
`144-148` rebus) are **out of scope** — native is `✗` in §10, not a serializer gap.

### Model fields carrying each feature
`Cell` dataclass — `xword/src/xword/board.py:25-37`:
- rebus text → `Cell.letter` (multi-char) — `board.py:29`
- prefill → `Cell.prefilled: bool` — `board.py:36`
- shading → `Cell.style: Meta | None` ("ipuz style passthrough (verbatim)") — `board.py:37`
- circle/bars have dedicated flags `Cell.circle` / `Cell.bar_right` / `Cell.bar_below` — `board.py:33-35`

### Fixtures grounding the encodings
- `xword/tests/fixtures/sample.ipuz.json` — has the `"PH"` rebus (line 13) and a
  `{"cell": 2, "style": {"shapebg": "circle"}}` style object (line 12). No prefill, no shade.
- `xword/tests/fixtures/sample_decorated.exolve` — exercises `@` circle, `|`/`_`
  bars (lines 8-11); does **not** exercise rebus, `!` prefill, or `exolve-colour`.
- New fixtures are therefore required for rebus / prefill / shade (§4).

## 4. Authoritative target-format encodings (verified)

### 4.1 Exolve rebus (Gap 1)
Verified against the Exolve README. Requires the directive `exolve-option:
rebus-cells`, after which **every entry in every grid row is space-separated**:

```
exolve-width: 3
exolve-height: 3
exolve-option: rebus-cells
exolve-grid:
  RAN G E
   DO . A
    M E T
```

Row 1 = cells `RAN`, `G`, `E`; row 2 = `DO`, block, `A`. **Not** square brackets.
Decorators follow the cell content (`RAN@` = rebus + circle). This is a whole-grid
mode switch — incompatible with the current whitespace-*skipping* tokeniser
(`exolve.py:112` `if main.isspace(): continue`), which must gain a rebus-mode
whitespace-*splitting* branch.

### 4.2 ipuz prefilled (Gap 2)
Verified against libipuz / go-ipuz / puzzazz references. The given letter lives in
the puzzle cell's `value` field:

```json
{ "cell": 1, "value": "A", "style": { ... } }
```

It appears in **both** `puzzle` (as `value`) and `solution`. This mirrors the
parser's existing solution-side `sc.get("value")` read (`ipuz.py:98`), just applied
to the puzzle cell object.

### 4.3 Exolve shaded cells (Gap 3)
Verified against the Exolve README. A top-level directive, **not** a grid decorator:

```
exolve-colour: <html-colour-no-spaces> <cell-spec> [<cell-spec>...]
exolve-color:  <html-colour-no-spaces> <cell-spec> [<cell-spec>...]
```

Cell-specs are chessboard notation (`rNcM`, `a3`, `c10r9`) or clue indices
(`A12`, `12a`). The ipuz source is a StyleSpec background-shade key (`color`,
`highlight`) carried verbatim in `Cell.style` (`ipuz.py:84`). Exolve's colour model
(one colour + cell list) is **strictly coarser** than the ipuz StyleSpec, so the
reverse direction is inherently lossy — see the scoping recommendation in §5.3.

## 5. The scoping fork (decide per gap)

Two ways to "close" each gap:

- **(a) Serialize-only** — emit the native target representation and add a
  *convert-level structural assertion* (the output contains the feature and is
  valid in the target). This closes the §10 fail-strict gap the tracker actually
  measures (`docs/xword-status.md:84`) at ~half the effort. No reverse reader.
- **(b) Serialize + parse-back round-trip** — also add the target-format *reader*
  so `parse(serialize(x))` is Board-equal — full symmetric fidelity.

**Recommendation, per gap:**

| Gap | Recommendation | Rationale |
|---|---|---|
| 2 prefill → ipuz | **(b)** full round-trip | Reader is ~5 lines; `value` is a clean symmetric encoding. |
| 1 rebus → Exolve | **(b)** full round-trip | Round-trip matters; the rebus-mode tokeniser branch is well-defined. |
| 3 shaded → Exolve | **(a) serialize-only + documented lossiness** | Exolve colour is strictly coarser than ipuz StyleSpec; a full reverse coordinate parser buys an *inherently lossy* round-trip. Not worth the complexity in v1. |

Per-gap detail and both-option effort are in §6.

## 6. Per-gap plan

Each gap is independent (distinct guard, fixture, test class) and can land as its
own commit. Sequence and risk are in §7.

### Gap 2 — prefilled → ipuz  (recommend option (b))

- **Serializer** (`ipuz.serialize`, cell loop `ipuz.py:177-194`): remove guard
  `ipuz.py:182-186`. For a prefilled cell, emit the puzzle cell as an object
  carrying `value` (the letter) alongside `cell`/`style` — merge with the existing
  `{"cell", "style"}` object path at `ipuz.py:187-189`. Keep the letter in
  `solution` (already done at `ipuz.py:191-194`).
- **Parser reader (option (b))** (`ipuz.parse`, puzzle-cell loop `ipuz.py:63-87`):
  when a puzzle cell dict carries `value`, set `Cell.prefilled=True` and
  `Cell.letter=value`; reconcile with the solution letter (`ipuz.py:90-103`).
- **Guard**: remove `ipuz.py:182-186`.
- **Tests**:
  - Invert `TestStructuralFailures::test_prefilled_cell_blocks_ipuz`
    (`xword/tests/test_convert.py:97-100`) — it must now *succeed*.
  - Keep `test_prefilled_cell_blocks_native` (`test_convert.py:87-90`) unchanged.
  - Option (a): convert `mini_exolve("C!AB", "A.U", "BUS")` (`test_convert.py:38`)
    → ipuz; assert the emitted ipuz has `value` on the given cell + is valid.
  - Option (b): also parse the ipuz back and assert Board-equal (prefill + letter
    survive), patterned on `TestIpuz::test_sample_board_round_trip`
    (`test_formats.py:64-68`); plus convert-level exolve(`!`)→ipuz→exolve.
- **Fixture**: none strictly required (reuse `mini_exolve`); optionally add a tiny
  `xword/tests/fixtures/sample_prefilled.exolve`.
- **Spec/status**: drop "prefilled cells → ipuz" from the gap bullet
  `docs/xword-spec.md:364-369`; update the Phase-2 item-2 row `docs/xword-status.md:84`
  (and, for option (b), the ipuz-parser row `docs/xword-status.md:60`).
- **Effort**: (a) ~0.25 day · (b) ~0.5 day.

### Gap 3 — shaded/styled → Exolve  (recommend option (a) serialize-only)

- **Serializer** (`exolve.serialize`, `exolve.py:162-207`): narrow the guard
  `exolve.py:182-186` and widen `_style_expressible` (`exolve.py:148-156`). Collect
  the background-shade style keys — **canonical key `color`**, plus `highlight`
  (confirm the exact ipuz StyleSpec set at implementation) — and emit deterministic
  `exolve-colour: <value> <cellspec>` lines. Keys Exolve genuinely can't hold
  (`colortext`, `colorborder`, `imagebg`, `mark`, …) **still fail-strict** — this is
  a *narrow*, not a removal.
- **Pre-work (do before starting this gap)**: pin the canonical key to `color` and
  update the non-standard `{"fill": "DDDDDD"}` case in
  `TestStructuralFailures::test_foreign_styling_blocks_native_and_exolve`
  (`test_convert.py:106-112`). After Gap 3, `color` must *serialize* to Exolve, so
  that test's exolve target must move to a genuinely-unmappable style key; its
  native half stays (native is still `✗`).
- **Parser reader (option (b), NOT recommended)**: would add an
  `exolve-colour`/`exolve-color` directive handler to `exolve.parse`
  (`exolve.py:32-96`) — a chessboard-notation coordinate parser (`rNcM`, `aN`,
  clue-index) back to `[row,col]`, re-populating `Cell.style` (mirror of
  `ipuz._synth_style`, `ipuz.py:159-167`) with the canonical `{"color": <value>}`.
  Because Exolve colour is coarser than the StyleSpec, even this round-trip is
  puzzle-lossless, not payload-lossless.
- **Guard**: narrow `exolve.py:182-186` + widen `_style_expressible`
  (`exolve.py:148-156`).
- **Tests**:
  - Update `test_foreign_styling_blocks_native_and_exolve` (`test_convert.py:106-112`)
    per the pre-work above.
  - Option (a): convert an ipuz with `{"style": {"color": "pink"}}` → exolve; assert
    a deterministic `exolve-colour: pink rNcM` line appears + non-styled cells
    unchanged.
  - Option (b) only: also parse back and assert Board-equal (canonical `color`).
- **Fixture**: add `xword/tests/fixtures/sample_shaded.ipuz.json` (a cell with a
  background-colour StyleSpec).
- **Spec/status**: drop "shaded/arbitrary-styled cells → Exolve" from
  `docs/xword-spec.md:364-369`; update `docs/xword-status.md:84`.
- **Residual lossiness (option (a))**: `exolve→ipuz` cannot reconstruct the exact
  StyleSpec — document this in §10 as puzzle-lossless-only, consistent with
  `docs/xword-spec.md:360-363`. It is **structure**, not metadata, so it does not
  warn (warnings are metadata-only, `xword/src/xword/convert.py:19-38`); it either
  serializes the mapped shade key or fails-strict on an unmapped one.
- **Effort**: (a) ~0.75-1 day · (b) ~1.5-2 days.

### Gap 1 — rebus → Exolve  (recommend option (b); HIGHEST risk to D6 goldens)

- **Serializer** (`exolve.serialize`, grid loop `exolve.py:171-196`): remove guard
  `exolve.py:177-181`. Detect any multi-char `cell.letter`; when present, emit
  `exolve-option: rebus-cells` (positioned after the height line `exolve.py:165`,
  before `exolve-grid:` line 170) and switch the whole grid to **space-joined
  tokens** (multi-char letters verbatim, decorators appended). When no rebus cell
  exists, keep the current compact non-spaced output **byte-for-byte** (see §8).
- **Parser reader (option (b))** (`exolve.parse` `exolve.py:32-96`; `_tokenize_row`
  `exolve.py:106-127`): recognise `exolve-option: rebus-cells` (currently dropped at
  `exolve.py:58`) and, in rebus mode, tokenise rows by **splitting on whitespace**
  (each token = multi-char main + trailing decorator run) instead of the char-by-char
  scan. `letter = main if main.isalpha()` (`exolve.py:79`) already generalises to
  multi-char.
- **Guard**: remove `exolve.py:177-181`.
- **Tests**:
  - Invert `TestStructuralFailures::test_rebus_blocks_exolve`
    (`test_convert.py:102-104`) — must now *succeed*.
  - Keep `TestIpuz::test_rebus_blocks_native_serialization` (`test_formats.py:82-93`).
  - Option (a): convert `sample.ipuz.json` (has `"PH"`) → exolve; assert
    `exolve-option: rebus-cells` present and grid space-separated.
  - Option (b): also parse back → Board-equal; plus a `TestExolve`
    parse→serialize→parse over a new rebus fixture.
- **Fixture**: add `xword/tests/fixtures/sample_rebus.exolve` (small grid, one
  multi-char cell + a decorator, exercising the space-separated tokeniser).
  `sample.ipuz.json` is the convert source.
- **Spec/status**: drop "rebus → Exolve" from `docs/xword-spec.md:364-369`; update
  `docs/xword-status.md:84` (and the exolve-parser row `docs/xword-status.md:61`
  for option (b)).
- **Why highest risk**: this rewrites the grid tokeniser's contract (compact vs
  space-separated) and is the change most likely to move the D6 byte-goldens. It
  must be gated behind the hardest "non-rebus output stays byte-identical" test (§8).
- **Effort**: (a) ~0.75 day · (b) ~1.25-1.5 days.

## 7. Sequence & risk ordering

Gap 2 first (lowest risk, isolated, ipuz-only) is correct. Reordering the rest:
with Gap 3 taken as **serialize-only** (no reverse coordinate parser) it becomes a
contained, additive new-directive emission; Gap 1 as **serialize + parse-back**
rewrites the Exolve grid-tokeniser contract and is the biggest threat to the D6
byte-goldens. So tackle the tokeniser-contract change **last**, behind the strongest
gate:

1. **Gap 2 — prefill → ipuz** — option (b). Lowest risk; ipuz only; `dump_json`
   keeps output sorted/idempotent.
2. **Gap 3 — shaded → Exolve** — option (a). Additive `exolve-colour` lines; gate
   the no-style path byte-identical; do the canonical-`color` pre-work first (§6).
3. **Gap 1 — rebus → Exolve** — option (b). Rewrites the grid tokeniser contract;
   strongest byte-identity gate; do last.

Total recommended effort (2b + 3a + 1b): **~2.5-3 days**.

## 8. D6 byte-identity acceptance gate (per gap)

D6 (`docs/xword-spec.md:94-96`) requires deterministic output. The **no-feature
path must not move**: `bundled_17` carries none of rebus/prefill/shading, so its
convert output must be **byte-identical before and after each gap**. Make this an
explicit acceptance gate:

- **Before touching a serializer**, capture the current bytes of
  `convert_text(bundled_17.native.json, <target>)` for the affected target(s)
  (ipuz for Gap 2; exolve for Gaps 1 & 3) as a golden string.
- **After the change**, assert byte-equality against that golden. The gap's change
  must keep it identical.
- Existing structural cross-checks stay green:
  `TestExolve::test_board_round_trip` (`test_formats.py:97-101`),
  `TestIpuz::test_round_trip_identity_vs_engine` (`test_formats.py:47-51`),
  `TestEngineCrossCheck::test_exolve_matches_engine_export` /
  `test_ipuz_matches_engine_export` (`test_convert.py:197-203`).
- Determinism of the *new* output:
  - ipuz `value` — `dump_json` sorts keys (`ipuz.py:220`); stays sorted & idempotent.
  - Exolve `exolve-option` line — fixed position; emitted **only** when a multi-char
    cell exists.
  - Exolve `exolve-colour` lines — emitted **only** when styled cells exist, in a
    deterministic order (sort by `(colour, row, col)`), one canonical cell-spec form
    (`rNcM`, 1-indexed).

**Render/view byte-goldens are unaffected** — this is a convert-only change; the
SVG renderer already draws circle/bar/**rebus** (`docs/xword-status.md:104`), so no
presentation work is needed, and the render/view goldens (all over `bundled_17`,
`xword/tests/test_render.py`, `xword/tests/test_view.py`) do not move.

## 9. Spec & tracker updates (summary)

- `docs/xword-spec.md`:
  - §10 gap bullet **`docs/xword-spec.md:364-369`** — remove each gap as it closes;
    when all three land, reword the bullet to note the gaps are closed (and, for
    Gap 3 option (a), that `exolve→ipuz` style is puzzle-lossless-only per the §10
    round-trip note at `docs/xword-spec.md:360-363`).
  - §10 **table** (`docs/xword-spec.md:322-330`) — **no change**; ipuz/Exolve are
    already `✓` for these properties (the *format* could always hold them).
  - **D7** wording (`docs/xword-spec.md:97-108`) — **no change** (the structure vs
    metadata split is unchanged).
- `docs/xword-status.md`:
  - Phase-2 item-2 row **`docs/xword-status.md:84`** — drop the closed serializer
    gaps from "v1 serializer gaps (rebus→exolve, prefill→ipuz, styling→exolve) fail
    likewise".
  - For option (b) gaps, update the parser rows **`docs/xword-status.md:60`** (ipuz)
    / **`docs/xword-status.md:61`** (exolve) to note the new readers, and the
    round-trip row **`docs/xword-status.md:87`** for any new round-trips.

## 10. Test-inversion checklist (existing tests that must change)

| Test | File:line | Change |
|---|---|---|
| `test_prefilled_cell_blocks_ipuz` | `test_convert.py:97-100` | Invert → now succeeds (Gap 2) |
| `test_rebus_blocks_exolve` | `test_convert.py:102-104` | Invert → now succeeds (Gap 1) |
| `test_foreign_styling_blocks_native_and_exolve` | `test_convert.py:106-112` | Narrow: keep native `✗`; move exolve target to an unmappable style key (Gap 3) |
| `test_prefilled_cell_blocks_native` | `test_convert.py:87-90` | **Unchanged** (native still `✗`) |
| `test_rebus_blocks_native_serialization` | `test_formats.py:82-93` | **Unchanged** (native still `✗`) |

## 11. Out of scope

- Native serializer/parser (genuine `✗` in §10, not a serializer gap).
- Presentation/render of these features (SVG already renders rebus/circle/bar;
  terminal/HTML/PNG/PDF unchanged).
- `.puz` and other deferred formats (D5, `docs/xword-spec.md:92-93`).
