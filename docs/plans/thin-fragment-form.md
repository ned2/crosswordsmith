# Plan: thin fragment convenience form (AC-FRAG-4)

Status: **done** (this branch; 2026-07-06). Spec: design-spec §6.6; acceptance
AC-FRAG-4 ("Thin convenience form and canonical form for the same fragment
produce identical results"), deferred by DP-1/OD-9, now implemented.

## Problem

`arrange --fragment` (and `fill --seeds`) accept only the **canonical**
fragment JSON — the emit format made partial: a `{gridLength, words:
[{answer, direction, cells:[[r,c]...]}...]}` object. Hand-authoring the
full `cells` run for every pin is tedious and error-prone; §6.6 always
envisaged a thin, hand-authorable convenience form that **desugars into**
the canonical one.

## Shape of the change

**Thin form** = a top-level JSON **list** of `{answer, row, col, dir}`
entries (`row`/`col` 0-based start cell, `dir` = `"across"|"down"`). The
top-level type is the form discriminator (object = canonical, list = thin);
anything else is a shaped error. Extra keys in an entry are ignored, like
canonical's `number`/`meta`.

The thin form carries **no `gridLength`**, so the consumer's size frame
provides N at desugar time:

- `arrange --fragment`: an explicit `--size N` / `--max-size N` frames the
  desugar; neither given → the default 15 (the same default as a bare
  invocation, `resolve_size/2`). The desugared fragment then reports
  `FragGridLen = N`, so the existing `reconcile_fragment_size/3` agrees
  trivially and the pin-and-search path runs unchanged.
- `fill --seeds`: **canonical-only (v1)**, rejected with a precise shaped
  error. Fill matches seeds to slots by raw cell numbers on its own grid, so
  a thin seed desugared at any other width would silently miss every slot
  (`fill_seed_no_slot`) — worse than a clear rejection. Threading the fill
  grid's size through `apply_seeds/3` would also change an internal arity
  that `benchmarks/fill_subjects.pl` (and `tests/fill.plt`) reach white-box;
  the benchmark tree is owned by a concurrent work item, so the fill
  extension is deferred, documented, and cheap to add later (pass the grid
  Size as the new `SizeCtx`).

**Desugar site**: `arrange.pl`'s fragment parse boundary. Both forms
converge on the same `frag(Answer, Dir, Start, CellNums)` terms before
`seed_from_fragment/6`; the legality/pin core and everything downstream is
form-agnostic and untouched.

## Mechanics

`prolog/crosswordsmith/arrange.pl`:

- `load_fragment/4` (new export): `load_fragment(+File, +SizeCtx,
  -FragGridLen, -Frags) is det`. `SizeCtx` ∈ `none` (no size flag → thin
  defaults to 15) | positive integer (the `--size`/`--max-size` value) |
  `canonical_only` (thin rejected — the fill boundary).
- `load_fragment/3` (kept export, unchanged contract for canonical files):
  now delegates to `/4` with `SizeCtx = canonical_only`. Existing importer
  (`fill.pl`) needs no change.
- `fragment_json/4` (internal): dispatch dict → existing
  `fragment_dict_words/3`; list → `thin_size/2` + `maplist
  fragment_thin_word/3`; other → throw `fragment_bad_toplevel`.
- `fragment_thin_word/3` (internal): validate `answer` (reuse
  `fragment_invalid_answer`), `dir` (`fragment_thin_invalid_dir`),
  `row`/`col` integers within [0, N) (`fragment_thin_invalid_position`);
  `Start is Row*N + Col + 1`; length = the answer's separator-stripped
  footprint (`word_letters/3`, same as the pin path); `fits_on_grid/4`
  else `fragment_thin_off_grid`; `word_cells/5` **generates** the run —
  a thin entry cannot be internally inconsistent, only off-grid.
- New `prolog:error_message//1` clauses for the four new shaped errors.

`crosswordsmith` (CLI script): the `--fragment` clause of `arrange_action`
computes `SizeFrame` (`none` | N) *before* `load_fragment` and passes it as
`SizeCtx`; `reconcile_fragment_size/3` keeps running after (a canonical
fragment still errors on a disagreeing `--size`). `--fragment` help text
mentions both forms.

## Equivalence argument (AC-FRAG-4)

For the same word pinned at the same place, the canonical entry's declared
`cells` (validated later by `expected_run/5` against the input answer's
footprint) are exactly `word_cells(Start, Dir, WLen, N)` — the run the thin
desugar generates. Identical `frag/4` terms in, identical deterministic
search out (INV-2). The golden check makes this byte-level: the thin fixture
is compared against the **same** golden file as the canonical fixture.

## Tests

- `tests/arrange.plt` (Phase 5 block): thin desugars to the identical
  `frag/4` terms as the canonical dict (the bundled two-word fragment);
  default-15 frame when `SizeCtx = none`; shaped errors for a `direction:`
  (not `dir:`) entry, out-of-range `row`/`col`, an answer that falls off
  the grid, a non-object non-list top level, and thin under
  `canonical_only`.
- Golden: new fixture `fixtures/bundled_17_fragment_thin.json` (the thin
  spelling of `fixtures/bundled_17_fragment.json`), checked with
  `--size 17` against the **existing** `tests/golden/
  arrange_bundled_17_fragment.json` in both `run_tests.sh` and the
  Makefile `golden` target. No golden file changes.

## Docs (same commit)

- design-spec §6.6: the DP-1 parenthetical → thin form implemented for
  `arrange --fragment` (AC-FRAG-4); `fill --seeds` canonical-only. §9 OD-9
  row annotated likewise.
- docs/STATUS.md: arrange Phase-5 row + "Fragment-grid primitive (§6.6)"
  row: AC-FRAG-4 no longer deferred; test counts refreshed.
- README: `--fragment` row + example gains the thin form; `--seeds` row
  notes canonical-only.
- docs/json-input-spec.md: **unchanged** — it specifies the `--input` clue
  file only; the fragment schema lives in §6.6/README.

## Verification

`make unit`, `make golden`, `make test`; determinism tags probed with
`deterministic/1` on the new `load_fragment/4` (thin + canonical paths).

## Non-goals

Letter-level pins (ninas) stay future §6.6 scope; `fill --seeds` thin
support deferred as above; no change to the canonical schema or emit.
