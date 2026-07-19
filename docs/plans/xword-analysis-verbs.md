# xword analysis verbs implementation plan

**Status:** COMPLETE (2026-07-19). Implements xword-spec D10-D12 and
§§6.4-6.5. The strategic rationale remains in
[`xword-breadth-expansion.md`](xword-breadth-expansion.md).

## Scope

Add two deterministic, read-only commands over the existing normalised `Board`:

- `stats`: fixed descriptive metrics in human or canonical JSON form.
- `diff`: structural comparison of two independently parsed layouts, with human
  or canonical JSON reports and diff-style result exits.

Explicitly out of scope: quality/lint verdicts, dictionary scoring, format
validation, barred-entry derivation, `.puz`/`.xd`, rendering changes, Textual,
and engine changes. `stats` rejects bars until the Board has truthful bar-aware
entry derivation; `diff` still compares bar flags.

## Architecture

- `xword/stats.py`: pure metric calculation plus human/JSON formatting.
- `xword/diff.py`: explicit semantic projection, stable granular comparison,
  plus human/JSON formatting.
- `xword/cli.py`: one-input `stats`; two-input `diff` with independent format
  overrides and one optional stdin (`-`) operand.
- Existing `parse_board`, `Board`, and `dump_json` remain the parse/model/JSON
  seams. The ipuz parser canonicalises left/top bars onto Board's right/below
  edge model, and native malformed-schema guards keep parse failures in the
  documented error channel. No serializer, dependency, or engine change is
  needed.

The comparison projection deliberately does not use dataclass equality. It
normalises missing clues/enumerations, removes duplicate circle/bar StyleSpec
keys, and excludes payload-only fields that interchange formats cannot share.

## Delivery

1. Lock D10-D12, §§6.4-6.5, this plan, and the status row.
2. Implement pure `stats` and `diff` modules.
3. Wire CLI streams, format overrides, output, and result exits.
4. Add focused unit and subprocess coverage, including cross-format fixtures.
5. Run `make test-xword`, `make xword-parity`, and the native suite; mark the
   tracker and this plan done in the same change.

Completed in one delivery. Verification: `make test-xword` = 155 passed, 3
optional raster skips; `make test` = 542 plunit checks plus all goldens and CLI
contracts passed; `make xword-parity` retained the documented baseline (Exolve
byte-identical, ipuz serializer-formatting differences only).

## Acceptance

- Human and JSON outputs are byte-deterministic and newline-terminated.
- Stats categories are internally consistent and use cell, not character,
  lengths for rebuses.
- Semantic equality survives representation-only format changes and metadata
  loss. Real structure loss remains visible: Exolve colour parsing is currently
  one-way, and the sanctioned title-less→Exolve `Untitled` addition is a title
  difference rather than a false equality.
- Every compared structural field yields a stable, localised difference.
- `diff` exits 0 equal, 1 different after emitting its report, and 2 on an
  operational failure; read/parse failures do not create `--out`.
- Existing xword output and engine parity behavior remain unchanged.
