# Spec: JSON solution output with per-word link metadata

Status: **implemented** — `crossword.pl`, `fixtures/bundled_17_clues.pl`, the
tests, and the README now reflect this design.

Companion: [`json-input-spec.md`](json-input-spec.md) defines the JSON
*input* format, designed as the natural mirror of this output so a solution
payload reduces back to a valid input.

## 1. Problem statement

`crossword.pl` emits a solved crossword as an ad-hoc, delimiter-based text
format (see `crossword/3`, `crossword.pl:105-111`, and the printing
predicates `print_grid/2` / `print_clues/1`, `crossword.pl:444-486`). The
output is meant to be consumed by a UI that renders the grid and turns each
word into a hyperlink. The format has three sections:

1. A grid: `grid_length` rows of cells. An empty cell is `*`; a filled cell
   is five comma-joined fields `Across,Down,Letter,Marker,Link`.
2. A line containing a single `@`.
3. Clue lists under `Across Clues` / `Down Clues`, one per line as
   `Word|ClueNum|Clue|Link`.

### Issues with the current format

1. **Ad-hoc string format.** Deserialising requires bespoke parsing logic:
   split on newlines, then on spaces, then on commas for the grid and on
   pipes for the clues. It is also latently fragile — the grid uses `,` as a
   delimiter and the clues use `|`, so any clue text, answer, or URL
   containing those characters would mis-parse. Nothing currently escapes
   them.

2. **Links are tracked per cell, and crossing cells lose optionality.** Each
   cell of a word repeats that word's URL. For a cell where an across and a
   down word cross, only one link can be stored, so the solver must pick one
   *upstream* — in the solving/annotation layer (`letters_to_grid/8`,
   `crossword.pl:307-341`). This is the wrong layer to make that call: the UI
   is what knows the user's current direction (across vs down) and should
   decide which link to follow. The current pick is also lossy in a subtle
   way. Tracing `letters_to_grid/8` (`crossword.pl:331-338`): the first word
   to reach a cell writes its own link, then the across pass *keeps* the
   existing link while the down pass *overwrites* it — so at any across/down
   crossing the **down** word's link wins, deterministically (confirmed at
   cell 4 of `tests/golden/grid_17_topleft_across.txt`, where the crossing of
   OMEGA POINT and GNOSTIC GOSPELS stores the Gnostic Gospels link). The
   across word's link survives only on its *non-crossing* cells. The practical
   consequence: the link a UI gets by clicking a cell depends on *which* cell
   of the word was clicked.

## 2. Goals

- **G1.** Emit the solution as **JSON**, eliminating custom parsing and the
  delimiter-collision fragility.
- **G2.** Preserve the association between every cell/letter and the word(s)
  it belongs to, so that for a crossing cell **both** the across and the down
  word (and therefore both links) remain available. The UI chooses which to
  follow.
- **G3.** Stop attaching links to cells. A link is a property of a **word**,
  stored once per word.
- **G4.** Make the solver **metadata-agnostic**: it needs only the answer
  string to solve. Clue text, link, and any other per-word metadata are
  opaque and passed straight through into the output payload, so the UI gets
  one self-contained artifact (no separate merge step).

## 3. Non-goals

- Changing the solver algorithm, placement rules, or clue-numbering logic.
- Inferring grid size, deduplicating solutions, or any search behaviour.
- Defining the UI's rendering or interaction model beyond what the payload
  must support.

## 4. Decisions taken

These were settled during design discussion:

- **Metadata strategy: emit-time join, single payload.** The solver is given
  only the answer string and never carries metadata. At output time the
  emitter rejoins each placed word to its metadata (keyed on the answer
  string, looked up in the in-scope input list) and writes one self-contained
  JSON payload. What we reject is pushing that merge *downstream* onto the UI
  as a separate artifact — the join happens in-process, so the UI still
  receives a single document. This is the purest form of G4: the solver does
  not merely ignore metadata, it never sees it. The one assumption it buys is
  that answers are unique (see §8.1 and §9).
- **Input format: answer + metadata dict.** The bundled Prolog fixture moves
  from positional triples `[Word, Clue, Link]` to `[Answer, Metadata]`, where
  `Metadata` is a dict of arbitrary passthrough keys (conventionally `clue`
  and `link`). The solver consumes only `Answer`; it does not special-case any
  metadata key.

## 5. Proposed input format

`clues/1` becomes a list of `[Answer, Metadata]` pairs:

```prolog
clues([
    ['OMEGA POINT',
     _{clue: 'Transcending entropy',
       link: 'http://en.wikipedia.org/wiki/Omega_Point'}],
    ['GNOSTIC GOSPELS',
     _{clue: 'Some apocrypha',
       link: 'http://en.wikipedia.org/wiki/Gnostic_Gospels'}],
    % metadata is optional; an answer may be given with an empty dict ...
    ['FLOW', _{}],
    % ... and, for convenience, the metadata element may be omitted entirely:
    ['BIAS']
]).
```

Rules:

- **Answer** — required. Spaces are allowed for multi-word answers and are
  stripped before letters are placed on the grid (unchanged behaviour). The
  answer retains its spaces as a display form in the output.
- **Metadata** — optional dict. It is treated as opaque and copied verbatim
  into the corresponding word's output entry, under a `meta` key (§6.2).
  `clue` and `link` are conventions, not requirements; a word with no metadata
  has an empty `meta` object.

This satisfies G4: the solver only consumes `Answer`; the metadata is held
aside and reattached at output time (§8), never threaded through the solver.

## 6. Proposed output format

A single JSON object:

```json
{
  "gridLength": 17,
  "grid": [
    [ {"letter": "O", "number": 1,    "across": 1,    "down": null},
      {"letter": "M", "number": null, "across": 1,    "down": null},
      {"letter": "E", "number": null, "across": 1,    "down": null},
      {"letter": "G", "number": 2,    "across": 1,    "down": 2},
      null, null, "... (gridLength cells per row) ..." ],
    "... (gridLength rows) ..."
  ],
  "words": [
    {
      "number": 1,
      "direction": "across",
      "answer": "OMEGA POINT",
      "cells": [[0,0],[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7],[0,8],[0,9]],
      "meta": {
        "clue": "Transcending entropy",
        "link": "http://en.wikipedia.org/wiki/Omega_Point"
      }
    },
    {
      "number": 2,
      "direction": "down",
      "answer": "GNOSTIC GOSPELS",
      "cells": [[0,3],[1,3],"..."],
      "meta": {
        "clue": "Some apocrypha",
        "link": "http://en.wikipedia.org/wiki/Gnostic_Gospels"
      }
    }
  ]
}
```

Field order in the rendered output is not significant and will not match the
example above literally: `json_write_dict/2` emits object keys in sorted
(alphabetical) order. That order is **deterministic** for a given input, so a
byte-exact golden test remains viable; the field order shown in this document
is chosen for readability only. (If logical key order in the output ever
matters, build the payload as `json([...])` terms and use `json_write/2`,
which preserves insertion order.)

### 6.1 `grid`

A `gridLength` × `gridLength` array, row-major (row 0 is the top row).
Coordinates are 0-indexed `[row, col]`. Each entry is either:

- `null` — an empty (black) cell; or
- a cell object:

  | Field     | Type            | Meaning |
  |-----------|-----------------|---------|
  | `letter`  | string          | the cell's letter |
  | `number`  | integer \| null | the corner label shown in the UI; non-null **iff** a word starts in this cell. Across and down words starting in the same cell share this number. |
  | `across`  | integer \| null | clue number of the across word passing through this cell, or `null` if none |
  | `down`    | integer \| null | clue number of the down word passing through this cell, or `null` if none |

`across` / `down` replace the old `x` sentinel with `null`. To follow a link
from a cell, the UI reads `cell.across` (or `cell.down`) and looks up the word
whose `number` **and** `direction` match — number alone is ambiguous, since an
across and a down word starting in the same cell share a number. A crossing
cell carries both `across` and `down`, so both links are reachable (**G2**).

The old single-character `Marker` (`a`/`d`/`n`) is subsumed by `number`: a
cell shows a label iff `number` is non-null. If the UI needs to know whether
that start is across, down, or both, it is re-derivable from `words` (a word
starts at `cells[0]`).

### 6.2 `words`

One entry per placed word:

| Field       | Type    | Source | Meaning |
|-------------|---------|--------|---------|
| `number`    | integer | solver | clue number |
| `direction` | string  | solver | `"across"` or `"down"` |
| `answer`    | string  | input  | the answer, spaces preserved (display form) |
| `cells`     | array   | solver | ordered `[row,col]` pairs the word occupies |
| `meta`      | object  | input  | the input `Metadata` dict, copied verbatim (e.g. `clue`, `link`); `{}` if none |

The metadata — and therefore the link — lives here under `meta`, once per
word (**G3**). It is nested rather than spread across the word object so that
arbitrary passthrough keys can never collide with the structural fields
(`number`, `direction`, `answer`, `cells`). `cells` gives the UI word → cells
for highlighting; `cell.across` / `cell.down` give it cell → word for click
handling.

`answer` is a display label only. For a multi-word answer its spaces are not
placed on the grid, so `cells` has one entry per *letter* and is shorter than
`answer` (e.g. `"OMEGA POINT"` → 10 cells). A UI must therefore not map
`answer` characters to `cells` positionally; the authoritative per-cell letter
is `grid[row][col].letter`.

### 6.3 No solution

If no layout exists for the requested `grid_length` (the grid is too small, or
the start position admits no solution), `find_crossword/5` fails and
`crossword/3` produces **no output** and exits non-zero — the same outcome as
today, but now an explicit contract. The emitter runs only on success, so a
consumer can treat "empty stdout / non-zero exit" as "no solution" and any
stdout as a complete JSON document. (Emitting `{"error": "no solution"}` with
exit 0 was considered and rejected: a non-zero exit is the more Unix-friendly
signal and keeps every success payload unconditionally valid JSON.)

### 6.4 `diagnostics`

An **optional** top-level object carrying quality caveats about how the
payload was produced. Arrange output is best-effort by nature — most real
clue sets will not achieve perfect health scores — so compromises are data
for the consumer, not stderr noise (design-spec §5.1, INV-3). Verbs that
compromise nothing (e.g. `fill`) omit the property entirely; a consumer must
treat a missing `diagnostics` as "no caveats reported".

Each key names the producing verb, so the object can grow richer per-producer
metadata later (e.g. an embedded `lint` report) without collisions:

```json
"diagnostics": {
  "arrange": {
    "capInert": true,
    "dropped": ["NARRATIVE FALLACY"],
    "reward": 47
  }
}
```

The `arrange` sub-object (present on every arrange payload, including each
element of a `--candidates` array — caveats are per-layout):

| Field      | Type    | Meaning |
|------------|---------|---------|
| `capInert` | boolean | `true` when **no** placed word reaches its checking target — the capped objective has degenerated to plain total-crossings (design-spec §7.2); tune with `--check-target` |
| `dropped`  | array   | answers from the input that were not placed (best-effort drops, AC-ARR-2), in input order; `[]` when everything placed |
| `reward`   | integer | the engine's objective value for this layout (capped-checking reward, design-spec §7.2) |
| `seed`     | integer | **optional** — the RNG seed that produced this layout, present **only** when `--seed N` or `--shuffle` perturbed the search; reproduce the exact layout with `--seed <value>`. Absent (the default) means the deterministic search, so a missing `seed` = "not perturbed" |

Structural fields (`gridLength`, `grid`, `words`) remain the complete layout
contract: a consumer may ignore `diagnostics` wholesale and lose nothing but
provenance. Downstream converters that cannot represent it (ipuz, Exolve)
treat it as metadata-class content (drop-and-warn).

## 7. Mapping from the current format

| Current | New |
|---------|-----|
| grid cell `Across,Down` (`x` = none) | `cell.across` / `cell.down` (`null` = none) |
| grid cell `Letter` | `cell.letter` |
| grid cell `Marker` `a`/`d`/`n` | `cell.number` (non-null iff a start cell) |
| grid cell `Link` (repeated per cell) | **removed**; now `words[].meta.link` |
| `*` empty cell | `null` |
| `@` separator | removed (JSON is self-describing) |
| `Across Clues` / `Down Clues` headings | removed; `words[].direction` distinguishes |
| `Word\|ClueNum\|Clue\|Link` clue line | a `words[]` entry |

## 8. Implementation sketch

*Point-in-time record (pre-module-split): file/line references below name the
old layout — root `crossword.pl` is now `prolog/crosswordsmith/core.pl`, and
`tests/crossword.plt` is now `tests/core.plt` (legacy-surface dissolution,
2026-07-07). The sketch itself was executed; kept as the design record.*

The solver *algorithm* is unchanged — placement, intersection, and clue
numbering all stay as-is. What changes is the shape of the data the solver
carries (the placed-word tuple loses its metadata slots) and the I/O edges
(input format and output emitter). `find_crossword/5`'s logic is untouched;
`assign_words/8` changes only in how it destructures each input entry.

1. **Input normalisation.** Update the bundled Prolog fixture to the
   `[Answer, Metadata]` form, and have `assign_words/8` take only the answer
   (`[Answer | _]`, `crossword.pl:133`). Metadata is **not** carried by the solver, so the
   placed-word tuple drops its `Clue` and `Link` slots entirely (9 elements
   to 7). That arity change shifts every positional match on the tuple. Sites
   to update (or, for the old grid builder, remove):

   - `assign_word/11` — the tuple is built here (`crossword.pl:168`)
   - `find_intersecting_word/6` — `member([_,_,_,PLetters,_,PDir,_,PStart,_], …)`
     (`crossword.pl:152`); two leading underscores go away
   - `dir_is_across/1` (`crossword.pl:352`) and `start_is/2` (`crossword.pl:354`)
   - `add_clue_word/3` (`crossword.pl:357-359`)
   - `words_to_grid/3` + `letters_to_grid/8` (`crossword.pl:298-341`) —
     **removed** with the old grid builder (step 2), not updated
   - `x_print_clues/1` (`crossword.pl:482`) — superseded by the emitter
   - the literal tuples and helpers in `tests/crossword.plt`
     (`:110`, `:114`, `:145-146`, `:150-166`)

   Because every site is being touched anyway, this is the moment to replace
   the bare positional list with a named representation — an SWI dict
   (`word{answer:…, letters:…, …}`) or a tagged compound `pw/7` with accessor
   predicates — so future field changes stay local and an arity mismatch
   becomes a hard error rather than a silent underscore miscount. A dict
   composes neatly with the dict-based input and `json_write_dict` output; a
   tagged compound keeps the solver core ISO-portable (see §9).

   **Answer-uniqueness guard.** The emit-time join (step 2) keys on the answer
   string, so it is well-defined only if answers are distinct. Add a check
   (e.g. `sort` the answers and compare lengths, or `is_set/1`) that fails
   with a clear message if `clues/1` contains a duplicate answer. Crosswords
   do not repeat answers, so this is a guard, not a feature.

2. **Output emitter.** Replace `print_grid/2` + the `@` write +
   `print_clues/1` in `crossword/3` (`crossword.pl:108-111`) with a single
   emitter that builds a dict and writes it via `json_write_dict/2`
   (`library(http/json)`; canonically `library(json)` on SWI 10+).

   - Cell → coordinates: `Row is (Cell-1) // GridLen`,
     `Col is (Cell-1) mod GridLen`.
   - `grid`: build fresh from the numbered placed words. For each across word,
     every covered cell gets `across = number`; likewise `down`. A cell's
     `number` is the clue number iff it is the start cell of some word. This
     replaces `annotate_grid/3` + `letters_to_grid/8` (`crossword.pl:293-341`),
     which can be **deleted** — they bake the per-cell link in (the very bug
     G3 removes) and read the now-removed `Link` slot.
   - `words`: from the numbered placed words — `number`, `direction`,
     `answer`, `cells` (mapped to `[row,col]`). The emitter has the input list
     `Words` in scope (`crossword/3`'s second argument), so it builds an
     `Answer → Metadata` map from it and attaches each word's metadata under a
     `meta` key (`{}` when absent). The solver output and the metadata are
     joined here, in one place.

3. **Docs.** Update `README.md` ("Finishing and output", lines ~162-191, and
   the clues-file section, ~80-106) to describe the JSON schema and the new
   input format.

4. **Tests.** Both layers in `run_tests.sh` need updating:
   - The golden regression compares the raw CLI output of
     `crossword.pl --input fixtures/bundled_17_clues.pl 17 topleft_across` against
     `tests/golden/grid_17_topleft_across.txt`. Regenerate the golden file in
     the new JSON format. Consider parsing the JSON and asserting on structure
     rather than byte-exact match, since key ordering in `json_write_dict` is
     stable but the format is now structured.
   - The plunit suite (`tests/crossword.plt:132`) currently only checks that
     `crossword/3` runs without error via `with_output_to/2`. Add assertions
     that the emitted string parses as JSON and has the expected
     `gridLength` / `grid` / `words` shape, that a known crossing cell carries
     both `across` and `down`, and that each word's metadata appears nested
     under `meta`. The fixtures in the suite also move to the new input /
     placed-word shape (see the site list in §8.1).

## 9. Risks and considerations

- **Breaking change.** Both the Prolog clue fixture shape and output formats
  change. Any existing consumer of the text format must be updated. This is
  acceptable given the redesign is the point; note it in the README / commit.
- **Emit-time join assumes unique answers.** Metadata is reattached by
  matching the answer string, so two `clues/1` entries with the same answer
  cannot be told apart. The uniqueness guard (§8.1) turns this into a clear
  up-front error. Threading metadata through the solver would remove the
  assumption (each placement would carry its own copy), but at the cost of
  plumbing metadata through the entire placement pipeline; the join keeps the
  solver metadata-free, and the answer-uniqueness assumption is harmless for
  crosswords (which do not repeat answers).
- **Portability vs. ergonomics of the placed-word representation.** Switching
  the placed-word to an SWI dict is the most ergonomic option and composes
  with the dict input and `json_write_dict` output, but dicts are
  SWI-specific. A tagged compound (`pw/7`) with accessor predicates is
  ISO-portable and still gives named access. The JSON I/O path is SWI-only
  regardless (dict input syntax + `json_write_dict`), so dict-ifying the core
  adds little new dependency — but the README's "ports with minimal effort"
  claim should be narrowed to note the JSON path is SWI-specific.
- **Delimiter fragility resolved.** JSON encoding removes the `,` / `|`
  collision risk for arbitrary clue text and URLs (a motivating issue).
- **Metadata typing.** Values are passed through verbatim. If a metadata value
  is a Prolog term that does not map cleanly to JSON, `json_write_dict` will
  error; we assume metadata values are JSON-friendly (strings, numbers,
  booleans, nested dicts/lists). Worth stating as an input contract.
- **Redundancy is intentional.** `cells` appears on words while `across` /
  `down` appear on cells. This two-way linking is deliberate (highlight vs.
  click-to-link) and cheap.

## 10. Open questions

- Should the empty-cell representation be `null` in a dense 2D array (as
  specified) or a sparse list of only-filled cells? Dense is simplest for a
  grid renderer; revisit if payload size matters.
- Should the golden test assert byte-exact JSON or parse-and-assert? Output is
  deterministic (§6), so byte-exact works; still leaning parse-and-assert for
  robustness against incidental formatting changes.
- Do we want a top-level `version` field on the payload to ease future format
  evolution?
