crosswordsmith
===============================================================================

A crossword layout generator in Prolog. 

This is a crossword grid generator implemented in SWI Prolog. The solver
core is portable Prolog, though the clue-input and JSON-output paths use
SWI-specific features (see Requirements). The solver works
deterministically, laying down one word at a time by placing it such
that it intersects with an existing word, and also such that it
doesn't sit hard up against any letters from another word.

The solver is designed to run online, so it returns the first solution
rather than trying to find them all.  The list of input words is
shuffled first so repeated solutions on different runs are less
likely. The start location for laying down the first word can also be
shuffled or specified in advance. The algorithm does not try guess the
grid length, so this must be specified by the user. It is of course
possible to specify a length for which there are no solutions. Also
note that duplicate solutions exist. These consist of the same layout
of words, but found by laying out words in a different order.

The words and their accompanying solutions are stored in a separate
file, clues.pl, which also allow you to attach a URL to each word, as
this was originally intented for embedding in a web page, where the
clues/solved words can be links.


## Requirements

[SWI-Prolog](https://www.swi-prolog.org/). The solver core is plain
Prolog, but the input and output paths are SWI-specific: clue metadata and
the JSON output both use SWI dicts and `library(http/json)`, so porting to
another Prolog would mean replacing those. It is developed and tested
against SWI and has been confirmed to run on SWI-Prolog 9.2.x (the version
in the Ubuntu repositories) and 10.0.x (from the `ppa:swi-prolog/stable`
PPA).


## Usage

`crossword.pl` is an executable SWI-Prolog script (note the `#!`
shebang at the top). Run it directly:

    # Place the words using a fixed starting position.
    $ ./crossword.pl <grid_length> <start_loc>

    # Find a vaguely random layout by shuffling the input words and the
    # order in which start positions are tried.
    $ ./crossword.pl --shuffle <grid_length>

    # Count how many distinct solutions exist (see the note below).
    $ ./crossword.pl --all <grid_length> [<start_loc>]

    # Use an external JSON clue file instead of the bundled clues.pl, and/or
    # write the output to a file. --clues and --out compose with any of the
    # forms above and may appear in any order.
    $ ./crossword.pl --clues puzzle.json <grid_length> <start_loc>
    $ ./crossword.pl --clues puzzle.json --shuffle <grid_length>
    $ ./crossword.pl --out solution.json <grid_length> <start_loc>

    # Show all options.
    $ ./crossword.pl --help

Arguments:

- `grid_length` — an integer giving the side length of the (square)
  grid. The words in `clues.pl` are intended for a grid of length 17.
- `start_loc` — where the *first* word is placed. One of:
  `topleft_across`, `topleft_down`, `topright`, `bottomleft`.

Options:

- `--clues <file>` — load the word/clue set from an external **JSON file**
  instead of the bundled `clues.pl`. This is the interchange format: it
  mirrors the JSON *output*, so a solution payload can be reduced back to a
  valid input. Without the flag, the bundled `clues.pl` is used (unchanged).
  See [Clues file](#clues-file) and
  [`docs/json-input-spec.md`](docs/json-input-spec.md).
- `--out <file>` — write the output to `<file>` instead of stdout. The file
  is written only once a solution is found, so an unsolvable run leaves no
  empty file behind. Without the flag, output goes to stdout as before.
- `--help` / `-h` — print the usage summary and exit.

Flags are parsed with `library(optparse)`, so they compose in any order and
also accept the `--flag=value` form (e.g. `--clues=puzzle.json`).

If you do not have execute permission set, you can equivalently run it
as `swipl crossword.pl <args...>`.

Examples:

    $ ./crossword.pl 17 topleft_across
    $ ./crossword.pl --shuffle 17

> **Note on `--all`:** this enumerates *every* solution by
> backtracking, and the search space is large — many solutions are
> duplicates that differ only in the order words were laid down. With a
> fixed `start_loc` it is bounded but can still be slow; with no
> `start_loc` it additionally enumerates all four starting positions and
> can take a very long time. Use it for small grids / small clue sets.


## Clues file

The words and clues live in `clues.pl`, which is pulled in via
`:- include('clues.pl').` at the top of `crossword.pl`. It defines a
single predicate, `clues/1`, whose argument is a list of `[Answer,
Metadata]` pairs:

    clues([
           ['OMEGA POINT',
            _{clue: 'Transcending entropy',
              link: 'http://en.wikipedia.org/wiki/Omega_Point'}],
           ['FLOW',
            _{clue: 'Autotelic activity',
              link: 'http://en.wikipedia.org/wiki/Flow_(psychology)'}],
           ...
          ]).

- **Answer** — the word. Spaces are allowed for multi-word answers; they
  are stripped before the letters are placed on the grid (so
  `'OMEGA POINT'` occupies a single 10-cell run). The answer keeps its
  spaces as a display label in the output.
- **Metadata** — an optional dict of arbitrary data carried through to the
  output verbatim. By convention it holds a `clue` and a `link`, but the
  solver never inspects it, so any keys are fine. The metadata element may
  be omitted entirely (just `['BIAS']`) or given as an empty dict
  (`['FLOW', _{}]`); either yields an empty `meta` object in the output.

Answers must be unique: the output associates metadata with each placed
word by its answer string, so a repeated answer is rejected up front.

To use your own words, replace the contents of `clues/1` and pick a
`grid_length` large enough to lay them out.

### JSON clue file (`--clues`)

For programmatic use you can supply the clue set as a JSON file instead,
without editing any Prolog:

    $ ./crossword.pl --clues puzzle.json 17 topleft_across

The file is a JSON object with a `clues` array; each entry mirrors an
output `words[]` object minus the solver-computed fields:

    {
      "clues": [
        { "answer": "OMEGA POINT",
          "meta": { "clue": "Transcending entropy",
                    "link": "http://en.wikipedia.org/wiki/Omega_Point" } },
        { "answer": "FLOW" },
        { "answer": "BIAS", "meta": {} }
      ]
    }

- **`answer`** — required string. Spaces are allowed and stripped before
  placement, exactly as in `clues.pl`.
- **`meta`** — optional object (default `{}`), copied verbatim to the
  output entry's `meta`. By convention `clue` and `link`, but any keys are
  fine; the solver never inspects it. Omitting `meta` and giving `{}` are
  equivalent.

The two sources are otherwise identical: the JSON path produces the same
internal word list as `clues/1`, so uniqueness checking, solving, and output
are unchanged (the bundled `clues.pl` remains the default when `--clues` is
absent). A missing file, malformed JSON, or a schema violation (no `clues`
array, a non-string `answer`, a non-object `meta`) is reported and the
program exits non-zero. See
[`docs/json-input-spec.md`](docs/json-input-spec.md) for the full schema and
rationale; `tests/clues.json` is a worked example.


## How it works

### Data structures

The solver uses two structures (see the comments at the top of
`crossword.pl`):

1. **The grid** — an association list (`library(assoc)`) keyed by cell
   number. Cells are numbered `1 .. grid_length²` in row-major order
   (left to right, top to bottom). Moving **across** is `+1` (along a
   row); moving **down** is `+grid_length` (down a column). Each value
   is initially the atom `empty`.

2. **Placed words** — a list, one entry per word already laid down. Each
   entry is a dict:

       word{answer:Answer, letters:Letters, cells:Cells, dir:Dir,
            len:Len, start:Start, num:ClueNum}

   where `letters` is the space-stripped character list, `cells` is the
   list of grid-cell numbers it occupies, `dir` is `across`/`down`,
   `len` is the length, `start` is the starting cell, and `num` (the clue
   number) is added at the end by `assign_clue_numbers/2`. The solver
   carries no clue text or link here — metadata is rejoined to each word by
   its answer string only at output time.

### The algorithm

The driver is `assign_words/8`. It is a greedy, backtracking placement:

1. Pick a not-yet-placed word (via `member/2`, so the choice is a
   backtrack point).
2. `find_intersecting_word/6` chooses where to put it. For the very
   first word there is nothing to cross, so the grounded `start_loc` is
   used. For every word after that, it finds an already-placed word that
   shares a letter, computes the crossing cell, flips the direction
   (`swap_dir/2`), derives the candidate start cell, and checks the word
   fits on the grid (`fits_on_grid/4`).
3. `assign_word/9` lays the letters down with `assign_letters/7`,
   enforcing the two layout rules:
   - each cell must either already hold the *same* letter (a legal
     crossing) or be empty; and
   - an empty cell's perpendicular neighbours must be free
     (`adj_is_free/4`), plus the cells just before the start and just
     after the end must be empty (`check_prev_cell/4`,
     `check_next_cell/4`). Together these stop words being laid flush
     alongside one another — only genuine crossings are allowed.
4. Recurse with the word removed from the pool and added to the placed
   list. When the pool is empty, every word has been placed and the
   layout succeeds.

Because the words are tried in order and the solver returns the first
success, a `--shuffle` run randomises both the word order and the order
the start positions are tried (`shuffle/2`, built on `choose/2`).

### Finishing and output

Once all words are placed:

- `assign_clue_numbers/2` sorts the placed words by their start cell and
  assigns clue numbers in that order (cells that start both an across
  and a down word share a number).
- `emit_json/3` builds the grid and word list and writes the whole
  solution as a single JSON object via `json_write_dict/2` — to standard
  output, or to the named file when `--out <file>` is given.

The output is one JSON object with three top-level keys (`json_write_dict`
emits object keys in sorted order, so on the wire they appear as `grid`,
`gridLength`, `words`):

    {
      "gridLength": 17,
      "grid":  [ ...gridLength rows... ],
      "words": [ ...one entry per placed word... ]
    }

- **`grid`** — a dense `grid_length` × `grid_length` array, row-major and
  0-indexed `[row, col]`. Each entry is either `null` (an empty cell) or a
  cell object `{"letter", "number", "across", "down"}`:

  - `letter` — the cell's letter;
  - `number` — the corner-label clue number, or `null` unless a word
    starts in this cell (across and down words starting together share it);
  - `across` / `down` — the clue numbers of the across/down words passing
    through the cell, or `null` for none. A crossing cell carries both, so
    both words (and both links) stay reachable; the consumer decides which
    to follow.

- **`words`** — one object per placed word: `number`, `direction`
  (`"across"`/`"down"`), `answer` (spaces preserved), `cells` (the ordered
  `[row, col]` pairs it occupies), and `meta` (the verbatim metadata dict
  from `clues.pl`, `{}` if none). To follow a link from a cell, match the
  cell's `across`/`down` number **and** direction against `words`.

If no layout exists for the requested `grid_length`, the program produces
no output and exits non-zero (and with `--out`, no file is written).

The full schema and design rationale live in
[`docs/json-output-spec.md`](docs/json-output-spec.md).


## Implementation notes / limitations

- **Grid size is not inferred.** You must pass a `grid_length` that is
  large enough; too small and there is no solution, and the solver will
  simply fail to find one.
- **First-solution semantics.** The solver stops at the first valid
  layout. `--all` exists to enumerate them but is expensive (see the
  usage note).
- **Duplicate solutions.** Many enumerated solutions are the same
  physical layout reached by placing words in a different order.
