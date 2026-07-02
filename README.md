crosswordsmith
===============================================================================

A CLI-first crossword **layout, validation, and fill** toolkit in SWI-Prolog. It
covers the setter's grid workflow end to end — arrange a closed set of words into
a layout, validate a grid against house-style rules, fill a blocked grid from a
dictionary, and export to standard interchange formats — emitting deterministic,
diffable output throughout. It feeds clue-writing tools; it is not one (no
auto-cluing). The full product vision and acceptance criteria live in
[`docs/design-spec.md`](docs/design-spec.md).

Everything is **deterministic**: identical input and flags always produce
byte-identical output — no randomness, no shuffle. One CLI, one verb per
capability:

- **`arrange`** — lay out *these specific* words into an interlocking,
  crossword-shaped layout (*Flavour A* — a closed word set, aesthetic interlock).
  It places words so each crosses an already-placed word and never sits flush
  against another, maximising a capped interlock objective over an N×N grid.
- **`lint`** — validate a layout against a profile (`toc`, `blocked-uk`,
  `american`, `barred-ximenean`): PASS / WARN / FAIL per rule, per word.
- **`fill`** — take a legal blocked grid and fill every slot from a dictionary,
  your own words pinned as seeds (*Flavour B* — grid-first, open-dictionary).
- **`export`** — convert a layout to ipuz v2 or Exolve (onward to `.puz`/`.jpz`/
  PDF via off-the-shelf [kotwords](https://github.com/jpd236/kotwords)).

The four verbs are glued together by one **canonical layout JSON**: `arrange`
and `fill` emit it, `lint` and `export` consume it, and it round-trips back in
as a seed/fragment. A small bundled **stock-grid library** (`grids/`) of
lint-validated blocked templates feeds `lint` and `fill`. The core is portable
Prolog; the clue-input and JSON paths use SWI-specific features (see
Requirements).

The former `./crossword.pl …` interface has been replaced; the implementation
now lives under `prolog/crosswordsmith/` and the root `crossword.pl` is only a
migration-message shim (see
[Migration](#migration-from-the-old-crosswordpl-cli)).

Words — and optional per-word metadata (a clue, a link, anything) — are supplied
to `arrange` as a JSON file or a Prolog `clues/1` fixture. The bundled
`fixtures/bundled_17_clues.pl` attaches a URL to each word (it was originally for
embedding in a web page, where the clues/solved words are links).


## Requirements

[SWI-Prolog](https://www.swi-prolog.org/) 10.1.x. The solver core is plain
Prolog, but the input and output paths are SWI-specific: clue metadata and
the JSON output both use SWI dicts and `library(http/json)`, so porting to
another Prolog would mean replacing those. It is developed and tested against
the SWI development release currently on `PATH` as `swipl` (10.1.10 at the time
of writing).


## Usage

`crosswordsmith` is an executable SWI-Prolog script (note the `#!/usr/bin/env
swipl` shebang, so it follows your `PATH`). A bare invocation prints usage and
exits non-zero; every capability is a verb. If you do not have execute permission
set, run it as `swipl crosswordsmith <args…>`.

    $ ./crosswordsmith arrange --input <file> [options]

### `arrange` — lay out a closed set of words

    # Strict (place every word) on a tight square crop. The bundled set needs a
    # side of 17 (it has a 16-letter answer); --size is a ceiling under the
    # default --size-mode max, so the result auto-shrinks to fit.
    $ ./crosswordsmith arrange --size 17 --input fixtures/bundled_17_clues.pl

    # A fixed 17x17 canvas (uncovered cells emitted as blocks).
    $ ./crosswordsmith arrange --strict --size-mode fixed --size 17 \
        --input fixtures/bundled_17_clues.pl

    # Best-effort: place a maximal subset on a tight grid; report the dropped
    # words on stderr.
    $ ./crosswordsmith arrange --best-effort --size 11 \
        --input fixtures/bundled_17_clues.pl

    # Seed from a partial layout (anchors): pin some words, let the engine
    # finish. The fragment's gridLength sets the size.
    $ ./crosswordsmith arrange --fragment fixtures/bundled_17_fragment.json \
        --input fixtures/bundled_17_clues.pl

    # Up to K meaningfully-distinct layouts, emitted as a JSON array.
    $ ./crosswordsmith arrange --candidates 3 --size 17 \
        --input fixtures/bundled_17_clues.pl

    # Count every feasible full placement (the old `--all`).
    $ ./crosswordsmith arrange --enumerate --size 17 \
        --input fixtures/bundled_17_clues.pl

    # Write to a file instead of stdout (composes with any form above; the file
    # is written only on success, so a failed run leaves no partial file).
    $ ./crosswordsmith arrange --size 17 --input <file> --out solution.json

**Options** (parsed with `library(optparse)`, so they compose in any order and
accept the `--flag=value` form):

| flag | meaning |
| --- | --- |
| `--input <file>` | **required** — the word/clue set (`.json` or `.pl`). |
| `--strict` | fail unless every word is placed (the **default**). |
| `--best-effort` | place a maximal subset; report dropped words on stderr. |
| `--size <N>` | square grid side (default `15`; a *ceiling* under `--size-mode max`). |
| `--size-mode fixed\|max` | `fixed` = exact N×N (blocks for empty cells); `max` = tight enclosing-square crop (the **default**). |
| `--fragment <file>` | seed from a partial-layout fragment (JSON — the emit format made partial). Its `gridLength` sets `N`; `--size` is then redundant, and an error if it disagrees. |
| `--candidates <K>` | emit up to `K` diverse layouts as a JSON array. Returns fewer than `K` (reported on stderr) when fewer ≥τ-distinct layouts exist. |
| `--enumerate` | count every feasible full placement instead of emitting a layout. |
| `--out <file>` | write output to `<file>` instead of stdout. |
| `--help` / `-h` | print the arrange options. |

`--strict` and `--best-effort` are mutually exclusive; `--enumerate` does not
combine with `--candidates`/`--fragment`, and `--candidates` does not combine
with `--fragment` (v1). The output is **deterministic**: identical input + flags
⇒ byte-identical JSON. `--strict` resolves to exactly one of — a full legal
layout (exit 0); a non-zero exit naming a genuinely unplaceable word; or a
non-zero *"not proven within budget"*. On any failure with `--out`, no partial
file is written.

> **Note on `--enumerate`:** it counts *every* solution by backtracking, and
> many are the same physical layout reached by placing words in a different
> order, so the count is large and the search can be slow. Use it on small
> grids / small clue sets.

### Migration from the old `crossword.pl` CLI

The previous `./crossword.pl --input F <N> <loc>` interface is replaced by the
`crosswordsmith` subcommands — a one-time breaking change:

| old | new |
| --- | --- |
| `crossword.pl --input F <N> <loc>` | `crosswordsmith arrange --strict --size-mode fixed --size <N> --input F` |
| `crossword.pl --input F --quality` | `crosswordsmith arrange --best-effort --size-mode max --input F` |
| `crossword.pl --input F --all <N>` | `crosswordsmith arrange --enumerate --size <N> --input F` |

`--shuffle` is removed (output is now deterministic). `--strategy`/`--start_loc`
are gone (the engine uses the production strategy internally and sweeps start
corners itself). Running `./crossword.pl` directly — or `crosswordsmith` with
old-style arguments — prints this mapping.

### `lint` — validate a layout against a profile

`lint` consumes a **canonical layout** (exactly what `arrange` emits) and reports
**PASS / WARN / FAIL per rule, per word**, plus a summary verdict, under a named
profile. It is a validator — no engine — so it works on any canonical layout,
including hand-authored ones.

    # Validate a saved layout under the relaxed (advisory) default profile.
    $ ./crosswordsmith lint --profile toc layout.json

    # Pipe arrange straight into lint under a strict profile.
    $ ./crosswordsmith arrange --size-mode fixed --size 17 \
        --input fixtures/bundled_17_clues.pl --out layout.json
    $ ./crosswordsmith lint --profile blocked-uk layout.json

| flag | meaning |
| --- | --- |
| `--profile <name>` | **required** — `toc` (advisory-only), `blocked-uk`, `american`, or `barred-ximenean` (the Ximenean per-length unchecked-letter band). |
| `--allow-asymmetry` | symmetry never hard-FAILs (its deficit is downgraded to WARN). |
| `--out <file>` | write the report to `<file>` instead of stdout. |
| `--help` / `-h` | print the lint options. |

The report is deterministic JSON: a `verdict` (`PASS`/`WARN`/`FAIL`), a `summary`
count, grid-level results (`connectivity`, `symmetry`), and per-word `results`
(min length, checked fraction, longest unchecked run, double-unchecked ends,
odd/even balance — exactly the rule set the chosen profile defines). The report
is always emitted; the **exit code** is non-zero only when a FAIL-severity rule
trips under the chosen profile.

### `export` — convert a layout to a standard format

`export` transforms a canonical layout into a standard interchange format (it is
a transformation, not a new emitter):

    # ipuz v2 (JSON). From here, reach .puz/.jpz/PDF via off-the-shelf kotwords.
    $ ./crosswordsmith export --to ipuz   layout.json > puzzle.ipuz

    # Exolve plain text (git-diffable; round-trips to Exet).
    $ ./crosswordsmith export --to exolve layout.json > puzzle.exolve

| flag | meaning |
| --- | --- |
| `--to ipuz\|exolve` | **required** — the target format. |
| `--out <file>` | write to `<file>` instead of stdout. |
| `--help` / `-h` | print the export options. |

Enumerations (e.g. `(5,5)`, `(4-5)`) are derived from the spaces and hyphens in
each answer; clue text rides through from each word's `meta.clue`. Nothing is
invented — a word with no clue exports an empty clue. (Spec-valid ingestion by a
third party — kotwords for ipuz, Exet for Exolve — is the intended consumer; that
round-trip is a manual verification step.)

### `fill` — grid-first auto-fill

`fill` takes a **legal blocked grid** (a stock-grid mask under `grids/`, or your
own) and fills every slot with a dictionary word so that crossings agree, with
your own words pinned as **seeds**. Output is a canonical layout, so it composes
with `lint`/`export`.

    # Fill a grid from a word list (a tiny sample ships; real fills want UKACD18).
    $ ./crosswordsmith fill --grid grids/blocked_13a.json --dict UKACD18.txt

    # Pin some answers (a fragment, the same format arrange consumes) and fill
    # around them.
    $ ./crosswordsmith fill --grid grids/blocked_13a.json --seeds seeds.json \
        --dict UKACD18.txt

| flag | meaning |
| --- | --- |
| `--grid <file>` | **required** — the grid template (a `grids/` black-square mask). |
| `--seeds <file>` | seed words to pin (a fragment, §6.6); filled around as hard pins. |
| `--dict <file>` | word list, one per line (default: a small bundled sample; real fills: `--dict UKACD18`). |
| `--out <file>` | write to `<file>` instead of stdout. |

Each white cell is a shared logical variable, so crossings are consistent by
construction; the search is MRV backtracking (most-constrained slot first) over
an in-memory pattern index, and is **deterministic**. When no complete fill
exists, `fill` reports the unfillable slot(s) and exits non-zero (it never emits
a partial grid). The bundled lexicon is a tiny sample for the demo grids; supply
a real dictionary (UKACD18 — redistributable freeware; ship its license verbatim) with `--dict` for production fills.


## Word/clue input (for `arrange`)

`arrange`'s `--input` is a word/clue set. The bundled example lives in
`fixtures/bundled_17_clues.pl`. Prolog fixture files define a single `clues/1`
term whose argument is a list of `[Answer, Metadata]` pairs:

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

To use your own words with the main CLI, pass either a JSON clue file or a
Prolog fixture file with `--input`.

### JSON clue file

For programmatic use you can supply the clue set as a JSON file instead,
without editing any Prolog:

    $ ./crosswordsmith arrange --size 17 --input puzzle.json

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
  placement, exactly as in the bundled Prolog fixture.
- **`meta`** — optional object (default `{}`), copied verbatim to the
  output entry's `meta`. By convention `clue` and `link`, but any keys are
  fine; the solver never inspects it. Omitting `meta` and giving `{}` are
  equivalent.

The two sources are otherwise identical: JSON files and Prolog fixtures
produce the same internal word list as `clues/1`, so uniqueness checking,
solving, and output are unchanged. A missing file, malformed JSON, Prolog
fixture without `clues/1`, unsupported extension, or JSON schema violation
(no `clues` array, a non-string `answer`, a non-object `meta`) is reported
and the program exits non-zero. See
[`docs/json-input-spec.md`](docs/json-input-spec.md) for the full schema and
rationale; `tests/clues.json` is a worked example.


## How it works

### Modules

One shared substrate, two solver tops, with the Flavour-B tools as
transformations/validators over the canonical JSON:

The implementation lives under `prolog/crosswordsmith/`; the root `load.pl`
loads it all in the right order (the CLI, tests, and benchmarks go through it):

| file | role |
| --- | --- |
| `prolog/crosswordsmith/core.pl` | **shared substrate** — grid model, the free-canvas legality core, clue numbering, JSON emit + input loading. |
| `prolog/crosswordsmith/metrics.pl` | shared **metric predicates** (checked cells, unchecked runs, bbox, crossings) + the greedy density constructor that `arrange` reuses. |
| `prolog/crosswordsmith/arrange.pl` | **Flavour A** — the deterministic MRV-first layout engine: construct + rescore + emit, with fragment seeding and diverse candidates. |
| `prolog/crosswordsmith/lint.pl` | **Flavour B** — the profile-driven grid validator (consumes the canonical layout, reuses the metric predicates). |
| `prolog/crosswordsmith/export.pl` | **Flavour B** — ipuz v2 / Exolve transforms of the canonical layout. |
| `prolog/crosswordsmith/stockgrid.pl` + `grids/` | **Flavour B** — the bundled stock-grid library: black-square masks, slots derived on load, each validated by `lint --profile blocked-uk`. |
| `prolog/crosswordsmith/fill.pl` | **Flavour B** — grid-first auto-fill: each white cell a shared logical variable, MRV backtracking over an in-memory dictionary index, seeds pinned. |
| `load.pl` | loads the implementation in the known-good order; defines the `crosswordsmith` file-search alias. |
| `crosswordsmith` | the CLI: verb dispatch (`arrange`/`lint`/`export`/`fill`). |
| `crossword.pl` | migration-message shim for the old CLI (prints the mapping below, exits non-zero). |

The rest of this section describes the shared substrate and the `arrange`
engine; the `lint`/`export`/`fill` sections above describe their own behaviour.

### Data structures

The substrate uses two structures (see the comments at the top of
`prolog/crosswordsmith/core.pl`):

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

This describes the shared **legality core** (`prolog/crosswordsmith/core.pl`),
which `arrange`
reuses. `arrange` drives it with a deterministic most-constrained-first
(MRV) ordering, constructs over the four start corners, rescores each complete
layout by a capped interlock objective, and emits the best — so the output is
stable and shuffle-free (the old random `--shuffle` path was removed).

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
  from the input clue set, `{}` if none). To follow a link from a cell, match
  the cell's `across`/`down` number **and** direction against `words`.

If no layout exists for the requested `grid_length`, the program produces
no output and exits non-zero (and with `--out`, no file is written).

The full schema and design rationale live in
[`docs/json-output-spec.md`](docs/json-output-spec.md).


## Implementation notes / limitations

- **Grid size is not inferred.** You pass `--size N` (a ceiling under
  `--size-mode max`). Too small for the words and `--strict` reports the
  failure and exits non-zero rather than silently mangling the layout; under
  `--best-effort` the unplaceable words are dropped and reported.
- **Best-within-budget, not proven-optimal.** `arrange` constructs a strong
  layout and rescores it; it does not exhaustively prove the optimum (the
  search superstructure was descoped after measurement — see
  [`docs/arrange-implementation-plan.md`](docs/arrange-implementation-plan.md)).
- **`--enumerate` is expensive.** Many enumerated solutions are the same
  physical layout reached by placing words in a different order; the count is
  large and the search slow on big inputs.
- **`fill` needs a dictionary.** Only a tiny sample wordlist ships (enough for
  the demo grids); pass a real lexicon (UKACD18 — redistributable freeware,
  ship its license verbatim) with `--dict` for production fills. The full
  dictionary is not bundled.
- **`export`'s third-party round-trip is a manual step.** The output is
  spec-valid ipuz v2 / Exolve by construction, but actual ingestion by kotwords
  (ipuz) or Exet (Exolve) is verified by hand, not in CI.
- **The stock-grid library is a curated starter set.** A handful of
  lint-validated blocked grids under `grids/`, not a grid *generator*; it grows
  by adding more validated masks.
- **No barred bar/edge model.** The `barred-ximenean` lint profile applies the
  Ximenean per-length checking math to a canonical (blocked-shaped) layout's
  words; a dedicated barred-grid representation is out of scope.


## Development

Run the full regression suite with:

    $ make test

For local performance comparisons, run the benchmark harness:

    $ make bench

or directly:

    $ swipl -q benchmarks/run_benchmarks.pl -- fixtures/bundled_17_clues.pl

The harness reports min, median, and mean wall time, CPU time, and inference
counts for a Prolog fixture file that defines `clues/1`. These numbers are
machine-specific and are intended for comparing local branches, not as CI
thresholds or authoritative baselines. The default fixture is
`fixtures/bundled_17_clues.pl`; pass a different fixture path directly, or use
`make bench BENCH_FIXTURE=fixtures/other.pl BENCH_GRID=<n>`. Use `--grid <n>`,
`--start-loc <loc>`, `--iterations <n>`, `--warmup <n>`, and
`--format text|csv|json` to shape a direct run. With `make bench`, the matching
variables are `BENCH_GRID`, `BENCH_START_LOC`, `BENCH_ITERATIONS`,
`BENCH_WARMUP`, and `BENCH_FORMAT`. The default format is `text`.

Synthetic benchmark fixtures are included for deeper search workloads:

| fixture | words | grid |
| --- | ---: | ---: |
| `fixtures/benchmark_08_words.pl` | 8 | 13 |
| `fixtures/benchmark_14_words.pl` | 14 | 17 |
| `fixtures/benchmark_16_dense_words.pl` | 16 | 17 |
| `fixtures/benchmark_20_words.pl` | 20 | 37 |
| `fixtures/benchmark_26_words.pl` | 26 | 49 |

For example:

    $ make bench BENCH_FIXTURE=fixtures/benchmark_20_words.pl BENCH_GRID=37

The dense 16-word fixture is intentionally much slower than the others; start
with one measured iteration:

    $ make bench BENCH_FIXTURE=fixtures/benchmark_16_dense_words.pl BENCH_GRID=17 BENCH_ITERATIONS=1 BENCH_WARMUP=0
