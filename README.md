crosswordsmith
===============================================================================

A CLI-first crossword **layout, validation, and fill** toolkit in SWI-Prolog. It
covers the setter's grid workflow end to end — arrange a closed set of words into
a layout, validate a grid against house-style rules, fill a blocked grid from a
dictionary, and export to standard interchange formats — emitting deterministic,
diffable output throughout. It feeds clue-writing tools; it is not one (no
auto-cluing). The full product vision and acceptance criteria live in
[`docs/design-spec.md`](docs/design-spec.md).

Everything is **deterministic by default**: identical input and flags always
produce byte-identical output — no randomness, no shuffle. (`arrange` and
`fill` take opt-in `--seed N` / `--shuffle` flags for a reproducible or fresh
pseudo-random result; neither ever touches the default deterministic path —
see below.) One CLI, one verb per capability:

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

The implementation lives under `prolog/crosswordsmith/`, driven by the
`crosswordsmith` CLI.

The same engine also runs **client-side in the browser** (SWI-Prolog compiled
to WASM, in a Web Worker) behind a typed JS SDK — `arrange`, `lint`, and
`export` today, value-locked against the CLI's output. See
[`wasm/README.md`](wasm/README.md). The measured plan to reduce its browser
payload and startup cost is
[`docs/plans/wasm-payload-performance.md`](docs/plans/wasm-payload-performance.md).

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

    # Strict (place every word) on an exact 17x17 grid — the bundled set has a
    # 16-letter answer, so it needs a side of at least 17. Uncovered cells are
    # emitted as blocks. (--size is the default framing; --strict is the default.)
    $ ./crosswordsmith arrange --size 17 --input fixtures/bundled_17_clues.pl

    # --max-size instead builds up to NxN and crops to the tight enclosing
    # square, so the result auto-shrinks to fit the placed words.
    $ ./crosswordsmith arrange --max-size 17 \
        --input fixtures/bundled_17_clues.pl

    # Best-effort: place a maximal subset, cropped to a tight grid; the dropped
    # words are recorded in the output's diagnostics property.
    $ ./crosswordsmith arrange --best-effort --max-size 11 \
        --input fixtures/bundled_17_clues.pl

    # Seed from a partial layout (anchors): pin some words, let the engine
    # finish. The fragment's gridLength sets the size.
    $ ./crosswordsmith arrange --fragment fixtures/bundled_17_fragment.json \
        --input fixtures/bundled_17_clues.pl

    # The same pins in the thin hand-authorable form - a JSON list of
    # {answer, row, col, dir} (0-based start cell). No gridLength, so
    # --size/--max-size frames it (default 15). Identical output (AC-FRAG-4).
    $ ./crosswordsmith arrange --fragment fixtures/bundled_17_fragment_thin.json \
        --size 17 --input fixtures/bundled_17_clues.pl

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
| `--best-effort` | place a maximal subset; the dropped words are recorded in the output's `diagnostics` property. |
| `--size <N>` | exact grid side — an N×N grid, uncovered cells emitted as blocks (default `15`, a standard crossword side). Mutually exclusive with `--max-size`. |
| `--max-size <N>` | build up to N×N, then crop to the tight enclosing square (`≤ N`). Mutually exclusive with `--size`. |
| `--fragment <file>` | seed from a partial-layout fragment (JSON), in either §6.6 form. **Canonical** (the emit format made partial): its `gridLength` sets `N`; `--size`/`--max-size` are then redundant, and an error if they disagree. **Thin** (hand-authorable): a top-level list of `{answer, row, col, dir}` entries (0-based start cell, `dir` = `across`/`down`); it carries no `gridLength`, so `--size`/`--max-size` frames it (default 15). Both produce identical results for the same pins (AC-FRAG-4). |
| `--candidates <K>` | emit up to `K` diverse layouts as a JSON array. Returns fewer than `K` (reported on stderr) when fewer ≥τ-distinct layouts exist. |
| `--enumerate` | count every feasible full placement instead of emitting a layout. |
| `--out <file>` | write output to `<file>` instead of stdout. |
| `--verbose` | report the success summary (grid, placed, dropped, reward, cap status) on stderr; by default a clean success is silent there — quality caveats live in the output's `diagnostics` property instead. Fewer-than-K-candidates warnings and failures print regardless. |
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

### `lint` — validate a layout against a profile

`lint` consumes a **canonical layout** (exactly what `arrange` emits) and reports
**PASS / WARN / FAIL per rule, per word**, plus a summary verdict, under a named
profile. It is a validator — no engine — so it works on any canonical layout,
including hand-authored ones.

    # Validate a saved layout under the relaxed (advisory) default profile.
    $ ./crosswordsmith lint --profile toc layout.json

    # Pipe arrange straight into lint under a strict profile.
    $ ./crosswordsmith arrange --size 17 \
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

    # Scored fill out of the box: the bundled Collaborative Word List
    # clean-floor derivative (dicts/cwl50.dict, MIT — see dicts/README.md)
    # on the bundled American-style 11x11. High-score-first, ~6s.
    $ ./crosswordsmith fill --grid grids/amer11.json --dict dicts/cwl50.dict

    # UK-style grids want a UK lexicon (the bundled CWL is American-style
    # and does not fill the blocked_* stock grids): supply your own, e.g.
    # UKACD18 (redistributable freeware).
    $ ./crosswordsmith fill --grid grids/blocked_13a.json --dict UKACD18.txt

    # Pin some answers (a §6.6 fragment, canonical or thin form) and fill
    # around them. A thin seed file is just [{"answer": "COW", "row": 0,
    # "col": 0, "dir": "across"}, ...], framed by the grid itself.
    $ ./crosswordsmith fill --grid grids/blocked_13a.json --seeds seeds.json \
        --dict UKACD18.txt

    # Scored fill: point --dict at a word;score list (e.g. Spread the
    # Wordlist - CC BY-NC-SA, download it yourself; it is never bundled).
    # Candidates are tried high-score-first; --min-score 50 is the usual
    # "clean" floor for 0-100-scale lists.
    $ ./crosswordsmith fill --grid grids/blocked_13a.json \
        --dict spread-the-wordlist.txt --min-score 50 --report-json quality.json

| flag | meaning |
| --- | --- |
| `--grid <file>` | **required** — the grid template (a `grids/` black-square mask). |
| `--seeds <file>` | seed words to pin (a fragment, §6.6, canonical or thin form — a thin `[{answer,row,col,dir}]` list is framed by the grid itself, no `gridLength` needed); filled around as hard pins. |
| `--dict <file>` | word list, one per line, **UTF-8** (default: a small bundled sample; for real fills pass `dicts/cwl50.dict` — the bundled scored lexicon, opt-in by design — or your own list, e.g. UKACD18). Words are normalized to A–Z: accented Latin letters fold to their base (café → CAFE, Straße → STRASSE), punctuation/digits are squeezed; a word with letters that cannot be folded (Cyrillic, Greek, …) is dropped, and the drop count is reported on stderr (unconditionally — pure-ASCII lists load in silence). A file containing `;` is read as a **scored** list (`word;score`, integer scores in the list's **own** units — no normalisation): candidates are tried score-descending (band order strict; the equal-score tail order is deterministic — lexicographic for multi-band lists, the engine's pinned shuffle when every word scores the same), unscored lines score 1, malformed lines are dropped + reported, and a duplicate word keeps its highest score. |
| `--min-score <N>` | hard-prune every candidate scoring below `N` (native units) *before* search; default `1` drops only score-0 blocklist entries, so unscored dictionaries are unaffected. `50` is the documented clean floor for 0–100-scale lists (STW/Broda). Pruning shrinks slot domains, so it can only *reduce* feasibility — a hard grid may lose its fill (reported as the ordinary no-fill outcome). Any nonzero prune count is reported on stderr. Text `--dict` path only (index artifacts carry no scores). |
| `--report-json <file>` | write the fill-quality report for the produced fill to `<file>` as one sorted-key JSON object: `{"belowThreshold":…,"mean":…,"min":…,"n":…,"threshold":50}` (threshold 50 = the clean-floor convention; entries absent from the dict score 0). The stdout layout is byte-unchanged; nothing is written when no fill is produced. Text `--dict` path only. |
| `--budget <N>` | override the search's inference budget (default 800,000,000). Deterministic: a budget change never alters a *produced* fill, only fill-vs-not-proven — an **escape hatch** for marginal grids (the default covers every bundled benchmark rung and the reference hard 13×13; [`benchmarks/fill_quality/`](benchmarks/fill_quality/README.md)). Composes with every other flag and both index modes. |
| `--seed <N>` | (N ≥ 0) a *reproducible alternative fill* of the same grid: every restart attempt picks among the top-3 candidates (weighted 4:2:1) on a PRNG stream seeded by `N`, and equal-score candidate order is seed-perturbed at load — the score-descending band order is untouched, so quality stays score-first. Variety with real per-seed variance (a seed can complete a marginal grid faster *or* exhaust the budget where the default completes). Equal-score bands use an O(n log n) order-statistic replay that preserves the historical draw sequence and exact permutation. Text `--dict` path only; mutually exclusive with `--shuffle`. |
| `--shuffle` | a *fresh* random alternative fill each run (same scope as `--seed`). Recoverable: `--verbose` prints `fill: shuffle seed N (reproduce with --seed N)` on stderr (the fill payload carries no diagnostics). Text `--dict` path only. |
| `--out <file>` | write to `<file>` instead of stdout. |
| `--verbose` | report the success summary (grid, filled slots) and the fill-quality line (`n/mean/min/below50`) on stderr; a clean success is silent there by default. Failures print regardless. |

Each white cell is a shared logical variable, so crossings are consistent by
construction. The search (design-spec §8.4c) keeps per-slot candidate sets as
bitmask domains over the pattern index, propagates crossing constraints to an
arc-consistency fixpoint after every placement, picks the next slot by
dom/wdeg conflict-weight learning, and restarts under a growing node cap —
and it is **deterministic**: the default path is a pure function of the
input (its diversification runs on pinned engine constants; no OS entropy,
byte-identical CLI and browser), while opt-in `--seed`/`--shuffle` swap in a
user-seeded stream for reproducible/fresh variety. When no complete fill
exists, `fill` reports the unfillable slot(s) and exits non-zero (it never emits
a partial grid). The `--dict` default is a tiny sample for the demo grids;
production fills pass a real lexicon explicitly — the bundled scored
`dicts/cwl50.dict` (Collaborative Word List clean-floor derivative, MIT,
provenance in [dicts/README.md](dicts/README.md)) or your own (UKACD18 —
redistributable freeware — for UK-style grids).


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

The implementation lives under `prolog/crosswordsmith/`, one SWI-Prolog module
per file (`crosswordsmith_core`, `crosswordsmith_metrics`, ... — a
`crosswordsmith_` prefix because module names share a flat namespace) with
explicit export lists; the root `load.pl` loads it all (the CLI, tests, and
benchmarks go through it):

| file | role |
| --- | --- |
| `prolog/crosswordsmith/core.pl` | **shared substrate** — grid model, the free-canvas legality core, clue numbering, JSON emit + input loading. |
| `prolog/crosswordsmith/metrics.pl` | shared **metric predicates** (checked cells, unchecked runs, bbox) consumed by `arrange` as optimizer signals and by `lint` as validators. |
| `prolog/crosswordsmith/arrange.pl` | **Flavour A** — the deterministic MRV-first layout engine: construct + rescore + emit, with fragment seeding, diverse candidates, and the greedy density constructor. |
| `prolog/crosswordsmith/lint.pl` | **Flavour B** — the profile-driven grid validator (consumes the canonical layout, reuses the metric predicates). |
| `prolog/crosswordsmith/export.pl` | **Flavour B** — ipuz v2 / Exolve transforms of the canonical layout. |
| `prolog/crosswordsmith/stockgrid.pl` + `grids/` | **Flavour B** — the bundled stock-grid library: black-square masks, slots derived on load, each validated by `lint --profile blocked-uk`. |
| `prolog/crosswordsmith/fill.pl` | **Flavour B** — grid-first auto-fill: each white cell a shared logical variable, MAC + dom/wdeg + restarts over bitmask domains (design-spec §8.4c), seeds pinned. |
| `load.pl` | loads the implementation in the known-good order; defines the `crosswordsmith` file-search alias. |
| `crosswordsmith` | the CLI: verb dispatch (`arrange`/`lint`/`export`/`fill`). |

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
(MRV) ordering, constructs over two non-transpose start-corner
representatives (`topleft_across`, `topright`) under one shared operation-wide
inference budget, and emits the best — so by default the
output is stable and shuffle-free.

For variety, two opt-in flags perturb only the branch **ordering** of the same
MRV search (shuffling the seed word and shuffling within equal-constraint
buckets, so completeness and solution quality are unchanged):

- **`arrange --seed N`** (N ≥ 0) — a *reproducible* pseudo-random layout: a
  given `N` always reproduces the same grid.
- **`arrange --shuffle`** — a *fresh* random layout on every run (seeded from
  OS entropy). It stays recoverable: `--verbose` prints the seed it drew
  (`arrange: shuffle seed 508756889 (reproduce with --seed 508756889)`), so a
  layout you like is never lost. Mutually exclusive with `--seed`.

Both live entirely on this opt-in path — with neither flag, no RNG is seeded or
consulted, so the deterministic default is byte-for-byte untouched. They apply
to the strict search (with or without `--fragment`) and are rejected with
`--best-effort`, `--candidates`, or `--enumerate`, which don't route through
that seam. (Use `--candidates K` for *deterministic* diverse layouts.)

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

- **`diagnostics`** *(optional; `arrange` output only)* — quality caveats
  about how the layout was produced, under a per-producer key:
  `diagnostics.arrange` carries `capInert` (no placed word reached its
  checking target — the quality objective degenerated to plain
  total-crossings; tune with `--check-target`), `dropped` (best-effort's
  unplaced answers, `[]` when everything placed), and `reward` (the engine's
  objective value). Arrange output is best-effort by nature, so these ride
  the payload rather than stderr; consumers may ignore the property entirely.

If no layout exists for the requested `grid_length`, the program produces
no output and exits non-zero (and with `--out`, no file is written).

The full schema and design rationale live in
[`docs/json-output-spec.md`](docs/json-output-spec.md).


## Implementation notes / limitations

- **Grid size is not inferred.** You pass `--size N` (an exact N×N grid) or
  `--max-size N` (a ceiling, cropped to fit). Too small for the words and `--strict` reports the
  failure and exits non-zero rather than silently mangling the layout; under
  `--best-effort` the unplaceable words are dropped and recorded in the
  output's `diagnostics`.
- **Best-within-budget, not proven-optimal.** `arrange` constructs a strong
  layout and rescores it; it does not exhaustively prove the optimum (the
  search superstructure was descoped after measurement — see
  [`docs/arrange-implementation-plan.md`](docs/arrange-implementation-plan.md)).
- **`--enumerate` is expensive.** Many enumerated solutions are the same
  physical layout reached by placing words in a different order; the count is
  large and the search slow on big inputs.
- **`fill`'s default dictionary is a tiny sample** (enough for the demo
  grids — the default deliberately never changed, design-spec §10 DP-9);
  production fills pass a lexicon explicitly. A **scored list is bundled**:
  `dicts/cwl50.dict`, the MIT Collaborative Word List clean-floor derivative
  (252,200 `word;score` entries ≥ 50, snapshot pinned + regenerable —
  [dicts/README.md](dicts/README.md)). Caveats, documented not hidden: CWL
  is American-style (it does not fill the UK `blocked_*` stock grids — pair
  it with `grids/amer11.json`; UK fills want your own UKACD18) and its
  upstream data has been frozen since 2023-02. Lists much beyond the
  bundled size hit measured engine capacity: the full 567k-entry CWL
  exceeds SWI's default 1GB stack at dict load, and its ≥30 band (437k)
  can exceed it in search on full-length-slot grids. Since DP-10
  (design-spec §8.4 AC-FILL-15) both report as an ordinary one-line
  failure naming the limit and the remedies — shrink/`--min-score`, or
  relaunch as `swipl --stack-limit=4g ./crosswordsmith …` — instead of a
  raw stack dump; the capacity envelope itself is unchanged (measured
  bounds in [`benchmarks/fill_quality/`](benchmarks/fill_quality/README.md)).
- **`fill` quality needs a scored dictionary.** With a plain (unscored)
  list the search places any *legal* word — it can pick obscure/junk
  entries a scored filler would reject (measured: it fills a small open
  grid with `AAAAA` where a scored filler yields real words). A scored `--dict` + `--min-score`
  (design-spec §8.4a) closes that quality gap — the bundled
  `dicts/cwl50.dict` needs no flag at all (its ≥ 50 floor is baked in),
  and on it the §8.4c engine *beats* `ingrid_core` on every completable
  benchmark mask (means 78–83 vs 77–79). Stricter floors: `--min-score
  75`/`90`. Other scored sources stay user-supplied: Spread the Wordlist
  (CC BY-NC-SA, never bundled), etc. The once-documented
  completion gap against the closest competitor is **closed**: the §8.4c
  search core (MAC propagation + dom/wdeg + diversified restarts, adopted
  at design-spec §10 DP-8) fills the reference blocked 13×13 that used to
  budget-exhaust — `--min-score 30` in ~2½ min and `--min-score 1` in ~20s
  under the default budget, at a mean fill score *above* `ingrid_core`'s on
  the same row (45.0 vs 44.4, score-first order intact; measured with a
  user-supplied Spread the Wordlist — the row is not completable from the
  bundled CWL, whose loadable floors are too high for it). Two honest edges
  remain: a `--min-score` prune shrinks slot domains, so high thresholds
  still make hard grids *less* fillable, never faster; and the two
  hardest benchmark grids (`blocked_13b`/`blocked_15a`) stay out of reach —
  they defeat `ingrid_core` too, and completion there was never promised.
  Both axes are quantified in
  [`benchmarks/fill_quality/`](benchmarks/fill_quality/README.md).
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

The Python conversion companion under `xword/` has its own suite, not run by
`make test`:

    $ make test-xword

and `make xword-parity` byte-diffs the engine's ipuz/exolve exports against
xword's over one layout (best-effort parity — see `docs/xword-spec.md` §14).

There are two benchmarks, answering different questions. Both share the
measurement core in `benchmarks/bench_core.pl` (warmup, iterate, summarize;
one median definition). All numbers are machine-specific and reporting-only;
compare on **inferences**, which are deterministic and portable, not on wall time.

**Product benchmark** — how long `arrange` takes a user, for the workloads in
`benchmarks/workloads.pl`: a synthetic mesh **cost ladder** (the hill-climbing
instrument, 9×9/15×15/21×21) plus **real-word realism anchors** at the blocked
daily sizes 13/15 (`fixtures/real_*`, ENABLE dictionary words planted on a legal
witness layout by `benchmarks/gen_real_fixture.py`):

    $ make bench

For each workload it reports three layers over the same word set: `cmd_wall`
(end-to-end `crosswordsmith arrange` process latency), `search` (the in-process
two-representative strict search alone), and `rest` (the CLI-wrapper overhead
between them, `cmd_wall - search`). Core workloads run by default; the heavy
tail — the hard ladder rungs (subsecond to a few seconds each) plus one
budget-saturating **latency probe**
(`benchmark_16_dense`, ~30–40 s per layer; its inference count pins to the
budget constant, so it is reported latency-only and never gates) — is opt-in:

    $ make bench BENCH_ARGS=--heavy
    $ make bench BENCH_FORMAT=csv BENCH_ARGS="--fixture real"

`rss` is the whole-process peak footprint, not a search-memory metric. Note that
`arrange` latency is bimodal: strict mode tries two non-transpose start-corner
representatives under a shared inference budget, so an input whose non-placing
corners fail fast is ~0.1-1 s,
while one whose corner triggers deep search burns the whole budget (~26-40 s)
even after a valid layout was already found — the latency probe exists to keep
that cliff visible.

The tracked metric is **search inferences** (deterministic, machine-independent;
the same count native or under WASM). `benchmarks/baseline.json` records each
rung as a **ratchet**: `make bench-check` fails on a rise past tolerance, a drop
is a win you accept with `make bench-record` (which also appends the run to the
`benchmarks/history.jsonl` trend ledger; `make bench-history` renders it).
`make bench-arrange-verify` runs the native tests, diagnostics-bearing ladder
identity, and ratchet in one fail-fast flow. `make bench-arrange-promote` checks
and records one measurement (rather than rerunning an accepted result), then
read-verifies every persisted measured rung; add `BENCH_ARGS=--heavy` to either
target for the complete core+heavy selection.

**Strategy matrix** — algorithm research (comparing the solver strategies, not a
shipped path). One CSV row per (strategy, fixture) over `benchmarks/fixtures.pl`,
each on its manifest grid:

    $ make bench-matrix
    $ make bench-matrix BENCH_STRATEGIES="baseline mrv_capped"

The matrix measures the single-corner search directly, so its per-fixture numbers
are lighter than the product command (which runs both strict representatives).
The synthetic fixtures it drives:

| fixture | words | grid |
| --- | ---: | ---: |
| `fixtures/benchmark_08_words.pl` | 8 | 13 |
| `fixtures/benchmark_14_words.pl` | 14 | 17 |
| `fixtures/benchmark_16_dense_words.pl` | 16 | 17 |
| `fixtures/benchmark_20_words.pl` | 20 | 37 |
| `fixtures/benchmark_26_words.pl` | 26 | 49 |
