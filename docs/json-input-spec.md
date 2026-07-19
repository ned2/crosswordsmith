# Spec: JSON clue-input file

Status: **implemented** — `load_clues/2` and `doc_to_words/2` live in
`prolog/crosswordsmith/core.pl`; `tests/clues.json`, the `json_input` plunit
suite in `tests/core.plt`, and the README cover the contract.

Companion to [`json-output-spec.md`](json-output-spec.md), which is implemented.

## 1. Problem statement

The repository's example input lives in `fixtures/bundled_17_clues.pl`, a
Prolog fixture file defining `clues/1` as a list of `[Answer, Metadata]` pairs
(`Metadata` an SWI dict). `main/0` now requires `--input File` and loads either
a `.pl` fixture or a `.json` clue file at runtime.

Originally the only way to supply a puzzle was to author Prolog source
compiled into the script. The output is already JSON (`json-output-spec.md`),
and the program was built to feed a web page, so an external consumer that can
*read* the solution should also be able to *produce* a puzzle without writing
Prolog. We accept clues from JSON and from Prolog fixture files through the
same mandatory CLI input-file option.

## 2. Goals

- **G1.** Accept a clue set from an external **JSON file** at runtime.
- **G2.** Accept Prolog fixture files at runtime, so the in-repo dataset does
  not need to be compiled into the implementation.
- **G3.** Add **no new dependency** — reuse `library(http/json)`, already
  loaded for output.
- **G4.** Leave the solver and emitter untouched: both input sources converge
  on the existing internal representation before the pipeline runs.
- **G5.** Make the input schema the natural mirror of the output, so a
  solution payload can be reduced back to a valid input.

## 3. Non-goals

- Moving `grid_length` / `start_loc` into the file (they stay CLI arguments;
  see §11).
- Changing the solver, placement, clue-numbering, or output format.

## 4. Decisions taken

- **D1. Mandatory input file.** The main CLI has no compiled-in default clue
  set. `--input File` is required; `.json` files use the interchange schema
  and `.pl` files are Prolog fixtures containing a `clues/1` term. This is an
  option rather than the first positional because direct SWI shebang execution
  can treat a leading `.pl` argument as another Prolog source file before the
  script receives argv.
- **D2. Single internal representation.** Both sources produce the existing
  `Words = [[Answer, MetaDict], …]` list (atom `Answer`, dict `MetaDict`)
  before anything downstream runs. The loader lives in `core.pl`; the CLI
  applies the shared uniqueness guard and passes the same representation into
  the selected arrange path.
- **D3. Answers normalised to atoms — for cross-source uniformity, not
  because the solver demands it.** `json_read_dict` yields SWI *strings*; the
  loader converts each answer to an atom (`atom_string/2`). To be precise about
  *why*: the solver runs fine on strings (`atom_chars/2` accepts a string), and
  atom-vs-string answers emit byte-identical JSON, so neither the solver nor
  output fidelity forces the conversion. What forces it is the `answer_meta/3`
  join and `check_unique_answers/1`, which compare answers with `==`:
  `atom == string` *fails*, so the answer stored in the placed word and the
  answer in the lookup list must be the *same type*. Normalising to atom
  matches the Prolog path (which supplies atoms), giving one type across both
  sources so the join, the `msort` in uniqueness-checking, and the output all
  behave identically. Byte-identical output across sources holds **given equal
  metadata content**: `json_read_dict` returns string meta *values* (vs. the
  Prolog path's atoms) and tags its dicts `#` (vs. an unbound tag), but
  `json_write_dict` renders both value types to the same JSON token and ignores
  the tag, so the difference is invisible in output.

## 5. Architecture

```
   fixture.pl (read clues/1 term) ─┐
                                   ├─► Words = [[Answer, MetaDict], …] ─► unchanged
   puzzle.json (json_read) ────────┘                                      pipeline
```

The `arrange` CLI decides the source from the input-file extension, calls the
shared loader, and passes the resulting `Words` into the existing flow.

## 6. JSON input schema

The input mirrors the output: an entry is an output `words[]` object minus the
solver-computed fields (`number`, `direction`, `cells`).

```json
{
  "clues": [
    { "answer": "OMEGA POINT",
      "meta": { "clue": "Transcending entropy",
                "link": "http://en.wikipedia.org/wiki/Omega_Point" } },
    { "answer": "FLOW" },
    { "answer": "BIAS", "meta": {} }
  ]
}
```

- **Top level** — an object with a `clues` array. A wrapper object (rather than
  a bare top-level array) leaves room for future keys such as `version` or a
  puzzle-default `gridLength` without another breaking change.
- **`answer`** — required string. Spaces allowed (multi-word); stripped before
  placement, preserved as the display form in output (unchanged behaviour).
- **`meta`** — optional object, default `{}`. Opaque passthrough, copied
  verbatim to the output entry's `meta`. By convention `clue` and `link`, but
  any keys are allowed; the solver never inspects it.

Nesting metadata under `meta` (rather than a flat `{answer, clue, link}`) is
**decided** (§11): it keeps exact symmetry with the output `meta` object and
reuses the collision-safety rationale from the output spec (§6.2 there). As
with the output, `meta` values are re-emitted via `json_write_dict`, so they
inherit that spec's contract (§9 there) — values must be JSON-friendly. A flat
alias is rejected: it would double the validation surface, and hand-authoring
ergonomics are what the bundled Prolog fixture is for, not the interchange
format.

## 7. Mapping to the internal form

| JSON | Internal |
|------|----------|
| `clues` array entry | one `[Answer, MetaDict]` pair |
| entry `answer` (string) | `Answer` (atom, via `atom_string/2`) |
| entry `meta` (object) | `MetaDict` (dict; `_{}` if absent) |
| (whole file) | `Words = [[Answer, MetaDict], …]` |

The result is indistinguishable from what `clues/1` returns, so
`check_unique_answers/1` (the existing answer-uniqueness guard) applies
uniformly.

## 8. CLI

The input file is mandatory on `arrange`:

    ./crosswordsmith arrange --input puzzle.json --size 17
    ./crosswordsmith arrange --input fixtures/bundled_17_clues.pl --max-size 17
    ./crosswordsmith arrange --input puzzle.json --size 17 --seed 7

`library(optparse)` handles the verb's options. `--input FILE` composes with
the documented arrange framing, fragment, candidate, enumeration, output, and
variety options subject to their normal compatibility rules; see the README
for the complete current matrix.

The clue loader dispatches by extension:

```prolog
load_clues(File, Words) :-
    file_name_extension(_, Ext0, File),
    downcase_atom(Ext0, Ext),
    load_clues_by_extension(Ext, File, Words).

load_clues_by_extension(json, File, Words) :- read_clues_json(File, Words).
load_clues_by_extension(pl, File, Words)   :- read_clues_prolog(File, Words).

read_clues_json(File, Words) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Doc), close(S)),
    doc_to_words(Doc, Words).
```

`read_clues_prolog/2` reads terms until it finds `clues(Words)`, rather than
consulting the file, so input fixtures do not define or redefine global
predicates. `doc_to_words/2` validates JSON documents (§9) and maps each entry
to `[AtomAnswer, MetaDict]`.

**`--out File` (added alongside).** Writes the output to a file instead of
stdout. Because `emit_json/3` already writes to `current_output`, this needed
no change to the solver or emitter: `run/2` wraps the dispatched goal in
`with_output/2`, which for a file captures the output to a string and writes
it only on success — so an unsolvable run leaves no empty file — and for the
default (`''`) calls the goal directly, keeping the stdout path byte-identical.

## 9. Validation

Two classes of error, handled differently:

- **From `json_read_dict` (free, and already readable).** A missing file or
  malformed JSON throws standard ISO errors (`existence_error`,
  `syntax_error`) that SWI's default handler renders clearly. No work needed.

- **Schema violations (must be checked *explicitly*).** `json_read_dict` does
  not validate shape, and `atom_string/2` silently *coerces* — `atom_string(X,
  42)` gives `'42'`, and a JSON `null` arrives as the string `"null"`. So a
  numeric or null answer is **not** rejected unless the loader guards types
  *before* converting. `doc_to_words/2` must explicitly check and fail on:
  - `clues` missing or not a list (`get_dict/3` merely *fails* on a missing
    key; a non-list value parses happily);
  - an entry whose `answer` is absent or not a string (guard with `string/1`
    before `atom_string/2`);
  - a `meta` that is present but not an object (`is_dict/1`).

**Error rendering.** Throw standard `error(Formal, _)` terms with a
domain-specific `Formal`, and add a `prolog:error_message//1` clause per
`Formal` so the message reads cleanly. A bare custom `throw/1` term (e.g.
`throw(bad_clues(...))`) is **not** clear: under this program's harness
(`set_prolog_flag(verbose, silent)` + `initialization(main, main)`) the
default handler prints `Unknown message: …`. That was the model originally
named here; `check_unique_answers/1` has since been converted to the
`error/2` + `error_message//1` pattern, and the loader should follow it.

Answer uniqueness is deliberately not duplicated in the loader — the shared
`check_unique_answers/1` guard covers both sources before search.

## 10. Testing & docs

- **Fixture.** A small `tests/<name>.json` clue set.
- **plunit.** Assert it loads to the expected internal `Words`; that an
  arrange pipeline over it emits valid JSON with the right shape; and that
  malformed inputs (no `clues`, non-string answer, invalid JSON) throw.
- **Symmetry.** A test that takes the emitted output, reduces each `words[]`
  entry to `{answer, meta}`, and confirms it round-trips as valid input.
  Compare at the **JSON-text** level (re-emit and diff) or normalise both
  sides first: emitted `answer`/`meta` values are JSON *strings* while authored
  Prolog fixtures use *atoms*, so a direct term comparison against `clues/1`
  would see a spurious `"FLOW" \== 'FLOW'` mismatch.
- **Docs.** README usage section; this spec; cross-link from the output spec.

## 11. Resolved questions

- **Nested `meta`, not flat.** Resolved nested (see §6): it's the only shape
  under which output→input symmetry (G5) holds, and it inherits the output's
  collision-safety. No flat alias.
- **Required `--input File`.** There is no `--clues` flag and no compiled-in
  default. The input file is part of every CLI form, and extension dispatch
  selects JSON (`.json`) or Prolog fixture (`.pl`) loading.
- **Config in the file: deferred, but coupled.** `gridLength` / `startLoc`
  stay CLI args for v1; the wrapper-object shape leaves room to add them later.
  When that happens, CLI overrides file and a `--shuffle` run ignores a file
  `startLoc` — and note the coupling: making `gridLength` optional-in-file
  means changing `main/0`'s "`grid_length` is a required positional" grammar,
  so the positional grammar and the file config must be designed together.
- **Golden coverage: parse-and-assert for the JSON-input path.** Output is
  deterministic (a byte-exact golden is viable), but the existing shared
  byte-exact golden already covers determinism; the JSON-input tests should
  assert structure (consistent with output spec §10) rather than add a second
  byte-exact golden that doubles maintenance for the same guarantee.
