# Spec: JSON clue-input file

Status: **proposed** — not yet implemented; audited and revised.

Companion to [`json-output-spec.md`](json-output-spec.md), which is implemented.

## 1. Problem statement

The crossword's input lives in `clues.pl`, a Prolog source file pulled in at
compile time via `:- include('clues.pl').`. It defines `clues/1` as a list of
`[Answer, Metadata]` pairs (`Metadata` an SWI dict). `main/0` calls `clues/1`
to obtain the word list.

This means the *only* way to supply a puzzle is to author Prolog source and
have it compiled into the script. The output is already JSON
(`json-output-spec.md`), and the program was built to feed a web page, so an
external consumer that can *read* the solution cannot *produce* a puzzle
without writing Prolog. We want to accept clues from a JSON file so puzzles
can be generated programmatically and the input/output loop is symmetric.

## 2. Goals

- **G1.** Accept a clue set from an external **JSON file** at runtime.
- **G2.** Keep the bundled Prolog `clues/1` working unchanged as the default,
  so existing usage and the in-repo dataset need no migration.
- **G3.** Add **no new dependency** — reuse `library(http/json)`, already
  loaded for output.
- **G4.** Leave the solver and emitter untouched: both input sources converge
  on the existing internal representation before the pipeline runs.
- **G5.** Make the input schema the natural mirror of the output, so a
  solution payload can be reduced back to a valid input.

## 3. Non-goals

- Moving `grid_length` / `start_loc` into the file (they stay CLI arguments;
  see §11).
- Supporting arbitrary *external Prolog* clue files (see §4, decision D1).
- Changing the solver, placement, clue-numbering, or output format.

## 4. Decisions taken

- **D1. Two formats, distinct roles — not interchangeable peers.** The bundled
  Prolog `clues/1` remains the **zero-argument default** and the native
  authoring format for the in-repo dataset. JSON is the **external /
  interchange** format, selected explicitly. We are *not* replacing `clues.pl`
  with JSON (that forces a second input-format migration right after the
  output redesign and loses the compile-time `include`), and we are *not*
  supporting arbitrary external `.pl` clue files. The latter is a YAGNI call,
  not a technical one: editing the bundled `clues.pl` already covers Prolog
  authoring, so a second Prolog source path earns nothing. (It is also
  slightly awkward — `clues/1` is already defined via `include`, so a naive
  runtime `load_files/2` with a `module` option can clobber it; a `read_term/3`
  scan would sidestep that, but it isn't worth building.)
- **D2. Single internal representation.** Both sources produce the existing
  `Words = [[Answer, MetaDict], …]` list (atom `Answer`, dict `MetaDict`)
  before anything downstream runs. The change is a *loader* in `main/0`; the
  pipeline (`check_unique_answers/1` → `find_crossword/5` → `emit_json/3`) is
  unchanged.
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
   clues.pl  (include, default) ─┐
                                 ├─► Words = [[Answer, MetaDict], …] ─► unchanged
   --clues f.json (json_read) ───┘                                      pipeline
```

`main/0` decides the source from the CLI, calls a loader, and passes the
resulting `Words` into the existing flow. Everything after the loader is the
code that ships today.

## 6. Proposed JSON input schema

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
ergonomics are what the bundled `clues.pl` is for, not the interchange format.

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

## 8. CLI changes

A new flag selects the source:

    ./crossword.pl --clues puzzle.json 17 topleft_across
    ./crossword.pl --clues puzzle.json --shuffle 17
    ./crossword.pl 17 topleft_across          # unchanged: bundled clues

`--clues File` is a flag-with-value that must compose with the existing
`--shuffle` / `--all` flags and the positional `grid_length` / `start_loc`.
Today `main/0` matches `Argv` positionally, so the flag is stripped first and
the *remaining* args drive the existing patterns:

1. **Strip the flag, then match the rest (minimal).** Remove `--clues File`
   from `Argv` and run the existing `main/0` clauses on what's left.
   **Caution:** the extraction must preserve flag→value adjacency. A
   `select('--clues', Argv, R), select(File, R, _)` pair is *wrong* — the
   second `select/3` binds `File` to the first remaining token, which is the
   flag's value only when `--clues` happens to be first (it mis-binds
   `--shuffle` for `[--shuffle, 17, --clues, p.json]`). Use an adjacency-
   preserving extractor instead.
2. **`library(optparse)`.** With `--shuffle`, `--all`, `--clues`, and
   positionals, the flag set is at the point where a real parser pays its way;
   prefer this if touching the parsing at all.

The loader, sketched (extension sniffing dropped — JSON is the only external
format, so just parse and let JSON errors speak):

```prolog
% Strip --clues from Argv; RestArgv feeds the existing positional matching.
load_clues(Argv, RestArgv, Words) :-
    ( strip_clues_flag(Argv, File, RestArgv)
    -> read_clues_json(File, Words)
    ;  RestArgv = Argv, clues(Words)        % bundled default
    ).

% Adjacency-preserving: only the token *after* --clues is the file.
strip_clues_flag(['--clues', File | Rest], File, Rest) :- !.
strip_clues_flag([A | T], File, [A | Rest]) :- strip_clues_flag(T, File, Rest).
% (no clause for []: absence of --clues makes load_clues fall to the default)

read_clues_json(File, Words) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Doc), close(S)),
    doc_to_words(Doc, Words).
```

`doc_to_words/2` validates the document (§9) and maps each entry to
`[AtomAnswer, MetaDict]`.

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

Answer uniqueness is *not* re-checked in the loader — the existing
`check_unique_answers/1` in `crossword/3` already covers it for both sources.

## 10. Testing & docs

- **Fixture.** A small `tests/<name>.json` clue set.
- **plunit.** Assert it loads to the expected internal `Words`; that
  `crossword/3` over it emits valid JSON with the right shape; and that
  malformed inputs (no `clues`, non-string answer, invalid JSON) throw.
- **Symmetry.** A test that takes the emitted output, reduces each `words[]`
  entry to `{answer, meta}`, and confirms it round-trips as valid input.
  Compare at the **JSON-text** level (re-emit and diff) or normalise both
  sides first: emitted `answer`/`meta` values are JSON *strings* while authored
  `clues.pl` uses *atoms*, so a direct term comparison against `clues/1` would
  see a spurious `"FLOW" \== 'FLOW'` mismatch.
- **Docs.** README usage section; this spec; cross-link from the output spec.

## 11. Resolved questions

- **Nested `meta`, not flat.** Resolved nested (see §6): it's the only shape
  under which output→input symmetry (G5) holds, and it inherits the output's
  collision-safety. No flat alias.
- **Flag name `--clues <file>`.** Matches the predicate `clues/1`; `--input`
  is too generic, `--clues-file` redundant. Extension sniffing is dropped
  (§8) — JSON is the only external format.
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
