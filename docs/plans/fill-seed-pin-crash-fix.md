# Fix plan: `fill` seed-pin crash via duplicate answer at emit

**Status:** proposed (plan only — no source/test/golden/baseline changes yet)
**Component:** `prolog/crosswordsmith/fill.pl`, `prolog/crosswordsmith/core.pl`
**Contract touched:** AC-FILL / INV-3 (fill failure = "report unfillable slots + fail", never a raw internal error)
**Design-spec change:** none. §8.4's contract ("output is a canonical layout; exit is non-zero, naming the unfillable slot(s), when no complete fill exists") is already correct. This change makes the *implementation honor it* on a path that currently escapes with an uncaught throw.
**Origin:** latent product bug found (out of scope) during the fill perf campaign's Phase 0 fixture trap — `docs/experiments.md:1074-1081`, `benchmarks/fill_workloads.pl:40-47`.

---

## 1. Summary

A `crosswordsmith fill` whose seed pin's answer is **also reachable from the dictionary** can make the search re-place that same word in an unseeded slot. The two identical answers collide at emit time and throw an **uncaught** `domain_error(unique_key_pairs, ...)`, exiting the CLI with code 1 and a raw internal error message instead of the clean `fill:` diagnostic the failure contract requires.

There are **two independent triggers**:

- **A. Seed reused by the search** — a seed answer that the dictionary can also place in an unseeded slot. The search's no-duplicate set never includes seed answers, so it may pick the seed word again.
- **B. Two identical seed answers** — two seed pins with the same answer. `apply_seed` only checks *cell* clashes, not duplicate answers.

Trigger A is not merely a diagnostics bug: the current code **crashes on grids that have a valid solution** (see the control case in §3). The fix is therefore a *correctness* fix, not just a catch-and-report.

---

## 2. Root cause (precise, at file:line)

### 2.1 The throw site

The crash is **not** in `assign_clue_numbers/2` (the campaign note's attribution is loose). `assign_clue_numbers/2` (`core.pl:1186`) numbers by start cell via `keysort` + `group_pairs_by_key`, neither of which rejects duplicate keys — it succeeds fine on duplicate answers.

The real throw is at **`prolog/crosswordsmith/core.pl:1346`**, inside `answer_meta_assoc/2`:

```prolog
answer_meta_assoc(Words, Assoc) :-
    findall(A-Meta,
            ( member(Entry, Words), Entry = [A|_],
              ( Entry = [_, M] -> Meta = M ; Meta = _{} ) ),
            Pairs),
    list_to_assoc(Pairs, Assoc).        % <-- domain_error(unique_key_pairs, Pairs) on duplicate keys
```

`list_to_assoc/2` (SWI `library(assoc)`) throws `domain_error(unique_key_pairs, List)` when `List` has duplicate keys. When two placed words share an answer, `Pairs` has two entries with the same key. The `_{}` empty-dict values in the observed error term (`['COW'-_{},'COW'-_{}]`) are the `Meta = _{}` branch — decisive confirmation of this exact site.

**Call chain on the fill path:**
`fill_place_and_emit/7` (`fill.pl:569-582`) → `emit_fill/4` (`fill.pl:586-589`) → `emit_json/3` (`core.pl:1259`) → `build_words/4` (`core.pl:1331-1334`) → `answer_meta_assoc/2` (`core.pl:1341-1346`).

The docstring at `core.pl:1338-1340` already states the invariant this violates: *"answers are unique (`check_unique_answers/1`), so no key clash."* The fill path breaks that assumption because it never runs `check_unique_answers/1` (see §5, Fix C).

### 2.2 Why a duplicate reaches emit (trigger A)

The search's no-duplicate rule is at `fill.pl:369-370`:

```prolog
member(Word, Cands),
\+ memberchk(Word, Used),
```

But the search is started with an **empty** `Used` in both entry points:

- `fill.pl:491` (in `fill_attempt/8`): `fill_search(SearchSlots, DictByLen, Index, none, [])`
- `fill.pl:507` (in `fill_attempt_masked/9`): `fill_search(SearchSlots, DictByLen, Index, Masks, [])`

Seed pins are excluded from `SearchSlots` (`fill_prepare/5`, `fill.pl:555-560`: `exclude(seeded_slot(SeededKeys), Slots, SearchSlots)`), but their answers are **never seeded into `Used`**. So a searched slot can be filled with a word equal to a seed answer. `slots_to_layout/3` (`fill.pl:454-461`) then builds `InputWords = [[A1],[A2],...]` — one `[Answer]` per placed word, seeds included — and emit throws on the duplicate.

> Note: the campaign docs reference `fill.pl:196` as the dedup guard. That line has since moved — `fill.pl:196` is now `pair_ordset/2`; the live dedup is `fill.pl:369-370`. The described mechanism is accurate; only the line number drifted.

### 2.3 Why a duplicate reaches emit (trigger B)

`apply_seed/4` (`fill.pl:121-127`) only rejects a seed that matches no slot (`fill_seed_no_slot`) or clashes at a shared *cell* (`fill_seed_clash`). It does **not** reject two seeds with the same *answer*. Two identical seed answers land in `AllSlots`, both get emitted, and `answer_meta_assoc` throws. Fix A does **not** cover this (neither slot is searched) — it needs its own guard (Fix B).

---

## 3. Reproduction recipe

Minimal grid: mask `["...", "###", "..."]` → two independent 3-cell **across** slots (rows 0 and 2) and **no down slots** (a run must be length ≥ 2; the blocked middle row splits every column into two length-1 cells). Slots do not cross, so the search fills row 2 independently of the row-0 seed.

Create three files (any scratch dir):

`grid_dup.json`
```json
{ "name": "dup-3", "size": 3, "symmetry": "rot180", "mask": ["...", "###", "..."] }
```

`seed_dup.json`
```json
{ "gridLength": 3, "words": [ { "answer": "COW", "direction": "across", "cells": [[0,0],[0,1],[0,2]] } ] }
```

`dict_dup.txt`
```
COW
```

### Trigger A (seed reused by search)
```
$ ./crosswordsmith fill --grid grid_dup.json --seeds seed_dup.json --dict dict_dup.txt
ERROR: [Thread main] Domain error: `unique_key_pairs' expected, found `['COW'-_86{},'COW'-_102{}]'
$ echo $?
1
```

### Control (proves it fails a SOLVABLE grid) — this is the key evidence
Dict = `COW\nPIG`. A valid non-duplicate fill exists (seed `COW` in row 0, `PIG` in row 2), yet the current code **still crashes**, because the search picks `COW` first (alphabetical) and never tries `PIG`:
```
$ printf 'COW\nPIG\n' > dict_cowpig.txt
$ ./crosswordsmith fill --grid grid_dup.json --seeds seed_dup.json --dict dict_cowpig.txt
ERROR: [Thread main] Domain error: `unique_key_pairs' ...
$ echo $?  # 1
```

### Trigger B (two identical seed answers)
```json
// seed_dupseed.json
{ "gridLength": 3, "words": [
  { "answer": "COW", "direction": "across", "cells": [[0,0],[0,1],[0,2]] },
  { "answer": "COW", "direction": "across", "cells": [[2,0],[2,1],[2,2]] }
] }
```
```
$ ./crosswordsmith fill --grid grid_dup.json --seeds seed_dupseed.json --dict dict_dup.txt
ERROR: [Thread main] Domain error: `unique_key_pairs' ...   # exit 1
```

### Direct throw-site confirmation
```
?- crosswordsmith_core:answer_meta_assoc([['COW'],['COW']], _).
ERROR: Domain error: `unique_key_pairs' expected, found `['COW'-_{},'COW'-_{}]'
```

---

## 4. Correct contract (what should happen)

Existing infeasible/not_proven fills already behave correctly:
`fill_place_and_emit/7` (`fill.pl:580-581`) reports on stderr via `fill_report_failure/5` and then `fail`s; `main/0` maps a plain `fail` to exit 1 (`crosswordsmith:32-37`: `catch((cli(Argv) -> Code=0 ; Code=1), Error, (print_message(error,Error), Code=1)), halt(Code)`). Verified:

```
$ ./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict dict_with_no_fit.txt
fill: no complete fill exists for this grid + dictionary (+ seeds)   # stderr, no stdout
$ echo $?  # 1
```

Both the crash and the correct path exit 1. The difference is **diagnostic quality**:
- **Crash:** a raw `ERROR: [Thread main] Domain error: unique_key_pairs ...` rendered by `main/0`'s catch (an uncaught internal throw). No `fill:` framing; leaks an engine internal.
- **Correct:** a clean, `fill:`-prefixed message (INV-3: never silent, never a raw internal error). Seed *input* errors already do this via hooked messages (`fill_seed_no_slot` / `fill_seed_clash`, hooks at `fill.pl:782-785`), exit 1 — the pattern Fix B follows.

So the target behaviour for each trigger:
- **A, when a non-duplicate fill exists:** emit that valid fill (exit 0). *(correctness — the search must avoid the seed word)*
- **A, when no non-duplicate fill exists:** clean `infeasible` report + exit 1.
- **B:** clean, hooked seed-duplicate error + exit 1, reported before searching.

---

## 5. Proposed change

### Fix A — root cause (REQUIRED): seed the search's `Used` with the seed answers

Seed the initial `Used` set with the seed answers in **both** search entry points:
- `fill_attempt/8` — `fill.pl:491`
- `fill_attempt_masked/9` — `fill.pl:507`

The seeded slots are exactly `AllSlots \ SearchSlots`, and their `Vars` are **already** the bound word in the exact representation `Used` uses. Add a helper computed **before** `call_with_inference_limit/3`, so its cost is *not* charged to `search_inf`:

```prolog
% seed_used(+AllSlots, +SearchSlots, -Used0)
% The seeded slots (AllSlots minus the searched slots) are fully bound; their
% Vars ARE the placed word (a char list, the same shape Used holds), so the
% search's `\+ memberchk(Word, Used)` at fill.pl:369-370 dedups searched slots
% against the seed pins.
seed_used(AllSlots, SearchSlots, Used0) :-
    subtract(AllSlots, SearchSlots, SeededSlots),
    maplist(slot_word, SeededSlots, Used0).

slot_word(slot(_, _, _, Vars), Vars).
```

Then change both call sites to pass `Used0` instead of `[]`, e.g.:

```prolog
seed_used(AllSlots, SearchSlots, Used0),
call_with_inference_limit(
    ( fill_search(SearchSlots, DictByLen, Index, none, Used0) -> R = ok ; R = exhausted ),
    Budget, Limit),
```

`subtract/3` is **not** currently imported; add it to the `library(lists)` import list at `fill.pl:32`:
```prolog
:- use_module(library(lists), [append/3, member/2, nth0/3, select/3, subtract/3]).
```

**⚠ The char-list-vs-atom trap (validated the hard way).** `Used`, `Cands`, and dictionary words are **char lists** (e.g. `['C','O','W']`), *not* atoms. A first prototype seeded `Used` with `atom_chars`-collapsed atoms (`'COW'`) — `\+ memberchk(['C','O','W'], ['COW'])` never matches, so the dedup silently no-ops and the crash persists. Seed with the slot `Vars` **directly** (they are already the char list). Do **not** apply `atom_chars`.

**Prototype validation** (run against the live engine via internal predicates, no source edited):
- dict `[COW,PIG]`, seed `COW`: search skips `COW` (now in `Used`), places **PIG** → `Outcome = filled`, answers `[COW,PIG]` (unique) → valid JSON. *(the control grid from §3 now solves instead of crashing.)*
- dict `[COW]` only, seed `COW`: search **fails** → `fill_report_failure(infeasible, ...)` → `empty_slots/4` is `[]` (a slot's raw `candidate_count/4` ignores `Used`, so it is not "0 candidates") → clean generic message *"fill: no complete fill exists for this grid + dictionary (+ seeds)"* + exit 1.

**Why unseeded fills stay byte-identical:** when `SeedFile == none`, `SearchSlots == AllSlots`, so `subtract` yields `[]` → `Used0 = []` → `fill_search` is called exactly as today. Zero change to unseeded rungs/goldens.

### Fix B — seed-vs-seed duplicate (REQUIRED): reject in `apply_seeds`

Trigger B is not covered by Fix A. In `apply_seeds/3` (`fill.pl:117-119`), after the fragments are parsed, reject duplicate seed answers with a clean, hooked throw consistent with the existing seed errors:

```prolog
throw(error(fill_seed_duplicate(Answer), _))
```
Detect the duplicate the same way `check_unique_answers/1` does (sort the seed answers, `nextto(Dup, Dup, Sorted)`), and add an `error_message//1` clause beside the others at `fill.pl:782-785`:

```prolog
prolog:error_message(fill_seed_duplicate(A)) -->
    [ 'fill: seed ~q is pinned more than once; answers must be unique'-[A] ].
```

This is reported *before* searching, matching the documented seed-error contract (`fill.pl:471-472`: "a seed clash / no-slot throws"). Exit 1.

### Fix C — defense in depth (RECOMMENDED, nearly free): route fill's emit through the existing `check_unique_answers/1`

`check_unique_answers/1` **already exists, is exported, and is hooked** — it is not new code:
- exported at `core.pl:67`; definition `core.pl:1239`; docstring `core.pl:1231`
- clean hooked message `duplicate_answer(_)` at `core.pl:1250-1251`
- already called on the *solve* path in `crossword/4` at `core.pl:341`, immediately before `assign_clue_numbers`/`emit_json`.

The fill path bypasses `crossword/4` (it calls `emit_json/3` directly via `emit_fill/4`), so it never runs this check. Fix C = **route fill's emit boundary through the check that the solve path already uses.** Call `check_unique_answers/1` on `InputWords` in `slots_to_layout/3` (`fill.pl:454-461`), before emit. Any *future* duplicate (from an unforeseen path) then yields the clean, hooked `duplicate_answer` message + exit 1 instead of the raw `domain_error(unique_key_pairs)`. Cost: one predicate call over the already-materialized answer list; it reuses the existing message.

**Caveat (keep this explicit):** Fix C **cannot be the sole fix.** A naive catch-and-report at emit would wrongly fail *solvable* grids — in the §3 control (dict `[COW,PIG]`), `PIG` was available, so the correct outcome is a valid fill, not a failure. Fix A carries the correctness; Fix C is a safety net that turns any residual, genuinely-unsolvable duplicate into a clean message rather than an internal throw. Ship A + B; add C as insurance.

---

## 6. Regression tests

### 6.1 White-box plunit tests (`tests/fill.plt`, suite `fill`)

Follow the existing `do_fill` / `apply_seeds` / `fill_attempt` idioms (`tests/fill.plt:20-27, 88-97, 331`). Add small fixtures under `fixtures/`:
- `fill_grid_split3.json` — the `["...","###","..."]` grid (two independent across slots).
- `fill_seed_cow_top.json` — seed `COW` across row 0.
- `fill_seed_cow_both.json` — seeds `COW` across row 0 **and** row 2 (for trigger B).
- `dict_cow_pig.txt` — two lines `COW`, `PIG`.
- `dict_cow.txt` — one line `COW`.

Tests:

1. **`fill_seed_answer_reused_takes_alternative`** — split grid, seed `COW`, dict `[COW,PIG]`.
   Build `AllSlots`/`SearchSlots` via `fill_grid` + `apply_seeds` + `exclude(seeded_slot(...))`, `load_dict`, then `fill_attempt(SearchSlots, AllSlots, DBL, Idx, Outcome, Numbered, _)`.
   Assert `Outcome == filled` **and** the placed answers `sort` to `[COW,PIG]` (proves the searched slot took the non-duplicate `PIG` — this guards the *correctness* fix, not just no-throw).

2. **`fill_seed_answer_only_word_is_infeasible`** — same setup, dict `[COW]` only.
   Assert `Outcome == infeasible`. The body must **not** throw — plunit fails an uncaught error, so this is the direct crash-regression guard (it fails today with `domain_error(unique_key_pairs)`).

3. **`fill_duplicate_seed_answers_rejected`**, `[throws(error(fill_seed_duplicate('COW'), _))]` — apply two `COW` seed frags via `apply_seeds_frags/2` (`tests/fill.plt:331`) or `apply_seeds` on `fill_seed_cow_both.json`.

*(If Fix C is included, optionally also assert that a hand-built duplicate `InputWords` through `slots_to_layout/3` throws `error(duplicate_answer(_), _)` rather than `domain_error(unique_key_pairs, _)`.)*

### 6.2 CLI exit-code + clean-stderr test (REQUIRED — not optional)

A white-box test through `fill_attempt/8` cannot catch a regression in `main/0`'s error rendering: the contract under test is **CLI exit 1 + a clean `fill:` stderr line + no raw `Domain error`**. Because *both* the crash and the clean report exit 1, the stock `check_exit` helper (`run_tests.sh:45-55`, which discards stderr) is **insufficient on its own** — it would pass on the crash. Add a small helper that also inspects stderr, and use it.

Add near the other `check_*` helpers in `run_tests.sh`:
```sh
# check_fail_report <name> <command...>: assert a clean fill-failure (exit 1,
# a `fill:`-prefixed stderr line, and NO raw internal throw). Distinguishes the
# clean report from an uncaught engine error, which check_exit (code only) cannot.
check_fail_report() {
    local name="$1"; shift
    local err; err="$("$@" 2>&1 >/dev/null)"; local got=$?
    if [ "$got" -eq 1 ] \
       && printf '%s' "$err" | grep -q '^fill:' \
       && ! printf '%s' "$err" | grep -qi 'Domain error\|unique_key_pairs'; then
        echo "fail-report ($name): OK"
    else
        echo "fail-report ($name): FAILED (exit $got)"; printf '%s\n' "$err" | head -3; status=1
    fi
}
```
Then, in the CLI exit-code section (alongside the existing cases at `run_tests.sh:93-105`), using the committed fixtures:
```sh
# A seed answer the dictionary can also place must report cleanly, not throw a
# raw unique_key_pairs domain error (regression: fill-seed-pin-crash-fix).
check_fail_report "fill seed-reused unsolvable -> clean fail" \
    ./crosswordsmith fill --grid fixtures/fill_grid_split3.json \
        --seeds fixtures/fill_seed_cow_top.json --dict fixtures/dict_cow.txt
# Two identical seed answers must be rejected before searching (Fix B).
check_fail_report "fill duplicate seeds -> clean fail" \
    ./crosswordsmith fill --grid fixtures/fill_grid_split3.json \
        --seeds fixtures/fill_seed_cow_both.json --dict fixtures/dict_cow.txt
```
*(For trigger B, the stderr line is `fill: seed 'COW' is pinned more than once; ...` — also `fill:`-prefixed, so the same helper matches.)*

No new **golden** stdout files are needed: the failing cases emit nothing on stdout; the solvable case (test #1) is exercised at the predicate level.

---

## 7. Benchmark baseline re-record

### 7.1 Why it moves, and by how much

Fix A lengthens the hot `\+ memberchk(Word, Used)` (`fill.pl:369-370`) by the seed count. This affects **only the one seeded rung, `g11_full_seed`** (`benchmarks/fill_workloads.pl:74`; seeds `CYANO/TOMMY/READD`). Its recorded `search_inf` baseline is **`409191`** (`benchmarks/fill_baseline.json:55`), gated at **±0.5% relative** (`benchmarks/check_fill_baseline.pl:11-12`; classify at `:151-153` — a rise past tolerance is a hard FAIL when the SWI version matches).

The **filled grid is unchanged** (the search never wanted those seed words in the searched slots), so:
- The **identity oracle** (`benchmarks/fill_identity.sha256`) stays green — output is byte-identical; **do not** re-record it.
- Only `search_inf` ticks up, by ≈ `seed_count × (number of `\+ memberchk` scans)` = a small constant-factor rise (3 extra list cells per scan). This will very likely exceed the ~2046-inference (0.5% of 409191) tolerance, so the ratchet will FAIL until re-recorded.
- All **unseeded** rungs are byte-identical in `search_inf` (their `Used0 = []`).

### 7.2 Procedure

1. Land Fix A/B (and C).
2. Run the fill ratchet in CHECK mode; confirm the **only** `search_inf` movement is `g11_full_seed` and it is a **rise** of the expected small magnitude:
   ```
   swipl -q benchmarks/check_fill_baseline.pl        # diff + PASS/FAIL
   ```
3. **Sanity bound (do not skip):** the `g11_full_seed` rise should be ≈ `seed_count(3) × per-node scans` — a modest constant factor. If the rise is *substantially* larger than that (e.g. a multiple, or any other rung moves), Fix A introduced an **unintended traversal** (e.g. `subtract`/`seed_used` accidentally inside the inference-limited goal, or a wrong `Used` shape) — stop and investigate before recording.
4. Re-record and log:
   ```
   swipl -q benchmarks/check_fill_baseline.pl --record   # ratchet fill_baseline.json (+ read the diff it prints)
   swipl -q benchmarks/check_fill_baseline.pl --log      # append fill_history.jsonl
   ```
   Label the `fill_history.jsonl` entry (and the commit message) **"correctness-driven"** — this is a bug fix that changes seed-search dedup semantics, *not* a performance regression. Note the g11_full_seed delta explicitly.

---

## 8. Documentation updates

- `docs/experiments.md:1074-1081` — mark the "Fixture trap found (latent product bug)" note **resolved**, pointing at this plan / the fixing commit.
- `benchmarks/fill_workloads.pl:40-47` — the `Seeds` CAUTION comment describes the now-fixed behaviour and cites the stale `fill.pl:196`; update it to note the trap is closed (seed answers now participate in the dedup) and correct the line reference to `fill.pl:369-370`.
- `prolog/crosswordsmith/fill.pl:471-472` — extend the entry-point contract docstring to mention `fill_seed_duplicate` alongside `fill_seed_clash`/`fill_seed_no_slot`.
- If Fix C lands: note in `core.pl:1338-1340` (or the fill emit path) that the fill path now also runs `check_unique_answers/1`, so the "answers are unique" invariant holds on both the solve and fill paths.
- **No design-spec change.** §8.4's contract is already correct; this change only makes the implementation honor it.

---

## 9. Risks and edge cases

- **Char-list vs atom** (the validated trap): pre-seed `Used` with slot `Vars`, never `atom_chars` atoms — otherwise the dedup silently no-ops and the crash survives the "fix." The white-box test #1 (asserting `PIG` is chosen) catches this; a mere no-throw test would not.
- **`subtract/3` term identity:** correctness relies on `SearchSlots` being the *same* slot terms as in `AllSlots` (true — `exclude/3` at `fill.pl:560` preserves terms, and `subtract/3` is `\+ memberchk`-based). Slot counts are small (≤154) and this runs once, outside the inference limit — negligible cost, uncounted in `search_inf`.
- **Seeded slots fully bound:** a seed frag spans exactly its slot (`apply_seed` matches cells + direction, and answer length = slot length), so `Vars` is ground at this point — safe to read as the word.
- **Bench seam unchanged:** `fill_attempt/8`'s signature and arity are untouched, so the in-process bench call site (`benchmarks/fill_subjects.pl:77`) needs no edit; only the one seeded rung's number moves.
- **Fix C is not sufficient alone** (would mis-report solvable grids) — must ship with Fix A. Do not let a reviewer "simplify" the change down to just the emit guard.
- **Case/normalization of seeds:** the dedup compares in the search's own char-list representation, which is self-consistent with how seeds are pinned into cells; no new normalization is introduced or required here.

---

## 10. Effort estimate

**Small–medium, ~1–2 hours** including the bench re-record and doc touch-ups.

- Fix A: ~8 lines (2 call sites + `seed_used/3` + `slot_word/2` + 1 import). 
- Fix B: ~6 lines (dup check in `apply_seeds` + one `error_message//1` clause).
- Fix C: ~1-2 lines (route `slots_to_layout/3` through `check_unique_answers/1`).
- Tests: 3 plunit tests + 1 `run_tests.sh` helper + 2 `check_fail_report` lines + 5 tiny fixtures (~60 lines total).
- Baseline: run ratchet, sanity-check the `g11_full_seed` delta, `--record` + `--log`.

Low risk; changes are localized to `fill.pl` (+ one line in `core.pl`/`slots_to_layout` for Fix C), fully covered by the new white-box + CLI regression tests, with the unchanged identity oracle guarding output byte-identity.

---

## Appendix — file/line index (as of this plan)

| Item | Location |
|---|---|
| Throw site (`list_to_assoc` on dup keys) | `core.pl:1346` in `answer_meta_assoc/2` (`core.pl:1341-1346`) |
| Emit chain | `emit_json/3` `core.pl:1259` → `build_words/4` `core.pl:1331-1334` |
| `assign_clue_numbers/2` (NOT the culprit) | `core.pl:1186` |
| Fill emit boundary | `slots_to_layout/3` `fill.pl:454-461`; `emit_fill/4` `fill.pl:586-589` |
| Empty `Used` (root cause A) | `fill.pl:491` (`fill_attempt/8`), `fill.pl:507` (`fill_attempt_masked/9`) |
| Search dedup guard | `fill.pl:369-370` |
| Seed application | `apply_seeds/3` `fill.pl:117-119`; `apply_seed/4` `fill.pl:121-127` |
| SearchSlots exclusion | `fill_prepare/5` `fill.pl:555-560` |
| Failure report + fail | `fill_place_and_emit/7` `fill.pl:569-582`; `fill_report_failure/5` `fill.pl:591-603` |
| CLI exit mapping | `crosswordsmith:32-37` (`main/0`) |
| Existing seed error hooks | `fill.pl:782-785` |
| `check_unique_answers/1` (Fix C) | export `core.pl:67`; def `core.pl:1239`; hook `core.pl:1250-1251`; solve-path call `core.pl:341` |
| `lists` import (add `subtract/3`) | `fill.pl:32` |
| Seeded bench rung + baseline | `benchmarks/fill_workloads.pl:74`; `benchmarks/fill_baseline.json:44-60` (`search_inf:409191`) |
| Ratchet gate (±0.5% rise = FAIL) | `benchmarks/check_fill_baseline.pl:11-12, 151-153` |
| Test helpers | `tests/fill.plt:20-27` (`do_fill`), `:331` (`apply_seeds_frags`) |
| CLI test helpers | `run_tests.sh:28-55` (`check_golden`/`check_exit`), cases `:93-105` |
