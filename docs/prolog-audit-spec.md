# Spec: SWI-Prolog code audit (mission brief for the review workflow)

## Mission

Audit all Prolog source in this repo for three things, in priority order:

- **(a) Best-fit standard-library predicates/builtins** — are we using the right
  library predicate for the need, or hand-rolling something the stdlib already
  provides (better, clearer, safer)?
- **(b) Idiomatic, elegant SWI-Prolog** — all else equal on performance, is the
  code written the way an expert SWI user would write it (determinism, control
  constructs, data structures, naming, clause shape)?
- **(c) Gotchas / footguns** — are we falling afoul of SWI/Prolog traps
  (term-copying, dedup/stability, negation-as-failure, double_quotes, stream/
  resource handling, instantiation, unbounded search, etc.)?

This is a **read-only audit**: produce a prioritized, verified **findings
report** with concrete fixes. Do **not** edit source. Applying approved fixes is
a separate follow-up.

## Ground truth (mandatory)

The authority is the **local manual** at `docs/reference/swi-manual/` — it is
SWI-Prolog **10.0.2**, the exact version this repo runs against.

- Every "use X instead" or "this is the idiom" claim **must cite a specific
  manual file** (e.g. `swi-manual/lists.md`, `swi-manual/clpfd.md`,
  `swi-manual/apply.md`, `swi-manual/pairs.md`, `swi-manual/assoc.md`,
  `swi-manual/ordsets.md`), ideally naming the predicate/section. Grep the tree;
  see `swi-manual/INDEX.md` for the topic map.
- Reason about **SWI 10.0.2 semantics specifically** (e.g. `double_quotes` is
  `string` by default; `sort/4` dedups only with `@<`/`@>`; autoloading). Do not
  import ISO/SICStus/other-Prolog assumptions.
- A claim the manual does not support is **dropped or marked UNVERIFIED** — never
  asserted from model memory.

## Scope

- **Tier 1 (deep audit):** `crossword.pl` (1000 LOC — solver, strategies, CLI,
  JSON I/O, error hooks), `quality.pl` (309 — greedy quality engine).
- **Tier 2 (review):** `tests/crossword.plt`, `tests/run_tests.pl`,
  `benchmarks/run_matrix.pl`, `benchmarks/run_benchmarks.pl`,
  `benchmarks/start_sensitivity.pl`.
- **Tier 3 (light — bugs/footguns only, not style):** `fixtures/*.pl` (data) and
  any fixture generators.
- **Out of scope:** `docs/reference/swi-manual/**` (generated), `mercury/`
  (third-party), `scratch/` (user research).

## The three lenses — targeted checklists

Use these as starting points, not an exhaustive list. Examples reference real
patterns in this codebase.

### (a) Most appropriate predicates / builtins
- Hand-rolled recursions that a library predicate replaces (`library(lists)`,
  `library(apply)` — maplist/foldl/include/exclude/partition, `library(pairs)`,
  `library(aggregate)`).
- **Counting idiom:** `mrv_count/8` does `findall(t, capped(...), Ts),
  length(Ts, Count)` — is `aggregate_all(count, Goal, N)` (or a counting fold)
  the better fit, and does it compose with `limit/2` capping? Weigh against
  perf.
- **String/atom handling:** the recurring `atom_chars(Word, L), delete(L, ' ',
  L2), length(L2, N)` (strip spaces + measure) appears in several predicates —
  is `delete/3` the right tool, is there a cleaner SWI idiom (`exclude`,
  `split_string`, `atomic_list_concat`), and should it be factored into one
  helper?
- Right set/dict/assoc choices: `library(ordsets)` ops require sorted lists —
  confirm `ord_memberchk`/`ord_intersection` inputs are genuinely ordsets (e.g.
  `dir_cells/3` sorts first — verify every caller does). assoc vs pairs vs dict
  appropriateness.
- `library(solution_sequences)` (`distinct/2`, `limit/2`, `order_by/2`) where it
  would replace manual machinery.
- Randomisation: is `shuffle/2` a custom predicate where
  `random_permutation/2` exists?
- Option handling: `library(option)` for `memberchk(flag(X), Opts)` patterns in
  `run/2` / `build_floors/2`.

### (b) Idiomatic / elegant SWI-Prolog
- **Determinism & control:** appropriate `!` vs `->/;` vs `once/1` vs `\+`; green
  vs red cuts; predicates that should be deterministic leaving choicepoints (or
  vice-versa). NB the quality engine is **intentionally cut-free** (spec v1b.1) —
  respect that.
- Accumulator/tail-recursion vs `foldl`; `maplist` over manual recursion.
- `library(yall)` lambdas (`[X]>>Goal`) where they would clarify a maplist/foldl
  goal.
- Type checking: `must_be/2`, `is_of_type/2` at boundaries (CLI/JSON input).
- **Error/message idiom:** verify the `prolog:error_message//1` DCG hook used for
  custom errors is the correct, current SWI extension hook (vs
  `prolog:message//1`), citing the manual — and that `throw(error(Formal,
  Context))` uses standard formal terms (`type_error/3`, `domain_error/2`,
  `existence_error/2`).
- Module/script structure: `crossword.pl` is a non-module script loaded by tests;
  note (do not necessarily change) any cost of `user`-module predicates and
  clause grouping (`:- discontiguous` smells).
- First-argument indexing: are multi-clause predicates shaped to index well?
- `format/2,3` directive correctness; `with_output_to/2` usage.

### (c) Gotchas / footguns (SWI 10.0.2)
- **Term-copying:** `findall/3` copies bindings — already the cause of a fixed
  bug (`map_list_to_pairs` is used deliberately to keep original terms in
  `select_inc`/quality; **do not flag that as wrong**). Look for *other* places
  where a copy could bite (e.g. terms later compared with `==` or removed with an
  `==`-based op).
- **`sort/4` semantics:** `@>=` keeps duplicates and is stable; `@>` dedups.
  Confirm every `sort/2|4`/`keysort/2`/`predsort/3` uses the variant its caller
  actually needs (dup-keeping where ties matter; stability assumptions).
- **`double_quotes`:** default is `string` in SWI — any place a `"..."` literal
  is assumed to be codes/chars, or mixed with `atom_chars`/`atom_codes`?
- **Negation/instantiation:** `\+`/`once` over goals with unbound vars
  (floundering); `is/2` and arithmetic comparison on possibly-unbound terms.
- **Resource/stream handling:** `--out` file path (`with_output_to`, partial
  files on failure, closing/flushing), `format(user_error, ...)` for reports.
- **Unbounded search / memory:** the project has an OOM history; flag any
  `findall`/`aggregate_all` or non-tail recursion that can blow the stack on
  large inputs (cross-ref `docs/experiments.md` op-notes).
- **assoc:** `get_assoc/3` *fails* (not errors) on a missing key — confirm
  callers rely on the right behaviour.
- Benchmark/test harness: `call_with_time_limit/2`, `statistics/2`, any
  `nb_setval`/global state and its interaction with backtracking.

## Hard constraints (behaviour preservation)

- The **66 plunit tests + the golden output must stay green**, and the JSON
  **input/output contract** (`docs/json-input-spec.md`, `docs/json-output-spec.md`)
  must not change. Any finding whose fix *could* alter observable output is
  flagged **behaviour-risk: HIGH**.
- Solver **soundness/completeness** must be preserved (strategies only reorder
  the same search tree).
- **Performance is governed by `docs/experiments.md` + the benchmarks**, with
  **inferences** as the portable metric. "All else equal" means: a fix that may
  trade performance for elegance is **not auto-recommended** — it goes in a
  separate **NEEDS-BENCHMARK** bucket with an expected-impact note.

## Already decided — do NOT re-litigate (cite, don't re-open)

- `map_list_to_pairs` (not `findall`) in `select_inc`/quality is **intentional**
  (findall copies terms → breaks `==`-based `remove_x`). See code comments + E5.
- The quality engine is **deliberately cut-free** (spec v1b.1).
- `mrv_inc_deg` (degree tie-break) — tried, **rejected** (experiments.md I6);
  value ordering — **rejected** (I2); static length-order — **rejected** (E4).
- The I5 `no_word_merge` maximality rule is a real fix; the assoc-list grid is
  the chosen structure; `sort/4` with `@>=` is chosen deliberately where used.

If you believe one of these is genuinely wrong, you may raise it **once** as a
`CHALLENGE` finding with strong, manual-backed evidence — but the default is they
are settled.

## Finding schema (one record per finding)

```
id | lens(a/b/c) | file:line(s) | severity(high/med/low/nit) | category
title
current      : <short snippet>
issue        : <what is wrong/suboptimal, concretely, on SWI 10.0.2>
proposal     : <concrete minimal change; code if small>
manual_cite  : swi-manual/<file>.md  (predicate/section)
behaviour    : preserving | risk-HIGH | uncertain
perf         : none | likely-better | needs-benchmark | likely-worse
confidence   : high/med/low
steelman     : <why the current code might actually be fine / counter-argument>
```

## Severity rubric

- **High** — a real bug, or a soundness/termination/leak/footgun risk likely to
  bite.
- **Med** — a clearly better stdlib fit or idiom, low risk, no perf cost.
- **Low** — minor idiom/readability.
- **Nit** — cosmetic.

## Verification bar (adversarial — every finding before it lands)

1. The issue is **real on SWI 10.0.2** (checked against the manual; reproduced if
   feasible).
2. The fix is **correct and behaviour-preserving** (would keep tests + golden
   green; output contract intact).
3. It is **genuinely more appropriate/idiomatic per the manual** — not merely
   different.
4. A **skeptic** considered whether the current code is intentional (indexing,
   determinism, a documented design decision, or a settled-decision above).

Findings failing any check are **dropped or downgraded to a NOTE**. Prefer
perspective-diverse verification (correctness / idiom / perf-risk as distinct
checks) over repeating the same check.

## Recommended workflow shape (for the orchestrator)

- **Phase 0 — Map:** per file, inventory predicates, library imports, and
  recurring patterns; emit a worklist of `(lens × Tier-1 file)` cells plus a
  hotspot list (e.g. `mrv_count`, `select_inc`, the strip-spaces idiom, the
  error hook, `--out` handling).
- **Phase 1 — Find (fan out):** one agent per `(lens × Tier-1 file)`, plus a
  Tier-2 pass; each emits structured candidate findings, **each with a manual
  citation** (agents grep `docs/reference/swi-manual/`).
- **Phase 2 — Verify (fan out, adversarial):** each candidate independently
  checked against the bar above, ideally by a different agent and lens.
- **Phase 3 — Synthesize:** dedup, group by file & theme, sort by
  severity × confidence, drop unverified; emit the report.
- **Optional completeness critic:** which predicates/files/lenses went
  uncovered?

## Deliverable

A single report at **`docs/prolog-audit-findings.md`**, with:

- An **executive summary**: counts by severity, the **top quick wins** (Med+,
  behaviour-preserving, no perf cost), anything **HIGH**, and the
  **NEEDS-BENCHMARK** list.
- Findings **grouped by file**, each in the schema above.
- A **"checked and cleared"** list — patterns that were audited and found correct
  (so the next pass doesn't re-tread them), explicitly including the
  settled-decisions confirmed still sound.
