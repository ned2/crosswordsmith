# Implementation plan: §8.4b fill search levers (`--budget`, `--seed`/`--shuffle`)

**Status: EXECUTING (2026-07-15).** Builds the LOCKED §8.4b contract
([design-spec](../design-spec.md), DP-6) — AC-FILL-9..11. Scope is **control
and variety only**: DP-6 measured both levers as non-completion-fixes (×20
budget ladder and 0/8 ordering probe on `blocked_13a`,
[`benchmarks/fill_quality/`](../../benchmarks/fill_quality/README.md)), and
the docs must say so. **No search-power change of any kind** — crossing-aware
forward-checking is the separate DP-7-class pass (§8.5 row).

## Architecture — four small seams, no new machinery

1. **`--budget N` (AC-FILL-9).** `fill_solve/5` gains a `budget(N)` option
   (default: the existing `fill_budget/1` fact — it stays the single source
   of truth). `fill_place_and_emit` gains the resolved budget argument and
   calls the budget-explicit `fill_attempt/8` (already the tested/gated seam;
   its arity and body are untouched) and `fill_attempt_masked/9` as today.
   The default-budget delegate `fill_attempt/7` stays (tests call it).
2. **PRNG: reuse `core.pl`'s engine-internal splitmix64 wholesale.**
   `set_search_seed/1`, `set_shuffle_seed/1`, `current_search_seed/1` are
   already exported; `seeded_permutation/2` is white-box-callable. Nothing
   new is written — the wasm-parity guarantee (never SWI's native RNG) is
   inherited (see wasm/README PRNG section).
3. **Candidate-order seam — at LOAD, not in the search (AC-FILL-10).**
   Bucket order *is* candidate order (`keysort/2` is stable), so the
   equal-score tie shuffle lands where the order is born:
   - scored path: after `order_score_desc/2`'s msort, iff `search_seed(_)`,
     `seeded_permutation` **within each equal-score group** (score-descending
     stays primary);
   - plain path: after `sort(Ws0, Words)`, iff seeded, `seeded_permutation`
     of the whole list (uniform dict ⇒ the tie is the whole length bucket).
   Unseeded cost: one failed `search_seed(_)` lookup per load — well inside
   the ratchet's 0.5% (`load_inf` already carries +30 from the §8.4a sniff).
4. **Slot-tie seam — a seeded twin of the search loop (AC-FILL-10).** The
   F-H2 idiom (`fill_attempt_masked` precedent): a `fill_search_inc_seeded`
   twin whose slot selection collects the equal-min-`Count` slots and picks
   one via `prng_draw` (Start/Dir order today = the tie DP-6's contract lets
   the seed vary). Dispatch **once** at `fill_search` entry on a
   `search_seed(_)` check — zero added goals per node on the deterministic
   path, `fill_attempt/8` byte-identical for the ratchet and the
   `fill_identity.sha256` oracle.

## CLI (mirrors arrange §7.6 verbatim)

`fill_opts_spec` gains `seed` (integer, default −1 sentinel), `shuffle`
(boolean), `budget` (integer, default −1 sentinel → engine default).
Reuse/mirror `validate_seed`, `seed_mode`, `apply_seed_mode`. Guards:
`--seed`/`--shuffle` mutually exclusive; either is rejected with
`--index`/`--save-index` (artifacts are order-pinned, score-free caches —
AC-FILL-11); `--budget` must be ≥ 1 and composes with every mode.
`--shuffle` prints `fill: seed N` under `--verbose` (fill's payload carries
**no** diagnostics — the stderr-contract test pins that — so provenance is
verbose-only, unlike arrange's `diagnostics.arrange.seed`).

## Tests / goldens / fuzz (the FS-3(c)-style CI subset)

- plunit: budget determinism + tiny-budget `not_proven` + budget-never-
  changes-content (AC-FILL-9); per-seed reproducibility, seeded-vs-default
  bucket-order divergence (white-box), no-flag identity (AC-FILL-10);
  rejection matrix (AC-FILL-11).
- Goldens: one seeded fill golden (bundled fixture + `--seed`), checked in
  `run_tests.sh` + `make golden`. Existing goldens prove the no-flag path.
- Fuzz: `--seed N` 3× byte-identity; small `--budget` failure determinism;
  `--shuffle` excluded from identity checks by construction (fresh seed per
  run) — cover its seed-report line instead.
- Ratchet: `make bench-fill-check` (and `--heavy` once, pre-final-commit)
  must stay green; identity oracle untouched.

## Commits

- C1 `feat(fill): --budget N` (engine option + CLI + plunit).
- C2 `feat(fill): fill --seed/--shuffle` (load seam + seeded twin + CLI +
  plunit).
- C3 `test(fill)/docs`: seeded golden, fuzz cases, README flag rows, §5-
  usage "build pending" marker dropped in design-spec, STATUS sweep.
- Each commit: `make unit` while iterating, `make test` before landing;
  `make fuzz` after C2/C3.

## Out of scope

Forward-checking (DP-7); wasm wire exposure of the new flags (browser fill
is untouched; the PRNG is portable by construction so nothing blocks a later
exposure); internal restart scheduling (a user-level loop over `--seed` ×
`--budget` composes the same envelope — documented in the README row).
