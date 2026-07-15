# fill capacity clean-fail — implementation plan (DP-10)

Build plan for design-spec §10 **DP-10** (the §8.5 capacity-robustness row,
resolved 2026-07-16 as the clean-error slice). Failure-path only: no search,
load, default, or output change on any success path — goldens, identity, and
ratchets untouched by construction.

## The two modes (measured, DP-9 grounding)

- **Load/index**: >~500k-word dict → `resource_error(stack)` inside
  `load_dict`'s index build (reproducer: `scripts/fetch-cwl.sh --min-score 1
  --out /tmp/cwl_full.dict`, then any `fill --dict /tmp/cwl_full.dict`).
- **Search**: ~437k-word band on a full-13 grid → overflow in `mac_support`
  bignum masks ~17s in (reproducer: the same dict at `--min-score 30` on
  `grids/blocked_13a.json`).

Both currently exit 1 with SWI's raw multi-line stack dump (bignum spew on
the search mode). AC-FILL-15 replaces the dump with the house failure line.

## Deliverables

1. **`fill.pl` guard** — one helper (`fill_capacity_guard/2`-shaped: catch
   `error(resource_error(stack), _)`, emit the report, fail; everything
   else rethrows) wrapped around the two phases at the entry seams:
   - load/index: the `load_dict` call in `fill_solve/5`, `fill_load_index`
     in `fill_solve_index/6`, and `fill_save_index`'s dict load;
   - search+emit: the `fill_place_and_emit` calls on both paths.
   Report shape (one line + one remedy line, budget-exhaustion style):
   - load: `fill: capacity exceeded loading/indexing DICT (stack limit
     ~1.0Gb)` + `hint: shrink the list (e.g. --min-score 50), or relaunch
     as: swipl --stack-limit=4g ./crosswordsmith ...`
   - search: `fill: capacity exceeded during search on NxN grid (stack
     limit ~1.0Gb)` + the same hint with `raise --min-score` first.
   Limit read from `current_prolog_flag(stack_limit, _)` at report time.
2. **Tests (`tests/fill.plt`)** — (a) white-box mapping: the guard turns a
   thrown `resource_error(stack)` into the stderr report + plain failure
   (via `with_stderr_string`), and rethrows a non-capacity error unchanged;
   (b) a real-overflow test under a lowered per-thread `stack_limit`
   (setup/cleanup restores the flag) — keep only if stable across runs,
   else (a) plus the documented CLI reproducers carry the coverage.
3. **Docs** — README limitation bullet (crash → clean failure; envelope
   unchanged), `dicts/README.md` fetch-floor caveat wording,
   `benchmarks/fill_quality/README.md` finding 2 post-note (behaviour since
   DP-10; envelope numbers still current), `docs/STATUS.md` backlog note.

## Non-goals (explicit, per DP-10)

- No default stack-limit change, no streamed/segmented index build, no new
  CLI flag. The envelope moves nowhere; it just stops looking like a bug.

## Verification

`make test` (incl. new plunit; goldens/identity/ratchets diff-clean),
`make fuzz`; manual run of both reproducers above confirming the one-line
report + exit 1 + no partial `--out` file (AC-FILL-15).
