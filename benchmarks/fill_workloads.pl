% benchmarks/fill_workloads.pl - the PRODUCT / performance-tracking manifest for `fill`.
%
% A fixed set of dictionary-fill instances covering grid sizes, dictionary
% sizes, and a seeded rung. This is the fill product bench's single source of
% truth (run_fill.pl reads only this file), deliberately separate from the arrange
% manifest (benchmarks/workloads.pl) so the two ratchets never cross.
%
% WHY A FIXED SET: fill search and dictionary-load costs (inferences) are
% deterministic and machine-INDEPENDENT - the SAME counts native or under WASM,
% only the wall-per-inference constant changes. benchmarks/fill_baseline.json
% records both metrics and `make bench-fill-check` reports each per-rung delta
% (regressions fail, improvements are wins).
%
% WHY THESE RUNGS: every row completes deterministically under the current MAC,
% dom/wdeg, and restart engine with ample headroom under the benchmark's 2e9
% budget. The set keeps short-word grid-size coverage, two byte-frozen dictionary
% scales, one seeded row, and the historical core/heavy selection stable.
% Budget-saturating, not_proven, and knife-edge combinations remain excluded.
%
% Each rung was chosen by its MEASURED search cost (a low-budget triage probe
% found the completing regime; each promoted rung was confirmed at full budget).
% The dictionaries are byte-frozen, so the counts are exactly reproducible under
% the pinned SWI-Prolog version.
%
% fill_workload(RungId, GridFile, DictFile, Seeds, Iterations, Warmup, Expected, Tier, Budget)
%   RungId     - unique rung key (baseline key; grids/dicts are reused across rungs)
%   GridFile   - stock-grid mask, path relative to repo root
%   DictFile   - word list, path relative to repo root (byte-frozen)
%   Seeds      - none | a fragment file path (the seeds path; a rung must still
%                complete). The Phase 0 fixture trap (a pinned answer was EXEMPT
%                from the search's no-duplicate rule, so seeding a word the
%                search would place anyway - e.g. AAH - made it re-place the
%                word and the CLI emit threw unique_key_pairs, exit 1) is
%                CLOSED: seed answers now pre-seed the search's Used set
%                (fill.pl seed_used/3 feeds fill_search_inc/5's
%                `\+ memberchk(Word, Used)` dedup; see
%                docs/plans/fill-seed-pin-crash-fix.md), so such a rung fills
%                with the next candidate instead of crashing. The current seeds
%                (fill_seed_11a: CYANO/TOMMY/READD, distinctive later-alphabet
%                words from the rung's own unseeded solution) are kept - the
%                rung's tree and output are unchanged.
%   Iterations - measured repetitions (the gated metrics search_inf/load_inf are
%                deterministic - 1 sample is exact; extra samples confirm + median the wall)
%   Warmup     - unmeasured warmup repetitions. The first in-process fill_attempt
%                JIT-indexes fill's clauses; warmup>=1 measures the WARM count.
%                Cold/warm delta is a FIXED 319 search inferences (MEASURED,
%                Phase 0: identical on sq04/g21/g09 - a one-time first-call JIT
%                cost, independent of tree size). warmup>=1 means both the
%                recorded baseline and `bench-fill-check` measure the WARM count,
%                so the gate compares warm-to-warm (zero delta). Relative size of
%                the cold delta is 0.06% on the smallest rung, <0.001% on the
%                largest; see the results doc.
%   Expected   - the asserted Outcome (all rungs here: filled) / exit (0)
%   Tier       - core | heavy
%       core : routine rows; run by default (make bench-fill / -check).
%       heavy: additional dictionary/grid coverage; opt in with --heavy. Still
%              deterministic and ratcheted.
%   Budget     - per-operation inference budget for the SEARCH layer
%                (fill_attempt/8), raised above fill's shipped 8e8 default so a
%                completing rung runs to TRUE completion (a deterministic,
%                ratchetable count). 2e9 is ample (top rung ~35M, ~57x headroom).
%
% Order matters: all core rungs precede the heavy rungs, so a --heavy run has
% already warmed the process before any warmup-lean rung is measured.

% -- core: fast completers, run by default ------------------------------------
fill_workload(sq04_full,     'fixtures/fill_grid_04a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 4x4 square, ~3.30M
fill_workload(g11_full,      'fixtures/fill_grid_11a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 11x11, ~3.46M
fill_workload(g11_full_seed, 'fixtures/fill_grid_11a.json', 'fixtures/dict/enable1.txt',   'fixtures/fill_seed_11a.json', 3, 1, filled, core, 2_000_000_000). % 11x11, 3 seed pins, ~3.41M
fill_workload(sq05_full,     'fixtures/fill_grid_05a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 5x5 square, ~3.38M
fill_workload(g17_full,      'fixtures/fill_grid_17a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 17x17, ~3.77M
fill_workload(g21_full,      'fixtures/fill_grid_21a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 21x21, 154 slots, ~7.54M
fill_workload(g13_full,      'fixtures/fill_grid_13a.json', 'fixtures/dict/enable1.txt',   none, 3, 1, filled, core, 2_000_000_000). % 13x13, ~3.61M

% -- heavy: additional coverage; opt in with --heavy --------------------------
fill_workload(sq04_50k,      'fixtures/fill_grid_04a.json', 'fixtures/dict/enable_50k.txt', none, 3, 1, filled, heavy, 2_000_000_000). % 4x4 x 50k, ~1.15M
fill_workload(g15_full,      'fixtures/fill_grid_15a.json', 'fixtures/dict/enable1.txt',   none, 2, 1, filled, heavy, 2_000_000_000). % 15x15, ~3.70M
fill_workload(g17_50k,       'fixtures/fill_grid_17a.json', 'fixtures/dict/enable_50k.txt', none, 2, 1, filled, heavy, 2_000_000_000). % 17x17 x 50k, ~9.76M
fill_workload(g09_full,      'fixtures/fill_grid_09a.json', 'fixtures/dict/enable1.txt',   none, 2, 1, filled, heavy, 2_000_000_000). % 9x9, ~3.69M
