% benchmarks/workloads.pl - the PRODUCT benchmark manifest.
%
% One row per thing we ship `arrange` to do. This is the product bench's single
% source of truth (run_arrange.pl reads only this file). It is deliberately
% separate from benchmarks/fixtures.pl (the strategy-research matrix): the two
% benches answer different questions and must not share a manifest, so a research
% knob can never leak into the product numbers (plan §6). No strategy column
% exists here - arrange always runs its production mrv_inc 4-corner search.
%
% arrange_workload(Fixture, Size, Mode, Iterations, Warmup, Expected, Tier)
%   Fixture    - path relative to the repo root; read by BOTH layers (arrange
%                --input accepts .pl, and the search layer reads its clues/1)
%   Size       - the N in --size N / --max-size N (both frame an N x N canvas)
%   Mode       - size | max_size (which flag the command layer passes)
%   Iterations - measured repetitions (1 for the slow budget-saturated probes)
%   Warmup     - unmeasured warmup repetitions (>= 1 for multi-iteration rows so
%                the first cold run never lands in the summary; plan m5)
%   Expected   - placed | infeasible; asserted every iteration, both layers
%   Tier       - core | heavy
%       core  : sub-~1.1 s realistic inputs; run by default, safe for `make bench`.
%       heavy : budget-SATURATING latency probes (~26 s/call). Excluded by default
%               (they'd dominate wall-clock); opt in with --heavy, or target one
%               with --fixture. Single-shot, latency-only, NOT regression-gated.
%
% The tier reflects a MEASURED product finding, not input size: arrange latency is
% bimodal. `arrange` tries four start corners under a SHARED ~500M-inference
% budget; when a non-placing corner's search fails fast the command is ~0.1-1 s,
% but when one triggers deep search it burns the whole budget and the command
% takes ~26 s - even though a valid placement was already found from another
% corner. benchmark_14 (14 words!) and benchmark_16_dense both hit the cliff; the
% single-corner research matrix hides it entirely (it reports benchmark_14 at
% 277k inf / 12 ms vs the command's 500M inf / 26 s). See plan §6/§7.
%
% Only fixtures that are realistic arrange inputs at a realistic grid are product
% workloads. benchmark_20 (37x37) and benchmark_26 (49x49) are intentionally
% EXCLUDED: they place only on very sparse oversized grids, so they exercise the
% research matrix, not a shape a user would actually arrange (plan §6).

% -- core: run by default -----------------------------------------------------
arrange_workload('fixtures/bundled_17_clues.pl',         17, size, 20, 3, placed, core).  % ~90 ms, wrapper-bound (~42% search)
arrange_workload('fixtures/benchmark_08_words.pl',       13, size,  8, 2, placed, core).  % ~470 ms
arrange_workload('fixtures/benchmark_70_mesh_words.pl',  21, size,  3, 1, placed, core).  % ~1.05 s, search-dominated (~94%)

% -- heavy: budget-saturating latency probes (opt in with --heavy) ------------
arrange_workload('fixtures/benchmark_14_words.pl',       17, size,  1, 0, placed, heavy). % ~26 s (500M-inf cliff)
arrange_workload('fixtures/benchmark_16_dense_words.pl', 17, size,  1, 0, placed, heavy). % ~26 s (500M-inf cliff)
