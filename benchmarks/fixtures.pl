% benchmarks/fixtures.pl - machine-readable benchmark manifest.
%
% Single source of truth for each fixture's intended solve configuration. The
% human-facing table in fixtures/README.md must stay in sync with this; the
% matrix harness (benchmarks/run_matrix.pl) reads ONLY this file, so a result
% batch can never silently run a fixture on the wrong grid.
%
% bench_fixture(RelPath, Grid, StartLoc, Iterations, Warmup)
%   RelPath    - fixture path relative to the repo root
%   Grid       - grid side length the fixture is designed for
%   StartLoc   - first-word start location
%   Iterations - measured repetitions (1 for the very slow dense fixture)
%   Warmup     - unmeasured warmup repetitions

bench_fixture('fixtures/bundled_17_clues.pl',        17, topleft_across, 30, 3).
bench_fixture('fixtures/benchmark_08_words.pl',      13, topleft_across, 30, 3).
bench_fixture('fixtures/benchmark_14_words.pl',      17, topleft_across, 30, 3).
bench_fixture('fixtures/benchmark_16_dense_words.pl', 17, topleft_across,  1, 0).
bench_fixture('fixtures/benchmark_20_words.pl',      37, topleft_across, 30, 3).
bench_fixture('fixtures/benchmark_26_words.pl',      49, topleft_across, 30, 3).
% Hard mesh (idea I4): hard for baseline (does not finish in the matrix's
% per-cell limit), tractable for the MRV strategies. 1 iteration / 0 warmup
% because the baseline cell is slow. See docs/experiments.md (I4).
bench_fixture('fixtures/benchmark_70_mesh_words.pl', 21, topleft_across,  1, 0).
