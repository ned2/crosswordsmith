% benchmarks/greedy_workloads.pl - permanent greedy arrange benchmark manifest.
%
% greedy_workload(Id, Fixture, Size, Framing, Command, SeedAnswer, Corner,
%                 ExpectedConstruction, Iterations, Warmup, Tier).
%
% Framing is size|max_size. Command is candidates(strict,K)|best_effort.
% ExpectedConstruction reifies the pinned product construction's result as
% completed(PlacedCount,DroppedCount); setup failure is represented by
% setup_failed rather than by sampler failure.

greedy_workload(
    bundled_17_candidates,
    'fixtures/bundled_17_clues.pl', 17, size, candidates(strict, 3),
    'NARRATIVE FALLACY', topleft_across, completed(6, 0),
    20, 3, core).

% The three longest bundled seeds do not fit on an 11x11 canvas. OMEGA POINT is
% deliberately pinned here because it is the first completing construction seed.
greedy_workload(
    bundled_11_best_effort,
    'fixtures/bundled_17_clues.pl', 11, max_size, best_effort,
    'OMEGA POINT', topleft_across, completed(1, 5),
    20, 3, core).

greedy_workload(
    benchmark_08_candidates,
    'fixtures/benchmark_08_words.pl', 13, size, candidates(strict, 5),
    'AOBSCFDMJJJJV', topleft_across, completed(8, 0),
    20, 3, core).

greedy_workload(
    real_13x13_12w_best_effort,
    'fixtures/real_13x13_12w.pl', 13, size, best_effort,
    'KINKINESSES', topleft_across, completed(8, 4),
    10, 2, core).

greedy_workload(
    real_15x15_18w_best_effort,
    'fixtures/real_15x15_18w.pl', 15, size, best_effort,
    'KATZENJAMMERS', topleft_across, completed(15, 3),
    10, 2, core).

greedy_workload(
    ladder_15x15_32w_best_effort,
    'fixtures/ladder_15x15_32w.pl', 15, size, best_effort,
    'FAFB', topleft_across, completed(28, 4),
    7, 2, heavy).

greedy_workload(
    ladder_21x21_80w_best_effort,
    'fixtures/ladder_21x21_80w.pl', 21, size, best_effort,
    'EDDE', topleft_across, completed(9, 71),
    5, 2, heavy).
