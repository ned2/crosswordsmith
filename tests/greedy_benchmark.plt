% Focused subprocess test for the greedy ratchet's atomic promote/read-back path.

:- use_module(library(plunit)).
:- use_module(library(process)).
:- use_module('../benchmarks/greedy_subjects.pl').

:- begin_tests(greedy_benchmark).

test(record_readback_and_retention) :-
    process_create(path(swipl),
                   ['-q','benchmarks/check_greedy_baseline.pl','--self-test'],
                   [stdout(null),stderr(std),process(PID)]),
    process_wait(PID,exit(0)).

test(replay_equivalence_easy) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    greedy_subjects:replay_equivalent(
        Words, 13, topleft_across, 'AOBSCFDMJJJJV').

% The far-corner EDDE construction places 78/80 words, exercising the dense
% scorer/legality loop rather than only the easy top-left 9/80 construction.
test(replay_equivalence_heavy) :-
    load_clues('fixtures/ladder_21x21_80w.pl', Words),
    greedy_subjects:replay_equivalent(Words, 21, topright, 'EDDE').

:- end_tests(greedy_benchmark).
