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

test(candidate_phases_match_product) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    greedy_subjects:build_raw_pool(Words, 13, candidates(strict, 5), Raw),
    assertion(Raw = [_|_]),
    assertion(forall(member(_Score-Placed, Raw), is_list(Placed))),
    once(greedy_subjects:identity_raw_pool(
        Words, 13, candidates(strict, 5), IdentityRaw)),
    assertion(length(IdentityRaw, 4)),
    assertion(forall(member(Row, IdentityRaw),
                     ( get_dict(placed_signature, Row, _),
                       \+ get_dict(dropped_answers, Row, _) ))),
    once(greedy_subjects:postprocess(
        candidates(strict, 5), Raw, 13, 8, 5,
        candidates(TwinLayouts, TwinReturned))),
    once(crosswordsmith_arrange:arrange_candidates(
        Words, 13, strict, 5, ProductLayouts, ProductReturned)),
    assertion(TwinReturned == ProductReturned),
    assertion(TwinLayouts =@= ProductLayouts).

test(best_effort_phases_match_product) :-
    load_clues('fixtures/bundled_17_clues.pl', Words),
    greedy_subjects:build_raw_pool(Words, 11, best_effort, Raw),
    assertion(Raw = [_|_]),
    assertion(forall(member(_Score-pd(_Placed, Dropped), Raw),
                     maplist(atom, Dropped))),
    once(greedy_subjects:identity_raw_pool(
        Words, 11, best_effort, IdentityRaw)),
    assertion(length(IdentityRaw, 8)),
    assertion(forall(member(Row, IdentityRaw),
                     get_dict(dropped_answers, Row, _))),
    length(Words, Total),
    once(greedy_subjects:postprocess(
        best_effort, Raw, 11, Total, 1,
        best_effort(TwinLayout, TwinReward, TwinPlaced, TwinDropped))),
    once(crosswordsmith_arrange:arrange_best_effort(
        Words, 11, ProductLayout, ProductReward, ProductPlaced, ProductDropped)),
    assertion(TwinReward == ProductReward),
    assertion(TwinPlaced == ProductPlaced),
    assertion(TwinDropped == ProductDropped),
    assertion(TwinLayout =@= ProductLayout).

test(direct_attempt_slots_and_strict_omission) :-
    load_clues('fixtures/bundled_17_clues.pl', Words),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    start_locs(Corners),
    findall(Corner-Answer,
            ( member(Corner, Corners), member([Answer|_], Seeds) ),
            ExpectedSlots),
    greedy_subjects:direct_attempts(Words, 11, candidates(strict, 3), Attempts),
    maplist(attempt_slot, Attempts, ActualSlots),
    assertion(ActualSlots == ExpectedSlots),
    findall(A, (member(A, Attempts), get_dict(outcome, A, setup_failed)), Failed),
    findall(A, (member(A, Attempts), get_dict(outcome, A, completed),
                get_dict(eligibility, A, false)), Ineligible),
    assertion(length(Attempts, 20)),
    assertion(length(Failed, 12)),
    assertion(length(Ineligible, 8)),
    assertion(forall(member(A, Failed),
                     ( get_dict(score, A, null),
                       get_dict(placed_signature, A, null),
                       get_dict(dropped_signature, A, null) ))),
    greedy_subjects:identity_raw_pool(
        Words, 11, candidates(strict, 3), RawPool),
    assertion(RawPool == []).

attempt_slot(Attempt, Corner-SeedAnswer) :-
    get_dict(corner, Attempt, Corner),
    get_dict(seed_answer, Attempt, SeedAnswer).

:- end_tests(greedy_benchmark).
