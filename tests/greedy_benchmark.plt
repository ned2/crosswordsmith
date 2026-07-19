% Focused subprocess test for the greedy ratchet's atomic promote/read-back path.

:- use_module(library(plunit)).
:- use_module('../benchmarks/bench_process.pl', [capture_process/6]).
:- use_module('../benchmarks/greedy_subjects.pl').

:- begin_tests(greedy_benchmark).

test(record_readback_and_retention) :-
    capture_process(path(swipl),
                    ['-q', 'benchmarks/check_greedy_baseline.pl', '--self-test'],
                    inherit, _Stdout, _Stderr, Status),
    assertion(Status == exit(0)).

test(exported_samplers_are_deterministic) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    Command = candidates(strict, 5),
    greedy_subjects:construction_sampler(
        Words, 13, 'AOBSCFDMJJJJV', topleft_across, completed(8, 0),
        _Construction),
    deterministic(ConstructionDet),
    assertion(ConstructionDet == true),
    greedy_subjects:sweep_sampler(Words, 13, Command, _Sweep),
    deterministic(SweepDet),
    assertion(SweepDet == true),
    greedy_subjects:build_raw_pool(Words, 13, Command, Raw),
    deterministic(PoolDet),
    assertion(PoolDet == true),
    greedy_subjects:postprocess_sampler(Raw, 13, 8, Command, 5, _Postprocess),
    deterministic(PostprocessDet),
    assertion(PostprocessDet == true),
    greedy_subjects:command_sampler(
        '/bin/true', ignored, 13, size, Command, 0, _Command),
    deterministic(CommandDet),
    assertion(CommandDet == true).

test(exported_samplers_reject_invalid_inputs,
     [forall(greedy_invalid_sampler(Goal, ExpectedError))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = ExpectedError).

test(identity_row_is_deterministic_without_selected_layouts) :-
    greedy_subjects:identity_row(
        no_layout, 'fixtures/bundled_17_clues.pl',
        'fixtures/bundled_17_clues.pl', 11, size, candidates(strict, 3),
        'OMEGA POINT', topleft_across, completed(1, 5), './crosswordsmith', Row),
    deterministic(Det),
    assertion(Det == true),
    assertion(Row.selected == []),
    assertion(Row.candidate_count =:= 0).

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

greedy_invalid_sampler(Goal, ExpectedError) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    greedy_invalid_sampler_case(Words, Goal, ExpectedError).

greedy_invalid_sampler_case(
    Words,
    greedy_subjects:construction_sampler(
        Words, 13, 'MISSING', topleft_across, completed(8, 0), _),
    error(greedy_seed_answer_missing('MISSING'), _)).
greedy_invalid_sampler_case(
    Words,
    greedy_subjects:sweep_sampler(Words, 13, unsupported, _),
    error(greedy_benchmark_command(unsupported), _)).
greedy_invalid_sampler_case(
    Words,
    greedy_subjects:build_raw_pool(Words, 13, unsupported, _),
    error(greedy_benchmark_command(unsupported), _)).
greedy_invalid_sampler_case(
    _Words,
    greedy_subjects:postprocess_sampler([], 13, 8, best_effort, 1, _),
    error(greedy_empty_raw_pool, _)).
greedy_invalid_sampler_case(
    _Words,
    greedy_subjects:command_sampler(
        '/bin/true', ignored, 13, unsupported, best_effort, 0, _),
    error(greedy_benchmark_framing(unsupported), _)).
greedy_invalid_sampler_case(
    _Words,
    greedy_subjects:identity_row(
        id, fixture, ignored, 13, size, unsupported, seed, corner,
        setup_failed, '/bin/true', _),
    error(greedy_benchmark_command(unsupported), _)).

:- end_tests(greedy_benchmark).
