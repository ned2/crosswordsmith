% Focused subprocess test for the greedy ratchet's atomic promote/read-back path.

:- use_module(library(plunit)).
:- use_module(library(process)).
:- use_module('../benchmarks/greedy_subjects.pl').
:- use_module('../benchmarks/probe_arrange/g2_transpose.pl').

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

test(semantic_counter_partitions) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    greedy_subjects:semantic_counters(
        Words, 13, candidates(strict, 5), Counters),
    assertion(Counters.generated_crossing_descriptors =:=
              Counters.legality_probes),
    assertion(Counters.legality_successes + Counters.legality_rejects =:=
              Counters.legality_probes),
    assertion(Counters.legality_successes =:= Counters.scored_candidates),
    assertion(Counters.start_lt_one =< Counters.legality_rejects),
    assertion(Counters.generated_crossing_descriptors =:= 1308),
    assertion(Counters.legality_rejects =:= 1134),
    assertion(Counters.scored_candidates =:= 174),
    assertion(Counters.score_cell_visits =:= 4524),
    assertion(Counters.rejected_score_cell_visits_avoided =:= 3640),
    assertion(Counters.start_lt_one =:= 994),
    assertion(Counters.direct_attempted_constructions =:= 10),
    assertion(Counters.completed_constructions =:= 10).

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

test(g2_transpose_is_involution_with_fresh_clue_vars) :-
    Placed = [pw('CAT',[c,a,t],[1,2,3],across,3,1,3,SourceNum)],
    g2_transpose_probe:transpose_placed(Placed, 7, Once),
    g2_transpose_probe:transpose_placed(Once, 7, Twice),
    Once = [pw('CAT',[c,a,t],[1,8,15],down,3,1,15,OnceNum)],
    Twice = [pw('CAT',[c,a,t],[1,2,3],across,3,1,3,TwiceNum)],
    assertion(var(SourceNum)), assertion(var(OnceNum)), assertion(var(TwiceNum)),
    assertion(SourceNum \== OnceNum), assertion(SourceNum \== TwiceNum),
    assertion(OnceNum \== TwiceNum).

test(g2_product_transpose_geometry_and_fresh_clue_vars) :-
    Placed = [pw('CAT',[c,a,t],[1,2,3],across,3,1,3,SourceNum1),
              pw('CAR',[c,a,r],[1,8,15],down,3,1,15,SourceNum2)],
    crosswordsmith_arrange:transpose_placed(Placed, 7, Once),
    crosswordsmith_arrange:transpose_placed(Once, 7, Twice),
    Once = [pw('CAT',[c,a,t],[1,8,15],down,3,1,15,OnceNum1),
            pw('CAR',[c,a,r],[1,2,3],across,3,1,3,OnceNum2)],
    Twice = [pw('CAT',[c,a,t],[1,2,3],across,3,1,3,TwiceNum1),
             pw('CAR',[c,a,r],[1,8,15],down,3,1,15,TwiceNum2)],
    term_variables([SourceNum1,SourceNum2,OnceNum1,OnceNum2,
                    TwiceNum1,TwiceNum2], ClueVars),
    assertion(length(ClueVars, 6)).

test(g2_product_preserves_four_corner_block_order) :-
    Words = [['CAT', _{sample:only}]],
    crosswordsmith_arrange:greedy_constructions(Words, 7, Constructions),
    maplist(construction_seed_start_dir, Constructions, StartsDirs),
    assertion(StartsDirs == [1-across,1-down,7-down,43-across]).

test(g2_product_omits_symmetric_setup_failures) :-
    g2_transpose_probe:additional_sample(
        setup_failure, Words, 3, _Command, setup_failure),
    crosswordsmith_arrange:greedy_constructions(Words, 3, Constructions),
    assertion(length(Constructions, 8)),
    maplist(construction_seed_start_dir, Constructions, StartsDirs),
    assertion(StartsDirs == [1-across,1-across,1-down,1-down,
                             3-down,3-down,7-across,7-across]).

test(g2_product_preserves_dropped_terms_order_and_fresh_copies) :-
    Words = [['CAT',meta(CatTag)],['CAR',meta(CarTag)],['DOG',meta(DogTag)]],
    crosswordsmith_arrange:greedy_constructions(Words, 7, Constructions),
    Constructions = [gc(_,DroppedSource),gc(_,_)|_],
    nth0(3, Constructions, gc(_,DroppedPartner)),
    assertion(DroppedSource =@= [['DOG',meta(_)] ]),
    assertion(DroppedPartner =@= [['DOG',meta(_)] ]),
    DroppedSource = [['DOG',meta(SourceTag)]],
    DroppedPartner = [['DOG',meta(PartnerTag)]],
    assertion(SourceTag \== PartnerTag),
    assertion(SourceTag \== DogTag),
    assertion(PartnerTag \== DogTag),
    assertion(var(CatTag)), assertion(var(CarTag)).

test(g2_additional_samples_cover_required_behaviors) :-
    findall(Row,
            ( g2_transpose_probe:additional_sample(Id, Words, Size, Command, Behavior),
              crosswordsmith_arrange:seed_candidates(Words, Seeds),
              length(Seeds, SeedCount), ExpectedSlots is 4 * SeedCount,
              g2_transpose_probe:probe_case(
                  Id, Words, Size, Command, ExpectedSlots, Row0),
              put_dict(behavior, Row0, Behavior, Row) ),
            Rows),
    assertion(length(Rows, 3)),
    assertion(forall(member(R, Rows), get_dict(gate, R, pass))),
    once((member(Setup, Rows), Setup.behavior == setup_failure)),
    assertion(Setup.setup_failures > 0),
    once((member(Dropped, Rows), Dropped.behavior == dropped_word)),
    assertion(Dropped.completed_with_drops > 0),
    once((member(AllFit, Rows), AllFit.behavior == all_fit)),
    assertion(AllFit.all_fit =:= AllFit.completed).

attempt_slot(Attempt, Corner-SeedAnswer) :-
    get_dict(corner, Attempt, Corner),
    get_dict(seed_answer, Attempt, SeedAnswer).

construction_seed_start_dir(gc(Placed, _Dropped), Start-Dir) :-
    last(Placed, PW),
    pw_start(PW, Start),
    pw_dir(PW, Dir).

:- end_tests(greedy_benchmark).
