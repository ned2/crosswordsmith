% A-G2 measurement-only direct greedy transpose probe.
%
% This is deliberately benchmark code. It independently enumerates the old
% four corner-major, seed-major direct construction blocks and compares them to
% the product's synthesized raw pool.

:- module(g2_transpose_probe,
          [ probe_all/1,
            probe_case/6,
            transpose_placed/3,
            additional_sample/5
          ]).

:- use_module(library(apply), [maplist/3]).
:- use_module(library(assoc), [assoc_to_list/2]).
:- use_module(library(lists), [append/2, member/2]).

:- use_module('../greedy_subjects.pl').
:- ensure_loaded('../greedy_workloads.pl').

%!  probe_all(-Document:dict) is det.
%
%   Run the seven-rung manifest and three generated differential cases.
probe_all(_{tool:'crosswordsmith-arrange-g2-transpose-premise', schema:1,
            manifest:Manifest, additional:Additional, gate:Gate}) :-
    findall(Spec, manifest_spec(Spec), Specs),
    maplist(probe_manifest_spec, Specs, Manifest),
    findall(Sample, additional_sample_term(Sample), Samples),
    maplist(probe_additional_sample, Samples, Additional),
    rows_gate(Manifest, ManifestGate),
    rows_gate(Additional, AdditionalGate),
    ( ManifestGate == pass, AdditionalGate == pass -> Gate = pass ; Gate = reject ).

manifest_spec(spec(Id, Fixture, Size, Command, ExpectedSlots)) :-
    greedy_workload(Id, Fixture, Size, _Framing, Command,
                    _Seed, _Corner, _Expected, _Iterations, _Warmup, _Tier),
    expected_manifest_slots(Id, ExpectedSlots).

expected_manifest_slots(bundled_17_candidates, 20).
expected_manifest_slots(bundled_11_best_effort, 20).
expected_manifest_slots(benchmark_08_candidates, 20).
expected_manifest_slots(real_13x13_12w_best_effort, 20).
expected_manifest_slots(real_15x15_18w_best_effort, 16).
expected_manifest_slots(ladder_15x15_32w_best_effort, 8).
expected_manifest_slots(ladder_21x21_80w_best_effort, 4).

probe_manifest_spec(spec(Id, Fixture, Size, Command, ExpectedSlots), Row) :-
    format(user_error, 'heartbeat: A-G2 manifest ~w direct replay~n', [Id]),
    greedy_subjects:load_words(Fixture, Words),
    probe_case(Id, Words, Size, Command, ExpectedSlots, Row).

probe_additional_sample(sample(Id, Words, Size, Command, Behavior), Row) :-
    format(user_error, 'heartbeat: A-G2 differential ~w direct replay~n', [Id]),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    length(Seeds, SeedCount),
    ExpectedSlots is 4 * SeedCount,
    probe_case(Id, Words, Size, Command, ExpectedSlots, Row0),
    behavior_observed(Behavior, Row0, Observed),
    ( Row0.gate == pass, Observed == true -> Gate = pass ; Gate = reject ),
    put_dict(_{behavior:Behavior, behavior_observed:Observed, gate:Gate}, Row0, Row).

% Deterministic generated word sets. Metadata is intentionally retained so the
% dropped-term checks exercise original input terms rather than answer-only data.
%!  additional_sample(-Id:atom, -Words:list, -GridLen:integer, -Command,
%!                    -Behavior:atom) is multi.
%
%   The three campaign-only differential subjects and their expected behavior.
additional_sample(setup_failure,
                  [['ABCDE', _{sample:long}],
                   ['AXE', _{sample:short_a}],
                   ['BEE', _{sample:short_b}]],
                  3, candidates(strict, 3), setup_failure).
additional_sample(dropped_word,
                  [['CAT', _{sample:connected_a}],
                   ['CAR', _{sample:connected_b}],
                   ['DOG', _{sample:isolated}]],
                  7, best_effort, dropped_word).
additional_sample(all_fit,
                  [['CAT', _{sample:only}]],
                  7, candidates(strict, 2), all_fit).

additional_sample_term(sample(Id, Words, Size, Command, Behavior)) :-
    additional_sample(Id, Words, Size, Command, Behavior).

behavior_observed(setup_failure, Row, True) :-
    bool(Row.setup_failures > 0, True).
behavior_observed(dropped_word, Row, True) :-
    bool(Row.completed_with_drops > 0, True).
behavior_observed(all_fit, Row, True) :-
    bool(Row.all_fit =:= Row.completed, True).

bool(Goal, true) :- call(Goal), !.
bool(_Goal, false).

%!  probe_case(+Id:atom, +Words:list, +GridLen:integer, +Command,
%!             +ExpectedSlots:integer, -Row:dict) is det.
%
%   Compare all direct transpose-pair slots for one unchanged greedy sweep.
probe_case(Id, Words, GridLen, Command, ExpectedSlots, Row) :-
    crosswordsmith_core:reset_search_memos,
    tagged_direct_attempts(Words, GridLen, Command, Attempts, SeedCount),
    length(Attempts, DirectSlots),
    split_blocks(SeedCount, Attempts, TLA, TLD, TR, BL),
    compare_blocks(topleft_across, topleft_down, GridLen, Words,
                   TLA, TLD, PairA, MismatchesA),
    compare_blocks(topright, bottomleft, GridLen, Words,
                   TR, BL, PairB, MismatchesB),
    append(MismatchesA, MismatchesB, Mismatches),
    attempt_stats(Attempts, Stats),
    raw_pool_check(Words, GridLen, Command, Attempts, RawEntries, RawMatch),
    bool(DirectSlots =:= ExpectedSlots, SlotCountMatch),
    ( Mismatches == [], RawMatch == true, SlotCountMatch == true
    -> Gate = pass
    ;  Gate = reject
    ),
    command_json(Command, Mode),
    Row = _{id:Id, size:GridLen, words:Stats.words, seeds:SeedCount, mode:Mode,
            direct_slots:DirectSlots, expected_slots:ExpectedSlots,
            slot_count_match:SlotCountMatch, pairs:[PairA,PairB],
            setup_failures:Stats.setup_failures, completed:Stats.completed,
            eligible:Stats.eligible, ineligible:Stats.ineligible,
            completed_with_drops:Stats.completed_with_drops,
            all_fit:Stats.all_fit, raw_pool_entries:RawEntries,
            raw_pool_match:RawMatch, mismatches:Mismatches, gate:Gate}.

command_json(candidates(strict, K), _{kind:candidates,drop:strict,k:K}).
command_json(best_effort, _{kind:best_effort}).

tagged_direct_attempts(Words, GridLen, Command, Attempts, SeedCount) :-
    crosswordsmith_arrange:arrange_weights(WCap, WTail),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    length(Seeds, SeedCount),
    length(Words, Total),
    start_locs(Corners),
    direct_corners(Corners, Seeds, 1, Words, GridLen, Command, Total,
                   WCap, WTail, Attempts, _).

direct_corners([], _Seeds, Slot, _Words, _GridLen, _Command, _Total,
               _WCap, _WTail, [], Slot).
direct_corners([Corner|Corners], Seeds, Slot0, Words, GridLen, Command, Total,
               WCap, WTail, Attempts, Slot) :-
    direct_seeds(Seeds, 1, Slot0, Corner, Words, GridLen, Command, Total,
                 WCap, WTail, Here, Slot1),
    direct_corners(Corners, Seeds, Slot1, Words, GridLen, Command, Total,
                   WCap, WTail, Rest, Slot),
    append(Here, Rest, Attempts).

direct_seeds([], _SeedIndex, Slot, _Corner, _Words, _GridLen, _Command,
             _Total, _WCap, _WTail, [], Slot).
direct_seeds([Seed|Seeds], SeedIndex, Slot0, Corner, Words, GridLen, Command,
             Total, WCap, WTail,
             [attempt(Slot0,Corner,SeedIndex,Seed,Result)|Attempts], Slot) :-
    direct_result(Words, GridLen, Command, Total, WCap, WTail,
                  Corner, Seed, Result),
    Slot1 is Slot0 + 1,
    SeedIndex1 is SeedIndex + 1,
    direct_seeds(Seeds, SeedIndex1, Slot1, Corner, Words, GridLen, Command,
                 Total, WCap, WTail, Attempts, Slot).

direct_result(Words, GridLen, Command, Total, WCap, WTail,
              Corner, Seed, Result) :-
    (   crosswordsmith_arrange:greedy_construct(
            Words, GridLen, Corner, Seed, Placed, Dropped)
    ->  length(Placed, NP), length(Dropped, ND),
        crosswordsmith_arrange:layout_reward(WCap, WTail, Placed, Reward),
        eligibility(Command, NP, Total, Eligible),
        Result = completed(Placed,Dropped,NP,ND,Reward,Eligible)
    ;   Result = setup_failed
    ).

eligibility(candidates(strict, _), NP, Total, Eligible) :-
    ( NP =:= Total -> Eligible = true ; Eligible = false ).
eligibility(best_effort, _NP, _Total, true).

split_blocks(N, Attempts, A, B, C, D) :-
    length(A, N), append(A, R1, Attempts),
    length(B, N), append(B, R2, R1),
    length(C, N), append(C, D, R2),
    length(D, N).

compare_blocks(SourceCorner, PartnerCorner, GridLen, Words,
               Sources, Partners, Pair, Mismatches) :-
    compare_slots(Sources, Partners, GridLen, Words, Mismatches),
    length(Sources, Checked),
    length(Mismatches, MismatchCount),
    SlotsChecked is Checked * 2,
    Pair = _{source:SourceCorner, partner:PartnerCorner,
             pairs_checked:Checked, slots_checked:SlotsChecked,
             mismatches:MismatchCount}.

compare_slots([], [], _GridLen, _Words, []).
compare_slots([Source|Sources], [Partner|Partners], GridLen, Words, Mismatches) :-
    compare_slot(Source, Partner, GridLen, Words, Result),
    compare_slots(Sources, Partners, GridLen, Words, Rest),
    add_mismatch(Result, Rest, Mismatches).

add_mismatch(ok, Rest, Rest).
add_mismatch(mismatch(M), Rest, [M|Rest]).

compare_slot(attempt(SourceSlot,SourceCorner,SeedIndex,SourceSeed,SourceResult),
             attempt(PartnerSlot,PartnerCorner,SeedIndex2,PartnerSeed,PartnerResult),
             GridLen, Words, Result) :-
    SourceSeed = [SourceAnswer|_],
    PartnerSeed = [PartnerAnswer|_],
    (   SeedIndex =\= SeedIndex2
    ->  Field = seed_index(SeedIndex,SeedIndex2)
    ;   SourceSeed \== PartnerSeed
    ->  Field = seed_term(SourceAnswer,PartnerAnswer)
    ;   first_result_mismatch(SourceResult, PartnerResult, GridLen, Words, Field)
    ->  true
    ;   Field = none
    ),
    (   Field == none
    ->  Result = ok
    ;   Result = mismatch(_{source_slot:SourceSlot,partner_slot:PartnerSlot,
                            source_corner:SourceCorner,partner_corner:PartnerCorner,
                            seed_index:SeedIndex,seed_answer:SourceAnswer,
                            first_field:Field})
    ).

first_result_mismatch(setup_failed, setup_failed, _GridLen, _Words, _) :- !, fail.
first_result_mismatch(setup_failed, completed(_,_,_,_,_,_), _GridLen, _Words,
                      outcome(setup_failed,completed)) :- !.
first_result_mismatch(completed(_,_,_,_,_,_), setup_failed, _GridLen, _Words,
                      outcome(completed,setup_failed)) :- !.
first_result_mismatch(completed(P1,D1,NP1,ND1,R1,E1),
                      completed(P2,D2,NP2,ND2,R2,E2), GridLen, Words, Field) :-
    (   E1 \== E2 -> Field = eligibility(E1,E2)
    ;   NP1 =\= NP2 -> Field = placed_count(NP1,NP2)
    ;   ND1 =\= ND2 -> Field = dropped_count(ND1,ND2)
    ;   R1 =\= R2 -> Field = reward(R1,R2)
    ;   placed_first_mismatch(P1, P2, GridLen, Field0)
    ->  Field = Field0
    ;   D1 \== D2 -> Field = dropped_original_terms_or_order
    ;   \+ dropped_subsequence_eq(D1, Words) -> Field = source_dropped_not_original_subsequence
    ;   \+ dropped_subsequence_eq(D2, Words) -> Field = partner_dropped_not_original_subsequence
    ;   \+ normalized_transpose(P1, P2, GridLen) -> Field = normalized_assoc_direction
    ;   \+ fresh_clue_vars(P1, P2) -> Field = clue_number_variable_freshness
    ;   \+ transpose_involution(P1, GridLen) -> Field = transpose_involution
    ;   fail
    ).

placed_first_mismatch(P1, P2, _GridLen, placed_list_length(N1,N2)) :-
    length(P1, N1), length(P2, N2), N1 =\= N2, !.
placed_first_mismatch(P1, P2, GridLen, Field) :-
    placed_word_mismatch(P1, P2, GridLen, 1, Field).

placed_word_mismatch([], [], _GridLen, _Index, _) :- fail.
placed_word_mismatch([PW1|Rest1], [PW2|Rest2], GridLen, Index, Field) :-
    (   pw_field_mismatch(PW1, PW2, GridLen, Index, Field)
    ->  true
    ;   Index1 is Index + 1,
        placed_word_mismatch(Rest1, Rest2, GridLen, Index1, Field)
    ).

pw_field_mismatch(pw(A1,L1,C1,D1,Len1,S1,E1,N1),
                  pw(A2,L2,C2,D2,Len2,S2,E2,N2), GridLen, Index, Field) :-
    maplist(transpose_cell(GridLen), C1, ExpectedCells),
    transpose_cell(GridLen, S1, ExpectedStart),
    transpose_cell(GridLen, E1, ExpectedEnd),
    transpose_dir(D1, ExpectedDir),
    (   A1 \== A2 -> Field = placed(Index,answer,A1,A2)
    ;   L1 \== L2 -> Field = placed(Index,letters,L1,L2)
    ;   C2 \== ExpectedCells -> Field = placed(Index,cells,ExpectedCells,C2)
    ;   D2 \== ExpectedDir -> Field = placed(Index,direction,ExpectedDir,D2)
    ;   Len1 =\= Len2 -> Field = placed(Index,length,Len1,Len2)
    ;   S2 =\= ExpectedStart -> Field = placed(Index,start,ExpectedStart,S2)
    ;   E2 =\= ExpectedEnd -> Field = placed(Index,end,ExpectedEnd,E2)
    ;   \+ var(N1) -> Field = placed(Index,source_clue_number_not_var)
    ;   \+ var(N2) -> Field = placed(Index,partner_clue_number_not_var)
    ;   N1 == N2 -> Field = placed(Index,clue_number_alias)
    ).

%!  transpose_placed(+Placed:list, +GridLen:integer, -Transposed:list) is det.
%
%   Apply the independent square-grid transpose, allocating fresh clue vars.
transpose_placed([], _GridLen, []).
transpose_placed([PW|PWs], GridLen, [TPW|TPWs]) :-
    transpose_pw(PW, GridLen, TPW),
    transpose_placed(PWs, GridLen, TPWs).

transpose_pw(pw(A,L,Cells,Dir,Len,Start,End,_SourceNum), GridLen,
             pw(A,L,TCells,TDir,Len,TStart,TEnd,_FreshNum)) :-
    maplist(transpose_cell(GridLen), Cells, TCells),
    transpose_dir(Dir, TDir),
    transpose_cell(GridLen, Start, TStart),
    transpose_cell(GridLen, End, TEnd).

transpose_cell(GridLen, Cell, Transposed) :-
    Row is (Cell - 1) // GridLen,
    Col is (Cell - 1) mod GridLen,
    Transposed is Col * GridLen + Row + 1.

transpose_dir(across, down).
transpose_dir(down, across).

transpose_involution(Placed, GridLen) :-
    transpose_placed(Placed, GridLen, Once),
    transpose_placed(Once, GridLen, Twice),
    geometry_equal(Placed, Twice),
    fresh_clue_vars(Placed, Once),
    fresh_clue_vars(Placed, Twice),
    fresh_clue_vars(Once, Twice).

geometry_equal([], []).
geometry_equal([pw(A,L,C,D,Len,S,E,_)|Ps],
               [pw(A2,L2,C2,D2,Len2,S2,E2,_)|Qs]) :-
    A == A2, L == L2, C == C2, D == D2,
    Len =:= Len2, S =:= S2, E =:= E2,
    geometry_equal(Ps, Qs).

fresh_clue_vars(P1, P2) :-
    clue_vars(P1, N1), clue_vars(P2, N2),
    append(N1, N2, All),
    all_fresh_vars(All).

clue_vars([], []).
clue_vars([pw(_,_,_,_,_,_,_,N)|PWs], [N|Ns]) :- clue_vars(PWs, Ns).

all_fresh_vars([]).
all_fresh_vars([N|Ns]) :-
    var(N), none_identical(N, Ns), all_fresh_vars(Ns).

none_identical(_N, []).
none_identical(N, [M|Ms]) :- N \== M, none_identical(N, Ms).

normalized_transpose(P1, P2, GridLen) :-
    crosswordsmith_arrange:placement_assoc(P1, GridLen, A1),
    crosswordsmith_arrange:placement_assoc(P2, GridLen, A2),
    assoc_to_list(A1, L1), assoc_to_list(A2, L2),
    maplist(transpose_assoc_pair, L1, Expected),
    Expected == L2.

transpose_assoc_pair(Answer-(Row-Col-Dir), Answer-(Col-Row-TDir)) :-
    transpose_dir(Dir, TDir).

dropped_subsequence_eq([], _Words).
dropped_subsequence_eq([Dropped|DroppedRest], [Word|Words]) :-
    (   Dropped == Word
    ->  dropped_subsequence_eq(DroppedRest, Words)
    ;   dropped_subsequence_eq([Dropped|DroppedRest], Words)
    ).

attempt_stats(Attempts,
              _{words:Words,setup_failures:Setup,completed:Completed,
                eligible:Eligible,ineligible:Ineligible,
                completed_with_drops:WithDrops,all_fit:AllFit}) :-
    Attempts = [attempt(_,_,_,_,FirstResult)|_],
    result_total_words(FirstResult, Attempts, Words),
    count_attempts(Attempts, setup, Setup),
    count_attempts(Attempts, completed, Completed),
    count_attempts(Attempts, eligible, Eligible),
    count_attempts(Attempts, ineligible, Ineligible),
    count_attempts(Attempts, with_drops, WithDrops),
    count_attempts(Attempts, all_fit(Words), AllFit).

result_total_words(completed(_,_,NP,ND,_,_), _Attempts, Words) :- !,
    Words is NP + ND.
result_total_words(setup_failed, Attempts, Words) :-
    member(attempt(_,_,_,_,completed(_,_,NP,ND,_,_)), Attempts), !,
    Words is NP + ND.

count_attempts(Attempts, Kind, Count) :-
    findall(1, (member(Attempt, Attempts), attempt_kind(Kind, Attempt)), Ones),
    length(Ones, Count).

attempt_kind(setup, attempt(_,_,_,_,setup_failed)).
attempt_kind(completed, attempt(_,_,_,_,completed(_,_,_,_,_,_))).
attempt_kind(eligible, attempt(_,_,_,_,completed(_,_,_,_,_,true))).
attempt_kind(ineligible, attempt(_,_,_,_,completed(_,_,_,_,_,false))).
attempt_kind(with_drops, attempt(_,_,_,_,completed(_,_,_,ND,_,_))) :- ND > 0.
attempt_kind(all_fit(Total), attempt(_,_,_,_,completed(_,_,NP,0,_,_))) :- NP =:= Total.

raw_pool_check(Words, GridLen, Command, Attempts, RawEntries, Match) :-
    expected_raw(Command, Attempts, Expected),
    greedy_subjects:build_raw_pool(Words, GridLen, Command, Actual),
    length(Actual, RawEntries),
    ( Expected =@= Actual -> Match = true ; Match = false ).

expected_raw(candidates(_, _), Attempts, Raw) :-
    findall(score(NP,Reward)-Placed,
            member(attempt(_,_,_,_,completed(Placed,_Dropped,NP,_ND,Reward,true)),
                   Attempts),
            Raw).
expected_raw(best_effort, Attempts, Raw) :-
    findall(score(NP,Reward)-pd(Placed,DroppedAnswers),
            ( member(attempt(_,_,_,_,completed(Placed,Dropped,NP,_ND,Reward,true)),
                     Attempts),
              maplist(entry_answer, Dropped, DroppedAnswers) ),
            Raw).

entry_answer([Answer|_], Answer).

rows_gate([], pass).
rows_gate([Row|Rows], Gate) :-
    ( Row.gate == pass -> rows_gate(Rows, Gate) ; Gate = reject ).
