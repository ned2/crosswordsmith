% A-D2 benchmark-only differential seam. It drives the product delta refresh,
% validates every changed bucket/residue against the unchanged full recount,
% and records exact old/full/delta decision traces. Observer mutation never
% enters search decisions.
:- module(probe_arrange_ad2,
          [ ad2_corner/8,
            ad2_operation/3
          ]).

:- use_module(library(assoc), [empty_assoc/1, get_assoc/3, put_assoc/4]).
:- use_module(library(lists), [member/2, reverse/2]).
:- use_module('./d0_support.pl', []).

%!  ad2_corner(+Words:list, +Grid:integer, +Corner:atom, +Seed,
%!             -Result, -Summary:dict, -CountTrace:list,
%!             -DecisionTrace:list) is det.
%
%   Replay one completing strict corner without an inference limit. Every
%   product delta result is checked immediately against mrv_count/8 and the
%   full proof recount before its bucket can influence ordering.
ad2_corner(Words, Grid, Corner, Seed, Result, Summary, CountTrace,
           DecisionTrace) :-
    Holder = holder([], [], [], Counters, []),
    empty_assoc(Counters),
    setup_call_cleanup(
        install_seed(Seed),
        ( crosswordsmith_core:reset_search_memos,
          crosswordsmith_core:init_gs(Grid, GS),
          crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
          crosswordsmith_core:tag_words(Words, Tagged),
          length(Words, N), Arity is N + 1,
          functor(Buckets, buckets, Arity),
          functor(Residues, residues, Arity),
          crosswordsmith_core:state_perturb(none, Perturb),
          once(ad2_assign(Tagged, [],
                          direct(Buckets, Residues, none, Perturb),
                          Grid, Start, Dir, GS, _, Placed, Holder)) ),
        crosswordsmith_core:set_search_seed(-1)),
    crosswordsmith_core:assign_clue_numbers(Placed, Numbered),
    crosswordsmith_arrange:layout_reward(5, 1, Numbered, Reward),
    Result = result(Numbered, Reward),
    holder_output(Holder, Summary, CountTrace, DecisionTrace).

install_seed(null) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(none) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(Seed) :- crosswordsmith_core:set_search_seed(Seed).

%!  ad2_operation(+Words:list, +Grid:integer, -Result) is det.
%
%   Unobserved product-delta two-representative operation for paired wall
%   measurement. It owns one memo reset shared across both corners.
ad2_operation(Words, Grid, Result) :-
    setup_call_cleanup(
        crosswordsmith_core:set_search_seed(-1),
        ( crosswordsmith_core:reset_search_memos,
          ad2_corners([topleft_across,topright], Words, Grid, Tagged),
          ad2_best(Tagged, Result) ),
        crosswordsmith_core:set_search_seed(-1)).

ad2_corners([], _, _, []).
ad2_corners([Corner|Corners], Words, Grid, [Tagged|Rest]) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    (   once(crosswordsmith_core:assign_words_direct(
                 Words, [], Grid, Start, Dir, GS, _, Placed))
    ->  crosswordsmith_arrange:layout_reward(5, 1, Placed, Reward),
        Tagged = ok(Reward, Placed)
    ;   Tagged = exhausted
    ),
    ad2_corners(Corners, Words, Grid, Rest).

ad2_best(Tagged, result(Numbered, Reward)) :-
    findall(R-P, member(ok(R, P), Tagged), OKs),
    (   OKs = [_|_]
    ->  sort(1, @>=, OKs, [Reward-Best|_]),
        crosswordsmith_core:assign_clue_numbers(Best, Numbered)
    ;   Numbered = [], Reward = -1
    ).

ad2_assign([], Placed, _State, _, _, _, GS, GS, Placed, _).
ad2_assign([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
           Holder) :-
    Words = [W|Ws],
    ad2_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
               Tagged, Remaining, State, Holder),
    Tagged = wid(ID, Entry), Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    record(Holder, 3, d(ID, Start, Dir)),
    ad2_assign(Remaining, [PW|Placed0], State, Grid, _S, _D, GS1, GS,
               Placed, Holder).

ad2_select(Words, [], direct(Buckets, Residues, _Last, Perturb),
           _Grid, _Start, _Dir, _GS,
           Tagged, Remaining,
           direct(Buckets, Residues, EntryLetters, Perturb), Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Tagged, Ordered), Tagged = wid(ID, Entry),
    record(Holder, 2, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).
ad2_select(Words, [P|Ps],
           direct(Buckets, Residues, LastLetters, Perturb),
           Grid, Start, Dir, GS,
           Tagged, Remaining,
           direct(Buckets, Residues, EntryLetters, Perturb), Holder) :-
    Placed = [P|Ps],
    verify_refreshes(Words, Buckets, Residues, LastLetters, Placed,
                     Grid, Start, Dir, GS, Holder),
    crosswordsmith_core:direct_partitions(Words, Buckets, Ones, Twos),
    crosswordsmith_core:direct_order(Perturb, Ones, Twos, Ordered),
    snapshot(Words, Buckets, Counts),
    tagged_ids(Ordered, Order),
    record(Holder, 1, node(Counts, Order)),
    member(Tagged, Ordered), Tagged = wid(ID, Entry),
    record(Holder, 2, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).

verify_refreshes([], _, _, _, _, _, _, _, _, _).
verify_refreshes([Tagged|Words], Buckets, Residues, LastLetters, Placed,
                 Grid, Start, Dir, GS, Holder) :-
    Tagged = wid(ID, Entry),
    arg(ID, Buckets, Previous),
    (   var(Previous)
    ->  crosswordsmith_core:refresh_delta_counts(
            [Tagged], Buckets, Residues, LastLetters, Placed,
            Grid, Start, Dir, GS),
        verify_full(ID, Entry, Buckets, Residues, Placed,
                    Grid, Start, Dir, GS, _Count, _Residue)
    ;   crosswordsmith_core:entry_letters(Entry, Letters),
        (   crosswordsmith_core:shares_letter(Letters, LastLetters)
        ->  arg(ID, Residues, OldResidue),
            semantic_work(Previous, Letters, Placed, Grid, GS,
                          ExactChecks, ProposedChecks),
            crosswordsmith_core:refresh_delta_counts(
                [Tagged], Buckets, Residues, LastLetters, Placed,
                Grid, Start, Dir, GS),
            verify_full(ID, Entry, Buckets, Residues, Placed,
                        Grid, Start, Dir, GS, Exact, FullResidue),
            arg(ID, Buckets, Delta), arg(ID, Residues, DeltaResidue),
            refresh_mode(Previous, Mode),
            counter_add(Holder, refreshes, 1),
            counter_add(Holder, ExactChecks, ProposedChecks, Mode),
            record(Holder, 5,
                   refresh(ID, Previous, Exact, Delta, Mode,
                           OldResidue, FullResidue, DeltaResidue))
        ;   crosswordsmith_core:refresh_delta_counts(
                [Tagged], Buckets, Residues, LastLetters, Placed,
                Grid, Start, Dir, GS),
            arg(ID, Buckets, Retained),
            Retained =:= Previous,
            counter_add(Holder, stale_retained, 1)
        )
    ),
    verify_refreshes(Words, Buckets, Residues, LastLetters, Placed,
                     Grid, Start, Dir, GS, Holder).

verify_full(ID, Entry, Buckets, Residues, Placed, Grid, Start, Dir, GS,
            Count, FullResidue) :-
    crosswordsmith_core:mrv_count(
        2, Placed, Grid, Start, Dir, GS, Entry, Authority),
    crosswordsmith_core:direct_count_proofs(
        Entry, Placed, Grid, Start, Dir, GS, Count, FullResidue),
    Authority =:= Count,
    arg(ID, Buckets, Product), Product =:= Authority,
    arg(ID, Residues, ProductResidue),
    ProductResidue == FullResidue.

semantic_work(Previous, Letters, Placed, Grid, GS,
              ExactChecks, ProposedChecks) :-
    Placed = [Newest|_],
    length(Letters, Length),
    probe_arrange_d0:scan_support(
        Letters, Length, Placed, Grid, GS,
        _Exact, ExactChecks, _ExactProofs, _ExactUnique),
    (   Previous =:= 2
    ->  ProposedChecks = ExactChecks
    ;   probe_arrange_d0:scan_support(
            Letters, Length, [Newest], Grid, GS,
            NewestCount, NewestChecks, _NewestProofs, _NewestUnique),
        ( Previous =:= 1, NewestCount < 2 -> Extra = 1 ; Extra = 0 ),
        ProposedChecks is NewestChecks + Extra
    ).

refresh_mode(2, fallback) :- !.
refresh_mode(_, classified).

snapshot([], _, []).
snapshot([wid(ID, _)|Words], Buckets, [ID-Count|Counts]) :-
    arg(ID, Buckets, Count), snapshot(Words, Buckets, Counts).

tagged_ids([], []).
tagged_ids([wid(ID, _)|Words], [ID|IDs]) :- tagged_ids(Words, IDs).

record(Holder, Arg, Event) :-
    arg(Arg, Holder, Events), nb_setarg(Arg, Holder, [Event|Events]).

counter_add(Holder, Key, Add) :-
    arg(4, Holder, C0),
    ( get_assoc(Key, C0, N0) -> true ; N0 = 0 ),
    N is N0 + Add, put_assoc(Key, C0, N, C), nb_setarg(4, Holder, C).

counter_add(Holder, ExactChecks, ProposedChecks, Mode) :-
    counter_add(Holder, exact_candidate_checks, ExactChecks),
    counter_add(Holder, proposed_candidate_checks, ProposedChecks),
    counter_add(Holder, verified_refreshes, 1),
    ( Mode == classified
    -> counter_add(Holder, classified_refreshes, 1)
    ;  counter_add(Holder, fallback_refreshes, 1)
    ).

counter(Holder, Key, Value) :-
    arg(4, Holder, Counters),
    ( get_assoc(Key, Counters, Value) -> true ; Value = 0 ).

holder_output(Holder, Summary, CountTrace, DecisionTrace) :-
    arg(1, Holder, CountRev), reverse(CountRev, Nodes),
    arg(5, Holder, RefreshRev), reverse(RefreshRev, RefreshTrace),
    CountTrace = trace(Nodes, RefreshTrace),
    arg(2, Holder, SelectRev), reverse(SelectRev, Selections),
    arg(3, Holder, DecisionRev), reverse(DecisionRev, Decisions),
    DecisionTrace = trace(Selections, Decisions),
    counter(Holder, refreshes, Refreshes),
    counter(Holder, verified_refreshes, Verified),
    counter(Holder, classified_refreshes, Classified),
    counter(Holder, fallback_refreshes, Fallback),
    counter(Holder, exact_candidate_checks, ExactChecks),
    counter(Holder, proposed_candidate_checks, ProposedChecks),
    counter(Holder, stale_retained, Stale),
    Summary = _{refreshes:Refreshes,verified_refreshes:Verified,
                classified_refreshes:Classified,
                fallback_refreshes:Fallback,
                exact_candidate_checks:ExactChecks,
                proposed_candidate_checks:ProposedChecks,
                stale_retained:Stale}.
