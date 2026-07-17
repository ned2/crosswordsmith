% A-D1 benchmark-only differential replay: unchanged assoc state versus the
% strict direct-bucket representation. Observer holders never enter a decision.
:- module(probe_arrange_ad1,
          [ ad1_corner/8,
            ad1_operation/4
          ]).

:- use_module(library(apply), [include/3]).
:- use_module(library(assoc), [empty_assoc/1, get_assoc/3, put_assoc/4]).
:- use_module(library(lists), [member/2, reverse/2]).
:- use_module(library(pairs), [map_list_to_pairs/3]).

%!  ad1_corner(+Mode:oneof([assoc,direct]), +Words:list, +Grid:integer,
%!             +Corner:atom, +Seed, -Result, -CountTrace:list,
%!             -DecisionTrace:list) is det.
%
%   Replay one completing strict corner without an inference limit. CountTrace
%   records every interior node's remaining ID->visible-bucket correspondence,
%   stable candidate order, and selected IDs. DecisionTrace records every legal
%   placement reached, including attempts later undone by backtracking.
ad1_corner(Mode, Words, Grid, Corner, Seed, Result, CountTrace,
           DecisionTrace) :-
    stable_ids(Words, IDs),
    Holder = holder([], [], []),
    setup_call_cleanup(
        install_seed(Seed),
        ( crosswordsmith_core:reset_search_memos,
          ad1_goal(Mode, Words, Grid, Corner, IDs, Holder, Result) ),
        crosswordsmith_core:set_search_seed(-1)),
    arg(1, Holder, CountRev), reverse(CountRev, CountTrace),
    arg(2, Holder, SelectRev), reverse(SelectRev, Selections),
    arg(3, Holder, DecisionRev), reverse(DecisionRev, Decisions),
    DecisionTrace = trace(Selections, Decisions).

install_seed(null) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(none) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(Seed) :- crosswordsmith_core:set_search_seed(Seed).

%!  ad1_operation(+Mode:oneof([assoc,direct]), +Words:list, +Grid:integer,
%!                -Result) is det.
%
%   Unobserved two-representative operation for secondary paired-wall
%   measurement. Both modes own one memo reset and share it across corners.
ad1_operation(Mode, Words, Grid, Result) :-
    setup_call_cleanup(
        crosswordsmith_core:set_search_seed(-1),
        ( crosswordsmith_core:reset_search_memos,
          plain_corners([topleft_across,topright], Mode, Words, Grid, Tagged),
          plain_best(Tagged, Result) ),
        crosswordsmith_core:set_search_seed(-1)).

plain_corners([], _, _, _, []).
plain_corners([Corner|Corners], Mode, Words, Grid, [Tagged|Rest]) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    (   once(plain_assign(Mode, Words, Grid, Start, Dir, GS, Placed))
    ->  crosswordsmith_arrange:layout_reward(5, 1, Placed, Reward),
        Tagged = ok(Reward, Placed)
    ;   Tagged = exhausted
    ),
    plain_corners(Corners, Mode, Words, Grid, Rest).

plain_assign(assoc, Words, Grid, Start, Dir, GS, Placed) :-
    crosswordsmith_core:assign_words_inc(
        Words, [], none, Grid, Start, Dir, GS, _, Placed).
plain_assign(direct, Words, Grid, Start, Dir, GS, Placed) :-
    crosswordsmith_core:assign_words_direct(
        Words, [], Grid, Start, Dir, GS, _, Placed).

plain_best(Tagged, result(Numbered, Reward)) :-
    findall(R-P, member(ok(R, P), Tagged), OKs),
    (   OKs = [_|_]
    ->  sort(1, @>=, OKs, [Reward-Best|_]),
        crosswordsmith_core:assign_clue_numbers(Best, Numbered)
    ;   Numbered = [], Reward = -1
    ).

ad1_goal(assoc, Words, Grid, Corner, IDs, Holder, Result) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    once(assoc_assign(Words, [], none, Grid, Start, Dir, GS, _, Placed,
                      IDs, Holder)),
    result(Placed, Result).
ad1_goal(direct, Words, Grid, Corner, IDs, Holder, Result) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    crosswordsmith_core:tag_words(Words, Tagged),
    length(Words, N), Arity is N + 1, functor(Buckets, buckets, Arity),
    crosswordsmith_core:state_perturb(none, Perturb),
    once(direct_assign(Tagged, [], direct(Buckets, none, Perturb), Grid,
                       Start, Dir, GS, _, Placed, IDs, Holder)),
    result(Placed, Result).

result(Placed, result(Numbered, Reward)) :-
    crosswordsmith_core:assign_clue_numbers(Placed, Numbered),
    crosswordsmith_arrange:layout_reward(5, 1, Numbered, Reward).

assoc_assign([], Placed, _State, _, _, _, GS, GS, Placed, _, _).
assoc_assign([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
             IDs, Holder) :-
    Words = [W|Ws],
    assoc_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
                 Entry, Remaining, State, IDs, Holder),
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters), length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    id_for(IDs, Answer, ID), record_decision(Holder, d(ID, Start, Dir)),
    assoc_assign(Remaining, [PW|Placed0], State, Grid, _S, _D, GS1, GS,
                 Placed, IDs, Holder).

assoc_select(Words, [], _State, _Grid, _Start, _Dir, _GS,
             Entry, Remaining, none, IDs, Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Entry, Ordered),
    Entry = [Answer|_], id_for(IDs, Answer, ID), record_selection(Holder, ID),
    crosswordsmith_core:remove_x(Entry, Words, Remaining).
assoc_select(Words, [P|Ps], State0, Grid, Start, Dir, GS,
             Entry, Remaining, State, IDs, Holder) :-
    Placed = [P|Ps],
    crosswordsmith_core:state_perturb(State0, Perturb),
    crosswordsmith_core:inc_counts(
        State0, Words, Placed, Grid, Start, Dir, GS, CountMap),
    map_list_to_pairs(crosswordsmith_core:count_of(CountMap), Words, Pairs),
    include(crosswordsmith_core:positive_key, Pairs, Placeable),
    keysort(Placeable, Sorted),
    crosswordsmith_core:order_candidates(Perturb, Sorted, Ordered),
    assoc_snapshot(Words, CountMap, IDs, Counts),
    entry_ids(Ordered, IDs, OrderIDs),
    record_counts(Holder, node(Counts, OrderIDs)),
    member(Entry, Ordered),
    Entry = [Answer|_], id_for(IDs, Answer, ID), record_selection(Holder, ID),
    crosswordsmith_core:remove_x(Entry, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters),
    State = state(CountMap, EntryLetters, Perturb).

direct_assign([], Placed, _State, _, _, _, GS, GS, Placed, _, _).
direct_assign([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
              IDs, Holder) :-
    Words = [W|Ws],
    direct_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
                  Tagged, Remaining, State, Holder),
    Tagged = wid(ID, Entry), Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters), length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    record_decision(Holder, d(ID, Start, Dir)),
    direct_assign(Remaining, [PW|Placed0], State, Grid, _S, _D, GS1, GS,
                  Placed, IDs, Holder).

direct_select(Words, [], direct(Buckets, _Last, Perturb),
              _Grid, _Start, _Dir, _GS,
              Tagged, Remaining, direct(Buckets, EntryLetters, Perturb), Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Tagged, Ordered), Tagged = wid(ID, Entry),
    record_selection(Holder, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).
direct_select(Words, [P|Ps], direct(Buckets, LastLetters, Perturb),
              Grid, Start, Dir, GS,
              Tagged, Remaining, direct(Buckets, EntryLetters, Perturb), Holder) :-
    Placed = [P|Ps],
    crosswordsmith_core:refresh_direct_counts(
        Words, Buckets, LastLetters, Placed, Grid, Start, Dir, GS),
    crosswordsmith_core:direct_partitions(Words, Buckets, Ones, Twos),
    crosswordsmith_core:direct_order(Perturb, Ones, Twos, Ordered),
    direct_snapshot(Words, Buckets, Counts),
    tagged_ids(Ordered, OrderIDs),
    record_counts(Holder, node(Counts, OrderIDs)),
    member(Tagged, Ordered), Tagged = wid(ID, Entry),
    record_selection(Holder, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).

stable_ids(Words, IDs) :-
    empty_assoc(Empty), stable_ids_(Words, 1, Empty, IDs).

stable_ids_([], _, IDs, IDs).
stable_ids_([[Answer|_]|Words], ID, IDs0, IDs) :-
    put_assoc(Answer, IDs0, ID, IDs1),
    ID1 is ID + 1,
    stable_ids_(Words, ID1, IDs1, IDs).

id_for(IDs, Answer, ID) :- get_assoc(Answer, IDs, ID).

assoc_snapshot([], _, _, []).
assoc_snapshot([[Answer|_]|Words], CountMap, IDs, [ID-Count|Counts]) :-
    id_for(IDs, Answer, ID), get_assoc(Answer, CountMap, Count),
    assoc_snapshot(Words, CountMap, IDs, Counts).

direct_snapshot([], _, []).
direct_snapshot([wid(ID, _)|Words], Buckets, [ID-Count|Counts]) :-
    arg(ID, Buckets, Count),
    direct_snapshot(Words, Buckets, Counts).

entry_ids([], _, []).
entry_ids([[Answer|_]|Entries], IDs, [ID|Rest]) :-
    id_for(IDs, Answer, ID), entry_ids(Entries, IDs, Rest).

tagged_ids([], []).
tagged_ids([wid(ID, _)|Tagged], [ID|IDs]) :- tagged_ids(Tagged, IDs).

record_counts(Holder, Event) :-
    arg(1, Holder, Events), nb_setarg(1, Holder, [Event|Events]).
record_selection(Holder, ID) :-
    arg(2, Holder, Events), nb_setarg(2, Holder, [ID|Events]).
record_decision(Holder, Event) :-
    arg(3, Holder, Events), nb_setarg(3, Holder, [Event|Events]).
