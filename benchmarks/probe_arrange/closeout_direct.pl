% Final close-out attribution for the current A-D2 strict direct driver.
% This benchmark-only twin observes events outside product inference limits.
:- module(probe_arrange_closeout_direct,
          [ closeout_corner/7,
            product_corner/4,
            trace_digest/2
          ]).

:- use_module(library(lists), [member/2, reverse/2]).

%!  closeout_corner(+Words:list, +Grid:integer, +Corner:atom, +Seed,
%!                  -Result, -Stats:dict, -Trace) is det.
%
%   Replay one current direct corner without an inference limit. Product
%   refresh_delta_counts/9, direct_partitions/4, direct_order/4, crossing
%   generation, and placement legality make every decision. The holder records
%   semantic events only and cannot steer search.
closeout_corner(Words, Grid, Corner, Seed, Result, Stats, Trace) :-
    Holder = holder(0, 0, 0, 0, 0, 0, [], []),
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
          once(closeout_assign(Tagged, [],
                               direct(Buckets, Residues, none, Perturb),
                               Grid, Start, Dir, GS, _, Placed, 0, Holder)) ),
        crosswordsmith_core:set_search_seed(-1)),
    result(Placed, Result),
    holder_output(Holder, Stats, Trace).

%!  product_corner(+Words:list, +Grid:integer, +Corner:atom, -Result) is det.
%
%   Run the unchanged current direct product driver for one unseeded completing
%   corner, with the same one-reset request lifecycle as attribution rows.
product_corner(Words, Grid, Corner, Result) :-
    setup_call_cleanup(
        crosswordsmith_core:set_search_seed(-1),
        ( crosswordsmith_core:reset_search_memos,
          crosswordsmith_core:init_gs(Grid, GS),
          crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
          once(crosswordsmith_core:assign_words_direct(
                   Words, [], Grid, Start, Dir, GS, _, Placed)) ),
        crosswordsmith_core:set_search_seed(-1)),
    result(Placed, Result).

%!  trace_digest(+Trace, -Digest:atom) is det.
%
%   Stable SHA-1 variant digest of the ground trace used in compact result rows.
trace_digest(Trace, Digest) :-
    variant_sha1(Trace, Digest).

install_seed(null) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(none) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(Seed) :- crosswordsmith_core:set_search_seed(Seed).

result(Placed, result(Numbered, Reward)) :-
    crosswordsmith_core:assign_clue_numbers(Placed, Numbered),
    crosswordsmith_arrange:layout_reward(5, 1, Numbered, Reward).

closeout_assign([], Placed, _State, _, _, _, GS, GS, Placed, Depth, Holder) :-
    update_depth(Holder, Depth).
closeout_assign([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
                Depth, Holder) :-
    Words = [W|Ws],
    bump(Holder, 1),
    update_depth(Holder, Depth),
    closeout_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
                    Tagged, Remaining, State, Holder),
    Tagged = wid(ID, Entry),
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    bump(Holder, 2),
    bump(Holder, 3),
    record(Holder, 8, d(ID, Start, Dir)),
    Depth1 is Depth + 1,
    ( closeout_assign(Remaining, [PW|Placed0], State, Grid, _S, _D,
                      GS1, GS, Placed, Depth1, Holder)
    ; bump(Holder, 4), fail
    ).

closeout_select(Words, [], direct(Buckets, Residues, _Last, Perturb),
                _Grid, _Start, _Dir, _GS,
                Tagged, Remaining,
                direct(Buckets, Residues, EntryLetters, Perturb), Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Tagged, Ordered),
    Tagged = wid(ID, Entry),
    record(Holder, 7, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).
closeout_select(Words, [P|Ps],
                direct(Buckets, Residues, LastLetters, Perturb),
                Grid, Start, Dir, GS,
                Tagged, Remaining,
                direct(Buckets, Residues, EntryLetters, Perturb), Holder) :-
    Placed = [P|Ps],
    crosswordsmith_core:refresh_delta_counts(
        Words, Buckets, Residues, LastLetters, Placed,
        Grid, Start, Dir, GS),
    crosswordsmith_core:direct_partitions(Words, Buckets, Ones, Twos),
    crosswordsmith_core:direct_order(Perturb, Ones, Twos, Ordered),
    ( Ordered == [] -> bump(Holder, 5), fail ; true ),
    member(Tagged, Ordered),
    Tagged = wid(ID, Entry),
    record(Holder, 7, ID),
    crosswordsmith_core:remove_word_id(ID, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters).

bump(Holder, Arg) :-
    arg(Arg, Holder, N0), N is N0 + 1, nb_setarg(Arg, Holder, N).

update_depth(Holder, Depth) :-
    arg(6, Holder, Old),
    ( Depth > Old -> nb_setarg(6, Holder, Depth) ; true ).

record(Holder, Arg, Event) :-
    arg(Arg, Holder, Events), nb_setarg(Arg, Holder, [Event|Events]).

holder_output(Holder, Stats, trace(Selections, Decisions)) :-
    arg(1, Holder, Nodes),
    arg(2, Holder, LegalDecisions),
    arg(3, Holder, Places),
    arg(4, Holder, Unplaces),
    arg(5, Holder, Wipeouts),
    arg(6, Holder, MaxDepth),
    arg(7, Holder, SelectionRev), reverse(SelectionRev, Selections),
    arg(8, Holder, DecisionRev), reverse(DecisionRev, Decisions),
    Stats = stats{nodes:Nodes,legal_decisions:LegalDecisions,places:Places,
                  unplaces:Unplaces,wipeouts:Wipeouts,max_depth:MaxDepth}.
