% Benchmark-only arrange authority and exact-replay mechanism rigs.
% Product modules are called unchanged; no probe hook exists in product code.
:- module(probe_arrange,
          [ load_fixture/3,
            authority_operation/6,
            authority_corner/7,
            twin_operation/7,
            twin_corner/9,
            layout_signature/2,
            trace_row/8
          ]).

:- use_module(library(apply), [foldl/4, include/3, maplist/3]).
:- use_module(library(assoc),
              [assoc_to_list/2, empty_assoc/1, get_assoc/3, put_assoc/4]).
:- use_module(library(lists), [member/2]).
:- use_module(library(pairs), [map_list_to_pairs/3]).
:- use_module(library(readutil), [read_file_to_terms/3]).
:- use_module(library(statistics), [call_time/2]).
% library(lists):list_to_set/2 uses must_be/2 through autoload. The project
% disables autoload in some test roots, so make that library-internal import
% explicit before product pair_crossings/3 reaches it.
:- lists:use_module(library(error), [must_be/2]).

% load_fixture(+File, +ExpectedCount, -Words)
load_fixture(File, ExpectedCount, Words) :-
    read_file_to_terms(File, Terms, []),
    ( memberchk(clues(Words0), Terms) -> Words = Words0
    ; throw(error(fixture_missing_clues(File), _))
    ),
    length(Words, Actual),
    ( Actual =:= ExpectedCount -> true
    ; throw(error(fixture_word_count(File, ExpectedCount, Actual), _))
    ).

% authority_operation(+Words,+Grid,+Budget,+Seed,-Result,-Timing)
% Only unchanged arrange_best_layout/6 executes in the measured call. Its own
% construct_one/7 owns the inference limit, exactly as in product.
authority_operation(Words, Grid, Budget, Seed, Result, Timing) :-
    with_search_seed(Seed,
        call_time(
            crosswordsmith_arrange:arrange_best_layout(
                Words, Grid, Budget, Numbered, Reward, Outcome),
            T)),
    authority_result(Outcome, Numbered, Reward, T.inferences, Result),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

% authority_corner defines a standalone request lifecycle: reset memos once,
% install the seed once, then call the unchanged product construct_one/7.
authority_corner(Words, Grid, Corner, Budget, Seed, Result, Timing) :-
    with_search_seed(Seed,
        ( crosswordsmith_core:reset_search_memos,
          call_time(
              crosswordsmith_arrange:construct_one(
                  Corner, Words, Grid, 5, 1, Budget, Tagged),
              T) )),
    corner_authority_result(Tagged, T.inferences, Result),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

authority_result(placed, Numbered, Reward, Inferences,
                 result(placed, ok, false, Inferences, Numbered, Reward)).
authority_result(not_proven, _, _, _,
                 result(not_proven, budget, true, null, [], null)).
authority_result(infeasible, _, _, _,
                 result(infeasible, exhausted, false, null, [], null)).

corner_authority_result(ok(Reward, Placed), Inferences,
                        result(placed, ok, false, Inferences, Numbered, Reward)) :-
    crosswordsmith_core:assign_clue_numbers(Placed, Numbered).
corner_authority_result(budget, _,
                        result(not_proven, budget, true, null, [], null)).
corner_authority_result(exhausted, _,
                        result(infeasible, exhausted, false, null, [], null)).

with_search_seed(Seed, Goal) :-
    setup_call_cleanup(install_seed(Seed), Goal,
                       crosswordsmith_core:set_search_seed(-1)).

install_seed(null) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(none) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(Seed) :- crosswordsmith_core:set_search_seed(Seed).

% twin_operation runs the same two representatives with one reset and one
% mutable seed stream. It has no inference cutoff and is therefore only an
% exact-replay control, never a substitute for authority at 500M.
twin_operation(Words, Grid, Seed, Mode, Result, Stats, Timing) :-
    init_holder(Mode, none, null, Holder),
    with_search_seed(Seed,
        ( crosswordsmith_core:reset_search_memos,
          call_time(twin_corners([topleft_across,topright], Words, Grid,
                                 Holder, Tagged), T) )),
    twin_best_result(Tagged, Result),
    holder_stats(Holder, Stats),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

twin_corners([], _, _, _, []).
twin_corners([Corner|Corners], Words, Grid, Holder, [Tagged|Rest]) :-
    twin_one(Words, Grid, Corner, Holder, Tagged),
    twin_corners(Corners, Words, Grid, Holder, Rest).

twin_one(Words, Grid, Corner, Holder, Tagged) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    ( once(twin_assign(Words, [], none, Grid, Start, Dir, GS, _, Placed,
                       0, Holder))
    -> crosswordsmith_arrange:layout_reward(5, 1, Placed, Reward),
       Tagged = ok(Reward, Placed)
    ;  Tagged = exhausted
    ).

twin_best_result(Tagged, result(placed, ok, false, null, Numbered, Reward)) :-
    findall(R-P, member(ok(R, P), Tagged), OKs),
    OKs = [_|_], !,
    sort(1, @>=, OKs, [Reward-Best|_]),
    crosswordsmith_core:assign_clue_numbers(Best, Numbered).
twin_best_result(_, result(infeasible, exhausted, false, null, [], null)).

% twin_corner(+Words,+Grid,+Corner,+Seed,+Mode,+LimitKind,+Limit,-Result,-Run)
% LimitKind is none|nodes|decisions. Instrumented limits are semantic work
% limits outside any product inference cutoff.
twin_corner(Words, Grid, Corner, Seed, Mode, LimitKind, Limit, Result,
            run(Stats, Timing)) :-
    init_holder(Mode, LimitKind, Limit, Holder),
    with_search_seed(Seed,
        ( crosswordsmith_core:reset_search_memos,
          call_time(
              catch(twin_corner_goal(Words, Grid, Corner, Holder, Result0),
                    probe_limit(_Kind),
                    Result0 = result(not_proven, budget, true, null, [], null)),
              T) )),
    Result = Result0,
    holder_stats(Holder, Stats),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

twin_corner_goal(Words, Grid, Corner, Holder, Result) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    ( once(twin_assign(Words, [], none, Grid, Start, Dir, GS, _, Placed,
                       0, Holder))
    -> crosswordsmith_core:assign_clue_numbers(Placed, Numbered),
       crosswordsmith_arrange:layout_reward(5, 1, Numbered, Reward),
       Result = result(placed, ok, false, null, Numbered, Reward)
    ;  Result = result(infeasible, exhausted, false, null, [], null)
    ).

% This is a local replay of assign_words_inc/9. The ordering and support logic
% below call product helpers verbatim; only event calls are added around branch
% boundaries. Grid and boundary behavior remain product assign_word/9 behavior.
twin_assign([], Placed, _State, _, _, _, GS, GS, Placed, Depth, Holder) :-
    event_leaf(Holder, Depth).
twin_assign([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
            Depth, Holder) :-
    Words = [W|Ws],
    event_node(Holder, Depth, Words, Placed0, GS0, Parent),
    twin_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
                Entry, Remaining, State, Holder),
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    event_place(Holder, Parent, Answer, Start, Dir),
    Depth1 is Depth + 1,
    ( twin_assign(Remaining, [PW|Placed0], State, Grid, _S, _D,
                  GS1, GS, Placed, Depth1, Holder)
    ; event_unplace(Holder), fail
    ).

% Exact local replay of select_inc/10. Equal-count bucket order, crossing proof
% multiplicity, and capped stale-overcount semantics all remain product code.
twin_select(Words, [], _State0, _Grid, _Start, _Dir, _GS,
            Entry, Remaining, none, _Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Entry, Ordered),
    crosswordsmith_core:remove_x(Entry, Words, Remaining).
twin_select(Words, [P|Ps], State0, Grid, Start, Dir, GS,
            Entry, Remaining, State, Holder) :-
    Placed = [P|Ps],
    crosswordsmith_core:state_perturb(State0, Perturb),
    crosswordsmith_core:inc_counts(
        State0, Words, Placed, Grid, Start, Dir, GS, CountMap),
    event_support(Holder, State0, Words, CountMap),
    map_list_to_pairs(crosswordsmith_core:count_of(CountMap), Words, Pairs),
    include(crosswordsmith_core:positive_key, Pairs, Placeable),
    keysort(Placeable, Sorted),
    crosswordsmith_core:order_candidates(Perturb, Sorted, Ordered),
    ( Ordered == [] -> event_wipeout(Holder), fail ; true ),
    member(Entry, Ordered),
    crosswordsmith_core:remove_x(Entry, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters),
    State = state(CountMap, EntryLetters, Perturb).

% Counter holder fields are intentionally fixed and extensible. Mode=none has
% no event work; lean has depth/churn/decision counters; full adds support,
% duplicate-state, and explicit state-size observations.
init_holder(Mode, Kind, Limit,
            h(Mode,Kind,Limit,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Children,States)) :-
    empty_assoc(Children), empty_assoc(States).

event_leaf(none_holder, _) :- !.
event_leaf(Holder, Depth) :- holder_mode(Holder, Mode), event_leaf_mode(Mode, Holder, Depth).
event_leaf_mode(none, _, _).
event_leaf_mode(lean, Holder, Depth) :- update_depth(Holder, Depth).
event_leaf_mode(full, Holder, Depth) :- update_depth(Holder, Depth).

holder_mode(Holder, Mode) :- arg(1, Holder, Mode).

event_node(Holder, Depth, Words, Placed, GS, Parent) :-
    holder_mode(Holder, Mode),
    event_node_mode(Mode, Holder, Depth, Words, Placed, GS, Parent).
event_node_mode(none, _, _, _, _, _, 0).
event_node_mode(lean, Holder, Depth, _, _, _, Parent) :-
    check_limit(Holder, nodes, 4),
    bump_arg(Holder, 4), update_depth(Holder, Depth), next_parent(Holder, Parent).
event_node_mode(full, Holder, Depth, Words, Placed, GS, Parent) :-
    event_node_mode(lean, Holder, Depth, Words, Placed, GS, Parent),
    state_sizes(Holder, Words, GS),
    state_duplicate(Holder, Words, Placed).

event_place(Holder, Parent, Answer, Start, Dir) :-
    holder_mode(Holder, Mode),
    event_place_mode(Mode, Holder, Parent, Answer, Start, Dir).
event_place_mode(none, _, _, _, _, _).
event_place_mode(lean, Holder, _, _, _, _) :-
    check_limit(Holder, decisions, 5), bump_arg(Holder, 5), bump_arg(Holder, 6).
event_place_mode(full, Holder, Parent, Answer, Start, Dir) :-
    event_place_mode(lean, Holder, Parent, Answer, Start, Dir),
    child_duplicate(Holder, Parent-Answer-Start-Dir).

event_unplace(Holder) :-
    holder_mode(Holder, Mode),
    ( Mode == none -> true ; bump_arg(Holder, 7) ).
event_wipeout(Holder) :-
    holder_mode(Holder, Mode),
    ( Mode == none -> true ; bump_arg(Holder, 8) ).

event_support(Holder, State0, Words, Current) :-
    holder_mode(Holder, Mode),
    ( Mode == full -> support_changes(State0, Words, Current, N), add_arg(Holder, 14, N)
    ; true
    ).
support_changes(none, _, _, 0).
support_changes(state(Previous, _, _), Words, Current, N) :-
    foldl(count_transition(Previous, Current), Words, 0, N).
count_transition(Previous, Current, [Answer|_], N0, N) :-
    ( get_assoc(Answer, Previous, A), get_assoc(Answer, Current, B), A =\= B
    -> N is N0 + 1
    ; N = N0
    ).

check_limit(Holder, Kind, CounterArg) :-
    arg(2, Holder, Configured),
    ( Configured == Kind
    -> arg(3, Holder, Limit), arg(CounterArg, Holder, Used),
       ( Used >= Limit -> throw(probe_limit(Kind)) ; true )
    ; true
    ).

bump_arg(Holder, I) :- arg(I, Holder, A), B is A + 1, nb_setarg(I, Holder, B).
add_arg(Holder, I, Add) :- arg(I, Holder, A), B is A + Add, nb_setarg(I, Holder, B).
max_arg(Holder, I, Value) :-
    arg(I, Holder, Old), ( Value > Old -> nb_setarg(I, Holder, Value) ; true ).
update_depth(Holder, Depth) :- max_arg(Holder, 9, Depth).
next_parent(Holder, Parent) :- bump_arg(Holder, 18), arg(18, Holder, Parent).

state_sizes(Holder, Words, gs(Letters, Boundaries)) :-
    length(Words, Entries), max_arg(Holder, 10, Entries),
    bound_args(Letters, LetterCount), max_arg(Holder, 11, LetterCount),
    bound_args(Boundaries, BoundaryCount), max_arg(Holder, 12, BoundaryCount).
bound_args(Term, Count) :-
    functor(Term, _, Arity), bound_args(1, Arity, Term, 0, Count).
bound_args(I, Arity, _, Count, Count) :- I > Arity, !.
bound_args(I, Arity, Term, Count0, Count) :-
    arg(I, Term, V), ( nonvar(V) -> Count1 is Count0 + 1 ; Count1 = Count0 ),
    I1 is I + 1, bound_args(I1, Arity, Term, Count1, Count).

child_duplicate(Holder, Key) :-
    arg(19, Holder, Seen0),
    ( get_assoc(Key, Seen0, _) -> bump_arg(Holder, 15)
    ; put_assoc(Key, Seen0, true, Seen), nb_setarg(19, Holder, Seen)
    ).

state_duplicate(Holder, Words, Placed) :-
    maplist(entry_answer, Words, Remaining0), msort(Remaining0, Remaining),
    maplist(placed_key, Placed, Placed0), msort(Placed0, PlacedKeys),
    Key = state(Remaining, PlacedKeys),
    arg(20, Holder, Seen0),
    ( get_assoc(Key, Seen0, _) -> bump_arg(Holder, 16)
    ; put_assoc(Key, Seen0, true, Seen), nb_setarg(20, Holder, Seen)
    ).
entry_answer([Answer|_], Answer).
placed_key(PW, Answer-Start-Dir) :-
    crosswordsmith_core:pw_answer(PW, Answer),
    crosswordsmith_core:pw_start(PW, Start),
    crosswordsmith_core:pw_dir(PW, Dir).

holder_stats(Holder, Stats) :-
    holder_mode(Holder, Mode),
    holder_stats_mode(Mode, Holder, Stats).
holder_stats_mode(none, _,
    _{nodes:null,decisions:null,places:null,unplaces:null,wipeouts:null,
      max_depth:null,state_entries_max:null,letter_cells_max:null,
      boundary_cells_max:null,support_transitions:null,
      duplicate_children:null,duplicate_states:null}).
holder_stats_mode(Mode, Holder, Stats) :-
    Mode \== none,
    arg(4,Holder,Nodes), arg(5,Holder,Decisions), arg(6,Holder,Places),
    arg(7,Holder,Unplaces), arg(8,Holder,Wipeouts), arg(9,Holder,Depth),
    ( Mode == full
    -> arg(10,Holder,Entries), arg(11,Holder,Letters), arg(12,Holder,Bounds),
       arg(14,Holder,Transitions), arg(15,Holder,DupChildren),
       arg(16,Holder,DupStates)
    ; Entries=null, Letters=null, Bounds=null, Transitions=null,
      DupChildren=null, DupStates=null
    ),
    Stats = _{nodes:Nodes,decisions:Decisions,places:Places,
              unplaces:Unplaces,wipeouts:Wipeouts,max_depth:Depth,
              state_entries_max:Entries,letter_cells_max:Letters,
              boundary_cells_max:Bounds,support_transitions:Transitions,
              duplicate_children:DupChildren,duplicate_states:DupStates}.

layout_signature([], null).
layout_signature(Numbered, Signature) :-
    Numbered = [_|_],
    maplist(placed_key, Numbered, Keys0), msort(Keys0, Keys),
    term_string(Keys, Signature, [quoted(true), numbervars(true)]).

% trace_row(+Rig,+Meta,+Result,+Stats,+Timing,+Commit,+Swi,-Row)
trace_row(Rig, Meta,
          result(Outcome,Termination,Censored,SuccessInf,Numbered,Reward),
          Stats, Timing, Commit, Swi, Row) :-
    layout_signature(Numbered, Signature),
    Base = _{rig:Rig,limit_kind:Meta.limit_kind,
             operation_id:Meta.operation_id,attempt_index:Meta.attempt_index,
             fixture:Meta.fixture,fixture_seed:Meta.fixture_seed,
             search_seed:Meta.search_seed,corner:Meta.corner,arm:Meta.arm,
             cutoff:Meta.cutoff,termination:Termination,outcome:Outcome,
             success_inferences:SuccessInf,
             censored:Censored,max_depth:Stats.max_depth,places:Stats.places,
             unplaces:Stats.unplaces,wipeouts:Stats.wipeouts,reward:Reward,
             layout_signature:Signature,swi_version:Swi,commit:Commit,
             measured_inferences:Timing.inferences,wall_seconds:Timing.wall,
             cpu_seconds:Timing.cpu,nodes:Stats.nodes,decisions:Stats.decisions,
             state_entries_max:Stats.state_entries_max,
             letter_cells_max:Stats.letter_cells_max,
             boundary_cells_max:Stats.boundary_cells_max,
             support_transitions:Stats.support_transitions,
             duplicate_children:Stats.duplicate_children,
             duplicate_states:Stats.duplicate_states},
    Row = Base.

:- multifile prolog:error_message//1.
prolog:error_message(fixture_missing_clues(File)) -->
    ['probe fixture ~q does not define clues/1'-[File]].
prolog:error_message(fixture_word_count(File, Expected, Actual)) -->
    ['probe fixture ~q has ~d words; expected exactly ~d'-[File,Actual,Expected]].
