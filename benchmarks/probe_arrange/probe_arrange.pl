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

:- dynamic pc0_exact/5.
:- dynamic pc0_simple/4.
:- dynamic pc0_decision/3.

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
    init_holder(Mode, none, null, Words, Holder),
    setup_call_cleanup(
        true,
        ( with_search_seed(Seed,
              ( crosswordsmith_core:reset_search_memos,
                call_time(twin_corners([topleft_across,topright], Words, Grid,
                                       Holder, Tagged), T) )),
          twin_best_result(Tagged, Result),
          holder_stats(Holder, Stats) ),
        cleanup_holder(Holder)),
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
    init_holder(Mode, LimitKind, Limit, Words, Holder),
    setup_call_cleanup(
        true,
        ( with_search_seed(Seed,
              ( crosswordsmith_core:reset_search_memos,
                call_time(
                    catch(twin_corner_goal(Words, Grid, Corner, Holder, Result0),
                          probe_limit(_Kind),
                          Result0 = result(not_proven,budget,true,null,[],null)),
                    T) )),
          Result = Result0,
          holder_stats(Holder, Stats) ),
        cleanup_holder(Holder)),
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
twin_assign(Words, Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
            Depth, Holder) :-
    twin_assign_(Words, Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
                 Depth, Holder, false, false).

twin_assign_([], Placed, _State, _, _, _, GS, GS, Placed, Depth, Holder,
             _StateRepeat, _ParentRepeat) :-
    event_leaf(Holder, Depth).
twin_assign_([W|Ws], Placed0, State0, Grid, Start, Dir, GS0, GS, Placed,
             Depth, Holder, StateRepeat0, ParentRepeat0) :-
    Words = [W|Ws],
    event_node(Holder, Depth, Grid, Words, Placed0, GS0, Parent,
               StateKey, StateKnownDead, NodeBefore),
    repeat_context(StateRepeat0, StateKnownDead, StateRepeat, StateRepeatRoot),
    ( twin_assign_body(Words, Placed0, State0, Grid, Start, Dir, GS0,
                       GS, Placed, Depth, Holder, Parent, StateRepeat,
                       ParentRepeat0)
    ; event_node_exhausted(Holder, StateKey, StateRepeatRoot, NodeBefore), fail
    ).

twin_assign_body(Words, Placed0, State0, Grid, Start, Dir, GS0,
                 GS, Placed, Depth, Holder, Parent, StateRepeat,
                 ParentRepeat0) :-
    twin_select(Words, Placed0, State0, Grid, Start, Dir, GS0,
                Entry, Remaining, State, Holder),
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    word_id(Holder, Answer, WordID),
    event_proof(Holder, Parent, WordID, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    event_child(Holder, Parent, WordID, Start, Dir, ParentRepeat0,
                ParentRepeat, ParentRepeatRoot, ChildBefore),
    Depth1 is Depth + 1,
    ( twin_assign_(Remaining, [PW|Placed0], State, Grid, _S, _D,
                   GS1, GS, Placed, Depth1, Holder, StateRepeat, ParentRepeat)
    ; event_child_exhausted(Holder, Parent, WordID, Start, Dir,
                            ParentRepeatRoot, ChildBefore),
      event_unplace(Holder), fail
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
% no event work; lean has depth/churn/decision counters; full adds exact
% duplicate/dead-state observation. Fingerprints only select a bucket: every
% possible hit is verified against the bucket's canonical keys with (==).
init_holder(Mode, Kind, Limit, Words,
            h(Mode,Kind,Limit,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              Run,none,0,0,0,0,0,0,0,0,0,none,none,IDs,none,0,none)) :-
    flag(pc0_probe_run, Run, Run+1), stable_word_ids(Words, IDs).

cleanup_holder(Holder) :-
    arg(19, Holder, Run),
    retractall(pc0_exact(Run, _, _, _, _)),
    retractall(pc0_simple(Run, _, _, _)),
    retractall(pc0_decision(Run, _, _)).

stable_word_ids(Words, IDs) :-
    empty_assoc(Empty), stable_word_ids(Words, 1, Empty, IDs).
stable_word_ids([], _, IDs, IDs).
stable_word_ids([[Answer|_]|Words], ID, IDs0, IDs) :-
    put_assoc(Answer, IDs0, ID, IDs1),
    ID1 is ID + 1,
    stable_word_ids(Words, ID1, IDs1, IDs).

word_id(Holder, Answer, ID) :-
    arg(32, Holder, IDs), get_assoc(Answer, IDs, ID).

event_leaf(none_holder, _) :- !.
event_leaf(Holder, Depth) :- holder_mode(Holder, Mode), event_leaf_mode(Mode, Holder, Depth).
event_leaf_mode(none, _, _).
event_leaf_mode(lean, Holder, Depth) :- update_depth(Holder, Depth).
event_leaf_mode(full, Holder, Depth) :- update_depth(Holder, Depth).

holder_mode(Holder, Mode) :- arg(1, Holder, Mode).

event_node(Holder, Depth, Grid, Words, Placed, GS, Parent,
           StateKey, KnownDead, NodeBefore) :-
    holder_mode(Holder, Mode),
    event_node_mode(Mode, Holder, Depth, Grid, Words, Placed, GS, Parent,
                    StateKey, KnownDead, NodeBefore).
event_node_mode(none, _, _, _, _, _, _, 0, none, false, 0).
event_node_mode(lean, Holder, Depth, _, _, _, _, Parent,
                none, false, NodeBefore) :-
    arg(4, Holder, NodeBefore),
    check_limit(Holder, nodes, 4),
    bump_arg(Holder, 4), update_depth(Holder, Depth), next_parent(Holder, Parent).
event_node_mode(full, Holder, Depth, Grid, Words, Placed, GS, Parent,
                StateKey, KnownDead, NodeBefore) :-
    event_node_mode(lean, Holder, Depth, Grid, Words, Placed, GS, Parent,
                    _LeanKey, _LeanDead, NodeBefore),
    state_sizes(Holder, Words, GS),
    canonical_state_key(Holder, Grid, Words, Placed, StateKey),
    state_observe(Holder, StateKey, KnownDead).

event_proof(Holder, Parent, WordID, Start, Dir) :-
    holder_mode(Holder, Mode),
    event_proof_mode(Mode, Holder, Parent-WordID-Start-Dir).
event_proof_mode(none, _, _).
event_proof_mode(lean, _, _).
event_proof_mode(full, Holder, Key) :-
    bump_arg(Holder, 21), remember_or_duplicate(Holder, proof, Key, 22).

event_child(Holder, Parent, WordID, Start, Dir, ParentRepeat0,
            ParentRepeat, ParentRepeatRoot, NodeBefore) :-
    holder_mode(Holder, Mode),
    event_child_mode(Mode, Holder, Parent-WordID-Start-Dir, ParentRepeat0,
                     ParentRepeat, ParentRepeatRoot, NodeBefore).
event_child_mode(none, _, _, ParentRepeat, ParentRepeat, false, 0).
event_child_mode(lean, Holder, Key, ParentRepeat, ParentRepeat, false,
                 NodeBefore) :-
    arg(4, Holder, NodeBefore),
    check_limit(Holder, decisions, 5), bump_arg(Holder, 5), bump_arg(Holder, 6),
    arg(5, Holder, Index), arg(19, Holder, Run),
    assertz(pc0_decision(Run, Index, Key)).
event_child_mode(full, Holder, Key, ParentRepeat0,
                 ParentRepeat, ParentRepeatRoot, NodeBefore) :-
    event_child_mode(lean, Holder, Key, ParentRepeat0, _, _, NodeBefore),
    remember_or_duplicate(Holder, child, Key, 15),
    ( simple_lookup(Holder, failed, Key)
    -> bump_arg(Holder, 25), KnownFailed = true
    ;  KnownFailed = false
    ),
    repeat_context(ParentRepeat0, KnownFailed, ParentRepeat, ParentRepeatRoot).

event_child_exhausted(Holder, Parent, WordID, Start, Dir,
                      ParentRepeatRoot, NodeBefore) :-
    holder_mode(Holder, Mode),
    event_child_exhausted_mode(Mode, Holder, Parent-WordID-Start-Dir,
                               ParentRepeatRoot, NodeBefore).
event_child_exhausted_mode(none, _, _, _, _).
event_child_exhausted_mode(lean, _, _, _, _).
event_child_exhausted_mode(full, Holder, Key, ParentRepeatRoot, NodeBefore) :-
    remember_key(Holder, failed, Key),
    ( ParentRepeatRoot == true
    -> arg(4, Holder, Nodes), Repeated is Nodes - NodeBefore,
       add_arg(Holder, 26, Repeated)
    ;  true
    ).

event_node_exhausted(Holder, StateKey, StateRepeatRoot, NodeBefore) :-
    holder_mode(Holder, Mode),
    event_node_exhausted_mode(Mode, Holder, StateKey, StateRepeatRoot, NodeBefore).
event_node_exhausted_mode(none, _, _, _, _).
event_node_exhausted_mode(lean, _, _, _, _).
event_node_exhausted_mode(full, Holder, StateKey, StateRepeatRoot, NodeBefore) :-
    exact_table_add(Holder, dead, StateKey, Added),
    ( Added == true -> bump_arg(Holder, 27) ; true ),
    ( StateRepeatRoot == true
    -> arg(4, Holder, Nodes), Repeated is Nodes - NodeBefore,
       add_arg(Holder, 24, Repeated)
    ;  true
    ).

repeat_context(Inside0, Known, Inside, Root) :-
    ( Known == true
    -> Inside = true, ( Inside0 == true -> Root = false ; Root = true )
    ;  Inside = Inside0, Root = false
    ).

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

remember_or_duplicate(Holder, Table, Key, DuplicateArg) :-
    ( simple_lookup(Holder, Table, Key) -> bump_arg(Holder, DuplicateArg)
    ; arg(19, Holder, Run), term_hash(Key, Hash),
      assertz(pc0_simple(Run, Table, Hash, Key))
    ).

remember_key(Holder, Table, Key) :-
    ( simple_lookup(Holder, Table, Key) -> true
    ; arg(19, Holder, Run), term_hash(Key, Hash),
      assertz(pc0_simple(Run, Table, Hash, Key))
    ).

simple_lookup(Holder, Table, Key) :-
    arg(19, Holder, Run), term_hash(Key, Hash),
    pc0_simple(Run, Table, Hash, Stored), Stored == Key, !.

canonical_state_key(Holder, Grid, Words, Placed,
                    state(Grid, Remaining, PlacedKeys)) :-
    maplist(entry_id(Holder), Words, Remaining0), msort(Remaining0, Remaining),
    maplist(placed_id_key(Holder), Placed, Placed0), msort(Placed0, PlacedKeys).
entry_id(Holder, [Answer|_], ID) :- word_id(Holder, Answer, ID).
placed_id_key(Holder, PW, p(ID,Start,Dir)) :-
    crosswordsmith_core:pw_answer(PW, Answer), word_id(Holder, Answer, ID),
    crosswordsmith_core:pw_start(PW, Start),
    crosswordsmith_core:pw_dir(PW, Dir).

state_observe(Holder, Key, KnownDead) :-
    exact_table_lookup(Holder, seen, Key, Seen),
    ( Seen == true -> bump_arg(Holder, 16)
    ; exact_table_add(Holder, seen, Key, _), bump_arg(Holder, 34)
    ),
    exact_table_lookup(Holder, dead, Key, KnownDead),
    ( KnownDead == true -> bump_arg(Holder, 23) ; true ).

exact_table_lookup(Holder, Table, Key, Found) :-
    state_fingerprint(Key, fp(H1,H2)), arg(19, Holder, Run),
    findall(Stored, pc0_exact(Run, Table, H1, H2, Stored), Bucket),
    ( Bucket = [_|_]
    -> ( exact_member(Key, Bucket)
       -> Found = true, bump_arg(Holder, 29)
       ;  Found = false, bump_arg(Holder, 28)
       )
    ; Found = false
    ).

exact_table_add(Holder, Table, Key, Added) :-
    exact_table_lookup(Holder, Table, Key, Found),
    ( Found == true -> Added = false
    ; state_fingerprint(Key, fp(H1,H2)), arg(19, Holder, Run),
      assertz(pc0_exact(Run, Table, H1, H2, Key)), Added = true
    ).

state_fingerprint(Key, fp(H1,H2)) :-
    ( ground(Key) -> true ; throw(error(non_ground_pc0_state(Key), _)) ),
    term_hash(Key, -1, 2147483647, H1),
    term_hash(pc0-Key, -1, 2147483629, H2).

exact_member(Key, [Stored|_]) :- Stored == Key, !.
exact_member(Key, [_|Stored]) :- exact_member(Key, Stored).

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
      crossing_proofs:null,duplicate_proofs:null,duplicate_children:null,
      duplicate_states:null,repeated_dead_entries:null,
      repeated_dead_nodes:null,parent_dedup_hits:null,
      parent_dedup_nodes:null,dead_states:null,canonical_states:null,
      exact_hash_hits:null,hash_collisions:null,decision_trace:null}).
holder_stats_mode(Mode, Holder, Stats) :-
    Mode \== none,
    arg(4,Holder,Nodes), arg(5,Holder,Decisions), arg(6,Holder,Places),
    arg(7,Holder,Unplaces), arg(8,Holder,Wipeouts), arg(9,Holder,Depth),
    arg(19,Holder,Run), findall(Key, pc0_decision(Run, _, Key), DecisionTrace),
    ( Mode == full
    -> arg(10,Holder,Entries), arg(11,Holder,Letters), arg(12,Holder,Bounds),
        arg(14,Holder,Transitions), arg(15,Holder,DupChildren),
       arg(16,Holder,DupStates), arg(21,Holder,Proofs),
       arg(22,Holder,DupProofs), arg(23,Holder,DeadEntries),
       arg(24,Holder,DeadNodes), arg(25,Holder,ParentHits),
       arg(26,Holder,ParentNodes), arg(27,Holder,DeadStates),
       arg(28,Holder,HashCollisions), arg(29,Holder,ExactHashHits),
       arg(34,Holder,CanonicalStates)
    ; Entries=null, Letters=null, Bounds=null, Transitions=null,
      Proofs=null, DupProofs=null, DupChildren=null, DupStates=null,
      DeadEntries=null, DeadNodes=null, ParentHits=null, ParentNodes=null,
      DeadStates=null, HashCollisions=null, ExactHashHits=null,
      CanonicalStates=null
    ),
    Stats = _{nodes:Nodes,decisions:Decisions,places:Places,
              unplaces:Unplaces,wipeouts:Wipeouts,max_depth:Depth,
              state_entries_max:Entries,letter_cells_max:Letters,
              boundary_cells_max:Bounds,support_transitions:Transitions,
              crossing_proofs:Proofs,duplicate_proofs:DupProofs,
              duplicate_children:DupChildren,duplicate_states:DupStates,
              repeated_dead_entries:DeadEntries,repeated_dead_nodes:DeadNodes,
              parent_dedup_hits:ParentHits,parent_dedup_nodes:ParentNodes,
              dead_states:DeadStates,canonical_states:CanonicalStates,
              exact_hash_hits:ExactHashHits,hash_collisions:HashCollisions,
              decision_trace:DecisionTrace}.

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
              crossing_proofs:Stats.crossing_proofs,
              duplicate_proofs:Stats.duplicate_proofs,
              duplicate_children:Stats.duplicate_children,
              duplicate_states:Stats.duplicate_states,
              repeated_dead_entries:Stats.repeated_dead_entries,
              repeated_dead_nodes:Stats.repeated_dead_nodes,
              parent_dedup_hits:Stats.parent_dedup_hits,
              parent_dedup_nodes:Stats.parent_dedup_nodes,
              dead_states:Stats.dead_states,
              canonical_states:Stats.canonical_states,
              exact_hash_hits:Stats.exact_hash_hits,
              hash_collisions:Stats.hash_collisions},
    Row = Base.

:- multifile prolog:error_message//1.
prolog:error_message(fixture_missing_clues(File)) -->
    ['probe fixture ~q does not define clues/1'-[File]].
prolog:error_message(fixture_word_count(File, Expected, Actual)) -->
    ['probe fixture ~q has ~d words; expected exactly ~d'-[File,Actual,Expected]].
