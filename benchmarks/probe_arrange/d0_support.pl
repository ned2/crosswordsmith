% P-D0 benchmark-only exact replay and support-delta observer.
% Product predicates are unchanged. Shadow residues are ordinary threaded
% terms; mutable holders contain counters only and are never read by search.
:- module(probe_arrange_d0,
          [ d0_corner/9,
            d0_operation/8
          ]).

:- use_module(library(apply), [foldl/4, include/3, maplist/3]).
:- use_module(library(assoc),
              [assoc_to_list/2, empty_assoc/1, get_assoc/3, put_assoc/4]).
:- use_module(library(lists),
              [append/3, intersection/3, list_to_set/2, member/2, nth1/3,
               reverse/2]).
:- use_module(library(ordsets),
              [ord_intersection/3, ord_subtract/3, ord_union/3]).
:- use_module(library(pairs), [map_list_to_pairs/3]).
:- use_module(library(statistics), [call_time/2]).

%!  d0_corner(+Words:list, +Grid:integer, +Corner:atom, +Seed,
%!             +Observe:boolean, -Result, -Summary:dict,
%!             -DecisionOrder:list, -Timing:dict) is det.
%
%   Replay one strict corner without a product inference limit. Observe=true
%   records P-D0 shadow evidence; false is the decision-order control.
d0_corner(Words, Grid, Corner, Seed, Observe, Result, Summary,
          DecisionOrder, Timing) :-
    init_holder(Holder),
    with_search_seed(Seed,
        ( crosswordsmith_core:reset_search_memos,
          call_time(d0_corner_goal(Words, Grid, Corner, Observe, Holder, Result),
                    T) )),
    holder_output(Observe, Holder, Summary, DecisionOrder),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

%!  d0_operation(+Words:list, +Grid:integer, +Seed, +Observe:boolean,
%!                -Result, -Summary:dict, -DecisionOrder:list,
%!                -Timing:dict) is det.
%
%   Replay the strict two-corner operation without a product inference limit.
% The strict operation uses exactly one reset and one seed stream across the
% two non-transpose representatives, as arrange_best_layout/6 does.
d0_operation(Words, Grid, Seed, Observe, Result, Summary, DecisionOrder,
             Timing) :-
    init_holder(Holder),
    with_search_seed(Seed,
        ( crosswordsmith_core:reset_search_memos,
          call_time(d0_corners([topleft_across,topright], Words, Grid,
                               Observe, Holder, Tagged), T) )),
    best_result(Tagged, Result),
    holder_output(Observe, Holder, Summary, DecisionOrder),
    Timing = _{wall:T.wall, cpu:T.cpu, inferences:T.inferences}.

with_search_seed(Seed, Goal) :-
    setup_call_cleanup(install_seed(Seed), Goal,
                       crosswordsmith_core:set_search_seed(-1)).

install_seed(null) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(none) :- !, crosswordsmith_core:set_search_seed(-1).
install_seed(Seed) :- crosswordsmith_core:set_search_seed(Seed).

d0_corners([], _, _, _, _, []).
d0_corners([Corner|Corners], Words, Grid, Observe, Holder, [Tagged|Rest]) :-
    d0_one(Words, Grid, Corner, Observe, Holder, Tagged),
    d0_corners(Corners, Words, Grid, Observe, Holder, Rest).

d0_corner_goal(Words, Grid, Corner, Observe, Holder, Result) :-
    d0_one(Words, Grid, Corner, Observe, Holder, Tagged),
    corner_result(Tagged, Result).

d0_one(Words, Grid, Corner, Observe, Holder, Tagged) :-
    crosswordsmith_core:init_gs(Grid, GS),
    crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
    ( once(d0_assign(Words, [], none, none, Grid, Start, Dir, GS, _, Placed,
                     Observe, Holder))
    -> crosswordsmith_arrange:layout_reward(5, 1, Placed, Reward),
       Tagged = ok(Reward, Placed)
    ;  Tagged = exhausted
    ).

corner_result(ok(Reward, Placed),
              result(placed, ok, false, null, Numbered, Reward)) :-
    crosswordsmith_core:assign_clue_numbers(Placed, Numbered).
corner_result(exhausted,
              result(infeasible, exhausted, false, null, [], null)).

best_result(Tagged, result(placed, ok, false, null, Numbered, Reward)) :-
    findall(R-P, member(ok(R, P), Tagged), OKs),
    OKs = [_|_], !,
    sort(1, @>=, OKs, [Reward-Best|_]),
    crosswordsmith_core:assign_clue_numbers(Best, Numbered).
best_result(_, result(infeasible, exhausted, false, null, [], null)).

% Exact local replay of assign_words_inc/9. PState is the unchanged product
% count state. RState is shadow-only last-recount evidence and never enters a
% branch decision.
d0_assign([], Placed, _PState, _RState, _, _, _, GS, GS, Placed, _, Holder) :-
    counter_add(Holder, leaves, 1).
d0_assign([W|Ws], Placed0, PState0, RState0, Grid, Start, Dir, GS0, GS,
           Placed, Observe, Holder) :-
    Words = [W|Ws],
    counter_add(Holder, nodes, 1),
    d0_select(Words, Placed0, PState0, RState0, Grid, Start, Dir, GS0,
              Entry, Remaining, PState, RState, Observe, Holder),
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    crosswordsmith_core:find_intersecting_word(
        Letters, Length, Placed0, Grid, Start, Dir),
    crosswordsmith_core:assign_word(
        Answer, Letters, Length, Start, Dir, Grid, GS0, PW, GS1),
    record_decision(Holder, Answer, Start, Dir),
    d0_assign(Remaining, [PW|Placed0], PState, RState, Grid, _S, _D,
              GS1, GS, Placed, Observe, Holder).

d0_select(Words, [], _PState0, _RState0, _Grid, _Start, _Dir, _GS,
          Entry, Remaining, none, none, _Observe, _Holder) :-
    crosswordsmith_core:seed_word_order(Words, Ordered),
    member(Entry, Ordered),
    crosswordsmith_core:remove_x(Entry, Words, Remaining).
d0_select(Words, [P|Ps], PState0, RState0, Grid, Start, Dir, GS,
          Entry, Remaining, PState, RState, Observe, Holder) :-
    Placed = [P|Ps],
    crosswordsmith_core:state_perturb(PState0, Perturb),
    crosswordsmith_core:inc_counts(
        PState0, Words, Placed, Grid, Start, Dir, GS, CountMap),
    structure_work(Observe, Holder, PState0, Words, CountMap),
    shadow_counts(Observe, PState0, RState0, Words, Placed, Grid, Start, Dir,
                  GS, CountMap, RState, Holder),
    map_list_to_pairs(crosswordsmith_core:count_of(CountMap), Words, Pairs),
    include(crosswordsmith_core:positive_key, Pairs, Placeable),
    length(Placeable, SortItems),
    counter_if(Observe, Holder, sort_items, SortItems),
    keysort(Placeable, Sorted),
    crosswordsmith_core:order_candidates(Perturb, Sorted, Ordered),
    member(Entry, Ordered),
    crosswordsmith_core:remove_x(Entry, Words, Remaining),
    crosswordsmith_core:entry_letters(Entry, EntryLetters),
    PState = state(CountMap, EntryLetters, Perturb).

structure_work(false, _, _, _, _).
structure_work(true, Holder, PState0, Words, CountMap) :-
    length(Words, N),
    counter_add(Holder, assoc_rebuilds, 1),
    counter_add(Holder, assoc_updates, N),
    counter_add(Holder, ordering_assoc_gets, N),
    counter_add(Holder, sort_calls, 1),
    ( PState0 == none
    -> counter_add(Holder, initial_full_counts, N)
    ;  PState0 = state(_Previous, LastLetters, _),
       sharing_partition(Words, LastLetters, Sharing, Carried),
       counter_add(Holder, refreshes, Sharing),
       counter_add(Holder, carry_assoc_gets, Carried)
    ),
    % Force the map argument to be observed here, not optimized away in this
    % benchmark twin; no value is changed.
    nonvar(CountMap).

sharing_partition([], _, 0, 0).
sharing_partition([Entry|Entries], LastLetters, Sharing, Carried) :-
    crosswordsmith_core:entry_letters(Entry, Letters),
    sharing_partition(Entries, LastLetters, S0, C0),
    ( crosswordsmith_core:shares_letter(Letters, LastLetters)
    -> Sharing is S0 + 1, Carried = C0
    ;  Sharing = S0, Carried is C0 + 1
    ).

shadow_counts(false, _PState, RState, _Words, _Placed, _Grid, _Start, _Dir,
              _GS, _CountMap, RState, _Holder).
shadow_counts(true, none, none, Words, Placed, Grid, _Start, _Dir, GS,
              CountMap, Residues, Holder) :-
    empty_assoc(R0),
    foldl(initial_residue(Placed, Grid, GS, CountMap, Holder),
          Words, R0, Residues).
shadow_counts(true, state(Previous, LastLetters, _), Residues0, Words,
              Placed, Grid, _Start, _Dir, GS, CountMap, Residues, Holder) :-
    empty_assoc(R0),
    placement_dirty(Placed, Grid, Dirty),
    foldl(refresh_or_carry(Previous, LastLetters, Residues0, Placed, Grid, GS,
                           CountMap, Dirty, Holder),
          Words, R0, Residues).

initial_residue(Placed, Grid, GS, CountMap, Holder, Entry, R0, R) :-
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    length(Letters, Length),
    scan_support(Letters, Length, Placed, Grid, GS,
                 Count, _Checks, RawProofs, _Unique),
    get_assoc(Answer, CountMap, ProductCount),
    require_equal_count(Answer, ProductCount, Count),
    maplist(add_proof_watches(Letters, Length, Grid, GS), RawProofs, Proofs),
    put_assoc(Answer, R0, residue(Count, Proofs), R),
    counter_add(Holder, initial_shadow_counts, 1).

refresh_or_carry(Previous, LastLetters, Residues0, Placed, Grid, GS,
                 CountMap, Dirty, Holder, Entry, R0, R) :-
    Entry = [Answer|_],
    crosswordsmith_core:entry_letters(Entry, Letters),
    ( crosswordsmith_core:shares_letter(Letters, LastLetters)
    -> get_assoc(Answer, Previous, PreviousBucket),
       get_assoc(Answer, Residues0, OldResidue),
       measure_refresh(Answer, Letters, PreviousBucket, OldResidue, Placed,
                       Grid, GS, CountMap, Dirty, Holder, NewResidue),
       put_assoc(Answer, R0, NewResidue, R)
    ;  get_assoc(Answer, Residues0, Residue),
       put_assoc(Answer, R0, Residue, R)
    ).

measure_refresh(Answer, Letters, PreviousBucket,
                residue(ResidueBucket, OldProofs), Placed,
                Grid, GS, CountMap, Dirty, Holder,
                residue(ExactBucket, Proofs)) :-
    Placed = [Newest|_],
    require_residue_bucket(Answer, PreviousBucket, ResidueBucket),
    length(Letters, Length),
    scan_support(Letters, Length, Placed, Grid, GS,
                 ExactBucket, ExactChecks, RawProofs, ExactUnique),
    get_assoc(Answer, CountMap, ProductCount),
    require_equal_count(Answer, ProductCount, ExactBucket),
    scan_support(Letters, Length, [Newest], Grid, GS,
                 NewestCount, NewestChecks, _NewestProofs, NewestUnique),
    surviving_residues(OldProofs, Letters, Grid, GS, Surviving),
    watch_dirty(OldProofs, Dirty, Holder),
    record_refresh(Holder, PreviousBucket, ExactBucket,
                   ExactChecks, ExactUnique, NewestCount, NewestChecks,
                   NewestUnique, OldProofs, Surviving),
    classify_refresh(PreviousBucket, OldProofs, Surviving,
                     NewestCount, NewestChecks, ExactBucket, ExactChecks,
                     Holder),
    maplist(add_proof_watches(Letters, Length, Grid, GS), RawProofs, Proofs).

require_equal_count(_Answer, Count, Count) :- !.
require_equal_count(Answer, Product, Shadow) :-
    throw(error(d0_count_mismatch(Answer, Product, Shadow), _)).

require_residue_bucket(_Answer, Bucket, Bucket) :- !.
require_residue_bucket(Answer, Visible, Residue) :-
    throw(error(d0_residue_bucket_mismatch(Answer, Visible, Residue), _)).

% The independent generator is deliberately not tabled and does not call the
% product pair_crossings table. It reproduces proof order and multiplicity
% without changing the product memo residue seen by later search nodes.
scan_support(Letters, Length, Placed, Grid, GS,
             Count, Checks, Proofs, Unique) :-
    Scan = scan(0, 0, []),
    ( \+ ( support_candidate(Letters, Length, Placed, Grid,
                              Start, Dir, RawProof),
             scan_bump(Scan, 1),
             crosswordsmith_core:check_word_fits(
                 Letters, Start, Dir, Grid, GS),
             scan_add_proof(Scan, RawProof),
             arg(2, Scan, N),
             N >= 2 )
    -> true
    ;  true
    ),
    arg(1, Scan, Checks),
    arg(2, Scan, Count),
    arg(3, Scan, Reversed),
    reverse(Reversed, Proofs),
    maplist(raw_geometry, Proofs, Geometries0),
    sort(Geometries0, Geometries),
    length(Geometries, Unique).

scan_bump(Scan, I) :-
    arg(I, Scan, N0), N is N0 + 1, nb_setarg(I, Scan, N).
scan_add_proof(Scan, Proof) :-
    scan_bump(Scan, 2),
    arg(3, Scan, Proofs), nb_setarg(3, Scan, [Proof|Proofs]).

support_candidate(Letters, Length, Placed, Grid, Start, Dir,
                  raw(SourceAnswer, SourceStart, SourceDir, Start, Dir)) :-
    member(PW, Placed),
    crosswordsmith_core:pw_answer(PW, SourceAnswer),
    crosswordsmith_core:pw_letters(PW, SourceLetters),
    crosswordsmith_core:pw_dir(PW, SourceDir),
    crosswordsmith_core:pw_start(PW, SourceStart),
    local_pair_crossings(Letters, SourceLetters, Crossings),
    member(x(SourcePos, Pos), Crossings),
    crosswordsmith_core:calc_num(
        SourceDir, Grid, SourcePos, SourceStart, CrossCell),
    crosswordsmith_core:swap_dir(SourceDir, Dir),
    crosswordsmith_core:calc_start(Dir, Grid, Pos, CrossCell, Start),
    crosswordsmith_core:fits_on_grid(Dir, Start, Length, Grid).

local_pair_crossings(Letters, SourceLetters, Crossings) :-
    findall(x(SourcePos, Pos),
            ( intersection(Letters, SourceLetters, Values),
              list_to_set(Values, Distinct),
              member(Value, Distinct),
              nth1(SourcePos, SourceLetters, Value),
              nth1(Pos, Letters, Value) ),
            Crossings).

raw_geometry(raw(_, _, _, Start, Dir), Start-Dir).

add_proof_watches(Letters, Length, Grid, GS,
                  raw(SA, SS, SD, Start, Dir),
                  proof(SA, SS, SD, Start, Dir, Watches)) :-
    proof_watches(Letters, Length, Start, Dir, Grid, GS, Watches).

proof_watches(Letters, Length, Start, Dir, Grid, gs(LGrid, _),
              watches(LetterWatches, BoundaryWatches, AdjacencyWatches)) :-
    run_cells(Start, Dir, Length, Grid, Cells, After),
    endpoint_watches(Start, After, Dir, Grid, Endpoints),
    append(Cells, Endpoints, Letter0),
    sort(Letter0, LetterWatches),
    sort(Cells, BoundaryWatches),
    adjacency_watches(Cells, Letters, Dir, Grid, LGrid, Adj0),
    sort(Adj0, AdjacencyWatches).

run_cells(Start, Dir, Length, Grid, Cells, After) :-
    run_cells_(Length, Start, Dir, Grid, Cells, After).
run_cells_(0, After, _, _, [], After).
run_cells_(N, Cell, Dir, Grid, [Cell|Cells], After) :-
    N > 0,
    crosswordsmith_core:next_cell(Dir, Cell, Grid, Next),
    N1 is N - 1,
    run_cells_(N1, Next, Dir, Grid, Cells, After).

endpoint_watches(Start, After, Dir, Grid, Endpoints) :-
    ( crosswordsmith_core:is_start_cell(Dir, Start, Grid)
    -> BeforeList = []
    ;  crosswordsmith_core:prev_cell(Dir, Start, Grid, Before),
       BeforeList = [Before]
    ),
    crosswordsmith_core:prev_cell(Dir, After, Grid, End),
    ( crosswordsmith_core:is_end_cell(Dir, End, Grid)
    -> AfterList = []
    ;  AfterList = [After]
    ),
    append(BeforeList, AfterList, Endpoints).

adjacency_watches([], [], _, _, _, []).
adjacency_watches([Cell|Cells], [_|Letters], Dir, Grid, LGrid, Watches) :-
    arg(Cell, LGrid, Value),
    ( var(Value) -> adjacent_cells(Dir, Cell, Grid, Here) ; Here = [] ),
    adjacency_watches(Cells, Letters, Dir, Grid, LGrid, Rest),
    append(Here, Rest, Watches).

adjacent_cells(down, Cell, Grid, Adjacent) :-
    Left is Cell - 1, Right is Cell + 1, M is Cell mod Grid,
    ( M == 0 -> Adjacent = [Left]
    ; M == 1 -> Adjacent = [Right]
    ; Adjacent = [Left,Right]
    ).
adjacent_cells(across, Cell, Grid, Adjacent) :-
    Above is Cell - Grid, Below is Cell + Grid, Last is Grid * Grid,
    ( Above =< 0 -> Adjacent = [Below]
    ; Below > Last -> Adjacent = [Above]
    ; Adjacent = [Above,Below]
    ).

surviving_residues([], _, _, _, 0).
surviving_residues([proof(_,_,_,Start,Dir,_)|Proofs], Letters, Grid, GS,
                    Count) :-
    surviving_residues(Proofs, Letters, Grid, GS, Count0),
    ( crosswordsmith_core:check_word_fits(Letters, Start, Dir, Grid, GS)
    -> Count is Count0 + 1
    ;  Count = Count0
    ).

placement_dirty([Newest|Older], Grid,
                dirty(LetterDirty, BoundaryDirty, LetterDirty)) :-
    crosswordsmith_core:pw_cells(Newest, NewCells0),
    sort(NewCells0, NewCells),
    placed_cells(Older, OldCells),
    ord_subtract(NewCells, OldCells, LetterDirty),
    word_boundaries(Newest, Grid, NewBoundaries),
    placed_boundaries(Older, Grid, OldBoundaries),
    ord_subtract(NewBoundaries, OldBoundaries, BoundaryDirty).

placed_cells(Placed, Cells) :-
    findall(Cell,
            ( member(PW, Placed),
              crosswordsmith_core:pw_cells(PW, PWCells),
              member(Cell, PWCells) ),
            Cells0),
    sort(Cells0, Cells).

placed_boundaries(Placed, Grid, Boundaries) :-
    foldl(add_word_boundaries(Grid), Placed, [], Boundaries).

add_word_boundaries(Grid, PW, B0, B) :-
    word_boundaries(PW, Grid, WB),
    ord_union(B0, WB, B).

word_boundaries(PW, Grid, Boundaries) :-
    crosswordsmith_core:pw_start(PW, Start),
    crosswordsmith_core:pw_end(PW, End),
    crosswordsmith_core:pw_dir(PW, Dir),
    ( crosswordsmith_core:is_start_cell(Dir, Start, Grid)
    -> Before = []
    ;  crosswordsmith_core:prev_cell(Dir, Start, Grid, Cell0), Before = [Cell0]
    ),
    ( crosswordsmith_core:is_end_cell(Dir, End, Grid)
    -> After = []
    ;  crosswordsmith_core:next_cell(Dir, End, Grid, Cell1), After = [Cell1]
    ),
    append(Before, After, B0),
    sort(B0, Boundaries).

watch_dirty([], _, _).
watch_dirty([proof(_,_,_,_,_,watches(Letter,Boundary,Adjacency))|Proofs],
            dirty(LetterDirty,BoundaryDirty,AdjacencyDirty), Holder) :-
    record_watch(Holder, letter, Letter, LetterDirty),
    record_watch(Holder, boundary, Boundary, BoundaryDirty),
    record_watch(Holder, adjacency, Adjacency, AdjacencyDirty),
    counter_add(Holder, watched_residues, 1),
    watch_dirty(Proofs, dirty(LetterDirty,BoundaryDirty,AdjacencyDirty), Holder).

record_watch(Holder, Type, Watches, DirtyCells) :-
    length(Watches, N),
    ord_intersection(Watches, DirtyCells, Dirty),
    length(Dirty, ND),
    watch_key(Type, watched, WatchKey),
    watch_key(Type, dirty, DirtyKey),
    watch_key(Type, events, EventsKey),
    watch_key(Type, dirty_events, DirtyEventsKey),
    counter_add(Holder, WatchKey, N),
    counter_add(Holder, DirtyKey, ND),
    counter_add(Holder, EventsKey, 1),
    ( ND > 0 -> counter_add(Holder, DirtyEventsKey, 1) ; true ).

watch_key(Type, Suffix, Key) :- atomic_list_concat([watch,Type,Suffix], '_', Key).

record_refresh(Holder, Previous, Exact, ExactChecks, ExactUnique,
               NewestCount, NewestChecks, NewestUnique, OldProofs, Surviving) :-
    transition_add(Holder, Previous, Exact),
    counter_add(Holder, exact_candidate_checks, ExactChecks),
    counter_add(Holder, newest_candidate_checks, NewestChecks),
    counter_add(Holder, exact_proofs, Exact),
    counter_add(Holder, exact_unique_geometries, ExactUnique),
    counter_add(Holder, newest_proofs, NewestCount),
    counter_add(Holder, newest_unique_geometries, NewestUnique),
    candidate_bucket_add(Holder, Exact, ExactChecks),
    ( Exact =\= ExactUnique
    -> counter_add(Holder, proof_geometry_divergences, 1)
    ;  true
    ),
    ( NewestCount =\= NewestUnique
    -> counter_add(Holder, newest_proof_geometry_divergences, 1)
    ;  true
    ),
    length(OldProofs, Residues),
    counter_add(Holder, residue_proofs, Residues),
    counter_add(Holder, residue_survivors, Surviving),
    ( Residues > 0 -> counter_add(Holder, residue_events, 1) ; true ),
    ( Exact =:= 2 -> counter_add(Holder, saturated_refreshes, 1) ; true ).

classify_refresh(0, _OldProofs, _Surviving, NewestCount, NewestChecks,
                 Exact, _ExactChecks, Holder) :- !,
    verify_classification(0, NewestCount, Exact, Holder),
    counter_add(Holder, proposed_candidate_checks, NewestChecks).
classify_refresh(1, OldProofs, Surviving, NewestCount, NewestChecks,
                 Exact, _ExactChecks, Holder) :- !,
    require_single_residue(OldProofs),
    ( NewestCount =:= 2
    -> Predicted = 2, ProposedChecks = NewestChecks
    ;  Predicted is min(2, NewestCount + Surviving),
       ProposedChecks is NewestChecks + 1
    ),
    verify_classification(1, Predicted, Exact, Holder),
    counter_add(Holder, proposed_candidate_checks, ProposedChecks).
classify_refresh(2, _OldProofs, _Surviving, _NewestCount, _NewestChecks,
                 _Exact, ExactChecks, Holder) :-
    % Bucket 2 has unknown unobserved older proofs once its two residues fail.
    % A-D2 therefore admits no precheck and retains the current full recount.
    counter_add(Holder, proposed_candidate_checks, ExactChecks),
    counter_add(Holder, fallback_refreshes, 1).

require_single_residue([_]) :- !.
require_single_residue(Proofs) :-
    length(Proofs, Count),
    throw(error(d0_bucket_one_residue_count(Count), _)).

verify_classification(Previous, Predicted, Exact, Holder) :-
    ( Predicted =:= Exact
    -> counter_add(Holder, classified_refreshes, 1),
       counter_add(Holder, verified_classifications, 1),
       classification_add(Holder, Previous)
    ;  throw(error(d0_classification_mismatch(Previous, Predicted, Exact), _))
    ).

candidate_bucket_add(Holder, Bucket, Checks) :-
    format(atom(EventsKey), 'bucket_~d_events', [Bucket]),
    format(atom(ChecksKey), 'bucket_~d_candidate_checks', [Bucket]),
    counter_add(Holder, EventsKey, 1),
    counter_add(Holder, ChecksKey, Checks).

% Holder: counters assoc, transition assoc, classification assoc, reverse
% chronological decision list. Counters are observational only.
init_holder(holder(Counters, Transitions, Classifications, [])) :-
    empty_assoc(Counters), empty_assoc(Transitions), empty_assoc(Classifications).

counter_if(false, _, _, _).
counter_if(true, Holder, Key, Add) :- counter_add(Holder, Key, Add).

counter_add(Holder, Key, Add) :-
    arg(1, Holder, C0),
    ( get_assoc(Key, C0, N0) -> true ; N0 = 0 ),
    N is N0 + Add,
    put_assoc(Key, C0, N, C),
    nb_setarg(1, Holder, C).

counter_value(Holder, Key, Value) :-
    arg(1, Holder, Counters),
    ( get_assoc(Key, Counters, Value) -> true ; Value = 0 ).

transition_add(Holder, Previous, Exact) :-
    arg(2, Holder, T0), Key = t(Previous, Exact),
    ( get_assoc(Key, T0, N0) -> true ; N0 = 0 ),
    N is N0 + 1, put_assoc(Key, T0, N, T), nb_setarg(2, Holder, T).

classification_add(Holder, Previous) :-
    arg(3, Holder, C0),
    ( get_assoc(Previous, C0, N0) -> true ; N0 = 0 ),
    N is N0 + 1, put_assoc(Previous, C0, N, C), nb_setarg(3, Holder, C).

record_decision(Holder, Answer, Start, Dir) :-
    counter_add(Holder, decisions, 1),
    arg(4, Holder, D0), nb_setarg(4, Holder, [d(Answer,Start,Dir)|D0]).

holder_output(Observe, Holder, Summary, DecisionOrder) :-
    arg(4, Holder, Reversed), reverse(Reversed, DecisionOrder),
    counter_value(Holder, nodes, Nodes),
    counter_value(Holder, decisions, Decisions),
    ( Observe == true
    -> observed_summary(Holder, Nodes, Decisions, Summary)
    ;  Summary = _{observed:false,nodes:Nodes,decisions:Decisions}
    ).

observed_summary(Holder, Nodes, Decisions, Summary) :-
    counter_value(Holder, refreshes, Refreshes),
    counter_value(Holder, classified_refreshes, Classified),
    counter_value(Holder, verified_classifications, Verified),
    counter_value(Holder, fallback_refreshes, Fallback),
    counter_value(Holder, exact_candidate_checks, ExactChecks),
    counter_value(Holder, proposed_candidate_checks, ProposedChecks),
    counter_value(Holder, newest_candidate_checks, NewestChecks),
    counter_value(Holder, exact_proofs, ExactProofs),
    counter_value(Holder, exact_unique_geometries, ExactUnique),
    counter_value(Holder, newest_proofs, NewestProofs),
    counter_value(Holder, newest_unique_geometries, NewestUnique),
    counter_value(Holder, proof_geometry_divergences, Divergences),
    counter_value(Holder, newest_proof_geometry_divergences, NewDivergences),
    counter_value(Holder, residue_events, ResidueEvents),
    counter_value(Holder, residue_proofs, ResidueProofs),
    counter_value(Holder, residue_survivors, ResidueSurvivors),
    counter_value(Holder, saturated_refreshes, Saturated),
    counter_value(Holder, assoc_rebuilds, AssocRebuilds),
    counter_value(Holder, assoc_updates, AssocUpdates),
    counter_value(Holder, ordering_assoc_gets, OrderingGets),
    counter_value(Holder, carry_assoc_gets, CarryGets),
    counter_value(Holder, sort_calls, SortCalls),
    counter_value(Holder, sort_items, SortItems),
    counter_value(Holder, initial_full_counts, InitialFull),
    counter_value(Holder, initial_shadow_counts, InitialShadow),
    bucket_summaries(Holder, Buckets),
    watch_summary(Holder, Watch),
    transition_summary(Holder, Transitions),
    classification_summary(Holder, Classifications),
    Summary = _{observed:true,nodes:Nodes,decisions:Decisions,
                refreshes:Refreshes,classified_refreshes:Classified,
                verified_classifications:Verified,
                fallback_refreshes:Fallback,
                exact_candidate_checks:ExactChecks,
                proposed_candidate_checks:ProposedChecks,
                newest_candidate_checks:NewestChecks,
                exact_proofs:ExactProofs,
                exact_unique_geometries:ExactUnique,
                newest_proofs:NewestProofs,
                newest_unique_geometries:NewestUnique,
                proof_geometry_divergences:Divergences,
                newest_proof_geometry_divergences:NewDivergences,
                residue_events:ResidueEvents,residue_proofs:ResidueProofs,
                residue_survivors:ResidueSurvivors,
                saturated_refreshes:Saturated,
                assoc_rebuilds:AssocRebuilds,assoc_updates:AssocUpdates,
                ordering_assoc_gets:OrderingGets,carry_assoc_gets:CarryGets,
                sort_calls:SortCalls,sort_items:SortItems,
                initial_full_counts:InitialFull,
                initial_shadow_counts:InitialShadow,
                buckets:Buckets,watch:Watch,transitions:Transitions,
                classifications:Classifications}.

bucket_summaries(Holder, Buckets) :-
    maplist(bucket_summary(Holder), [0,1,2], Buckets).
bucket_summary(Holder, Bucket,
               _{bucket:Bucket,events:Events,candidate_checks:Checks}) :-
    format(atom(EventsKey), 'bucket_~d_events', [Bucket]),
    format(atom(ChecksKey), 'bucket_~d_candidate_checks', [Bucket]),
    counter_value(Holder, EventsKey, Events),
    counter_value(Holder, ChecksKey, Checks).

watch_summary(Holder, Watch) :-
    maplist(watch_type_summary(Holder), [letter,boundary,adjacency], Watch).
watch_type_summary(Holder, Type,
                   _{type:Type,watches:Watches,dirty:Dirty,
                     events:Events,dirty_events:DirtyEvents}) :-
    watch_key(Type, watched, WatchKey), watch_key(Type, dirty, DirtyKey),
    watch_key(Type, events, EventsKey),
    watch_key(Type, dirty_events, DirtyEventsKey),
    counter_value(Holder, WatchKey, Watches),
    counter_value(Holder, DirtyKey, Dirty),
    counter_value(Holder, EventsKey, Events),
    counter_value(Holder, DirtyEventsKey, DirtyEvents).

transition_summary(Holder, Summary) :-
    arg(2, Holder, Assoc), assoc_to_list(Assoc, Pairs),
    maplist(transition_row, Pairs, Summary).
transition_row(t(Previous,Exact)-Count,
               _{previous:Previous,exact:Exact,count:Count}).

classification_summary(Holder, Summary) :-
    arg(3, Holder, Assoc), assoc_to_list(Assoc, Pairs),
    maplist(classification_row, Pairs, Summary).
classification_row(Previous-Count, _{previous:Previous,count:Count}).

:- multifile prolog:error_message//1.
prolog:error_message(d0_count_mismatch(Answer, Product, Shadow)) -->
    ['P-D0 exact recount for ~q returned ~d; product returned ~d'-
     [Answer,Shadow,Product]].
prolog:error_message(d0_residue_bucket_mismatch(Answer, Visible, Residue)) -->
    ['P-D0 residue for ~q has bucket ~d; visible bucket is ~d'-
     [Answer,Residue,Visible]].
prolog:error_message(d0_bucket_one_residue_count(Count)) -->
    ['P-D0 visible bucket 1 has ~d stored residues; expected exactly 1'-[Count]].
prolog:error_message(d0_classification_mismatch(Previous, Predicted, Exact)) -->
    ['P-D0 bucket ~d classifier predicted ~d; exact recount returned ~d'-
     [Previous,Predicted,Exact]].
