% B0-I measurement rig for the shipped MAC + dom/wdeg fill policy.
%
% This is benchmark-only code.  It calls unchanged fill internals for setup,
% selection, candidate enumeration, support masks, and layout conversion.  The
% local twins exist only where instrumentation must observe a hot predicate.
:- module(probe_mac_b0,
          [ compare_product/5,
            run_row/9,
            run_profile/7,
            run_cli/0
          ]).

:- use_module(library(assoc)).
:- use_module(library(json)).
:- use_module(library(lists)).
:- use_module(library(pairs)).
:- use_module(library(prolog_profile)).
:- use_module(library(sha)).
:- use_module(library(time)).

% compare_product(+Grid,+Dict,+MinScore,+Seeds,+Budget)
% Replay the product and the counter-free twin from fresh grids and require the
% same outcome, complete layout, input words, and final-attempt node count.
% Product internals expose the final-attempt node counter but not the attempt
% number, so the twin's attempt count is reported rather than compared.
compare_product(Grid, Dict, MinScore, Seeds, Budget) :-
    crosswordsmith_fill:load_dict(Dict, [min_score(MinScore)], DBL, Index, _),
    fresh_slots(Grid, Seeds, ProductSearch, ProductAll),
    crosswordsmith_fill:fill_attempt(ProductSearch, ProductAll, DBL, Index,
                                     Budget, ProductOutcome, ProductNumbered,
                                     ProductWords),
    ( nb_current(cw_mac_attempt_nodes, ProductNodes) -> true ; ProductNodes = 0 ),
    fresh_slots(Grid, Seeds, TwinSearch, TwinAll),
    run_loaded(TwinSearch, TwinAll, DBL, Index, Budget, 0, off,
               TwinOutcome, TwinNumbered, TwinWords, TwinReport),
    TwinAttempts = TwinReport.attempts_count,
    TwinNodes = TwinReport.last_attempt_nodes,
    (   ProductOutcome == TwinOutcome,
        ProductNumbered =@= TwinNumbered,
        ProductWords =@= TwinWords,
        ProductNodes =:= TwinNodes
    ->  Verdict = pass
    ;   Verdict = fail
    ),
    truth_value(ProductNumbered =@= TwinNumbered, LayoutEqual),
    truth_value(ProductWords =@= TwinWords, WordsEqual),
    Result = _{verdict:Verdict, outcome:ProductOutcome,
               product_nodes:ProductNodes, twin_nodes:TwinNodes,
               twin_attempts:TwinAttempts,
               layout_equal:LayoutEqual, words_equal:WordsEqual},
    json_write_dict(current_output, Result, [width(0)]), nl,
    Verdict == pass.

% run_row(+Label,+Grid,+Dict,+MinScore,+Seeds,+Budget,+Timeout,+Mode,-Report)
% Mode is off, counters, or timing. Timeout=0 disables the wall timeout.
run_row(Label, Grid, Dict, MinScore, Seeds, Budget, Timeout, Mode, Report) :-
    format(user_error, 'b0: row ~w setup~n', [Label]),
    get_time(T0),
    crosswordsmith_fill:load_dict(Dict, [min_score(MinScore)], DBL, Index, Scores),
    fresh_slots(Grid, Seeds, Search, All),
    get_time(T1),
    SetupSeconds is T1-T0,
    run_loaded(Search, All, DBL, Index, Budget, Timeout, Mode,
               Outcome, Numbered, _InputWords, Core),
    get_time(T2),
    TotalSeconds is T2-T0,
    SetupTotalSeconds is SetupSeconds+Core.mac_setup_seconds,
    score_layout(Numbered, Scores, DBL, Quality),
    layout_digest(Numbered, Digest),
    Report = Core.put(_{label:Label, grid:Grid, dict:Dict,
                        min_score:MinScore, budget:Budget, timeout:Timeout,
                        mode:Mode, outcome:Outcome,
                        input_setup_seconds:SetupSeconds,
                        setup_seconds:SetupTotalSeconds,
                        total_seconds:TotalSeconds,
                        fill_digest:Digest, quality:Quality}),
    format(user_error, 'b0: row ~w ~w after ~2f s~n',
           [Label, Outcome, TotalSeconds]).

% Statistical profiling is deliberately separate from decision counters.
% Wrapper cumulative ticks attribute selection/candidate, placement, and
% propagation/support without placing a clock inside every revision.
run_profile(Label, Grid, Dict, MinScore, Seeds, Budget, Report) :-
    reset_profiler,
    with_output_to(string(_ProfileText),
        profile(run_row(Label, Grid, Dict, MinScore, Seeds, Budget, 0, off, Core),
                [time(wall), sample_rate(1000), ports(false), top(0)])),
    profile_data(Data),
    phase_profile(Data, Profile),
    Report = Core.put(profile, Profile).

phase_profile(Data, Profile) :-
    profile_ticks(Data.nodes, probe_mac_b0:b0_select_candidate/9, Selection),
    profile_ticks(Data.nodes, probe_mac_b0:b0_place_phase/9, Placement),
    profile_ticks(Data.nodes, probe_mac_b0:b0_propagate_phase/7, Propagation),
    TotalTicks = Data.summary.ticks,
    Profile = _{total_ticks:TotalTicks,
                selection_candidate_ticks:Selection,
                placement_ticks:Placement,
                propagation_support_ticks:Propagation,
                sample_period_us:Data.summary.sample_period,
                accounting_ticks:Data.summary.accounting}.

profile_ticks(Nodes, PI, Ticks) :-
    ( member(N, Nodes), N.predicate == PI
    -> Ticks is N.ticks_self + N.ticks_siblings
    ;  Ticks = 0
    ).

fresh_slots(Grid, none, Search, All) :- !,
    crosswordsmith_fill:fill_grid(Grid, _, All, _),
    Search = All.
fresh_slots(Grid, Seeds, Search, All) :-
    crosswordsmith_fill:fill_grid(Grid, Size, All, _),
    crosswordsmith_fill:apply_seeds(Seeds, Size, All, Keys),
    exclude(crosswordsmith_fill:seeded_slot(Keys), All, Search).

run_loaded(Search, All, DBL, Index, Budget, Timeout, Mode,
           Outcome, Numbered, InputWords, Report) :-
    % Match fill_attempt/8 exactly: seed_used/3 is outside the inference
    % budget, while every operation corresponding to fill_search/5 is inside.
    crosswordsmith_fill:seed_used(All, Search, Used),
    init_run_state(Mode),
    run_limited(Timeout, Budget,
                b0_fill_search(Search, DBL, Index, Used, State, Assign), Stop),
    run_times(MacSetupSeconds, SearchSeconds),
    finish_active_attempt(Stop),
    map_stop(Stop, Outcome0),
    (   Outcome0 == filled
    ->  crosswordsmith_fill:mac_bind_fill(Assign, State.id_slots),
        crosswordsmith_fill:slots_to_layout(All, Numbered, InputWords),
        Outcome = filled
    ;   Outcome = Outcome0, Numbered = [], InputWords = []
    ),
    core_report(SearchSeconds, Report0),
    Report=Report0.put(mac_setup_seconds,MacSetupSeconds).

init_run_state(Mode) :-
    nb_setval(cw_b0_mode, Mode),
    nb_setval(cw_b0_setup_complete, false),
    nb_setval(cw_b0_mac_setup_seconds, 0.0),
    nb_setval(cw_b0_search_started, 0.0),
    nb_setval(cw_b0_attempt_count, 0),
    nb_setval(cw_b0_attempt_records, []),
    nb_setval(cw_b0_attempt_active, none),
    nb_setval(cw_b0_root_record, _{}).

% Local equivalent of fill_search/5's non-empty arm. Setup timing is itself
% inside the limited goal, and all report state needed after an inference/time
% interrupt is persisted before root propagation begins.
b0_fill_search(Search, DBL, Index, Used, State, Assign) :-
    get_time(TSetup0),
    setup_core(Search, DBL, Index, Used, State),
    get_time(TSetup1),
    MacSetupSeconds is TSetup1-TSetup0,
    nb_setval(cw_b0_mac_setup_seconds, MacSetupSeconds),
    nb_setval(cw_b0_edge_a, State.edge_a),
    nb_setval(cw_b0_setup_complete, true),
    nb_setval(cw_b0_search_started, TSetup1),
    b0_core(State, Assign).

run_times(MacSetupSeconds, SearchSeconds) :-
    nb_getval(cw_b0_mac_setup_seconds, MacSetupSeconds),
    (   nb_getval(cw_b0_setup_complete, true)
    ->  nb_getval(cw_b0_search_started, SearchStart),
        get_time(Now), SearchSeconds is Now-SearchStart
    ;   SearchSeconds = 0.0
    ).

run_limited(0, Budget, Goal, Stop) :- !,
    call_with_inference_limit((Goal -> R=ok ; R=exhausted), Budget, Limit),
    ( Limit == inference_limit_exceeded -> Stop=budget ; Stop=R ).
run_limited(Timeout, Budget, Goal, Stop) :-
    catch(call_with_time_limit(
              Timeout,
              ( call_with_inference_limit((Goal -> R=ok ; R=exhausted),
                                          Budget, Limit),
                ( Limit == inference_limit_exceeded -> Stop=budget ; Stop=R ) )),
          time_limit_exceeded, Stop=timeout).

map_stop(ok, filled).
map_stop(exhausted, infeasible).
map_stop(budget, not_proven).
map_stop(timeout, timeout).

setup_core(Slots, DBL, Index, Used, State) :-
    crosswordsmith_fill:build_masks(Index, MA),
    crosswordsmith_fill:mac_flat_masks(MA, LmA),
    crosswordsmith_fill:mac_buckets(DBL, BucketA),
    length(Slots, N), NMax is N-1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, Slots),
    maplist(crosswordsmith_fill:mac_slot_len, Slots, Lens),
    pairs_keys_values(IdLens, Ids, Lens), list_to_assoc(IdLens, LenA),
    crosswordsmith_fill:mac_edges(IdSlots, EdgeA, NEdges),
    crosswordsmith_fill:mac_init_weights(NEdges),
    maplist(crosswordsmith_fill:mac_init_domain(DBL, LmA), IdSlots, IdDoms),
    list_to_assoc(IdDoms, DomA0),
    init_probe(N, NEdges),
    State = state{ids:Ids,id_slots:IdSlots,len_a:LenA,bucket_a:BucketA,
                   edge_a:EdgeA,lm_a:LmA,dom_a:DomA0,used:Used}.

init_probe(NSlots, NEdges) :-
    nb_setval(cw_b0_attempt_count, 0),
    nb_setval(cw_b0_attempt_records, []),
    nb_setval(cw_b0_attempt_active, none),
    nb_setval(cw_b0_root_record, _{}),
    C = c(0,0,0,0,0,0,0), nb_setval(cw_b0_counts, C),
    T = t(0.0,0.0,0.0), nb_setval(cw_b0_times, T),
    functor(Q, q, NSlots), zero_args(Q, NSlots), nb_setval(cw_b0_queue, Q),
    functor(B, b, NEdges), zero_args(B, NEdges), nb_setval(cw_b0_bumps, B),
    functor(P, p, NEdges), ones_args(P, NEdges), nb_setval(cw_b0_peaks, P).

zero_args(_, 0) :- !.
zero_args(T, I) :- nb_setarg(I, T, 0), I1 is I-1, zero_args(T, I1).
ones_args(_, 0) :- !.
ones_args(T, I) :- nb_setarg(I, T, 1.0), I1 is I-1, ones_args(T, I1).

b0_core(State, Assign) :-
    \+ ( gen_assoc(_, State.dom_a, D), D =:= 0 ),
    counter_snapshot(RootStart),
    (   b0_propagate_phase(State.ids, State.dom_a, State.len_a, State.edge_a,
                           State.lm_a, State.ids, DomA1)
    ->  counter_delta(RootStart, Root),
        nb_setval(cw_b0_root_record, Root),
        ( crosswordsmith_fill:current_search_seed(_)
        -> SearchMode=seeded
        ;  SearchMode=default ),
        b0_restart_loop(1, 500, State.ids, DomA1, State.len_a,
                        State.bucket_a, State.edge_a, State.lm_a, SearchMode,
                        State.used, Assign)
    ;   counter_delta(RootStart, Root),
        nb_setval(cw_b0_root_record, Root),
        fail
    ).

b0_restart_loop(Attempt, Cap, Ids, DomA, LenA, BucketA, EdgeA, LmA, Mode,
                Used, Assign) :-
    crosswordsmith_fill:mac_attempt_pick(Mode, Attempt, Pick),
    nb_setval(cw_mac_attempt_nodes, 0),
    nb_setval(cw_b0_attempt_count, Attempt),
    counter_snapshot(Start),
    bump_vector(BumpStart),
    weight_vector(AttemptPeak0), list_compound(ap,AttemptPeak0,AttemptPeaks),
    nb_setval(cw_b0_attempt_peaks,AttemptPeaks),
    nb_setval(cw_b0_attempt_active, active(Attempt,Cap,Start,BumpStart)),
    format(user_error, 'b0: attempt ~w cap ~w~n', [Attempt, Cap]),
    catch(( b0_search(Ids, DomA, LenA, BucketA, EdgeA, LmA, Pick, Cap,
                      Used, Assign)
          -> AttemptStop=filled
          ;  AttemptStop=exhausted ),
          mac_cap, AttemptStop=cap),
    record_attempt(AttemptStop),
    (   AttemptStop == cap
    ->  crosswordsmith_fill:mac_age_weights(0.99),
        Cap1 is (Cap*3+1)//2, A1 is Attempt+1,
        b0_restart_loop(A1, Cap1, Ids, DomA, LenA, BucketA, EdgeA, LmA,
                        Mode, Used, Assign)
    ;   AttemptStop == filled
    ).

b0_search([], _DomA, _LenA, _BucketA, _EdgeA, _LmA, _Pick, _Cap, _Used, []).
b0_search(Unfilled, DomA, LenA, BucketA, EdgeA, LmA, Pick, Cap, Used,
          [Best-Word|Assign]) :-
    Unfilled=[_|_],
    b0_select_candidate(Unfilled, DomA, LenA, BucketA, EdgeA, Pick,
                        Best, Rest, Word),
    \+ memberchk(Word, Used),
    nb_getval(cw_mac_attempt_nodes, N0), N1 is N0+1,
    nb_setval(cw_mac_attempt_nodes, N1), count_node,
    ( N1 > Cap -> throw(mac_cap) ; true ),
    ( N1 mod 50000 =:= 0 -> format(user_error, 'b0: heartbeat node ~w~n',[N1]) ; true ),
    crosswordsmith_fill:mac_edges_of(EdgeA, Best, Es),
    b0_place_phase(Es, Word, Rest, LenA, LmA, DomA, DomA1, [], Dirty),
    b0_propagate_phase(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2),
    b0_search(Rest, DomA2, LenA, BucketA, EdgeA, LmA, Pick, Cap,
              [Word|Used], Assign).

% Wrapper kept non-recursive so profiler cumulative ticks form one phase.
b0_select_candidate(Unfilled, DomA, LenA, BucketA, EdgeA, Pick,
                    Best, Rest, Word) :-
    phase_start(T0),
    crosswordsmith_fill:mac_select_dwd(Unfilled, DomA, EdgeA, Best),
    selectchk(Best, Unfilled, Rest),
    get_assoc(Best, DomA, Dom), get_assoc(Best, LenA, Len),
    get_assoc(Len, BucketA, Bucket),
    % End before the nondet lazy enumerator: on redo Prolog resumes inside
    % mac_candidate/4, so timing after it would repeatedly include descendant
    % search time. Its self-time remains covered by the profiler cross-check.
    phase_end(1, T0),
    crosswordsmith_fill:mac_candidate(Pick, Dom, Bucket, Word).

b0_place_phase(Es, Word, Rest, LenA, LmA, DomA0, DomA, Dirty0, Dirty) :-
    phase_start(T0),
    (   b0_place(Es, Word, Rest, LenA, LmA, DomA0, DomA, Dirty0, Dirty)
    ->  phase_end(2,T0)
    ;   phase_end(2,T0), fail
    ).

b0_place([], _Word, _Unfilled, _LenA, _LmA, DomA, DomA, Dirty, Dirty).
b0_place([e(Pos,T,TPos,EId)|Es], Word, Unfilled, LenA, LmA, DomA0, DomA,
         Dirty0, Dirty) :-
    (   memberchk(T, Unfilled)
    ->  nth0(Pos, Word, Ch), get_assoc(T, LenA, TLen),
        crosswordsmith_fill:mac_lm_of(LmA, TLen, TPos, LmT),
        crosswordsmith_fill:mac_letter_arg(Ch, ChI), arg(ChI, LmT, ChMask),
        get_assoc(T, DomA0, TDom), TDom1 is TDom /\ ChMask,
        ( TDom1 =:= 0 -> count_place_dwo, b0_bump_edge(EId), fail ; true ),
        ( TDom1 =:= TDom
        -> DomA1=DomA0, Dirty1=Dirty0
        ;  put_assoc(T, DomA0, TDom1, DomA1), Dirty1=[T|Dirty0] )
    ;   DomA1=DomA0, Dirty1=Dirty0
    ),
    b0_place(Es, Word, Unfilled, LenA, LmA, DomA1, DomA, Dirty1, Dirty).

b0_propagate_phase(Queue, DomA0, LenA, EdgeA, LmA, Unfilled, DomA) :-
    phase_start(T0),
    (   b0_propagate(Queue, DomA0, LenA, EdgeA, LmA, Unfilled, DomA)
    ->  phase_end(3,T0)
    ;   phase_end(3,T0), fail
    ).

b0_propagate([], DomA, _LenA, _EdgeA, _LmA, _Unfilled, DomA).
b0_propagate([S|Queue], DomA0, LenA, EdgeA, LmA, Unfilled, DomA) :-
    count_queue([S|Queue]),
    crosswordsmith_fill:mac_edges_of(EdgeA, S, Es),
    get_assoc(S, DomA0, SDom), get_assoc(S, LenA, SLen),
    b0_revise(Es, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA1,
              Queue, Queue1),
    b0_propagate(Queue1, DomA1, LenA, EdgeA, LmA, Unfilled, DomA).

b0_revise([], _SDom, _SLen, _LenA, _LmA, _Unfilled, DomA, DomA, Q, Q).
b0_revise([e(Pos,T,TPos,EId)|Es], SDom, SLen, LenA, LmA, Unfilled,
          DomA0, DomA, Q0, Q) :-
    (   memberchk(T, Unfilled)
    ->  get_assoc(T, LenA, TLen),
        crosswordsmith_fill:mac_lm_of(LmA, SLen, Pos, LmS),
        crosswordsmith_fill:mac_lm_of(LmA, TLen, TPos, LmT),
        crosswordsmith_fill:mac_support(1, LmS, SDom, LmT, 0, Supp),
        get_assoc(T, DomA0, TDom), TDom1 is TDom /\ Supp,
        Deleted is popcount(TDom)-popcount(TDom1),
        ( TDom1 =:= 0
        -> count_revision(dwo, Deleted), b0_bump_edge(EId), fail
        ; TDom1 =:= TDom
        -> count_revision(redundant, 0), DomA1=DomA0, Q1=Q0
        ;  count_revision(fruitful, Deleted),
           put_assoc(T, DomA0, TDom1, DomA1),
           ( memberchk(T,Q0) -> Q1=Q0 ; Q1=[T|Q0] ) )
    ;   DomA1=DomA0, Q1=Q0
    ),
    b0_revise(Es, SDom, SLen, LenA, LmA, Unfilled, DomA1, DomA, Q1, Q).

b0_bump_edge(EId) :-
    nb_getval(cw_mac_weights, WT), arg(EId, WT, V), V1 is V+1.0,
    nb_setarg(EId, WT, V1),
    ( probe_on
    -> nb_getval(cw_b0_bumps, B), arg(EId,B,N), N1 is N+1, nb_setarg(EId,B,N1),
       nb_getval(cw_b0_peaks, P), arg(EId,P,Old),
       ( V1 > Old -> nb_setarg(EId,P,V1) ; true ),
       nb_getval(cw_b0_attempt_peaks,AP), arg(EId,AP,AOld),
       ( V1 > AOld -> nb_setarg(EId,AP,V1) ; true )
    ;  true ).

probe_on :- nb_getval(cw_b0_mode, counters).
probe_timing :- nb_getval(cw_b0_mode, timing).

phase_start(T) :- ( probe_timing -> get_time(T) ; T=0.0 ).
phase_end(_, _) :- \+ probe_timing, !.
phase_end(I,T0) :-
    get_time(T1), D is T1-T0, nb_getval(cw_b0_times,T),
    arg(I,T,V), V1 is V+D, nb_setarg(I,T,V1).

count_node :- count_arg(6, 1).
count_place_dwo :- count_arg(5, 1).
count_revision(redundant, _) :- count_arg(1, 1).
count_revision(fruitful, Deleted) :- count_arg(2, 1), count_arg(4, Deleted).
count_revision(dwo, Deleted) :- count_arg(3, 1), count_arg(4, Deleted).

count_arg(_, _) :- \+ probe_on, !.
count_arg(I, D) :-
    nb_getval(cw_b0_counts, C), arg(I,C,V), V1 is V+D, nb_setarg(I,C,V1).

count_queue(_) :- \+ probe_on, !.
count_queue(Q0) :-
    length(Q0, N), nb_getval(cw_b0_queue,Q), arg(N,Q,V), V1 is V+1,
    nb_setarg(N,Q,V1), count_arg(7,1).

counter_snapshot(S) :-
    nb_getval(cw_b0_counts,C), compound_args(C,Counts),
    nb_getval(cw_b0_queue,Q), compound_args(Q,Queue),
    S = snap(Counts,Queue).

counter_delta(snap(C0,Q0), D) :-
    counter_snapshot(snap(C1,Q1)), maplist(minus,C1,C0,CD),
    maplist(minus,Q1,Q0,QD),
    CD=[Red,Fruit,Dwo,Deleted,PlaceDwo,Nodes,Pops],
    D=_{revision_redundant:Red,revision_fruitful:Fruit,
        revision_dwo:Dwo,deleted_popcount:Deleted,
        placement_dwo:PlaceDwo,nodes:Nodes,queue_pops:Pops,
        queue_lengths:QD}.

minus(A,B,D) :- D is A-B.
compound_args(T, L) :- functor(T,_,N), compound_args(1,N,T,L).
compound_args(I,N,_,[]) :- I>N, !.
compound_args(I,N,T,[V|Vs]) :- arg(I,T,V), I1 is I+1,
    compound_args(I1,N,T,Vs).

bump_vector(V) :- nb_getval(cw_b0_bumps,B), compound_args(B,V).
weight_vector(V) :- nb_getval(cw_mac_weights,W), compound_args(W,V).
peak_vector(V) :- nb_getval(cw_b0_peaks,P), compound_args(P,V).
attempt_peak_vector(V) :-
    nb_getval(cw_b0_attempt_peaks,P), compound_args(P,V).

list_compound(Name,Args,T) :- T=..[Name|Args].

record_attempt(Stop) :-
    nb_getval(cw_b0_attempt_active, active(A,Cap,Start,B0)),
    counter_delta(Start, Delta), bump_vector(B1), maplist(minus,B1,B0,BD),
    sum_list(BD, Learned), weight_vector(W), attempt_peak_vector(P),
    nb_getval(cw_mac_attempt_nodes, Nodes),
    Rec=Delta.put(_{attempt:A,cap:Cap,stop:Stop,attempt_nodes:Nodes,
                    learned_bumps:Learned,bump_counts:BD,
                    final_weights:W,peak_weights:P}),
    nb_getval(cw_b0_attempt_records,R0), append(R0,[Rec],R),
    nb_setval(cw_b0_attempt_records,R), nb_setval(cw_b0_attempt_active,none).

finish_active_attempt(Stop) :-
    ( nb_getval(cw_b0_attempt_active, active(_,_,_,_))
    -> record_attempt(Stop)
    ;  true ).

core_report(SearchSeconds, Report) :-
    (   nb_getval(cw_b0_setup_complete, true)
    ->  core_report_ready(SearchSeconds, Report)
    ;   Report=_{search_seconds:SearchSeconds,setup_complete:false,root:_{},
                 attempts:[],attempts_count:0,last_attempt_nodes:0,
                 bump_counts:[],final_weights:[],peak_weights:[],edges:[],
                 weight_concentration:_{learned_bumps:0,
                                        top_quartile_share:0.0,gini:0.0},
                 phase_seconds:_{selection_candidate:0.0,placement:0.0,
                                 propagation_support:0.0,other_search:0.0}}
    ).

core_report_ready(SearchSeconds, Report) :-
    nb_getval(cw_b0_root_record,Root),
    nb_getval(cw_b0_attempt_records,Attempts),
    nb_getval(cw_b0_attempt_count,AttemptCount),
    ( Attempts=[] -> LastNodes=0 ; last(Attempts,Last), LastNodes=Last.attempt_nodes ),
    bump_vector(Bumps), weight_vector(Weights), peak_vector(Peaks),
    concentration(Bumps, Concentration),
    nb_getval(cw_b0_edge_a,EdgeA),
    edge_rows(EdgeA, Bumps, Weights, Peaks, Edges),
    nb_getval(cw_b0_times,Times), compound_args(Times,PhaseValues),
    PhaseValues=[SelectionSeconds,PlacementSeconds,PropagationSeconds],
    OtherSeconds is max(0.0,SearchSeconds-SelectionSeconds-PlacementSeconds-
                             PropagationSeconds),
    Report=_{search_seconds:SearchSeconds,setup_complete:true,
             root:Root,attempts:Attempts,
             attempts_count:AttemptCount,last_attempt_nodes:LastNodes,
             bump_counts:Bumps,final_weights:Weights,peak_weights:Peaks,
             weight_concentration:Concentration,edges:Edges,
             phase_seconds:_{selection_candidate:SelectionSeconds,
                             placement:PlacementSeconds,
                             propagation_support:PropagationSeconds,
                             other_search:OtherSeconds}}.

edge_rows(EdgeA,Bumps,Weights,Peaks,Rows) :-
    assoc_to_list(EdgeA,Pairs),
    findall(EId-edge{edge:EId,a:A,b:B,bumps:N,final:W,peak:P},
            ( member(A-Es,Pairs), member(e(_,B,_,EId),Es), A<B,
              nth1(EId,Bumps,N),nth1(EId,Weights,W),nth1(EId,Peaks,P) ),
            Keyed), keysort(Keyed,Sorted), pairs_values(Sorted,Rows).

concentration([], _{learned_bumps:0,top_quartile_share:0.0,gini:0.0}).
concentration(V, C) :-
    sum_list(V,Sum), length(V,N),
    ( Sum=:=0 -> Share=0.0,Gini=0.0
    ; msort(V,Asc), reverse(Asc,Desc), K is max(1,(N+3)//4),
      length(Top,K), append(Top,_,Desc), sum_list(Top,TopSum), Share is TopSum/Sum,
      gini(Asc,Sum,N,Gini) ),
    C=_{learned_bumps:Sum,top_quartile_share:Share,gini:Gini}.

gini(Asc,Sum,N,G) :-
    weighted_sum(Asc,1,0,WS), G is (2.0*WS)/(N*Sum)-(N+1.0)/N.
weighted_sum([],_,S,S).
weighted_sum([X|Xs],I,S0,S) :- S1 is S0+I*X,I1 is I+1,
    weighted_sum(Xs,I1,S1,S).

score_layout([], _, _, _{n:0,mean:0.0,min:0,below50:0}).
score_layout(Numbered, Scores, DBL, Q) :-
    crosswordsmith_fill:fill_quality_stats(Numbered, Scores, DBL,
                                            stats(N,Mean,Min,Below)),
    Q=_{n:N,mean:Mean,min:Min,below50:Below}.

layout_digest([], none) :- !.
layout_digest(Numbered, Digest) :-
    with_output_to(string(S), write_canonical(Numbered)),
    sha_hash(S, Hash, [algorithm(sha256)]), hash_atom(Hash,Digest).

run_cli :-
    current_prolog_flag(argv, Argv),
    ( Argv=[compare,Grid,Dict,MinA,SeedsA,BudgetA]
    -> atom_number(MinA,Min),atom_number(BudgetA,Budget),seed_arg(SeedsA,Seeds),
       ( compare_product(Grid,Dict,Min,Seeds,Budget) -> halt(0);halt(1) )
    ; Argv=[run,Label,Grid,Dict,MinA,SeedsA,BudgetA,TimeoutA,Mode]
    -> atom_number(MinA,Min),atom_number(BudgetA,Budget),
       atom_number(TimeoutA,Timeout),seed_arg(SeedsA,Seeds),
       run_row(Label,Grid,Dict,Min,Seeds,Budget,Timeout,Mode,R),
       json_write_dict(current_output,R,[width(0)]),nl
    ; Argv=[profile,Label,Grid,Dict,MinA,SeedsA,BudgetA]
    -> atom_number(MinA,Min),atom_number(BudgetA,Budget),seed_arg(SeedsA,Seeds),
       run_profile(Label,Grid,Dict,Min,Seeds,Budget,R),
       json_write_dict(current_output,R,[width(0)]),nl
    ; format(user_error,
             'usage: ... -- compare GRID DICT MIN SEEDS BUDGET | run LABEL GRID DICT MIN SEEDS BUDGET TIMEOUT off|counters|timing | profile LABEL GRID DICT MIN SEEDS BUDGET~n',[]),
      halt(2) ).

seed_arg(none,none) :- !.
seed_arg(A,A).

truth_value(Goal, true) :- call(Goal), !.
truth_value(_, false).
