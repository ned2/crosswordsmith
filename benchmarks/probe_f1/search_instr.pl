% benchmarks/probe_f1/search_instr.pl - P-F1 probe: counter-instrumented COPY of
% fill.pl's search (fill_search/select_mrv/candidates/candidate_count/
% slot_bucket/index_intersection), PROBE-ONLY, beside - never inside - the
% measured engine. fill.pl is untouched.
%
% Instrumentation model:
%   - reg(Key, Goal): accumulate statistics(inferences) delta of a DETERMINISTIC
%     span into an nb_ global. Region deltas approximate the CLEAN cost of the
%     span (bookkeeping sits outside the delta except the exit statistics call +
%     the call/1 wrapper; calib/0 measures that per-entry constant d_in so it
%     can be subtracted).
%   - counters (cnt/1) + per-node trace (assertz) sit OUTSIDE all regions, so
%     they inflate only the instrumented TOTAL, never a region. Attribution is
%     therefore reported as region shares of the CLEAN total (regions ~= clean
%     spans), never of the instrumented total.
%   - leaf mode (pf1_leaf): the per-slot count path (select_mrv's full recount)
%     runs VERBATIM copies when off (V1: r_count_all is a clean superset
%     measure) and per-call leaf regions when on (V2: splits r_bound_c /
%     r_inter_c out of the recount; nested bookkeeping calibrated out).
%   - the try loop (member + \+memberchk + unify) is NOT region-wrapped (it is
%     nondeterministic); its clean cost is derived as the residual
%     clean_total - sum(V1 regions), with counters (tries/placements) sizing it.
%
% EQUIVALENCE: same predicates, same clause order, same standard-order sorts,
% counters never fail and add no choicepoints; the one structural addition is
% the ( Body ; note_fail, fail ) disjunction, which fires exactly once per
% FAILING activation (a succeeded activation's choicepoint is pruned by the
% caller's once/1 before any redo can reach it). Verified externally: the
% instrumented run's filled grid must be term-identical to the clean run's.

:- module(probe_f1_search,
          [ pf1_reset/1,            % pf1_reset(+LeafMode:boolean)
            pf1_search/5,           % pf1_search(+Slots,+DictByLen,+Index,+Used,+Depth0)
            pf1_counter/2,          % pf1_counter(?Key, ?Value)
            pf1_node/2,             % trace: pf1_node(Depth, BestCount) in order
            pf1_fail/1,             % trace: pf1_fail(Depth) in order
            pf1_calib/0
          ]).

:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(assoc)).
:- use_module(library(ordsets)).

:- dynamic pf1_node/2.
:- dynamic pf1_fail/1.

counter_key(pf1_nodes).      % select_mrv activations
counter_key(pf1_tries).      % member/2 yields (words tried)
counter_key(pf1_place).      % placements (passed memberchk; unification attempted)
counter_key(pf1_inter_c).    % count-path index_intersection calls (leaf mode only)
counter_key(r_count_all).    % region: select_mrv full recount (maplist slot_candidate_count)
counter_key(r_sortsel).      % region: select_mrv sort + once(select)
counter_key(r_bound_m).      % region: winner-path slot_bucket Bound findall
counter_key(r_inter_m).      % region: winner-path index_intersection
counter_key(r_mat).          % region: winner-path maplist(nth0_of) materialization
counter_key(r_bound_c).      % region: count-path Bound findall (leaf mode)
counter_key(r_inter_c).      % region: count-path index_intersection (leaf mode)
counter_key(calib_r).        % calibration scratch
counter_key(calib_c).        % calibration scratch

pf1_reset(Leaf) :-
    forall(counter_key(K), nb_setval(K, 0)),
    nb_setval(pf1_leaf, Leaf),
    retractall(pf1_node(_, _)),
    retractall(pf1_fail(_)).

pf1_counter(K, V) :- counter_key(K), nb_getval(K, V).

cnt(Key) :- nb_getval(Key, N0), N is N0 + 1, nb_setval(Key, N).

reg(Key, Goal) :-
    statistics(inferences, I0),
    call(Goal),
    statistics(inferences, I1),
    nb_getval(Key, A0), A is A0 + (I1 - I0), nb_setval(Key, A).

% --- instrumented copies (verbatim from fill.pl but for the marked lines) ----

pf1_search([], _DictByLen, _Index, _Used, _Depth) :- !.
pf1_search(Slots, DictByLen, Index, Used, Depth) :-
    p_select_mrv(Slots, DictByLen, Index, slot(_, _, _, Vars), Rest, Cands, Depth),
    Depth1 is Depth + 1,                                    % [instr]
    (   member(Word, Cands),
        cnt(pf1_tries),                                     % [instr]
        \+ memberchk(Word, Used),
        cnt(pf1_place),                                     % [instr]
        Vars = Word,
        pf1_search(Rest, DictByLen, Index, [Word|Used], Depth1)
    ;   assertz(pf1_fail(Depth)), fail                      % [instr]
    ).

p_select_mrv(Slots, DictByLen, Index, Best, Rest, BestCands, Depth) :-
    cnt(pf1_nodes),                                         % [instr]
    nb_getval(pf1_leaf, Leaf),                              % [instr]
    (   Leaf == true                                        % [instr]
    ->  reg(r_count_all, maplist(p_scc_leaf(DictByLen, Index), Slots, Counted))
    ;   reg(r_count_all, maplist(p_scc_plain(DictByLen, Index), Slots, Counted))
    ),
    reg(r_sortsel,
        ( sort(0, @=<, Counted, [c(BestCount, BestStart, BestDir)|_]),
          Best = slot(BestStart, BestDir, _, _),
          once(select(Best, Slots, Rest)) )),
    Best = slot(_, _, _, Vars),
    p_candidates(Vars, DictByLen, Index, BestCands),
    assertz(pf1_node(Depth, BestCount)).                    % [instr]

% count path, V1: VERBATIM slot_candidate_count/candidate_count/slot_bucket -
% zero instrumentation inside, so r_count_all measures the clean recount cost.
p_scc_plain(DictByLen, Index, slot(Start, Dir, _, Vars), c(Count, Start, Dir)) :-
    p_candidate_count_plain(Vars, DictByLen, Index, Count).

p_candidate_count_plain(Vars, DictByLen, Index, Count) :-
    p_slot_bucket_plain(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> length(Words, Count)
    ; Sel = idx(Indices), length(Indices, Count)
    ).

p_slot_bucket_plain(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    findall(P-V, ( nth0(P, Vars, V), nonvar(V) ), Bound),
    ( Bound == []
    ->  Sel = all
    ;   p_index_intersection(Bound, Len, Index, Indices), Sel = idx(Indices)
    ).

% count path, V2 (leaf mode): same code with per-call leaf regions.
p_scc_leaf(DictByLen, Index, slot(Start, Dir, _, Vars), c(Count, Start, Dir)) :-
    p_candidate_count_leaf(Vars, DictByLen, Index, Count).

p_candidate_count_leaf(Vars, DictByLen, Index, Count) :-
    p_slot_bucket_leaf(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> length(Words, Count)
    ; Sel = idx(Indices), length(Indices, Count)
    ).

p_slot_bucket_leaf(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    reg(r_bound_c, findall(P-V, ( nth0(P, Vars, V), nonvar(V) ), Bound)),
    ( Bound == []
    ->  Sel = all
    ;   cnt(pf1_inter_c),                                   % [instr]
        reg(r_inter_c, p_index_intersection(Bound, Len, Index, Indices)),
        Sel = idx(Indices)
    ).

% winner path: candidates/4 with leaf regions (entries = nodes, so overhead is
% per-node, negligible).
p_candidates(Vars, DictByLen, Index, Cands) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    reg(r_bound_m, findall(P-V, ( nth0(P, Vars, V), nonvar(V) ), Bound)),
    ( Bound == []
    ->  Cands = Words
    ;   reg(r_inter_m, p_index_intersection(Bound, Len, Index, Indices)),
        reg(r_mat, maplist(p_nth0_of(Words), Indices, Cands))
    ).

p_nth0_of(Words, I, W) :- nth0(I, Words, W).

p_index_intersection([P-Ch|Rest], Len, Index, Indices) :-
    p_index_set(Len, P, Ch, Index, S0),
    foldl(p_index_intersect(Len, Index), Rest, S0, Indices).
p_index_intersect(Len, Index, P-Ch, Acc, Acc1) :-
    p_index_set(Len, P, Ch, Index, S), ord_intersection(Acc, S, Acc1).
p_index_set(Len, P, Ch, Index, S) :-
    ( get_assoc(k(Len, P, Ch), Index, S0) -> S = S0 ; S = [] ).

% --- calibration: per-entry bookkeeping constants ----------------------------
% d_in    : inferences a reg/2 entry adds INSIDE its own delta (call/1 wrapper +
%           exit statistics call) - subtract entries*d_in from leaf regions.
% t_total : inferences one reg(k, true) entry adds inside an ENCLOSING region -
%           the inflation r_count_all suffers per nested leaf entry (leaf mode).
% c_cnt   : inferences one cnt/1 adds inside an enclosing region.
pf1_calib :-
    N = 100000,
    nb_setval(calib_r, 0),
    statistics(inferences, A0),
    forall(between(1, N, _), reg(calib_r, true)),
    statistics(inferences, A1),
    nb_getval(calib_r, DIn_Total),
    % forall/between loop baseline
    statistics(inferences, B0),
    forall(between(1, N, _), true),
    statistics(inferences, B1),
    Loop is B1 - B0,
    DIn is DIn_Total / N,
    TTotal is (A1 - A0 - Loop) / N,
    nb_setval(calib_c, 0),
    statistics(inferences, C0),
    forall(between(1, N, _), cnt(calib_c)),
    statistics(inferences, C1),
    CCnt is (C1 - C0 - Loop) / N,
    LoopPer is Loop / N,
    format("calib: d_in=~3f t_total=~3f c_cnt=~3f loop_per_iter=~3f~n",
           [DIn, TTotal, CCnt, LoopPer]).
