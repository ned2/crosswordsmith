#!/usr/bin/env swipl
% benchmarks/probe_f1/run_search.pl - P-F1 probe runner for the search layer.
% One rung per invocation (command hygiene). PROBE-ONLY.
%
%   swipl -q benchmarks/probe_f1/run_search.pl -- RUNG MODE
%     RUNG: a fill_workloads.pl rung id (e.g. g17_50k)
%     MODE: v1 (coarse regions) | v2 (leaf regions) | calib (constants only)
%
% Sequence per invocation:
%   1. dict loaded once; CLEAN reference: fresh slots, fill_attempt/8 at 2e9
%      under call_time -> must be `filled`; search_inf printed for the +0.00%
%      check against fill_baseline.json; solution term captured.
%   2. INSTRUMENTED run: fresh slots, pf1_search under call_time; solution term
%      compared == against the clean one (choice-identity gate).
%   3. counters/regions/trace dump (deciles, fail depths, derived work counts).

:- set_prolog_flag(verbose, silent).
%% [pf1 A/B] :- set_prolog_flag(stack_limit, 4_000_000_000).
:- use_module(library(lists)).
:- use_module(library(apply)).

:- dynamic repo_root/1.
:- prolog_load_context(directory, ProbeDir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(ProbeDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(RepoRoot, 'benchmarks/fill_workloads.pl', WL),
   consult(WL),
   directory_file_path(ProbeDir, 'search_instr.pl', SI),
   use_module(SI).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, [RungA, ModeA|_]),
    ( ModeA == calib -> pf1_calib, halt(0) ; true ),
    fill_workload(RungA, GridRel, DictRel, Seeds, _It, _Wu, Expected, _Tier, Budget),
    Seeds == none,   % the four probe rungs are unseeded; keep the harness simple
    repo_root(Root),
    directory_file_path(Root, GridRel, GridFile),
    directory_file_path(Root, DictRel, DictFile),
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    % --- clean reference (warm: one throwaway attempt first, Phase 0's 319-inf
    % cold delta) -------------------------------------------------------------
    build_slots(GridFile, WSlots0, _),
    crosswordsmith_fill:fill_attempt(WSlots0, WSlots0, DictByLen, Index, Budget,
                                     _WOutcome, _, _),
    build_slots(GridFile, CSlots, CAll),
    call_time(crosswordsmith_fill:fill_attempt(CSlots, CAll, DictByLen, Index,
                                               Budget, Outcome, _N, _IW), TClean),
    ( Outcome == Expected -> true ; throw(error(pf1_outcome(Outcome), _)) ),
    solution_sig(CAll, CleanSig),
    format("rung=~w clean_search_inf=~d clean_wall_ms=~2f~n",
           [RungA, TClean.inferences, TClean.wall*1000.0]),
    % --- instrumented run ----------------------------------------------------
    ( ModeA == v2 -> Leaf = true ; Leaf = false ),
    build_slots(GridFile, ISlots, IAll),
    pf1_reset(Leaf),
    length(ISlots, NSlots),
    call_time(once(pf1_search(ISlots, DictByLen, Index, [], 0)), TInstr),
    solution_sig(IAll, InstrSig),
    ( InstrSig == CleanSig
    ->  format("equivalence=IDENTICAL (~d slots)~n", [NSlots])
    ;   format("equivalence=MISMATCH~n"),
        format("clean=~q~ninstr=~q~n", [CleanSig, InstrSig]),
        halt(3)
    ),
    format("instr_total_inf=~d instr_wall_ms=~2f mode=~w~n",
           [TInstr.inferences, TInstr.wall*1000.0, ModeA]),
    % --- counters + regions --------------------------------------------------
    forall(( member(K, [pf1_nodes, pf1_tries, pf1_place, pf1_inter_c,
                        r_count_all, r_sortsel, r_bound_m, r_inter_m, r_mat,
                        r_bound_c, r_inter_c]),
             pf1_counter(K, V) ),
           format("~w=~d~n", [K, V])),
    % --- trace-derived stats -------------------------------------------------
    findall(D-C, pf1_node(D, C), Nodes),
    length(Nodes, NN),
    findall(D, member(D-_, Nodes), Depths),
    max_list(Depths, MaxD),
    % count_calls = sum over nodes of |slots remaining| = NSlots - depth
    findall(R, ( member(D-_, Nodes), R is NSlots - D ), Rs), sum_list(Rs, CountCalls),
    format("nodes_traced=~d n_search_slots=~d max_depth_reached=~d count_calls=~d~n",
           [NN, NSlots, MaxD, CountCalls]),
    findall(FD, pf1_fail(FD), FDs),
    length(FDs, NF),
    ( NF > 0
    ->  max_list(FDs, MaxFD), sum_list(FDs, SumFD), MeanFD is SumFD/NF,
        min_list(FDs, MinFD),
        format("node_failures=~d fail_depth_min=~d fail_depth_mean=~2f fail_depth_max=~d~n",
               [NF, MinFD, MeanFD, MaxFD])
    ;   format("node_failures=0~n") ),
    pf1_counter(pf1_place, Placed),
    Unwound is Placed - NSlots,
    format("placements=~d unwound_placements=~d~n", [Placed, Unwound]),
    % --- BestCount distribution per depth decile ------------------------------
    format("decile,nodes,min_cand,median_cand,max_cand~n"),
    forall(between(0, 9, Dec),
           ( findall(C, ( member(D-C, Nodes), Dec =:= min(9, D*10//NSlots) ), Cs),
             ( Cs == []
             ->  true
             ;   msort(Cs, S), length(S, L), min_list(S, Mn), max_list(S, Mx),
                 MI is L//2, nth0(MI, S, Med),
                 format("~d,~d,~d,~d,~d~n", [Dec, L, Mn, Med, Mx]) ) )).

build_slots(GridFile, Slots, Slots) :-
    crosswordsmith_fill:fill_grid(GridFile, _Size, Slots, _CV).

solution_sig(Slots, Sig) :-
    findall(Start-Dir-Word,
            ( member(slot(Start, Dir, _, Vars), Slots), atom_chars(Word, Vars) ),
            Sig).
