#!/usr/bin/env swipl
% benchmarks/probe_f1/load_stages.pl - P-F1 probe: decompose load_dict/3 stage
% by stage (campaign Phase 1, question A). PROBE-ONLY: never merged to master;
% lives beside, not inside, the measured path (fill.pl is untouched).
%
% Method: replicate the load_dict/3 + build_index/2 pipeline stage-by-stage,
% running each stage on the PREVIOUS stage's output under call_time/2, calling
% the module-internal originals (crosswordsmith_fill:read_file_lines etc.)
% wherever a stage is a named predicate, and verbatim copies of the inline
% goals otherwise. Reconcile: sum(stage inf) vs call_time(load_dict/3) whole.
%
% Usage: swipl -q benchmarks/probe_f1/load_stages.pl -- fixtures/dict/enable_10k.txt

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(assoc)).
:- use_module(library(pairs)).
:- use_module(library(ordsets)).

:- dynamic repo_root/1.
:- prolog_load_context(directory, ProbeDir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(ProbeDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, [DictRel|_]),
    repo_root(Root),
    directory_file_path(Root, DictRel, File),
    % Warm the load path once (first-call JIT; Phase 0 measured the cold/warm
    % delta at 2,383 inf on 10k, 0 on full ENABLE), then take the WARM whole.
    call_time(crosswordsmith_fill:load_dict(File, _, _), _TCold),
    call_time(crosswordsmith_fill:load_dict(File, _, _), TWhole),
    % --- stage-by-stage, each on the previous stage's output ---------------
    call_time(crosswordsmith_fill:read_file_lines(File, Lines), T1),
    call_time(findall(W, ( member(L, Lines),
                           crosswordsmith_fill:normalize_word(L, W),
                           W \== [] ), Ws0), T2),
    call_time(sort(Ws0, Words), T3),
    call_time(map_list_to_pairs(length, Words, LPairs), T4),
    call_time(keysort(LPairs, LSorted), T5),
    call_time(group_pairs_by_key(LSorted, LGroups), T6),
    call_time(list_to_assoc(LGroups, DictByLen), T7),
    % build_index/2 stages (verbatim copies of fill.pl:133-141 goals)
    call_time(findall(k(Len, P, Ch)-Idx,
                      ( gen_assoc(Len, DictByLen, WsL),
                        nth0(Idx, WsL, W), nth0(P, W, Ch) ),
                      Triples), T8),
    call_time(keysort(Triples, SortedT), T9),
    call_time(group_pairs_by_key(SortedT, Grouped), T10),
    call_time(maplist([K-Is, K-Set]>>list_to_ord_set(Is, Set), Grouped, GroupedSets), T11),
    call_time(list_to_assoc(GroupedSets, _Index), T12),
    % --- report -------------------------------------------------------------
    length(Lines, NLines), length(Words, NWords), length(Triples, NTriples),
    length(Grouped, NKeys),
    format("dict=~w lines=~d words=~d triples=~d index_keys=~d~n",
           [DictRel, NLines, NWords, NTriples, NKeys]),
    Stages = [ 'read_file_lines(read+split)'-T1,
               'normalize findall'-T2,
               'sort/2 dedupe'-T3,
               'map_list_to_pairs(length)'-T4,
               'keysort LPairs'-T5,
               'group_pairs_by_key len'-T6,
               'list_to_assoc DictByLen'-T7,
               'build_index findall triples'-T8,
               'build_index keysort'-T9,
               'build_index group_pairs'-T10,
               'build_index ord_set per key'-T11,
               'build_index list_to_assoc'-T12 ],
    WholeInf = TWhole.inferences,
    WholeWall is TWhole.wall * 1000.0,
    format("~w~t~30| ~t~w~12+ ~t~w~10+ ~t~w~8+~n", ['stage','inf','wall_ms','inf%']),
    forall(member(Name-T, Stages),
           ( I = T.inferences, Wms is T.wall*1000.0, Pct is 100.0*I/WholeInf,
             format("~w~t~30| ~t~d~12+ ~t~2f~10+ ~t~2f~8+~n", [Name, I, Wms, Pct]) )),
    findall(I, ( member(_-T, Stages), get_dict(inferences, T, I) ), Infs),
    sum_list(Infs, SumInf),
    findall(W, ( member(_-T, Stages), get_dict(wall, T, W0), W is W0*1000.0 ), Walls),
    sum_list(Walls, SumWall),
    ResInf is WholeInf - SumInf,
    ResPct is 100.0*ResInf/WholeInf,
    format("~w~t~30| ~t~d~12+ ~t~2f~10+~n", ['SUM(stages)', SumInf, SumWall]),
    format("~w~t~30| ~t~d~12+ ~t~2f~10+~n", ['WHOLE load_dict/3 (warm)', WholeInf, WholeWall]),
    format("residue (whole-sum): ~d inf = ~3f% of whole~n", [ResInf, ResPct]).
