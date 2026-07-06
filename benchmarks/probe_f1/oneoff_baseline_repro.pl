#!/usr/bin/env swipl
% P-F1 one-off: reproduce a rung's search_inf through the REAL bench pathway
% (bench_core:measure + fill_subjects:fill_search_sampler) inside a probe
% process, to locate the +1 offset seen from the probe runner's direct
% call_time. PROBE-ONLY.
%   swipl -q benchmarks/probe_f1/oneoff_baseline_repro.pl -- sq04_50k [stack]
:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).

:- dynamic repo_root/1.
:- prolog_load_context(directory, ProbeDir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(ProbeDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(RepoRoot, 'benchmarks/bench_core.pl', BC), use_module(BC),
   directory_file_path(RepoRoot, 'benchmarks/fill_subjects.pl', FS), use_module(FS),
   directory_file_path(RepoRoot, 'benchmarks/fill_workloads.pl', WL), consult(WL).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, [RungA|Rest]),
    ( Rest == [stack] -> set_prolog_flag(stack_limit, 4_000_000_000) ; true ),
    fill_workload(RungA, GridRel, DictRel, none, _, _, Expected, _, Budget),
    repo_root(Root),
    directory_file_path(Root, GridRel, GridFile),
    directory_file_path(Root, DictRel, DictFile),
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    measure(fill_search_sampler(GridFile, none, DictByLen, Index, Budget-Expected),
            _{warmup:1, iterations:1}, S),
    format("bench_pathway search_inf=~w (stack flag: ~w)~n",
           [S.stats.inferences.median, Rest]),
    % same process, direct call_time (the probe runner's pathway)
    fill_subjects:build_search_slots(GridFile, none, slots(SS, AS)),
    call_time(crosswordsmith_fill:fill_attempt(SS, AS, DictByLen, Index, Budget,
                                               O2, _, _), T2),
    format("direct call_time search_inf=~d outcome=~w~n", [T2.inferences, O2]),
    % and a bare non-qualified-context variant via a local wrapper
    fill_subjects:build_search_slots(GridFile, none, slots(SS3, AS3)),
    G3 = crosswordsmith_fill:fill_attempt(SS3, AS3, DictByLen, Index, Budget, O3, _, _),
    call_time(G3, T3),
    format("via-var call_time search_inf=~d outcome=~w~n", [T3.inferences, O3]).
