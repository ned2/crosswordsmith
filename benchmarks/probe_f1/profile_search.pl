#!/usr/bin/env swipl
% benchmarks/probe_f1/profile_search.pl - P-F1 probe: SWI profile/2 over the
% CLEAN fill_attempt/8 for one rung; per-predicate call/redo counts are exact
% (the reconciliation signal); time is sampled (treat with suspicion).
% PROBE-ONLY.  Usage: swipl -q benchmarks/probe_f1/profile_search.pl -- g17_50k

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(statistics)).

:- dynamic repo_root/1.
:- prolog_load_context(directory, ProbeDir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(ProbeDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(RepoRoot, 'benchmarks/fill_workloads.pl', WL),
   consult(WL).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, [RungA|_]),
    fill_workload(RungA, GridRel, DictRel, none, _, _, _, _, Budget),
    repo_root(Root),
    directory_file_path(Root, GridRel, GridFile),
    directory_file_path(Root, DictRel, DictFile),
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    % warm once
    crosswordsmith_fill:fill_grid(GridFile, _, W0, _),
    crosswordsmith_fill:fill_attempt(W0, W0, DictByLen, Index, Budget, _, _, _),
    % profiled run on fresh slots
    crosswordsmith_fill:fill_grid(GridFile, _, S1, _),
    profile(crosswordsmith_fill:fill_attempt(S1, S1, DictByLen, Index, Budget,
                                             _O, _, _),
            [top(0)]),
    profile_data(Data),
    Summary = Data.summary,
    format("summary: ~q~n", [Summary]),
    format("predicate,calls,redos,exits,ticks_self,ticks_siblings~n"),
    findall(Node, member(Node, Data.nodes), Nodes),
    forall(member(N, Nodes),
           ( P = N.predicate, C = N.call, R = N.redo, E = N.exit,
             T = N.ticks_self, TS = N.ticks_siblings,
             ( C + R > 0
             ->  format("~q,~d,~d,~d,~d,~d~n", [P, C, R, E, T, TS])
             ;   true ) )).
