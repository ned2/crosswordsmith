#!/usr/bin/env swipl
% benchmarks/start_sensitivity.pl - start-position sensitivity sweep.
%
% The regular matrix (run_matrix.pl) pins each fixture to its manifest start
% position. This sweeps EVERY named start position (start_locs/1) for a small
% set of strategies, to show how solvability and cost vary with where the first
% word is seeded. One measured solve per cell (compare inferences within one SWI
% version; wall is not reported) with a 30 s per-cell timeout, so a hard cell
% records `no` or `timeout` instead of hanging.
%
% Usage:  swipl -q benchmarks/start_sensitivity.pl
%
% See docs/experiments.md ("Suite property: start-position bias").

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists), [member/2]).
:- use_module(library(time), [call_with_time_limit/2]).
% call_time/2 is autoload-only (library(statistics)); explicit so this root
% also runs under autoload(false) (P11/C5).
:- use_module(library(statistics), [call_time/2]).
:- use_module('bench_paths.pl', [repo_path/2]).
:- repo_path('load.pl', Load), consult(Load).
:- use_module('bench_fixture.pl', [load_arrange_fixture/2]).
:- consult('fixtures.pl').

:- initialization(main, main).

cell_limit_seconds(30).
sweep_strategies([baseline, mrv_inc]).

main :-
    sweep_strategies(Strategies),
    start_locs(Starts),
    format("strategy,fixture,grid,start,result,inferences~n", []),
    forall(member(Strategy, Strategies),
      forall(bench_fixture(Rel, Grid, _ManifestStart, _, _),
        forall(member(Start, Starts),
          run_cell(Strategy, Rel, Grid, Start)))).

run_cell(Strategy, Rel, Grid, Start) :-
    repo_path(Rel, File),
    load_arrange_fixture(File, Words),
    file_base_name(Rel, Name),
    cell_limit_seconds(Limit),
    ( catch( call_with_time_limit(Limit,
               ( call_time(once_solve(Strategy, Grid, Words, Start), T)
                 -> Res = yes(T.inferences) ; Res = no )),
             time_limit_exceeded, Res = timeout )
    -> true ; Res = no ),
    ( Res = yes(Inf) -> format("~w,~w,~d,~w,yes,~d~n", [Strategy,Name,Grid,Start,Inf])
    ; Res = timeout  -> format("~w,~w,~d,~w,timeout,~n", [Strategy,Name,Grid,Start])
    ;                   format("~w,~w,~d,~w,no,~n", [Strategy,Name,Grid,Start]) ).

once_solve(Strategy, Grid, Words, Start) :-
    find_crossword(Strategy, Grid, Words, Start, _, _),
    !.
