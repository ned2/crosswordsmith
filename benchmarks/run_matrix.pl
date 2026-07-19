#!/usr/bin/env swipl
% benchmarks/run_matrix.pl - strategy x fixture benchmark matrix.
%
% Runs every solver strategy against every fixture in the manifest
% (benchmarks/fixtures.pl), each on its OWN configured grid/start/iterations,
% and emits one CSV row per (strategy, fixture) cell. This is the harness used
% to produce comparable result batches for the experiment log
% (docs/experiments.md): regenerate the whole matrix in one command whenever
% the fixtures or an algorithm change.
%
% Usage:
%   swipl -q benchmarks/run_matrix.pl                 % all strategies
%   swipl -q benchmarks/run_matrix.pl -- baseline mrv_capped
%
% Wall time is machine-dependent; compare INFERENCES only across runs using the
% same SWI-Prolog version.

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists), [member/2]).
:- use_module(library(time), [call_with_time_limit/2]).
:- use_module('bench_paths.pl', [repo_path/2]).
:- repo_path('load.pl', Load), consult(Load).
:- use_module('bench_core.pl', [inproc_sampler/2, measure/3]).
:- use_module('bench_fixture.pl', [load_arrange_fixture/2]).
:- use_module('bench_report.pl', [swi_version/1]).
:- consult('fixtures.pl').

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    chosen_strategies(Argv, Strategies),
    emit_metadata(Strategies),
    format("strategy,fixture,grid,start,iterations,solved,wall_min_ms,wall_median_ms,inferences_min,inferences_median~n", []),
    forall(member(Strategy, Strategies),
           forall(bench_fixture(Rel, Grid, Start, Iters, Warmup),
                  run_cell(Strategy, Rel, Grid, Start, Iters, Warmup))).

% Strategies from argv (already atoms), else every strategy the solver core
% defines. Each named strategy is validated; an unknown one throws.
chosen_strategies([], Strategies) :-
    !,
    strategies(Strategies).
chosen_strategies(Argv, Argv) :-
    forall(member(S, Argv), require_strategy(S)).

emit_metadata(Strategies) :-
    swi_version(Ver),
    format("# tool: crosswordsmith-matrix~n", []),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# strategies: ~w~n", [Strategies]),
    format("# metric_note: inferences are same-SWI regression signals; wall time is machine-dependent~n", []).

% Per-cell wall-clock guard. A strategy that cannot solve a fixture within this
% many seconds is recorded as `timeout` rather than hanging the whole matrix
% (e.g. baseline on a hard mesh). Keep fixtures' manifest iterations at 1 for
% such fixtures so a slow-but-solved cell is not multiplied.
cell_limit_seconds(60).

run_cell(Strategy, Rel, Grid, Start, Iters, Warmup) :-
    repo_path(Rel, File),
    load_arrange_fixture(File, Words),
    file_base_name(Rel, Name),
    cell_limit_seconds(Limit),
    solve_status(Strategy, Words, Grid, Start, Limit, Status),
    ( Status == solved
    ->  % Warmup + measured loop + summary now come from bench_core (shared with
        % run_arrange). solve_status above keeps the timeout/`no` gate OUTSIDE
        % measure/3, which is success-only (plan §4). Only solve_status runs
        % under call_with_time_limit; the measured solves are unguarded, safe
        % because the search is deterministic - a cell that solved within the
        % probe's limit cannot newly hang here.
        measure(strategy_sample(Strategy, Words, Grid, Start),
                _{warmup: Warmup, iterations: Iters}, Summary),
        WMin is Summary.stats.wall.min * 1000.0,
        WMed is Summary.stats.wall.median * 1000.0,
        IMin = Summary.stats.inferences.min,
        IMed = Summary.stats.inferences.median,
        format("~w,~w,~d,~w,~d,yes,~3f,~3f,~0f,~0f~n",
               [Strategy, Name, Grid, Start, Iters, WMin, WMed, IMin, IMed])
    ;   % Status is `no` (search exhausted) or `timeout` (exceeded the limit);
        % distinguished so a slow cell is never mistaken for unsatisfiable.
        format("~w,~w,~d,~w,~d,~w,,,,~n", [Strategy, Name, Grid, Start, Iters, Status])
    ).

% Sampler closure for bench_core: one measured single-corner search -> one sample
% (wall/cpu/inferences). measure/3 calls this as call(strategy_sample(...), Sample).
strategy_sample(Strategy, Words, Grid, Start, Sample) :-
    inproc_sampler(solve_once(Strategy, Words, Grid, Start), Sample).

% solved | no | timeout, never throwing the time-limit error to the caller.
solve_status(Strategy, Words, Grid, Start, Limit, Status) :-
    catch( ( call_with_time_limit(Limit, solve_once(Strategy, Words, Grid, Start))
           ->  Status = solved
           ;   Status = no ),
           time_limit_exceeded,
           Status = timeout ).

solve_once(Strategy, Words, Grid, Start) :-
    find_crossword(Strategy, Grid, Words, Start, _Grid, _Placed),
    !.
