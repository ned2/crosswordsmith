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
% Wall time is machine-dependent; treat INFERENCES as the primary, portable
% metric when comparing across runs or machines.

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(statistics)).
:- use_module(library(time)).        % call_with_time_limit/2 (per-cell guard)

:- dynamic repo_root/1.

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', RepoRoot,
                      [ relative_to(BenchDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(RepoRoot, 'benchmarks/fixtures.pl', Manifest),
   consult(Manifest).

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
    current_prolog_flag(version_data, V),
    ( V = swi(Ma, Mi, Pa, _) -> format(atom(Ver), '~d.~d.~d', [Ma, Mi, Pa]) ; Ver = V ),
    format("# tool: crosswordsmith-matrix~n", []),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# strategies: ~w~n", [Strategies]),
    format("# metric_note: wall time is machine-dependent; inferences are the portable metric~n", []).

% Per-cell wall-clock guard. A strategy that cannot solve a fixture within this
% many seconds is recorded as `timeout` rather than hanging the whole matrix
% (e.g. baseline on a hard mesh). Keep fixtures' manifest iterations at 1 for
% such fixtures so a slow-but-solved cell is not multiplied.
cell_limit_seconds(60).

run_cell(Strategy, Rel, Grid, Start, Iters, Warmup) :-
    repo_file(Rel, File),
    read_clues(File, Words),
    file_base_name(Rel, Name),
    cell_limit_seconds(Limit),
    solve_status(Strategy, Words, Grid, Start, Limit, Status),
    ( Status == solved
    ->  % Only solve_status (below) runs under call_with_time_limit; the warmup
        % and measured solves here are unguarded. Safe because the search is
        % deterministic - inference counts are identical across iterations, so a
        % cell that solved within the probe's limit cannot newly hang here.
        forall(between(1, Warmup, _), ignore(solve_once(Strategy, Words, Grid, Start))),
        findall(W-I,
                ( between(1, Iters, _),
                  call_time(solve_once(Strategy, Words, Grid, Start), T),
                  W is T.wall * 1000.0, I = T.inferences ),
                Pairs),
        pairs_keys_values(Pairs, Walls, Infs),
        col(Walls, WMin, WMed),
        col(Infs, IMin, IMed),
        format("~w,~w,~d,~w,~d,yes,~3f,~3f,~0f,~0f~n",
               [Strategy, Name, Grid, Start, Iters, WMin, WMed, IMin, IMed])
    ;   % Status is `no` (search exhausted) or `timeout` (exceeded the limit);
        % distinguished so a slow cell is never mistaken for unsatisfiable.
        format("~w,~w,~d,~w,~d,~w,,,,~n", [Strategy, Name, Grid, Start, Iters, Status])
    ).

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

% Min and "median" of a value list. For even-length lists this reports the upper
% of the two central values (no averaging), unlike run_benchmarks.pl's median/3.
% Immaterial for the authoritative inference metric (deterministic, so
% upper-median == true median); only wall-clock (machine noise) would differ.
col(Values, Min, Median) :-
    msort(Values, Sorted),
    Sorted = [Min|_],
    length(Sorted, Len),
    Mid is Len // 2,
    nth0(Mid, Sorted, Median).

repo_file(Rel, File) :-
    repo_root(Root),
    directory_file_path(Root, Rel, File).

read_clues(File, Words) :-
    setup_call_cleanup(open(File, read, S), read_loop(S, File, Words), close(S)).
% A fixture with no clues/1 term (typo, truncation, wrong path) must be a hard
% error: silently unifying Words=[] makes an empty puzzle "solve" trivially and
% records a bogus measured cell, corrupting the batch docs/experiments.md reads.
% Mirror run_benchmarks.pl's read_fixture_clues/3 (throw fixture_missing_clues/1).
read_loop(S, File, Words) :-
    read_term(S, T, []),
    ( T == end_of_file -> throw(error(fixture_missing_clues(File), _))
    ; T = clues(Words)  -> true
    ; read_loop(S, File, Words) ).

:- multifile prolog:error_message//1.
prolog:error_message(fixture_missing_clues(Fixture)) -->
    [ 'benchmark fixture ~q does not define clues/1'-[Fixture] ].
