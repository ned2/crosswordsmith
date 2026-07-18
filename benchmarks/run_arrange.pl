#!/usr/bin/env swipl
% benchmarks/run_arrange.pl - the PRODUCT benchmark for `arrange`.
%
% For every workload in benchmarks/workloads.pl, measures two layers over the
% SAME word set and reports both plus their difference:
%
%   command  - end-to-end `crosswordsmith arrange ...` process latency (the number
%              a user feels): SWI startup + load + parse + strict two-corner
%              representative search + emit.
%   search   - the in-process budget-explicit arrange_best_layout/6 search alone.
%   rest     - command - search: the fixed CLI-wrapper overhead, clamped at 0.
%
% Wall/rss are machine-dependent; SEARCH INFERENCES are the portable, deterministic
% metric - compare those across runs and machines. rss is whole-process peak RSS
% (a footprint, not a search-memory metric).
%
% Usage:
%   swipl -q benchmarks/run_arrange.pl                       % all workloads, text
%   swipl -q benchmarks/run_arrange.pl -- --format csv
%   swipl -q benchmarks/run_arrange.pl -- --fixture bundled  % one workload
%   swipl -q benchmarks/run_arrange.pl -- --iterations 5 --warmup 1

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(optparse)).
:- use_module(library(json)).
:- use_module('bench_paths.pl', [repo_path/2]).
:- repo_path('load.pl', Load), consult(Load).
:- use_module('bench_core.pl').
:- use_module('bench_cli.pl').
:- use_module('bench_fixture.pl', [load_arrange_fixture/2]).
:- use_module('bench_report.pl', [benchmark_report/3, swi_version/1]).
:- use_module('subjects.pl').
:- consult('workloads.pl').

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    catch(parse_runner_options(arrange, Spec, Argv, _Opts, Common),
          E, (print_message(error, E), halt(2))),
    Common = runner_options{format:Fmt, fixture:Filter, heavy:Heavy,
                            iterations:ItOv, warmup:WuOv, help:Help},
    ( Help == true -> usage, halt(0) ; true ),
    selected_workloads(Filter, Heavy, Workloads),
    catch(require_selected(arrange, Filter, Workloads),
          E, (print_message(error, E), halt(2))),
    maplist(run_selected_workload(ItOv, WuOv), Workloads, Rows),
    emit(Fmt, Rows).

selected_workloads(Filter, Heavy, Workloads) :-
    findall(workload(F, Size, Mode, It0, Wu0, Exp, Tier, Gate, Budget, Words),
            ( arrange_workload(F, Size, Mode, It0, Wu0, Exp, Tier, Gate, Budget, Words),
              file_base_name(F, SelectionId),
              workload_selected(Filter, Heavy, SelectionId, Tier) ),
            Workloads).

run_selected_workload(ItOv, WuOv,
        workload(F, Size, Mode, It0, Wu0, Exp, Tier, Gate, Budget, Words), Row) :-
    apply_override(ItOv, It0, It),
    apply_override(WuOv, Wu0, Wu),
    format(user_error, "bench arrange: START ~w (~w)~n", [F, Tier]),
    run_workload(F, Size, Mode, It, Wu, Exp, Budget, Words, Row0),
    format(user_error, "bench arrange: DONE  ~w~n", [F]),
    Row = Row0.put(_{tier: Tier, gate: Gate}).

opts_spec(
    [ [opt(help),       type(boolean), default(false),
       shortflags([h]), longflags([help]),       help('show this help')],
      [opt(format),     type(atom),    default(text),
       longflags([format]),   help('output format: text | csv | json')],
      [opt(fixture),    type(atom),    default(''),
       longflags([fixture]),  help('only workloads whose fixture basename contains this substring (any tier)')],
      [opt(heavy),      type(boolean), default(false),
       longflags([heavy]),    help('also run the heavy tail: hard ladder rungs (subsecond to a few seconds) and the budget-saturating latency probe (~1 min)')],
      [opt(iterations), type(integer), default(-1),
       longflags([iterations]), help('override measured iterations for every workload')],
      [opt(warmup),     type(integer), default(-1),
       longflags([warmup]),   help('override warmup iterations for every workload')]
    ]).

usage :-
    opts_spec(Spec),
    opt_help(Spec, Help),
    format(user_output, "Usage: swipl -q benchmarks/run_arrange.pl -- [options]~n~n~w", [Help]).

% --- one workload: measure both layers, derive the breakdown -----------------
run_workload(Fixture, Size, Mode, Iters, Warmup, Expected, Budget, ExpectedWords, Row) :-
    repo_path(Fixture, File),
    load_arrange_fixture(File, Words),
    length(Words, NumWords),
    ( NumWords =:= ExpectedWords -> true
    ; throw(error(fixture_word_count(Fixture, expected(ExpectedWords), got(NumWords)), _)) ),
    file_base_name(Fixture, Name),
    crosswordsmith_exe(Exe),
    Opts = _{warmup: Warmup, iterations: Iters},
    measure(arrange_command_sampler(Exe, File, Size, Mode, Expected), Opts, Cmd),
    measure(arrange_search_sampler(Words, Size, Budget, Expected),    Opts, Search),
    CmdWallMinMs    is Cmd.stats.wall.min      * 1000.0,
    CmdWallMedMs    is Cmd.stats.wall.median   * 1000.0,
    CmdRssMedKiB    is round(Cmd.stats.rss.median),
    SearchWallMinMs is Search.stats.wall.min    * 1000.0,
    SearchWallMedMs is Search.stats.wall.median * 1000.0,
    SearchInfMin     = Search.stats.inferences.min,
    SearchInfMed     = Search.stats.inferences.median,
    RestMedMs       is max(0.0, CmdWallMedMs - SearchWallMedMs),
    ( CmdWallMedMs > 0.0 -> Share is 100.0 * SearchWallMedMs / CmdWallMedMs ; Share = 0.0 ),
    % warmup/budget/words ride along (with tier + gate, added by the caller) so a
    % recorder can build a COMPLETE baseline spec for a rung it has never seen
    % (check_baseline --record adds new ladder rungs from these fields).
    Row = _{ fixture:Name, size:Size, mode:Mode, iterations:Iters, expected:Expected,
             warmup:Warmup, budget:Budget, words:NumWords,
             cmd_wall_min_ms:CmdWallMinMs, cmd_wall_med_ms:CmdWallMedMs,
             cmd_rss_med_kib:CmdRssMedKiB,
             search_wall_min_ms:SearchWallMinMs, search_wall_med_ms:SearchWallMedMs,
             search_inf_min:SearchInfMin, search_inf_med:SearchInfMed,
             rest_med_ms:RestMedMs, search_share_pct:Share }.

crosswordsmith_exe(Exe) :-
    repo_path(crosswordsmith, Exe).

% --- output formats ----------------------------------------------------------
emit(text, Rows) :- !,
    metadata_lines,
    format("~w~t~30|~t~w~6+~t~w~13+~t~w~11+~t~w~15+~t~w~13+~t~w~12+~t~w~9+~n",
           ['fixture','size','cmd_wall','cmd_rss','srch_wall','srch_inf','rest','srch%']),
    forall(member(R, Rows), print_text_row(R)).
emit(csv, Rows) :- !,
    metadata_lines,
    format("fixture,size,mode,iterations,expected,cmd_wall_min_ms,cmd_wall_med_ms,cmd_rss_med_kib,search_wall_min_ms,search_wall_med_ms,search_inf_min,search_inf_med,rest_med_ms,search_share_pct~n"),
    forall(member(R, Rows), print_csv_row(R)).
emit(json, Rows) :- !,
    benchmark_report('crosswordsmith-arrange-bench', Rows, Doc0),
    Doc = Doc0.put(metric_note,
        'wall/rss machine-dependent; search inferences are the portable metric; rss is whole-process footprint'),
    json_write_dict(user_output, Doc, [width(80)]), nl.
emit(Other, _) :-
    format(user_error, "run_arrange: unknown --format ~w (use text|csv|json)~n", [Other]),
    halt(2).

metadata_lines :-
    swi_version(Ver),
    format("# tool: crosswordsmith-arrange-bench~n"),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# metric_note: wall/rss machine-dependent; search inferences are the portable metric; rss is whole-process footprint~n").

print_text_row(R) :-
    RssMB is R.cmd_rss_med_kib / 1024.0,
    Fixture = R.fixture, Size = R.size,
    CmdWall = R.cmd_wall_med_ms, SearchWall = R.search_wall_med_ms,
    SearchInf = R.search_inf_med, Rest = R.rest_med_ms, Share = R.search_share_pct,
    format("~w~t~30|~t~d~6+~t~1f ms~13+~t~1f MB~11+~t~1f ms~15+~t~d~13+~t~1f ms~12+~t~1f%~9+~n",
           [Fixture, Size, CmdWall, RssMB, SearchWall, SearchInf, Rest, Share]).

print_csv_row(R) :-
    format("~w,~d,~w,~d,~w,~3f,~3f,~d,~3f,~3f,~d,~d,~3f,~2f~n",
           [ R.fixture, R.size, R.mode, R.iterations, R.expected,
             R.cmd_wall_min_ms, R.cmd_wall_med_ms, R.cmd_rss_med_kib,
             R.search_wall_min_ms, R.search_wall_med_ms,
             R.search_inf_min, R.search_inf_med, R.rest_med_ms, R.search_share_pct ]).

:- multifile prolog:error_message//1.
prolog:error_message(fixture_word_count(Fixture, expected(Expected), got(Got))) -->
    [ 'benchmark fixture ~q has ~d words; manifest requires ~d'-[Fixture, Got, Expected] ].
