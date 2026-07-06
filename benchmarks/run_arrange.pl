#!/usr/bin/env swipl
% benchmarks/run_arrange.pl - the PRODUCT benchmark for `arrange`.
%
% For every workload in benchmarks/workloads.pl, measures two layers over the
% SAME word set and reports both plus their difference:
%
%   command  - end-to-end `crosswordsmith arrange ...` process latency (the number
%              a user feels): SWI startup + load + parse + 4-corner search + emit.
%   search   - the in-process arrange_best_layout/5 search alone.
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

% directory_file_path/3 is autoload-only (library(filesex)); explicit so this
% root also loads under autoload(false) (P11/C5, matching load.pl).
:- use_module(library(filesex), [directory_file_path/3]).

:- dynamic repo_root/1.

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', RepoRoot,
                      [ relative_to(BenchDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(BenchDir, 'bench_core.pl', BenchCore),
   use_module(BenchCore),
   directory_file_path(BenchDir, 'subjects.pl', Subjects),
   use_module(Subjects),
   directory_file_path(BenchDir, 'workloads.pl', Workloads),
   consult(Workloads).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    catch(opt_parse(Spec, Argv, Opts, _Pos), E, (print_message(error, E), halt(2))),
    ( memberchk(help(true), Opts) -> usage, halt(0) ; true ),
    memberchk(format(Fmt), Opts),
    memberchk(fixture(Filter), Opts),
    memberchk(heavy(Heavy), Opts),
    memberchk(iterations(ItOv), Opts),
    memberchk(warmup(WuOv), Opts),
    findall(Row,
            ( arrange_workload(F, Size, Mode, It0, Wu0, Exp, Tier, Gate, Budget),
              workload_selected(Filter, Heavy, F, Tier),
              apply_override(ItOv, It0, It),
              apply_override(WuOv, Wu0, Wu),
              run_workload(F, Size, Mode, It, Wu, Exp, Budget, Row0),
              Row = Row0.put(_{tier: Tier, gate: Gate}) ),
            Rows),
    emit(Fmt, Rows).

opts_spec(
    [ [opt(help),       type(boolean), default(false),
       shortflags([h]), longflags([help]),       help('show this help')],
      [opt(format),     type(atom),    default(text),
       longflags([format]),   help('output format: text | csv | json')],
      [opt(fixture),    type(atom),    default(''),
       longflags([fixture]),  help('only workloads whose fixture basename contains this substring (any tier)')],
      [opt(heavy),      type(boolean), default(false),
       longflags([heavy]),    help('also run the heavy tail: the hard ladder rungs (~0.6-9s each) and the budget-saturating latency probe (~1 min)')],
      [opt(iterations), type(integer), default(-1),
       longflags([iterations]), help('override measured iterations for every workload')],
      [opt(warmup),     type(integer), default(-1),
       longflags([warmup]),   help('override warmup iterations for every workload')]
    ]).

usage :-
    opts_spec(Spec),
    opt_help(Spec, Help),
    format(user_output, "Usage: swipl -q benchmarks/run_arrange.pl -- [options]~n~n~w", [Help]).

% An explicit --fixture filter selects across ANY tier (you asked for it by name);
% otherwise core always runs and heavy only under --heavy.
workload_selected(Filter, _Heavy, Fixture, _Tier) :-
    Filter \== '', !,
    file_base_name(Fixture, Base),
    sub_atom(Base, _, _, _, Filter).
workload_selected('', _Heavy, _Fixture, core) :- !.
workload_selected('', true,   _Fixture, heavy).

apply_override(-1, Default, Default) :- !.
apply_override(Override, _, Override).

% --- one workload: measure both layers, derive the breakdown -----------------
run_workload(Fixture, Size, Mode, Iters, Warmup, Expected, Budget, Row) :-
    repo_file(Fixture, File),
    read_clues(File, Words),
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
    length(Words, NumWords),
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
    repo_root(Root),
    directory_file_path(Root, crosswordsmith, Exe).

repo_file(Rel, File) :-
    repo_root(Root),
    directory_file_path(Root, Rel, File).

% Same strict clue reader as run_matrix: a fixture with no clues/1 is a hard error
% (an empty word set would "arrange" trivially and record a bogus row).
read_clues(File, Words) :-
    read_file_to_terms(File, Terms, []),
    (   memberchk(clues(Words0), Terms)
    ->  Words = Words0
    ;   throw(error(fixture_missing_clues(File), _))
    ).

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
    swi_version(Ver),
    Doc = _{ tool: 'crosswordsmith-arrange-bench', swi_prolog: Ver,
             metric_note: 'wall/rss machine-dependent; search inferences are the portable metric; rss is whole-process footprint',
             results: Rows },
    json_write_dict(user_output, Doc, [width(80)]), nl.
emit(Other, _) :-
    format(user_error, "run_arrange: unknown --format ~w (use text|csv|json)~n", [Other]),
    halt(2).

metadata_lines :-
    swi_version(Ver),
    format("# tool: crosswordsmith-arrange-bench~n"),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# metric_note: wall/rss machine-dependent; search inferences are the portable metric; rss is whole-process footprint~n").

swi_version(Ver) :-
    current_prolog_flag(version_data, V),
    ( V = swi(Ma, Mi, Pa, _) -> format(atom(Ver), '~d.~d.~d', [Ma, Mi, Pa]) ; Ver = V ).

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
prolog:error_message(fixture_missing_clues(Fixture)) -->
    [ 'benchmark fixture ~q does not define clues/1'-[Fixture] ].
