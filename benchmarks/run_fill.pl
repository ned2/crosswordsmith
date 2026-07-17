#!/usr/bin/env swipl
% benchmarks/run_fill.pl - the PRODUCT benchmark for `fill`.
%
% For every rung in benchmarks/fill_workloads.pl, measures FOUR attribution
% buckets over the same (grid, dict, seeds) and reports them plus a derived rest:
%
%   command   - end-to-end `crosswordsmith fill ...` process latency (SWI startup
%               + load.pl + dict load + slot derivation + search + emit).
%   dict_load - the in-process load_dict/3 alone. load_inf (inferences) is the
%               deterministic, REPORTED metric (not gated - Phase 3 owns it).
%   grid      - the in-process fill_grid/4 slot derivation alone (its own bucket).
%   search    - the in-process budget-explicit fill_attempt/8 (FRESH slots per
%               sample, pre-loaded dict). search_inf is the GATED metric of record.
%   rest      - command - (dict_load + grid + search) walls, clamped at 0.
%
% Wall/rss are machine-dependent; SEARCH/LOAD INFERENCES are the portable,
% deterministic metrics - compare those across runs and machines. rss is
% whole-process peak RSS (a footprint, not a search-memory metric).
%
% Usage:
%   swipl -q benchmarks/run_fill.pl                          % core rungs, text
%   swipl -q benchmarks/run_fill.pl -- --format csv
%   swipl -q benchmarks/run_fill.pl -- --heavy               % + the heavy tail
%   swipl -q benchmarks/run_fill.pl -- --fixture g11         % rungs whose id contains g11
%   swipl -q benchmarks/run_fill.pl -- --iterations 3 --warmup 1

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(optparse)).
:- use_module(library(json)).
:- use_module(library(assoc)).
% yall: dict_word_count's metadata fold uses a lambda; explicit import so this
% root also runs under autoload(false) (P11/C5). NB importing yall also
% compile-expands the lambda — harness metadata only, outside every sampler
% (the bench-fill-check ratchet locks the measured counts bit-identical).
:- use_module(library(yall)).

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
   directory_file_path(BenchDir, 'bench_cli.pl', BenchCli),
   use_module(BenchCli),
   directory_file_path(BenchDir, 'fill_subjects.pl', Subjects),
   use_module(Subjects),
   directory_file_path(BenchDir, 'fill_workloads.pl', Workloads),
   consult(Workloads).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    catch(opt_parse(Spec, Argv, Opts, Pos), E, (print_message(error, E), halt(2))),
    memberchk(format(Fmt), Opts),
    memberchk(fixture(Filter), Opts),
    memberchk(heavy(Heavy), Opts),
    memberchk(iterations(ItOv), Opts),
    memberchk(warmup(WuOv), Opts),
    catch(validate_runner_options(fill, Pos, Fmt, ItOv, WuOv),
          E, (print_message(error, E), halt(2))),
    ( memberchk(help(true), Opts) -> usage, halt(0) ; true ),
    selected_workloads(Filter, Heavy, Workloads),
    catch(require_selected(fill, Filter, Workloads),
          E, (print_message(error, E), halt(2))),
    maplist(run_selected_workload(ItOv, WuOv), Workloads, Rows),
    emit(Fmt, Rows).

selected_workloads(Filter, Heavy, Workloads) :-
    findall(workload(Id, Grid, Dict, Seeds, It0, Wu0, Exp, Tier, Budget),
            ( fill_workload(Id, Grid, Dict, Seeds, It0, Wu0, Exp, Tier, Budget),
              workload_selected(Filter, Heavy, Id, Tier) ),
            Workloads).

run_selected_workload(ItOv, WuOv,
        workload(Id, Grid, Dict, Seeds, It0, Wu0, Exp, Tier, Budget), Row) :-
    apply_override(ItOv, It0, It),
    apply_override(WuOv, Wu0, Wu),
    run_workload(Id, Grid, Dict, Seeds, It, Wu, Exp, Budget, Row0),
    Row = Row0.put(_{tier: Tier}).

opts_spec(
    [ [opt(help),       type(boolean), default(false),
       shortflags([h]), longflags([help]),       help('show this help')],
      [opt(format),     type(atom),    default(text),
       longflags([format]),   help('output format: text | csv | json')],
      [opt(fixture),    type(atom),    default(''),
       longflags([fixture]),  help('only rungs whose id contains this substring (any tier)')],
      [opt(heavy),      type(boolean), default(false),
       longflags([heavy]),    help('also run the heavy (10-35M-inference) tail rungs')],
      [opt(iterations), type(integer), default(-1),
       longflags([iterations]), help('override measured iterations for every rung')],
      [opt(warmup),     type(integer), default(-1),
       longflags([warmup]),   help('override warmup iterations for every rung')]
    ]).

usage :-
    opts_spec(Spec),
    opt_help(Spec, Help),
    format(user_output, "Usage: swipl -q benchmarks/run_fill.pl -- [options]~n~n~w", [Help]).

% An explicit --fixture filter selects across ANY tier; otherwise core always
% runs and heavy only under --heavy.
workload_selected(Filter, _Heavy, Id, _Tier) :-
    Filter \== '', !,
    sub_atom(Id, _, _, _, Filter).
workload_selected('', _Heavy, _Id, core) :- !.
workload_selected('', true,   _Id, heavy).

apply_override(-1, Default, Default) :- !.
apply_override(Override, _, Override).

% --- one rung: measure all four buckets, derive the breakdown ----------------
run_workload(Id, GridRel, DictRel, SeedsRel, Iters, Warmup, Expected, Budget, Row) :-
    repo_file(GridRel, GridFile),
    repo_file(DictRel, DictFile),
    ( SeedsRel == none -> SeedsFile = none ; repo_file(SeedsRel, SeedsFile) ),
    crosswordsmith_exe(Exe),
    % pre-load the dict ONCE so the search layer isolates search cost; also grab
    % the grid size + dict word count for the row metadata.
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    dict_word_count(DictByLen, NWords),
    crosswordsmith_fill:fill_grid(GridFile, Size, _Slots, _CV),
    Opts = _{warmup: Warmup, iterations: Iters},
    measure(fill_command_sampler(Exe, GridFile, DictFile, SeedsFile, Expected), Opts, Cmd),
    measure(fill_load_sampler(DictFile),                                        Opts, Load),
    measure(fill_grid_sampler(GridFile),                                        Opts, Grid),
    measure(fill_search_sampler(GridFile, SeedsFile, DictByLen, Index, Budget-Expected), Opts, Search),
    CmdWallMinMs   is Cmd.stats.wall.min      * 1000.0,
    CmdWallMedMs   is Cmd.stats.wall.median   * 1000.0,
    CmdRssMedKiB   is round(Cmd.stats.rss.median),
    LoadInfMed      = Load.stats.inferences.median,
    LoadWallMedMs  is Load.stats.wall.median  * 1000.0,
    GridInfMed      = Grid.stats.inferences.median,
    GridWallMedMs  is Grid.stats.wall.median  * 1000.0,
    SearchInfMin    = Search.stats.inferences.min,
    SearchInfMed    = Search.stats.inferences.median,
    SearchWallMinMs is Search.stats.wall.min    * 1000.0,
    SearchWallMedMs is Search.stats.wall.median * 1000.0,
    RestMedMs      is max(0.0, CmdWallMedMs - (LoadWallMedMs + GridWallMedMs + SearchWallMedMs)),
    share(LoadWallMedMs,   CmdWallMedMs, LoadShare),
    share(GridWallMedMs,   CmdWallMedMs, GridShare),
    share(SearchWallMedMs, CmdWallMedMs, SearchShare),
    share(RestMedMs,       CmdWallMedMs, RestShare),
    seeds_label(SeedsRel, SeedsLabel),
    % metadata rides along so the recorder can build a COMPLETE first-seen spec.
    Row = _{ rung:Id, grid_file:GridRel, dict_file:DictRel, seeds:SeedsLabel,
             size:Size, words:NWords, iterations:Iters, warmup:Warmup,
             budget:Budget, expected:Expected,
             cmd_wall_min_ms:CmdWallMinMs, cmd_wall_med_ms:CmdWallMedMs,
             cmd_rss_med_kib:CmdRssMedKiB,
             load_inf:LoadInfMed, load_wall_med_ms:LoadWallMedMs,
             grid_inf:GridInfMed, grid_wall_med_ms:GridWallMedMs,
             search_inf_min:SearchInfMin, search_inf_med:SearchInfMed,
             search_wall_min_ms:SearchWallMinMs, search_wall_med_ms:SearchWallMedMs,
             rest_med_ms:RestMedMs,
             load_share_pct:LoadShare, grid_share_pct:GridShare,
             search_share_pct:SearchShare, rest_share_pct:RestShare }.

share(_, Cmd, 0.0) :- Cmd =< 0.0, !.
share(Part, Cmd, Pct) :- Pct is 100.0 * Part / Cmd.

seeds_label(none, none) :- !.
seeds_label(Path, Base) :- file_base_name(Path, Base).

dict_word_count(DictByLen, N) :-
    assoc_to_values(DictByLen, Buckets),
    foldl([B, A0, A1]>>(length(B, L), A1 is A0 + L), Buckets, 0, N).

crosswordsmith_exe(Exe) :-
    repo_root(Root),
    directory_file_path(Root, crosswordsmith, Exe).

repo_file(Rel, File) :-
    repo_root(Root),
    directory_file_path(Root, Rel, File).

% --- output formats ----------------------------------------------------------
emit(text, Rows) :- !,
    metadata_lines,
    format("~w~t~16|~t~w~6+~t~w~9+~t~w~13+~t~w~13+~t~w~13+~t~w~11+~t~w~9+~n",
           ['rung','size','words','load_inf','grid_inf','search_inf','cmd_ms','srch%']),
    forall(member(R, Rows), print_text_row(R)).
emit(csv, Rows) :- !,
    metadata_lines,
    format("rung,size,words,seeds,dict,load_inf,grid_inf,search_inf_min,search_inf_med,cmd_wall_med_ms,cmd_rss_med_kib,load_share_pct,grid_share_pct,search_share_pct,rest_share_pct~n"),
    forall(member(R, Rows), print_csv_row(R)).
emit(json, Rows) :- !,
    swi_version(Ver),
    Doc = _{ tool: 'crosswordsmith-fill-bench', swi_prolog: Ver,
             metric_note: 'search_inf gated; load_inf reported (Phase 3); wall/rss machine-dependent; rss is whole-process footprint',
             results: Rows },
    json_write_dict(user_output, Doc, [width(80)]), nl.
emit(Other, _) :-
    format(user_error, "run_fill: unknown --format ~w (use text|csv|json)~n", [Other]),
    halt(2).

metadata_lines :-
    swi_version(Ver),
    format("# tool: crosswordsmith-fill-bench~n"),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# metric_note: search_inf gated; load_inf reported; wall/rss machine-dependent~n").

swi_version(Ver) :-
    current_prolog_flag(version_data, V),
    ( V = swi(Ma, Mi, Pa, _) -> format(atom(Ver), '~d.~d.~d', [Ma, Mi, Pa]) ; Ver = V ).

print_text_row(R) :-
    format("~w~t~16|~t~d~6+~t~d~9+~t~D~13+~t~D~13+~t~D~13+~t~1f~11+~t~1f%~9+~n",
           [R.rung, R.size, R.words, R.load_inf, R.grid_inf, R.search_inf_med,
            R.cmd_wall_med_ms, R.search_share_pct]).

print_csv_row(R) :-
    format("~w,~d,~d,~w,~w,~d,~d,~d,~d,~3f,~d,~2f,~2f,~2f,~2f~n",
           [ R.rung, R.size, R.words, R.seeds, R.dict_file,
             R.load_inf, R.grid_inf, R.search_inf_min, R.search_inf_med,
             R.cmd_wall_med_ms, R.cmd_rss_med_kib,
             R.load_share_pct, R.grid_share_pct, R.search_share_pct, R.rest_share_pct ]).
