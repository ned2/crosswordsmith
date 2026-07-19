#!/usr/bin/env swipl
% benchmarks/run_fill.pl - the PRODUCT benchmark for `fill`.
%
% For every rung in benchmarks/fill_workloads.pl, measures FOUR attribution
% buckets over the same (grid, dict, seeds) and reports them plus a derived rest:
%
%   command   - end-to-end `crosswordsmith fill ...` process latency (SWI startup
%               + load.pl + dict load + slot derivation + search + emit).
%   dict_load - the in-process load_dict/3 alone. load_inf (inferences) is a
%               deterministic GATED metric.
%   grid      - the in-process fill_grid/4 slot derivation alone (its own bucket).
%   search    - the in-process budget-explicit fill_attempt/8 (FRESH slots per
%               sample, pre-loaded dict). search_inf is the GATED metric of record.
%   rest      - command - (dict_load + grid + search) walls, clamped at 0.
%
% Wall/rss are machine-dependent. SEARCH/LOAD INFERENCES are compared only
% against baselines recorded by the same SWI-Prolog version. rss is whole-process
% peak RSS (a footprint, not a search-memory metric).
%
% Usage:
%   swipl -q benchmarks/run_fill.pl                          % core rungs, text
%   swipl -q benchmarks/run_fill.pl -- --format csv
%   swipl -q benchmarks/run_fill.pl -- --heavy               % + the heavy tail
%   swipl -q benchmarks/run_fill.pl -- --fixture g11         % rungs whose id contains g11
%   swipl -q benchmarks/run_fill.pl -- --iterations 3 --warmup 1

:- set_prolog_flag(verbose, silent).
:- use_module(library(apply), [foldl/4, maplist/3]).
:- use_module(library(assoc), [assoc_to_values/2]).
:- use_module(library(json), [json_write_dict/3]).
:- use_module(library(lists), [member/2]).
:- use_module(library(optparse), [opt_help/2]).
:- use_module('bench_paths.pl', [repo_path/2]).
:- repo_path('load.pl', Load), consult(Load).
:- use_module('bench_core.pl', [measure/3]).
:- use_module('bench_cli.pl',
              [ apply_override/3,
                parse_runner_options/5,
                require_unique_ids/2,
                require_selected/3,
                workload_selected/4
              ]).
:- use_module('bench_report.pl', [benchmark_report/3, swi_version/1]).
:- use_module('fill_subjects.pl',
              [ fill_command_sampler/6,
                fill_load_sampler/2,
                fill_grid_sampler/2,
                fill_search_sampler/6
              ]).
:- consult('fill_workloads.pl').

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    catch(parse_runner_options(fill, Spec, Argv, _Opts, Common),
          E, (print_message(error, E), halt(2))),
    Common = runner_options{format:Fmt, fixture:Filter, heavy:Heavy,
                            iterations:ItOv, warmup:WuOv, help:Help},
    ( Help == true -> usage, halt(0) ; true ),
    selected_workloads(Filter, Heavy, Workloads),
    catch(require_selected(fill, Filter, Workloads),
          E, (print_message(error, E), halt(2))),
    findall(Id,
            member(workload(Id, _, _, _, _, _, _, _, _), Workloads),
            Ids),
    catch(require_unique_ids(fill, Ids),
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

% --- one rung: measure all four buckets, derive the breakdown ----------------
run_workload(Id, GridRel, DictRel, SeedsRel, Iters, Warmup, Expected, Budget, Row) :-
    repo_path(GridRel, GridFile),
    repo_path(DictRel, DictFile),
    ( SeedsRel == none -> SeedsFile = none ; repo_path(SeedsRel, SeedsFile) ),
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
    foldl(add_bucket_length, Buckets, 0, N).

add_bucket_length(Bucket, N0, N) :-
    length(Bucket, Length),
    N is N0 + Length.

crosswordsmith_exe(Exe) :-
    repo_path(crosswordsmith, Exe).

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
    benchmark_report('crosswordsmith-fill-bench', Rows, Doc0),
    Doc = Doc0.put(metric_note,
        'search_inf/load_inf are same-SWI gated metrics; wall/rss machine-dependent; rss is whole-process footprint'),
    json_write_dict(user_output, Doc, [width(80)]), nl.
emit(Other, _) :-
    format(user_error, "run_fill: unknown --format ~w (use text|csv|json)~n", [Other]),
    halt(2).

metadata_lines :-
    swi_version(Ver),
    format("# tool: crosswordsmith-fill-bench~n"),
    format("# swi_prolog: ~w~n", [Ver]),
    format("# metric_note: search_inf/load_inf are same-SWI gated metrics; wall/rss machine-dependent~n").

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
