% Focused arrange baseline promotion safety tests. The recorder is exercised
% only against temporary baseline/history files; committed benchmark state is
% never opened for writing.

:- use_module(library(plunit)).
:- use_module(library(json)).
:- use_module(library(readutil)).
:- use_module(library(time), [call_with_time_limit/2]).
:- use_module(library(process), [process_create/3, process_wait/2]).
:- use_module(library(apply), [exclude/3]).
:- use_module('../benchmarks/check_baseline.pl', []).
:- use_module('../benchmarks/check_fill_baseline.pl', []).
:- use_module('../benchmarks/bench_process.pl', [capture_process/6]).
:- use_module('../benchmarks/bench_cli.pl', [checker_mode/3]).
:- use_module('bench_core_caller.pl', []).

:- begin_tests(arrange_benchmark_promotion).

test(core_only_promotion_retains_heavy) :-
    baseline_doc(['core.pl'-spec(core, inf, 100),
                  'heavy.pl'-spec(heavy, inf, 500)], Baseline),
    result_row('core.pl', core, inf, 90, CoreRow),
    run_temp_promotion(Baseline, [CoreRow], Written, HistoryLines, Fails),
    assertion(Fails =:= 0),
    baseline_inf(Written, 'core.pl', 90),
    baseline_inf(Written, 'heavy.pl', 500),
    assertion(HistoryLines =:= 1).

test(new_rung_promotion_persists_complete_row) :-
    baseline_doc(['core.pl'-spec(core, inf, 100)], Baseline),
    result_row('new.pl', heavy, inf, 321, NewRow),
    run_temp_promotion(Baseline, [NewRow], Written, _HistoryLines, Fails),
    assertion(Fails =:= 0),
    baseline_spec(Written, 'new.pl', Spec),
    assertion(Spec.search_inf =:= 321),
    assertion(Spec.cmd_wall_med_ms =:= 12.5),
    assertion(Spec.cmd_rss_med_kib =:= 12345),
    assertion(Spec.grid =:= 17),
    assertion(Spec.words =:= 16),
    assertion(Spec.tier == "heavy"),
    assertion(Spec.gate == "inf"),
    assertion(Spec.iterations =:= 1),
    assertion(Spec.warmup =:= 0),
    assertion(Spec.budget =:= 500000000).

test(latency_only_promotion_persists_without_corruption) :-
    baseline_doc(['core.pl'-spec(core, inf, 100),
                  'latency.pl'-spec(heavy, latency, 500)], Baseline),
    result_row('latency.pl', heavy, latency, 900, LatencyRow),
    run_temp_promotion(Baseline, [LatencyRow], Written, _HistoryLines, Fails),
    assertion(Fails =:= 0),
    baseline_inf(Written, 'latency.pl', 900),
    baseline_inf(Written, 'core.pl', 100).

test(regression_rejected_without_writes) :-
    baseline_doc(['core.pl'-spec(core, inf, 100)], Baseline),
    result_row('core.pl', core, inf, 110, RegressedRow),
    with_temp_files(Baseline, BaselinePath, HistoryPath,
        ( check_arrange_baseline:promote_doc(
              BaselinePath, HistoryPath, Baseline,
              _{swi_prolog:"10.1.10", results:[RegressedRow]}, [], Fails, _Wins),
          check_arrange_baseline:load_baseline(BaselinePath, Written),
          baseline_inf(Written, 'core.pl', 100),
          history_line_count(HistoryPath, HistoryLines),
          assertion(Fails =:= 1),
          assertion(HistoryLines =:= 0) )).

test(meta_closures_resolve_in_caller_module) :-
    bench_core_test_caller:measure_local(Summary),
    assertion(Summary.stats.value.median =:= 7),
    bench_core_test_caller:inproc_local(Sample),
    assertion(number(Sample.wall)),
    assertion(number(Sample.cpu)),
    assertion(integer(Sample.inferences)).

test(dual_capture_drains_large_stderr_before_stdout_close) :-
    Goal = 'forall(between(1,1048576,_),put_char(user_error,x)),format("done"),close(user_output)',
    call_with_time_limit(
        10,
        capture_process(path(swipl), ['-q', '-g', Goal, '-t', halt], capture,
                        Stdout, Stderr, Status)),
    assertion(Status == exit(0)),
    assertion(Stdout == "done"),
    string_length(Stderr, StderrBytes),
    assertion(StderrBytes =:= 1048576).

test(process_cleanup_waits_after_goal_exception) :-
    tmp_file_stream(text, Path, Out),
    setup_call_cleanup(
        true,
        catch(setup_call_cleanup(
                  bench_process:start_child(
                      path(swipl), ['-q', '-g', true, '-t', halt], Out, null,
                      Child),
                  throw(capture_test_exception),
                  bench_process:reap_child(Child)),
              capture_test_exception, true),
        ( close(Out), delete_file(Path) )),
    Child = child(PID, pending),
    catch(process_wait(PID, _), ReapedError, true),
    assertion(nonvar(ReapedError)).

test(fill_core_record_retains_heavy_and_verifies_readback) :-
    fill_baseline_doc(['core'-fill_spec(core, 100, 200, 300),
                       'heavy'-fill_spec(heavy, 500, 600, 700)], Baseline),
    fill_result_row('core', core, 90, 190, 290, CoreRow),
    run_temp_fill_record(Baseline, [CoreRow], Written),
    fill_baseline_spec(Written, 'core', CoreSpec),
    assertion(CoreSpec.search_inf =:= 90),
    assertion(CoreSpec.load_inf =:= 190),
    assertion(CoreSpec.grid_inf =:= 290),
    fill_baseline_spec(Written, 'heavy', HeavySpec),
    assertion(HeavySpec.search_inf =:= 500),
    assertion(HeavySpec.load_inf =:= 600),
    assertion(HeavySpec.grid_inf =:= 700).

test(fill_new_rung_record_is_complete) :-
    fill_baseline_doc([], Baseline),
    fill_result_row('new', heavy, 321, 654, 987, NewRow),
    run_temp_fill_record(Baseline, [NewRow], Written),
    fill_baseline_spec(Written, 'new', Spec),
    assertion(Spec.search_inf =:= 321),
    assertion(Spec.load_inf =:= 654),
    assertion(Spec.grid_inf =:= 987),
    assertion(Spec.cmd_wall_med_ms =:= 12.5),
    assertion(Spec.cmd_rss_med_kib =:= 12345),
    assertion(Spec.grid_file == "grids/test.json"),
    assertion(Spec.dict_file == "dicts/test.txt"),
    assertion(Spec.seeds == "none"),
    assertion(Spec.grid =:= 15),
    assertion(Spec.words =:= 10),
    assertion(Spec.tier == "heavy"),
    assertion(Spec.iterations =:= 1),
    assertion(Spec.warmup =:= 0),
    assertion(Spec.budget =:= 800000000).

test(fill_failed_write_leaves_live_baseline_unchanged) :-
    fill_baseline_doc(['core'-fill_spec(core, 100, 200, 300)], Baseline),
    fill_result_row('core', core, 90, 190, 290, Row0),
    Row = Row0.put(cmd_wall_med_ms, not_json(foo)),
    with_temp_fill_baseline(Baseline, Path,
        ( catch(check_fill_baseline:do_record(
                    Path, Baseline, _{swi_prolog:"10.1.10", results:[Row]}),
                _, true),
          check_fill_baseline:load_baseline(Path, Written),
          fill_baseline_spec(Written, 'core', Spec),
          get_dict(search_inf, Spec, SearchInf),
          assertion(SearchInf =:= 100),
          atom_concat(Path, '.tmp', TempPath),
          assertion(\+ exists_file(TempPath)) )).

test(runners_reject_invalid_arguments_before_measurement,
     [forall((runner_script(Script), invalid_runner_args(Args)))]) :-
    append(['-q', Script, '--'], Args, Argv),
    capture_process(path(swipl), Argv, capture, _Stdout, _Stderr, Status),
    assertion(Status == exit(2)).

test(checker_mode_rejects_conflicts,
     [throws(error(bench_conflicting_modes([record,log]), _))]) :-
    checker_mode(['--record', '--log'], _, _).

test(checkers_reject_conflicting_modes_before_measurement,
     [forall(checker_script(Script))]) :-
    capture_process(path(swipl),
                    ['-q', Script, '--record', '--log'], capture,
                    _Stdout, _Stderr, Status),
    assertion(Status == exit(2)).

test(checkers_reject_history_runner_arguments,
     [forall(checker_script(Script))]) :-
    capture_process(path(swipl),
                    ['-q', Script, '--history', '--heavy'], capture,
                    _Stdout, _Stderr, Status),
    assertion(Status == exit(2)).

test(identity_scripts_reject_unknown_options,
     [forall(identity_script(Script))]) :-
    capture_process(path(bash), [Script, '--definitely-invalid'], capture,
                    _Stdout, _Stderr, Status),
    assertion(Status == exit(2)).

test(fill_identity_partial_record_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_fill_identity_fixture(ManifestPath, WorkloadsPath, CliPath, Original),
        ( format(atom(ManifestEnv), 'FILL_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'FILL_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(CliEnv), 'FILL_IDENTITY_CLI=~w', [CliPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, CliEnv,
                           'benchmarks/check_fill_identity.sh', '--record'],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(1)),
          read_file_to_string(ManifestPath, Written, []),
          assertion(Written == Original) ),
        maplist(delete_if_exists, [ManifestPath, WorkloadsPath, CliPath])).

:- end_tests(arrange_benchmark_promotion).

baseline_doc(Pairs0, _{host:"test-host", swi_prolog:"10.1.10",
                       regression_tolerance_pct:0.5, workloads:Workloads}) :-
    maplist(baseline_pair, Pairs0, Pairs),
    dict_pairs(Workloads, workloads, Pairs).

baseline_pair(Name-spec(Tier, Gate, Inf), Name-Spec) :-
    Spec = _{search_inf:Inf, cmd_wall_med_ms:10.0, cmd_rss_med_kib:12000,
             tier:Tier, gate:Gate, grid:15, words:10, iterations:1,
             warmup:0, budget:2000000000,
             info_only:["cmd_wall_med_ms", "cmd_rss_med_kib"]}.

result_row(Fixture, Tier, Gate, Inf,
           _{fixture:Fixture, search_inf_med:Inf, cmd_wall_med_ms:12.5,
             cmd_rss_med_kib:12345, tier:Tier, gate:Gate, size:17, words:16,
             iterations:1, warmup:0, budget:500000000}).

run_temp_promotion(Baseline, Rows, Written, HistoryLines, Fails) :-
    with_temp_files(Baseline, BaselinePath, HistoryPath,
        ( check_arrange_baseline:promote_doc(
              BaselinePath, HistoryPath, Baseline,
              _{swi_prolog:"10.1.10", results:Rows}, [], Fails, _Wins),
          check_arrange_baseline:load_baseline(BaselinePath, Written),
          history_line_count(HistoryPath, HistoryLines) )).

with_temp_files(Baseline, BaselinePath, HistoryPath, Goal) :-
    tmp_file_stream(text, BaselinePath, BaselineStream),
    json_write_dict(BaselineStream, Baseline),
    close(BaselineStream),
    tmp_file_stream(text, HistoryPath, HistoryStream),
    close(HistoryStream),
    setup_call_cleanup(
        true,
        Goal,
        ( delete_file(BaselinePath), delete_file(HistoryPath) )).

baseline_spec(Baseline, Fixture, Spec) :-
    check_arrange_baseline:find_baseline(Baseline.workloads, Fixture, Spec).

baseline_inf(Baseline, Fixture, Inf) :-
    baseline_spec(Baseline, Fixture, Spec),
    assertion(Spec.search_inf =:= Inf).

history_line_count(Path, Count) :-
    read_file_to_string(Path, Text, []),
    split_string(Text, "\n", "\n", Lines0),
    exclude(==(""), Lines0, Lines),
    length(Lines, Count).

fill_baseline_doc(Pairs0, _{host:"test-host", swi_prolog:"10.1.10",
                             regression_tolerance_pct:0.5,
                             workloads:Workloads}) :-
    maplist(fill_baseline_pair, Pairs0, Pairs),
    dict_pairs(Workloads, workloads, Pairs).

fill_baseline_pair(Name-fill_spec(Tier, Search, Load, Grid), Name-Spec) :-
    Spec = _{search_inf:Search, load_inf:Load, grid_inf:Grid,
             cmd_wall_med_ms:10.0, cmd_rss_med_kib:12000,
             grid_file:"grids/original.json", dict_file:"dicts/original.txt",
             seeds:"none", grid:15, words:10, tier:Tier,
             iterations:1, warmup:0, budget:800000000,
             info_only:["grid_inf", "cmd_wall_med_ms", "cmd_rss_med_kib"]}.

fill_result_row(Rung, Tier, Search, Load, Grid,
                _{rung:Rung, search_inf_med:Search, load_inf:Load,
                  grid_inf:Grid, cmd_wall_med_ms:12.5,
                  cmd_rss_med_kib:12345, grid_file:"grids/test.json",
                  dict_file:"dicts/test.txt", seeds:"none", size:15,
                  words:10, tier:Tier, iterations:1, warmup:0,
                  budget:800000000}).

run_temp_fill_record(Baseline, Rows, Written) :-
    with_temp_fill_baseline(Baseline, Path,
        ( check_fill_baseline:do_record(
              Path, Baseline, _{swi_prolog:"10.1.10", results:Rows}),
          check_fill_baseline:load_baseline(Path, Written) )).

with_temp_fill_baseline(Baseline, Path, Goal) :-
    tmp_file_stream(text, Path, Stream),
    json_write_dict(Stream, Baseline),
    close(Stream),
    setup_call_cleanup(true, Goal, delete_file(Path)).

fill_baseline_spec(Baseline, Rung, Spec) :-
    check_fill_baseline:find_baseline(Baseline.workloads, Rung, Spec).

runner_script('benchmarks/run_arrange.pl').
runner_script('benchmarks/run_fill.pl').
runner_script('benchmarks/run_arrange_greedy.pl').

invalid_runner_args([positional]).
invalid_runner_args(['--format', bogus]).
invalid_runner_args(['--iterations', '0']).
invalid_runner_args(['--warmup', '-2']).
invalid_runner_args(['--fixture', 'definitely-no-such-workload']).

checker_script('benchmarks/check_baseline.pl').
checker_script('benchmarks/check_fill_baseline.pl').
checker_script('benchmarks/check_greedy_baseline.pl').

identity_script('benchmarks/check_arrange_identity.sh').
identity_script('benchmarks/check_fill_identity.sh').
identity_script('benchmarks/check_fill_identity_artifact.sh').
identity_script('benchmarks/check_greedy_identity.sh').

make_fill_identity_fixture(ManifestPath, WorkloadsPath, CliPath, Original) :-
    Original = "original\tmanifest\n",
    tmp_file_stream(text, ManifestPath, ManifestStream),
    format(ManifestStream, '~s', [Original]),
    close(ManifestStream),
    tmp_file_stream(text, WorkloadsPath, WorkloadsStream),
    format(WorkloadsStream,
           "fill_workload(ok,'ok-grid',dict,none,1,0,filled,core,1).~n", []),
    format(WorkloadsStream,
           "fill_workload(fail,'fail-grid',dict,none,1,0,filled,core,1).~n", []),
    close(WorkloadsStream),
    tmp_file_stream(text, CliPath, CliStream),
    format(CliStream, "#!/usr/bin/env bash~n", []),
    format(CliStream, "for arg in \"$@\"; do~n", []),
    format(CliStream, "    [ \"$arg\" = fail-grid ] && exit 1~n", []),
    format(CliStream, "done~n", []),
    format(CliStream, "printf 'ok\\n'~n", []),
    close(CliStream),
    process_create(path(chmod), ['u+x', CliPath], [process(PID)]),
    process_wait(PID, exit(0)).

delete_if_exists(Path) :-
    ( exists_file(Path) -> delete_file(Path) ; true ).
