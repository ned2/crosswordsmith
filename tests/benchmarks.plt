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
:- use_module('../benchmarks/check_greedy_baseline.pl', []).
:- use_module('../benchmarks/bench_process.pl', [capture_process/6]).
:- use_module('../benchmarks/bench_cli.pl',
              [ apply_override/3,
                 checker_mode/3,
                 exact_runner_args/2,
                 parse_runner_options/5,
                 require_persistence_args/1,
                 require_selected/3,
                 workload_selected/4
              ]).
:- use_module('../benchmarks/bench_fixture.pl', [load_arrange_fixture/2]).
:- use_module('../benchmarks/bench_paths.pl', [benchmark_path/2, repo_path/2]).
:- use_module('../benchmarks/bench_report.pl',
              [benchmark_report/3, swi_version/1]).
:- use_module('../benchmarks/subjects.pl', []).
:- use_module('../benchmarks/fill_subjects.pl', []).
:- use_module('../benchmarks/bench_exact.pl', [exact_metric/4, exact_presence/5]).
:- use_module('../benchmarks/bench_store.pl',
              [append_history/4, read_history/2, read_json_dict/2,
               render_history/3]).
:- use_module('bench_core_caller.pl', []).
:- use_module('bench_store_caller.pl', []).

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
    deterministic(MeasureDet),
    assertion(MeasureDet == true),
    assertion(Summary.stats.value.median =:= 7),
    bench_core_test_caller:inproc_local(Sample),
    deterministic(InprocDet),
    assertion(InprocDet == true),
    assertion(number(Sample.wall)),
    assertion(number(Sample.cpu)),
    assertion(integer(Sample.inferences)).

test(sample_schema_rejects_missing_numeric_metric,
     [throws(error(bench_sample_schema_mismatch([cpu,wall], [wall]), _))]) :-
    bench_core:summarize_samples([_{cpu:1, wall:2}, _{wall:3}], _).

test(sample_schema_rejects_added_numeric_metric,
     [throws(error(bench_sample_schema_mismatch([wall], [rss,wall]), _))]) :-
    bench_core:summarize_samples([_{tag:ok, wall:1},
                                  _{rss:2, tag:ok, wall:3}], _).

test(subject_modules_resolve_standalone_dependencies,
     [forall(standalone_subject_goal(Goal))]) :-
    capture_process(path(swipl), ['-q', '-g', Goal, '-t', halt], capture,
                    _Stdout, _Stderr, Status),
    assertion(Status == exit(0)).

test(subject_sampler_exports_are_deterministic) :-
    bench_subjects:arrange_command_sampler(
        '/bin/true', ignored, 3, size, placed, _ArrangeCommand),
    deterministic(ArrangeCommandDet),
    assertion(ArrangeCommandDet == true),
    load_clues('fixtures/benchmark_08_words.pl', ArrangeWords),
    bench_subjects:arrange_search_sampler(
        ArrangeWords, 13, 2_000_000_000, placed, _ArrangeSearch),
    deterministic(ArrangeSearchDet),
    assertion(ArrangeSearchDet == true),
    fill_subjects:fill_command_sampler(
        '/bin/true', ignored_grid, ignored_dict, none, filled, _FillCommand),
    deterministic(FillCommandDet),
    assertion(FillCommandDet == true),
    fill_subjects:fill_load_sampler(
        'fixtures/wordlist_sample.txt', _FillLoad),
    deterministic(FillLoadDet),
    assertion(FillLoadDet == true),
    fill_subjects:fill_grid_sampler('fixtures/fill_grid_3.json', _FillGrid),
    deterministic(FillGridDet),
    assertion(FillGridDet == true),
    crosswordsmith_fill:load_dict(
        'fixtures/wordlist_sample.txt', DictByLen, Index),
    fill_subjects:fill_search_sampler(
        'fixtures/fill_grid_3.json', none, DictByLen, Index,
        800_000_000-filled, _FillSearch),
    deterministic(FillSearchDet),
    assertion(FillSearchDet == true).

test(arrange_command_sampler_is_semidet_for_invalid_mode, [fail]) :-
    bench_subjects:arrange_command_sampler(
        '/bin/true', ignored, 3, unsupported, placed, _).

test(exact_mode_rejects_one_inference_increase) :-
    exact_arrange_fails(101, "10.1.10", Fails),
    assertion(Fails =:= 1),
    exact_metric(100, 101, increase, 1).

test(exact_mode_rejects_one_inference_decrease) :-
    exact_arrange_fails(99, "10.1.10", Fails),
    assertion(Fails =:= 1),
    exact_metric(100, 99, decrease, 1).

test(exact_mode_accepts_identical_inferences) :-
    exact_arrange_fails(100, "10.1.10", Fails),
    assertion(Fails =:= 0).

test(exact_presence_rejects_duplicate_measured_ids,
     [throws(error(bench_exact_duplicate_ids(["same","same"]), _))]) :-
    exact_presence([same], [same, "same"], _Missing, _Unexpected, _Fails).

test(exact_mode_rejects_measured_gate_drift,
     [throws(error(
         bench_record_protocol_mismatch('core.pl', gate, inf, latency), _))]) :-
    baseline_doc(['core.pl'-spec(core, inf, 100)], Baseline),
    result_row('core.pl', core, latency, 101, Row),
    with_output_to(string(_),
                   check_arrange_baseline:do_exact_check(
                       Baseline, _{swi_prolog:"10.1.10", results:[Row]}, _)).

test(exact_mode_skips_reference_latency_metric) :-
    baseline_doc(['latency.pl'-spec(heavy, latency, 100)], Baseline),
    result_row('latency.pl', heavy, latency, 101, Row),
    with_output_to(string(_),
                   check_arrange_baseline:do_exact_check(
                       Baseline, _{swi_prolog:"10.1.10", results:[Row]}, Fails)),
    assertion(Fails =:= 0).

test(exact_mode_requires_same_swi) :-
    exact_arrange_fails(100, "99.0.0", Fails),
    assertion(Fails =:= 1).

test(cross_version_decreases_are_informational,
     [forall(cross_version_classifier(Goal, Kind))]) :-
    call(Goal),
    assertion(Kind == version_warn).

test(promotions_reject_different_swi,
     [forall(promotion_version_goal(Goal))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = error(bench_swi_version_mismatch("10.1.10", "10.2.0"), _)).

test(cross_version_records_require_complete_rows,
     [forall(record_version_goal(incomplete, Goal))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = error(bench_version_migration_incomplete(["heavy"], []), _)).

test(cross_version_records_accept_complete_rows,
      [forall(record_version_goal(complete, Goal))]) :-
    call(Goal),
    deterministic(Det),
    assertion(Det == true).

test(records_reject_measurement_protocol_changes,
     [forall(record_protocol_mismatch_goal(Goal))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = error(
        bench_record_protocol_mismatch(_, warmup, 0, 1), _)).

test(comparisons_reject_measurement_protocol_changes,
     [forall(comparison_protocol_mismatch_goal(Goal))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = error(
        bench_record_protocol_mismatch(_, warmup, 0, 1), _)).

test(normal_checks_reject_duplicate_measured_ids,
     [forall(duplicate_check_goal(Tool, Goal))]) :-
    catch((call(Goal), Result = succeeded), Error, Result = threw(Error)),
    assertion(Result = threw(Caught)),
    assertion(Caught = error(bench_duplicate_ids(Tool, ["core","core"]), _)).

test(exact_mode_requires_every_core_and_heavy_row) :-
    baseline_doc(['core.pl'-spec(core, inf, 100),
                  'heavy.pl'-spec(heavy, inf, 200)], Baseline),
    result_row('core.pl', core, inf, 100, CoreRow),
    with_output_to(string(_),
                   check_arrange_baseline:do_exact_check(
                       Baseline, _{swi_prolog:"10.1.10", results:[CoreRow]}, Fails)),
    assertion(Fails =:= 1).

test(exact_mode_forces_full_ladder) :-
    checker_mode(['--exact'], exact, []),
    exact_runner_args([], ['--heavy']).

test(exact_mode_rejects_partial_selection,
     [throws(error(bench_exact_args(['--fixture', core]), _))]) :-
    exact_runner_args(['--fixture', core], _).

test(persistence_rejects_sampling_overrides,
     [forall(member(Args,
                    [['--warmup', '0'], ['--iterations', '1'],
                     ['--warmup=0'], ['--iterations=1']]))]) :-
    catch((require_persistence_args(Args), Result = succeeded),
          Error, Result = threw(Error)),
    assertion(Result = threw(error(bench_persistence_sampling_override(_), _))).

test(workload_filter_is_semidet_with_repeated_substrings) :-
    workload_selected(a, false, banana, core),
    deterministic(Det),
    assertion(Det == true),
    findall(selected, workload_selected(a, false, banana, core), Selected),
    assertion(Selected == [selected]).

test(workload_tier_defaults_match_existing_runner_policy) :-
    workload_selected('', false, core_rung, core),
    \+ workload_selected('', false, heavy_rung, heavy),
    workload_selected('', true, heavy_rung, heavy),
    deterministic(Det),
    assertion(Det == true).

test(default_empty_selection_is_rejected,
     [throws(error(bench_empty_selection(arrange, ''), _))]) :-
    require_selected(arrange, '', []).

test(nonempty_selection_is_deterministic) :-
    require_selected(arrange, '', [workload]),
    deterministic(Det),
    assertion(Det == true).

test(common_runner_options_are_parsed_once) :-
    common_runner_spec(Spec),
    parse_runner_options(test, Spec,
                         ['--fixture', banana, '--heavy',
                          '--iterations', '3', '--warmup', '1'],
                         Options, Common),
    deterministic(Det),
    assertion(Det == true),
    assertion(memberchk(heavy(true), Options)),
    assertion(Common == runner_options{
        help:false, format:text, fixture:banana, heavy:true,
        iterations:3, warmup:1
    }).

test(runner_override_uses_manifest_or_explicit_value) :-
    apply_override(-1, 7, 7),
    apply_override(3, 7, 3),
    deterministic(Det),
    assertion(Det == true).

test(shared_paths_resolve_module_and_repository_files) :-
    benchmark_path('bench_cli.pl', BenchCli),
    deterministic(BenchmarkPathDet),
    assertion(BenchmarkPathDet == true),
    once(source_file(bench_cli:checker_mode(_, _, _), BenchCliSource)),
    assertion(BenchCli == BenchCliSource),
    repo_path('README.md', Readme),
    deterministic(RepoPathDet),
    assertion(RepoPathDet == true),
    assertion(exists_file(Readme)).

test(arrange_fixture_loader_matches_existing_fixture_terms) :-
    repo_path('fixtures/benchmark_08_words.pl', File),
    load_arrange_fixture(File, Words),
    deterministic(Det),
    assertion(Det == true),
    read_file_to_terms(File, Terms, []),
    memberchk(clues(Expected), Terms),
    assertion(Words == Expected).

test(arrange_fixture_loader_rejects_missing_clues,
     [throws(error(bench_fixture_missing_clues(_), _))]) :-
    with_temp_text("not_clues([]).\n", Path,
                   load_arrange_fixture(Path, _)).

test(arrange_fixture_loader_rejects_empty_clues,
     [throws(error(bench_fixture_empty_clues(_), _))]) :-
    with_temp_text("clues([]).\n", Path,
                   load_arrange_fixture(Path, _)).

test(benchmark_report_has_shared_envelope) :-
    Rows = [_{rung:one}],
    benchmark_report('crosswordsmith-test-bench', Rows, Doc),
    deterministic(ReportDet),
    assertion(ReportDet == true),
    swi_version(Version),
    deterministic(VersionDet),
    assertion(VersionDet == true),
    assertion(Doc =@= _{tool:'crosswordsmith-test-bench',
                        swi_prolog:Version, results:Rows}).

test(fill_exact_mode_checks_load_inferences) :-
    fill_baseline_doc(['core'-fill_spec(core, 100, 200, 300)], Baseline),
    fill_result_row('core', core, 100, 199, 300, Row),
    with_output_to(string(_),
                   check_fill_baseline:do_exact_check(
                       Baseline, _{swi_prolog:"10.1.10", results:[Row]}, Fails)),
    assertion(Fails =:= 1).

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

test(process_capture_modes_are_deterministic,
     [forall(member(Mode, [capture, inherit, null]))]) :-
    capture_process(path(swipl),
                    ['-q', '-g', 'format(ok)', '-t', halt], Mode,
                    Stdout, Stderr, Status),
    deterministic(Det),
    assertion(Det == true),
    assertion(Status == exit(0)),
    assertion(Stdout == "ok"),
    assertion(Stderr == "").

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

test(store_callbacks_resolve_and_retain_unmeasured_rows) :-
    Baseline = _{host:"old", swi_prolog:"10.1.10",
                 workloads:_{kept:_{value:7}, old:_{value:1, tier:"core"}}},
    Run = _{swi_prolog:"10.1.10",
            results:[_{id:"old", value:2, tier:"core"},
                     _{id:"new", value:3, tier:"heavy"}]},
    bench_store_test_caller:build_recorded(Baseline, Run, Recorded),
    deterministic(Det),
    assertion(Det == true),
    assertion(Recorded.workloads.kept.value =:= 7),
    assertion(Recorded.workloads.old.value =:= 2),
    assertion(Recorded.workloads.new.value =:= 3),
    assertion(Recorded.workloads.new.tier == "heavy").

test(store_rejects_duplicate_measured_keys,
     [throws(error(bench_record_duplicate_key(same), _))]) :-
    Baseline = _{host:"old", swi_prolog:"10.1.10", workloads:_{}},
    Run = _{swi_prolog:"10.1.10",
            results:[_{id:"same", value:1, tier:"core"},
                     _{id:"same", value:2, tier:"core"}]},
    bench_store_test_caller:build_recorded(Baseline, Run, _).

test(store_failed_verification_preserves_live_file) :-
    with_temp_json(_{value:"original"}, Path,
        ( catch(bench_store_test_caller:replace_rejected(
                    Path, _{value:"replacement"}),
                error(bench_record_verification_failed(Path), _), true),
          read_json_dict(Path, Written),
          get_dict(value, Written, Value),
          assertion(Value == "original"),
          atom_concat(Path, '.tmp', TempPath),
          assertion(\+ exists_file(TempPath)) )).

test(store_jsonl_preserves_string_whitespace) :-
    with_temp_text("", Path,
        ( with_output_to(string(_),
              append_history(Path, _{swi_prolog:"10.1.10"}, [],
                             [row-_{inf:1, note:"a  b"}])),
          read_history(Path, Entries),
          deterministic(Det),
          assertion(Det == true),
          Entries = [Entry],
          get_dict(rungs, Entry, Rungs),
          get_dict(row, Rungs, Cell),
          get_dict(note, Cell, Note),
          assertion(Note == "a  b") )).

test(store_jsonl_repairs_missing_final_newline) :-
    with_temp_text("{\"commit\":\"old\",\"date\":\"d\",\"dirty\":false,\"host\":\"h\",\"rungs\":{\"old\":{\"inf\":1}},\"swi\":\"10.1.10\",\"tiers\":\"core\"}", Path,
        ( with_output_to(string(_),
              append_history(Path, _{swi_prolog:"10.1.10"}, [],
                             [new-_{inf:2}])),
          read_history(Path, Entries),
          assertion(length(Entries, 2)) )).

test(store_jsonl_rejects_empty_rungs,
     [throws(error(bench_history_bad_rungs, _))]) :-
    with_temp_text("", Path,
        append_history(Path, _{swi_prolog:"10.1.10"}, [], [])).

test(store_history_reports_malformed_line_number) :-
    with_temp_text("{\"commit\":\"old\",\"date\":\"d\",\"dirty\":false,\"host\":\"h\",\"rungs\":{\"row\":{\"inf\":1}},\"swi\":\"10.1.10\",\"tiers\":\"core\"}\nnot-json\n", Path,
        ( catch(read_history(Path, _), Error, true),
          assertion(nonvar(Error)),
          assertion(Error = error(bench_history_line(Path, 2, _), _)) )).

test(store_history_reports_bad_envelope_line_number) :-
    with_temp_text("{\"rungs\":{}}\n", Path,
        ( catch(read_history(Path, _), Error, true),
          assertion(nonvar(Error)),
          assertion(Error = error(bench_history_line(Path, 1, _), _)) )).

test(store_history_renderer_emits_each_metric) :-
    Entries = [_{date:"2026-07-18T00:00:00", commit:"abc1234",
                 dirty:false, host:"host", swi:"10.1.10", tiers:core,
                 rungs:_{row:_{inf:10, load_inf:20}}}],
    with_output_to(string(Text),
                   render_history('test benchmark',
                                  [metric(search_inf, inf),
                                   metric(load_inf, load_inf)], Entries)),
    deterministic(Det),
    assertion(Det == true),
    assertion(sub_string(Text, _, _, _, "test benchmark history")),
    assertion(sub_string(Text, _, _, _, "search_inf per rung")),
    assertion(sub_string(Text, _, _, _, "load_inf per rung")).

test(store_history_partitions_swi_versions) :-
    Entries = [_{date:"2026-07-18T00:00:00", commit:"old",
                  dirty:false, host:"host", swi:"10.1.10", tiers:core,
                  rungs:_{row:_{inf:10}}},
               _{date:"2026-07-19T00:00:00", commit:"new",
                  dirty:false, host:"host", swi:"10.2.0", tiers:core,
                  rungs:_{row:_{inf:20}}}],
    with_output_to(string(Text),
                   render_history('test benchmark',
                                  [metric(search_inf, inf)], Entries)),
    assertion(sub_string(Text, _, _, _, "SWI-Prolog 10.1.10 trend")),
    assertion(sub_string(Text, _, _, _, "SWI-Prolog 10.2.0 trend")),
    assertion(\+ sub_string(Text, _, _, _, "+100.00%")).

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

test(checkers_reject_persistence_sampling_overrides_before_measurement,
     [forall(checker_script(Script))]) :-
    capture_process(path(swipl),
                    ['-q', Script, '--record', '--warmup', '0',
                     '--fixture', 'definitely-no-such-workload'], capture,
                    _Stdout, Stderr, Status),
    assertion(Status == exit(2)),
    assertion(sub_string(Stderr, _, _, _,
                         "persistence requires manifest sampling")).

test(checker_history_commands_render,
     [forall(checker_script(Script))]) :-
    capture_process(path(swipl), ['-q', Script, '--history'], capture,
                    Stdout, _Stderr, Status),
    assertion(Status == exit(0)),
    assertion(Stdout \== "").

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

test(fill_identity_partial_enumeration_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_partial_fill_identity_fixture(
            ManifestPath, WorkloadsPath, CliPath, Original),
        ( format(atom(ManifestEnv), 'FILL_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'FILL_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(CliEnv), 'FILL_IDENTITY_CLI=~w', [CliPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, CliEnv,
                           'benchmarks/check_fill_identity.sh', '--record'],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(2)),
          read_file_to_string(ManifestPath, Written, []),
          assertion(Written == Original) ),
        maplist(delete_if_exists, [ManifestPath, WorkloadsPath, CliPath])).

test(arrange_identity_partial_enumeration_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_partial_arrange_identity_fixture(
            ManifestPath, WorkloadsPath, RunnerPath, Original),
        ( format(atom(ManifestEnv), 'ARRANGE_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'ARRANGE_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(RunnerEnv), 'ARRANGE_IDENTITY_RUNNER=~w', [RunnerPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, RunnerEnv,
                           'benchmarks/check_arrange_identity.sh',
                           '--record', '--heavy'],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(2)),
          read_file_to_string(ManifestPath, Written, []),
          assertion(Written == Original) ),
        maplist(delete_if_exists,
                [ManifestPath, WorkloadsPath, RunnerPath])).

test(fill_identity_duplicate_ids_leave_manifest_unchanged) :-
    setup_call_cleanup(
        make_duplicate_fill_identity_fixture(
            ManifestPath, WorkloadsPath, CliPath, Original),
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

test(arrange_identity_duplicate_ids_leave_manifest_unchanged) :-
    setup_call_cleanup(
        make_duplicate_arrange_identity_fixture(
            ManifestPath, WorkloadsPath, RunnerPath, Original),
        ( format(atom(ManifestEnv), 'ARRANGE_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'ARRANGE_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(RunnerEnv), 'ARRANGE_IDENTITY_RUNNER=~w', [RunnerPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, RunnerEnv,
                           'benchmarks/check_arrange_identity.sh',
                           '--record', '--heavy'],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(1)),
          read_file_to_string(ManifestPath, Written, []),
          assertion(Written == Original) ),
        maplist(delete_if_exists,
                [ManifestPath, WorkloadsPath, RunnerPath])).

test(arrange_identity_rejects_duplicate_manifest_entries) :-
    setup_call_cleanup(
        make_single_arrange_identity_fixture(
            ManifestPath, WorkloadsPath, RunnerPath, _Original),
        ( write_duplicate_manifest(ManifestPath, 'dup.pl'),
          format(atom(ManifestEnv), 'ARRANGE_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'ARRANGE_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(RunnerEnv), 'ARRANGE_IDENTITY_RUNNER=~w', [RunnerPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, RunnerEnv,
                           'benchmarks/check_arrange_identity.sh'],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(2)) ),
        maplist(delete_if_exists,
                [ManifestPath, WorkloadsPath, RunnerPath])).

test(fill_identity_scripts_reject_duplicate_manifest_entries,
     [forall(member(Script,
                    ['benchmarks/check_fill_identity.sh',
                     'benchmarks/check_fill_identity_artifact.sh']))]) :-
    setup_call_cleanup(
        make_single_fill_identity_fixture(
            ManifestPath, WorkloadsPath, CliPath, _Original),
        ( write_duplicate_manifest(ManifestPath, ok),
          format(atom(ManifestEnv), 'FILL_IDENTITY_MANIFEST=~w', [ManifestPath]),
          format(atom(WorkloadsEnv), 'FILL_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
          format(atom(CliEnv), 'FILL_IDENTITY_CLI=~w', [CliPath]),
          capture_process(path(env),
                          [ManifestEnv, WorkloadsEnv, CliEnv, Script],
                          capture, _Stdout, _Stderr, Status),
          assertion(Status == exit(2)) ),
        maplist(delete_if_exists, [ManifestPath, WorkloadsPath, CliPath])).

test(fill_identity_write_failure_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_single_fill_identity_fixture(
            ManifestPath, WorkloadsPath, CliPath, Original),
        setup_call_cleanup(
            make_bash_override(printf_failure, BashEnvPath),
            ( format(atom(ManifestEnv), 'FILL_IDENTITY_MANIFEST=~w', [ManifestPath]),
              format(atom(WorkloadsEnv), 'FILL_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
              format(atom(CliEnv), 'FILL_IDENTITY_CLI=~w', [CliPath]),
              format(atom(BashEnv), 'BASH_ENV=~w', [BashEnvPath]),
              capture_process(path(env),
                              [BashEnv, ManifestEnv, WorkloadsEnv, CliEnv,
                               'benchmarks/check_fill_identity.sh', '--record'],
                              capture, _Stdout, _Stderr, Status),
              assertion(Status == exit(1)),
              read_file_to_string(ManifestPath, Written, []),
              assertion(Written == Original) ),
            delete_if_exists(BashEnvPath)),
        maplist(delete_if_exists, [ManifestPath, WorkloadsPath, CliPath])).

test(arrange_identity_write_failure_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_single_arrange_identity_fixture(
            ManifestPath, WorkloadsPath, RunnerPath, Original),
        setup_call_cleanup(
            make_bash_override(printf_failure, BashEnvPath),
            ( format(atom(ManifestEnv), 'ARRANGE_IDENTITY_MANIFEST=~w', [ManifestPath]),
              format(atom(WorkloadsEnv), 'ARRANGE_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
              format(atom(RunnerEnv), 'ARRANGE_IDENTITY_RUNNER=~w', [RunnerPath]),
              format(atom(BashEnv), 'BASH_ENV=~w', [BashEnvPath]),
              capture_process(path(env),
                              [BashEnv, ManifestEnv, WorkloadsEnv, RunnerEnv,
                               'benchmarks/check_arrange_identity.sh',
                               '--record', '--heavy'],
                              capture, _Stdout, _Stderr, Status),
              assertion(Status == exit(1)),
              read_file_to_string(ManifestPath, Written, []),
              assertion(Written == Original) ),
            delete_if_exists(BashEnvPath)),
        maplist(delete_if_exists,
                [ManifestPath, WorkloadsPath, RunnerPath])).

test(arrange_identity_hash_failure_leaves_manifest_unchanged) :-
    setup_call_cleanup(
        make_single_arrange_identity_fixture(
            ManifestPath, WorkloadsPath, RunnerPath, Original),
        setup_call_cleanup(
            make_bash_override(sha256_failure, BashEnvPath),
            ( format(atom(ManifestEnv), 'ARRANGE_IDENTITY_MANIFEST=~w', [ManifestPath]),
              format(atom(WorkloadsEnv), 'ARRANGE_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
              format(atom(RunnerEnv), 'ARRANGE_IDENTITY_RUNNER=~w', [RunnerPath]),
              format(atom(BashEnv), 'BASH_ENV=~w', [BashEnvPath]),
              capture_process(path(env),
                              [BashEnv, ManifestEnv, WorkloadsEnv, RunnerEnv,
                               'benchmarks/check_arrange_identity.sh',
                               '--record', '--heavy'],
                              capture, _Stdout, _Stderr, Status),
              assertion(Status == exit(1)),
              read_file_to_string(ManifestPath, Written, []),
              assertion(Written == Original) ),
            delete_if_exists(BashEnvPath)),
        maplist(delete_if_exists,
                [ManifestPath, WorkloadsPath, RunnerPath])).

test(fill_artifact_identity_scratch_failure_is_closed) :-
    setup_call_cleanup(
        make_single_fill_identity_fixture(
            ManifestPath, WorkloadsPath, CliPath, Original),
        setup_call_cleanup(
            make_bash_override(mktemp_failure, BashEnvPath),
            ( format(atom(ManifestEnv), 'FILL_IDENTITY_MANIFEST=~w', [ManifestPath]),
              format(atom(WorkloadsEnv), 'FILL_IDENTITY_WORKLOADS=~w', [WorkloadsPath]),
              format(atom(CliEnv), 'FILL_IDENTITY_CLI=~w', [CliPath]),
              format(atom(BashEnv), 'BASH_ENV=~w', [BashEnvPath]),
              capture_process(path(env),
                              [BashEnv, ManifestEnv, WorkloadsEnv, CliEnv,
                               'benchmarks/check_fill_identity_artifact.sh'],
                              capture, _Stdout, _Stderr, Status),
              assertion(Status == exit(2)),
              read_file_to_string(ManifestPath, Written, []),
              assertion(Written == Original) ),
            delete_if_exists(BashEnvPath)),
        maplist(delete_if_exists, [ManifestPath, WorkloadsPath, CliPath])).

:- end_tests(arrange_benchmark_promotion).

baseline_doc(Pairs0, _{host:"test-host", swi_prolog:"10.1.10",
                       regression_tolerance_pct:0.5, workloads:Workloads}) :-
    maplist(baseline_pair, Pairs0, Pairs),
    dict_pairs(Workloads, workloads, Pairs).

baseline_pair(Name-spec(Tier, Gate, Inf), Name-Spec) :-
    Spec = _{search_inf:Inf, cmd_wall_med_ms:10.0, cmd_rss_med_kib:12000,
              tier:Tier, gate:Gate, mode:size, expected:placed,
              grid:17, words:16, iterations:1,
              warmup:0, budget:500000000,
             info_only:["cmd_wall_med_ms", "cmd_rss_med_kib"]}.

result_row(Fixture, Tier, Gate, Inf,
           _{fixture:Fixture, search_inf_med:Inf, cmd_wall_med_ms:12.5,
              cmd_rss_med_kib:12345, tier:Tier, gate:Gate, mode:size,
              expected:placed, size:17, words:16,
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
              grid_file:"grids/test.json", dict_file:"dicts/test.txt",
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

exact_arrange_fails(Measured, Swi, Fails) :-
    baseline_doc(['core.pl'-spec(core, inf, 100)], Baseline),
    result_row('core.pl', core, inf, Measured, Row),
    with_output_to(string(_),
                    check_arrange_baseline:do_exact_check(
                        Baseline, _{swi_prolog:Swi, results:[Row]}, Fails)).

cross_version_classifier(
    check_arrange_baseline:classify(-10.0, 0.5, false, Kind), Kind).
cross_version_classifier(
    check_fill_baseline:classify(-10.0, 0.5, false, Kind), Kind).
cross_version_classifier(
    check_greedy_baseline:classify(-10.0, 0.5, false, Kind), Kind).

promotion_version_goal(
    check_arrange_baseline:require_promotion_version(
        _{swi_prolog:"10.1.10"}, _{swi_prolog:"10.2.0"})).
promotion_version_goal(
    check_fill_baseline:require_promotion_version(
        _{swi_prolog:"10.1.10"}, _{swi_prolog:"10.2.0"})).
promotion_version_goal(
    check_greedy_baseline:require_promotion_version(
        _{swi_prolog:"10.1.10"}, _{swi_prolog:"10.2.0"})).

record_version_goal(Completeness,
                    check_arrange_baseline:require_record_version(Baseline, Doc)) :-
    arrange_migration_fixture(Completeness, Baseline, Doc).
record_version_goal(Completeness,
                    check_fill_baseline:require_record_version(Baseline, Doc)) :-
    fill_migration_fixture(Completeness, Baseline, Doc).
record_version_goal(Completeness,
                    check_greedy_baseline:require_record_version(Baseline, Doc)) :-
    greedy_migration_fixture(Completeness, Baseline, Doc).

record_protocol_mismatch_goal(
        check_arrange_baseline:require_record_version(Baseline, Changed)) :-
    arrange_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
record_protocol_mismatch_goal(
        check_fill_baseline:require_record_version(Baseline, Changed)) :-
    fill_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
record_protocol_mismatch_goal(
        check_greedy_baseline:require_record_version(Baseline, Changed)) :-
    greedy_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).

comparison_protocol_mismatch_goal(
        check_arrange_baseline:do_check(Baseline, Changed, _, _)) :-
    arrange_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
comparison_protocol_mismatch_goal(
        check_arrange_baseline:do_exact_check(Baseline, Changed, _)) :-
    arrange_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
comparison_protocol_mismatch_goal(
        check_fill_baseline:do_check(Baseline, Changed, _, _)) :-
    fill_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
comparison_protocol_mismatch_goal(
        check_fill_baseline:do_exact_check(Baseline, Changed, _)) :-
    fill_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
comparison_protocol_mismatch_goal(
        check_greedy_baseline:do_check(Baseline, Changed, _, _)) :-
    greedy_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).
comparison_protocol_mismatch_goal(
        check_greedy_baseline:do_exact_check(Baseline, Changed, _)) :-
    greedy_migration_fixture(complete, Baseline, Doc),
    protocol_changed_doc(Doc, Changed).

duplicate_check_goal(arrange,
        check_arrange_baseline:do_check(Baseline, DuplicateDoc, _, _)) :-
    arrange_migration_fixture(complete, Baseline, Doc),
    duplicate_first_row(Doc, DuplicateDoc).
duplicate_check_goal(fill,
        check_fill_baseline:do_check(Baseline, DuplicateDoc, _, _)) :-
    fill_migration_fixture(complete, Baseline, Doc),
    duplicate_first_row(Doc, DuplicateDoc).
duplicate_check_goal(greedy,
        check_greedy_baseline:do_check(Baseline, DuplicateDoc, _, _)) :-
    greedy_migration_fixture(complete, Baseline, Doc),
    duplicate_first_row(Doc, DuplicateDoc).

duplicate_first_row(Doc, DuplicateDoc) :-
    Doc.results = [First|_],
    DuplicateDoc = Doc.put(_{swi_prolog:"10.1.10", results:[First, First]}).

protocol_changed_doc(Doc, Changed) :-
    Doc.results = [First|Rest],
    ChangedFirst = First.put(warmup, 1),
    Changed = Doc.put(_{swi_prolog:"10.1.10",
                        results:[ChangedFirst|Rest]}).

arrange_migration_fixture(Completeness, Baseline, Doc) :-
    baseline_doc(['core'-spec(core, inf, 100),
                  'heavy'-spec(heavy, inf, 200)], Baseline),
    result_row(core, core, inf, 100, Core),
    result_row(heavy, heavy, inf, 200, Heavy),
    migration_rows(Completeness, Core, Heavy, Rows),
    Doc = _{swi_prolog:"10.2.0", results:Rows}.

fill_migration_fixture(Completeness, Baseline, Doc) :-
    fill_baseline_doc([core-fill_spec(core, 100, 200, 300),
                       heavy-fill_spec(heavy, 400, 500, 600)], Baseline),
    fill_result_row(core, core, 100, 200, 300, Core),
    fill_result_row(heavy, heavy, 400, 500, 600, Heavy),
    migration_rows(Completeness, Core, Heavy, Rows),
    Doc = _{swi_prolog:"10.2.0", results:Rows}.

greedy_migration_fixture(Completeness, Baseline, Doc) :-
    greedy_result_row(core, core, Core),
    greedy_result_row(heavy, heavy, Heavy),
    check_greedy_baseline:record_spec(Core, new, CoreSpec),
    check_greedy_baseline:record_spec(Heavy, new, HeavySpec),
    Baseline = _{swi_prolog:"10.1.10",
                 workloads:_{core:CoreSpec, heavy:HeavySpec}},
    migration_rows(Completeness, Core, Heavy, Rows),
    Doc = _{swi_prolog:"10.2.0", results:Rows}.

greedy_result_row(Rung, Tier,
                  _{rung:Rung, fixture:"fixtures/test.pl", size:15,
                    framing:"size", mode:_{kind:"best_effort"},
                    construction:_{seed_answer:"SEED", corner:"topleft_across",
                                   expected:_{outcome:"completed", placed:1,
                                              dropped:0}},
                    words:1, iterations:1, warmup:0, tier:Tier,
                    metrics:_{construction_inf_med:1, sweep_inf_med:2,
                              postprocess_inf_med:3,
                              command_wall_med_ms:4.0,
                              command_rss_med_kib:5}}).

migration_rows(complete, Core, Heavy, [Core, Heavy]).
migration_rows(incomplete, Core, _Heavy, [Core]).

runner_script('benchmarks/run_arrange.pl').
runner_script('benchmarks/run_fill.pl').
runner_script('benchmarks/run_arrange_greedy.pl').

invalid_runner_args([positional]).
invalid_runner_args(['--format', bogus]).
invalid_runner_args(['--iterations', '0']).
invalid_runner_args(['--warmup', '-2']).
invalid_runner_args(['--fixture', 'definitely-no-such-workload']).

common_runner_spec([
    [opt(help), type(boolean), default(false), longflags([help])],
    [opt(format), type(atom), default(text), longflags([format])],
    [opt(fixture), type(atom), default(''), longflags([fixture])],
    [opt(heavy), type(boolean), default(false), longflags([heavy])],
    [opt(iterations), type(integer), default(-1), longflags([iterations])],
    [opt(warmup), type(integer), default(-1), longflags([warmup])]
]).

standalone_subject_goal(
    "consult('load.pl'),use_module('benchmarks/subjects.pl'),bench_subjects:arrange_command_sampler('/bin/true',ignored,3,size,placed,S),is_dict(S)").
standalone_subject_goal(
    "consult('load.pl'),use_module('benchmarks/fill_subjects.pl'),fill_subjects:fill_command_sampler('/bin/true',grid,dict,none,filled,S),is_dict(S)").

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

make_partial_fill_identity_fixture(ManifestPath, WorkloadsPath, CliPath,
                                   Original) :-
    Original = "original\tmanifest\n",
    tmp_file_stream(text, ManifestPath, ManifestStream),
    format(ManifestStream, '~s', [Original]),
    close(ManifestStream),
    tmp_file_stream(text, WorkloadsPath, WorkloadsStream),
    format(WorkloadsStream,
           "fill_workload(ok,'ok-grid',dict,none,1,0,filled,core,1).~n", []),
    format(WorkloadsStream,
           "fill_workload(_,_,_,_,_,_,_,_,_) :- throw(partial_enumeration).~n", []),
    close(WorkloadsStream),
    tmp_file_stream(text, CliPath, CliStream),
    format(CliStream, "#!/usr/bin/env bash~nprintf 'ok\\n'~n", []),
    close(CliStream),
    make_executable(CliPath).

make_partial_arrange_identity_fixture(ManifestPath, WorkloadsPath, RunnerPath,
                                      Original) :-
    Original = "original\tmanifest\n",
    tmp_file_stream(text, ManifestPath, ManifestStream),
    format(ManifestStream, '~s', [Original]),
    close(ManifestStream),
    tmp_file_stream(text, WorkloadsPath, WorkloadsStream),
    format(WorkloadsStream,
           "arrange_workload('fixtures/ok.pl',3,size,1,0,placed,core,inf,100,1).~n", []),
    format(WorkloadsStream,
           "arrange_workload(_,_,_,_,_,_,_,_,_,_) :- throw(partial_enumeration).~n", []),
    close(WorkloadsStream),
    tmp_file_stream(text, RunnerPath, RunnerStream),
    format(RunnerStream, ":- initialization(main, main).~n", []),
    format(RunnerStream, "main :- format(ok).~n", []),
    close(RunnerStream).

make_single_fill_identity_fixture(ManifestPath, WorkloadsPath, CliPath,
                                  Original) :-
    Original = "original\tmanifest\n",
    tmp_file_stream(text, ManifestPath, ManifestStream),
    format(ManifestStream, '~s', [Original]),
    close(ManifestStream),
    tmp_file_stream(text, WorkloadsPath, WorkloadsStream),
    format(WorkloadsStream,
           "fill_workload(ok,'ok-grid',dict,none,1,0,filled,core,1).~n", []),
    close(WorkloadsStream),
    tmp_file_stream(text, CliPath, CliStream),
    format(CliStream, "#!/usr/bin/env bash~nprintf 'ok\\n'~n", []),
    close(CliStream),
    make_executable(CliPath).

make_duplicate_fill_identity_fixture(ManifestPath, WorkloadsPath, CliPath,
                                     Original) :-
    make_single_fill_identity_fixture(
        ManifestPath, WorkloadsPath, CliPath, Original),
    setup_call_cleanup(
        open(WorkloadsPath, append, Stream),
        format(Stream,
               "fill_workload(ok,'other-grid',dict,none,1,0,filled,core,1).~n",
               []),
        close(Stream)).

make_single_arrange_identity_fixture(ManifestPath, WorkloadsPath, RunnerPath,
                                     Original) :-
    Original = "original\tmanifest\n",
    tmp_file_stream(text, ManifestPath, ManifestStream),
    format(ManifestStream, '~s', [Original]),
    close(ManifestStream),
    tmp_file_stream(text, WorkloadsPath, WorkloadsStream),
    format(WorkloadsStream,
           "arrange_workload('one/dup.pl',3,size,1,0,placed,core,inf,100,1).~n",
           []),
    close(WorkloadsStream),
    tmp_file_stream(text, RunnerPath, RunnerStream),
    format(RunnerStream, ":- initialization(main, main).~n", []),
    format(RunnerStream, "main :- format(ok).~n", []),
    close(RunnerStream).

make_duplicate_arrange_identity_fixture(ManifestPath, WorkloadsPath,
                                        RunnerPath, Original) :-
    make_single_arrange_identity_fixture(
        ManifestPath, WorkloadsPath, RunnerPath, Original),
    setup_call_cleanup(
        open(WorkloadsPath, append, Stream),
        format(Stream,
               "arrange_workload('two/dup.pl',3,size,1,0,placed,heavy,inf,100,1).~n",
               []),
        close(Stream)).

make_bash_override(printf_failure, Path) :-
    tmp_file_stream(text, Path, Stream),
    format(Stream, "printf() {~n", []),
    format(Stream, "    if [ \"$1\" = '%s\\t%s\\n' ]; then return 1; fi~n", []),
    format(Stream, "    builtin printf \"$@\"~n", []),
    format(Stream, "}~n", []),
    close(Stream).
make_bash_override(mktemp_failure, Path) :-
    tmp_file_stream(text, Path, Stream),
    format(Stream, "mktemp() { return 1; }~n", []),
    close(Stream).
make_bash_override(sha256_failure, Path) :-
    tmp_file_stream(text, Path, Stream),
    format(Stream, "sha256sum() {~n", []),
    format(Stream, "    builtin printf '0000000000000000000000000000000000000000000000000000000000000000  %s\\n' \"$1\"~n", []),
    format(Stream, "    return 1~n", []),
    format(Stream, "}~n", []),
    close(Stream).

write_duplicate_manifest(Path, Id) :-
    setup_call_cleanup(
        open(Path, write, Stream),
        ( format(Stream,
                 '~w\t0000000000000000000000000000000000000000000000000000000000000000~n',
                 [Id]),
          format(Stream, '~w\t', [Id]) ),
        close(Stream)).

make_executable(Path) :-
    process_create(path(chmod), ['u+x', Path], [process(PID)]),
    process_wait(PID, exit(0)).

delete_if_exists(Path) :-
    ( exists_file(Path) -> delete_file(Path) ; true ).

with_temp_json(Dict, Path, Goal) :-
    tmp_file_stream(text, Path, Stream),
    json_write_dict(Stream, Dict),
    close(Stream),
    setup_call_cleanup(true, Goal, delete_if_exists(Path)).

with_temp_text(Text, Path, Goal) :-
    tmp_file_stream(text, Path, Stream),
    format(Stream, '~s', [Text]),
    close(Stream),
    setup_call_cleanup(true, Goal, delete_if_exists(Path)).
