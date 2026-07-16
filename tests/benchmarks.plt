% Focused arrange baseline promotion safety tests. The recorder is exercised
% only against temporary baseline/history files; committed benchmark state is
% never opened for writing.

:- use_module(library(plunit)).
:- use_module(library(json)).
:- use_module(library(readutil)).
:- use_module(library(apply), [exclude/3]).
:- use_module('../benchmarks/check_baseline.pl', []).

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
