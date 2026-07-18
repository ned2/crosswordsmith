#!/usr/bin/env swipl
% benchmarks/check_baseline.pl - performance ratchet for the strict arrange ladder.
%
% Runs the product bench (benchmarks/run_arrange.pl --format json) and diffs the
% measured search-inference count of every rung against benchmarks/baseline.json,
% reporting the delta so you can HILL-CLIMB the arrange algorithm:
%
%   * a DROP past the tolerance is a WIN (print it, prompt `make bench-record`);
%   * a RISE past regression_tolerance_pct is a REGRESSION -> the check fails (exit 1);
%   * anything within +/-tolerance is unchanged.
%
% search inferences are deterministic and machine-independent (the SAME count native
% or under WASM), so the delta is a real algorithm signal, not machine noise. wall
% and rss are printed as informational deltas but NEVER fail the check. If the
% running SWI-Prolog differs from the baseline's, a regression is downgraded to a
% WARN (the version, not the code, moved the count) - regenerate the baseline.
%
% Modes:
%   swipl -q benchmarks/check_baseline.pl            % CHECK: diff + PASS/FAIL (exit 0/1)
%   swipl -q benchmarks/check_baseline.pl --exact    % EXACT: same-SWI equality over the
%                                                    %   complete core + heavy ladder
%   swipl -q benchmarks/check_baseline.pl --record   % RECORD: ratchet baseline.json to
%                                                    %   the measured numbers (accept wins),
%                                                    %   and append the run to history.jsonl
%   swipl -q benchmarks/check_baseline.pl --log      % LOG: diff + append to history.jsonl,
%                                                    %   WITHOUT moving the baseline
%   swipl -q benchmarks/check_baseline.pl --promote  % PROMOTE: check one run, then if clean
%                                                    %   record that exact run + history
%   swipl -q benchmarks/check_baseline.pl --history  % HISTORY: render the recorded trend
%                                                    %   (no bench run)
% Extra args pass through to run_arrange, so add --heavy to include the tail rungs:
%   make bench-check  BENCH_ARGS=--heavy
%   make bench-record BENCH_ARGS=--heavy
%
% The baseline is a MOVING reference (record ratchets it down), so it can't show a
% trajectory. benchmarks/history.jsonl is the append-only, git-tracked ledger: one
% JSON line per recorded/logged run, stamped with the git commit + timestamp, so the
% per-rung search_inf is comparable over time regardless of how the baseline moved.

:- module(check_arrange_baseline, []).

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(json)).
:- use_module(library(pairs), [pairs_keys/2]).
:- use_module('bench_process.pl', [capture_process/6]).
:- use_module('bench_cli.pl', [checker_mode/3, exact_runner_args/2]).
:- use_module('bench_exact.pl',
              [exact_version/4, exact_metric/4, exact_presence/5]).
:- use_module('bench_store.pl',
              [ read_json_dict/2,
                build_recorded_baseline/5,
                replace_json_dict/4,
                append_history/4,
                read_history/2,
                render_history/3,
                current_host/1
              ]).

% directory_file_path/3 is autoload-only (library(filesex)); explicit so this
% root also runs under autoload(false) (P11/C5, matching load.pl).
:- use_module(library(filesex), [directory_file_path/3]).

:- dynamic bench_dir/1.
:- prolog_load_context(directory, D), asserta(bench_dir(D)).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv0),
    catch(checker_mode(Argv0, Mode, Extra),
          E, (print_message(error, E), halt(2))),
    bench_dir(BenchDir),
    ( Mode == history -> show_history(BenchDir), halt(0) ; true ),
    directory_file_path(BenchDir, 'baseline.json', BaselinePath),
    history_path(BenchDir, HistoryPath),
    load_baseline(BaselinePath, Baseline),
    run_note(Extra, Note),
    ( Mode == record
    ->  format("crosswordsmith arrange - RECORD baseline~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_record(BaselinePath, Baseline, Doc),
        append_history(BenchDir, Doc, Extra),
        halt(0)
    ; Mode == log
    ->  format("crosswordsmith arrange - LOG run to history~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_check(Baseline, Doc, _Fails, _Wins),
        append_history(BenchDir, Doc, Extra),
        halt(0)
    ; Mode == promote
    ->  format("crosswordsmith arrange - PROMOTE (check, then ratchet this run)~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        promote_doc(BaselinePath, HistoryPath, Baseline, Doc, Extra, Fails, Wins),
        report_result(Fails, Wins),
        ( Fails =:= 0
        ->  format("~npromote: clean run recorded from the single measurement~n"),
            halt(0)
        ;   format("~npromote: regressions present - baseline and history NOT changed~n"),
            halt(1) )
    ; Mode == exact
    ->  catch(exact_runner_args(Extra, ExactArgs),
              E, (print_message(error, E), halt(2))),
        format("crosswordsmith arrange - exact inference check (+heavy tail)~n~n"),
        run_product_bench(BenchDir, ExactArgs, Doc),
        do_exact_check(Baseline, Doc, Fails),
        report_exact_result(Fails),
        ( Fails =:= 0 -> halt(0) ; halt(1) )
    ;   format("crosswordsmith arrange - performance ratchet~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_check(Baseline, Doc, Fails, Wins),
        report_result(Fails, Wins),
        ( Fails =:= 0 -> halt(0) ; halt(1) ) ).

run_note(Extra, " (+heavy tail)") :- memberchk('--heavy', Extra), !.
run_note(_, "").

% --- load + run --------------------------------------------------------------

load_baseline(Path, Baseline) :-
    ( exists_file(Path) -> true ; throw(error(baseline_missing(Path), _)) ),
    read_json_dict(Path, Baseline).

% Spawn run_arrange.pl, capture its JSON stdout, parse it. Extra args (e.g. --heavy)
% pass through. stderr flows to ours; we gate on the child exit code.
run_product_bench(BenchDir, Extra, Doc) :-
    directory_file_path(BenchDir, 'run_arrange.pl', RunArrange),
    append(['-q', RunArrange, '--', '--format', json], Extra, Args),
    capture_process(path(swipl), Args, inherit, JsonText, _Stderr, Status),
    ( Status == exit(0) -> true ; throw(error(bench_run_failed(Status), _)) ),
    atom_json_dict(JsonText, Doc, [default_tag(json)]).

% --- CHECK -------------------------------------------------------------------

do_check(Baseline, Doc, Fails, Wins) :-
    get_dict(results, Doc, Results),
    get_dict(swi_prolog, Doc, RunSwi),
    baseline_meta(Baseline, BaseHost, BaseSwi, Tol),
    current_host(RunHost),
    ( same_text(RunSwi, BaseSwi) -> VMatch = true ; VMatch = false ),
    ( same_text(RunHost, BaseHost) -> HMatch = true ; HMatch = false ),
    print_env(BaseHost, RunHost, HMatch, BaseSwi, RunSwi, VMatch, Tol),
    get_dict(workloads, Baseline, WL),
    format("~w~t~30|~t~w~14+~t~w~14+~t~w~11+   ~w~n",
           ['rung', 'baseline', 'measured', 'delta', 'status']),
    foldl(check_row(WL, Tol, VMatch), Results, 0-0, Fails-Wins),
    report_unmeasured(WL, Results),
    info_section(WL, Results, HMatch).

do_exact_check(Baseline, Doc, Fails) :-
    get_dict(swi_prolog, Baseline, BaseSwi),
    get_dict(swi_prolog, Doc, RunSwi),
    exact_version(BaseSwi, RunSwi, VersionStatus, VersionFails),
    format("exact SWI: ~w -> ~w (~w)~n", [BaseSwi, RunSwi, VersionStatus]),
    get_dict(workloads, Baseline, WL),
    get_dict(results, Doc, Results),
    exact_row_presence(WL, Results, Missing, Unexpected, PresenceFails),
    format("~w~t~30|~t~w~14+~t~w~14+   ~w~n",
           ['rung', 'baseline', 'measured', 'exact status']),
    foldl(exact_row(WL), Results, 0, MetricFails),
    report_exact_presence(Missing, Unexpected),
    Fails is VersionFails + PresenceFails + MetricFails.

exact_row(WL, Row, F0, F1) :-
    get_dict(fixture, Row, Fixture),
    get_dict(search_inf_med, Row, Measured),
    (   find_baseline(WL, Fixture, Spec)
    ->  get_dict(search_inf, Spec, Base),
        (   baseline_latency_gated(Spec)
        ->  Status = not_gated, MetricFails = 0
        ;   exact_metric(Base, Measured, Status, MetricFails)
        ),
        format("~w~t~30|~t~D~14+~t~D~14+   ~w~n",
               [Fixture, Base, Measured, Status])
    ;   MetricFails = 0,
        format("~w~t~30|~t~w~14+~t~D~14+   no_reference~n",
               [Fixture, '(none)', Measured])
    ),
    F1 is F0 + MetricFails.

exact_row_presence(WL, Results, Missing, Unexpected, Fails) :-
    dict_pairs(WL, _, Pairs),
    pairs_keys(Pairs, ReferenceIds),
    findall(Fixture,
            ( member(Row, Results), get_dict(fixture, Row, Fixture) ),
            MeasuredIds),
    exact_presence(ReferenceIds, MeasuredIds, Missing, Unexpected, Fails).

report_exact_presence([], []) :- !.
report_exact_presence(Missing, Unexpected) :-
    ( Missing == [] -> true ; format("missing exact rows: ~w~n", [Missing]) ),
    ( Unexpected == [] -> true ; format("rows without reference: ~w~n", [Unexpected]) ).

report_exact_result(0) :-
    format("~nEXACT RESULT: PASS (all gated inferences identical)~n").
report_exact_result(Fails) :-
    Fails > 0,
    format("~nEXACT RESULT: FAIL (~d exact comparison failure(s))~n", [Fails]).

% Exact mode trusts the committed reference policy, not a measured row that the
% code under review produced.
baseline_latency_gated(Spec) :-
    get_dict(gate, Spec, Gate),
    same_text(Gate, latency).

check_row(WL, Tol, VMatch, Row, F0-W0, F1-W1) :-
    get_dict(fixture, Row, Fx),
    get_dict(search_inf_med, Row, Meas),
    ( find_baseline(WL, Fx, Spec)
    ->  get_dict(search_inf, Spec, Base),
        ( Base =:= 0 -> Delta = 0.0 ; Delta is (Meas - Base) / Base * 100.0 ),
        ( row_latency_gated(Row)
        ->  Kind = latency_info
        ;   classify(Delta, Tol, VMatch, Kind) ),
        kind_counts(Kind, DF, DW), F1 is F0 + DF, W1 is W0 + DW,
        signed(Delta, DStr), kind_label(Kind, Label),
        format("~w~t~30|~t~D~14+~t~D~14+~t~w%~11+   ~w~n",
               [Fx, Base, Meas, DStr, Label])
    ;   F1 = F0, W1 = W0,
        format("~w~t~30|~t~w~14+~t~D~14+~t~w~11+   ~w~n",
               [Fx, '(none)', Meas, '-', 'NEW (not in baseline)']) ).

% A gate:latency workload (workloads.pl) is measured for its wall-latency only:
% its search count is pinned to the budget constant (budget-saturated), so the
% inference delta is not an algorithm signal and can neither fail nor win the
% ratchet. The manifest (via the run_arrange row) is authoritative; rows from
% pre-gate JSON default to inference-gated.
row_latency_gated(Row) :-
    get_dict(gate, Row, Gate),
    same_text(Gate, latency).

% win: dropped past tolerance. ok: within +/-tolerance. regression: rose past it
% (a hard fail when the SWI version matches; a WARN when it differs - version, not
% code, likely moved the count).
classify(Delta, Tol, _, win)  :- Delta =< -Tol, !.
classify(Delta, Tol, _, ok)   :- Delta =<  Tol, !.
classify(_, _, true,  regression) :- !.
classify(_, _, false, regression_warn).

kind_counts(win,             0, 1).
kind_counts(ok,              0, 0).
kind_counts(regression,      1, 0).
kind_counts(regression_warn, 0, 0).
kind_counts(latency_info,    0, 0).

kind_label(win,             'WIN (improvement)').
kind_label(ok,              'ok').
kind_label(regression,      'REGRESSION').
kind_label(regression_warn, 'regress? (swi-ver)').
kind_label(latency_info,    'info (latency-only)').

% A rung in the baseline that this run did not measure (e.g. heavy rungs on a
% core-only run). Informational - never a failure.
report_unmeasured(WL, Results) :-
    dict_pairs(WL, _, Pairs),
    findall(K, ( member(K-_, Pairs), \+ ( member(R, Results), get_dict(fixture, R, F), same_text(F, K) ) ), Missing),
    ( Missing == [] -> true
    ; format("~nnot run this pass (add --heavy to include): ~w~n", [Missing]) ).

% --- info (wall/rss), never gates --------------------------------------------

info_section(WL, Results, HMatch) :-
    format("~ninfo (host-specific, reporting-only):~n"),
    format("~w~t~30|~w~t~28+~w~n", ['rung', 'cmd_wall_med_ms', 'cmd_rss_med_kib']),
    forall(member(Row, Results), info_row(WL, HMatch, Row)).

info_row(WL, _HMatch, Row) :-
    get_dict(fixture, Row, Fx),
    ( find_baseline(WL, Fx, Spec)
    ->  metric_delta(Spec, Row, cmd_wall_med_ms, WallStr),
        metric_delta(Spec, Row, cmd_rss_med_kib, RssStr),
        format("~w~t~30|~w~t~28+~w~n", [Fx, WallStr, RssStr])
    ;   format("~w~t~30|~w~n", [Fx, '(new)']) ).

metric_delta(Spec, Row, Key, Str) :-
    ( get_dict(Key, Spec, Base), get_dict(Key, Row, Meas)
    ->  ( Base =:= 0 -> Pct = 0.0 ; Pct is 100.0 * (Meas - Base) / Base ),
        signed(Pct, PStr),
        format(atom(Str), "~w -> ~w (~w%)", [Base, Meas, PStr])
    ;   Str = 'n/a' ).

% --- RECORD (ratchet the baseline to the measured numbers) -------------------

% Check the in-memory measurement before any write. On a clean result, the
% exact same Doc is recorded and logged; no second benchmark run can drift from
% the accepted measurement.
promote_doc(BaselinePath, HistoryPath, Baseline, Doc, Extra, Fails, Wins) :-
    do_check(Baseline, Doc, Fails, Wins),
    (   Fails =:= 0
    ->  write_record(BaselinePath, Baseline, Doc),
        append_history_to(HistoryPath, Doc, Extra)
    ;   true
    ).

do_record(BaselinePath, Baseline, Doc) :-
    write_record(BaselinePath, Baseline, Doc).

write_record(BaselinePath, Baseline, Doc) :-
    build_recorded_baseline(Baseline, Doc, Recorded),
    replace_json_dict(BaselinePath, Recorded, 90, verify_recorded_doc(Doc)),
    get_dict(results, Doc, Results),
    get_dict(workloads, Baseline, WL0),
    format("baseline updated and read-back verified: ~w~n~n", [BaselinePath]),
    forall(member(Row, Results), report_recorded(WL0, Row)),
    unmeasured_note(WL0, Results).

build_recorded_baseline(Baseline, Doc, Recorded) :-
    build_recorded_baseline(Baseline, Doc, record_key, record_spec, Recorded).

record_key(Row, Row.fixture).

record_spec(Row, existing(Old), Spec) :-
    Spec = Old.put(_{ search_inf:      Row.search_inf_med,
                      cmd_wall_med_ms: Row.cmd_wall_med_ms,
                      cmd_rss_med_kib: Row.cmd_rss_med_kib }).
record_spec(Row, new, Spec) :-
    new_rung_spec(Row, Spec).

verify_recorded_doc(Doc, Written) :-
    get_dict(results, Doc, Results),
    verify_recorded_results(Written, Results).

% A success banner is not evidence that persistence worked. Read the file back
% and verify every measured rung carries every value the writer records.
verify_recorded_results(Written, Results) :-
    get_dict(workloads, Written, Workloads),
    forall(member(Row, Results), verify_recorded_row(Workloads, Row)).

verify_recorded_row(Workloads, Row) :-
    get_dict(fixture, Row, Fixture),
    (   find_baseline(Workloads, Fixture, Spec)
    ->  verify_recorded_value(Fixture, search_inf, Spec, Row.search_inf_med),
        verify_recorded_value(Fixture, cmd_wall_med_ms, Spec, Row.cmd_wall_med_ms),
        verify_recorded_value(Fixture, cmd_rss_med_kib, Spec, Row.cmd_rss_med_kib)
    ;   throw(error(record_readback_missing(Fixture), _))
    ).

verify_recorded_value(Fixture, Key, Spec, Expected) :-
    (   get_dict(Key, Spec, Actual), Actual =:= Expected
    ->  true
    ;   throw(error(record_readback_mismatch(Fixture, Key, Expected), _))
    ).

% A measured rung with no baseline entry joins the baseline as a COMPLETE spec
% (the first --record after a rung is added to workloads.pl). The spec metadata
% (tier/warmup/budget/words) comes from the run_arrange result row.
new_rung_spec(Row, _{ search_inf:      Row.search_inf_med,
                      cmd_wall_med_ms: Row.cmd_wall_med_ms,
                      cmd_rss_med_kib: Row.cmd_rss_med_kib,
                      info_only:       ["cmd_wall_med_ms", "cmd_rss_med_kib"],
                      grid:            Row.size,
                      words:           Row.words,
                      tier:            Row.tier,
                      gate:            Row.get(gate, "inf"),
                      iterations:      Row.iterations,
                      warmup:          Row.warmup,
                      budget:          Row.budget }).

report_recorded(WL0, Row) :-
    get_dict(fixture, Row, Fx),
    get_dict(search_inf_med, Row, New),
    ( find_baseline(WL0, Fx, Spec), get_dict(search_inf, Spec, Old)
    ->  ( Old =:= 0 -> D = 0.0 ; D is (New - Old) / Old * 100.0 ), signed(D, DStr),
        format("  ~w~t~30|~D -> ~D  (~w%)~n", [Fx, Old, New, DStr])
    ;   format("  ~w~t~30|new rung, added to baseline at ~D~n", [Fx, New]) ).

unmeasured_note(WL0, Results) :-
    dict_pairs(WL0, _, Pairs),
    findall(K, ( member(K-_, Pairs), \+ result_for(Results, K, _) ), Kept),
    ( Kept == [] -> true
    ; format("~nkept unchanged (not measured this pass; use --heavy): ~w~n", [Kept]) ).

result_for(Results, Key, Row) :-
    member(Row, Results), get_dict(fixture, Row, F), same_text(F, Key), !.

% --- result banner -----------------------------------------------------------

report_result(Fails, Wins) :-
    ( Fails > 0
    ->  format("~nRESULT: FAIL  (~d regression(s), ~d win(s))~n", [Fails, Wins]),
        format("A gated rung rose past tolerance. If the change was intentional,~n"),
        format("re-baseline with `make bench-record` and review the diff.~n")
    ;   Wins > 0
    ->  format("~nRESULT: PASS  (~d improvement(s), 0 regressions)~n", [Wins]),
        format("Lock the win(s) in with `make bench-record`.~n")
    ;   format("~nRESULT: PASS  (no change; 0 regressions)~n") ).

% --- env banner --------------------------------------------------------------

print_env(BaseHost, RunHost, HMatch, BaseSwi, RunSwi, VMatch, Tol) :-
    ( HMatch == true
    ->  HostNote = 'same host -> wall/rss deltas meaningful'
    ;   HostNote = 'DIFFERENT host -> wall/rss NOT comparable (search_inf still portable)' ),
    ( VMatch == true
    ->  VerNote = 'same swi  -> regressions gate'
    ;   VerNote = 'DIFFERENT swi -> regressions downgraded to WARN (regenerate baseline)' ),
    format("baseline:  host ~w,  swi ~w~n", [BaseHost, BaseSwi]),
    format("this run:  host ~w,  swi ~w~n", [RunHost, RunSwi]),
    format("  ~w~n", [HostNote]),
    format("  ~w   (regression tolerance +~w%)~n~n", [VerNote, Tol]).

% --- helpers -----------------------------------------------------------------

baseline_meta(Baseline, Host, Swi, Tol) :-
    get_dict(host, Baseline, Host),
    get_dict(swi_prolog, Baseline, Swi),
    ( get_dict(regression_tolerance_pct, Baseline, Tol) -> true ; Tol = 0.5 ).

find_baseline(WL, Fixture, Spec) :-
    dict_pairs(WL, _, Pairs),
    member(K-Spec, Pairs), same_text(K, Fixture), !.

signed(X, S) :- ( X >= 0 -> format(atom(S), "+~2f", [X]) ; format(atom(S), "~2f", [X]) ).

same_text(A, B) :- text_to_string(A, S), text_to_string(B, S).

% --- HISTORY (append-only trend ledger) --------------------------------------
% One JSON object per line in benchmarks/history.jsonl. search_inf is the portable,
% comparable-over-time number; wall/rss ride along for same-host trend but stay
% informational. Each entry is stamped with the git commit + local timestamp so a
% row can be attributed to the change that produced it.

history_path(BenchDir, Path) :-
    directory_file_path(BenchDir, 'history.jsonl', Path).

append_history(BenchDir, Doc, Extra) :-
    history_path(BenchDir, Path),
    append_history_to(Path, Doc, Extra).

append_history_to(Path, Doc, Extra) :-
    get_dict(results, Doc, Results),
    findall(F-Cell,
            ( member(Row, Results),
              get_dict(fixture, Row, F0), atom_string(F, F0),
              Cell = _{ inf:     Row.search_inf_med,
                        wall_ms: Row.cmd_wall_med_ms,
                         rss_kib: Row.cmd_rss_med_kib } ),
            RungPairs),
    append_history(Path, Doc, Extra, RungPairs).

show_history(BenchDir) :-
    history_path(BenchDir, Path),
    read_history(Path, Entries),
    ( Entries == []
    ->  format("no benchmark history yet.~n"),
        format("record one with:  make bench-record   (or  make bench-log)~n")
    ;   render_history('benchmark', [metric(search_inf, inf)], Entries) ).

:- multifile prolog:error_message//1.
prolog:error_message(baseline_missing(Path)) -->
    [ 'check_baseline: ~w not found (regenerate with: make bench-record)'-[Path] ].
prolog:error_message(bench_run_failed(Status)) -->
    [ 'check_baseline: the product bench (run_arrange.pl) failed: ~q'-[Status] ].
prolog:error_message(record_readback_missing(Fixture)) -->
    [ 'check_baseline: recorded rung ~w is absent after baseline read-back'-[Fixture] ].
prolog:error_message(record_readback_mismatch(Fixture, Key, Expected)) -->
    [ 'check_baseline: recorded rung ~w field ~w does not equal measured value ~w after read-back'-[Fixture, Key, Expected] ].
