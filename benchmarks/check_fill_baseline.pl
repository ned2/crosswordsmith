#!/usr/bin/env swipl
% benchmarks/check_fill_baseline.pl - performance ratchet for the FILL ladder.
%
% Forked from check_baseline.pl (the arrange ratchet) per the campaign plan's
% "fork, don't parameterize" call: the two ratchets share a schema shape but
% gate different rung semantics and must never couple. This file owns
% fill_baseline.json / fill_history.jsonl; the arrange files are untouched.
%
% TWO measured layers per rung:
%
%   search_inf - the fill_attempt/8 search-inference count. GATED: a rise past
%                regression_tolerance_pct (relative, default 0.5%) FAILS the
%                check (exit 1); a drop past it is a WIN. Deterministic and
%                machine-independent (same count native or under WASM).
%   load_inf   - the load_dict/3 dictionary-load inference count. GATED since
%                the F-L1 acceptance (2026-07-05, baseline recorded at 16452de):
%                Phase 3 took the scheduled promotion decision - load is the
%                dominant end-to-end term (P-F1) and now carries recorded wins
%                worth defending (F-L1 -59.7%), and F-H2's measured mask-
%                construction tax lands exactly here. Same tolerance and
%                swi-version-downgrade semantics as search_inf; either layer
%                regressing fails the check independently.
%
% wall/rss are printed as informational deltas and never fail the check. If the
% running SWI-Prolog differs from the baseline's, a search_inf regression is
% downgraded to a WARN (the version, not the code, moved the count).
%
% Modes:
%   swipl -q benchmarks/check_fill_baseline.pl            % CHECK: diff + PASS/FAIL (exit 0/1)
%   swipl -q benchmarks/check_fill_baseline.pl --exact    % EXACT: same-SWI equality over the
%                                                         %   complete core + heavy ladder
%   swipl -q benchmarks/check_fill_baseline.pl --record   % RECORD: ratchet fill_baseline.json,
%                                                         %   append to fill_history.jsonl
%   swipl -q benchmarks/check_fill_baseline.pl --log      % LOG: diff + append to history,
%                                                         %   WITHOUT moving the baseline
%   swipl -q benchmarks/check_fill_baseline.pl --promote  % PROMOTE: check, then if CLEAN
%                                                         %   (zero regressions) record THIS
%                                                         %   run as the new baseline + history
%                                                         %   entry - one measurement instead of
%                                                         %   the check-then-record double run
%   swipl -q benchmarks/check_fill_baseline.pl --history  % HISTORY: render the trend (no run)
% Extra args pass through to run_fill, so add --heavy for the tail rungs:
%   make bench-fill-check  BENCH_ARGS=--heavy
%
% NEW-RUNG RECORDING (the arrange campaign's silent-drop incident, carried fix):
% --record must persist rungs it has never seen. run_fill rows carry the full
% spec metadata (grid/dict/seeds/tier/warmup/budget/words/size), and
% new_rung_pairs/3 builds a COMPLETE baseline spec for any measured rung absent
% from the baseline - never silently dropping it. After any --record, READ THE
% FILE BACK and verify every rung is present (the tool's own success message is
% not evidence; the state is).

:- module(check_fill_baseline, []).

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
    directory_file_path(BenchDir, 'fill_baseline.json', BaselinePath),
    load_baseline(BaselinePath, Baseline),
    run_note(Extra, Note),
    ( Mode == record
    ->  format("crosswordsmith fill - RECORD baseline~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_record(BaselinePath, Baseline, Doc),
        append_history(BenchDir, Doc, Extra),
        halt(0)
    ; Mode == log
    ->  format("crosswordsmith fill - LOG run to history~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_check(Baseline, Doc, _Fails, _Wins),
        append_history(BenchDir, Doc, Extra),
        halt(0)
    ; Mode == promote
    ->  format("crosswordsmith fill - PROMOTE (check, then ratchet this run)~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_check(Baseline, Doc, Fails, Wins),
        report_result(Fails, Wins),
        ( Fails =:= 0
        ->  format("~npromote: clean run - recording it as the new baseline~n~n"),
            do_record(BaselinePath, Baseline, Doc),
            append_history(BenchDir, Doc, Extra),
            halt(0)
        ;   format("~npromote: regressions present - baseline NOT recorded~n"),
            halt(1) )
    ; Mode == exact
    ->  catch(exact_runner_args(Extra, ExactArgs),
              E, (print_message(error, E), halt(2))),
        format("crosswordsmith fill - exact inference check (+heavy tail)~n~n"),
        run_product_bench(BenchDir, ExactArgs, Doc),
        do_exact_check(Baseline, Doc, Fails),
        report_exact_result(Fails),
        ( Fails =:= 0 -> halt(0) ; halt(1) )
    ;   format("crosswordsmith fill - performance ratchet~s~n~n", [Note]),
        run_product_bench(BenchDir, Extra, Doc),
        do_check(Baseline, Doc, Fails, Wins),
        report_result(Fails, Wins),
        ( Fails =:= 0 -> halt(0) ; halt(1) ) ).

run_note(Extra, " (+heavy tail)") :- memberchk('--heavy', Extra), !.
run_note(_, "").

% --- load + run --------------------------------------------------------------

load_baseline(Path, Baseline) :-
    ( exists_file(Path) -> true ; throw(error(fill_baseline_missing(Path), _)) ),
    read_json_dict(Path, Baseline).

% Spawn run_fill.pl, capture its JSON stdout, parse it. Extra args (e.g. --heavy)
% pass through. stderr flows to ours; we gate on the child exit code.
run_product_bench(BenchDir, Extra, Doc) :-
    directory_file_path(BenchDir, 'run_fill.pl', RunFill),
    append(['-q', RunFill, '--', '--format', json], Extra, Args),
    capture_process(path(swipl), Args, inherit, JsonText, _Stderr, Status),
    ( Status == exit(0) -> true ; throw(error(fill_bench_run_failed(Status), _)) ),
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
    format("search_inf (GATED):~n"),
    format("~w~t~18|~t~w~14+~t~w~14+~t~w~11+   ~w~n",
           ['rung', 'baseline', 'measured', 'delta', 'status']),
    foldl(check_row(WL, Tol, VMatch), Results, 0-0, F0-W0),
    load_section(WL, Tol, VMatch, Results, F0-W0, Fails-Wins),
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
    format("~w~t~18|~w~t~20+~t~w~14+~t~w~14+   ~w~n",
           ['rung', 'metric', 'baseline', 'measured', 'exact status']),
    foldl(exact_fill_row(WL), Results, 0, MetricFails),
    report_exact_presence(Missing, Unexpected),
    Fails is VersionFails + PresenceFails + MetricFails.

exact_fill_row(WL, Row, F0, F2) :-
    get_dict(rung, Row, Rung),
    (   find_baseline(WL, Rung, Spec)
    ->  exact_fill_metric(Rung, search_inf, Spec.search_inf,
                          Row.search_inf_med, F0, F1),
        exact_fill_metric(Rung, load_inf, Spec.load_inf,
                          Row.load_inf, F1, F2)
    ;   F2 = F0,
        format("~w~t~18|~w~t~20+~t~w~14+~t~w~14+   no_reference~n",
               [Rung, gated_metrics, '(none)', '(n/a)'])
    ).

exact_fill_metric(Rung, Metric, Base, Measured, F0, F1) :-
    exact_metric(Base, Measured, Status, MetricFails),
    F1 is F0 + MetricFails,
    format("~w~t~18|~w~t~20+~t~D~14+~t~D~14+   ~w~n",
           [Rung, Metric, Base, Measured, Status]).

exact_row_presence(WL, Results, Missing, Unexpected, Fails) :-
    dict_pairs(WL, _, Pairs),
    pairs_keys(Pairs, ReferenceIds),
    findall(Rung, (member(Row, Results), get_dict(rung, Row, Rung)), MeasuredIds),
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

check_row(WL, Tol, VMatch, Row, F0-W0, F1-W1) :-
    get_dict(rung, Row, Rung),
    get_dict(search_inf_med, Row, Meas),
    ( find_baseline(WL, Rung, Spec)
    ->  get_dict(search_inf, Spec, Base),
        ( Base =:= 0 -> Delta = 0.0 ; Delta is (Meas - Base) / Base * 100.0 ),
        classify(Delta, Tol, VMatch, Kind),
        kind_counts(Kind, DF, DW), F1 is F0 + DF, W1 is W0 + DW,
        signed(Delta, DStr), kind_label(Kind, Label),
        format("~w~t~18|~t~D~14+~t~D~14+~t~w%~11+   ~w~n",
               [Rung, Base, Meas, DStr, Label])
    ;   F1 = F0, W1 = W0,
        format("~w~t~18|~t~w~14+~t~D~14+~t~w~11+   ~w~n",
               [Rung, '(none)', Meas, '-', 'NEW (not in baseline)']) ).

% win: dropped past tolerance. ok: within +/-tolerance. regression: rose past it
% (a hard fail when the SWI version matches; a WARN when it differs).
classify(Delta, Tol, _, win)  :- Delta =< -Tol, !.
classify(Delta, Tol, _, ok)   :- Delta =<  Tol, !.
classify(_, _, true,  regression) :- !.
classify(_, _, false, regression_warn).

kind_counts(win,             0, 1).
kind_counts(ok,              0, 0).
kind_counts(regression,      1, 0).
kind_counts(regression_warn, 0, 0).

kind_label(win,             'WIN (improvement)').
kind_label(ok,              'ok').
kind_label(regression,      'REGRESSION').
kind_label(regression_warn, 'regress? (swi-ver)').

% load_inf: GATED since the F-L1 acceptance (Phase 3's scheduled promotion
% decision) - same classify/tolerance machinery as search_inf, accumulating
% into the same fail/win totals so either layer defends independently.
load_section(WL, Tol, VMatch, Results, Acc0, Acc) :-
    format("~nload_inf (GATED; promoted at F-L1 acceptance):~n"),
    format("~w~t~18|~t~w~14+~t~w~14+~t~w~11+   ~w~n",
           ['rung', 'baseline', 'measured', 'delta', 'status']),
    foldl(load_row(WL, Tol, VMatch), Results, Acc0, Acc).

load_row(WL, Tol, VMatch, Row, F0-W0, F1-W1) :-
    get_dict(rung, Row, Rung),
    get_dict(load_inf, Row, Meas),
    ( find_baseline(WL, Rung, Spec), get_dict(load_inf, Spec, Base)
    ->  ( Base =:= 0 -> Delta = 0.0 ; Delta is (Meas - Base) / Base * 100.0 ),
        classify(Delta, Tol, VMatch, Kind),
        kind_counts(Kind, DF, DW), F1 is F0 + DF, W1 is W0 + DW,
        signed(Delta, DStr), kind_label(Kind, Label),
        format("~w~t~18|~t~D~14+~t~D~14+~t~w%~11+   ~w~n",
               [Rung, Base, Meas, DStr, Label])
    ;   F1 = F0, W1 = W0,
        format("~w~t~18|~t~w~14+~t~D~14+~t~w~11+   ~w~n",
               [Rung, '(none)', Meas, '-', 'NEW (not in baseline)']) ).

% A rung in the baseline that this run did not measure (e.g. heavy rungs on a
% core-only run). Informational - never a failure.
report_unmeasured(WL, Results) :-
    dict_pairs(WL, _, Pairs),
    findall(K, ( member(K-_, Pairs), \+ ( member(R, Results), get_dict(rung, R, F), same_text(F, K) ) ), Missing),
    ( Missing == [] -> true
    ; format("~nnot run this pass (add --heavy to include): ~w~n", [Missing]) ).

% --- info (wall/rss), never gates --------------------------------------------

info_section(WL, Results, _HMatch) :-
    format("~ninfo (host-specific, reporting-only):~n"),
    format("~w~t~18|~w~t~30+~w~n", ['rung', 'cmd_wall_med_ms', 'cmd_rss_med_kib']),
    forall(member(Row, Results), info_row(WL, Row)).

info_row(WL, Row) :-
    get_dict(rung, Row, Rung),
    ( find_baseline(WL, Rung, Spec)
    ->  metric_delta(Spec, Row, cmd_wall_med_ms, WallStr),
        metric_delta(Spec, Row, cmd_rss_med_kib, RssStr),
        format("~w~t~18|~w~t~30+~w~n", [Rung, WallStr, RssStr])
    ;   format("~w~t~18|~w~n", [Rung, '(new)']) ).

metric_delta(Spec, Row, Key, Str) :-
    ( get_dict(Key, Spec, Base), get_dict(Key, Row, Meas)
    ->  ( Base =:= 0 -> Pct = 0.0 ; Pct is 100.0 * (Meas - Base) / Base ),
        signed(Pct, PStr),
        format(atom(Str), "~w -> ~w (~w%)", [Base, Meas, PStr])
    ;   Str = 'n/a' ).

% --- RECORD (ratchet the baseline to the measured numbers) -------------------

do_record(BaselinePath, Baseline, Doc) :-
    build_recorded_baseline(Baseline, Doc, Recorded),
    replace_json_dict(BaselinePath, Recorded, 90, verify_recorded_doc(Doc)),
    get_dict(results, Doc, Results),
    get_dict(workloads, Baseline, WL0),
    format("baseline updated and read-back verified: ~w~n~n", [BaselinePath]),
    forall(member(Row, Results), report_recorded(WL0, Row)),
    unmeasured_note(WL0, Results).

build_recorded_baseline(Baseline, Doc, Recorded) :-
    build_recorded_baseline(Baseline, Doc, record_key, record_spec, Recorded).

record_key(Row, Row.rung).

record_spec(Row, existing(Old), Spec) :-
    Spec = Old.put(_{ search_inf:      Row.search_inf_med,
                      load_inf:        Row.load_inf,
                      grid_inf:        Row.grid_inf,
                      cmd_wall_med_ms: Row.cmd_wall_med_ms,
                      cmd_rss_med_kib: Row.cmd_rss_med_kib }).
record_spec(Row, new, Spec) :-
    new_rung_spec(Row, Spec).

verify_recorded_doc(Doc, Written) :-
    get_dict(results, Doc, Results),
    verify_recorded_results(Written, Results).

verify_recorded_results(Written, Results) :-
    get_dict(workloads, Written, Workloads),
    forall(member(Row, Results), verify_recorded_row(Workloads, Row)).

verify_recorded_row(Workloads, Row) :-
    get_dict(rung, Row, Rung),
    (   find_baseline(Workloads, Rung, Spec)
    ->  verify_recorded_value(Rung, search_inf, Spec, Row.search_inf_med),
        verify_recorded_value(Rung, load_inf, Spec, Row.load_inf),
        verify_recorded_value(Rung, grid_inf, Spec, Row.grid_inf),
        verify_recorded_value(Rung, cmd_wall_med_ms, Spec, Row.cmd_wall_med_ms),
        verify_recorded_value(Rung, cmd_rss_med_kib, Spec, Row.cmd_rss_med_kib)
    ;   throw(error(fill_record_readback_missing(Rung), _))
    ).

verify_recorded_value(Rung, Key, Spec, Expected) :-
    (   get_dict(Key, Spec, Actual), Actual =:= Expected
    ->  true
    ;   throw(error(fill_record_readback_mismatch(Rung, Key, Expected), _))
    ).

% A measured rung with no baseline entry joins the baseline as a COMPLETE spec
% (the first --record after a rung is added to fill_workloads.pl). The spec
% metadata comes from the run_fill result row - this is the carried fix for the
% arrange campaign's --record silent-drop bug.
new_rung_spec(Row, _{ search_inf:      Row.search_inf_med,
                      load_inf:        Row.load_inf,
                      grid_inf:        Row.grid_inf,
                      cmd_wall_med_ms: Row.cmd_wall_med_ms,
                      cmd_rss_med_kib: Row.cmd_rss_med_kib,
                      info_only:       ["grid_inf", "cmd_wall_med_ms", "cmd_rss_med_kib"],
                      grid_file:       Row.grid_file,
                      dict_file:       Row.dict_file,
                      seeds:           Row.seeds,
                      grid:            Row.size,
                      words:           Row.words,
                      tier:            Row.tier,
                      iterations:      Row.iterations,
                      warmup:          Row.warmup,
                      budget:          Row.budget }).

report_recorded(WL0, Row) :-
    get_dict(rung, Row, Rung),
    get_dict(search_inf_med, Row, New),
    ( find_baseline(WL0, Rung, Spec), get_dict(search_inf, Spec, Old)
    ->  ( Old =:= 0 -> D = 0.0 ; D is (New - Old) / Old * 100.0 ), signed(D, DStr),
        format("  ~w~t~18|~D -> ~D  (~w%)~n", [Rung, Old, New, DStr])
    ;   format("  ~w~t~18|new rung, added to baseline at ~D~n", [Rung, New]) ).

unmeasured_note(WL0, Results) :-
    dict_pairs(WL0, _, Pairs),
    findall(K, ( member(K-_, Pairs), \+ result_for(Results, K, _) ), Kept),
    ( Kept == [] -> true
    ; format("~nkept unchanged (not measured this pass; use --heavy): ~w~n", [Kept]) ).

result_for(Results, Key, Row) :-
    member(Row, Results), get_dict(rung, Row, F), same_text(F, Key), !.

% --- result banner -----------------------------------------------------------

report_result(Fails, Wins) :-
    ( Fails > 0
    ->  format("~nRESULT: FAIL  (~d regression(s), ~d win(s))~n", [Fails, Wins]),
        format("A gated rung's search_inf or load_inf rose past tolerance. If intentional,~n"),
        format("re-baseline with `make bench-fill-record` and review the diff.~n")
    ;   Wins > 0
    ->  format("~nRESULT: PASS  (~d improvement(s), 0 regressions)~n", [Wins]),
        format("Lock the win(s) in with `make bench-fill-record`.~n")
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

find_baseline(WL, Rung, Spec) :-
    dict_pairs(WL, _, Pairs),
    member(K-Spec, Pairs), same_text(K, Rung), !.

signed(X, S) :- ( X >= 0 -> format(atom(S), "+~2f", [X]) ; format(atom(S), "~2f", [X]) ).

same_text(A, B) :- text_to_string(A, S), text_to_string(B, S).

% --- HISTORY (append-only trend ledger) --------------------------------------
% One JSON object per line in benchmarks/fill_history.jsonl. search_inf and
% load_inf are the portable, comparable-over-time numbers; wall/rss ride along
% for same-host trend but stay informational. Each entry is stamped with the git
% commit + local timestamp.

history_path(BenchDir, Path) :-
    directory_file_path(BenchDir, 'fill_history.jsonl', Path).

append_history(BenchDir, Doc, Extra) :-
    history_path(BenchDir, Path),
    get_dict(results, Doc, Results),
    findall(F-Cell,
            ( member(Row, Results),
              get_dict(rung, Row, F0), atom_string(F, F0),
              Cell = _{ inf:      Row.search_inf_med,
                        load_inf: Row.load_inf,
                        wall_ms:  Row.cmd_wall_med_ms,
                         rss_kib:  Row.cmd_rss_med_kib } ),
            RungPairs),
    append_history(Path, Doc, Extra, RungPairs).

show_history(BenchDir) :-
    history_path(BenchDir, Path),
    read_history(Path, Entries),
    ( Entries == []
    ->  format("no fill benchmark history yet.~n"),
        format("record one with:  make bench-fill-record   (or  make bench-fill-log)~n")
    ;   render_history('fill benchmark',
                       [metric(search_inf, inf), metric(load_inf, load_inf)],
                       Entries) ).

:- multifile prolog:error_message//1.
prolog:error_message(fill_baseline_missing(Path)) -->
    [ 'check_fill_baseline: ~w not found (regenerate with: make bench-fill-record)'-[Path] ].
prolog:error_message(fill_bench_run_failed(Status)) -->
    [ 'check_fill_baseline: the fill product bench (run_fill.pl) failed: ~q'-[Status] ].
prolog:error_message(fill_record_readback_missing(Rung)) -->
    [ 'check_fill_baseline: recorded rung ~w is absent after baseline read-back'-[Rung] ].
prolog:error_message(fill_record_readback_mismatch(Rung, Key, Expected)) -->
    [ 'check_fill_baseline: recorded rung ~w field ~w does not equal measured value ~w after read-back'-[Rung, Key, Expected] ].
