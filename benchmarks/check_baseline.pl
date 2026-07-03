#!/usr/bin/env swipl
% benchmarks/check_baseline.pl - regression gate for the PRODUCT arrange bench.
%
% Runs the core product bench (benchmarks/run_arrange.pl --format json, core tier
% only) and diffs the result against the committed benchmarks/baseline.json.
%
% The GATED metric is `search_inf` on the warm core cells: arrange search
% inferences are deterministic and machine-independent, so any change is a real
% search regression (or an intended algorithm change that needs the baseline
% regenerated). A mismatch on a gate:"exact" cell fails the check (halt(1)).
%
% Everything else is INFORMATIONAL and NEVER fails the check:
%   - wall/rss are machine-dependent (and `time %e` is only 10ms-resolution): the
%     deltas are printed, and large SAME-HOST drift is WARNed, but never gated.
%   - Inference counts are portable across MACHINES but not across SWI-Prolog
%     VERSIONS. If the running SWI differs from the baseline's, a search_inf
%     mismatch is downgraded to a WARN (the version, not a regression, explains
%     it) and the baseline should be regenerated on this version.
%
% Usage:
%   swipl -q benchmarks/check_baseline.pl     % run core bench + compare; exit 0/1
%   make bench-check

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(process)).
:- use_module(library(readutil)).
:- use_module(library(http/json)).

:- dynamic bench_dir/1.
:- prolog_load_context(directory, BenchDir), asserta(bench_dir(BenchDir)).

:- initialization(main, main).

main :-
    bench_dir(BenchDir),
    directory_file_path(BenchDir, 'baseline.json', BaselinePath),
    load_baseline(BaselinePath, Baseline),
    format("crosswordsmith arrange - regression check~n"),
    format("running core product bench (this takes ~~20s)...~n~n"),
    run_product_bench(BenchDir, Doc),
    check_all(Baseline, Doc, Fails, Warns),
    ( Fails =:= 0
    ->  format("~nRESULT: PASS  (all exact-gate cells match; ~d warning(s))~n", [Warns]),
        halt(0)
    ;   format("~nRESULT: FAIL  (~d exact-gate cell(s) drifted; ~d warning(s))~n", [Fails, Warns]),
        format("If this was an INTENTIONAL search change, regenerate benchmarks/baseline.json~n"),
        format("(swipl -q benchmarks/run_arrange.pl -- --format json) and review the diff.~n"),
        halt(1) ).

% --- load + run --------------------------------------------------------------

load_baseline(Path, Baseline) :-
    ( exists_file(Path) -> true ; throw(error(baseline_missing(Path), _)) ),
    setup_call_cleanup(open(Path, read, S),
                       json_read_dict(S, Baseline, [default_tag(json)]),
                       close(S)).

% Spawn run_arrange.pl as a child, capture its JSON stdout, parse it. run_arrange
% with no --heavy measures exactly the core tier (the gated cells). Its stderr is
% passed through so any load warning is visible; we gate on the child exit code.
run_product_bench(BenchDir, Doc) :-
    directory_file_path(BenchDir, 'run_arrange.pl', RunArrange),
    process_create(path(swipl), ['-q', RunArrange, '--', '--format', 'json'],
                   [ stdout(pipe(Out)), stderr(std), process(PID) ]),
    read_string(Out, _, JsonText),
    close(Out),
    process_wait(PID, Status),
    ( Status == exit(0) -> true ; throw(error(bench_run_failed(Status), _)) ),
    setup_call_cleanup(open_string(JsonText, S),
                       json_read_dict(S, Doc, [default_tag(json)]),
                       close(S)).

% --- compare -----------------------------------------------------------------

check_all(Baseline, Doc, Fails, Warns) :-
    get_dict(results, Doc, Results),
    get_dict(swi_prolog, Doc, RunSwi),
    get_dict(swi_prolog, Baseline, BaseSwi),
    get_dict(host, Baseline, BaseHost),
    current_host(RunHost),
    ( same_text(RunSwi, BaseSwi) -> VersionMatch = true ; VersionMatch = false ),
    ( same_text(RunHost, BaseHost) -> HostMatch = true ; HostMatch = false ),
    print_env(BaseHost, RunHost, HostMatch, BaseSwi, RunSwi, VersionMatch),
    get_dict(workloads, Baseline, Workloads),
    dict_pairs(Workloads, _, Pairs),
    % gate table
    format("~w~t~34|~t~w~16+~t~w~16+  ~w~n",
           ['workload', 'baseline_inf', 'measured_inf', 'gate']),
    foldl(check_gate(Results, VersionMatch), Pairs, 0-0, Fails-GateWarns),
    % info table (wall/rss) - never gates
    format("~ninfo (host-specific, reporting-only):~n"),
    format("~w~t~34|~w~t~26+~w~n", ['workload', 'cmd_wall_med_ms', 'cmd_rss_med_kib']),
    foldl(check_info(Results, HostMatch), Pairs, 0, InfoWarns),
    Warns is GateWarns + InfoWarns.

check_gate(Results, VersionMatch, Fixture-Spec, F0-W0, F1-W1) :-
    ( get_dict(gate, Spec, G0), same_text(G0, exact) -> Gate = exact ; Gate = other ),
    get_dict(search_inf, Spec, Base),
    ( find_row(Results, Fixture, Row)
    ->  get_dict(search_inf_med, Row, Meas),
        classify_gate(Gate, VersionMatch, Base, Meas, Status, DF, DW),
        F1 is F0 + DF, W1 is W0 + DW
    ;   Meas = missing, Status = 'MISSING!', F1 is F0 + 1, W1 = W0 ),
    print_gate_row(Fixture, Base, Meas, Status).

% exact gate: equal -> PASS; unequal -> FAIL if the SWI version matches the
% baseline (a real regression), else a WARN (version, not code, moved the count).
classify_gate(exact, _, Base, Meas, 'PASS', 0, 0) :- Base =:= Meas, !.
classify_gate(exact, true,  _Base, _Meas, 'FAIL', 1, 0) :- !.
classify_gate(exact, false, _Base, _Meas, 'WARN(swi-ver)', 0, 1) :- !.
classify_gate(_Other, _, _Base, _Meas, 'info', 0, 0).

print_gate_row(Fixture, Base, missing, Status) :- !,
    format("~w~t~34|~t~d~16+~t~w~16+  ~w~n", [Fixture, Base, '-', Status]).
print_gate_row(Fixture, Base, Meas, Status) :-
    format("~w~t~34|~t~d~16+~t~d~16+  ~w~n", [Fixture, Base, Meas, Status]).

check_info(Results, HostMatch, Fixture-Spec, W0, W1) :-
    ( find_row(Results, Fixture, Row)
    ->  metric_delta(Spec, Row, cmd_wall_med_ms, WallStr, WWarn),
        metric_delta(Spec, Row, cmd_rss_med_kib, RssStr, RWarn),
        ( HostMatch == true -> HW is WWarn + RWarn ; HW = 0 ),
        W1 is W0 + HW,
        format("~w~t~34|~w~t~26+~w~n", [Fixture, WallStr, RssStr])
    ;   W1 = W0,
        format("~w~t~34|~w~n", [Fixture, '(not measured)']) ).

% Format "base -> meas (+d.d%)"; flag a same-host warn when |drift| exceeds 25%.
metric_delta(Spec, Row, Key, Str, Warn) :-
    get_dict(Key, Spec, Base),
    get_dict(Key, Row, Meas),
    ( Base =:= 0 -> Pct = 0.0 ; Pct is 100.0 * (Meas - Base) / Base ),
    ( abs(Pct) > 25.0 -> Warn = 1, Flag = ' !' ; Warn = 0, Flag = '' ),
    format(atom(Str), "~w -> ~w (~1f%)~w", [Base, Meas, Pct, Flag]).

% --- env banner --------------------------------------------------------------

print_env(BaseHost, RunHost, HostMatch, BaseSwi, RunSwi, VersionMatch) :-
    ( HostMatch == true
    ->  HostNote = 'same host -> wall/rss deltas meaningful'
    ;   HostNote = 'DIFFERENT host -> wall/rss NOT comparable (search_inf still portable)' ),
    ( VersionMatch == true
    ->  VerNote = 'same swi  -> search_inf gate ACTIVE'
    ;   VerNote = 'DIFFERENT swi -> search_inf mismatches downgraded to WARN' ),
    format("baseline:  host ~w,  swi ~w~n", [BaseHost, BaseSwi]),
    format("this run:  host ~w,  swi ~w~n", [RunHost, RunSwi]),
    format("  ~w~n", [HostNote]),
    format("  ~w~n~n", [VerNote]).

% --- helpers -----------------------------------------------------------------

% Baseline keys are atoms; run_arrange JSON strings parse to Prolog strings.
find_row(Results, Fixture, Row) :-
    member(Row, Results),
    get_dict(fixture, Row, F),
    same_text(F, Fixture), !.

same_text(A, B) :- text_to_string(A, S), text_to_string(B, S).

% Host only tags wall/rss comparability (never the gate), so a missing/odd uname
% degrades to 'unknown' rather than aborting the whole check.
current_host(Host) :-
    catch(uname_nm(Raw), _, fail),
    normalize_space(atom(Host), Raw), Host \== '', !.
current_host(unknown).

uname_nm(Raw) :-
    process_create(path(uname), ['-nm'],
                   [ stdout(pipe(Out)), stderr(null), process(PID) ]),
    read_string(Out, _, Raw), close(Out),
    process_wait(PID, _).

:- multifile prolog:error_message//1.
prolog:error_message(baseline_missing(Path)) -->
    [ 'check_baseline: ~w not found (regenerate from run_arrange --format json)'-[Path] ].
prolog:error_message(bench_run_failed(Status)) -->
    [ 'check_baseline: the product bench (run_arrange.pl) failed: ~q'-[Status] ].
