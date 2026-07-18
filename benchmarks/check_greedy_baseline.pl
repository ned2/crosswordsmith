#!/usr/bin/env swipl
% benchmarks/check_greedy_baseline.pl - independent greedy inference ratchet.
% Construction, sweep, and postprocess inferences gate at 0.5%; sweep is the
% primary metric. Command wall/RSS are informational.
% --exact requires same-SWI equality over every core and heavy gated metric.

:- module(check_greedy_baseline, []).

:- set_prolog_flag(verbose, silent).
:- use_module(library(apply), [foldl/4]).
:- use_module(library(filesex), [directory_file_path/3]).
:- use_module(library(json)).
:- use_module(library(lists)).
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
                render_history/3
              ]).

:- prolog_load_context(directory, D), assertz(bench_dir(D)).
:- dynamic bench_dir/1.
:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv0),
    ( Argv0 == ['--self-test'] -> self_test, halt(0) ; true ),
    catch(checker_mode(Argv0, Mode, Extra),
          E, (print_message(error, E), halt(2))),
    bench_dir(BenchDir),
    ( Mode == history -> show_history(BenchDir), halt(0) ; true ),
    directory_file_path(BenchDir, 'greedy_baseline.json', Path),
    load_json(Path, Baseline),
    ( Mode == exact
    -> catch(exact_runner_args(Extra, ExactArgs),
             E, (print_message(error, E), halt(2))),
       run_bench(BenchDir, ExactArgs, Doc),
       do_exact_check(Baseline, Doc, Fails),
       report_exact_result(Fails), halt_if(Fails)
    ;  run_bench(BenchDir, Extra, Doc),
       do_check(Baseline, Doc, Fails, Wins),
       ( Mode == check
       -> report(Fails,Wins), halt_if(Fails)
       ; Mode == log
       -> append_history(BenchDir,Doc,Extra), halt(0)
       ; Mode == record
       -> record_checked(Path,Baseline,Doc), append_history(BenchDir,Doc,Extra), halt(0)
       ; Mode == promote
       -> report(Fails,Wins),
          ( Fails =:= 0
          -> record_checked(Path,Baseline,Doc), append_history(BenchDir,Doc,Extra), halt(0)
          ;  format("promote: regressions present; baseline unchanged~n"), halt(1) )
       )
    ).

halt_if(0) :- halt(0).
halt_if(_) :- halt(1).

load_json(Path, Dict) :-
    read_json_dict(Path, Dict).

run_bench(BenchDir, Extra, Doc) :-
    directory_file_path(BenchDir, 'run_arrange_greedy.pl', Runner),
    append(['-q',Runner,'--','--format',json], Extra, Args),
    capture_process(path(swipl),Args,inherit,Text,_Stderr,Status),
    ( Status == exit(0) -> atom_json_dict(Text,Doc,[default_tag(json)])
    ; throw(error(greedy_bench_failed(Status),_)) ).

do_check(Baseline, Doc, Fails, Wins) :-
    get_dict(workloads,Baseline,WL), get_dict(results,Doc,Rows),
    get_dict(regression_tolerance_pct,Baseline,Tol),
    get_dict(swi_prolog,Baseline,BaseSwi), get_dict(swi_prolog,Doc,RunSwi),
    same_text(BaseSwi,RunSwi,VMatch),
    format("greedy inference ratchet (sweep primary; tolerance +~w%)~n",[Tol]),
    foldl(check_result(WL,Tol,VMatch),Rows,0-0,Fails-Wins),
    report_unmeasured(WL,Rows).

do_exact_check(Baseline, Doc, Fails) :-
    get_dict(swi_prolog, Baseline, BaseSwi),
    get_dict(swi_prolog, Doc, RunSwi),
    exact_version(BaseSwi, RunSwi, VersionStatus, VersionFails),
    format("greedy exact inference check (+heavy tail)~n"),
    format("exact SWI: ~w -> ~w (~w)~n", [BaseSwi, RunSwi, VersionStatus]),
    get_dict(workloads, Baseline, WL),
    get_dict(results, Doc, Rows),
    exact_row_presence(WL, Rows, Missing, Unexpected, PresenceFails),
    foldl(exact_greedy_row(WL), Rows, 0, MetricFails),
    report_exact_presence(Missing, Unexpected),
    Fails is VersionFails + PresenceFails + MetricFails.

exact_greedy_row(WL, Row, F0, F3) :-
    get_dict(rung, Row, Rung),
    format("~n~w~n", [Rung]),
    (   baseline_spec(WL, Rung, Spec)
    ->  exact_greedy_metric(construction_inf, Spec.construction_inf,
                            Row.metrics.construction_inf_med, F0, F1),
        exact_greedy_metric(sweep_inf, Spec.sweep_inf,
                            Row.metrics.sweep_inf_med, F1, F2),
        exact_greedy_metric(postprocess_inf, Spec.postprocess_inf,
                            Row.metrics.postprocess_inf_med, F2, F3)
    ;   F3 = F0,
        format("  gated metrics: no_reference~n")
    ).

exact_greedy_metric(Metric, Base, Measured, F0, F1) :-
    exact_metric(Base, Measured, Status, MetricFails),
    F1 is F0 + MetricFails,
    format("  ~w: ~D -> ~D ~w~n", [Metric, Base, Measured, Status]).

exact_row_presence(WL, Rows, Missing, Unexpected, Fails) :-
    dict_pairs(WL, _, Pairs),
    pairs_keys(Pairs, ReferenceIds),
    findall(Rung, (member(Row, Rows), get_dict(rung, Row, Rung)), MeasuredIds),
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

check_result(WL,Tol,VMatch,Row,F0-W0,F-W) :-
    get_dict(rung,Row,Rung),
    format("~n~w~n",[Rung]),
    check_metric(WL,Rung,construction_inf,Row.metrics.construction_inf_med,
                 Tol,VMatch,F0-W0,F1-W1),
    check_metric(WL,Rung,sweep_inf,Row.metrics.sweep_inf_med,
                 Tol,VMatch,F1-W1,F2-W2),
    check_metric(WL,Rung,postprocess_inf,Row.metrics.postprocess_inf_med,
                 Tol,VMatch,F2-W2,F-W),
    format("  info command wall ~1f ms, rss ~0f KiB~n",
           [Row.metrics.command_wall_med_ms,Row.metrics.command_rss_med_kib]).

check_metric(WL,Rung,Key,Measured,Tol,VMatch,F0-W0,F-W) :-
    ( baseline_spec(WL,Rung,Spec), get_dict(Key,Spec,Base)
    -> delta(Base,Measured,Pct), classify(Pct,Tol,VMatch,Kind),
       counts(Kind,DF,DW), F is F0+DF, W is W0+DW,
       format("  ~w: ~D -> ~D (~2f%) ~w~n",[Key,Base,Measured,Pct,Kind])
    ;  F=F0, W=W0,
       format("  ~w: NEW at ~D~n",[Key,Measured])
    ).

delta(0,_M,0.0) :- !.
delta(Base,M,Pct) :- Pct is (M-Base)/Base*100.0.

classify(P,T,_V,win) :- P =< -T, !.
classify(P,T,_V,ok) :- P =< T, !.
classify(_P,_T,true,regression) :- !.
classify(_P,_T,false,regression_warn).
counts(win,0,1). counts(ok,0,0). counts(regression,1,0). counts(regression_warn,0,0).

report(F,W) :-
    ( F > 0 -> format("~nRESULT: FAIL (~d regressions, ~d wins)~n",[F,W])
    ; W > 0 -> format("~nRESULT: PASS (~d wins, 0 regressions)~n",[W])
    ; format("~nRESULT: PASS (all deterministic metrics +0.00%)~n") ).

report_unmeasured(WL,Rows) :-
    dict_pairs(WL,_,Pairs),
    findall(K,(member(K-_,Pairs),\+ row_for(Rows,K,_)),Missing),
    ( Missing=[] -> true ; format("not measured this pass: ~w~n",[Missing]) ).

baseline_spec(WL,Rung,Spec) :-
    dict_pairs(WL,_,Pairs), member(K-Spec,Pairs), same_text(K,Rung,true), !.
row_for(Rows,Rung,Row) :-
    member(Row,Rows), same_text(Row.rung,Rung,true), !.

same_text(A,B,Same) :-
    text_to_string(A,S1), text_to_string(B,S2),
    ( S1==S2 -> Same=true ; Same=false ).

% Atomic one-measurement recording: write and verify a sibling temporary file,
% then rename it over the target. A failed write/read-back leaves the committed
% baseline intact. Unmeasured heavy rows are retained.
record_checked(Path, Baseline, Doc) :-
    build_recorded(Baseline,Doc,Recorded),
    replace_json_dict(Path, Recorded, 100, verify_readback(Doc)),
    format("baseline recorded and read-back verified: ~w~n",[Path]).

build_recorded(Baseline,Doc,Recorded) :-
    build_recorded_baseline(Baseline, Doc, record_key, record_spec, BaseRecorded),
    Recorded = BaseRecorded.put(generated_note,
        'construction/sweep/postprocess inferences gate independently at 0.5%; sweep is primary. Wall/RSS are informational. Regenerate with make bench-greedy-record BENCH_ARGS=--heavy.').

record_key(Row, Row.rung).

record_spec(Row, _Prior, Spec) :-
    Metrics = Row.metrics,
    Spec = _{construction_inf:Metrics.construction_inf_med,
             sweep_inf:Metrics.sweep_inf_med,
             postprocess_inf:Metrics.postprocess_inf_med,
             command_wall_med_ms:Metrics.command_wall_med_ms,
             command_rss_med_kib:Metrics.command_rss_med_kib,
             tier:Row.tier, fixture:Row.fixture, grid:Row.size,
             framing:Row.framing, mode:Row.mode,
             construction:Row.construction, words:Row.words,
             iterations:Row.iterations, warmup:Row.warmup,
             info_only:["command_wall_med_ms", "command_rss_med_kib"]}.

text_to_atom(Text,Atom) :- text_to_string(Text,S), atom_string(Atom,S).

verify_readback(Doc, ReadBack) :-
    get_dict(workloads,ReadBack,WL), get_dict(results,Doc,Rows),
    forall(member(Row,Rows),
           ( baseline_spec(WL,Row.rung,Spec),
             Spec.construction_inf =:= Row.metrics.construction_inf_med,
             Spec.sweep_inf =:= Row.metrics.sweep_inf_med,
             Spec.postprocess_inf =:= Row.metrics.postprocess_inf_med )),
    length(Rows,Expected),
    findall(R,(member(R,Rows),baseline_spec(WL,R.rung,_)),Found),
    length(Found,Expected).

history_path(BenchDir,Path) :- directory_file_path(BenchDir,'greedy_history.jsonl',Path).

append_history(BenchDir,Doc,Extra) :-
    history_path(BenchDir,Path), get_dict(results,Doc,Rows),
    findall(K-Cell,
            ( member(R,Rows), text_to_atom(R.rung,K), M=R.metrics,
              Cell=_{construction_inf:M.construction_inf_med,
                     sweep_inf:M.sweep_inf_med,
                     postprocess_inf:M.postprocess_inf_med,
                     wall_ms:M.command_wall_med_ms,
                      rss_kib:M.command_rss_med_kib} ),
            Pairs),
    append_history(Path, Doc, Extra, Pairs).

show_history(BenchDir) :-
    history_path(BenchDir,Path),
    read_history(Path, Entries),
    ( Entries == []
    -> format("no greedy benchmark history yet.~n")
    ; render_history('greedy benchmark',
                     [metric(construction_inf, construction_inf),
                      metric(sweep_inf, sweep_inf),
                      metric(postprocess_inf, postprocess_inf)],
                     Entries)
    ).

% Fast permanent self-test for promotion safety. It proves core-only recording
% retains a heavy row, updates an existing row, inserts a new row, and verifies
% all three measured metrics from the exact one Doc used to write.
self_test :-
    tmp_file_stream(text,Path,S), close(S),
    Initial=_{generated_note:"test",host:"h",metric:"sweep_inf_med",mode:"ratchet",
              regression_tolerance_pct:0.5,swi_prolog:"10.1.10",tool:"t",
              workloads:_{core:_{construction_inf:1,sweep_inf:2,postprocess_inf:3},
                           heavy:_{construction_inf:7,sweep_inf:8,postprocess_inf:9}}},
    replace_json_dict(Path, Initial, 100, accept_readback),
    R1=_{rung:"core",fixture:"f",size:1,framing:"size",mode:_{kind:"best_effort"},
         construction:_{},words:1,iterations:1,warmup:0,tier:"core",
         metrics:_{construction_inf_med:10,sweep_inf_med:20,postprocess_inf_med:30,
                   command_wall_med_ms:1.0,command_rss_med_kib:2}},
    R2=R1.put(_{rung:"new"}), Doc=_{swi_prolog:"10.1.10",results:[R1,R2]},
    record_checked(Path,Initial,Doc), load_json(Path,Back),
    Back.workloads.heavy.sweep_inf =:= 8,
    Back.workloads.core.sweep_inf =:= 20,
    Back.workloads.new.postprocess_inf =:= 30,
    delete_file(Path),
    format("greedy baseline self-test: PASS~n").

accept_readback(_).

:- multifile prolog:error_message//1.
prolog:error_message(greedy_bench_failed(Status)) -->
    ['greedy benchmark runner failed: ~q'-[Status]].
