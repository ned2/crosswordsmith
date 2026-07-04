#!/usr/bin/env swipl
% benchmarks/check_baseline.pl - performance ratchet for the arrange 15x15 ladder.
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
%   swipl -q benchmarks/check_baseline.pl --record   % RECORD: ratchet baseline.json to
%                                                    %   the measured numbers (accept wins),
%                                                    %   and append the run to history.jsonl
%   swipl -q benchmarks/check_baseline.pl --log      % LOG: diff + append to history.jsonl,
%                                                    %   WITHOUT moving the baseline
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

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(process)).
:- use_module(library(readutil)).
:- use_module(library(http/json)).

:- dynamic bench_dir/1.
:- prolog_load_context(directory, D), asserta(bench_dir(D)).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv0),
    ( select('--history', Argv0, _)      -> Mode = history, Extra = []
    ; select('--log', Argv0, Extra)      -> Mode = log
    ; select('--record', Argv0, Extra)   -> Mode = record
    ;                                       Mode = check, Extra = Argv0 ),
    bench_dir(BenchDir),
    ( Mode == history -> show_history(BenchDir), halt(0) ; true ),
    directory_file_path(BenchDir, 'baseline.json', BaselinePath),
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
    setup_call_cleanup(open(Path, read, S),
                       json_read_dict(S, Baseline, [default_tag(json)]),
                       close(S)).

% Spawn run_arrange.pl, capture its JSON stdout, parse it. Extra args (e.g. --heavy)
% pass through. stderr flows to ours; we gate on the child exit code.
run_product_bench(BenchDir, Extra, Doc) :-
    directory_file_path(BenchDir, 'run_arrange.pl', RunArrange),
    append(['-q', RunArrange, '--', '--format', json], Extra, Args),
    process_create(path(swipl), Args,
                   [ stdout(pipe(Out)), stderr(std), process(PID) ]),
    read_string(Out, _, JsonText),
    close(Out),
    process_wait(PID, Status),
    ( Status == exit(0) -> true ; throw(error(bench_run_failed(Status), _)) ),
    setup_call_cleanup(open_string(JsonText, S),
                       json_read_dict(S, Doc, [default_tag(json)]),
                       close(S)).

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

check_row(WL, Tol, VMatch, Row, F0-W0, F1-W1) :-
    get_dict(fixture, Row, Fx),
    get_dict(search_inf_med, Row, Meas),
    ( find_baseline(WL, Fx, Spec)
    ->  get_dict(search_inf, Spec, Base),
        ( Base =:= 0 -> Delta = 0.0 ; Delta is (Meas - Base) / Base * 100.0 ),
        classify(Delta, Tol, VMatch, Kind),
        kind_counts(Kind, DF, DW), F1 is F0 + DF, W1 is W0 + DW,
        signed(Delta, DStr), kind_label(Kind, Label),
        format("~w~t~30|~t~D~14+~t~D~14+~t~w%~11+   ~w~n",
               [Fx, Base, Meas, DStr, Label])
    ;   F1 = F0, W1 = W0,
        format("~w~t~30|~t~w~14+~t~D~14+~t~w~11+   ~w~n",
               [Fx, '(none)', Meas, '-', 'NEW (not in baseline)']) ).

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

kind_label(win,             'WIN (improvement)').
kind_label(ok,              'ok').
kind_label(regression,      'REGRESSION').
kind_label(regression_warn, 'regress? (swi-ver)').

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

do_record(BaselinePath, Baseline, Doc) :-
    get_dict(results, Doc, Results),
    get_dict(swi_prolog, Doc, RunSwi),
    current_host(RunHost),
    get_dict(workloads, Baseline, WL0),
    dict_pairs(WL0, Tag, Pairs0),
    maplist(record_pair(Results), Pairs0, Pairs1),
    new_rung_pairs(WL0, Results, NewPairs),
    append(Pairs1, NewPairs, Pairs2),
    dict_pairs(WL1, Tag, Pairs2),
    B1 = Baseline.put(_{host: RunHost, swi_prolog: RunSwi, workloads: WL1}),
    setup_call_cleanup(open(BaselinePath, write, S),
                       json_write_dict(S, B1, [width(90)]),
                       close(S)),
    format("baseline updated: ~w~n~n", [BaselinePath]),
    forall(member(Row, Results), report_recorded(WL0, Row)),
    unmeasured_note(WL0, Results).

record_pair(Results, K-V0, K-V1) :-
    ( result_for(Results, K, Row)
    ->  V1 = V0.put(_{ search_inf:      Row.search_inf_med,
                       cmd_wall_med_ms:  Row.cmd_wall_med_ms,
                       cmd_rss_med_kib:  Row.cmd_rss_med_kib })
    ;   V1 = V0 ).

% A measured rung with no baseline entry joins the baseline as a COMPLETE spec
% (the first --record after a rung is added to workloads.pl). The spec metadata
% (tier/warmup/budget/words) comes from the run_arrange result row.
new_rung_pairs(WL0, Results, Pairs) :-
    findall(K-Spec,
            ( member(Row, Results),
              get_dict(fixture, Row, F),
              \+ find_baseline(WL0, F, _),
              text_to_string(F, FS), atom_string(K, FS),
              new_rung_spec(Row, Spec) ),
            Pairs).

new_rung_spec(Row, _{ search_inf:      Row.search_inf_med,
                      cmd_wall_med_ms: Row.cmd_wall_med_ms,
                      cmd_rss_med_kib: Row.cmd_rss_med_kib,
                      info_only:       ["cmd_wall_med_ms", "cmd_rss_med_kib"],
                      grid:            Row.size,
                      words:           Row.words,
                      tier:            Row.tier,
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

% --- HISTORY (append-only trend ledger) --------------------------------------
% One JSON object per line in benchmarks/history.jsonl. search_inf is the portable,
% comparable-over-time number; wall/rss ride along for same-host trend but stay
% informational. Each entry is stamped with the git commit + local timestamp so a
% row can be attributed to the change that produced it.

history_path(BenchDir, Path) :-
    directory_file_path(BenchDir, 'history.jsonl', Path).

append_history(BenchDir, Doc, Extra) :-
    history_path(BenchDir, Path),
    get_dict(results, Doc, Results),
    get_dict(swi_prolog, Doc, Swi),
    current_host(Host),
    git_commit(Commit),
    git_dirty(Dirty),
    get_time(EpochF), EpochI is round(EpochF),
    format_time(atom(Date), '%Y-%m-%dT%H:%M:%S', EpochF),
    ( memberchk('--heavy', Extra) -> Tiers = 'core+heavy' ; Tiers = core ),
    findall(F-Cell,
            ( member(Row, Results),
              get_dict(fixture, Row, F0), atom_string(F, F0),
              Cell = _{ inf:     Row.search_inf_med,
                        wall_ms: Row.cmd_wall_med_ms,
                        rss_kib: Row.cmd_rss_med_kib } ),
            RungPairs),
    dict_pairs(Rungs, rungs, RungPairs),
    Entry = _{ epoch:EpochI, date:Date, commit:Commit, dirty:Dirty,
               host:Host, swi:Swi, tiers:Tiers, rungs:Rungs },
    with_output_to(string(Raw), json_write_dict(current_output, Entry, [width(0)])),
    normalize_space(string(Line), Raw),
    setup_call_cleanup(open(Path, append, S),
                       ( write(S, Line), nl(S) ),
                       close(S)),
    ( Dirty == true -> DMark = ' +dirty' ; DMark = '' ),
    format("~nlogged to ~w  (commit ~w~w, ~w)~n", [Path, Commit, DMark, Tiers]).

show_history(BenchDir) :-
    history_path(BenchDir, Path),
    ( exists_file(Path)
    ->  read_history(Path, Entries)
    ;   Entries = [] ),
    ( Entries == []
    ->  format("no benchmark history yet.~n"),
        format("record one with:  make bench-record   (or  make bench-log)~n")
    ;   render_history(Entries) ).

read_history(Path, Entries) :-
    read_file_to_string(Path, Str, []),
    split_string(Str, "\n", "", Lines0),
    exclude([L]>>(normalize_space(string(""), L)), Lines0, Lines),
    findall(E, ( member(L, Lines), parse_history_line(L, E) ), Entries).

parse_history_line(Line, Entry) :-
    catch(setup_call_cleanup(open_string(Line, S),
                             json_read_dict(S, Entry, [default_tag(json)]),
                             close(S)),
          _, fail).

render_history(Entries) :-
    length(Entries, N),
    format("benchmark history  (~d run(s), most recent last)~n~n", [N]),
    forall(nth1(I, Entries, E), print_entry_meta(I, E)),
    all_rung_keys(Entries, Keys),
    format("~nsearch_inf per rung  (vs prev = last step, vs first = cumulative):~n~n"),
    format("~w~t~30|~t~w~16+~t~w~12+~t~w~12+~t~w~7+~n",
           ['rung', 'latest', 'vs prev', 'vs first', 'runs']),
    forall(member(K, Keys), print_rung_trend(K, Entries)).

print_entry_meta(I, E) :-
    get_dict(date, E, Date), get_dict(commit, E, Commit),
    get_dict(tiers, E, Tiers), get_dict(host, E, Host),
    ( get_dict(dirty, E, true) -> DMark = '+dirty' ; DMark = '' ),
    format("  [~d] ~w  ~w~w  ~w  ~w~n", [I, Date, Commit, DMark, Tiers, Host]).

all_rung_keys(Entries, Keys) :-
    findall(K,
            ( member(E, Entries), get_dict(rungs, E, R),
              dict_pairs(R, _, Ps), member(K-_, Ps) ),
            Ks),
    sort(Ks, Keys).

% Trajectory for one rung across the runs that measured it: latest absolute value,
% step delta (vs the previous such run) and cumulative delta (vs the first).
print_rung_trend(K, Entries) :-
    rung_series(K, Entries, Series),
    ( Series == [] -> true
    ;   last(Series, Latest),
        Series = [First|_],
        ( append(_, [Prev, Latest], Series) -> true ; Prev = Latest ),
        length(Series, Cnt),
        pct_str(Latest, Prev, DPrev),
        pct_str(Latest, First, DFirst),
        format("~w~t~30|~t~D~16+~t~w~12+~t~w~12+~t~d~7+~n",
               [K, Latest, DPrev, DFirst, Cnt]) ).

rung_series(K, Entries, Series) :-
    findall(V,
            ( member(E, Entries), get_dict(rungs, E, R),
              get_dict(K, R, Cell), get_dict(inf, Cell, V) ),
            Series).

pct_str(_New, Old, 'n/a')  :- Old =:= 0, !.
pct_str(New, Old, Str) :- P is (New - Old) / Old * 100.0, signed(P, S), atom_concat(S, '%', Str).

% git provenance for a history entry. Both degrade gracefully outside a checkout.
git_commit(Commit) :-
    ( catch(run_capture(path(git), ['rev-parse', '--short', 'HEAD'], Raw), _, fail),
      normalize_space(atom(C), Raw), C \== ''
    -> Commit = C ; Commit = unknown ).

% dirty := any tracked modification (untracked files, e.g. scratch renders, are
% ignored so they don't perpetually flag the tree).
git_dirty(Dirty) :-
    ( catch(run_capture(path(git),
              ['status', '--porcelain', '--untracked-files=no'], Raw), _, fail),
      normalize_space(string(T), Raw), T \== ""
    -> Dirty = true ; Dirty = false ).

run_capture(Spec, Args, Out) :-
    process_create(Spec, Args, [stdout(pipe(S)), stderr(null), process(PID)]),
    read_string(S, _, Out), close(S), process_wait(PID, _).

:- multifile prolog:error_message//1.
prolog:error_message(baseline_missing(Path)) -->
    [ 'check_baseline: ~w not found (regenerate with: make bench-record)'-[Path] ].
prolog:error_message(bench_run_failed(Status)) -->
    [ 'check_baseline: the product bench (run_arrange.pl) failed: ~q'-[Status] ].
