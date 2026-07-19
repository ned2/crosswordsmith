% benchmarks/bench_core.pl - shared measurement engine for the crosswordsmith
% benchmarks (design: docs/benchmark-rework-plan.md §4).
%
% It measures an opaque SAMPLER: a closure that, invoked once, produces one
% sample as a dict of metric -> number (e.g. _{wall:0.03, cpu:0.03,
% inferences:319183}). It knows nothing about crosswords, strategies, or the CLI
% - the same loop serves the strategy matrix and the product bench. Adding a new
% thing to benchmark means writing a sampler, never editing this loop.

:- module(bench_core, [ measure/3, inproc_sampler/2, process_sampler/5 ]).

:- meta_predicate measure(1, +, -).
:- meta_predicate inproc_sampler(0, -).

:- use_module(library(apply), [exclude/3, foldl/4, maplist/2]).
:- use_module(library(lists), [append/3, member/2, nth0/3, sum_list/2]).
:- use_module(library(readutil), [read_file_to_string/3]).
:- use_module('bench_process.pl', [capture_process/6]).
% call_time/2 (inproc_sampler's measurement wrapper) is autoload-only
% (library(statistics)); explicit so the bench roots also run under
% autoload(false) (P11/C5).
:- use_module(library(statistics), [call_time/2]).

%!  measure(:Sampler, +Opts:dict, -Summary:dict) is det.
%
%   Sampler : called as call(Sampler, Sample); Sample is a dict metric->number.
%             Non-numeric keys (e.g. an outcome tag the caller attaches) are
%             carried in the sample but IGNORED by the stats here.
%   Opts    : _{warmup:W, iterations:N}  (W >= 0, N >= 1)
%   Summary : _{iterations:N, warmup:W, stats:_{Metric:_{min,median,mean}, ...}}
%
% Runs W warmup samples (discarded), then N measured samples; for each NUMERIC
% metric key present it computes {min, median, mean}. Metric sets may differ
% between samplers (command has wall+rss; search has wall+cpu+inferences), but
% every measured sample from one sampler must carry the same numeric keys.
%
% SUCCESS-ONLY contract (plan §4, stress-test m1): measure/3 assumes the caller
% has already established that the sampler will succeed - run_matrix gates each
% cell through solve_status/6 (solved|no|timeout) BEFORE calling here. A measured
% sample that fails is a caller-gating bug, not a data point, so it throws rather
% than record a bogus zero-work sample. Warmup tolerates failure (discarded
% regardless).
measure(Sampler, Opts, Summary) :-
    get_dict(warmup, Opts, Warmup),
    get_dict(iterations, Opts, Iters),
    ( integer(Warmup), Warmup >= 0 -> true ; throw(error(bench_bad_opt(warmup, Warmup), _)) ),
    ( integer(Iters),  Iters  >= 1 -> true ; throw(error(bench_bad_opt(iterations, Iters), _)) ),
    forall(between(1, Warmup, _), ignore(call(Sampler, _))),
    findall(Sample,
            ( between(1, Iters, _), measured_sample(Sampler, Sample) ),
            Samples),
    summarize_samples(Samples, Stats),
    Summary = _{iterations:Iters, warmup:Warmup, stats:Stats}.

measured_sample(Sampler, Sample) :-
    ( call(Sampler, Sample) -> true
    ; throw(error(bench_sampler_failed(Sampler), _)) ).

%!  inproc_sampler(:Goal, -Sample:dict) is det.
%
% In-process Prolog goal - the search layer. call_time/2 returns
% time{cpu,inferences,wall}; inferences is the SWI-version-locked metric of record.
% The Goal MUST be deterministic (the CALLER pins the first solution with its own
% cut - solve_once/4 cuts, arrange_best_layout/5 is single-valued). We do NOT wrap
% in once/1: that would add one inference per call, breaking byte-for-byte
% reproducibility of the recorded inference counts against pre-rework batches.
inproc_sampler(Goal, _{wall:W, cpu:C, inferences:I}) :-
    call_time(Goal, T),
    W = T.wall, C = T.cpu, I = T.inferences.

%!  process_sampler(+Exe, +Argv:list, -WallSeconds:number,
%!                  -MaxRssKiB:number, -ExitCode:integer) is det.
%
%   The COMMAND layer: run `Exe Argv...` as a fresh child under
%   `/usr/bin/time -f "%e %M"`, capturing end-to-end wall (%e, seconds) and peak
%   resident set (%M, KiB). This is the honest top-level latency a user feels -
%   SWI startup, load.pl, arg parsing, search, emit, all of it. The timing goes
%   to a scratch file via `-o` so it never mixes with the command's own streams
%   (both discarded; a successful arrange is quiet with --out /dev/null).
%
%   Deliberately NOT a measure/3 sampler on its own: it returns the raw exit code
%   so a caller (subjects.pl) can assert it against the workload's Expected before
%   wrapping wall+rss into a sample dict. rss is a whole-process footprint, not a
%   search-memory metric (plan §7).
process_sampler(Exe, Argv, Wall, Rss, Exit) :-
    tmp_file_stream(text, TimeFile, TStream), close(TStream),
    setup_call_cleanup(
        true,
        run_under_time(Exe, Argv, TimeFile, Wall, Rss, Exit),
        delete_file(TimeFile)).

run_under_time(Exe, Argv, TimeFile, Wall, Rss, Exit) :-
    append(['-o', TimeFile, '-f', '%e %M', Exe], Argv, TimeArgv),
    capture_process('/usr/bin/time', TimeArgv, null, _Stdout, _Stderr, Status),
    ( Status = exit(Exit) -> true ; Exit = -1 ),
    read_time_file(TimeFile, Wall, Rss).

% Parse the `%e %M` line GNU time wrote. Take the FIRST line carrying two numeric
% tokens: a failed command makes time append a "Command exited..." line, which we
% skip rather than misparse.
read_time_file(TimeFile, Wall, Rss) :-
    read_file_to_string(TimeFile, Str, []),
    split_string(Str, "\n", "", Lines),
    ( member(Line, Lines),
      split_string(Line, " ", " ", Toks0),
      exclude(==(""), Toks0, [WS, RS | _]),
      catch((number_string(Wall, WS), number_string(Rss, RS)), _, fail)
    ->  true
    ;   throw(error(bench_time_parse_failed(Str), _)) ).

% Per-metric {min,median,mean} over the samples, for every NUMERIC-valued key.
% Non-numeric annotations are skipped, but numeric keys must be stable across
% samples so a missing metric cannot silently reduce its sample count.
summarize_samples(Samples, Stats) :-
    Samples = [First|_],
    numeric_keys(First, NumKeys),
    maplist(require_numeric_keys(NumKeys), Samples),
    foldl(summarize_key(Samples), NumKeys, _{}, Stats).

numeric_keys(Sample, Keys) :-
    findall(Key, ( get_dict(Key, Sample, Value), number(Value) ), Keys0),
    sort(Keys0, Keys).

require_numeric_keys(Expected, Sample) :-
    numeric_keys(Sample, Actual),
    ( Actual == Expected -> true
    ; throw(error(bench_sample_schema_mismatch(Expected, Actual), _)) ).

summarize_key(Samples, Key, SIn, SOut) :-
    findall(V, ( member(S, Samples), get_dict(Key, S, V) ), Vs),
    metric_stat(Vs, Stat),
    put_dict(Key, SIn, Stat, SOut).

% min / mean / median over a non-empty value list. Median averages the two
% central values for even N - the SINGLE canonical definition, unifying the old
% run_benchmarks (averaged) / run_matrix (upper-of-two) disagreement
% (docs/prolog-audit-findings.md F025). Immaterial for the deterministic
% inferences metric (all samples equal => min == median == mean); a sub-noise
% shift for wall.
metric_stat(Vs, _{min:Min, median:Median, mean:Mean}) :-
    Vs = [_|_],
    msort(Vs, Sorted),
    Sorted = [Min|_],
    length(Sorted, Len),
    sum_list(Sorted, Sum),
    Mean is Sum / Len,
    median_of(Sorted, Len, Median).

median_of(Sorted, Len, Median) :-
    ( 1 is Len mod 2
    ->  Mid is Len // 2, nth0(Mid, Sorted, Median)
    ;   Hi is Len // 2, Lo is Hi - 1,
        nth0(Lo, Sorted, A), nth0(Hi, Sorted, B),
        Median is (A + B) / 2 ).

:- multifile prolog:error_message//1.
prolog:error_message(bench_sampler_failed(_)) -->
    [ 'bench_core: a measured sampler failed; the caller must gate for success (solve_status) before measure/3' ].
prolog:error_message(bench_bad_opt(Key, Val)) -->
    [ 'bench_core: measure/3 option ~w must be a valid count, got ~q'-[Key, Val] ].
prolog:error_message(bench_sample_schema_mismatch(Expected, Actual)) -->
    [ 'bench_core: measured sample numeric keys changed; expected ~q, got ~q'-
      [Expected, Actual] ].
