% benchmarks/subjects.pl - the crosswordsmith-specific measurement subjects.
%
% This is the ONLY benchmark file that knows what an "arrange" is. It adapts the
% generic samplers in bench_core into measure/3-ready closures for the two layers
% the product bench reports (plan §5):
%
%   - COMMAND layer: the whole `crosswordsmith arrange` process, timed end-to-end
%     (SWI startup + load + parse + search + emit) - the latency a user feels.
%   - SEARCH layer: the in-process 4-corner arrange search alone (arrange_best_
%     layout/5), so the CLI wrapper's fixed cost can be subtracted out (rest = A-B).
%
% Both layers read the SAME .pl fixture (arrange --input accepts .pl), so they
% measure the same word set. Each layer asserts its result against the workload's
% Expected before recording, so a workload that quietly stopped placing can never
% masquerade as a fast run.

:- module(bench_subjects,
          [ size_flag/4,
            arrange_command_sampler/6,
            arrange_search_sampler/5 ]).

:- use_module(library(lists)).

% --size N and --max-size N both frame the SAME N x N canvas (GridLen = N); only
% the emit framing differs (crop vs full). So the search layer takes GridLen = N
% for either mode, and the command layer just picks the matching flag.
size_flag(size,     N, '--size',     A) :- atom_number(A, N).
size_flag(max_size, N, '--max-size', A) :- atom_number(A, N).

% COMMAND layer sampler. process_sampler returns the raw exit code; we assert it
% against Expected here (not inside bench_core) and record only wall+rss.
% rss (peak RSS, KiB) is a whole-process footprint, not a search-memory metric.
arrange_command_sampler(Exe, File, Size, Mode, Expected, _{wall:Wall, rss:Rss}) :-
    size_flag(Mode, Size, Flag, Val),
    bench_core:process_sampler(Exe,
        ['arrange', '--input', File, Flag, Val, '--out', '/dev/null'],
        Wall, Rss, Exit),
    expected_exit(Expected, Exit).

expected_exit(placed, 0) :- !.
expected_exit(infeasible, 1) :- !.
expected_exit(Expected, Got) :-
    throw(error(bench_command_exit(expected(Expected), exit(Got)), _)).

% SEARCH layer sampler. Uses the budget-EXPLICIT arrange_best_layout/6 (arrange.pl)
% so the bench can raise the inference budget above the shipped 500M default: a
% pathological workload that would otherwise SATURATE the budget (its count pinned
% to the budget constant, useless as a hill-climbing signal) instead runs to true
% completion, yielding a deterministic count that reflects the search's real cost.
%
% WHITE-BOX REACH, deliberate (C25): arrange_best_layout/6 is NOT on
% crosswordsmith_arrange's export list - the module-qualified call below is
% this harness reaching an internal on purpose. It must stay on exactly this
% predicate: every recorded baseline.json/history.jsonl search_inf count is
% defined against it, so re-pointing at a public seam (e.g. arrange_outcome/5)
% would silently move the gated counts. If arrange.pl ever exports /6, switch
% to a plain import; until then this is the sanctioned exception.
% Any Budget >= a completing fixture's true cost gives the SAME count, so this does
% not perturb the fast workloads. arrange_best_layout/6 is single-valued and always
% succeeds with an Outcome (placed | not_proven | infeasible); we assert it matches
% Expected per iteration (plan m2) then DROP it - the sample carries only the
% numeric wall/cpu/inferences that measure/3 summarizes.
arrange_search_sampler(Words, GridLen, Budget, Expected, Sample) :-
    bench_core:inproc_sampler(
        crosswordsmith_arrange:arrange_best_layout(Words, GridLen, Budget, _, _, Outcome),
        Sample),
    ( outcome_ok(Expected, Outcome) -> true
    ; throw(error(bench_search_outcome(expected(Expected), got(Outcome)), _)) ).

outcome_ok(placed, placed).
outcome_ok(infeasible, infeasible).
outcome_ok(infeasible, not_proven).   % budget-limited non-placement still "didn't place"

:- multifile prolog:error_message//1.
prolog:error_message(bench_command_exit(expected(E), exit(G))) -->
    [ 'arrange command exit ~w does not match expected workload outcome ~w'-[G, E] ].
prolog:error_message(bench_search_outcome(expected(E), got(G))) -->
    [ 'arrange search outcome ~w does not match expected workload outcome ~w'-[G, E] ].
