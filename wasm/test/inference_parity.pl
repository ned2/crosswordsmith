% inference_parity.pl - gate #2: native <-> WASM arrange inference-count parity.
%
% This is the empirical native/WASM parity gate for the pinned runtimes. It
% measures each ladder rung's search-layer inference count and prints it. Run it
% under both native swipl and wasm/node swipl and diff the outputs; equality
% confirms that the current same-SWI ratchet is also a valid proxy for this WASM
% build. Re-run the gate after changing either runtime rather than assuming
% inference counts are portable.
%
% It measures the IDENTICAL quantity the product bench records: call_time/2 around
% arrange_best_layout/6 (the budget-explicit search), with NO once/1 wrap - exactly
% benchmarks/bench_core.pl:inproc_sampler/2, so a native run reproduces
% benchmarks/baseline.json byte-for-byte. It deliberately does NOT load bench_core /
% subjects, because those pull library(process) (the command layer) which is
% unavailable under wasm.
%
% Run:
%   swipl -q wasm/test/inference_parity.pl -- [--heavy]
%   node ~/src/swipl-devel/build.wasm/src/swipl.js -q wasm/test/inference_parity.pl -- [--heavy]
% Output (stdout): CSV rows `fixture,grid,words,warmup,inferences,outcome`, one per
% rung, in manifest order. Project-load warnings (the cosmetic http/json noise
% under wasm) go to stderr; redirect 2>/dev/null for a clean diff.

:- set_prolog_flag(verbose, silent).
:- dynamic repo_root/1.

:- use_module(library(filesex), [directory_file_path/3]).
:- use_module(library(lists), [memberchk/2]).
:- use_module(library(readutil), [read_file_to_terms/3]).
:- use_module(library(statistics), [call_time/2]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir), file_type(directory)]),
   asserta(repo_root(Root)),
   directory_file_path(Root, 'load.pl', Load),
   ensure_loaded(Load),
   directory_file_path(Root, 'benchmarks/workloads.pl', Wl),
   consult(Wl).

:- initialization(main, main).

main :-
    % Heavy rungs reach tens of millions of inferences; give the (heap-backed)
    % Prolog stacks room under both VMs. The measured goal is unchanged.
    set_prolog_flag(stack_limit, 2_000_000_000),
    current_prolog_flag(argv, Argv),
    ( memberchk('--heavy', Argv) -> Tiers = [core, heavy] ; Tiers = [core] ),
    % gate=latency rungs are skipped: their count pins to the budget constant
    % (no ratchet-comparable search signal to diff) at a ~5e8-inference cost.
    forall( ( arrange_workload(F, Grid, _Mode, _It, Wu, Exp, Tier, inf, Budget, _Words),
              memberchk(Tier, Tiers) ),
            run_rung(F, Grid, Wu, Exp, Budget) ).

% One rung: read its word set, run `Warmup` discarded searches, then measure the
% inference count of the next search - the same warm/cold protocol the rung's
% baseline row used, so the number is directly comparable.
run_rung(F, Grid, Warmup, Expected, Budget) :-
    repo_root(Root),
    directory_file_path(Root, F, File),
    read_clues(File, Words),
    length(Words, NumWords),
    forall(between(1, Warmup, _),
           ignore(call_time(
               crosswordsmith_arrange:arrange_best_layout(
                   Words, Grid, Budget, _, _, _), _))),
    call_time(crosswordsmith_arrange:arrange_best_layout(Words, Grid, Budget, _, _, Outcome), T),
    Inf = T.inferences,
    file_base_name(F, Base),
    ( outcome_ok(Expected, Outcome) -> Tag = Outcome ; Tag = mismatch(Expected, Outcome) ),
    format("~w,~w,~w,~w,~w,~w~n", [Base, Grid, NumWords, Warmup, Inf, Tag]).

outcome_ok(placed, placed).
outcome_ok(infeasible, infeasible).
outcome_ok(infeasible, not_proven).

% Strict clue reader (same contract as benchmarks/run_arrange.pl): a fixture with
% no clues/1 is a hard error, never a bogus empty-word run. read_file_to_terms/3
% resolves under the WASM image too (library/readutil.pl ships in its tree).
read_clues(File, Words) :-
    read_file_to_terms(File, Terms, []),
    (   memberchk(clues(Words0), Terms)
    ->  Words = Words0
    ;   throw(error(fixture_missing_clues(File), _))
    ).
