% benchmarks/fill_subjects.pl - the fill-specific measurement subjects.
%
% This is the ONLY fill-bench file that knows what a "fill" is. It adapts the
% generic samplers in bench_core into measure/3-ready closures for the FOUR
% attribution buckets the fill product bench reports (campaign plan Phase 0):
%
%   - command  : the whole `crosswordsmith fill` process, end-to-end (SWI startup
%                + load.pl + dict load + slot derivation + search + emit) - the
%                latency a user feels. process layer; asserts exit vs Expected.
%   - dict_load: the in-process load_dict/3 alone (the every-invocation startup
%                tax: read + normalize + positional index over all lengths).
%                load_inf (inferences) is deterministic and REPORTED (not gated).
%   - grid     : the in-process fill_grid/4 alone (mask -> shared-cell-variable
%                slots). Its own bucket, never silently folded into "rest".
%   - search   : the in-process budget-explicit fill_attempt/8 with a PRE-LOADED
%                dict and FRESH slots rebuilt per sample. search_inf is the metric
%                of record (GATED). Asserts Outcome vs Expected every iteration.
%
% CRITICAL (the fresh-slots trap): fill_attempt destructively unifies the shared
% cell variables of its slots. Timing N iterations against ONE slot set would run
% iterations 2..N on an already-solved grid - garbage counts. So the search
% sampler REBUILDS fill_grid (+ seeds) OUTSIDE the timed goal on every sample, and
% times ONLY fill_attempt/8. The dict (DictByLen+Index) is loaded ONCE by the
% caller and threaded in, so search cost is isolated from load cost. The slots are
% passed to fill_attempt as plain arguments (never through findall/yall) so the
% crossing cells stay genuinely shared (fill.pl:57-64 comment).

:- module(fill_subjects,
          [ fill_command_sampler/6,
            fill_load_sampler/2,
            fill_grid_sampler/2,
            fill_search_sampler/6,
            build_search_slots/3 ]).

:- use_module(library(lists)).
:- use_module(library(apply)).

% COMMAND layer. process_sampler returns the raw exit code; we assert it against
% Expected here (not inside bench_core) and record only wall+rss. rss (peak RSS,
% KiB) is a whole-process footprint, not a search-memory metric.
fill_command_sampler(Exe, GridFile, DictFile, Seeds, Expected, _{wall:Wall, rss:Rss}) :-
    seed_args(Seeds, SeedArgs),
    append(['fill', '--grid', GridFile, '--dict', DictFile], SeedArgs, A0),
    append(A0, ['--out', '/dev/null'], Argv),
    bench_core:process_sampler(Exe, Argv, Wall, Rss, Exit),
    expected_exit(Expected, Exit).

seed_args(none, []) :- !.
seed_args(SeedFile, ['--seeds', SeedFile]).

% fill exits 0 on a full fill; non-zero on any non-filled outcome (INV-3). The
% ladder rungs all expect `filled`.
expected_exit(filled, 0) :- !.
expected_exit(Expected, Got) :-
    throw(error(bench_fill_exit(expected(Expected), exit(Got)), _)).

% DICT-LOAD layer. Each sample reloads the whole dictionary (fresh index); the
% inference count is deterministic (the reported load_inf). We do NOT wrap in
% once/1 (adds an inference; bench_core's inproc_sampler note).
fill_load_sampler(DictFile, Sample) :-
    bench_core:inproc_sampler(
        crosswordsmith_fill:load_dict(DictFile, _DictByLen, _Index),
        Sample).

% GRID layer. Each sample re-derives the slots from the mask (fresh shared cell
% variables); grid_inf is deterministic.
fill_grid_sampler(GridFile, Sample) :-
    bench_core:inproc_sampler(
        crosswordsmith_fill:fill_grid(GridFile, _Size, _Slots, _CellVar),
        Sample).

% SEARCH layer. Rebuild FRESH slots (+ seeds) OUTSIDE the timed goal, then time
% ONLY fill_attempt/8 against the caller's PRE-LOADED dict. Assert the outcome.
fill_search_sampler(GridFile, Seeds, DictByLen, Index, Budget-Expected, Sample) :-
    build_search_slots(GridFile, Seeds, slots(SearchSlots, AllSlots)),
    bench_core:inproc_sampler(
        crosswordsmith_fill:fill_attempt(SearchSlots, AllSlots, DictByLen, Index,
                                         Budget, Outcome, _Numbered, _InputWords),
        Sample),
    ( Outcome == Expected -> true
    ; throw(error(bench_fill_outcome(expected(Expected), got(Outcome)), _)) ).

% Fresh slot set from the mask (+ seed pins). AllSlots are emitted (seed pins
% included); SearchSlots are what the engine fills (AllSlots minus seed pins).
% Slots are passed as plain args downstream, so the crossing cells stay shared.
build_search_slots(GridFile, Seeds, slots(SearchSlots, AllSlots)) :-
    crosswordsmith_fill:fill_grid(GridFile, Size, AllSlots, _CellVar),
    ( Seeds == none
    ->  SearchSlots = AllSlots
    ;   crosswordsmith_fill:apply_seeds(Seeds, Size, AllSlots, SeededKeys),
        exclude(crosswordsmith_fill:seeded_slot(SeededKeys), AllSlots, SearchSlots) ).

:- multifile prolog:error_message//1.
prolog:error_message(bench_fill_exit(expected(E), exit(G))) -->
    [ 'fill command exit ~w does not match expected workload outcome ~w'-[G, E] ].
prolog:error_message(bench_fill_outcome(expected(E), got(G))) -->
    [ 'fill search outcome ~w does not match expected workload outcome ~w'-[G, E] ].
