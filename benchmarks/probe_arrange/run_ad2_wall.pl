#!/usr/bin/env swipl
% Serialized alternating paired wall evidence for A-D2. Reporting-only.
:- set_prolog_flag(verbose, silent).
:- use_module(library(readutil), [read_file_to_terms/3]).
:- use_module(library(statistics), [call_time/2]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'ad1_buckets.pl', AD1), use_module(AD1),
   directory_file_path(Dir, 'd0_support.pl', D0), use_module(D0),
   directory_file_path(Dir, 'ad2_delta.pl', AD2), use_module(AD2).

:- initialization(main, main).

control(ladder_15x15_32w, 'fixtures/ladder_15x15_32w.pl', 15, 32).
control(ladder_21x21_80w, 'fixtures/ladder_21x21_80w.pl', 21, 80).

main :-
    forall(control(Name, Fixture, Grid, Count),
           run_control(Name, Fixture, Grid, Count)).

run_control(Name, Fixture, Grid, ExpectedCount) :-
    load_words(Fixture, ExpectedCount, Words),
    format(user_error, 'A-D2 wall ~w: warmup~n', [Name]),
    paired_sample(Words, Grid, full_first, _),
    paired_sample(Words, Grid, delta_first, _),
    findall(Pair,
            ( between(1, 21, Trial),
              ( 1 is Trial mod 2 -> Order = full_first ; Order = delta_first ),
              paired_sample(Words, Grid, Order, Pair),
              format(user_error, 'A-D2 wall ~w: trial ~d/21~n', [Name, Trial]) ),
            Pairs),
    findall(Full, member(pair(Full, _), Pairs), FullTimes),
    findall(Delta, member(pair(_, Delta), Pairs), DeltaTimes),
    findall(Ratio,
            ( member(pair(Full, Delta), Pairs), Ratio is Delta / Full ),
            Ratios),
    median(FullTimes, FullMedian),
    median(DeltaTimes, DeltaMedian),
    median(Ratios, RatioMedian),
    format('~w full_median=~6f delta_median=~6f paired_ratio=~6f~n',
           [Name, FullMedian, DeltaMedian, RatioMedian]).

paired_sample(Words, Grid, full_first, pair(FullWall, DeltaWall)) :-
    timed_full(Words, Grid, FullResult, FullWall),
    timed_delta(Words, Grid, DeltaResult, DeltaWall),
    FullResult == DeltaResult.
paired_sample(Words, Grid, delta_first, pair(FullWall, DeltaWall)) :-
    timed_delta(Words, Grid, DeltaResult, DeltaWall),
    timed_full(Words, Grid, FullResult, FullWall),
    FullResult == DeltaResult.

timed_full(Words, Grid, Result, Wall) :-
    call_time(probe_arrange_ad1:ad1_operation(direct, Words, Grid, Result), T),
    Wall = T.wall.
timed_delta(Words, Grid, Result, Wall) :-
    call_time(probe_arrange_ad2:ad2_operation(Words, Grid, Result), T),
    Wall = T.wall.

median(Values, Median) :-
    msort(Values, Sorted), length(Sorted, N), Middle is N // 2,
    nth0(Middle, Sorted, Median).

load_words(Fixture, ExpectedCount, Words) :-
    read_file_to_terms(Fixture, Terms, []), memberchk(clues(Words), Terms),
    length(Words, ExpectedCount).
