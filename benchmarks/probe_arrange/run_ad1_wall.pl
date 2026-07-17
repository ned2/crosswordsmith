#!/usr/bin/env swipl
% Serialized alternating paired wall evidence for A-D1. Reporting-only.
:- set_prolog_flag(verbose, silent).
:- use_module(library(readutil), [read_file_to_terms/3]).
:- use_module(library(statistics), [call_time/2]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'ad1_buckets.pl', AD1), use_module(AD1).

:- initialization(main, main).

control(ladder_15x15_32w, 'fixtures/ladder_15x15_32w.pl', 15, 32).
control(ladder_21x21_80w, 'fixtures/ladder_21x21_80w.pl', 21, 80).

main :-
    forall(control(Name, Fixture, Grid, Count),
           run_control(Name, Fixture, Grid, Count)).

run_control(Name, Fixture, Grid, ExpectedCount) :-
    load_words(Fixture, ExpectedCount, Words),
    format(user_error, 'A-D1 wall ~w: warmup~n', [Name]),
    paired_sample(Words, Grid, assoc_first, _),
    paired_sample(Words, Grid, direct_first, _),
    findall(Pair,
            ( between(1, 21, Trial),
              ( 1 is Trial mod 2 -> Order = assoc_first ; Order = direct_first ),
              paired_sample(Words, Grid, Order, Pair),
              format(user_error, 'A-D1 wall ~w: trial ~d/21~n', [Name, Trial]) ),
            Pairs),
    findall(Assoc, member(pair(Assoc, _), Pairs), AssocTimes),
    findall(Direct, member(pair(_, Direct), Pairs), DirectTimes),
    findall(Ratio,
            ( member(pair(Assoc, Direct), Pairs), Ratio is Direct / Assoc ),
            Ratios),
    median(AssocTimes, AssocMedian),
    median(DirectTimes, DirectMedian),
    median(Ratios, RatioMedian),
    format('~w assoc_median=~6f direct_median=~6f paired_ratio=~6f~n',
           [Name, AssocMedian, DirectMedian, RatioMedian]).

paired_sample(Words, Grid, assoc_first, pair(AssocWall, DirectWall)) :-
    timed_operation(assoc, Words, Grid, AssocResult, AssocWall),
    timed_operation(direct, Words, Grid, DirectResult, DirectWall),
    AssocResult == DirectResult.
paired_sample(Words, Grid, direct_first, pair(AssocWall, DirectWall)) :-
    timed_operation(direct, Words, Grid, DirectResult, DirectWall),
    timed_operation(assoc, Words, Grid, AssocResult, AssocWall),
    AssocResult == DirectResult.

timed_operation(Mode, Words, Grid, Result, Wall) :-
    call_time(probe_arrange_ad1:ad1_operation(Mode, Words, Grid, Result), Time),
    Wall = Time.wall.

median(Values, Median) :-
    msort(Values, Sorted),
    length(Sorted, Length),
    Middle is Length // 2,
    nth0(Middle, Sorted, Median).

load_words(Fixture, ExpectedCount, Words) :-
    read_file_to_terms(Fixture, Terms, []),
    memberchk(clues(Words), Terms),
    length(Words, ExpectedCount).
