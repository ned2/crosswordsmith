:- module(bench_fixture, [load_arrange_fixture/2]).

:- use_module(library(error), [must_be/2]).
:- use_module(library(lists), [memberchk/2]).
:- use_module(library(readutil), [read_file_to_terms/3]).

%!  load_arrange_fixture(+File:atom, -Words:list) is det.
%
%   Read a Prolog arrange fixture containing one nonempty clues/1 payload.
%   Reading terms rather than consulting the file prevents fixture predicates
%   from entering the harness module.
load_arrange_fixture(File, Words) :-
    read_file_to_terms(File, Terms, []),
    (   memberchk(clues(Words0), Terms)
    ->  must_be(list, Words0),
        ( Words0 = [_|_] -> Words = Words0
        ; throw(error(bench_fixture_empty_clues(File), _))
        )
    ;   throw(error(bench_fixture_missing_clues(File), _))
    ).

:- multifile prolog:error_message//1.
prolog:error_message(bench_fixture_missing_clues(File)) -->
    [ 'benchmark fixture ~q does not define clues/1'-[File] ].
prolog:error_message(bench_fixture_empty_clues(File)) -->
    [ 'benchmark fixture ~q defines an empty clues/1 payload'-[File] ].
