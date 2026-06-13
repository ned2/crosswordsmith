% tests/run_tests.pl - plunit test runner for crossword.pl
%
% Loads the code under test and the test suite, runs all plunit tests, and
% exits with a non-zero code if any test fails or times out. Intended to be
% invoked from the repository root (run_tests.sh / make test handle this):
%
%     swipl -q tests/run_tests.pl
%
% Using initialization(main, main) means consulting crossword.pl below does
% NOT trigger its own main/0 (see the comment at the top of crossword.pl).

:- use_module(library(plunit)).
% Named run_suite rather than main to avoid clashing with crossword.pl's
% own main/0 when we consult it below.
:- initialization(run_suite, main).

run_suite :-
    consult('crossword.pl'),
    consult('tests/crossword.plt'),
    % summary(S) gives us the counts AND stops run_tests/2 from failing on
    % its own, so we control the exit code explicitly from the dict.
    run_tests(all, [summary(S)]),
    format("~n== summary: ~w passed, ~w failed, ~w timed out, ~w blocked (of ~w) ==~n",
           [S.passed, S.failed, S.timeout, S.blocked, S.total]),
    ( S.failed =:= 0, S.timeout =:= 0
    -> halt(0)
    ;  halt(1)
    ).
