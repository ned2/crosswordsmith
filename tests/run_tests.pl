% tests/run_tests.pl - plunit test runner for the crosswordsmith library.
%
% Loads the code under test (via load.pl, the single owner of load order) and
% the test suite, runs all plunit tests, and exits with a non-zero code if any
% test fails or times out. Intended to be invoked from the repository root
% (run_tests.sh / make test handle this):
%
%     swipl -q tests/run_tests.pl

:- use_module(library(plunit)).
% Named run_suite rather than main so it can never clash with a main/0 in
% anything we consult below.
:- initialization(run_suite, main).

run_suite :-
    consult('load.pl'),                % the whole implementation, known-good order
    consult('tests/crossword.plt'),
    consult('tests/arrange.plt'),
    consult('tests/lint.plt'),
    consult('tests/export.plt'),
    consult('tests/stockgrid.plt'),
    consult('tests/fill.plt'),
    consult('tests/browser.plt'),
    % summary(S) gives us the counts AND stops run_tests/2 from failing on
    % its own, so we control the exit code explicitly from the dict. Note:
    % run_tests/2 and the summary/1 option are undocumented in the SWI 10.0.2
    % manual (which lists only run_tests/0,1) but are verified working here.
    run_tests(all, [summary(S)]),
    format("~n== summary: ~w passed, ~w failed, ~w timed out, ~w blocked (of ~w) ==~n",
           [S.passed, S.failed, S.timeout, S.blocked, S.total]),
    ( S.failed =:= 0, S.timeout =:= 0
    -> halt(0)
    ;  halt(1)
    ).
