#!/usr/bin/env swipl
% Fresh-process strict identity row generator. Search timing ends before the
% diagnostics-bearing layout is serialized, so identity work cannot consume
% the inference budget or alter the ratchet's search count.

:- set_prolog_flag(verbose, silent).
:- use_module(library(filesex), [directory_file_path/3]).
:- use_module(library(json), [json_write_dict/2]).
:- use_module(library(statistics), [call_time/2]).

:- dynamic repo_root/1.

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', RepoRoot,
                      [relative_to(BenchDir), file_type(directory), access(read)]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load),
   directory_file_path(BenchDir, 'workloads.pl', Workloads),
   consult(Workloads).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    (   Argv = [Fixture]
    ->  catch(identity_row(Fixture), E, (print_message(error, E), halt(1))),
        halt(0)
    ;   format(user_error, "usage: run_arrange_identity.pl FIXTURE~n", []),
        halt(2)
    ).

identity_row(Fixture0) :-
    text_to_atom(Fixture0, Fixture),
    findall(row(Size, Mode, Expected, Tier, Gate, Budget, ExpectedWords),
            arrange_workload(Fixture, Size, Mode, _Iters, _Warmup, Expected,
                             Tier, Gate, Budget, ExpectedWords),
            Rows),
    (   Rows = [row(Size, Mode, Expected, Tier, Gate, Budget, ExpectedWords)]
    ->  true
    ;   throw(error(identity_workload_not_unique(Fixture, Rows), _))
    ),
    mode_framing(Mode, Framing),
    setup_call_cleanup(
        clear_identity_state,
        identity_with_clear_state(Fixture, Size, Mode, Framing, Expected, Tier,
                                  Gate, Budget, ExpectedWords),
        clear_identity_state).

identity_with_clear_state(Fixture, Size, Mode, Framing, Expected, Tier, Gate,
                          Budget, ExpectedWords) :-
    repo_root(Root),
    directory_file_path(Root, Fixture, File),
    read_fixture(File, Words),
    length(Words, ActualWords),
    (   ActualWords =:= ExpectedWords
    ->  true
    ;   throw(error(identity_word_count(Fixture, ExpectedWords, ActualWords), _))
    ),
    measured_identity(Fixture, Size, Mode, Framing, Expected, Tier, Gate,
                      Budget, ExpectedWords, Words).

clear_identity_state :-
    set_search_seed(-1),
    set_check_target(-1).

measured_identity(Fixture, Size, Mode, Framing, Expected, Tier, Gate, Budget,
                  WordCount, Words) :-
    call_time(arrange_best_layout(Words, Size, Budget, Numbered, Reward, Outcome), _),
    expected_outcome(Expected, Outcome),
    format("fixture\t~w~nmode\t~w~nframing\t~w~ngrid\t~d~nwords\t~d~n",
           [Fixture, Mode, Framing, Size, WordCount]),
    format("tier\t~w~ngate\t~w~nbudget\t~d~noutcome\t~w~nreward\t~d~n",
           [Tier, Gate, Budget, Outcome, Reward]),
    emit_identity_payload(Outcome, Numbered, Words, Size, Framing).

emit_identity_payload(placed, Numbered, Words, Size, Framing) :-
    !,
    once(crosswordsmith_arrange:arrange_diag_layout_dict(
             Numbered, Words, Size, Framing, Payload)),
    format("payload\tjson~n", []),
    json_write_dict(user_output, Payload),
    nl.
emit_identity_payload(Outcome, _Numbered, _Words, _Size, _Framing) :-
    format("payload\tnon-placement:~w~n", [Outcome]).

mode_framing(size, fixed).
mode_framing(max_size, max).

text_to_atom(Text, Atom) :-
    ( atom(Text) -> Atom = Text ; atom_string(Atom, Text) ).

expected_outcome(placed, placed).
expected_outcome(infeasible, infeasible).
expected_outcome(infeasible, not_proven).

read_fixture(File, Words) :-
    read_file_to_terms(File, Terms, []),
    (   memberchk(clues(Words), Terms)
    ->  true
    ;   throw(error(identity_fixture_missing_clues(File), _))
    ).

:- multifile prolog:error_message//1.
prolog:error_message(identity_workload_not_unique(Fixture, Rows)) -->
    [ 'identity workload ~q resolved to ~q (expected exactly one row)'-[Fixture, Rows] ].
prolog:error_message(identity_word_count(Fixture, Expected, Actual)) -->
    [ 'identity fixture ~q has ~d words; manifest requires ~d'-[Fixture, Actual, Expected] ].
prolog:error_message(identity_fixture_missing_clues(File)) -->
    [ 'identity fixture ~q does not define clues/1'-[File] ].
