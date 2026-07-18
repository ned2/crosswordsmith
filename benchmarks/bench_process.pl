:- module(bench_process, [capture_process/6]).

:- use_module(library(error), [must_be/2]).
:- use_module(library(process), [process_create/3, process_wait/2]).
:- use_module(library(readutil), [read_file_to_string/3]).

%!  capture_process(+Exe, +Argv, +StderrMode, -Stdout, -Stderr, -Status) is det.
%
%   Run Exe with Argv and return its text output and process Status. StderrMode
%   is one of `capture`, `inherit`, or `null`; Stderr is empty unless capture is
%   selected. Standard input is inherited. Nonzero status is returned rather
%   than thrown. File-backed capture avoids bounded-pipe deadlocks while
%   setup_call_cleanup/3 closes streams, reaps the requested PID, and removes
%   temporary files. This predicate does not impose a timeout or kill a child.
capture_process(Exe, Argv, StderrMode, Stdout, Stderr, Status) :-
    must_be(oneof([capture, inherit, null]), StderrMode),
    setup_call_cleanup(
        capture_file(OutFile),
        capture_with_stderr(StderrMode, Exe, Argv, OutFile,
                            Stdout, Stderr, Status),
        delete_capture_file(OutFile)).

capture_with_stderr(capture, Exe, Argv, OutFile, Stdout, Stderr, Status) :-
    setup_call_cleanup(
        capture_file(ErrFile),
        ( run_capture(Exe, Argv, OutFile, file(ErrFile), Status),
          read_file_to_string(OutFile, Stdout, []),
          read_file_to_string(ErrFile, Stderr, [])
        ),
        delete_capture_file(ErrFile)).
capture_with_stderr(inherit, Exe, Argv, OutFile, Stdout, "", Status) :-
    run_capture(Exe, Argv, OutFile, inherit, Status),
    read_file_to_string(OutFile, Stdout, []).
capture_with_stderr(null, Exe, Argv, OutFile, Stdout, "", Status) :-
    run_capture(Exe, Argv, OutFile, null, Status),
    read_file_to_string(OutFile, Stdout, []).

capture_file(File) :-
    tmp_file_stream(text, File, Stream),
    close(Stream).

delete_capture_file(File) :-
    ( exists_file(File) -> delete_file(File) ; true ).

run_capture(Exe, Argv, OutFile, ErrMode, Status) :-
    setup_call_cleanup(
        open(OutFile, write, Out),
        run_with_stderr(ErrMode, Exe, Argv, Out, Status),
        close_quietly(Out)).

run_with_stderr(file(ErrFile), Exe, Argv, Out, Status) :-
    setup_call_cleanup(
        open(ErrFile, write, Err),
        run_child(Exe, Argv, Out, stream(Err), Status),
        close_quietly(Err)).
run_with_stderr(inherit, Exe, Argv, Out, Status) :-
    run_child(Exe, Argv, Out, std, Status).
run_with_stderr(null, Exe, Argv, Out, Status) :-
    run_child(Exe, Argv, Out, null, Status).

run_child(Exe, Argv, Out, ErrSpec, Status) :-
    setup_call_cleanup(
        start_child(Exe, Argv, Out, ErrSpec, Child),
        wait_child(Child, Status),
        reap_child(Child)).

start_child(Exe, Argv, Out, ErrSpec, child(PID, pending)) :-
    process_create(Exe, Argv,
                   [stdout(stream(Out)), stderr(ErrSpec), process(PID)]).

wait_child(Child, Status) :-
    arg(1, Child, PID),
    process_wait(PID, Actual),
    nb_setarg(2, Child, waited),
    Status = Actual.

reap_child(child(PID, pending)) :-
    !,
    catch(process_wait(PID, _), _, true).
reap_child(_).

close_quietly(Stream) :-
    catch(close(Stream), _, true).
