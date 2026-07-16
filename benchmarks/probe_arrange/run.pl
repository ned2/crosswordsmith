#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/2]).
:- use_module(library(process), [process_create/3, process_wait/2]).
:- use_module(library(time), [call_with_time_limit/2]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   nb_setval(probe_root, Root).

:- initialization(main, main).

% Usage is intentionally positional so batch launchers can construct commands
% without an option parser entering the measured process.
% authority-operation FIXTURE GRID COUNT FIXTURE_SEED SEARCH_SEED BUDGET OP ATT ARM OUTER
% authority-corner    ... BUDGET CORNER OP ATT ARM OUTER
% instrumented        ... CORNER MODE LIMIT_KIND LIMIT OP ATT ARM OUTER
main :-
    current_prolog_flag(argv, Argv),
    catch(dispatch(Argv), E, fatal(E)).

dispatch([Command|Args]) :-
    command_request(Command, Args, Request, Outer),
    heartbeat(Request, start),
    catch(call_with_time_limit(Outer, execute(Request, Row)),
          time_limit_exceeded,
          interrupted_row(Request, Row)),
    json_write_dict(user_output, Row, [width(0)]), nl,
    heartbeat(Request, done).
dispatch(_) :- throw(error(usage, _)).

command_request('authority-operation',
    [Fixture,G0,N0,FS0,SS0,B0,Op,A0,Arm,O0], Request, Outer) :-
    common(Fixture,G0,N0,FS0,SS0,Op,A0,Arm,O0,Request0,Outer),
    number_atom(B0, Budget),
    Request = Request0.put(_{command:authority_operation,budget:Budget,
                             corner:null,limit_kind:inferences,mode:none}).
command_request('authority-corner',
    [Fixture,G0,N0,FS0,SS0,B0,Corner,Op,A0,Arm,O0], Request, Outer) :-
    common(Fixture,G0,N0,FS0,SS0,Op,A0,Arm,O0,Request0,Outer),
    number_atom(B0, Budget),
    Request = Request0.put(_{command:authority_corner,budget:Budget,
                             corner:Corner,limit_kind:inferences,mode:none}).
command_request(instrumented,
    [Fixture,G0,N0,FS0,SS0,Corner,Mode,Kind,L0,Op,A0,Arm,O0], Request, Outer) :-
    common(Fixture,G0,N0,FS0,SS0,Op,A0,Arm,O0,Request0,Outer),
    parse_limit(Kind, L0, Limit),
    Request = Request0.put(_{command:instrumented,corner:Corner,mode:Mode,
                             limit_kind:Kind,limit:Limit}).

common(Fixture,G0,N0,FS0,SS0,Op,A0,Arm,O0,Request,Outer) :-
    number_atom(G0, Grid), number_atom(N0, Count), number_atom(FS0, FixtureSeed),
    search_seed(SS0, SearchSeed), number_atom(A0, Attempt), number_atom(O0, Outer),
    Request = _{fixture_path:Fixture,grid:Grid,count:Count,
                fixture_seed:FixtureSeed,search_seed:SearchSeed,
                operation_id:Op,attempt_index:Attempt,arm:Arm}.

parse_limit(none, _, null).
parse_limit(nodes, A, N) :- number_atom(A, N).
parse_limit(decisions, A, N) :- number_atom(A, N).
number_atom(A, N) :- atom_number(A, N).
search_seed(none, null) :- !.
search_seed(A, N) :- number_atom(A, N).

execute(Request, Row) :-
    probe_arrange:load_fixture(Request.fixture_path, Request.count, Words),
    execute_probe(Request, Words, Result, Stats, Timing, Rig),
    row_context(Request, Rig, Meta, Commit, Swi),
    probe_arrange:trace_row(Rig, Meta, Result, Stats, Timing, Commit, Swi, Row).

execute_probe(R, Words, Result, Stats, Timing, authority) :-
    R.command == authority_operation, !,
    probe_arrange:authority_operation(Words,R.grid,R.budget,R.search_seed,Result,Timing),
    null_stats(Stats).
execute_probe(R, Words, Result, Stats, Timing, authority) :-
    R.command == authority_corner, !,
    probe_arrange:authority_corner(Words,R.grid,R.corner,R.budget,R.search_seed,Result,Timing),
    null_stats(Stats).
execute_probe(R, Words, Result, Stats, Timing, instrumented) :-
    probe_arrange:twin_corner(Words,R.grid,R.corner,R.search_seed,R.mode,
                              R.limit_kind,R.limit,Result,run(Stats,Timing)).

null_stats(_{nodes:null,decisions:null,places:null,unplaces:null,wipeouts:null,
             max_depth:null,state_entries_max:null,letter_cells_max:null,
             boundary_cells_max:null,support_transitions:null,
             duplicate_children:null,duplicate_states:null}).

row_context(R, Rig, Meta, Commit, Swi) :-
    file_base_name(R.fixture_path, Fixture),
    Meta = _{limit_kind:R.limit_kind,operation_id:R.operation_id,
             attempt_index:R.attempt_index,fixture:Fixture,
             fixture_seed:R.fixture_seed,search_seed:R.search_seed,
             corner:R.corner,arm:R.arm},
    git_commit(Commit), swi_version(Swi),
    ( Rig == authority -> true ; true ).

git_commit(Commit) :-
    nb_getval(probe_root, Root),
    process_create(path(git), ['-C',Root,'rev-parse','HEAD'],
                   [stdout(pipe(S)), process(P)]),
    read_string(S, _, Raw), close(S), process_wait(P, exit(0)),
    normalize_space(string(Commit), Raw).
swi_version(Version) :-
    current_prolog_flag(version_data, swi(Major,Minor,Patch,_)),
    format(string(Version), '~d.~d.~d', [Major,Minor,Patch]).

interrupted_row(R, Row) :-
    interrupted_rig(R.command, Rig),
    row_context(R, Rig, Meta, Commit, Swi),
    Result = result(interrupted, interrupted, true, null, [], null),
    null_stats(Stats), Timing=_{inferences:null,wall:null,cpu:null},
    probe_arrange:trace_row(Rig, Meta, Result, Stats, Timing, Commit, Swi, Row).

interrupted_rig(instrumented, instrumented) :- !.
interrupted_rig(_, authority).

heartbeat(R, Phase) :-
    format(user_error, 'probe-arrange heartbeat operation=~w attempt=~w phase=~w~n',
           [R.operation_id,R.attempt_index,Phase]), flush_output(user_error).

fatal(E) :- print_message(error, E), halt(2).

:- multifile prolog:error_message//1.
prolog:error_message(usage) -->
    ['usage: run.pl authority-operation|authority-corner|instrumented ... (see README.md)'].
