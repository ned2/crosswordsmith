:- module(bench_cli,
          [ validate_runner_options/5,
            require_selected/3,
            checker_mode/3
          ]).

:- use_module(library(apply), [exclude/3]).

%!  validate_runner_options(+Tool, +Positional, +Format, +Iterations, +Warmup) is det.
%
%   Validate common benchmark runner options before any workload is measured.
validate_runner_options(Tool, Positional, Format, Iterations, Warmup) :-
    ( Positional == [] -> true
    ; throw(error(bench_positional_args(Tool, Positional), _)) ),
    ( memberchk(Format, [text, csv, json]) -> true
    ; throw(error(bench_bad_format(Tool, Format), _)) ),
    valid_override(Tool, iterations, Iterations, 1),
    valid_override(Tool, warmup, Warmup, 0).

valid_override(_Tool, _Key, -1, _Min) :- !.
valid_override(Tool, Key, Value, Min) :-
    ( integer(Value), Value >= Min -> true
    ; throw(error(bench_bad_count(Tool, Key, Value, Min), _)) ).

%!  require_selected(+Tool, +ExplicitFilter, +Selected) is det.
%
%   Reject an explicit filter that selects no workload.
require_selected(_Tool, '', _Selected) :- !.
require_selected(Tool, Filter, Selected) :-
    ( Selected = [_|_] -> true
    ; throw(error(bench_empty_selection(Tool, Filter), _)) ).

%!  checker_mode(+Argv, -Mode, -RunnerArgs) is det.
%
%   Separate one optional checker mode from arguments forwarded to the runner.
%   Multiple mode flags conflict. History accepts no runner arguments because it
%   performs no measurement.
checker_mode(Argv, Mode, RunnerArgs) :-
    findall(M, (member(Arg, Argv), checker_mode_arg(Arg, M)), Modes),
    ( Modes == [] -> Mode = check
    ; Modes = [Mode] -> true
    ; throw(error(bench_conflicting_modes(Modes), _)) ),
    exclude(is_checker_mode_arg, Argv, RunnerArgs),
    ( Mode == history, RunnerArgs \== []
    -> throw(error(bench_history_args(RunnerArgs), _))
    ; true ).

checker_mode_arg('--history', history).
checker_mode_arg('--log', log).
checker_mode_arg('--record', record).
checker_mode_arg('--promote', promote).

is_checker_mode_arg(Arg) :-
    checker_mode_arg(Arg, _).

:- multifile prolog:error_message//1.
prolog:error_message(bench_positional_args(Tool, Args)) -->
    [ '~w benchmark does not accept positional arguments: ~q'-[Tool, Args] ].
prolog:error_message(bench_bad_format(Tool, Format)) -->
    [ '~w benchmark format must be text, csv, or json; got ~q'-[Tool, Format] ].
prolog:error_message(bench_bad_count(Tool, Key, Value, Min)) -->
    [ '~w benchmark ~w must be at least ~d; got ~q'-[Tool, Key, Min, Value] ].
prolog:error_message(bench_empty_selection(Tool, Filter)) -->
    [ '~w benchmark filter ~q selected no workloads'-[Tool, Filter] ].
prolog:error_message(bench_conflicting_modes(Modes)) -->
    [ 'benchmark checker modes conflict: ~q'-[Modes] ].
prolog:error_message(bench_history_args(Args)) -->
    [ 'benchmark history mode accepts no additional arguments: ~q'-[Args] ].
