:- module(bench_cli,
           [ parse_runner_options/5,
             require_selected/3,
             workload_selected/4,
              apply_override/3,
              checker_mode/3,
              exact_runner_args/2,
              require_persistence_args/1,
              require_unique_ids/2
           ]).

:- use_module(library(apply), [exclude/3, maplist/3]).
:- use_module(library(error), [must_be/2]).
:- use_module(library(lists), [member/2, memberchk/2]).
:- use_module(library(optparse), [opt_parse/4]).

%!  parse_runner_options(+Tool, +Spec:list, +Argv:list, -Options:list,
%!                       -Common:dict) is det.
%
%   Parse a runner's domain-owned option Spec and extract its common options.
%   Spec must define help, format, fixture, heavy, iterations, and warmup.
%   Unknown options and invalid common values throw before measurement.
parse_runner_options(Tool, Spec, Argv, Options, Common) :-
    opt_parse(Spec, Argv, Options, Positional),
    memberchk(help(Help), Options),
    memberchk(format(Format), Options),
    memberchk(fixture(Filter), Options),
    memberchk(heavy(Heavy), Options),
    memberchk(iterations(Iterations), Options),
    memberchk(warmup(Warmup), Options),
    validate_runner_options(Tool, Positional, Format, Iterations, Warmup),
    Common = runner_options{
        help:Help, format:Format, fixture:Filter, heavy:Heavy,
        iterations:Iterations, warmup:Warmup
    }.

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

%!  require_selected(+Tool, +Filter, +Selected:list) is det.
%
%   Reject any runner selection with no workload, including a malformed
%   manifest whose default core selection is empty.
require_selected(Tool, Filter, Selected) :-
    must_be(list, Selected),
    ( Selected = [_|_] -> true
    ; throw(error(bench_empty_selection(Tool, Filter), _)) ).

%!  workload_selected(+Filter:atom, +Heavy:boolean, +Id:atom, +Tier:atom) is semidet.
%
%   Apply the common runner selection policy. An explicit substring filter
%   selects matching IDs from any tier; otherwise core is selected by default
%   and heavy is added only when Heavy is true. At most one proof is produced
%   even when Filter occurs repeatedly in Id.
workload_selected(Filter, _Heavy, Id, _Tier) :-
    Filter \== '',
    !,
    once(sub_atom(Id, _, _, _, Filter)).
workload_selected('', _Heavy, _Id, core) :-
    !.
workload_selected('', true, _Id, heavy).

%!  apply_override(+Override:integer, +Default:integer, -Value:integer) is det.
%
%   Resolve the runners' -1 sentinel to the manifest Default; any validated
%   explicit override replaces it.
apply_override(-1, Default, Default) :-
    !.
apply_override(Override, _Default, Override).

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
checker_mode_arg('--exact', exact).

is_checker_mode_arg(Arg) :-
    checker_mode_arg(Arg, _).

%!  exact_runner_args(+RunnerArgs, -FullRunnerArgs) is det.
%
%   Exact comparison always covers the complete core and heavy ladder. Reject
%   selection and sampling overrides that would make the proof partial.
exact_runner_args(RunnerArgs, ['--heavy']) :-
    exclude(==('--heavy'), RunnerArgs, OtherArgs),
    ( OtherArgs == [] -> true
    ; throw(error(bench_exact_args(OtherArgs), _)) ).

%!  require_persistence_args(+RunnerArgs:list) is det.
%
%   Reject warmup and iteration overrides before a checker records a baseline or
%   history row. Selection flags remain valid, but persisted measurements must
%   use the workload manifest's sampling protocol, including for first-seen rows.
require_persistence_args(RunnerArgs) :-
    ( member(Arg, RunnerArgs), sampling_override_arg(Arg)
    -> throw(error(bench_persistence_sampling_override(Arg), _))
    ; true
    ).

%!  require_unique_ids(+Tool:atom, +Ids:list) is det.
%
%   Reject duplicate selected or measured row IDs, accepting atom/string spelling
%   differences. Duplicate rows make both comparison totals and persistence
%   provenance ambiguous.
require_unique_ids(Tool, Ids) :-
    maplist(text_to_string, Ids, Strings),
    sort(Strings, Unique),
    length(Strings, Count),
    length(Unique, UniqueCount),
    ( Count =:= UniqueCount
    -> true
    ; throw(error(bench_duplicate_ids(Tool, Strings), _))
    ).

sampling_override_arg('--iterations').
sampling_override_arg('--warmup').
sampling_override_arg(Arg) :-
    atom(Arg),
    ( sub_atom(Arg, 0, _, _, '--iterations=')
    ; sub_atom(Arg, 0, _, _, '--warmup=')
    ).

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
prolog:error_message(bench_exact_args(Args)) -->
    [ 'benchmark exact mode accepts no selection or sampling overrides: ~q'-[Args] ].
prolog:error_message(bench_persistence_sampling_override(Arg)) -->
    [ 'benchmark persistence requires manifest sampling; remove ~w'-[Arg] ].
prolog:error_message(bench_duplicate_ids(Tool, Ids)) -->
    [ '~w benchmark contains duplicate row IDs: ~q'-[Tool, Ids] ].
