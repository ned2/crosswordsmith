:- module(bench_report,
          [ benchmark_report/3,
            swi_version/1
          ]).

:- use_module(library(error), [must_be/2]).

%!  benchmark_report(+Tool:atom, +Rows:list, -Report:dict) is det.
%
%   Construct the common JSON report envelope for a nonempty benchmark run.
%   Domain runners add their own metric metadata and retain their row schemas.
benchmark_report(Tool, Rows, Report) :-
    must_be(atom, Tool),
    must_be(list, Rows),
    ( Rows = [_|_] -> true
    ; throw(error(bench_report_empty_results(Tool), _))
    ),
    swi_version(Version),
    Report = _{tool:Tool, swi_prolog:Version, results:Rows}.

%!  swi_version(-Version:atom) is det.
%
%   Return the running SWI-Prolog semantic version used in benchmark metadata.
swi_version(Version) :-
    current_prolog_flag(version_data, Value),
    ( Value = swi(Major, Minor, Patch, _)
    -> format(atom(Version), '~d.~d.~d', [Major, Minor, Patch])
    ;  Version = Value
    ).

:- multifile prolog:error_message//1.
prolog:error_message(bench_report_empty_results(Tool)) -->
    [ '~w benchmark report has no result rows'-[Tool] ].
