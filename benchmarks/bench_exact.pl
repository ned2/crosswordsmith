:- module(bench_exact,
          [ exact_version/4,
            exact_metric/4,
            exact_presence/5
          ]).

:- use_module(library(apply), [maplist/3]).
:- use_module(library(ordsets), [ord_subtract/3]).

%!  exact_version(+Reference, +Measured, -Status, -Failures) is det.
%
%   Require the same textual SWI-Prolog version for an exact comparison.
exact_version(Reference, Measured, Status, Failures) :-
    text_to_string(Reference, ReferenceString),
    text_to_string(Measured, MeasuredString),
    ( ReferenceString == MeasuredString
    -> Status = same, Failures = 0
    ;  Status = version_mismatch, Failures = 1
    ).

%!  exact_metric(+Reference:number, +Measured:number, -Status, -Failures) is det.
%
%   Compare one gated metric. Both increases and decreases are failures.
exact_metric(Reference, Measured, Status, Failures) :-
    ( Measured =:= Reference
    -> Status = same, Failures = 0
    ; Measured > Reference
    -> Status = increase, Failures = 1
    ;  Status = decrease, Failures = 1
    ).

%!  exact_presence(+ReferenceIds, +MeasuredIds, -Missing, -Unexpected, -Failures) is det.
%
%   Compare complete row-id sets, accepting atom/string spelling differences.
exact_presence(ReferenceIds, MeasuredIds, Missing, Unexpected, Failures) :-
    maplist(text_to_string, ReferenceIds, ReferenceStrings0),
    maplist(text_to_string, MeasuredIds, MeasuredStrings0),
    sort(ReferenceStrings0, ReferenceStrings),
    sort(MeasuredStrings0, MeasuredStrings),
    ord_subtract(ReferenceStrings, MeasuredStrings, Missing),
    ord_subtract(MeasuredStrings, ReferenceStrings, Unexpected),
    length(Missing, MissingCount),
    length(Unexpected, UnexpectedCount),
    Failures is MissingCount + UnexpectedCount.
