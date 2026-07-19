:- module(bench_exact,
          [ exact_version/4,
             require_same_version/2,
             require_complete_migration/4,
             require_protocol_value/4,
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

%!  require_same_version(+Reference, +Measured) is det.
%
%   Succeed when the textual SWI-Prolog versions match; otherwise throw before
%   a caller can persist incomparable measurements.
require_same_version(Reference, Measured) :-
    exact_version(Reference, Measured, _Status, Failures),
    ( Failures =:= 0 -> true
    ; throw(error(bench_swi_version_mismatch(Reference, Measured), _))
    ).

%!  require_complete_migration(+ReferenceVersion, +MeasuredVersion,
%!                             +ReferenceIds:list, +MeasuredIds:list) is det.
%
%   A same-version record may retain unmeasured rows. When the SWI-Prolog version
%   changes, require the complete reference row set so retained rows cannot be
%   relabeled with counts from the old runtime.
require_complete_migration(ReferenceVersion, MeasuredVersion,
                           ReferenceIds, MeasuredIds) :-
    exact_version(ReferenceVersion, MeasuredVersion, _Status, VersionFailures),
    ( VersionFailures =:= 0
    -> true
    ; exact_presence(ReferenceIds, MeasuredIds,
                     Missing, Unexpected, PresenceFailures),
      ( PresenceFailures =:= 0 -> true
      ; throw(error(bench_version_migration_incomplete(Missing, Unexpected), _))
      )
    ).

%!  require_protocol_value(+RowId, +Field:atom, +Reference, +Measured) is det.
%
%   Require one persisted measurement-protocol value to match its measured row.
%   Numbers compare arithmetically and text accepts atom/string spelling; all
%   other values compare structurally. A mismatch throws before callers
%   can record metrics gathered under a different protocol.
require_protocol_value(RowId, Field, Reference, Measured) :-
    ( protocol_value_equal(Reference, Measured)
    -> true
    ; throw(error(bench_record_protocol_mismatch(
                      RowId, Field, Reference, Measured), _))
    ).

protocol_value_equal(Reference, Measured) :-
    number(Reference),
    number(Measured),
    !,
    Reference =:= Measured.
protocol_value_equal(Reference, Measured) :-
    text(Reference),
    text(Measured),
    !,
    text_to_string(Reference, ReferenceString),
    text_to_string(Measured, MeasuredString),
    ReferenceString == MeasuredString.
protocol_value_equal(Reference, Measured) :-
    is_dict(Reference),
    is_dict(Measured),
    !,
    dict_pairs(Reference, _, ReferencePairs),
    dict_pairs(Measured, _, MeasuredPairs),
    protocol_pairs_equal(ReferencePairs, MeasuredPairs).
protocol_value_equal(Reference, Measured) :-
    is_list(Reference),
    is_list(Measured),
    !,
    protocol_lists_equal(Reference, Measured).
protocol_value_equal(Reference, Measured) :-
    Reference =@= Measured.

protocol_pairs_equal([], []).
protocol_pairs_equal([Key-Reference|ReferencePairs],
                     [Key-Measured|MeasuredPairs]) :-
    protocol_value_equal(Reference, Measured),
    protocol_pairs_equal(ReferencePairs, MeasuredPairs).

protocol_lists_equal([], []).
protocol_lists_equal([Reference|References], [Measured|Measureds]) :-
    protocol_value_equal(Reference, Measured),
    protocol_lists_equal(References, Measureds).

text(Value) :- atom(Value), !.
text(Value) :- string(Value).

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
    length(MeasuredStrings0, MeasuredCount),
    length(MeasuredStrings, UniqueMeasuredCount),
    ( MeasuredCount =:= UniqueMeasuredCount -> true
    ; throw(error(bench_exact_duplicate_ids(MeasuredStrings0), _))
    ),
    ord_subtract(ReferenceStrings, MeasuredStrings, Missing),
    ord_subtract(MeasuredStrings, ReferenceStrings, Unexpected),
    length(Missing, MissingCount),
    length(Unexpected, UnexpectedCount),
    Failures is MissingCount + UnexpectedCount.

:- multifile prolog:error_message//1.
prolog:error_message(bench_exact_duplicate_ids(Ids)) -->
    [ 'benchmark exact comparison received duplicate row IDs: ~q'-[Ids] ].
prolog:error_message(bench_swi_version_mismatch(Reference, Measured)) -->
    [ 'benchmark reference SWI-Prolog version ~w differs from measured version ~w'-[Reference, Measured] ].
prolog:error_message(bench_version_migration_incomplete(Missing, Unexpected)) -->
    [ 'benchmark SWI-version migration is incomplete (missing ~q, unexpected ~q)'-[Missing, Unexpected] ].
prolog:error_message(bench_record_protocol_mismatch(RowId, Field, Reference, Measured)) -->
    [ 'benchmark row ~w changes recording protocol ~w from ~q to ~q'-
      [RowId, Field, Reference, Measured] ].
