:- module(bench_store,
          [ read_json_dict/2,
            build_recorded_baseline/5,
            replace_json_dict/4,
            append_history/4,
            read_history/2,
            render_history/3,
            current_host/1
          ]).

:- use_module(library(apply), [foldl/4]).
:- use_module(library(error), [must_be/2]).
:- use_module(library(json), [json_read_dict/3, json_write_dict/3]).
:- use_module(library(lists),
              [ append/3,
                last/2,
                member/2,
                memberchk/2,
                nth1/3,
                reverse/2
              ]).
:- use_module(library(readutil), [read_line_to_string/2]).
:- use_module('bench_process.pl', [capture_process/6]).

:- meta_predicate
    build_recorded_baseline(+, +, 2, 3, -),
    replace_json_dict(+, +, +, 1).

%!  read_json_dict(+Path:atom, -Dict:dict) is det.
%
%   Read exactly one UTF-8 JSON object. Missing, malformed, non-object, or
%   trailing non-whitespace input throws a contextual exception.
read_json_dict(Path, Dict) :-
    ( exists_file(Path) -> true
    ; throw(error(bench_json_missing(Path), _))
    ),
    setup_call_cleanup(
        open(Path, read, Stream, [encoding(utf8)]),
        read_json_stream(Path, Stream, Dict),
        close(Stream)).

read_json_stream(Context, Stream, Dict) :-
    catch(json_read_dict(Stream, Term, [default_tag(json)]),
          Error,
          throw(error(bench_json_invalid(Context, Error), _))),
    ( is_dict(Term) -> Dict = Term
    ; throw(error(bench_json_not_object(Context), _))
    ),
    read_string(Stream, _, Rest),
    ( blank_string(Rest) -> true
    ; throw(error(bench_json_trailing_data(Context), _))
    ).

%!  build_recorded_baseline(+Baseline:dict, +Run:dict, :RowKey,
%!                          :RowSpec, -Recorded:dict) is det.
%
%   Merge Run.results into Baseline.workloads while retaining every unmeasured
%   workload. RowKey(+Row, -Key) identifies a row. RowSpec(+Row, +Prior, -Spec)
%   builds its persisted spec, where Prior is `existing(OldSpec)` or `new`.
%   Empty runs and duplicate measured keys are rejected. The run's SWI version
%   and the current host replace their baseline metadata counterparts.
build_recorded_baseline(Baseline, Run, RowKey, RowSpec, Recorded) :-
    must_be(dict, Baseline),
    must_be(dict, Run),
    get_dict(results, Run, Rows),
    must_be(list, Rows),
    ( Rows = [_|_] -> true
    ; throw(error(bench_record_empty_results, _))
    ),
    get_dict(workloads, Baseline, Workloads0),
    must_be(dict, Workloads0),
    foldl(merge_recorded_row(RowKey, RowSpec), Rows,
          state(Workloads0, []), state(Workloads, _Seen)),
    get_dict(swi_prolog, Run, Swi),
    current_host(Host),
    Recorded = Baseline.put(_{host:Host, swi_prolog:Swi,
                              workloads:Workloads}).

merge_recorded_row(RowKey, RowSpec, Row,
                   state(Workloads0, Seen0), state(Workloads, [Key|Seen0])) :-
    call(RowKey, Row, Key0),
    text_to_atom(Key0, Key),
    ( memberchk(Key, Seen0)
    -> throw(error(bench_record_duplicate_key(Key), _))
    ; true
    ),
    ( get_dict(Key, Workloads0, Old) -> Prior = existing(Old) ; Prior = new ),
    call(RowSpec, Row, Prior, Spec),
    must_be(dict, Spec),
    put_dict(Key, Workloads0, Spec, Workloads).

%!  replace_json_dict(+Path:atom, +Dict:dict, +Width:integer,
%!                    :VerifyReadBack) is det.
%
%   Write Dict to a sibling temporary file, parse it completely, call
%   VerifyReadBack(+Parsed), and atomically rename it over Path only when all
%   steps succeed. The temporary file is removed on failure or exception.
replace_json_dict(Path, Dict, Width, VerifyReadBack) :-
    must_be(dict, Dict),
    must_be(integer, Width),
    atom_concat(Path, '.tmp', TempPath),
    setup_call_cleanup(
        true,
        ( setup_call_cleanup(
              open(TempPath, write, Stream, [encoding(utf8)]),
              json_write_dict(Stream, Dict, [width(Width)]),
              close(Stream)),
          read_json_dict(TempPath, Parsed),
          ( call(VerifyReadBack, Parsed) -> true
          ; throw(error(bench_record_verification_failed(Path), _))
          ),
          rename_file(TempPath, Path)
        ),
        delete_if_exists(TempPath)).

%!  append_history(+Path:atom, +Run:dict, +RunnerArgs:list,
%!                 +RungPairs:list) is det.
%
%   Append one compact JSONL history entry. RungPairs is the domain-owned list
%   of `Key-Cell` pairs; this predicate adds common host, Git, SWI, timestamp,
%   and tier provenance. JSON string whitespace is preserved verbatim.
append_history(Path, Run, RunnerArgs, RungPairs) :-
    must_be(list, RunnerArgs),
    must_be(list, RungPairs),
    dict_pairs(Rungs, rungs, RungPairs),
    ( valid_rungs(Rungs) -> true
    ; throw(error(bench_history_bad_rungs, _))
    ),
    get_dict(swi_prolog, Run, Swi),
    current_host(Host),
    git_commit(Commit),
    git_dirty(Dirty),
    get_time(Now),
    Epoch is round(Now),
    format_time(atom(Date), '%Y-%m-%dT%H:%M:%S', Now),
    ( memberchk('--heavy', RunnerArgs) -> Tiers = 'core+heavy' ; Tiers = core ),
    Entry = _{epoch:Epoch, date:Date, commit:Commit, dirty:Dirty,
              host:Host, swi:Swi, tiers:Tiers, rungs:Rungs},
    with_output_to(string(Line),
                   json_write_dict(current_output, Entry, [width(0)])),
    ( history_needs_newline(Path) -> Prefix = "\n" ; Prefix = "" ),
    setup_call_cleanup(
        open(Path, append, Stream, [encoding(utf8)]),
        format(Stream, '~s~s~n', [Prefix, Line]),
        close(Stream)),
    ( Dirty == true -> DirtyMark = ' +dirty' ; DirtyMark = '' ),
    format("history appended: ~w  (commit ~w~w, ~w)~n",
           [Path, Commit, DirtyMark, Tiers]).

%!  read_history(+Path:atom, -Entries:list(dict)) is det.
%
%   Read a JSONL ledger. A missing file is an empty ledger. Blank lines are
%   ignored; malformed or structurally invalid lines throw with a one-based line
%   number instead of disappearing from the trend.
read_history(Path, Entries) :-
    ( exists_file(Path)
    -> setup_call_cleanup(
           open(Path, read, Stream, [encoding(utf8)]),
           read_history_lines(Stream, Path, 1, [], Reversed),
           close(Stream)),
       reverse(Reversed, Entries)
    ; Entries = []
    ).

read_history_lines(Stream, Path, LineNo, Entries0, Entries) :-
    read_line_to_string(Stream, Line),
    ( Line == end_of_file
    -> Entries = Entries0
    ; blank_string(Line)
    -> Next is LineNo + 1,
       read_history_lines(Stream, Path, Next, Entries0, Entries)
    ; parse_history_line(Path, LineNo, Line, Entry),
      Next is LineNo + 1,
      read_history_lines(Stream, Path, Next, [Entry|Entries0], Entries)
    ).

parse_history_line(Path, LineNo, Line, Entry) :-
    catch(setup_call_cleanup(
              open_string(Line, Stream),
              read_json_stream(history(Path, LineNo), Stream, Parsed),
              close(Stream)),
          Error,
          throw(error(bench_history_line(Path, LineNo, Error), _))),
    ( valid_history_entry(Parsed)
    -> Entry = Parsed
    ; throw(error(bench_history_line(
                      Path, LineNo, error(bench_history_bad_envelope, _)), _))
    ).

valid_history_entry(Entry) :-
    get_dict(date, Entry, _),
    get_dict(commit, Entry, _),
    get_dict(dirty, Entry, _),
    get_dict(host, Entry, _),
    get_dict(swi, Entry, _),
    get_dict(tiers, Entry, _),
    get_dict(rungs, Entry, Rungs),
    valid_rungs(Rungs).

valid_rungs(Rungs) :-
    is_dict(Rungs),
    dict_pairs(Rungs, _, [_-First|Rest]),
    is_dict(First),
    forall(member(_-Cell, Rest), is_dict(Cell)).

%!  render_history(+Title, +Metrics:list, +Entries:list(dict)) is det.
%
%   Render common entry provenance and one trend table per
%   `metric(Label, CellKey)` in Metrics. Domain callers own the title, metric
%   selection, and empty-ledger guidance.
render_history(Title, Metrics, Entries) :-
    must_be(list, Metrics),
    must_be(list, Entries),
    length(Entries, Count),
    format("~w history  (~d run(s), most recent last)~n~n", [Title, Count]),
    forall(nth1(Index, Entries, Entry), print_entry_meta(Index, Entry)),
    all_rung_keys(Entries, Keys),
    forall(member(metric(Label, CellKey), Metrics),
           render_metric(Label, CellKey, Keys, Entries)).

print_entry_meta(Index, Entry) :-
    get_dict(date, Entry, Date),
    get_dict(commit, Entry, Commit),
    get_dict(tiers, Entry, Tiers),
    get_dict(host, Entry, Host),
    ( get_dict(dirty, Entry, true) -> DirtyMark = '+dirty' ; DirtyMark = '' ),
    format("  [~d] ~w  ~w~w  ~w  ~w~n",
           [Index, Date, Commit, DirtyMark, Tiers, Host]).

all_rung_keys(Entries, Keys) :-
    findall(Key,
            ( member(Entry, Entries),
              get_dict(rungs, Entry, Rungs),
              dict_pairs(Rungs, _, Pairs),
              member(Key-_, Pairs)
            ),
            Keys0),
    sort(Keys0, Keys).

render_metric(Label, CellKey, Keys, Entries) :-
    format("~n~w per rung  (vs prev = last step, vs first = cumulative):~n~n",
           [Label]),
    format("~w~t~30|~t~w~16+~t~w~12+~t~w~12+~t~w~7+~n",
           ['rung', 'latest', 'vs prev', 'vs first', 'runs']),
    forall(member(Key, Keys), print_rung_trend(Key, CellKey, Entries)).

print_rung_trend(Key, CellKey, Entries) :-
    rung_series(Key, CellKey, Entries, Series),
    ( Series == [] -> true
    ; last(Series, Latest),
      Series = [First|_],
      ( append(_, [Previous, Latest], Series) -> true ; Previous = Latest ),
      length(Series, Count),
      pct_string(Latest, Previous, PreviousDelta),
      pct_string(Latest, First, FirstDelta),
      format("~w~t~30|~t~D~16+~t~w~12+~t~w~12+~t~d~7+~n",
             [Key, Latest, PreviousDelta, FirstDelta, Count])
    ).

rung_series(Key, CellKey, Entries, Series) :-
    findall(Value,
            ( member(Entry, Entries),
              get_dict(rungs, Entry, Rungs),
              get_dict(Key, Rungs, Cell),
              get_dict(CellKey, Cell, Value)
            ),
            Series).

pct_string(_New, 0, 'n/a') :- !.
pct_string(New, Old, String) :-
    Percent is (New - Old) / Old * 100.0,
    ( Percent >= 0
    -> format(atom(Signed), '+~2f', [Percent])
    ; format(atom(Signed), '~2f', [Percent])
    ),
    atom_concat(Signed, '%', String).

%!  current_host(-Host:atom) is det.
%
%   Return `uname -nm` normalized to an atom, or `unknown` when unavailable.
%   Host is provenance only and never a benchmark gate.
current_host(Host) :-
    catch(capture_process(path(uname), ['-nm'], null,
                          Raw, _Stderr, Status), _, fail),
    Status == exit(0),
    normalize_space(atom(Candidate), Raw),
    Candidate \== '',
    !,
    Host = Candidate.
current_host(unknown).

git_commit(Commit) :-
    catch(capture_process(path(git), ['rev-parse', '--short', 'HEAD'], null,
                          Raw, _Stderr, Status), _, fail),
    Status == exit(0),
    normalize_space(atom(Candidate), Raw),
    Candidate \== '',
    !,
    Commit = Candidate.
git_commit(unknown).

git_dirty(Dirty) :-
    catch(capture_process(
              path(git), ['status', '--porcelain', '--untracked-files=no'], null,
              Raw, _Stderr, Status), _, fail),
    Status == exit(0),
    !,
    normalize_space(string(Text), Raw),
    ( Text == "" -> Dirty = false ; Dirty = true ).
git_dirty(false).

history_needs_newline(Path) :-
    exists_file(Path),
    size_file(Path, Size),
    Size > 0,
    setup_call_cleanup(
        open(Path, read, Stream, [type(binary)]),
        ( seek(Stream, -1, eof, _), get_byte(Stream, Last) ),
        close(Stream)),
    Last =\= 10.

blank_string(String) :-
    normalize_space(string(""), String).

text_to_atom(Text, Atom) :-
    text_to_string(Text, String),
    atom_string(Atom, String).

delete_if_exists(Path) :-
    ( exists_file(Path) -> delete_file(Path) ; true ).

:- multifile prolog:error_message//1.
prolog:error_message(bench_json_missing(Path)) -->
    [ 'benchmark JSON file does not exist: ~w'-[Path] ].
prolog:error_message(bench_json_invalid(Context, Error)) -->
    [ 'invalid benchmark JSON in ~w: ~p'-[Context, Error] ].
prolog:error_message(bench_json_not_object(Context)) -->
    [ 'benchmark JSON in ~w is not an object'-[Context] ].
prolog:error_message(bench_json_trailing_data(Context)) -->
    [ 'benchmark JSON in ~w has trailing non-whitespace data'-[Context] ].
prolog:error_message(bench_record_empty_results) -->
    [ 'benchmark recording requires at least one measured row' ].
prolog:error_message(bench_record_duplicate_key(Key)) -->
    [ 'benchmark recording contains duplicate row key ~w'-[Key] ].
prolog:error_message(bench_record_verification_failed(Path)) -->
    [ 'benchmark baseline read-back verification failed for ~w'-[Path] ].
prolog:error_message(bench_history_line(Path, LineNo, Error)) -->
    [ 'invalid benchmark history at ~w line ~d: ~p'-[Path, LineNo, Error] ].
prolog:error_message(bench_history_bad_envelope) -->
    [ 'benchmark history entry has invalid provenance or rungs' ].
prolog:error_message(bench_history_bad_rungs) -->
    [ 'benchmark history append requires nonempty dict-valued rungs' ].
