#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/2]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   directory_file_path(Dir, 'ad2_delta.pl', AD2), use_module(AD2),
   directory_file_path(Dir, 'closeout_direct.pl', Direct), use_module(Direct).

:- initialization(main, main).

control('15x15/34w', 'fixtures/ladder_15x15_34w.pl', 15, 34).
control('15x15/36w', 'fixtures/ladder_15x15_36w.pl', 15, 36).
control('21x21/80w', 'fixtures/ladder_21x21_80w.pl', 21, 80).

corner(topleft_across).
corner(topright).

main :- catch(run, Error, (print_message(error, Error), halt(2))).

run :-
    findall(Row, row(Row), Rows),
    json_write_dict(user_output,
                    _{rig:current_direct_attribution,rows:Rows}, [width(0)]),
    nl.

row(Row) :-
    control(Name, Fixture, Grid, Count),
    corner(Corner),
    format(user_error,
           'closeout-direct heartbeat control=~w corner=~w phase=start~n',
           [Name,Corner]),
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange_closeout_direct:product_corner(
        Words, Grid, Corner, Product),
    probe_arrange_closeout_direct:closeout_corner(
        Words, Grid, Corner, none, First, Stats, Trace),
    probe_arrange_closeout_direct:closeout_corner(
        Words, Grid, Corner, none, Second, Stats2, Trace2),
    probe_arrange_ad2:ad2_corner(
        Words, Grid, Corner, null, Differential, _Summary, _Counts, AD2Trace),
    require_identical(product_result, Product, First),
    require_identical(duplicate_result, First, Second),
    require_identical(differential_result, First, Differential),
    require_identical(duplicate_stats, Stats, Stats2),
    require_identical(duplicate_trace, Trace, Trace2),
    require_identical(differential_trace, Trace, AD2Trace),
    First = result(Numbered, Reward),
    Trace = trace(Selections, Decisions),
    length(Selections, SelectionCount),
    length(Decisions, DecisionCount),
    DecisionCount =:= Stats.legal_decisions,
    probe_arrange:layout_signature(Numbered, LayoutSignature),
    probe_arrange_closeout_direct:trace_digest(LayoutSignature, LayoutDigest),
    probe_arrange_closeout_direct:trace_digest(Trace, TraceDigest),
    Row = Stats.put(_{control:Name,fixture:Fixture,grid:Grid,words:Count,
                      corner:Corner,outcome:placed,reward:Reward,
                      layout_signature:LayoutSignature,layout_sha1:LayoutDigest,
                      selections:SelectionCount,trace_sha1:TraceDigest,
                      product_result_exact:true,ad2_trace_exact:true,
                      duplicate_counters_exact:true}),
    format(user_error,
           'closeout-direct heartbeat control=~w corner=~w phase=done~n',
           [Name,Corner]).

require_identical(_, A, B) :- A == B, !.
require_identical(Label, A, B) :-
    variant_sha1(A, HashA),
    variant_sha1(B, HashB),
    throw(error(closeout_direct_mismatch(Label,HashA,HashB), _)).
