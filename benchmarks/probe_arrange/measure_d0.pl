#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/3]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   directory_file_path(Dir, 'd0_support.pl', D0), use_module(D0).

:- initialization(main, main).

control(light_09x09_08w, light, 'fixtures/ladder_09x09_08w.pl', 9, 8).
control(light_15x15_12w, light, 'fixtures/ladder_15x15_12w.pl', 15, 12).
control(dense_15x15_32w, dense, 'fixtures/ladder_15x15_32w.pl', 15, 32).
control(dense_21x21_80w, dense, 'fixtures/ladder_21x21_80w.pl', 21, 80).

corner(topleft_across).
corner(topright).

operation_control(light_09x09_08w).
operation_control(dense_21x21_80w).

main :-
    catch(run, Error, fatal(Error)).

run :-
    findall(Row, measurement_row(Row), Rows),
    findall(Row, operation_row(Row), OperationRows),
    Result = _{probe:'P-D0',seed:null,budget:500000000,
               controls:Rows,operation_controls:OperationRows},
    json_write_dict(user_output, Result, [width(0)]), nl.

measurement_row(Row) :-
    control(Name, Class, Fixture, Grid, Count),
    corner(Corner),
    heartbeat(Name, Corner, warmup),
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange_d0:d0_corner(
        Words, Grid, Corner, null, true, _WarmResult, _WarmSummary,
        _WarmDecisions, _WarmTiming),
    heartbeat(Name, Corner, authority),
    probe_arrange:authority_corner(
        Words, Grid, Corner, 500000000, null, Authority, _AuthorityTiming),
    probe_arrange_d0:d0_corner(
        Words, Grid, Corner, null, false, Base, BaseSummary,
        BaseDecisions, _BaseTiming),
    heartbeat(Name, Corner, pass1),
    probe_arrange_d0:d0_corner(
        Words, Grid, Corner, null, true, Result1, Summary1, Decisions1, Timing1),
    heartbeat(Name, Corner, pass2),
    probe_arrange_d0:d0_corner(
        Words, Grid, Corner, null, true, Result2, Summary2, Decisions2, Timing2),
    equivalent_result(Authority, Base, Signature, Reward),
    equivalent_result(Authority, Result1, Signature, Reward),
    equivalent_result(Authority, Result2, Signature, Reward),
    BaseDecisions == Decisions1,
    Decisions1 == Decisions2,
    Summary1 =@= Summary2,
    Result1 == Result2,
    get_dict(inferences, Timing1, Inferences1),
    get_dict(inferences, Timing2, Inferences2),
    Inferences1 =:= Inferences2,
    length(Decisions1, DecisionCount),
    BaseSummary.decisions =:= DecisionCount,
    Summary1.decisions =:= DecisionCount,
    Row = _{control:Name,class:Class,fixture:Fixture,grid:Grid,
            words:Count,corner:Corner,outcome:placed,reward:Reward,
            layout_signature:Signature,decision_count:DecisionCount,
            duplicate_nonwall_exact:true,
            measured_inferences:Inferences1,
            duplicate_inferences:Inferences2,summary:Summary1},
    heartbeat(Name, Corner, done).

operation_row(Row) :-
    operation_control(Name),
    control(Name, Class, Fixture, Grid, Count),
    heartbeat(Name, operation, start),
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange:authority_operation(
        Words, Grid, 500000000, null, Authority, _AuthorityTiming),
    probe_arrange_d0:d0_operation(
        Words, Grid, null, false, Base, _BaseSummary, BaseDecisions, _BaseTiming),
    probe_arrange_d0:d0_operation(
        Words, Grid, null, true, Observed, Summary, ObservedDecisions, _Timing),
    equivalent_result(Authority, Base, Signature, Reward),
    equivalent_result(Authority, Observed, Signature, Reward),
    BaseDecisions == ObservedDecisions,
    length(ObservedDecisions, DecisionCount),
    Row = _{control:Name,class:Class,fixture:Fixture,grid:Grid,words:Count,
            corners:[topleft_across,topright],outcome:placed,reward:Reward,
            layout_signature:Signature,decision_count:DecisionCount,
            decision_order_exact:true,summary:Summary},
    heartbeat(Name, operation, done).

equivalent_result(result(placed,ok,false,_,NumberedA,Reward),
                  result(placed,ok,false,_,NumberedB,Reward),
                  Signature, Reward) :-
    probe_arrange:layout_signature(NumberedA, Signature),
    probe_arrange:layout_signature(NumberedB, Signature).

heartbeat(Name, Corner, Phase) :-
    format(user_error,
           'P-D0 heartbeat control=~w corner=~w phase=~w~n',
           [Name,Corner,Phase]),
    flush_output(user_error).

fatal(Error) :- print_message(error, Error), halt(2).
