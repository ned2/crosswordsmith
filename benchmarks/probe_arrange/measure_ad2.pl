#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/3]).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   directory_file_path(Dir, 'd0_support.pl', D0), use_module(D0),
   directory_file_path(Dir, 'ad1_buckets.pl', AD1), use_module(AD1),
   directory_file_path(Dir, 'ad2_delta.pl', AD2), use_module(AD2).

:- initialization(main, main).

control(light_09x09_08w, light, 'fixtures/ladder_09x09_08w.pl', 9, 8).
control(light_15x15_12w, light, 'fixtures/ladder_15x15_12w.pl', 15, 12).
control(dense_15x15_32w, dense, 'fixtures/ladder_15x15_32w.pl', 15, 32).
control(dense_21x21_80w, dense, 'fixtures/ladder_21x21_80w.pl', 21, 80).
corner(topleft_across).
corner(topright).

main :- catch(run, Error, (print_message(error, Error), halt(2))).

run :-
    findall(Row, row(Row), Rows),
    json_write_dict(user_output, _{experiment:'A-D2',controls:Rows},
                    [width(0)]), nl.

row(Row) :-
    control(Name, Class, Fixture, Grid, Count), corner(Corner),
    format(user_error, 'A-D2 control=~w corner=~w START~n', [Name,Corner]),
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange_ad1:ad1_corner(
        direct, Words, Grid, Corner, null,
        Full, FullCounts, FullDecisions),
    probe_arrange_ad2:ad2_corner(
        Words, Grid, Corner, null,
        Delta, Summary, trace(DeltaCounts, _RefreshTrace), DeltaDecisions),
    Full == Delta,
    FullCounts == DeltaCounts,
    FullDecisions == DeltaDecisions,
    Full = result(_Placed, Reward),
    Row = _{control:Name,class:Class,fixture:Fixture,grid:Grid,words:Count,
            corner:Corner,reward:Reward,equivalent:true,summary:Summary},
    format(user_error, 'A-D2 control=~w corner=~w DONE~n', [Name,Corner]).
