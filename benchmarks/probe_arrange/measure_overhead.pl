#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/2]).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir),file_type(directory)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe).
:- initialization(main, main).

% Serialized alternating protocol: two discarded warmup pairs followed by 21
% authority/lean pairs in one process. Median wall is the adjudication ratio;
% inference overhead is reported independently.
main :-
    current_prolog_flag(argv, [Fixture,G0,N0,Seed0,Corner|Rest]),
    atom_number(G0, Grid), atom_number(N0, Count),
    ( Seed0 == none -> Seed=null ; atom_number(Seed0, Seed) ),
    overhead_args(Rest, Mode, Trials),
    probe_arrange:load_fixture(Fixture, Count, Words),
    forall(between(1,2,_), pair(Mode,Words,Grid,Seed,Corner,_,_)),
    findall(A-I,
            ( between(1,Trials,T),
              format(user_error,'probe-arrange heartbeat overhead trial=~d/~d~n',[T,Trials]),
              pair(Mode,Words,Grid,Seed,Corner,A,I) ),
            Pairs),
    findall(W, member(sample(W,_)-_,Pairs), AWalls),
    findall(N, member(sample(_,N)-_,Pairs), AInfs),
    findall(W, member(_-sample(W,_),Pairs), IWalls),
    findall(N, member(_-sample(_,N),Pairs), IInfs),
    median(AWalls, AMW), median(IWalls, IMW),
    median(AInfs, AMI), median(IInfs, IMI),
    paired_ratios(AWalls, IWalls, Ratios),
    percentile(Ratios, 0.10, RatioP10),
    median(Ratios, RatioMedian),
    percentile(Ratios, 0.90, RatioP90),
    WallRatio is IMW / AMW, InfRatio is IMI / AMI,
    WallPct is (WallRatio-1.0)*100.0,
    InfPct is (InfRatio-1.0)*100.0,
    Row = _{fixture:Fixture,corner:Corner,search_seed:Seed,mode:Mode,
            warmup_pairs:2,trials:Trials,authority_wall_median:AMW,
            instrumented_wall_median:IMW,wall_ratio:WallRatio,
            wall_overhead_pct:WallPct,
            paired_wall_ratio_p10:RatioP10,
            paired_wall_ratio_median:RatioMedian,
            paired_wall_ratio_p90:RatioP90,
            authority_inferences_median:AMI,instrumented_inferences_median:IMI,
            inference_ratio:InfRatio,inference_overhead_pct:InfPct,
            authority_walls:AWalls,instrumented_walls:IWalls},
    json_write_dict(user_output, Row, [width(0)]), nl.

overhead_args([], lean, 21).
overhead_args([Trials0], lean, Trials) :- atom_number(Trials0, Trials).
overhead_args([Mode,Trials0|_], Mode, Trials) :-
    memberchk(Mode, [lean,full]), atom_number(Trials0, Trials).

pair(Mode,Words,Grid,Seed,Corner,sample(AW,AI),sample(IW,II)) :-
    probe_arrange:authority_corner(Words,Grid,Corner,500000000,Seed,
                                   result(placed,ok,false,_,_,_),AT),
    probe_arrange:twin_corner(Words,Grid,Corner,Seed,Mode,none,null,
                              result(placed,ok,false,_,_,_),run(_,IT)),
    AW=AT.wall, AI=AT.inferences, IW=IT.wall, II=IT.inferences.

median(Values, Median) :-
    msort(Values, Sorted), length(Sorted,N), I is N//2, nth0(I,Sorted,Median).

paired_ratios([], [], []).
paired_ratios([A|As], [I|Is], [R|Rs]) :-
    R is I/A, paired_ratios(As, Is, Rs).

percentile(Values, Fraction, Value) :-
    msort(Values, Sorted), length(Sorted,N),
    I is round(Fraction*(N-1)), nth0(I,Sorted,Value).
