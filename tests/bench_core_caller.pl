:- module(bench_core_test_caller, []).

:- use_module('../benchmarks/bench_core.pl').

measure_local(Summary) :-
    measure(local_sampler, _{warmup:0, iterations:1}, Summary).

local_sampler(_{value:7}).

inproc_local(Sample) :-
    inproc_sampler(local_goal, Sample).

local_goal.
