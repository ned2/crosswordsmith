#!/usr/bin/env swipl
% Serialized A-G2 wall comparison: two discarded warmup pairs, then 21 pairs.

:- set_prolog_flag(verbose, silent).
:- use_module(library(filesex), [directory_file_path/3]).
:- use_module(library(json), [json_write_dict/3]).
:- use_module(library(process), [process_create/3, process_wait/2]).

:- prolog_load_context(file, Script), assertz(script_path(Script)).
:- dynamic script_path/1.
:- initialization(main, main).

main :-
    current_prolog_flag(argv, Args),
    (   Args = ['--sample',Root,Id]
    ->  sample(Root, Id)
    ;   Args = [Baseline,Candidate,Id]
    ->  measure(Baseline, Candidate, Id)
    ;   format(user_error,
               'usage: run_arrange_g2_wall.pl BASELINE CANDIDATE RUNG~n', []),
        halt(2)
    ).

sample(Root, Id) :-
    directory_file_path(Root, 'load.pl', Load), consult(Load),
    directory_file_path(Root, 'benchmarks/greedy_subjects.pl', Subjects),
    use_module(Subjects),
    directory_file_path(Root, 'benchmarks/greedy_workloads.pl', Workloads),
    consult(Workloads),
    greedy_workload(Id, Fixture, Size, _Framing, Command,
                    _Seed, _Corner, _Expected, _Iterations, _Warmup, _Tier),
    directory_file_path(Root, Fixture, File),
    greedy_subjects:load_words(File, Words),
    greedy_subjects:sweep_operation(Words, Size, Command, _),
    greedy_subjects:sweep_sampler(Words, Size, Command, Sample),
    format('~16g ~d~n', [Sample.wall, Sample.inferences]).

measure(Baseline, Candidate, Id) :-
    forall(between(1, 2, Pair),
           ( heartbeat(warmup, Pair, 2, baseline, Id), run_sample(Baseline, Id, _),
             heartbeat(warmup, Pair, 2, candidate, Id), run_sample(Candidate, Id, _) )),
    findall(B-C,
            ( between(1, 21, Pair),
              heartbeat(measured, Pair, 21, baseline, Id), run_sample(Baseline, Id, B),
              heartbeat(measured, Pair, 21, candidate, Id), run_sample(Candidate, Id, C) ),
            Pairs),
    summarize(Id, Baseline, Candidate, Pairs, Row),
    json_write_dict(current_output, Row, [width(0)]), nl.

heartbeat(Kind, Pair, Total, Side, Id) :-
    format(user_error, 'heartbeat: A-G2 wall ~w ~w pair ~d/~d ~w~n',
           [Id, Kind, Pair, Total, Side]).

run_sample(Root, Id, sample(Wall, Inferences)) :-
    script_path(Script),
    process_create(path(swipl), ['-q',Script,'--','--sample',Root,Id],
                   [stdout(pipe(Out)),stderr(std),process(PID)]),
    read_string(Out, _, Text), close(Out),
    process_wait(PID, exit(0)),
    split_string(Text, ' ', ' \n', [WallText,InfText]),
    number_string(Wall, WallText), number_string(Inferences, InfText).

summarize(Id, Baseline, Candidate, Pairs,
          _{rung:Id,baseline:Baseline,candidate:Candidate,warmup_pairs:2,
            measured_pairs:21,baseline_wall_median_s:BaselineMedian,
            candidate_wall_median_s:CandidateMedian,median_ratio:MedianRatio,
            paired_ratio_p10:P10,paired_ratio_p90:P90,
            baseline_inferences:BaselineInferences,
            candidate_inferences:CandidateInferences}) :-
    findall(W, member(sample(W,_)-_, Pairs), BaselineWalls),
    findall(W, member(_-sample(W,_), Pairs), CandidateWalls),
    findall(I, member(sample(_,I)-_, Pairs), BaselineInferences),
    findall(I, member(_-sample(_,I), Pairs), CandidateInferences),
    paired_ratios(BaselineWalls, CandidateWalls, Ratios),
    median(BaselineWalls, BaselineMedian),
    median(CandidateWalls, CandidateMedian),
    median(Ratios, MedianRatio),
    percentile(Ratios, 0.10, P10),
    percentile(Ratios, 0.90, P90).

paired_ratios([], [], []).
paired_ratios([Baseline|Bs], [Candidate|Cs], [Ratio|Ratios]) :-
    Ratio is Candidate / Baseline,
    paired_ratios(Bs, Cs, Ratios).

median(Values, Median) :-
    msort(Values, Sorted), length(Sorted, N), Index is N // 2,
    nth0(Index, Sorted, Median).

percentile(Values, Fraction, Value) :-
    msort(Values, Sorted), length(Sorted, N),
    Index is round(Fraction * (N - 1)), nth0(Index, Sorted, Value).
