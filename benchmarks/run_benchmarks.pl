#!/usr/bin/env swipl

:- set_prolog_flag(verbose, silent).

:- use_module(library(http/json)).
:- use_module(library(lists)).
:- use_module(library(optparse)).
:- use_module(library(statistics)).

:- dynamic repo_root/1.

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', RepoRoot,
                      [ relative_to(BenchDir),
                        file_type(directory),
                        access(read)
                      ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'crossword.pl', Crossword),
   consult(Crossword).

:- initialization(benchmark_main, main).


benchmark_main :-
    current_prolog_flag(argv, Argv),
    benchmark_opts_spec(Spec),
    opt_parse(Spec, Argv, Opts, Positional),
    validate_options(Opts, Positional, Config),
    run_benchmark(Config, Report),
    output_report(Config.format, Report).


benchmark_opts_spec(
    [ [opt(help),       type(boolean), default(false),
       shortflags([h]), longflags([help]),
       help('show this help and exit')],
      [opt(fixture),    type(atom),    default(''), meta('FILE'),
       longflags([fixture]),
       help('Prolog fixture file defining clues/1')],
      [opt(grid),       type(integer), default(17), meta('N'),
       longflags([grid]),
       help('grid side length')],
      [opt(start_loc),  type(atom),    default(topleft_across), meta('LOC'),
       longflags(['start-loc']),
       help('first-word start location')],
      [opt(iterations), type(integer), default(30), meta('N'),
       longflags([iterations]),
       help('measured repetitions')],
      [opt(warmup),     type(integer), default(3), meta('N'),
       longflags([warmup]),
       help('unreported warmup repetitions')],
      [opt(format),     type(atom),    default(both), meta('FORMAT'),
       longflags([format]),
       help('output format: table, json, or both')]
    ]).


validate_options(Opts, _Positional, _Config) :-
    memberchk(help(true), Opts),
    !,
    print_benchmark_usage,
    halt(0).
validate_options(Opts, Positional, Config) :-
    memberchk(fixture(FixtureOpt), Opts),
    fixture_arg(FixtureOpt, Positional, FixtureArg),
    resolve_fixture_path(FixtureArg, FixtureFile),
    memberchk(grid(GridLen), Opts),
    memberchk(start_loc(StartLoc), Opts),
    memberchk(iterations(Iterations), Opts),
    memberchk(warmup(Warmup), Opts),
    memberchk(format(Format), Opts),
    validate_grid(GridLen),
    validate_start_loc(StartLoc),
    validate_iterations(Iterations),
    validate_warmup(Warmup),
    validate_format(Format),
    Config = _{fixture: FixtureFile,
               grid: GridLen,
               startLoc: StartLoc,
               iterations: Iterations,
               warmup: Warmup,
               format: Format}.


print_benchmark_usage :-
    benchmark_opts_spec(Spec),
    opt_help(Spec, Help),
    format("Usage: run_benchmarks.pl [options] [fixture.pl]~n~n~w", [Help]).


fixture_arg('', [], Fixture) :-
    !,
    default_fixture(Fixture).
fixture_arg('', [Fixture], Fixture) :-
    !.
fixture_arg(Fixture, [], Fixture) :-
    Fixture \== '',
    !.
fixture_arg(Fixture, Positional, _Selected) :-
    Fixture \== '',
    Positional \== [],
    !,
    throw(error(duplicate_fixture_argument(Fixture, Positional), _)).
fixture_arg(_Fixture, Positional, _Selected) :-
    throw(error(unexpected_positional_args(Positional), _)).


default_fixture('fixtures/bundled_17_clues.pl').


resolve_fixture_path(FixtureArg, FixtureFile) :-
    absolute_file_name(FixtureArg, Candidate, [file_errors(fail)]),
    exists_file(Candidate),
    !,
    FixtureFile = Candidate.
resolve_fixture_path(FixtureArg, FixtureFile) :-
    repo_file(FixtureArg, Candidate),
    exists_file(Candidate),
    !,
    FixtureFile = Candidate.
resolve_fixture_path(FixtureArg, _FixtureFile) :-
    throw(error(fixture_not_found(FixtureArg), _)).


validate_grid(GridLen) :-
    integer(GridLen),
    GridLen > 0,
    !.
validate_grid(GridLen) :-
    throw(error(invalid_grid(GridLen), _)).


validate_start_loc(StartLoc) :-
    valid_loc(StartLoc),
    !.
validate_start_loc(StartLoc) :-
    throw(error(invalid_start_loc(StartLoc), _)).


validate_iterations(Iterations) :-
    integer(Iterations),
    Iterations > 0,
    !.
validate_iterations(Iterations) :-
    throw(error(invalid_iterations(Iterations), _)).


validate_warmup(Warmup) :-
    integer(Warmup),
    Warmup >= 0,
    !.
validate_warmup(Warmup) :-
    throw(error(invalid_warmup(Warmup), _)).


validate_format(Format) :-
    memberchk(Format, [table, json, both]),
    !.
validate_format(Format) :-
    throw(error(invalid_format(Format), _)).


run_benchmark(Config, Report) :-
    load_fixture(Config.fixture, Words),
    warmup(Words, Config.grid, Config.startLoc, Config.warmup),
    measure(Words, Config, Result),
    report(Config, Result, Report).


load_fixture(FixtureFile, Words) :-
    setup_call_cleanup(open(FixtureFile, read, Stream),
                       read_fixture_clues(Stream, FixtureFile, Words),
                       close(Stream)).


read_fixture_clues(Stream, FixtureFile, Words) :-
    read_term(Stream, Term, []),
    (   Term == end_of_file
    ->  throw(error(fixture_missing_clues(FixtureFile), _))
    ;   Term = clues(Words)
    ->  true
    ;   read_fixture_clues(Stream, FixtureFile, Words)
    ).


warmup(_Words, _GridLen, _StartLoc, 0) :-
    !.
warmup(Words, GridLen, StartLoc, Warmup) :-
    forall(between(1, Warmup, _), must_solve(Words, GridLen, StartLoc)).


measure(Words, Config, Result) :-
    findall(Time,
            ( between(1, Config.iterations, _),
              timed_solve(Words, Config.grid, Config.startLoc, Time)
            ),
            Times),
    summarize_times(Times, Summary),
    file_base_name(Config.fixture, FixtureName),
    Result = _{fixture: Config.fixture,
               fixtureName: FixtureName,
               gridLength: Config.grid,
               startLoc: Config.startLoc,
               iterations: Config.iterations,
               wall: Summary.wall,
               cpu: Summary.cpu,
               inferences: Summary.inferences}.


timed_solve(Words, GridLen, StartLoc, Time) :-
    call_time(must_solve(Words, GridLen, StartLoc), Time).


must_solve(Words, GridLen, StartLoc) :-
    find_crossword(GridLen, Words, StartLoc, _Grid, _PlacedWords),
    !.
must_solve(_Words, GridLen, StartLoc) :-
    throw(error(fixture_unsolved(GridLen, StartLoc), _)).


repo_file(Relative, File) :-
    repo_root(Root),
    directory_file_path(Root, Relative, File).


summarize_times(Times, _{wall: Wall, cpu: Cpu, inferences: Inferences}) :-
    findall(W, (member(Time, Times), get_dict(wall, Time, W)), Walls),
    findall(C, (member(Time, Times), get_dict(cpu, Time, C)), Cpus),
    findall(I, (member(Time, Times), get_dict(inferences, Time, I)), InferenceCounts),
    stats(Walls, Wall),
    stats(Cpus, Cpu),
    stats(InferenceCounts, Inferences).


stats(Values, _{min: Min, median: Median, mean: Mean}) :-
    Values = [_|_],
    msort(Values, Sorted),
    Sorted = [Min|_],
    length(Sorted, Len),
    sum_list(Sorted, Sum),
    Mean is Sum / Len,
    median(Sorted, Len, Median).


median(Sorted, Len, Median) :-
    1 is Len mod 2,
    !,
    Index is Len // 2,
    nth0(Index, Sorted, Median).
median(Sorted, Len, Median) :-
    HiIndex is Len // 2,
    LoIndex is HiIndex - 1,
    nth0(LoIndex, Sorted, Lo),
    nth0(HiIndex, Sorted, Hi),
    Median is (Lo + Hi) / 2.


report(Config, Result, Report) :-
    current_prolog_flag(version_data, VersionData),
    version_atom(VersionData, Version),
    get_time(EpochSeconds),
    Report = _{metadata: _{tool: 'crosswordsmith-benchmarks',
                           fixture: Config.fixture,
                           gridLength: Config.grid,
                           startLoc: Config.startLoc,
                           iterations: Config.iterations,
                           warmup: Config.warmup,
                           epochSeconds: EpochSeconds,
                           swiProlog: Version},
               results: [Result]}.


version_atom(swi(Major, Minor, Patch, _Extra), Version) :-
    !,
    format(atom(Version), '~d.~d.~d', [Major, Minor, Patch]).
version_atom(VersionData, Version) :-
    term_to_atom(VersionData, Version).


output_report(table, Report) :-
    !,
    emit_table(Report).
output_report(json, Report) :-
    !,
    emit_json_report(Report).
output_report(both, Report) :-
    emit_table(Report),
    nl,
    emit_json_report(Report).


emit_table(Report) :-
    get_dict(metadata, Report, Metadata),
    format("crosswordsmith benchmarks~n", []),
    format("fixture: ~w, grid: ~d, start: ~w, iterations: ~d, warmup: ~d, swipl: ~w~n~n",
           [Metadata.fixture, Metadata.gridLength, Metadata.startLoc,
            Metadata.iterations, Metadata.warmup, Metadata.swiProlog]),
    format("fixture\tgrid\tstart\titerations\twall_min_s\twall_median_s\twall_mean_s\tcpu_min_s\tcpu_median_s\tcpu_mean_s\tinferences_min\tinferences_median\tinferences_mean~n",
           []),
    get_dict(results, Report, Results),
    forall(member(Result, Results), emit_table_row(Result)).


emit_table_row(Result) :-
    format("~w\t~d\t~w\t~d\t~6f\t~6f\t~6f\t~6f\t~6f\t~6f\t~0f\t~0f\t~0f~n",
           [ Result.fixtureName,
             Result.gridLength,
             Result.startLoc,
             Result.iterations,
             Result.wall.min,
             Result.wall.median,
             Result.wall.mean,
             Result.cpu.min,
             Result.cpu.median,
             Result.cpu.mean,
             Result.inferences.min,
             Result.inferences.median,
             Result.inferences.mean
           ]).


emit_json_report(Report) :-
    current_output(Out),
    json_write_dict(Out, Report),
    nl(Out).


:- multifile prolog:error_message//1.
prolog:error_message(duplicate_fixture_argument(Fixture, Positional)) -->
    [ 'benchmark fixture was given twice: --fixture ~q and positional ~q'-
      [Fixture, Positional] ].
prolog:error_message(unexpected_positional_args(Positional)) -->
    [ 'benchmark command accepts at most one fixture path, got ~q'-[Positional] ].
prolog:error_message(fixture_not_found(Fixture)) -->
    [ 'benchmark fixture not found: ~q'-[Fixture] ].
prolog:error_message(fixture_missing_clues(Fixture)) -->
    [ 'benchmark fixture ~q does not define clues/1'-[Fixture] ].
prolog:error_message(invalid_grid(GridLen)) -->
    [ 'benchmark grid length must be a positive integer, got ~q'-[GridLen] ].
prolog:error_message(invalid_start_loc(StartLoc)) -->
    [ 'benchmark start location is not valid: ~q'-[StartLoc] ].
prolog:error_message(invalid_iterations(Iterations)) -->
    [ 'benchmark iterations must be a positive integer, got ~q'-[Iterations] ].
prolog:error_message(invalid_warmup(Warmup)) -->
    [ 'benchmark warmup must be a non-negative integer, got ~q'-[Warmup] ].
prolog:error_message(invalid_format(Format)) -->
    [ 'benchmark format must be one of table, json, or both, got ~q'-[Format] ].
prolog:error_message(fixture_unsolved(GridLen, StartLoc)) -->
    [ 'benchmark fixture did not solve for grid length ~d and start location ~q'-
      [GridLen, StartLoc] ].
