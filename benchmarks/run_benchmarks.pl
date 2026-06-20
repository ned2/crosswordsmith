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
      [opt(format),     type(atom),    default(text), meta('FORMAT'),
       longflags([format]),
       help('output format: text, csv, or json')],
      [opt(strategy),   type(atom),    default(''), meta('STRAT'),
       longflags([strategy]),
       help('solver strategy (default: the production default_strategy/1)')]
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
    memberchk(strategy(StrategyOpt), Opts),
    resolve_strategy(StrategyOpt, Strategy),
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
               format: Format,
               strategy: Strategy}.

% '' means "use the production default"; otherwise validate against the
% strategies crossword.pl defines.
resolve_strategy('', Strategy) :-
    !,
    default_strategy(Strategy).
resolve_strategy(Strategy, Strategy) :-
    (   valid_strategy(Strategy)
    ->  true
    ;   throw(error(invalid_strategy(Strategy), _))
    ).


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
    memberchk(Format, [text, csv, json]),
    !.
validate_format(Format) :-
    throw(error(invalid_format(Format), _)).


run_benchmark(Config, Report) :-
    load_fixture(Config.fixture, Words),
    warmup(Config.strategy, Words, Config.grid, Config.startLoc, Config.warmup),
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


warmup(_Strategy, _Words, _GridLen, _StartLoc, 0) :-
    !.
warmup(Strategy, Words, GridLen, StartLoc, Warmup) :-
    forall(between(1, Warmup, _), must_solve(Strategy, Words, GridLen, StartLoc)).


measure(Words, Config, Result) :-
    findall(Time,
            ( between(1, Config.iterations, _),
              timed_solve(Config.strategy, Words, Config.grid, Config.startLoc, Time)
            ),
            Times),
    summarize_times(Times, Summary),
    file_base_name(Config.fixture, FixtureName),
    Result = _{fixture: Config.fixture,
               fixtureName: FixtureName,
               gridLength: Config.grid,
               startLoc: Config.startLoc,
               strategy: Config.strategy,
               iterations: Config.iterations,
               wall: Summary.wall,
               cpu: Summary.cpu,
               inferences: Summary.inferences}.


timed_solve(Strategy, Words, GridLen, StartLoc, Time) :-
    call_time(must_solve(Strategy, Words, GridLen, StartLoc), Time).


must_solve(Strategy, Words, GridLen, StartLoc) :-
    find_crossword(Strategy, GridLen, Words, StartLoc, _Grid, _PlacedWords),
    !.
must_solve(_Strategy, _Words, GridLen, StartLoc) :-
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


output_report(text, Report) :-
    !,
    emit_text(Report).
output_report(csv, Report) :-
    !,
    emit_csv(Report).
output_report(json, Report) :-
    !,
    emit_json_report(Report).


emit_text(Report) :-
    get_dict(metadata, Report, Metadata),
    get_dict(results, Report, [Result]),
    emit_text_atom(tool, Metadata.tool),
    emit_text_atom(fixture, Result.fixture),
    emit_text_atom(fixture_name, Result.fixtureName),
    emit_text_integer(grid_length, Result.gridLength),
    emit_text_atom(start_loc, Result.startLoc),
    emit_text_atom(strategy, Result.strategy),
    emit_text_integer(iterations, Result.iterations),
    emit_text_integer(warmup, Metadata.warmup),
    emit_text_atom(swi_prolog, Metadata.swiProlog),
    emit_text_float(epoch_seconds, Metadata.epochSeconds),
    emit_text_float(wall_min_s, Result.wall.min),
    emit_text_float(wall_median_s, Result.wall.median),
    emit_text_float(wall_mean_s, Result.wall.mean),
    emit_text_float(cpu_min_s, Result.cpu.min),
    emit_text_float(cpu_median_s, Result.cpu.median),
    emit_text_float(cpu_mean_s, Result.cpu.mean),
    emit_text_count(inferences_min, Result.inferences.min),
    emit_text_count(inferences_median, Result.inferences.median),
    emit_text_count(inferences_mean, Result.inferences.mean).


emit_text_atom(Key, Value) :-
    emit_text_key(Key),
    format("~w~n", [Value]).


emit_text_integer(Key, Value) :-
    emit_text_key(Key),
    format("~d~n", [Value]).


emit_text_float(Key, Value) :-
    emit_text_key(Key),
    format("~6f~n", [Value]).


emit_text_count(Key, Value) :-
    emit_text_key(Key),
    format("~0f~n", [Value]).


emit_text_key(Key) :-
    format("~w:~t~22|", [Key]).


emit_csv(Report) :-
    get_dict(results, Report, Results),
    format("fixture,grid,start,strategy,iterations,wall_min_s,wall_median_s,wall_mean_s,cpu_min_s,cpu_median_s,cpu_mean_s,inferences_min,inferences_median,inferences_mean~n",
           []),
    forall(member(Result, Results), emit_csv_row(Result)).


emit_csv_row(Result) :-
    format("~w,~d,~w,~w,~d,~6f,~6f,~6f,~6f,~6f,~6f,~0f,~0f,~0f~n",
           [ Result.fixtureName,
             Result.gridLength,
             Result.startLoc,
             Result.strategy,
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
    [ 'benchmark format must be text, csv, or json, got ~q'-[Format] ].
prolog:error_message(invalid_strategy(Strategy)) -->
    [ 'benchmark strategy is not a known solver strategy: ~q'-[Strategy] ].
prolog:error_message(fixture_unsolved(GridLen, StartLoc)) -->
    [ 'benchmark fixture did not solve for grid length ~d and start location ~q'-
      [GridLen, StartLoc] ].
