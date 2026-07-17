#!/usr/bin/env swipl
% benchmarks/run_arrange_greedy.pl - four-layer greedy arrange benchmark.

:- set_prolog_flag(verbose, silent).
:- use_module(library(json)).
:- use_module(library(optparse)).
:- use_module(library(filesex), [directory_file_path/3]).

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', Root,
                      [relative_to(BenchDir), file_type(directory), access(read)]),
   assertz(repo_root(Root)),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(BenchDir, 'bench_core.pl', Core), use_module(Core),
   directory_file_path(BenchDir, 'bench_cli.pl', BenchCli), use_module(BenchCli),
   directory_file_path(BenchDir, 'greedy_subjects.pl', Subjects), use_module(Subjects),
   directory_file_path(BenchDir, 'greedy_workloads.pl', Workloads), consult(Workloads).

:- dynamic repo_root/1.
:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    catch(opt_parse(Spec, Argv, Opts, Pos), E, (print_message(error, E), halt(2))),
    memberchk(identity(Identity), Opts),
    memberchk(format(Format), Opts),
    memberchk(fixture(Filter), Opts),
    memberchk(heavy(Heavy), Opts),
    memberchk(iterations(ItOverride), Opts),
    memberchk(warmup(WOverride), Opts),
    catch(validate_runner_options(greedy, Pos, Format, ItOverride, WOverride),
          E, (print_message(error, E), halt(2))),
    ( memberchk(help(true), Opts) -> usage, halt(0) ; true ),
    selected_specs(Filter, Heavy, Specs),
    catch(require_selected(greedy, Filter, Specs),
          E, (print_message(error, E), halt(2))),
    ( Identity == true
    -> identity_document(Specs, Doc), json_write_dict(current_output, Doc, [width(0)]), nl
    ;  maplist(run_spec(ItOverride, WOverride), Specs, Rows),
       emit(Format, Rows) ).

opts_spec([
    [opt(help), type(boolean), default(false), shortflags([h]), longflags([help])],
    [opt(format), type(atom), default(text), longflags([format])],
    [opt(fixture), type(atom), default(''), longflags([fixture])],
    [opt(heavy), type(boolean), default(false), longflags([heavy])],
    [opt(iterations), type(integer), default(-1), longflags([iterations])],
    [opt(warmup), type(integer), default(-1), longflags([warmup])],
    [opt(identity), type(boolean), default(false), longflags([identity])]
]).

usage :-
    opts_spec(Spec), opt_help(Spec, Help),
    format("Usage: swipl -q benchmarks/run_arrange_greedy.pl -- [options]~n~n~w", [Help]).

selected_specs(Filter, Heavy, Specs) :-
    findall(spec(Id,F,S,Fr,M,Seed,C,E,I,W,T),
            ( greedy_workload(Id,F,S,Fr,M,Seed,C,E,I,W,T),
              selected(Filter, Heavy, Id, T) ),
            Specs).

selected(Filter, _Heavy, Id, _Tier) :-
    Filter \== '', !, sub_atom(Id, _, _, _, Filter).
selected('', _Heavy, _Id, core) :- !.
selected('', true, _Id, heavy).

override(-1, Default, Default) :- !.
override(Value, _, Value).

run_spec(ItOverride, WOverride,
         spec(Id,Fixture,Size,Framing,Command,Seed,Corner,Expected,I0,W0,Tier), Row) :-
    override(ItOverride, I0, Iterations), override(WOverride, W0, Warmup),
    repo_path(Fixture, File), greedy_subjects:load_words(File, Words),
    length(Words, Total), command_k(Command, K),
    Opts = _{iterations:Iterations,warmup:Warmup},
    layer_measure(Id, construction,
                  greedy_subjects:construction_sampler(
                      Words,Size,Seed,Corner,Expected), Opts, Construction),
    layer_measure(Id, sweep,
                  greedy_subjects:sweep_sampler(Words,Size,Command), Opts, Sweep),
    format(user_error, "heartbeat: ~w postprocess setup (outside timing)~n", [Id]),
    greedy_subjects:build_raw_pool(Words, Size, Command, Raw),
    layer_measure(Id, postprocess,
                  greedy_subjects:postprocess_sampler(
                      Raw,Size,Total,Command,K), Opts, Post),
    repo_path(crosswordsmith, Exe),
    layer_measure(Id, command,
                  greedy_subjects:command_sampler(
                      Exe,File,Size,Framing,Command,0), Opts, Cmd),
    format(user_error, "heartbeat: ~w semantic replay/counters~n", [Id]),
    greedy_subjects:semantic_counters(Words, Size, Command, Counters),
    expected_json(Expected, ExpectedJson), command_json(Command, CommandJson),
    metrics_fields(Construction, construction, CF),
    metrics_fields(Sweep, sweep, SF),
    metrics_fields(Post, postprocess, PF),
    command_fields(Cmd, CmdF),
    put_dict(CF, _{}, D1), put_dict(SF, D1, D2), put_dict(PF, D2, D3),
    put_dict(CmdF, D3, Metrics),
    Row = _{rung:Id,fixture:Fixture,size:Size,framing:Framing,mode:CommandJson,
            construction:_{seed_answer:Seed,corner:Corner,expected:ExpectedJson},
            words:Total,iterations:Iterations,warmup:Warmup,tier:Tier,
            metrics:Metrics,semantic_counters:Counters}.

% One explicit first sample characterizes cold-vs-warm behavior. bench_core then
% performs the manifest warmup and measured repetitions. Heartbeats are outside
% each sampler's own call_time/process timing.
layer_measure(Id, Layer, Sampler, Opts, Result) :-
    format(user_error, "heartbeat: ~w ~w cold sample~n", [Id, Layer]),
    call(Sampler, Cold),
    flag(greedy_trial, _, 0),
    measure(heartbeat_sampler(Id, Layer, Sampler), Opts, Summary),
    Result = result(Cold, Summary).

heartbeat_sampler(Id, Layer, Sampler, Sample) :-
    flag(greedy_trial, N, N+1),
    format(user_error, "heartbeat: ~w ~w trial ~d~n", [Id, Layer, N]),
    call(Sampler, Sample).

metrics_fields(result(Cold, Summary), Prefix, Dict) :-
    Inf = Summary.stats.inferences,
    Wall = Summary.stats.wall,
    atom_concat(Prefix, '_inf_cold', KCold),
    atom_concat(Prefix, '_inf_min', KMin),
    atom_concat(Prefix, '_inf_med', KMed),
    atom_concat(Prefix, '_inf_mean', KMean),
    atom_concat(Prefix, '_wall_min_ms', KWMin),
    atom_concat(Prefix, '_wall_med_ms', KWMed),
    atom_concat(Prefix, '_wall_mean_ms', KWMean),
    WallMin is Wall.min*1000, WallMed is Wall.median*1000,
    WallMean is Wall.mean*1000,
    dict_create(Dict, metrics,
                [KCold-Cold.inferences,KMin-Inf.min,KMed-Inf.median,KMean-Inf.mean,
                 KWMin-WallMin,KWMed-WallMed,KWMean-WallMean]).

command_fields(result(Cold, Summary), Dict) :-
    Wall = Summary.stats.wall, Rss = Summary.stats.rss,
    ColdMs is Cold.wall*1000,
    MinMs is Wall.min*1000, MedMs is Wall.median*1000,
    MeanMs is Wall.mean*1000,
    Dict = _{command_wall_cold_ms:ColdMs,
             command_wall_min_ms:MinMs,
             command_wall_med_ms:MedMs,
             command_wall_mean_ms:MeanMs,
             command_rss_med_kib:Rss.median}.

command_k(candidates(_, K), K).
command_k(best_effort, 1).

expected_json(completed(P,D), _{outcome:completed,placed:P,dropped:D}).
expected_json(setup_failed, _{outcome:setup_failed}).
command_json(candidates(strict,K), _{kind:candidates,drop:strict,k:K}).
command_json(best_effort, _{kind:best_effort}).

repo_path(Rel, Path) :- repo_root(Root), directory_file_path(Root, Rel, Path).

identity_document(Specs, _{tool:'crosswordsmith-arrange-greedy-identity',schema:2,
                           results:Rows}) :-
    repo_path(crosswordsmith, Exe),
    maplist(identity_spec(Exe), Specs, Rows).

identity_spec(Exe, spec(Id,F0,S,Fr,M,Seed,C,E,_I,_W,_T), Row) :-
    repo_path(F0, F),
    greedy_subjects:identity_row(Id,F0,S,Fr,M,Seed,C,E,Exe,Row0),
    % identity_row needs the CLI-visible relative fixture path while loading is
    % rooted here; replace its load by running from the repo root via absolute F.
    ( exists_file(F) -> Row = Row0 ; throw(error(existence_error(file,F),_)) ).

emit(json, Rows) :- !,
    swi_version(Version),
    json_write_dict(current_output,
                    _{tool:'crosswordsmith-arrange-greedy-bench',swi_prolog:Version,
                      primary_metric:sweep_inf_med,results:Rows}, [width(100)]), nl.
emit(text, Rows) :- !,
    format("~w~t~38|~t~w~14+~t~w~14+~t~w~14+~t~w~12+~n",
           [rung,construction_inf,sweep_inf,postprocess_inf,command_ms]),
    forall(member(R,Rows), print_row(R)).
emit(csv, Rows) :- !,
    format("rung,tier,construction_inf,sweep_inf,postprocess_inf,command_wall_med_ms,command_rss_med_kib~n"),
    forall(member(R,Rows), print_csv(R)).
emit(Other, _) :- throw(error(domain_error(greedy_bench_format,Other),_)).

print_row(R) :-
    M=R.metrics,
    format("~w~t~38|~t~D~14+~t~D~14+~t~D~14+~t~1f~12+~n",
           [R.rung,M.construction_inf_med,M.sweep_inf_med,
            M.postprocess_inf_med,M.command_wall_med_ms]).
print_csv(R) :-
    M=R.metrics,
    format("~w,~w,~d,~d,~d,~3f,~0f~n",
           [R.rung,R.tier,M.construction_inf_med,M.sweep_inf_med,
            M.postprocess_inf_med,M.command_wall_med_ms,M.command_rss_med_kib]).

swi_version(Version) :-
    current_prolog_flag(version_data, swi(Ma,Mi,Pa,_)),
    format(atom(Version), '~d.~d.~d', [Ma,Mi,Pa]).
