#!/usr/bin/env swipl
% Run the measurement-only A-G2 direct transpose premise probe.

:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/3]).
:- use_module(library(filesex), [directory_file_path/3]).

:- prolog_load_context(directory, BenchDir),
   absolute_file_name('..', Root,
                      [relative_to(BenchDir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(BenchDir, 'probe_arrange/g2_transpose.pl', Probe),
   use_module(Probe).

:- initialization(main, main).

main :-
    g2_transpose_probe:probe_all(Doc),
    json_write_dict(current_output, Doc, [width(0)]), nl,
    ( Doc.gate == pass -> halt(0) ; halt(1) ).
