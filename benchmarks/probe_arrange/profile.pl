#!/usr/bin/env swipl
:- use_module(library(prolog_profile), [profile/1]).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir),file_type(directory)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe).
:- initialization(main, main).

main :-
    current_prolog_flag(argv, [Fixture,G0,N0,Seed0,Corner|_]),
    atom_number(G0, Grid), atom_number(N0, Count),
    ( Seed0 == none -> Seed=null ; atom_number(Seed0, Seed) ),
    probe_arrange:load_fixture(Fixture, Count, Words),
    profile(probe_arrange:twin_corner(Words,Grid,Corner,Seed,none,none,null,_,_)).
