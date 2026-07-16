#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/2]).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir),file_type(directory)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe).
:- initialization(main, main).

control(easy_topleft, 'fixtures/bundled_17_clues.pl', 17, 6, null, topleft_across).
control(hard_topleft, 'fixtures/ladder_09x09_16w.pl', 9, 16, null, topleft_across).
control(seeded_topleft, 'fixtures/bundled_17_clues.pl', 17, 6, 7, topleft_across).
control(corner_topright, 'fixtures/bundled_17_clues.pl', 17, 6, null, topright).

main :-
    findall(Row, control_row(Row), Rows),
    json_write_dict(user_output, _{controls:Rows}, [width(0)]), nl.

control_row(Row) :-
    control(Name, Fixture, Grid, Count, Seed, Corner),
    format(user_error, 'probe-arrange heartbeat control=~w phase=start~n', [Name]),
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange:authority_corner(Words,Grid,Corner,500000000,Seed,Authority,_),
    probe_arrange:twin_corner(Words,Grid,Corner,Seed,none,none,null,None,run(_,NT)),
    probe_arrange:twin_corner(Words,Grid,Corner,Seed,lean,none,null,Lean,run(LS,LT)),
    probe_arrange:twin_corner(Words,Grid,Corner,Seed,full,none,null,Full,run(FS,FT)),
    equivalent(Authority,None,Sig,Reward),
    equivalent(Authority,Lean,Sig,Reward),
    equivalent(Authority,Full,Sig,Reward),
    LS.decisions =:= FS.decisions,
    Row = _{control:Name,fixture:Fixture,corner:Corner,search_seed:Seed,
            outcome:placed,reward:Reward,layout_signature:Sig,
            lean_decisions:LS.decisions,full_decisions:FS.decisions,
            lean_nodes:LS.nodes,full_nodes:FS.nodes,
            max_depth:LS.max_depth,places:LS.places,
            unplaces:LS.unplaces,wipeouts:LS.wipeouts,
            lean_inferences:LT.inferences,full_inferences:FT.inferences,
            uninstrumented_twin_inferences:NT.inferences},
    format(user_error, 'probe-arrange heartbeat control=~w phase=done~n', [Name]).

equivalent(result(placed,ok,false,_,AN,R),
           result(placed,ok,false,_,BN,R), Sig, R) :-
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(BN, Sig).
