#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(json), [json_write_dict/2]).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir),file_type(directory)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe).
:- initialization(main, main).

subject(hard_15x15_34_topright, hard, 'fixtures/ladder_15x15_34w.pl', 15, 34,
        topright).
subject(hard_15x15_36_topright, hard, 'fixtures/ladder_15x15_36w.pl', 15, 36,
        topright).
subject(light_09x09_08_topright, light, 'fixtures/ladder_09x09_08w.pl', 9, 8,
        topright).
subject(light_15x15_12_topright, light, 'fixtures/ladder_15x15_12w.pl', 15, 12,
        topright).
subject(light_21x21_25_topright, light, 'fixtures/ladder_21x21_25w.pl', 21, 25,
        topright).
subject(easy_bundled_topleft, light, 'fixtures/bundled_17_clues.pl', 17, 6,
        topleft_across).

main :-
    findall(Row, subject_row(Row), Rows),
    json_write_dict(user_output,
                    _{probe:'P-C0',repeat_verified:true,rows:Rows}, [width(0)]),
    nl.

subject_row(Row) :-
    subject(Name, Class, Fixture, Grid, Count, Corner),
    probe_arrange:load_fixture(Fixture, Count, Words),
    run(Name, warmup, Words, Grid, Corner, _, _, _),
    run(Name, first, Words, Grid, Corner, Result1, Stats1, Timing1),
    run(Name, repeat, Words, Grid, Corner, Result2, Stats2, Timing2),
    result_fields(Result1, Outcome, Reward, Signature),
    result_fields(Result2, Outcome2, Reward2, Signature2),
    require_equal(Name-result, Outcome-Reward-Signature,
                  Outcome2-Reward2-Signature2),
    del_dict(decision_trace, Stats1, Trace1, ComparedStats1),
    del_dict(decision_trace, Stats2, Trace2, ComparedStats2),
    require_equal(Name-stats, ComparedStats1, ComparedStats2),
    require_trace_equal(Name-decision_trace, Trace1, Trace2),
    require_equal(Name-inferences, Timing1.inferences, Timing2.inferences),
    ratio(Stats1.repeated_dead_entries, Stats1.nodes, DeadEntryPct),
    ratio(Stats1.repeated_dead_nodes, Stats1.nodes, DeadWorkPct),
    ratio(Stats1.parent_dedup_nodes, Stats1.repeated_dead_nodes, CapturePct),
    Row = _{name:Name,class:Class,fixture:Fixture,grid:Grid,corner:Corner,
            outcome:Outcome,reward:Reward,layout_signature:Signature,
            nodes:Stats1.nodes,decisions:Stats1.decisions,
            crossing_proofs:Stats1.crossing_proofs,
            duplicate_proofs:Stats1.duplicate_proofs,
            duplicate_children:Stats1.duplicate_children,
            canonical_states:Stats1.canonical_states,
            duplicate_states:Stats1.duplicate_states,
            dead_states:Stats1.dead_states,
            repeated_dead_entries:Stats1.repeated_dead_entries,
            repeated_dead_nodes:Stats1.repeated_dead_nodes,
            parent_dedup_hits:Stats1.parent_dedup_hits,
            parent_dedup_nodes:Stats1.parent_dedup_nodes,
            dead_entry_pct:DeadEntryPct,dead_work_pct:DeadWorkPct,
            parent_capture_pct:CapturePct,
            exact_hash_hits:Stats1.exact_hash_hits,
            hash_collisions:Stats1.hash_collisions,
            measured_inferences:Timing1.inferences,
            repeat_exact:true}.

run(Name, Pass, Words, Grid, Corner, Result, Stats, Timing) :-
    format(user_error,
           'probe-arrange heartbeat probe=P-C0 control=~w pass=~w phase=start~n',
           [Name,Pass]), flush_output(user_error),
    probe_arrange:twin_corner(Words,Grid,Corner,null,full,none,null,
                              Result,run(Stats,Timing)),
    format(user_error,
           'probe-arrange heartbeat probe=P-C0 control=~w pass=~w phase=done~n',
           [Name,Pass]), flush_output(user_error).

result_fields(result(placed,ok,false,null,Numbered,Reward), placed, Reward,
              Signature) :-
    probe_arrange:layout_signature(Numbered, Signature).
result_fields(result(infeasible,exhausted,false,null,[],null), infeasible, null,
              null).

ratio(_, 0, null) :- !.
ratio(Numerator, Denominator, Percent) :-
    Percent is Numerator * 100.0 / Denominator.

require_equal(_Label, A, B) :- A =@= B, !.
require_equal(Label, A, B) :-
    throw(error(pc0_repeat_mismatch(Label,A,B), _)).

require_trace_equal(_Label, A, B) :- A == B, !.
require_trace_equal(Label, A, B) :-
    variant_sha1(A, HashA), variant_sha1(B, HashB),
    throw(error(pc0_repeat_mismatch(Label,HashA,HashB), _)).
