:- use_module(library(plunit)).
:- use_module('../benchmarks/probe_arrange/probe_arrange.pl').
:- use_module('../benchmarks/probe_arrange/ad1_buckets.pl').

:- begin_tests(probe_arrange).

test(fixture_exact_count) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    length(Words, 6).

test(fixture_wrong_count_throws,
     [throws(error(fixture_word_count(_, 7, 6), _))]) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 7, _).

test(authority_tiny_budget_is_censored_not_proven) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:authority_operation(Words, 17, 1000, null,
        result(not_proven,budget,true,null,[],null), _).

test(authority_matches_direct_product_tiny_and_500m) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    forall(member(Budget, [1000,500000000]),
           ( probe_arrange:authority_operation(Words,17,Budget,null,Authority,_),
             crosswordsmith_core:set_search_seed(-1),
             crosswordsmith_arrange:arrange_best_layout(
                 Words,17,Budget,DirectN,DirectR,DirectOutcome),
             direct_result(DirectOutcome,DirectN,DirectR,Direct),
             same_product_result(Authority,Direct) )).

test(twin_replays_both_nontranspose_corners,
     [forall(member(Corner, [topleft_across,topright]))]) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:authority_corner(Words,17,Corner,500000000,null,AR,_),
    probe_arrange:twin_corner(Words,17,Corner,null,lean,none,null,TR,run(Stats,_)),
    AR = result(placed,ok,false,_,AN,Reward),
    TR = result(placed,ok,false,null,TN,Reward),
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(TN, Sig),
    Stats.decisions =:= Stats.places.

test(twin_replays_seeded_layout) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:authority_corner(Words,17,topright,500000000,7,AR,_),
    probe_arrange:twin_corner(Words,17,topright,7,full,none,null,TR,run(Stats,_)),
    AR = result(placed,ok,false,_,AN,Reward),
    TR = result(placed,ok,false,null,TN,Reward),
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(TN, Sig),
    integer(Stats.support_transitions).

test(full_observer_preserves_exact_decision_order) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:twin_corner(Words,17,topleft_across,null,lean,none,null,
                              Lean,run(LS,_)),
    probe_arrange:twin_corner(Words,17,topleft_across,null,full,none,null,
                              Full,run(FS,_)),
    same_twin_result(Lean, Full),
    LS.decision_trace == FS.decision_trace,
    length(FS.decision_trace, FS.decisions).

test(ad1_assoc_and_direct_count_selection_decision_traces_match,
     [forall((member(Fixture-Grid-Count,
                     ['fixtures/ladder_09x09_08w.pl'-9-8,
                      'fixtures/ladder_15x15_12w.pl'-15-12,
                      'fixtures/ladder_15x15_32w.pl'-15-32,
                      'fixtures/ladder_21x21_80w.pl'-21-80]),
              member(Corner, [topleft_across,topright])))]) :-
    probe_arrange:load_fixture(Fixture, Count, Words),
    probe_arrange_ad1:ad1_corner(
        assoc, Words, Grid, Corner, null, Assoc, AssocCounts, AssocDecisions),
    probe_arrange_ad1:ad1_corner(
        direct, Words, Grid, Corner, null, Direct, DirectCounts, DirectDecisions),
    AssocCounts == DirectCounts,
    AssocDecisions == DirectDecisions,
    Assoc = result(AssocPlaced, Reward),
    Direct = result(DirectPlaced, Reward),
    probe_arrange:layout_signature(AssocPlaced, Signature),
    probe_arrange:layout_signature(DirectPlaced, Signature).

test(canonical_state_key_is_absolute_and_orientation_sensitive) :-
    Words = [['ONE'],['TWO']],
    probe_arrange:init_holder(full,none,null,Words,H),
    setup_call_cleanup(
        true,
        ( probe_arrange:canonical_state_key(
              H, 5, [['TWO']], [pw('ONE',[],[],across,3,2,4,_)], Across),
          probe_arrange:canonical_state_key(
              H, 5, [['TWO']], [pw('ONE',[],[],down,3,2,12,_)], Down),
          Across == state(5,[2],[p(1,2,across)]),
          Down == state(5,[2],[p(1,2,down)]),
          Across \== Down ),
        probe_arrange:cleanup_holder(H)).

test(fingerprint_bucket_never_establishes_equality) :-
    Words = [['ONE']], Key = state(5,[],[p(1,1,across)]),
    Other = state(5,[],[p(1,1,down)]),
    probe_arrange:init_holder(full,none,null,Words,H),
    setup_call_cleanup(
        true,
        ( probe_arrange:state_fingerprint(Key, fp(H1,H2)),
          arg(19,H,Run),
          assertz(probe_arrange:pc0_exact(Run,seen,H1,H2,Other)),
          probe_arrange:exact_table_lookup(H,seen,Key,false),
          arg(28,H,1) ),
        probe_arrange:cleanup_holder(H)).

test(semantic_cutoff_never_marks_dead) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:twin_corner(Words,17,topleft_across,null,full,nodes,0,
        result(not_proven,budget,true,null,[],null),run(Stats,_)),
    Stats.nodes =:= 0,
    Stats.dead_states =:= 0,
    Stats.repeated_dead_entries =:= 0.

test(twin_operation_preserves_shared_seed_stream) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:authority_operation(Words,17,500000000,7,Authority,_),
    once(probe_arrange:twin_operation(Words,17,7,none,Twin,_,_)),
    Authority = result(placed,ok,false,_,AN,Reward),
    Twin = result(placed,ok,false,null,TN,Reward),
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(TN, Sig).

test(seed_cleanup_on_failure_exception_and_interruption) :-
    \+ probe_arrange:with_search_seed(7, fail),
    \+ crosswordsmith_core:current_search_seed(_),
    catch(probe_arrange:with_search_seed(7, throw(probe_error)), probe_error, true),
    \+ crosswordsmith_core:current_search_seed(_),
    catch(probe_arrange:with_search_seed(7, throw(time_limit_exceeded)),
          time_limit_exceeded, true),
    \+ crosswordsmith_core:current_search_seed(_).

test(decision_limit_is_censored) :-
    probe_arrange:load_fixture('fixtures/bundled_17_clues.pl', 6, Words),
    probe_arrange:twin_corner(Words,17,topleft_across,null,lean,decisions,1,
        result(not_proven,budget,true,null,[],null),run(Stats,_)),
    Stats.decisions =:= 1.

test(trace_row_keeps_configured_cutoff_separate_from_termination) :-
    Meta = _{limit_kind:inferences,cutoff:500000000,operation_id:op,
             attempt_index:0,fixture:'f.pl',fixture_seed:11,search_seed:null,
             corner:null,arm:control},
    null_probe_stats(Stats),
    Timing = _{inferences:10,wall:0.1,cpu:0.1},
    probe_arrange:trace_row(authority,Meta,
        result(placed,ok,false,10,[pw('A',[],[],across,0,1,1,1)],12),
        Stats,Timing,"commit","10.1.10",Row),
    Row.cutoff =:= 500000000,
    Row.termination == ok,
    Row.outcome == placed.

test(trace_row_interruption_preserves_configured_cutoff) :-
    Meta = _{limit_kind:inferences,cutoff:1000,operation_id:op,
             attempt_index:0,fixture:'f.pl',fixture_seed:11,search_seed:null,
             corner:null,arm:control},
    null_probe_stats(Stats),
    Timing = _{inferences:null,wall:null,cpu:null},
    probe_arrange:trace_row(authority,Meta,
        result(interrupted,interrupted,true,null,[],null),
        Stats,Timing,"commit","10.1.10",Row),
    Row.cutoff =:= 1000,
    Row.termination == interrupted,
    Row.censored == true.

null_probe_stats(_{nodes:null,decisions:null,places:null,unplaces:null,
                    wipeouts:null,max_depth:null,state_entries_max:null,
                    letter_cells_max:null,boundary_cells_max:null,
                    support_transitions:null,crossing_proofs:null,
                    duplicate_proofs:null,duplicate_children:null,
                    duplicate_states:null,repeated_dead_entries:null,
                    repeated_dead_nodes:null,parent_dedup_hits:null,
                    parent_dedup_nodes:null,dead_states:null,
                    canonical_states:null,exact_hash_hits:null,
                    hash_collisions:null}).

direct_result(placed, N, R, result(placed,N,R)).
direct_result(not_proven, _, _, result(not_proven,[],null)).
direct_result(infeasible, _, _, result(infeasible,[],null)).

same_product_result(result(Outcome,_,_,_,AN,AR), result(Outcome,DN,DR)) :-
    AR == DR,
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(DN, Sig).

same_twin_result(result(Outcome,Termination,Censored,_,AN,Reward),
                 result(Outcome,Termination,Censored,_,BN,Reward)) :-
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(BN, Sig).

:- end_tests(probe_arrange).
