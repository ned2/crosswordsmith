:- use_module(library(plunit)).
:- use_module('../benchmarks/probe_arrange/probe_arrange.pl').

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

direct_result(placed, N, R, result(placed,N,R)).
direct_result(not_proven, _, _, result(not_proven,[],null)).
direct_result(infeasible, _, _, result(infeasible,[],null)).

same_product_result(result(Outcome,_,_,_,AN,AR), result(Outcome,DN,DR)) :-
    AR == DR,
    probe_arrange:layout_signature(AN, Sig),
    probe_arrange:layout_signature(DN, Sig).

:- end_tests(probe_arrange).
