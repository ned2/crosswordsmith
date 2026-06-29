% tests/arrange.plt - plunit test suite for arrange.pl (Flavour-A engine).
%
% Assumes crossword.pl and arrange.pl have been consulted into `user` by the
% runner (tests/run_tests.pl) before this file is loaded. Run via:
%
%     ./run_tests.sh        (or)        make test
%
% Covers Phase 1 (the layout_reward oracle), Phase 1.5's keeper scoring
% helpers, and Phase 2-3 (strict construct + rescore + emit, size framing,
% the 3-outcome strict contract, and emit determinism). The byte-exact emit
% regression lives in the golden files (tests/golden/arrange_*.json).

:- use_module(library(plunit)).
:- use_module(library(http/json)).

arrange_bundled(Words) :- load_clues('fixtures/bundled_17_clues.pl', Words).
arrange_toc(Words)     :- load_clues('fixtures/toc_demo.pl', Words).


:- begin_tests(arrange).

% --- Phase 1: scoring oracle -------------------------------------------------

% check_target is ceil(L/2): the per-word checking cap.
test(check_target_odd)  :- check_target(5, 3).
test(check_target_even) :- check_target(4, 2).
test(check_target_one)  :- check_target(1, 1).
test(check_target_two)  :- check_target(2, 1).

% The carried reward equals the from-scratch oracle on the final (numbered)
% layout (AC-ARR-9): numbering adds only `num`, so layout_reward is unchanged.
test(layout_reward_matches_after_numbering) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, R, Outcome),
    Outcome == placed,
    layout_reward(5, 1, Numbered, R2),
    R =:= R2.

% Regression pin: bundled_17 on grid 17 places all 6, reward 60 (cap inert
% here, so reward == 6*Sum-checked). Update deliberately if construction
% changes (mirror the golden update).
test(bundled_strict_reward, [true(R =:= 60)]) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, _Numbered, R, placed).

% cap_binding_count on bundled_17 is 0 (the reachability caveat in miniature).
test(bundled_cap_binding_zero, [true(CB =:= 0)]) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, _R, placed),
    cap_binding_count(Numbered, CB).

% --- Phase 2: strict construct + rescore (place-all-or-fail) ------------------

test(strict_places_all_words) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, _R, Outcome),
    Outcome == placed,
    length(Numbered, 6).

% AC-ARR-4: best-of-corners reward is never below a single first-solution.
test(reward_ge_single_corner) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, _N, RBest, placed),
    once(find_crossword(mrv_inc, 17, Words, topleft_across, _G, P0)),
    layout_reward(5, 1, P0, R0),
    RBest >= R0.

% AC-ARR-1 outcome (c)/(b): too-small grid -> infeasible, no stdout layout.
test(strict_grid_too_small_infeasible) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 3, Numbered, _R, Outcome),
    Outcome == infeasible,
    Numbered == [].

% AC-ARR-1 outcome (b): words sharing no letters -> infeasible, named.
test(strict_isolated_words_infeasible) :-
    arrange_best_layout([['ABC'], ['DEF']], 9, _N, _R, Outcome),
    Outcome == infeasible.

test(unplaceable_names_isolated, [true(Bad == ['ABC', 'DEF'])]) :-
    unplaceable_words([['ABC'], ['DEF']], Bad).

test(unplaceable_empty_when_interlocking, [true(Bad == [])]) :-
    arrange_bundled(Words),
    unplaceable_words(Words, Bad).

% --- Phase 3: emit framing ---------------------------------------------------

% fixed: exactly N x N, all words present, valid JSON.
test(emit_fixed_is_n_by_n) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S), emit_arrange(Numbered, Words, 17, fixed)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Rows), length(Rows, 17),
    get_dict(words, Dict, Ws), length(Ws, 6).

% max: cropped square (<= N), content anchored at (0,0).
test(emit_max_shrinks_and_anchors) :-
    arrange_toc(Words),
    arrange_best_layout(Words, 25, Numbered, _R, placed),
    with_output_to(string(S), emit_arrange(Numbered, Words, 25, max)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, GL),
    GL < 25,                                   % shrunk from the N=25 ceiling
    get_dict(words, Dict, Ws),
    once(( member(W,  Ws), get_dict(cells, W,  Cs),  member([0, _], Cs) )),  % touches row 0
    once(( member(W2, Ws), get_dict(cells, W2, Cs2), member([_, 0], Cs2) )). % touches col 0

% The crop translation arithmetic: original (row 5, col 7) on a 17-grid,
% cropped with origin (2,3), lands at [3,4].
test(cropped_coord_translates, [true(Coord == [3, 4])]) :-
    Cell is 5 * 17 + 7 + 1,
    cropped_coord(17, 2, 3, Cell, Coord).

% AC-ARR-6: identical input -> byte-identical emit (no nondeterminism).
test(emit_is_deterministic) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S1), emit_arrange(Numbered, Words, 17, max)),
    with_output_to(string(S2), emit_arrange(Numbered, Words, 17, max)),
    S1 == S2.

% --- Phase 4: best-effort (drop) ---------------------------------------------

% When the grid fits everything, best-effort places all and drops nothing.
test(best_effort_places_all_when_fits) :-
    arrange_bundled(Words),
    arrange_best_effort(Words, 17, Numbered, _R, NP, Dropped),
    NP =:= 6, Dropped == [], length(Numbered, 6).

% AC-ARR-2: best-effort succeeds where strict fails (isolated words sharing no
% letters), placing a maximal subset and reporting the dropped remainder.
test(best_effort_succeeds_on_isolated) :-
    arrange_best_effort([['ABC'], ['DEF']], 9, Numbered, _R, NP, Dropped),
    NP =:= 1, length(Dropped, 1), length(Numbered, 1).

% On a too-tight grid best-effort drops, and every word is placed XOR dropped.
test(best_effort_drops_and_partitions) :-
    arrange_bundled(Words), length(Words, Total),
    arrange_best_effort(Words, 11, _N, _R, NP, Dropped),
    length(Dropped, ND),
    NP >= 1, NP < Total, NP + ND =:= Total.

:- end_tests(arrange).
