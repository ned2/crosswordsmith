% tests/arrange.plt - plunit test suite for arrange.pl (Flavour-A engine).
%
% Assumes core.pl and arrange.pl have been consulted into `user` by the
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

% R8 (revamp audit): --check-target N (check_target_override/1) lowers the
% per-word target to min(ceil(L/2), N); N above ceil(L/2) is a no-op (only lowers).
% (check_target_override/1 is asserted into `user` - where arrange.pl's
% check_target/2 reads it - because plunit runs these in module plunit_arrange.)
test(check_target_override_lowers, [true(T =:= 2)]) :-
    setup_call_cleanup(assertz(user:check_target_override(2)),
                       check_target(9, T),                 % min(ceil(9/2)=5, 2) = 2
                       retractall(user:check_target_override(_))).
test(check_target_override_only_lowers, [true(T =:= 3)]) :-
    setup_call_cleanup(assertz(user:check_target_override(10)),
                       check_target(6, T),                 % min(ceil(6/2)=3, 10) = 3
                       retractall(user:check_target_override(_))).
% End to end: bundled_17's cap is inert by default (cap_binding_count 0); with
% the target lowered to 1, every interlocking placed word reaches it.
test(check_target_makes_cap_bind) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered0, _R0, placed),
    cap_binding_count(Numbered0, CB0), CB0 =:= 0,
    setup_call_cleanup(assertz(user:check_target_override(1)),
                       ( arrange_best_layout(Words, 17, Numbered1, _R1, placed),
                         cap_binding_count(Numbered1, CB1) ),
                       retractall(user:check_target_override(_))),
    CB1 > 0.

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

% AC-ARR-1 outcome (b): too-small grid -> infeasible, no stdout layout.
test(strict_grid_too_small_infeasible) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 3, Numbered, _R, Outcome),
    Outcome == infeasible,
    Numbered == [].

% AC-ARR-1 outcome (c) + AC-ARR-10 (R4, revamp audit): when the per-corner
% inference budget is exhausted before any corner completes, the outcome is
% `not_proven` - DISTINCT from `infeasible`. A budget too small to finish even
% one corner (1000 inferences << a 6-word grid-17 construction) forces it.
test(strict_budget_exhausted_not_proven) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, 1000, Numbered, Reward, Outcome),
    Outcome == not_proven,
    Numbered == [],
    Reward =:= -1.

% The same input completes (placed) under the real budget - so `not_proven`
% above is the budget biting, not genuine infeasibility (the (b) vs (c) split).
test(strict_budget_not_proven_is_not_infeasible) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, _N, _R, Outcome),
    Outcome == placed.

% R7 (revamp audit): the 4-corner sweep shares ONE inference budget, so a HARD
% input (toc_demo at size 15 - the finding's repro, which the full budget cannot
% finish) costs ~Budget total, not Budget x 4. With a 5M budget too small for
% any corner to complete, the total inferences consumed stay near one budget,
% well under the 4x the per-corner version would burn.
test(arrange_budget_shared_across_corners) :-
    arrange_toc(Words),
    B = 5_000_000,
    statistics(inferences, I0),
    arrange_best_layout(Words, 15, B, _N, _R, Outcome),
    statistics(inferences, I1),
    Used is I1 - I0,
    Outcome == not_proven,
    Used < 2 * B.

% AC-ARR-1 outcome (b): words sharing no letters -> infeasible, named.
test(strict_isolated_words_infeasible) :-
    arrange_best_layout([['ABC'], ['DEF']], 9, _N, _R, Outcome),
    Outcome == infeasible.

% Regression for P2: construct_one/7 must NOT swallow a genuine error from the
% search and report it as a normal outcome. call_with_inference_limit/3 handles
% the budget itself (binds inference_limit_exceeded and succeeds) and re-throws
% real exceptions; infeasibility is a search FAILURE, never a throw. So a real
% error (here a non-integer GridLen forcing is/2 to raise) must propagate, not
% be reclassified as `exhausted`. The old broad catch/3 returned `exhausted`.
% (construct_fragment_one/6 shares the identical no-catch pattern.)
test(construct_one_propagates_genuine_error,
     [throws(error(type_error(evaluable, not_an_int/0), _))]) :-
    construct_one(topleft_across, [['CAT', _{}]], not_an_int, 5, 3, 1_000_000, _Res).

test(unplaceable_names_isolated, [true(Bad == ['ABC', 'DEF'])]) :-
    unplaceable_words([['ABC'], ['DEF']], Bad).

test(unplaceable_empty_when_interlocking, [true(Bad == [])]) :-
    arrange_bundled(Words),
    unplaceable_words(Words, Bad).

% R1 (revamp audit): a multi-word/hyphenated answer is placed by its letters
% ONLY - the separator (space or hyphen) is an enumeration marker, not a grid
% cell. WELL-BEING is a 9-cell run (not 10), no placed cell is '-', and the
% original answer atom (with the hyphen) is preserved for export enumeration.
test(hyphenated_answer_strips_separator, [nondet]) :-
    arrange_best_layout([['WELL-BEING'], ['BELOW']], 15, Numbered, _R, placed),
    member(W, Numbered), get_dict(answer, W, 'WELL-BEING'), !,
    get_dict(len, W, Len), Len =:= 9,
    get_dict(cells, W, Cells), length(Cells, 9),
    get_dict(letters, W, Ls), \+ memberchk('-', Ls), \+ memberchk(' ', Ls).

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

% --- Phase 5: fragment-grid seeding (anchors) --------------------------------

% Two long words of bundled_17 pinned (an across + a down crossing at a 'T'),
% used by the partial-seed and best-effort tests below.
bundled_two_word_fragment(Frags) :-
    fragment_dict_words(
        _{gridLength:17, words:[
            _{answer:"GNOSTIC GOSPELS", direction:"down",
              cells:[[0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[7,3],
                     [8,3],[9,3],[10,3],[11,3],[12,3],[13,3]]},
            _{answer:"ETERNAL RETURN", direction:"across",
              cells:[[4,2],[4,3],[4,4],[4,5],[4,6],[4,7],[4,8],
                     [4,9],[4,10],[4,11],[4,12],[4,13],[4,14]]} ]},
        17, Frags).

% Parse: a 4-cell across word at [0,3] on a 17-grid -> dir across, start cell 4
% (0*17+3+1), and the declared run as 1-based cell numbers.
test(fragment_parse_start_and_run,
     [true(Frags == [frag('BIAS', across, 4, [4, 5, 6, 7])])]) :-
    fragment_dict_words(
        _{gridLength:17, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,3],[0,4],[0,5],[0,6]]}]},
        17, Frags).

% AC-EMIT-2 / AC-ARR-3 (load-bearing): an emitted layout IS a valid fragment;
% re-ingesting it and re-emitting reproduces byte-identical JSON (all words
% pinned -> the remainder is empty, so the seed is the whole layout).
test(fragment_roundtrip_byte_identical) :-
    arrange_bundled(Words),
    arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S1), emit_arrange(Numbered, Words, 17, fixed)),
    atom_json_dict(S1, Dict, []),
    fragment_dict_words(Dict, GL, Frags),
    GL =:= 17,
    arrange_fragment_strict(Words, Frags, 17, Numbered2, _R2, placed),
    with_output_to(string(S2), emit_arrange(Numbered2, Words, 17, fixed)),
    S1 == S2.

% R3 (revamp audit): the round-trip above passes on bundled_17 only because its
% shared-start pairs happen not to flip. benchmark_08 at size 13 DOES flip (the
% strict DFS and the fragment re-solve build a start-1 across/down pair in
% opposite order), so without the canonical (number, direction) emit ordering
% the re-ingest diverges. This is the finding's exact repro.
test(fragment_roundtrip_byte_identical_reordering) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    arrange_best_layout(Words, 13, Numbered, _R, placed),
    with_output_to(string(S1), emit_arrange(Numbered, Words, 13, fixed)),
    atom_json_dict(S1, Dict, []),
    fragment_dict_words(Dict, GL, Frags),
    GL =:= 13,
    arrange_fragment_strict(Words, Frags, 13, Numbered2, _R2, placed),
    with_output_to(string(S2), emit_arrange(Numbered2, Words, 13, fixed)),
    S1 == S2.

% AC-FRAG-3: pinned words appear at exactly their fragment cells; the engine
% places the remaining four around them (all six end up placed).
test(fragment_pins_preserved_partial) :-
    arrange_bundled(Words),
    bundled_two_word_fragment(Frags),
    arrange_fragment_strict(Words, Frags, 17, Numbered, _R, placed),
    length(Numbered, 6),
    once(( member(Wg, Numbered), get_dict(answer, Wg, 'GNOSTIC GOSPELS') )),
    get_dict(cells, Wg, Cg), maplist(cell_coord(17), Cg, RCg),
    RCg == [[0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[7,3],
            [8,3],[9,3],[10,3],[11,3],[12,3],[13,3]],
    once(( member(We, Numbered), get_dict(answer, We, 'ETERNAL RETURN') )),
    get_dict(cells, We, Ce), maplist(cell_coord(17), Ce, RCe),
    RCe == [[4,2],[4,3],[4,4],[4,5],[4,6],[4,7],[4,8],
            [4,9],[4,10],[4,11],[4,12],[4,13],[4,14]].

% AC-FRAG-1: a fragment word not present in --input is rejected up front.
test(fragment_word_not_in_input_throws,
     [throws(error(fragment_word_not_in_input('ZZZTOP'), _))]) :-
    arrange_bundled(Words),
    fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"ZZZTOP", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3],[0,4],[0,5]]}]},
        9, Frags),
    arrange_fragment_strict(Words, Frags, 9, _, _, _).

% AC-FRAG-2 (overlap): two pinned words demanding different letters in a shared
% cell are reported before any search, naming the clashing cell and letters.
test(fragment_letter_clash_throws,
     [throws(error(fragment_letter_clash('FLOW', [0,0], 'B', 'F'), _))]) :-
    arrange_bundled(Words),
    fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3]]},
            _{answer:"FLOW", direction:"down",   cells:[[0,0],[1,0],[2,0],[3,0]]}]},
        9, Frags),
    arrange_fragment_strict(Words, Frags, 9, _, _, _).

% AC-FRAG-2 (self-inconsistent geometry): cells that are not a straight run of
% the answer's length are rejected up front.
test(fragment_cells_inconsistent_throws,
     [throws(error(fragment_cells_inconsistent('BIAS'), _))]) :-
    arrange_bundled(Words),
    fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2]]}]},  % 3 cells, 4 letters
        9, Frags),
    arrange_fragment_strict(Words, Frags, 9, _, _, _).

% A fragment may not pin the same answer twice.
test(fragment_duplicate_answer_throws,
     [throws(error(fragment_duplicate_answer('BIAS'), _))]) :-
    arrange_bundled(Words),
    fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3]]},
            _{answer:"BIAS", direction:"down",   cells:[[0,0],[1,0],[2,0],[3,0]]}]},
        9, Frags),
    arrange_fragment_strict(Words, Frags, 9, _, _, _).

% An off-grid cell is a parse-time error.
test(fragment_invalid_cell_throws,
     [throws(error(fragment_invalid_cell('BIAS', [0,9]), _))]) :-
    fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,9],[0,10],[0,11],[0,12]]}]},
        _, _).

% Best-effort with a fragment: the seed is pinned, the rest greedily placed;
% on a roomy grid nothing is dropped and the pins are kept.
test(fragment_best_effort_places_all) :-
    arrange_bundled(Words),
    bundled_two_word_fragment(Frags),
    arrange_fragment_best_effort(Words, Frags, 17, Numbered, _R, NP, Dropped),
    NP =:= 6, Dropped == [], length(Numbered, 6).

% The fragment's gridLength sets N; an explicit size is redundant unless it
% disagrees, which is an error (design-spec §6.6).
test(fragment_size_reconcile) :-
    reconcile_fragment_size(17, none, 17),
    reconcile_fragment_size(17, 17, 17),
    catch(reconcile_fragment_size(17, 15, _),
          error(fragment_size_mismatch(17, 15), _), true).

% --- Phase 6: candidates (diverse layouts) -----------------------------------

% AC-ARR-7: --candidates K returns up to K layouts; bundled_17 yields 3 distinct
% full layouts.
test(candidates_returns_k_full_layouts) :-
    arrange_bundled(Words),
    arrange_candidates(Words, 17, strict, 3, Layouts, Returned),
    Returned =:= 3,
    length(Layouts, 3),
    forall(member(L, Layouts), length(L, 6)).

% AC-ARR-7: the returned layouts are pairwise >= tau apart (translation-invariant
% placement distance).
test(candidates_pairwise_tau_apart) :-
    arrange_bundled(Words), length(Words, Total),
    candidate_tau_pct(TauPct),
    arrange_candidates(Words, 17, strict, 3, Layouts, _),
    maplist([L, A]>>placement_assoc(L, 17, A), Layouts, Assocs),
    forall( ( nth0(I, Assocs, A1), nth0(J, Assocs, A2), I < J ),
            ( pos_diff_count(A1, A2, Diff), Diff * 100 >= TauPct * Total ) ).

% AC-ARR-7: fewer than K only when fewer >= tau-distinct layouts exist; on the
% 8-word benchmark the greedy breadth yields just two distinct full layouts.
test(candidates_short_return_when_fewer_distinct) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    arrange_candidates(Words, 13, strict, 5, Layouts, Returned),
    Returned < 5, Returned >= 1, length(Layouts, Returned).

% --candidates 1 yields exactly the single best layout.
test(candidates_one_is_single_best) :-
    arrange_bundled(Words),
    arrange_candidates(Words, 17, strict, 1, Layouts, 1),
    Layouts = [L], length(L, 6).

% AC-ARR-6 / INV-2: candidate emission is byte-identical across runs.
test(candidates_emit_deterministic) :-
    arrange_bundled(Words),
    arrange_candidates(Words, 17, strict, 3, Layouts, _),
    with_output_to(string(S1), emit_candidates(Layouts, Words, 17, fixed)),
    with_output_to(string(S2), emit_candidates(Layouts, Words, 17, fixed)),
    S1 == S2.

% The emitted candidates form a JSON array, each element a standalone layout.
test(candidates_emit_is_json_array) :-
    arrange_bundled(Words),
    arrange_candidates(Words, 17, strict, 3, Layouts, _),
    with_output_to(string(S), emit_candidates(Layouts, Words, 17, fixed)),
    atom_json_dict(S, Arr, []),
    is_list(Arr), length(Arr, 3),
    forall(member(D, Arr),
           ( get_dict(gridLength, D, 17), get_dict(words, D, Ws), length(Ws, 6) )).

% Placement distance is translation-invariant: a layout and its one-row shift
% have distance 0 (the same crossword, just moved).
test(candidates_distance_translation_invariant) :-
    L1 = [ word{answer:'ABC', dir:across, start:1,  cells:[1,2,3]},
           word{answer:'ADE', dir:down,   start:1,  cells:[1,18,35]} ],
    L2 = [ word{answer:'ABC', dir:across, start:18, cells:[18,19,20]},
           word{answer:'ADE', dir:down,   start:18, cells:[18,35,52]} ],
    placement_assoc(L1, 17, A1),
    placement_assoc(L2, 17, A2),
    pos_diff_count(A1, A2, 0).

% --- Phase 7: enumerate seam -------------------------------------------------

% AC-ARR-8: --enumerate counts every feasible full placement, matching the old
% --all count (it is all_crossword/5 under the production default strategy, over
% all start corners).
test(enumerate_matches_all_crossword) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    arrange_enumerate(Words, 17, N1),
    default_strategy(Strat),
    all_crossword(Strat, 17, Words, _StartLoc, N2),
    N1 =:= N2,
    N1 > 0.

% --- greedy constructor (moved here from metrics.pl in Phase 3) ---------------

% seed_candidates/2 returns the longest words first (restart diversity). With
% fewer words than the cap K, all are returned, longest-first.
test(seed_candidates_longest_first) :-
    seed_candidates([['AAAA',_{}],['BB',_{}],['CCC',_{}]], Seeds),
    Seeds = [[S1|_], [S2|_], [S3|_]],
    S1 == 'AAAA', S2 == 'CCC', S3 == 'BB'.

:- end_tests(arrange).
