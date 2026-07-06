% tests/arrange.plt - plunit test suite for arrange.pl (Flavour-A engine).
%
% Assumes the project is loaded (via load.pl) by the runner
% (tests/run_tests.pl) before this file. arrange's white-box helpers are not
% exported — they are reached as crosswordsmith_arrange:Pred(...) (migration
% plan, Resolved Decision 3). Run via:
%
%     ./run_tests.sh        (or)        make test
%
% Covers Phase 1 (the layout_reward oracle), Phase 1.5's keeper scoring
% helpers, and Phase 2-3 (strict construct + rescore + emit, size framing,
% the 3-outcome strict contract, and emit determinism). The byte-exact emit
% regression lives in the golden files (tests/golden/arrange_*.json).

:- use_module(library(plunit)).
:- use_module(library(json)).
% lists/apply: test-body helpers; explicit so the suite also runs under
% autoload(false) (P11/C5).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(yall)).    % test-body lambdas
:- use_module(library(debug)).   % assertion/1

arrange_bundled(Words) :- load_clues('fixtures/bundled_17_clues.pl', Words).
arrange_toc(Words)     :- load_clues('fixtures/toc_demo.pl', Words).

% statistics(inferences) delta around a deterministically-succeeding goal
% (once/1-pinned so measuring can never leave the goal's choicepoint behind).
arrange_inf_delta(Goal, Delta) :-
    statistics(inferences, I0),
    once(Goal),
    statistics(inferences, I1),
    Delta is I1 - I0.

% One greedy / one strict solve on the bundled fixture with FRESH output
% variables per call, measured (the C1 history-independence lock's probes).
greedy_inf(Words, Delta) :-
    arrange_inf_delta(
        crosswordsmith_arrange:arrange_best_effort(Words, 17, _, _, _, _),
        Delta).
strict_inf(Words, Delta) :-
    arrange_inf_delta(
        crosswordsmith_arrange:arrange_best_layout(Words, 17, _, _, _),
        Delta).


:- begin_tests(arrange).

% --- Phase 1: scoring oracle -------------------------------------------------

% check_target is ceil(L/2): the per-word checking cap.
test(check_target_odd)  :- crosswordsmith_arrange:check_target(5, 3).
test(check_target_even) :- crosswordsmith_arrange:check_target(4, 2).
test(check_target_one)  :- crosswordsmith_arrange:check_target(1, 1).
test(check_target_two)  :- crosswordsmith_arrange:check_target(2, 1).

% R8 (revamp audit): --check-target N (check_target_override/1) lowers the
% per-word target to min(ceil(L/2), N); N above ceil(L/2) is a no-op (only lowers).
% (The override is installed via the exported set_check_target/1 - the only
% sanctioned writer of the module-private dynamic; -1 clears it.)
test(check_target_override_lowers, [true(T =:= 2)]) :-
    setup_call_cleanup(set_check_target(2),
                       crosswordsmith_arrange:check_target(9, T),                 % min(ceil(9/2)=5, 2) = 2
                       set_check_target(-1)).
test(check_target_override_only_lowers, [true(T =:= 3)]) :-
    setup_call_cleanup(set_check_target(10),
                       crosswordsmith_arrange:check_target(6, T),                 % min(ceil(6/2)=3, 10) = 3
                       set_check_target(-1)).
% End to end: bundled_17's cap is inert by default (cap_binding_count 0); with
% the target lowered to 1, every interlocking placed word reaches it.
test(check_target_makes_cap_bind) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered0, _R0, placed),
    crosswordsmith_arrange:cap_binding_count(Numbered0, CB0), CB0 =:= 0,
    setup_call_cleanup(set_check_target(1),
                       ( crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered1, _R1, placed),
                         crosswordsmith_arrange:cap_binding_count(Numbered1, CB1) ),
                       set_check_target(-1)),
    CB1 > 0.

% The carried reward equals the from-scratch oracle on the final (numbered)
% layout (AC-ARR-9): numbering adds only `num`, so layout_reward is unchanged.
test(layout_reward_matches_after_numbering) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, R, Outcome),
    Outcome == placed,
    crosswordsmith_arrange:layout_reward(5, 1, Numbered, R2),
    R =:= R2.

% Regression pin: bundled_17 on grid 17 places all 6, reward 60 (cap inert
% here, so reward == 6*Sum-checked). Update deliberately if construction
% changes (mirror the golden update).
test(bundled_strict_reward, [true(R =:= 60)]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, _Numbered, R, placed).

% cap_binding_count on bundled_17 is 0 (the reachability caveat in miniature).
test(bundled_cap_binding_zero, [true(CB =:= 0)]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, placed),
    crosswordsmith_arrange:cap_binding_count(Numbered, CB).

% --- Phase 2: strict construct + rescore (place-all-or-fail) ------------------

test(strict_places_all_words) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, Outcome),
    Outcome == placed,
    length(Numbered, 6).

% AC-ARR-4: best-of-corners reward is never below a single first-solution.
test(reward_ge_single_corner) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, _N, RBest, placed),
    once(find_crossword(mrv_inc, 17, Words, topleft_across, _G, P0)),
    crosswordsmith_arrange:layout_reward(5, 1, P0, R0),
    RBest >= R0.

% AC-ARR-1 outcome (b): too-small grid -> infeasible, no stdout layout.
test(strict_grid_too_small_infeasible) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 3, Numbered, _R, Outcome),
    Outcome == infeasible,
    Numbered == [].

% AC-ARR-1 outcome (c) + AC-ARR-10 (R4, revamp audit): when the per-corner
% inference budget is exhausted before any corner completes, the outcome is
% `not_proven` - DISTINCT from `infeasible`. A budget too small to finish even
% one corner (1000 inferences << a 6-word grid-17 construction) forces it.
test(strict_budget_exhausted_not_proven) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, 1000, Numbered, Reward, Outcome),
    Outcome == not_proven,
    Numbered == [],
    Reward =:= -1.

% The same input completes (placed) under the real budget - so `not_proven`
% above is the budget biting, not genuine infeasibility (the (b) vs (c) split).
test(strict_budget_not_proven_is_not_infeasible) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, _N, _R, Outcome),
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
    crosswordsmith_arrange:arrange_best_layout(Words, 15, B, _N, _R, Outcome),
    statistics(inferences, I1),
    Used is I1 - I0,
    Outcome == not_proven,
    Used < 2 * B.

% AC-ARR-1 outcome (b): words sharing no letters -> infeasible, named.
test(strict_isolated_words_infeasible) :-
    crosswordsmith_arrange:arrange_best_layout([['ABC'], ['DEF']], 9, _N, _R, Outcome),
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
    crosswordsmith_arrange:construct_one(topleft_across, [['CAT', _{}]], not_an_int, 5, 3, 1_000_000, _Res).

test(unplaceable_names_isolated, [true(Bad == ['ABC', 'DEF'])]) :-
    crosswordsmith_arrange:unplaceable_words([['ABC'], ['DEF']], Bad).

test(unplaceable_empty_when_interlocking, [true(Bad == [])]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:unplaceable_words(Words, Bad).

% R1 (revamp audit): a multi-word/hyphenated answer is placed by its letters
% ONLY - the separator (space or hyphen) is an enumeration marker, not a grid
% cell. WELL-BEING is a 9-cell run (not 10), no placed cell is '-', and the
% original answer atom (with the hyphen) is preserved for export enumeration.
test(hyphenated_answer_strips_separator, [nondet]) :-
    crosswordsmith_arrange:arrange_best_layout([['WELL-BEING'], ['BELOW']], 15, Numbered, _R, placed),
    member(W, Numbered), pw_answer(W, 'WELL-BEING'), !,
    pw_len(W, Len), Len =:= 9,
    pw_cells(W, Cells), length(Cells, 9),
    pw_letters(W, Ls), \+ memberchk('-', Ls), \+ memberchk(' ', Ls).

% --- Phase 3: emit framing ---------------------------------------------------

% fixed: exactly N x N, all words present, valid JSON.
test(emit_fixed_is_n_by_n) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S), crosswordsmith_arrange:emit_arrange(Numbered, Words, 17, fixed)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Rows), length(Rows, 17),
    get_dict(words, Dict, Ws), length(Ws, 6).

% max: cropped square (<= N), content anchored at (0,0).
test(emit_max_shrinks_and_anchors) :-
    arrange_toc(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 25, Numbered, _R, placed),
    with_output_to(string(S), crosswordsmith_arrange:emit_arrange(Numbered, Words, 25, max)),
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
    crosswordsmith_arrange:cropped_coord(17, 2, 3, Cell, Coord).

% AC-ARR-6: identical input -> byte-identical emit (no nondeterminism).
test(emit_is_deterministic) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S1), crosswordsmith_arrange:emit_arrange(Numbered, Words, 17, max)),
    with_output_to(string(S2), crosswordsmith_arrange:emit_arrange(Numbered, Words, 17, max)),
    S1 == S2.

% --- Phase 4: best-effort (drop) ---------------------------------------------

% When the grid fits everything, best-effort places all and drops nothing.
test(best_effort_places_all_when_fits) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_effort(Words, 17, Numbered, _R, NP, Dropped),
    NP =:= 6, Dropped == [], length(Numbered, 6).

% AC-ARR-2: best-effort succeeds where strict fails (isolated words sharing no
% letters), placing a maximal subset and reporting the dropped remainder.
test(best_effort_succeeds_on_isolated) :-
    crosswordsmith_arrange:arrange_best_effort([['ABC'], ['DEF']], 9, Numbered, _R, NP, Dropped),
    NP =:= 1, length(Dropped, 1), length(Numbered, 1).

% On a too-tight grid best-effort drops, and every word is placed XOR dropped.
test(best_effort_drops_and_partitions) :-
    arrange_bundled(Words), length(Words, Total),
    crosswordsmith_arrange:arrange_best_effort(Words, 11, _N, _R, NP, Dropped),
    length(Dropped, ND),
    NP >= 1, NP < Total, NP + ND =:= Total.

% --- Phase 5: fragment-grid seeding (anchors) --------------------------------

% Two long words of bundled_17 pinned (an across + a down crossing at a 'T'),
% used by the partial-seed and best-effort tests below.
bundled_two_word_fragment(Frags) :-
    crosswordsmith_arrange:fragment_dict_words(
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
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:17, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,3],[0,4],[0,5],[0,6]]}]},
        17, Frags).

% AC-EMIT-2 / AC-ARR-3 (load-bearing): an emitted layout IS a valid fragment;
% re-ingesting it and re-emitting reproduces byte-identical JSON (all words
% pinned -> the remainder is empty, so the seed is the whole layout).
test(fragment_roundtrip_byte_identical) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, placed),
    with_output_to(string(S1), crosswordsmith_arrange:emit_arrange(Numbered, Words, 17, fixed)),
    atom_json_dict(S1, Dict, []),
    crosswordsmith_arrange:fragment_dict_words(Dict, GL, Frags),
    GL =:= 17,
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 17, Numbered2, _R2, placed),
    with_output_to(string(S2), crosswordsmith_arrange:emit_arrange(Numbered2, Words, 17, fixed)),
    S1 == S2.

% R3 (revamp audit): the round-trip above passes on bundled_17 only because its
% shared-start pairs happen not to flip. benchmark_08 at size 13 DOES flip (the
% strict DFS and the fragment re-solve build a start-1 across/down pair in
% opposite order), so without the canonical (number, direction) emit ordering
% the re-ingest diverges. This is the finding's exact repro.
test(fragment_roundtrip_byte_identical_reordering) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    crosswordsmith_arrange:arrange_best_layout(Words, 13, Numbered, _R, placed),
    with_output_to(string(S1), crosswordsmith_arrange:emit_arrange(Numbered, Words, 13, fixed)),
    atom_json_dict(S1, Dict, []),
    crosswordsmith_arrange:fragment_dict_words(Dict, GL, Frags),
    GL =:= 13,
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 13, Numbered2, _R2, placed),
    with_output_to(string(S2), crosswordsmith_arrange:emit_arrange(Numbered2, Words, 13, fixed)),
    S1 == S2.

% AC-FRAG-3: pinned words appear at exactly their fragment cells; the engine
% places the remaining four around them (all six end up placed).
test(fragment_pins_preserved_partial) :-
    arrange_bundled(Words),
    bundled_two_word_fragment(Frags),
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 17, Numbered, _R, placed),
    length(Numbered, 6),
    once(( member(Wg, Numbered), pw_answer(Wg, 'GNOSTIC GOSPELS') )),
    pw_cells(Wg, Cg), maplist(cell_coord(17), Cg, RCg),
    RCg == [[0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[7,3],
            [8,3],[9,3],[10,3],[11,3],[12,3],[13,3]],
    once(( member(We, Numbered), pw_answer(We, 'ETERNAL RETURN') )),
    pw_cells(We, Ce), maplist(cell_coord(17), Ce, RCe),
    RCe == [[4,2],[4,3],[4,4],[4,5],[4,6],[4,7],[4,8],
            [4,9],[4,10],[4,11],[4,12],[4,13],[4,14]].

% AC-FRAG-1: a fragment word not present in --input is rejected up front.
test(fragment_word_not_in_input_throws,
     [throws(error(fragment_word_not_in_input('ZZZTOP'), _))]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"ZZZTOP", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3],[0,4],[0,5]]}]},
        9, Frags),
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 9, _, _, _).

% AC-FRAG-2 (overlap): two pinned words demanding different letters in a shared
% cell are reported before any search, naming the clashing cell and letters.
test(fragment_letter_clash_throws,
     [throws(error(fragment_letter_clash('FLOW', [0,0], 'B', 'F'), _))]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3]]},
            _{answer:"FLOW", direction:"down",   cells:[[0,0],[1,0],[2,0],[3,0]]}]},
        9, Frags),
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 9, _, _, _).

% AC-FRAG-2 (self-inconsistent geometry): cells that are not a straight run of
% the answer's length are rejected up front.
test(fragment_cells_inconsistent_throws,
     [throws(error(fragment_cells_inconsistent('BIAS'), _))]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2]]}]},  % 3 cells, 4 letters
        9, Frags),
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 9, _, _, _).

% A fragment may not pin the same answer twice.
test(fragment_duplicate_answer_throws,
     [throws(error(fragment_duplicate_answer('BIAS'), _))]) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,0],[0,1],[0,2],[0,3]]},
            _{answer:"BIAS", direction:"down",   cells:[[0,0],[1,0],[2,0],[3,0]]}]},
        9, Frags),
    crosswordsmith_arrange:arrange_fragment_strict(Words, Frags, 9, _, _, _).

% An off-grid cell is a parse-time error.
test(fragment_invalid_cell_throws,
     [throws(error(fragment_invalid_cell('BIAS', [0,9]), _))]) :-
    crosswordsmith_arrange:fragment_dict_words(
        _{gridLength:9, words:[
            _{answer:"BIAS", direction:"across", cells:[[0,9],[0,10],[0,11],[0,12]]}]},
        _, _).

% Best-effort with a fragment: the seed is pinned, the rest greedily placed;
% on a roomy grid nothing is dropped and the pins are kept.
test(fragment_best_effort_places_all) :-
    arrange_bundled(Words),
    bundled_two_word_fragment(Frags),
    crosswordsmith_arrange:arrange_fragment_best_effort(Words, Frags, 17, Numbered, _R, NP, Dropped),
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
    crosswordsmith_arrange:arrange_candidates(Words, 17, strict, 3, Layouts, Returned),
    Returned =:= 3,
    length(Layouts, 3),
    forall(member(L, Layouts), length(L, 6)).

% AC-ARR-7: the returned layouts are pairwise >= tau apart (translation-invariant
% placement distance).
test(candidates_pairwise_tau_apart) :-
    arrange_bundled(Words), length(Words, Total),
    crosswordsmith_arrange:candidate_tau_pct(TauPct),
    crosswordsmith_arrange:arrange_candidates(Words, 17, strict, 3, Layouts, _),
    % note the parens: (:)/2 is priority 600, ABOVE (>>)/2's 400, so an
    % unparenthesized qualified body would parse as (Lambda:Goal), not a lambda
    maplist([L, A]>>(crosswordsmith_arrange:placement_assoc(L, 17, A)), Layouts, Assocs),
    forall( ( nth0(I, Assocs, A1), nth0(J, Assocs, A2), I < J ),
            ( crosswordsmith_arrange:pos_diff_count(A1, A2, Diff), Diff * 100 >= TauPct * Total ) ).

% AC-ARR-7: fewer than K only when fewer >= tau-distinct layouts exist; on the
% 8-word benchmark the greedy breadth yields just two distinct full layouts.
test(candidates_short_return_when_fewer_distinct) :-
    load_clues('fixtures/benchmark_08_words.pl', Words),
    crosswordsmith_arrange:arrange_candidates(Words, 13, strict, 5, Layouts, Returned),
    Returned < 5, Returned >= 1, length(Layouts, Returned).

% --candidates 1 yields exactly the single best layout.
test(candidates_one_is_single_best) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_candidates(Words, 17, strict, 1, Layouts, 1),
    Layouts = [L], length(L, 6).

% AC-ARR-6 / INV-2: candidate emission is byte-identical across runs.
test(candidates_emit_deterministic) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_candidates(Words, 17, strict, 3, Layouts, _),
    with_output_to(string(S1), crosswordsmith_arrange:emit_candidates(Layouts, Words, 17, fixed)),
    with_output_to(string(S2), crosswordsmith_arrange:emit_candidates(Layouts, Words, 17, fixed)),
    S1 == S2.

% The emitted candidates form a JSON array, each element a standalone layout.
test(candidates_emit_is_json_array) :-
    arrange_bundled(Words),
    crosswordsmith_arrange:arrange_candidates(Words, 17, strict, 3, Layouts, _),
    with_output_to(string(S), crosswordsmith_arrange:emit_candidates(Layouts, Words, 17, fixed)),
    atom_json_dict(S, Arr, []),
    is_list(Arr), length(Arr, 3),
    forall(member(D, Arr),
           ( get_dict(gridLength, D, 17), get_dict(words, D, Ws), length(Ws, 6) )).

% Placement distance is translation-invariant: a layout and its one-row shift
% have distance 0 (the same crossword, just moved).
test(candidates_distance_translation_invariant) :-
    L1 = [ pw('ABC', _, [1,2,3],    across, _, 1,  _, _),
           pw('ADE', _, [1,18,35],  down,   _, 1,  _, _) ],
    L2 = [ pw('ABC', _, [18,19,20], across, _, 18, _, _),
           pw('ADE', _, [18,35,52], down,   _, 18, _, _) ],
    crosswordsmith_arrange:placement_assoc(L1, 17, A1),
    crosswordsmith_arrange:placement_assoc(L2, 17, A2),
    crosswordsmith_arrange:pos_diff_count(A1, A2, 0).

% --- Phase 7: enumerate seam -------------------------------------------------

% AC-ARR-8: --enumerate counts every feasible full placement, matching the old
% --all count (it is all_crossword/5 under the production default strategy, over
% all start corners).
test(enumerate_matches_all_crossword) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    crosswordsmith_arrange:arrange_enumerate(Words, 17, N1),
    default_strategy(Strat),
    all_crossword(Strat, 17, Words, _StartLoc, N2),
    N1 =:= N2,
    N1 > 0.

% --- greedy constructor (moved here from metrics.pl in Phase 3) ---------------

% seed_candidates/2 returns the longest words first (restart diversity). With
% fewer words than the cap K, all are returned, longest-first.
test(seed_candidates_longest_first) :-
    crosswordsmith_arrange:seed_candidates([['AAAA',_{}],['BB',_{}],['CCC',_{}]], Seeds),
    Seeds = [[S1|_], [S2|_], [S3|_]],
    S1 == 'AAAA', S2 == 'CCC', S3 == 'BB'.

% --- --seed: opt-in pseudo-random search perturbation -----------------------
% Design contract: a fully deterministic default path alongside a reproducible
% seeded path, with the seed NEVER on the deterministic flow. set_search_seed/1
% (exported by core) is the sole writer; -1 clears it. The byte-exact
% deterministic emit is pinned by the golden (tests/golden/arrange_*.json);
% these cover reproducibility, quality preservation, and — critically — that a
% cleared seed restores the deterministic result even after the RNG was used.

% The PRNG is module-OWNED (portable splitmix64), so its output sequence is
% part of the reproducibility contract: same seed => same layout on EVERY
% build (native GMP, wasm LibBF — the VM RNG diverges between those backends,
% which is why we own the algorithm). Known-answer locks: the seed-0 values
% are the published splitmix64 test vectors; a refactor that shifts any draw
% silently re-maps every kept seed, and only these tests would notice.
test(splitmix64_known_answer) :-
    crosswordsmith_core:splitmix64(0, V1, S1),
    V1 =:= 16294208416658607535,
    crosswordsmith_core:splitmix64(S1, V2, S2),
    V2 =:= 7960286522194355700,
    crosswordsmith_core:splitmix64(S2, V3, _),
    V3 =:= 487617019471545679.

% ...and the full draw discipline (destructive state advance + V mod N
% selection shuffle), locked end-to-end through set_search_seed/1. Reference
% permutation computed with an independent Python implementation.
test(seeded_permutation_known_answer,
     [ setup(set_search_seed(42)), cleanup(set_search_seed(-1)) ]) :-
    crosswordsmith_core:seeded_permutation([a,b,c,d,e,f,g], P),
    P == [f,b,e,a,d,c,g].

% Placement fingerprint: the sorted (answer, start, dir) triples of a layout.
layout_sig(Numbered, Sig) :-
    findall(A-S-D,
            ( member(W, Numbered),
              pw_answer(W, A), pw_start(W, S), pw_dir(W, D) ),
            Sig0),
    sort(Sig0, Sig).

% Same seed => identical layout (reproducible). The seed is re-installed before
% the second run so the PRNG stream (module-owned splitmix64) restarts from the
% same point.
test(seed_reproducible) :-
    arrange_bundled(Words),
    setup_call_cleanup(
        set_search_seed(7),
        ( crosswordsmith_arrange:arrange_best_layout(Words, 17, N1, _, placed),
          set_search_seed(7),
          crosswordsmith_arrange:arrange_best_layout(Words, 17, N2, _, placed) ),
        set_search_seed(-1)),
    layout_sig(N1, S1), layout_sig(N2, S2),
    S1 == S2.

% A seeded search preserves completeness/quality: still a full placement.
test(seed_places_all_words) :-
    arrange_bundled(Words),
    setup_call_cleanup(
        set_search_seed(1),
        crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _R, Outcome),
        set_search_seed(-1)),
    Outcome == placed,
    length(Numbered, 6).

% At least one seed actually perturbs the layout away from the deterministic one
% (robust to RNG-impl detail: it suffices that some seed in the set differs).
test(seed_can_perturb_layout) :-
    arrange_bundled(Words),
    set_search_seed(-1),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, DetN, _, placed),
    layout_sig(DetN, DetSig),
    once(( member(Seed, [0, 1, 2, 3, 7, 42]),
           setup_call_cleanup(
               set_search_seed(Seed),
               crosswordsmith_arrange:arrange_best_layout(Words, 17, SeedN, _, placed),
               set_search_seed(-1)),
           layout_sig(SeedN, SeedSig),
           SeedSig \== DetSig )).

% The seed is OFF the deterministic path: after running a seeded search (which
% advances the PRNG state) and then clearing the seed, the deterministic layout
% is byte-for-byte (fingerprint) identical to the pre-seed run.
test(seed_cleared_restores_deterministic) :-
    arrange_bundled(Words),
    set_search_seed(-1),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Before, _, placed),
    layout_sig(Before, SigBefore),
    setup_call_cleanup(
        set_search_seed(5),
        crosswordsmith_arrange:arrange_best_layout(Words, 17, _, _, placed),
        set_search_seed(-1)),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, After, _, placed),
    layout_sig(After, SigAfter),
    SigAfter == SigBefore.

% --- C1: in-process inference counts are history-independent ------------------
% reset_search_memos/0 (core.pl) runs ONCE at every top-level search entry, so
% a search's inference count cannot inherit table work from whatever ran before
% it in the same process. Pre-fix, only the strict path abolished: a greedy
% run left warm pair_crossings/3 memos and later counts drifted by THOUSANDS
% of inferences (audit C1: 11,454 cold vs 4,237 warm on this very fixture).
% Locked here with statistics(inferences) deltas. Two irrelevant residuals are
% probe-documented (see reset_search_memos/0): first calls carry one-time
% JIT/autoload cost (absorbed by the warmup pair), and the reset's flush of
% the PREVIOUS entry's residue costs +-2 inferences depending on that
% residue's shape - so identical (entry, predecessor-type) pairs must agree
% EXACTLY, and cross-history comparisons must agree within 4 inferences
% (~0.002%; the pre-fix defect was ~3 orders of magnitude larger).
test(inference_counts_history_independent) :-
    arrange_bundled(Words),
    % warmup: one-time JIT/autoload, discarded (fresh output vars per call -
    % a reused goal term would carry the first run's bindings forward)
    greedy_inf(Words, _),
    strict_inf(Words, _),
    % measured sequence (predecessor type annotated)
    greedy_inf(Words, G1),   % greedy after strict
    strict_inf(Words, S1),   % strict after greedy
    greedy_inf(Words, G2),   % greedy after strict
    strict_inf(Words, S2),   % strict after greedy
    strict_inf(Words, S3),   % strict after strict ("cold": no greedy warmth)
    strict_inf(Words, S4),   % strict after strict
    greedy_inf(Words, G3),   % greedy after strict
    greedy_inf(Words, G4),   % greedy after greedy (pre-fix: warm memo, cheap)
    greedy_inf(Words, G5),   % greedy after greedy
    % identical in-process solves yield IDENTICAL deltas
    assertion(G1 =:= G2),
    assertion(G1 =:= G3),
    assertion(S1 =:= S2),
    assertion(S3 =:= S4),
    assertion(G4 =:= G5),
    % cross-history: strict-after-greedy == strict-cold, greedy-after-greedy ==
    % greedy-after-strict, up to the +-2 residue-flush walk
    assertion(abs(S1 - S3) =< 4),
    assertion(abs(G4 - G1) =< 4).

% --shuffle draws a fresh seed from OS entropy but stays RECOVERABLE: it picks a
% concrete integer N, and re-running with --seed N reproduces the same layout.
% (We assert recoverability, not that two shuffles differ — that is inherently
% probabilistic, and on a tiny word set only a few distinct layouts exist.)
test(shuffle_seed_is_recoverable) :-
    arrange_bundled(Words),
    setup_call_cleanup(
        set_shuffle_seed(N),
        crosswordsmith_arrange:arrange_best_layout(Words, 17, ShufN, _, placed),
        set_search_seed(-1)),
    integer(N), N >= 0,
    setup_call_cleanup(
        set_search_seed(N),
        crosswordsmith_arrange:arrange_best_layout(Words, 17, RepeatN, _, placed),
        set_search_seed(-1)),
    layout_sig(ShufN, S1), layout_sig(RepeatN, S2),
    S1 == S2.

% Provenance: a perturbed layout records its seed in diagnostics.arrange.seed,
% so the JSON artifact is self-documenting (reproduce with --seed N).
test(seed_recorded_in_diagnostics, [true(S == 7)]) :-
    arrange_bundled(Words),
    setup_call_cleanup(
        set_search_seed(7),
        ( crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _, placed),
          crosswordsmith_arrange:arrange_diag_dict(Numbered, Words, Diag) ),
        set_search_seed(-1)),
    get_dict(seed, Diag, S).

% The deterministic default omits the seed key entirely (missing = default),
% so its emitted diagnostics — and thus the golden — are byte-unchanged.
test(no_seed_absent_from_diagnostics) :-
    arrange_bundled(Words),
    set_search_seed(-1),
    crosswordsmith_arrange:arrange_best_layout(Words, 17, Numbered, _, placed),
    crosswordsmith_arrange:arrange_diag_dict(Numbered, Words, Diag),
    \+ get_dict(seed, Diag, _).

:- end_tests(arrange).
