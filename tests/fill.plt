% tests/fill.plt - plunit suite for fill.pl (grid-first auto-fill, §8.4).
%
% Assumes core.pl, metrics.pl, lint.pl, stockgrid.pl, arrange.pl and
% fill.pl are consulted by the runner before this file. Covers slot derivation
% + the SHARED cell-variable invariant (the crossing constraint), the fill
% search (AC-FILL-1), seeds (AC-FILL-2), determinism (AC-FILL-3), and the
% no-fill report. The byte-exact fill is pinned by tests/golden/fill_3.json.

:- use_module(library(plunit)).

% Run a fill and return the across/down answer lists (deterministic).
do_fill(GridFile, SeedFile, DictFile, Across, Down) :-
    fill_grid(GridFile, _Size, Slots, _CellVar),
    ( SeedFile == none -> SeededKeys = [] ; apply_seeds(SeedFile, Slots, SeededKeys) ),
    exclude(seeded_slot(SeededKeys), Slots, SearchSlots),
    load_dict(DictFile, DictByLen, Index),
    fill_attempt(SearchSlots, Slots, DictByLen, Index, filled, Numbered, _),
    findall(A, ( member(W, Numbered), get_dict(dir, W, across), get_dict(answer, W, A) ), Across),
    findall(A, ( member(W, Numbered), get_dict(dir, W, down),   get_dict(answer, W, A) ), Down).


:- begin_tests(fill).

% --- slots + the shared-variable invariant (the crossing constraint) ---------
% The 3x3 has 3 across + 3 down slots.
test(fill_derives_slots) :-
    fill_grid('fixtures/fill_grid_3.json', 3, Slots, _),
    length(Slots, 6).

% REGRESSION: a crossing cell's variable must be SHARED between its across and
% down slot (cell 1 starts both row-0-across and col-0-down). A findall-copied
% or yall-lambda-copied template breaks this and the fill is inconsistent.
test(fill_crossing_cells_are_shared) :-
    fill_grid('fixtures/fill_grid_3.json', 3, Slots, _),
    once(( member(slot(1, across, [1,2,3], [A1|_]), Slots) )),
    once(( member(slot(1, down,   [1,4,7], [D1|_]), Slots) )),
    A1 == D1.

% --- the fill (AC-FILL-1) ----------------------------------------------------
% A complete, consistent fill: every across is a row, every down is a column
% (deterministic - pinned to the exact result).
test(fill_produces_consistent_square) :-
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', Across, Down),
    Across == ['CAT', 'ORE', 'WED'],
    Down   == ['COW', 'ARE', 'TED'].

% No dictionary word of the slot length -> infeasible (reported, not silent).
test(fill_infeasible_when_no_matching_words) :-
    tmp_file_stream(text, F, S), write(S, "FOUR\nFIVE\n"), close(S),   % 4-letter only
    fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    load_dict(F, DictByLen, Index),
    fill_attempt(Slots, Slots, DictByLen, Index, Outcome, _, _),
    delete_file(F),
    Outcome == infeasible.

% Regression for P2: fill_attempt/8 must NOT swallow a genuine error from
% fill_search and report it as `infeasible`. call_with_inference_limit/3 handles
% the budget itself and re-throws real exceptions; infeasibility is a search
% FAILURE (R == exhausted), never a throw. Feeding a malformed pattern index (a
% non-assoc) makes get_assoc/3 raise a type_error, which must propagate rather
% than be masked as infeasibility. The old broad catch/3 returned `infeasible`.
test(fill_propagates_genuine_error,
     [throws(error(type_error(btree, _), _))]) :-
    fill_attempt([slot(a, across, 3, [_,_,_])], [slot(a, across, 3, [_,_,_])],
                 [], not_an_assoc, 1_000_000, _Outcome, _Numbered, _InputWords).

% --- seeds (AC-FILL-2) -------------------------------------------------------
% Pinning COW across row 0 forces the transpose square; the pin appears at its
% slot, and the unseeded default (CAT across) is overridden.
test(fill_respects_seed) :-
    do_fill('fixtures/fill_grid_3.json', 'fixtures/fill_seed_3.json',
            'fixtures/wordlist_sample.txt', Across, _Down),
    Across == ['COW', 'ARE', 'TED'],
    Across \== ['CAT', 'ORE', 'WED'].

% R2 (revamp audit): a seed is a HARD PIN, not required to be a dictionary word
% (a setter's theme/own word). With CAT pinned across row 0 and a dictionary
% that OMITS CAT, the grid still fills around the pin (CAT/ORE/WED across) - the
% old code reported the user's own pinned slot as unfillable and failed.
test(fill_seed_need_not_be_in_dict) :-
    tmp_file_stream(text, F, S),
    write(S, "ORE\nWED\nCOW\nARE\nTED\n"), close(S),   % the sample wordlist minus CAT
    fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    foldl(apply_seed(Slots), [frag('CAT', across, 1, [1, 2, 3])], [], SeededKeys),
    exclude(seeded_slot(SeededKeys), Slots, SearchSlots),
    load_dict(F, DictByLen, Index),
    fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, _),
    delete_file(F),
    Outcome == filled,
    findall(A, ( member(W, Numbered), get_dict(dir, W, across), get_dict(answer, W, A) ), Across),
    Across == ['CAT', 'ORE', 'WED'].

% A seed whose answer matches no slot of the grid is rejected.
test(fill_seed_no_slot_throws, [throws(error(fill_seed_no_slot('CAT'), _))]) :-
    fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    % CAT down at col 0 cells [1,4,7] is a real slot; ask for it as a 2-cell
    % run that no slot has -> no match.
    apply_seeds_frags(Slots, [frag('CAT', across, 2, [2, 3])]).

% R11 (revamp audit): AC-FILL-1 "not proven within budget". A budget too small
% to complete the search yields Outcome==not_proven - distinct from infeasible
% (a genuinely 0-candidate slot) and from filled.
test(fill_budget_exhausted_not_proven) :-
    fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    load_dict('fixtures/wordlist_sample.txt', DictByLen, Index),
    fill_attempt(Slots, Slots, DictByLen, Index, 100, Outcome, Numbered, _),
    Outcome == not_proven,
    Numbered == [].

% R11 (revamp audit): AC-FILL-4 - a produced fill re-validates as a legal
% layout. Linting the fill output gives a clean structural verdict: the 3x3
% full grid is every-cell-checked, connected, and symmetric, so it PASSes the
% strictest structural profile (american). Letters do not affect structure.
test(fill_revalidates_under_lint) :-
    lint_load('tests/golden/fill_3.json', GridLen, Placed),
    lint_run(Placed, GridLen, american, false, Report),
    get_dict(verdict, Report, V),
    V == 'PASS'.

% R6 (revamp audit): when a start cell begins both an across and a down slot,
% select_mrv must expand the slot that actually has the fewest candidates, not
% whichever shares the start and was generated first. Here start 1's down slot
% is fully constrained to CAT (1 candidate) while its across slot has 3 (C__);
% the down slot (the minimum) must be chosen.
test(select_mrv_recovers_correct_direction_on_tie) :-
    tmp_file_stream(text, F, S), write(S, "CAT\nCOW\nCUB\nDOG\n"), close(S),
    load_dict(F, DictByLen, Index), delete_file(F),
    Across = slot(1, across, [1, 2, 3], ['C', _, _]),       % C__ -> CAT/COW/CUB (3)
    Down   = slot(1, down,   [1, 4, 7], ['C', 'A', 'T']),   % CAT (1, the minimum)
    select_mrv([Across, Down], DictByLen, Index, Best, _Rest, Cands),
    Best = slot(1, down, _, _),
    Cands == [['C', 'A', 'T']].

% P13: select_mrv/6 must succeed DETERMINISTICALLY - the winning slot's Start+Dir
% is unique, so once(select/3) prunes the spurious choicepoint plain select/3
% would otherwise leave. The winner (Down, 1 candidate) is placed FIRST in the
% slot list so a stray choicepoint from select/3 (scanning the rest of the list)
% is possible - that is exactly what once/1 removes.
test(select_mrv_leaves_no_choicepoint) :-
    tmp_file_stream(text, F, S), write(S, "CAT\nCOW\nCUB\nDOG\n"), close(S),
    load_dict(F, DictByLen, Index), delete_file(F),
    Across = slot(1, across, [1, 2, 3], ['C', _, _]),
    Down   = slot(1, down,   [1, 4, 7], ['C', 'A', 'T']),
    select_mrv([Down, Across], DictByLen, Index, Best, _Rest, _Cands),
    Best = slot(1, down, _, _),
    deterministic(Det),
    Det == true.

% P3: candidate_count/4 (the count-only MRV metric, no word materialization) must
% agree with length(candidates/4) on BOTH branches - `all` (no cell bound) and
% idx (some cells bound). select_mrv orders slots by these counts, so if the two
% ever diverged the MRV choice would silently differ from the true candidate set.
test(candidate_count_matches_candidates) :-
    load_dict('fixtures/wordlist_sample.txt', DictByLen, Index),
    % all-unbound (`all` branch): count == number of length-3 words materialized
    length(Free, 3),
    candidates(Free, DictByLen, Index, CAll), length(CAll, NAll),
    candidate_count(Free, DictByLen, Index, NAll),
    NAll > 0,
    % one cell bound to 'C' (idx branch): only CAT, COW match
    Bound = ['C', _, _],
    candidates(Bound, DictByLen, Index, CB), length(CB, NB),
    candidate_count(Bound, DictByLen, Index, NB),
    NB =:= 2.

% --- determinism (AC-FILL-3) -------------------------------------------------
test(fill_deterministic) :-
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A1, D1),
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A2, D2),
    A1 == A2, D1 == D2.

:- end_tests(fill).

% Helper: apply a list of frag/4 seeds directly (mirrors apply_seeds/3 minus the
% file read), for the no-slot test; the seeded-key accumulator is discarded.
apply_seeds_frags(Slots, Frags) :- foldl(apply_seed(Slots), Frags, [], _).
