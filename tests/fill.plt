% tests/fill.plt - plunit suite for fill.pl (grid-first auto-fill, §8.4).
%
% Assumes crossword.pl, quality.pl, lint.pl, stockgrid.pl, arrange.pl and
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

% --- determinism (AC-FILL-3) -------------------------------------------------
test(fill_deterministic) :-
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A1, D1),
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A2, D2),
    A1 == A2, D1 == D2.

:- end_tests(fill).

% Helper: apply a list of frag/4 seeds directly (mirrors apply_seeds/3 minus the
% file read), for the no-slot test; the seeded-key accumulator is discarded.
apply_seeds_frags(Slots, Frags) :- foldl(apply_seed(Slots), Frags, [], _).
