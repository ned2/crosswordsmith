% tests/crossword.plt - plunit test suite for crossword.pl
%
% These tests assume crossword.pl has already been consulted into the
% `user` module (the runner, tests/run_tests.pl, does this before loading
% this file). Run them via:
%
%     ./run_tests.sh        (or)        make test
%
% The tests are grouped into units: pure grid geometry, generic list
% utilities, the solver, and clue numbering.

:- use_module(library(plunit)).


% Grid geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cells are numbered 1..GridLen*GridLen in row-major order. `across` steps
% by +1 along a row; `down` steps by +GridLen down a column.

:- begin_tests(geometry).

test(swap_dir_across) :- swap_dir(across, down).
test(swap_dir_down)   :- swap_dir(down, across).

% A word fits across only if there is room to its right in the row.
test(fits_across_full_row)      :- fits_on_grid(across, 1, 17, 17).
test(fits_across_too_long, [fail]) :- fits_on_grid(across, 1, 18, 17).
% Cell 10 is column 10 of 17, so only 8 cells remain (10..17).
test(fits_across_mid_row)       :- fits_on_grid(across, 10, 8, 17).
test(fits_across_overflow, [fail]) :- fits_on_grid(across, 10, 9, 17).
% A cell at the end of a row (mod == 0) cannot start an across word.
test(fits_across_row_end, [fail]) :- fits_on_grid(across, 17, 2, 17).

% A word fits down only if it does not run off the bottom of the grid.
test(fits_down_full_col)        :- fits_on_grid(down, 1, 17, 17).
test(fits_down_too_long, [fail]) :- fits_on_grid(down, 1, 18, 17).

% calc_num and calc_start are inverses: given a word starting at WStart,
% the letter at position P has cell number N; calc_start must recover
% WStart from N and P.
test(calc_num_across, [true(N =:= 7)])  :- calc_num(across, 17, 3, 5, N).
test(calc_num_down,   [true(N =:= 35)]) :- calc_num(down, 17, 3, 1, N).
test(calc_roundtrip_across, [true(S =:= 5)]) :-
    calc_num(across, 17, 3, 5, N), calc_start(across, 17, 3, N, S).
test(calc_roundtrip_down, [true(S =:= 1)]) :-
    calc_num(down, 17, 4, 1, N), calc_start(down, 17, 4, N, S).

test(next_cell_across, [true(C =:= 6)])  :- next_cell(across, 5, 17, C).
test(next_cell_down,   [true(C =:= 22)]) :- next_cell(down, 5, 17, C).
test(prev_cell_across, [true(C =:= 4)])  :- prev_cell(across, 5, 17, C).
test(prev_cell_down,   [true(C =:= 5)])  :- prev_cell(down, 22, 17, C).

% First column starts across rows; first row starts down columns.
test(is_start_across)            :- is_start_cell(across, 18, 17).   % col 1, row 2
test(is_start_across_no, [fail]) :- is_start_cell(across, 19, 17).
test(is_start_down)              :- is_start_cell(down, 5, 17).      % row 1
test(is_start_down_no, [fail])   :- is_start_cell(down, 18, 17).
test(is_end_across)              :- is_end_cell(across, 17, 17).     % last col
test(is_end_down)                :- is_end_cell(down, 289, 17).      % last row

% The four named start locations resolve to the expected cell + direction.
test(start_topleft_across, [true(N-D == 1-across)])  :- start_loc(topleft_across, 17, N, D).
test(start_topleft_down,   [true(N-D == 1-down)])    :- start_loc(topleft_down, 17, N, D).
test(start_topright,       [true(N-D == 17-down)])   :- start_loc(topright, 17, N, D).
test(start_bottomleft,     [true(N-D == 273-across)]):- start_loc(bottomleft, 17, N, D).

:- end_tests(geometry).


% Generic list utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(utilities).

% position/3 yields every 1-based index of X, on backtracking.
test(position_first, [true(P == 1)])  :- once(position(a, [a,b,c], P)).
test(position_all,   [all(P == [1,3])]) :- position(a, [a,b,a], P).
test(position_none,  [fail])          :- position(z, [a,b,c], _).

% remove_x removes only the first occurrence. (nondet: it leaves a harmless
% choicepoint; we only care about the first answer.)
test(remove_x_first, [true(R == [a,c,b]), nondet]) :- remove_x(b, [a,b,c,b], R).
test(remove_x_absent, [true(R == [a,b,c]), nondet]) :- remove_x(z, [a,b,c], R).

% shuffle is a permutation of its input (we cannot assert order: it is
% random, so just check the multiset of elements is preserved).
test(shuffle_is_permutation, [true(Sorted == [1,2,3,4,5])]) :-
    shuffle([1,2,3,4,5], Out), msort(Out, Sorted).

% init_grid builds a GridLen*GridLen assoc of `empty` cells.
test(init_grid_size, [true(Len =:= 9)]) :-
    init_grid(3, G), assoc_to_list(G, L), length(L, Len).
test(init_grid_empty) :-
    init_grid(3, G), get_assoc(1, G, empty), get_assoc(9, G, empty).
test(init_grid_no_extra, [fail]) :-
    init_grid(3, G), get_assoc(10, G, _).

:- end_tests(utilities).


% The solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(solver).

% Two words that cross on a shared letter are both placed. OMEGA POINT
% (across, starting cell 1) and GNOSTIC GOSPELS (down) cross on G.
% (nondet: find_crossword backtracks over many solutions; we assert one exists.)
test(places_two_crossing_words, [nondet]) :-
    Words = [['OMEGA POINT', 'c1', ''], ['GNOSTIC GOSPELS', 'c2', '']],
    find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 2),
    % the first word is laid across starting at cell 1
    memberchk([_,_,_,_,_,across,_,1,_], Placed).

% The full bundled clue set has a solution on a 17x17 grid.
test(bundled_clues_solve_at_17, [nondet]) :-
    clues(Words),
    find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 6).

% A grid that is too small for the words has no solution.
test(too_small_grid_fails, [fail]) :-
    clues(Words),
    find_crossword(3, Words, topleft_across, _Grid, _Placed).

% End-to-end: the top-level crossword/3 (which also numbers clues and
% prints) succeeds on the bundled clues. Output is captured to keep test
% logs clean.
test(full_pipeline_bundled) :-
    clues(Words),
    with_output_to(string(_), crossword(17, Words, topleft_across)).

:- end_tests(solver).


% Clue numbering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A placed word is [Word, Clue, Link, Letters, Cells, Dir, WLen, Start, ClueNum].
% assign_clue_numbers/2 sorts by Start and fills in ClueNum.

:- begin_tests(clue_numbering).

% Helper: clue number assigned to the word with the given text.
clue_num_of(Word, Placed, Num) :-
    member([Word,_,_,_,_,_,_,_,Num], Placed).

% Two words with distinct start cells get sequential numbers, in start order.
test(distinct_starts_numbered_in_order, [nondet]) :-
    W1 = ['CAT', 'c', '', [c,a,t], [1,2,3],   across, 3, 1, _],
    W2 = ['DOG', 'c', '', [d,o,g], [5,12,19], down,   3, 5, _],
    assign_clue_numbers([W2,W1], Placed),   % deliberately unsorted input
    clue_num_of('CAT', Placed, 1),
    clue_num_of('DOG', Placed, 2).

% Regression test for the (now-fixed) add_clue_nums/3 arity bug: when one
% cell is the start of BOTH an across and a down word, the two clues must
% share a clue number. Before the fix the buggy clause had the wrong arity
% and this scenario failed clue numbering entirely.
test(shared_start_cell_shares_number, [nondet]) :-
    Wa = ['CAT', 'c', '', [c,a,t], [1,2,3],   across, 3, 1, _],
    Wd = ['COW', 'c', '', [c,o,w], [1,18,35], down,   3, 1, _],
    assign_clue_numbers([Wa,Wd], Placed),
    clue_num_of('CAT', Placed, N),
    clue_num_of('COW', Placed, N),
    integer(N).

:- end_tests(clue_numbering).
