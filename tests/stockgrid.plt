% tests/stockgrid.plt - plunit suite for stockgrid.pl (bundled stock-grid
% library, §8.3; schema fixed by DP-1 / OD-5).
%
% Assumes core.pl, metrics.pl, lint.pl and stockgrid.pl are consulted by
% the runner before this file. Covers the mask loader, light derivation, the
% legality validator (which reuses lint blocked-uk), and - the load-bearing
% regression - that every bundled grid validates as legal.

:- use_module(library(plunit)).

% The shipped library (OD-6): a small curated, lint-validated starter set.
bundled_grid('grids/blocked_13a.json').
bundled_grid('grids/blocked_13b.json').
bundled_grid('grids/blocked_15a.json').


:- begin_tests(stockgrid).

% --- loader ------------------------------------------------------------------
test(stockgrid_loads_mask) :-
    stockgrid_load('grids/blocked_13a.json', grid(Name, Size, Mask)),
    Name == 'blocked-13a', Size == 13, length(Mask, 13),
    once(( member(Row, Mask), length(Row, 13) )).

% --- light derivation --------------------------------------------------------
% A 3x3 all-white grid has 3 across + 3 down lights, each length 3.
test(stockgrid_derives_lights) :-
    maplist(string_chars, ["...", "...", "..."], Mask),
    mask_white_cells(Mask, 3, WS),
    crosswordsmith_stockgrid:grid_lights(3, WS, Words),
    length(Words, 6),
    forall(member(W, Words), get_dict(len, W, 3)).

% A '#' splits a row into lights; runs shorter than 2 are not lights.
test(stockgrid_splits_on_blocks) :-
    maplist(string_chars, ["...#...", ".......", "...#...",
                           ".......", "...#...", ".......", "...#..."], Mask),
    mask_white_cells(Mask, 7, WS),
    crosswordsmith_stockgrid:grid_lights(7, WS, Words),
    once(( member(W, Words), get_dict(dir, W, across), get_dict(len, W, 3) )).

% --- validator catches illegal grids -----------------------------------------
% A grid whose interior columns carry no down light leaves a 3-long unchecked
% run + sub-half checking in the middle rows -> FAIL.
test(stockgrid_detects_illegal) :-
    maplist(string_chars, [".....", ".###.", ".....", ".###.", "....."], Mask),
    stockgrid_validate(grid('test-5', 5, Mask), Verdict, Fails),
    Verdict == 'FAIL',
    once(( member(fail(R, _, _), Fails), memberchk(R, [checked_half, max_unch_run]) )).

% An isolated white cell (in no light) is reported. Centre cell (2,2) = cell 13
% on a 5x5, ringed by blocks at its four orthogonal neighbours.
test(stockgrid_detects_isolated_cell) :-
    maplist(string_chars, [".....", "..#..", ".#.#.", "..#..", "....."], Mask),
    stockgrid_validate(grid('iso-5', 5, Mask), Verdict, Fails),
    Verdict == 'FAIL',
    once(( member(fail(isolated_cells, _, Cells), Fails), memberchk(13, Cells) )).

% --- the load-bearing regression: every bundled grid is legal ----------------
test(bundled_grids_all_legal) :-
    forall( bundled_grid(F),
            ( stockgrid_validate_file(F, V, Fails),
              ( V \== 'FAIL' -> true ; throw(grid_failed(F, V, Fails)) ) )).

% Each bundled grid is 180-degree symmetric (the symmetry rule passes).
test(bundled_grids_are_symmetric) :-
    forall( bundled_grid(F),
            ( stockgrid_load(F, grid(_, Size, Mask)),
              mask_white_cells(Mask, Size, WS),
              crosswordsmith_stockgrid:grid_lights(Size, WS, Words),
              lint_run(Words, Size, 'blocked-uk', false, R),
              get_dict(grid, R, G),
              once(( member(GR, G), get_dict(rule, GR, symmetry), get_dict(severity, GR, 'PASS') )) )).

:- end_tests(stockgrid).
