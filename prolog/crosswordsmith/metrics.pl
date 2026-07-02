% metrics.pl - the shared layout metric predicates (design-spec §6.4).
% (Formerly quality.pl; renamed in Phase 2 of the source-structure migration.
% The greedy density constructor that lived here moved to arrange.pl — its
% only consumer — in Phase 3.)
%
% Pure measurement over a placed layout (dir_cells, checked_cells,
% word_checked_count, word_meets_half, word_max_unch_run, placed_bbox,
% word_cells, ...), consulted AFTER core.pl (next_cell/4 is the only core
% predicate used, by word_cells/5). Consumed by arrange.pl as optimizer
% signals and by `lint` as validators. Per spec §4 these stay a separate
% module — lint depends only on the JSON contract plus this metric layer,
% never on the solver substrate.
%
% The former `--quality` CLI engine (quality_solve / quality_layout /
% grid_candidates + the floor subsystem) was retired at the Phase-7 CLI cutover
% (design-spec §4); arrange.pl supersedes it.
%
% Exports: only the predicates with real cross-module consumers (lint,
% arrange, fill — see the migration plan §4.5). checked_cells/2,
% word_meets_half/2, word_max_unch_run/3 and dir_cells/3 are test-only and
% deliberately NOT exported — tests reach them as
% crosswordsmith_metrics:Pred(...).

:- module(crosswordsmith_metrics,
          [ % lint's validator primitives
            layout_dir_cells/2,
            word_checked_bitmap/3,
            bits_checked_count/2,
            bits_max_unch_run/2,
            word_half_threshold/2,
            cell_rc/4,
            % arrange's optimizer signals (word_letters/3 also: fill)
            word_checked_count/3,
            word_letters/3,
            placed_bbox/4,
            word_cells/5
          ]).

:- use_module(library(ordsets)).
:- use_module(library(apply)).
% pairs_keys_values/3 lives here; imported explicitly (not via autoload) to
% match fill.pl's style and survive a qsave_program(..., [autoload(false)])
% (P11).
:- use_module(library(pairs)).

% next_cell/4 (word_cells/5's stepper) — metrics' ONLY core dependency,
% confirmed at the Phase-3 constructor extraction.
:- use_module(crosswordsmith(core), [next_cell/4]).

% --- shared answer/footprint helpers ----------------------------------------

% The "placement footprint" of an answer: its letters with word separators
% (spaces AND hyphens) removed - those are enumeration markers, not grid cells.
% The original answer atom is carried to emit (placed_to_word/4) so export still
% derives the enumeration; only the placed run/length drops the separators.
word_letters([Word|_], Letters, WLen) :-
    atom_chars(Word, L0), delete(L0, ' ', L1), delete(L1, '-', Letters),
    length(Letters, WLen).

% The run of cell numbers a word occupies from Start in Dir.
word_cells(_, _, 0, _, []) :- !.
word_cells(Num, Dir, K, GridLen, [Num|Rest]) :-
    K > 0, K1 is K - 1,
    next_cell(Dir, Num, GridLen, Num2),
    word_cells(Num2, Dir, K1, GridLen, Rest).

% --- geometry helpers -----------------------------------------------------

cell_rc(Cell, GridLen, R, C) :-
    R is (Cell - 1) // GridLen,
    C is (Cell - 1) mod GridLen.

placed_bbox(Placed, GridLen, bbox(MinR, MaxR, MinC, MaxC), Area) :-
    findall(R-C,
            ( member(PW, Placed), get_dict(cells, PW, Cells), member(Cell, Cells),
              cell_rc(Cell, GridLen, R, C) ),
            RCs),
    RCs = [_|_],
    pairs_keys_values(RCs, Rs, Cs),
    min_list(Rs, MinR), max_list(Rs, MaxR),
    min_list(Cs, MinC), max_list(Cs, MaxC),
    Area is (MaxR - MinR + 1) * (MaxC - MinC + 1).

% --- shared metric predicates over a placed layout (design-spec §6.4) ------
% Pure measurement: arrange.pl reads them as optimizer signals, `lint` as
% validators.

other_dir(across, down).
other_dir(down, across).

% checked cell = in an across word AND a down word.
checked_cells(Placed, Count) :-
    dir_cells(Placed, across, AcrossCells),
    dir_cells(Placed, down, DownCells),
    ord_intersection(AcrossCells, DownCells, Both),
    length(Both, Count).

dir_cells(Placed, Dir, Set) :-
    findall(Cell,
            ( member(PW, Placed), get_dict(dir, PW, Dir),
              get_dict(cells, PW, Cells), member(Cell, Cells) ),
            Cs),
    sort(Cs, Set).

% --- per-word checkedness: one canonical bitmap, everything derived from it ----
% A "checked" cell is covered by a perpendicular word. dir_cells/3 is a findall +
% sort over ALL placed words, so re-deriving it per rule per word is costly (P4).
% layout_dir_cells/2 computes both directions ONCE; lint threads the result
% through every per-word rule, and word_checked_bitmap/3 turns it into the word's
% per-cell 1/0 flags. Checked count (sum) and max-unchecked-run (a fold) derive
% from the bitmap, so "checkedness" has a single source of truth.

% Both directions' covered-cell ordsets over a layout, computed once.
layout_dir_cells(Placed, dircells(AcrossCells, DownCells)) :-
    dir_cells(Placed, across, AcrossCells),
    dir_cells(Placed, down, DownCells).

% A word's per-cell checked bitmap (1 = crossed by a perpendicular word) from a
% precomputed dircells/2: an across word's crossings live in the down set, and
% vice versa.
word_checked_bitmap(W, dircells(AcrossCells, DownCells), Bits) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    ( Dir == across -> Perp = DownCells ; Perp = AcrossCells ),
    checked_bits(Cells, Perp, Bits).

% Flag each cell 1/0 by membership in the perpendicular covered-cell ordset.
checked_bits(Cells, Perp, Bits) :- maplist(cell_checked(Perp), Cells, Bits).
cell_checked(Perp, C, B) :- ( ord_memberchk(C, Perp) -> B = 1 ; B = 0 ).

% Derived from the bitmap: number of checked cells; longest UNchecked run.
bits_checked_count(Bits, Count) :- sum_list(Bits, Count).
bits_max_unch_run(Bits, MaxRun) :- max_unch_run(Bits, 0, 0, MaxRun).
max_unch_run([], _, M, M).
max_unch_run([B|Bs], Cur, M0, M) :-
    ( B =:= 0 -> Cur1 is Cur + 1 ; Cur1 = 0 ),
    M1 is max(M0, Cur1),
    max_unch_run(Bs, Cur1, M1, M).

% The "half-checked" threshold: at least ceil(L/2) of a word's L cells must be
% crossings. Defined once here so lint's checked_half rule reuses it rather than
% re-deriving ceil(L/2) (P10).
word_half_threshold(L, T) :- T is (L + 1) // 2.

% A word is "half-checked" iff at least ceil(L/2) of its cells are crossings.
word_meets_half(W, Placed) :-
    get_dict(cells, W, Cells), length(Cells, L),
    word_checked_count(W, Placed, CC),
    word_half_threshold(L, T),
    CC >= T.

% Convenience (W, Placed) forms: derive from the same bitmap primitives, but over
% only THIS word's perpendicular direction. arrange.pl reads word_checked_count/3
% per word during rescore and needs just the one direction, so these compute
% dir_cells/3 once (not both, as lint's threaded path does).
word_checked_count(W, Placed, Count) :-
    word_perp_bits(W, Placed, Bits),
    bits_checked_count(Bits, Count).

word_max_unch_run(W, Placed, MaxRun) :-
    word_perp_bits(W, Placed, Bits),
    bits_max_unch_run(Bits, MaxRun).

word_perp_bits(W, Placed, Bits) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    checked_bits(Cells, ODCells, Bits).
