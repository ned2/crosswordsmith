% metrics.pl - the shared layout metric predicates (design-spec §6.4).
% (Formerly quality.pl; renamed in Phase 2 of the source-structure migration.
% The greedy density constructor that lived here moved to arrange.pl — its
% only consumer — in Phase 3.)
%
% Pure measurement over a placed layout (dir_cells, checked_cells,
% word_checked_count, word_meets_half, word_max_unch_run, placed_bbox,
% word_cells, ...); next_cell/4 is the only core predicate used (imported
% below, by word_cells/5). Consumed by arrange.pl as optimizer signals and by
% `lint` as validators. Per spec §4 these stay a separate module — lint
% depends only on the JSON contract plus this metric layer, never on the
% solver substrate.
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

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5).
:- use_module(library(ordsets), [ord_intersection/3, ord_memberchk/2]).
:- use_module(library(apply), [exclude/3, foldl/4, maplist/3]).
:- use_module(library(lists), [member/2, sum_list/2]).

% next_cell/4 (word_cells/5's stepper) plus the placed-word record (pw/8)
% field accessors — the metric predicates read cells/dir off each placed word.
:- use_module(crosswordsmith(core), [next_cell/4, pw_cells/2, pw_dir/2]).

% --- shared answer/footprint helpers ----------------------------------------

%!  word_letters(+Entry:list, -Letters:list, -WLen:integer) is det.
%
%   The "placement footprint" of an input word entry [Answer|_]: its letters
%   with word separators (spaces AND hyphens) removed - those are enumeration
%   markers, not grid cells - and its placed length. The original answer atom
%   is carried to emit (placed_to_word/4) so export still derives the
%   enumeration; only the placed run/length drops the separators.
word_letters([Word|_], Letters, WLen) :-
    atom_chars(Word, L0),
    exclude(separator_char, L0, Letters),
    length(Letters, WLen).

% Word-separator characters: enumeration markers, not grid cells. (delete/3,
% the previous strip, is deprecated in SWI 10 and took two passes.)
separator_char(' ').
separator_char('-').

%!  word_cells(+Start:integer, +Dir:oneof([across,down]), +Len:nonneg,
%!             +GridLen:integer, -Cells:list(integer)) is det.
%
%   The run of Len cell numbers a word occupies from Start in Dir.
word_cells(_, _, 0, _, []) :- !.
word_cells(Num, Dir, K, GridLen, [Num|Rest]) :-
    K > 0, K1 is K - 1,
    next_cell(Dir, Num, GridLen, Num2),
    word_cells(Num2, Dir, K1, GridLen, Rest).

% --- geometry helpers -----------------------------------------------------

%!  cell_rc(+Cell:integer, +GridLen:integer, -R:integer, -C:integer) is det.
%
%   Cell number (1-based, row-major) to 0-based row/column.
cell_rc(Cell, GridLen, R, C) :-
    R is (Cell - 1) // GridLen,
    C is (Cell - 1) mod GridLen.

%!  placed_bbox(+Placed:list, +GridLen:integer, -BBox, -Area:integer) is semidet.
%
%   BBox = bbox(MinR, MaxR, MinC, MaxC), the tight bounding box (0-based
%   rows/cols) of all cells covered by the placed pw/8 records, and its
%   Area. Fails for an empty placement set (no cells, no box).
%
%   One rc_extend/3 fold computes all four bounds in a single pass over the
%   R-C pairs (was pairs_keys_values + min_list/max_list x4 - five passes;
%   audit C34). The head match on RCs keeps the empty-set FAILURE contract.
placed_bbox(Placed, GridLen, bbox(MinR, MaxR, MinC, MaxC), Area) :-
    findall(R-C,
            ( member(PW, Placed), pw_cells(PW, Cells), member(Cell, Cells),
              cell_rc(Cell, GridLen, R, C) ),
            RCs),
    RCs = [R0-C0|Rest],                    % [] -> fail (no cells, no box)
    foldl(rc_extend, Rest, b(R0, R0, C0, C0), b(MinR, MaxR, MinC, MaxC)),
    Area is (MaxR - MinR + 1) * (MaxC - MinC + 1).

rc_extend(R-C, b(Rlo, Rhi, Clo, Chi), b(Rlo2, Rhi2, Clo2, Chi2)) :-
    Rlo2 is min(Rlo, R), Rhi2 is max(Rhi, R),
    Clo2 is min(Clo, C), Chi2 is max(Chi, C).

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
            ( member(PW, Placed), pw_dir(PW, Dir),
              pw_cells(PW, Cells), member(Cell, Cells) ),
            Cs),
    sort(Cs, Set).

% --- per-word checkedness: one canonical bitmap, everything derived from it ----
% A "checked" cell is covered by a perpendicular word. dir_cells/3 is a findall +
% sort over ALL placed words, so re-deriving it per rule per word is costly (P4).
% layout_dir_cells/2 computes both directions ONCE; lint threads the result
% through every per-word rule, and word_checked_bitmap/3 turns it into the word's
% per-cell 1/0 flags. Checked count (sum) and max-unchecked-run (a fold) derive
% from the bitmap, so "checkedness" has a single source of truth.

%!  layout_dir_cells(+Placed:list, -DirCells) is det.
%
%   Both directions' covered-cell ordsets over a layout, computed once:
%   DirCells = dircells(AcrossCells, DownCells). lint threads the result
%   through every per-word rule (P4).
layout_dir_cells(Placed, dircells(AcrossCells, DownCells)) :-
    dir_cells(Placed, across, AcrossCells),
    dir_cells(Placed, down, DownCells).

%!  word_checked_bitmap(+W, +DirCells, -Bits:list(integer)) is det.
%
%   W's per-cell checked bitmap (1 = crossed by a perpendicular word) from a
%   precomputed dircells/2 (layout_dir_cells/2): an across word's crossings
%   live in the down set, and vice versa.
word_checked_bitmap(W, dircells(AcrossCells, DownCells), Bits) :-
    pw_cells(W, Cells), pw_dir(W, Dir),
    ( Dir == across -> Perp = DownCells ; Perp = AcrossCells ),
    checked_bits(Cells, Perp, Bits).

% Flag each cell 1/0 by membership in the perpendicular covered-cell ordset.
checked_bits(Cells, Perp, Bits) :- maplist(cell_checked(Perp), Cells, Bits).
cell_checked(Perp, C, B) :- ( ord_memberchk(C, Perp) -> B = 1 ; B = 0 ).

%!  bits_checked_count(+Bits:list(integer), -Count:integer) is det.
%!  bits_max_unch_run(+Bits:list(integer), -MaxRun:integer) is det.
%
%   Derived from the checked bitmap: the number of checked cells; the
%   longest UNchecked run.
bits_checked_count(Bits, Count) :- sum_list(Bits, Count).
bits_max_unch_run(Bits, MaxRun) :- max_unch_run(Bits, 0, 0, MaxRun).
max_unch_run([], _, M, M).
max_unch_run([B|Bs], Cur, M0, M) :-
    ( B =:= 0 -> Cur1 is Cur + 1 ; Cur1 = 0 ),
    M1 is max(M0, Cur1),
    max_unch_run(Bs, Cur1, M1, M).

%!  word_half_threshold(+L:integer, -T:integer) is det.
%
%   The "half-checked" threshold: at least ceil(L/2) of a word's L cells must
%   be crossings. Defined once here so lint's checked_half rule reuses it
%   rather than re-deriving ceil(L/2) (P10).
word_half_threshold(L, T) :- T is (L + 1) // 2.

% A word is "half-checked" iff at least ceil(L/2) of its cells are crossings.
word_meets_half(W, Placed) :-
    pw_cells(W, Cells), length(Cells, L),
    word_checked_count(W, Placed, CC),
    word_half_threshold(L, T),
    CC >= T.

%!  word_checked_count(+W, +Placed:list, -Count:integer) is det.
%
%   The number of W's cells crossed by a perpendicular word of Placed - the
%   convenience (W, Placed) form of bits_checked_count/2. Derives from the
%   same bitmap primitives, but over only THIS word's perpendicular
%   direction: arrange.pl reads it per word during rescore and needs just the
%   one direction, so it computes dir_cells/3 once (not both, as lint's
%   threaded path does).
word_checked_count(W, Placed, Count) :-
    word_perp_bits(W, Placed, Bits),
    bits_checked_count(Bits, Count).

word_max_unch_run(W, Placed, MaxRun) :-
    word_perp_bits(W, Placed, Bits),
    bits_max_unch_run(Bits, MaxRun).

word_perp_bits(W, Placed, Bits) :-
    pw_cells(W, Cells), pw_dir(W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    checked_bits(Cells, ODCells, Bits).
