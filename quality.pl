% quality.pl - shared layout metrics + the Flavour-A greedy density constructor.
%
% Two things live here, both consulted AFTER crossword.pl (whose legality core
% they reuse: init_grid, start_loc, find_intersecting_word, assign_word incl. the
% I5 no_word_merge, fits_on_grid, next_cell):
%
%   1. Shared METRIC predicates - pure measurement over a placed layout
%      (dir_cells, checked_cells, word_checked_count, word_meets_half,
%      word_max_unch_run, crossing_count, placed_bbox, word_cells, ...). Consumed
%      by arrange.pl as optimizer signals and by `lint` as validators
%      (design-spec §6.4). (Lifting them into crossword.pl proper is a pending
%      tidy-up; functionally they are already the shared metric layer.)
%
%   2. The greedy density CONSTRUCTOR (greedy_construct/greedy_loop/next_move/
%      word_best_placement) - arrange.pl's best-effort and candidates path:
%      greedily place the highest-(new-crossings, then lowest-bbox-growth) legal
%      placement each step, dropping words that cannot be placed. One pass, no
%      backtracking, cut-free.
%
% The former `--quality` CLI engine (quality_solve / quality_layout /
% grid_candidates + the floor subsystem) was retired at the Phase-7 CLI cutover
% (design-spec §4); arrange.pl supersedes it. Only the reused primitives remain.

:- use_module(library(ordsets)).
:- use_module(library(apply)).

% --- seed selection (greedy restart diversity) ----------------------------

% Seed candidates: the K longest words (restart diversity to escape greedy
% local optima). K shrinks as the set grows, so the total start x seed sweep
% stays bounded (and deterministic).
seed_candidates(Words, Seeds) :-
    length(Words, N),
    K is max(1, min(5, 80 // N)),
    map_list_to_pairs(neg_answer_len, Words, Pairs),
    keysort(Pairs, Sorted),
    pairs_values(Sorted, ByLenDesc),
    ( length(Prefix, K), append(Prefix, _, ByLenDesc)
    ->  Seeds = Prefix
    ;   Seeds = ByLenDesc ).

% --- greedy construction from one start location --------------------------

greedy_construct(Words, GridLen, Loc, Seed, Placed, Dropped) :-
    init_grid(GridLen, G0),
    start_loc(Loc, GridLen, StartNum, StartDir),
    remove_x(Seed, Words, Rest),
    seed_word(Seed, StartNum, StartDir, GridLen, G0, SeedPW, G1),
    greedy_loop(Rest, [SeedPW], GridLen, G1, Placed, Dropped).

neg_answer_len(Entry, NL) :- word_letters(Entry, _, WLen), NL is -WLen.

% Seed: place the chosen word at the fixed start cell/direction (no crossings).
seed_word(Entry, Start, Dir, GridLen, GIn, PW, GOut) :-
    word_letters(Entry, Letters, WLen),
    fits_on_grid(Dir, Start, WLen, GridLen),
    Entry = [Word|_],
    assign_word(Word, Letters, WLen, Start, Dir, GridLen, [], GIn, PW, GOut).

% The "placement footprint" of an answer: its letters with word separators
% (spaces AND hyphens) removed - those are enumeration markers, not grid cells.
% The original answer atom is carried to emit (placed_to_word/4) so export still
% derives the enumeration; only the placed run/length drops the separators.
word_letters([Word|_], Letters, WLen) :-
    atom_chars(Word, L0), delete(L0, ' ', L1), delete(L1, '-', Letters),
    length(Letters, WLen).

% Place the globally best-scoring placeable word, repeat; drop the rest. The
% construction is cut-free: instead of an `( Best -> place ; stop )` if-then-else,
% the next move (or `none`) is reified as a term and dispatched on the mutually
% exclusive clause heads of best_move/2 ([] vs [_|_]). No !, ->, or \+ here or in
% word_best_placement below; only findall/sort/arithmetic remain. (This replaced
% an equivalent cut-based version; it is identical-output and faster - see
% docs/cryptic-layout-spec.md v1b.1.)
greedy_loop(Remaining, Placed, GridLen, Grid, FinalPlaced, Dropped) :-
    next_move(Remaining, Placed, GridLen, Grid, Move),
    apply_move(Move, Remaining, Placed, GridLen, FinalPlaced, Dropped).

apply_move(none, Remaining, Placed, _GridLen, Placed, Remaining).
apply_move(move(Entry, NewPW, NewGrid), Remaining, Placed, GridLen, FinalPlaced, Dropped) :-
    % next_move/5 collects candidates with findall, which COPIES each Entry, so
    % the move's Entry is NOT ==-identical to its term in Remaining when the
    % entry is non-ground (a .pl fixture's [Answer, _{...}] has an unbound dict
    % tag). remove_x (==-based) would then drop nothing and greedy_loop would
    % re-offer the just-placed word forever. Key the removal on the ground
    % answer atom instead (answers are unique, check_unique_answers/1); selectchk
    % keeps the original, uncopied tail. (Same term-copy footgun the MRV path
    % avoids via map_list_to_pairs - see crossword.pl select_inc.)
    Entry = [Answer|_],
    selectchk([Answer|_], Remaining, Remaining1),
    greedy_loop(Remaining1, [NewPW|Placed], GridLen, NewGrid, FinalPlaced, Dropped).

% The best placeable word as move(Entry,PW,Grid), or `none` when nothing fits.
% Among remaining words with a legal placement, the one whose best placement
% scores highest (most new crossings, then least bbox growth). The [] vs [_|_]
% dispatch in best_move/2 is what lets greedy_loop avoid an if-then-else.
next_move(Remaining, Placed, GridLen, Grid, Move) :-
    placed_bbox(Placed, GridLen, BBox, _),
    findall(Score-cand(Entry, PW, G1),
            ( member(Entry, Remaining),
              word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, PW, G1) ),
            Cands),
    best_move(Cands, Move).

best_move([], none).
best_move([C|Cs], move(Entry, PW, G1)) :-
    sort(1, @>=, [C|Cs], [_-cand(Entry, PW, G1)|_]).

% Best legal placement of one word on the current grid, cut-free: assign_word
% legality is folded INTO the findall generator so only legal placements are
% collected (keyed by the density score), and the head of the @>=-sorted list is
% the best - no first-solution `!`. Folding the legality filter ahead of scoring
% means placement_key runs only for legal candidates (not every crossing
% candidate), which is why this is faster than the cut version (see spec v1b.1).
word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, PW, GOut) :-
    word_letters(Entry, Letters, WLen),
    Entry = [Word|_],
    findall(Key-place(PW1, G1),
            ( find_intersecting_word(Letters, WLen, Placed, GridLen, Start, Dir),
              assign_word(Word, Letters, WLen, Start, Dir, GridLen, Placed, Grid, PW1, G1),
              placement_key(Letters, Start, Dir, WLen, GridLen, Grid, BBox, Key) ),
            Keyed),
    Keyed = [_|_],
    sort(1, @>=, Keyed, [Score-place(PW, GOut)|_]).

% Density score for a candidate placement: crossings dominate, bbox-growth
% breaks ties (smaller is better). 10000 > any plausible bbox area.
placement_key(Letters, Start, Dir, WLen, GridLen, Grid, BBox, Key) :-
    crossing_count(Letters, Start, Dir, GridLen, Grid, Crossings),
    word_cells(Start, Dir, WLen, GridLen, Cells),
    bbox_growth(BBox, Cells, GridLen, Growth),
    Key is Crossings * 10000 - Growth.

crossing_count(Letters, Start, Dir, GridLen, Grid, Count) :-
    cc_(Letters, Start, Dir, GridLen, Grid, 0, Count).
cc_([], _, _, _, _, A, A).
cc_([L|Ls], Num, Dir, GridLen, Grid, A0, Count) :-
    ( get_assoc(Num, Grid, L) -> A1 is A0 + 1 ; A1 = A0 ),
    next_cell(Dir, Num, GridLen, Num2),
    cc_(Ls, Num2, Dir, GridLen, Grid, A1, Count).

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

bbox_growth(bbox(MinR, MaxR, MinC, MaxC), NewCells, GridLen, Growth) :-
    OldArea is (MaxR - MinR + 1) * (MaxC - MinC + 1),
    foldl(extend_cell(GridLen), NewCells, b(MinR, MaxR, MinC, MaxC), b(R0, R1, C0, C1)),
    NewArea is (R1 - R0 + 1) * (C1 - C0 + 1),
    Growth is NewArea - OldArea.

extend_cell(GridLen, Cell, b(MinR, MaxR, MinC, MaxC), b(MinR2, MaxR2, MinC2, MaxC2)) :-
    cell_rc(Cell, GridLen, R, C),
    MinR2 is min(MinR, R), MaxR2 is max(MaxR, R),
    MinC2 is min(MinC, C), MaxC2 is max(MaxC, C).

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

% A word is "half-checked" iff at least ceil(L/2) of its cells are crossings.
word_meets_half(W, Placed) :-
    get_dict(cells, W, Cells), length(Cells, L),
    word_checked_count(W, Placed, CC),
    CC >= (L + 1) // 2.

% W's cells that are crossings = also covered by a perpendicular word.
word_checked_count(W, Placed, Count) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    aggregate_all(count, ( member(C, Cells), ord_memberchk(C, ODCells) ), Count).

% Longest run of consecutive UNchecked cells along the word.
word_max_unch_run(W, Placed, MaxRun) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    maxrun(Cells, ODCells, 0, 0, MaxRun).
maxrun([], _, _, M, M).
maxrun([C|Cs], ODCells, Cur, M0, MaxRun) :-
    ( ord_memberchk(C, ODCells) -> Cur1 = 0 ; Cur1 is Cur + 1 ),
    M1 is max(M0, Cur1),
    maxrun(Cs, ODCells, Cur1, M1, MaxRun).
