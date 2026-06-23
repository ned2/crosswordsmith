% quality.pl - greedy density-construction layout engine (v1b of
% docs/cryptic-layout-spec.md). Reuses crossword.pl's legality (init_grid,
% start_loc, find_intersecting_word, assign_word incl. I5 no_word_merge,
% fits_on_grid, next_cell), clue numbering and JSON emit. Consult AFTER
% crossword.pl.
%
% Unlike the strategies in crossword.pl (which return the first all-words layout
% or fail), this CONSTRUCTS toward density: greedily place the
% highest-(new-crossings, then lowest-bbox-growth) legal placement each step,
% over all four start locations, keeping the best-Q result and DROPPING words
% that cannot be placed (graceful degradation). One pass, no backtracking.

:- use_module(library(ordsets)).
:- use_module(library(apply)).

% quality_layout(+Words, -BestPlaced, -BestDropped, -BestGrid)
% Best greedy layout across candidate grid sizes x the four start locations.
% The engine picks the canvas (a tighter grid forces density; a looser one
% avoids drops), so the caller does not supply GridLen. Q is lexicographic:
% most words placed, then most checked cells, then smallest bounding box.
% quality_layout(+Words, +Floors, -BestPlaced, -BestDroppedAnswers, -BestGrid)
% Floors is a dict floors{min_half:M, max_unch_run:U, all_words:A} where each
% value is `off` (not imposed; Q handles it softly) or active (on / an integer).
% Floor groundedness gives the modes: all off = auto; some on = manual/partial;
% all on = strict. A bound per-word floor (min_half / max_unch_run) DROPS words
% that violate it (drop-to-satisfy); `all_words` rejects any layout with drops.
% The per-construction step (greedy_construct -> greedy_loop) is cut-free; the
% sweep, scoring and floors wrap it.
quality_layout(Words, Floors, BestPlaced, BestDropped, BestGrid) :-
    grid_candidates(Words, Grids),
    seed_candidates(Words, Seeds),
    findall(q(NP, Score)-pdg(FPlaced, AllDropped, GridLen),
            ( member(GridLen, Grids),
              start_locs(Locs), member(Loc, Locs),
              member(Seed, Seeds),
              greedy_construct(Words, GridLen, Loc, Seed, Placed0, Dropped0),
              floor_drop(Floors, Placed0, FPlaced, FloorDropped),
              findall(A, member([A|_], Dropped0), CantPlace),
              append(CantPlace, FloorDropped, AllDropped),
              floors_ok(Floors, AllDropped),
              length(FPlaced, NP),
              layout_score(FPlaced, GridLen, Score) ),
            Results),
    Results = [_|_],
    sort(1, @>=, Results, [_-pdg(BestPlaced, BestDropped, BestGrid)|_]).

% Q (secondary to words-placed): reward checked cells, penalise a large/elongated
% bounding box (a web table-of-contents must render compactly, so checking and
% compactness compete rather than checking strictly dominating). Weights are the
% editorial knob (quality_weights/3).
layout_score(Placed, GridLen, Score) :-
    checked_cells(Placed, Checked),
    placed_bbox(Placed, GridLen, bbox(MinR, MaxR, MinC, MaxC), Area),
    H is MaxR - MinR + 1, W is MaxC - MinC + 1,
    Elong is abs(H - W),
    quality_weights(WCheck, WArea, WElong),
    Score is WCheck * Checked - WArea * Area - WElong * Elong.

% w_check : w_area : w_elong. Checking leads, but bbox area and elongation pull
% layouts toward a compact square.
quality_weights(6, 1, 2).

% all_words floor: reject any candidate that dropped a word.
floors_ok(Floors, AllDropped) :-
    ( get_dict(all_words, Floors, on) -> AllDropped == [] ; true ).

% Seed candidates: the K longest words (restart diversity to escape greedy
% local optima). K shrinks as the set grows, so the total grid x start x seed
% sweep stays bounded (and deterministic).
seed_candidates(Words, Seeds) :-
    length(Words, N),
    K is max(1, min(5, 80 // N)),
    map_list_to_pairs(neg_answer_len, Words, Pairs),
    keysort(Pairs, Sorted),
    pairs_values(Sorted, ByLenDesc),
    ( length(Prefix, K), append(Prefix, _, ByLenDesc)
    ->  Seeds = Prefix
    ;   Seeds = ByLenDesc ).

% Candidate square sizes from the word set: hold the total letters at a few
% target densities, but never smaller than the longest word. A small, fixed,
% deterministic set (no wall-clock budget needed). On an empty word set max_list/2
% fails (it does not error), so this fails and quality_solve reports "no layout" -
% the intended fail-soft for a degenerate/empty puzzle.
grid_candidates(Words, Sizes) :-
    maplist(word_letter_count, Words, Lens),
    sum_list(Lens, Total),
    max_list(Lens, MaxL),
    % target densities from dense-interlock down to sparse: a dense word set
    % wins on a tight grid (smaller bbox), a sparse one needs the looser grids
    % to place all its words (Q prefers most-placed). Three candidates bound the
    % work (with the per-set seed cap) - the compactness-aware Q would not pick a
    % much larger canvas anyway, so we do not generate slow oversized grids.
    findall(S,
            ( member(D, [0.55, 0.38, 0.26]),
              S0 is ceiling(sqrt(Total / D)),
              S is max(S0, MaxL + 1) ),
            Raw),
    sort(Raw, Sizes).

word_letter_count(Entry, N) :- word_letters(Entry, _, N).

% --- greedy construction from one start location --------------------------

greedy_construct(Words, GridLen, Loc, Seed, Placed, Dropped) :-
    init_grid(GridLen, G0),
    start_loc(Loc, GridLen, StartNum, StartDir),
    remove_x(Seed, Words, Rest),
    seed_word(Seed, StartNum, StartDir, GridLen, G0, SeedPW, G1),
    greedy_loop(Rest, [SeedPW], GridLen, G1, Placed, Dropped).

neg_answer_len([A|_], NL) :- atom_length(A, L), NL is -L.

% Seed: place the chosen word at the fixed start cell/direction (no crossings).
seed_word(Entry, Start, Dir, GridLen, GIn, PW, GOut) :-
    word_letters(Entry, Letters, WLen),
    fits_on_grid(Dir, Start, WLen, GridLen),
    Entry = [Word|_],
    assign_word(Word, Letters, WLen, Start, Dir, GridLen, [], GIn, PW, GOut).

word_letters([Word|_], Letters, WLen) :-
    atom_chars(Word, L0), delete(L0, ' ', Letters), length(Letters, WLen).

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

% --- quality metrics over a placed layout ---------------------------------

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

% --- emit (reuse crossword.pl) + report -----------------------------------

% --- floors: drop-to-satisfy per-word floors ------------------------------

other_dir(across, down).
other_dir(down, across).

% Iteratively drop placed words that violate an active per-word floor (their
% removal can un-check a crosser, so re-evaluate after each drop). Returns the
% surviving placed words and the answers dropped to satisfy the floors.
floor_drop(Floors, Placed, FinalPlaced, FloorDropped) :-
    floor_drop_(Floors, Placed, [], FinalPlaced, FloorDropped).
floor_drop_(Floors, Placed, Acc, FinalPlaced, FloorDropped) :-
    ( select(W, Placed, Rest), violates_floor(W, Placed, Floors)
    ->  get_dict(answer, W, A),
        floor_drop_(Floors, Rest, [A|Acc], FinalPlaced, FloorDropped)
    ;   FinalPlaced = Placed, FloorDropped = Acc ).

violates_floor(W, Placed, Floors) :-
    get_dict(min_half, Floors, on),
    \+ word_meets_half(W, Placed).
violates_floor(W, Placed, Floors) :-
    get_dict(max_unch_run, Floors, K), integer(K),
    word_max_unch_run(W, Placed, R), R > K.

word_meets_half(W, Placed) :-
    get_dict(cells, W, Cells), length(Cells, L),
    word_checked_count(W, Placed, CC),
    CC >= (L + 1) // 2.

% W's cells that are crossings = also covered by a perpendicular word.
word_checked_count(W, Placed, Count) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    findall(x, ( member(C, Cells), ord_memberchk(C, ODCells) ), Xs),
    length(Xs, Count).

word_max_unch_run(W, Placed, MaxRun) :-
    get_dict(cells, W, Cells), get_dict(dir, W, Dir),
    other_dir(Dir, OD), dir_cells(Placed, OD, ODCells),
    maxrun(Cells, ODCells, 0, 0, MaxRun).
maxrun([], _, _, M, M).
maxrun([C|Cs], ODCells, Cur, M0, MaxRun) :-
    ( ord_memberchk(C, ODCells) -> Cur1 = 0 ; Cur1 is Cur + 1 ),
    M1 is max(M0, Cur1),
    maxrun(Cs, ODCells, Cur1, M1, MaxRun).

active_floors(Floors, Active) :-
    findall(K=V, ( member(K, [min_half, max_unch_run, all_words]),
                   get_dict(K, Floors, V), V \== off ),
            Active).

% --- entry point ----------------------------------------------------------

% quality_solve(+Words, +Floors): construct the best layout (engine picks the
% grid), emit it as JSON on stdout and a one-line report on stderr. Fails (no
% output) when active floors cannot be satisfied.
quality_solve(Words, Floors) :-
    ( quality_layout(Words, Floors, Placed, Dropped, GridLen)
    ->  assign_clue_numbers(Placed, Numbered),
        emit_json(Numbered, Words, GridLen),
        length(Placed, NP), length(Dropped, ND),
        format(user_error, "quality: grid ~w, placed ~w, dropped ~w ~w~n",
               [GridLen, NP, ND, Dropped])
    ;   active_floors(Floors, Active),
        format(user_error,
               "quality: no layout satisfies floors ~w (try relaxing)~n",
               [Active]),
        fail ).

% Auto mode: no floors imposed.
no_floors(floors{min_half:off, max_unch_run:off, all_words:off}).
