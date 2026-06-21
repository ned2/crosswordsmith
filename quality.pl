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
quality_layout(Words, BestPlaced, BestDropped, BestGrid) :-
    grid_candidates(Words, Grids),
    seed_candidates(Words, Seeds),
    findall(q(NP, Checked, NegArea)-pdg(Placed, Dropped, GridLen),
            ( member(GridLen, Grids),
              start_locs(Locs), member(Loc, Locs),
              member(Seed, Seeds),
              greedy_construct(Words, GridLen, Loc, Seed, Placed, Dropped),
              length(Placed, NP),
              checked_cells(Placed, Checked),
              placed_bbox(Placed, GridLen, _, Area),
              NegArea is -Area ),
            Results),
    Results = [_|_],
    sort(1, @>=, Results, [_-pdg(BestPlaced, BestDropped, BestGrid)|_]).

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
% deterministic set (no wall-clock budget needed).
grid_candidates(Words, Sizes) :-
    maplist(word_letter_count, Words, Lens),
    sum_list(Lens, Total),
    max_list(Lens, MaxL),
    % target densities from dense-interlock down to sparse: a dense word set
    % wins on a tight grid (smaller bbox), a sparse one needs the looser grids
    % to place all its words (Q prefers most-placed). Spanning both avoids
    % under-sizing sparse sets.
    findall(S,
            ( member(D, [0.55, 0.40, 0.28, 0.18]),
              S0 is ceiling(sqrt(Total / D)),
              S is max(S0, MaxL + 1) ),
            Raw),
    sort(Raw, Sizes).

word_letter_count([A|_], N) :-
    atom_chars(A, Cs), delete(Cs, ' ', Ls), length(Ls, N).

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

% Place the globally best-scoring placeable word, repeat; drop the rest.
greedy_loop(Remaining, Placed, GridLen, Grid, FinalPlaced, Dropped) :-
    ( best_global_placement(Remaining, Placed, GridLen, Grid, Entry, NewPW, NewGrid)
    ->  remove_x(Entry, Remaining, Remaining1),
        greedy_loop(Remaining1, [NewPW|Placed], GridLen, NewGrid, FinalPlaced, Dropped)
    ;   FinalPlaced = Placed, Dropped = Remaining ).

% Among remaining words with a legal placement, the one whose best placement
% scores highest (most new crossings, then least bbox growth).
best_global_placement(Remaining, Placed, GridLen, Grid, BestEntry, BestPW, BestGrid) :-
    placed_bbox(Placed, GridLen, BBox, _),
    findall(Score-cand(Entry, PW, G1),
            ( member(Entry, Remaining),
              word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, PW, G1) ),
            Cands),
    Cands = [_|_],
    sort(1, @>=, Cands, [_-cand(BestEntry, BestPW, BestGrid)|_]).

% Best legal placement of one word on the current grid.
word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, PW, GOut) :-
    word_letters(Entry, Letters, WLen),
    Entry = [Word|_],
    findall(Key-(Start-Dir),
            ( find_intersecting_word(Letters, WLen, Placed, GridLen, Start, Dir),
              placement_key(Letters, Start, Dir, WLen, GridLen, Grid, BBox, Key) ),
            Keyed),
    Keyed = [_|_],
    sort(1, @>=, Keyed, Sorted),
    member(Score-(Start-Dir), Sorted),
    assign_word(Word, Letters, WLen, Start, Dir, GridLen, Placed, Grid, PW, GOut),
    !.

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

% quality_solve(+Words): construct the best layout (engine picks the grid),
% emit it as JSON on stdout (placed subset) and a one-line report on stderr.
quality_solve(Words) :-
    quality_layout(Words, Placed, Dropped, GridLen),
    assign_clue_numbers(Placed, Numbered),
    emit_json(Numbered, Words, GridLen),
    length(Placed, NP), length(Dropped, ND),
    findall(A, member([A|_], Dropped), DroppedAnswers),
    format(user_error, "quality: grid ~w, placed ~w, dropped ~w ~w~n",
           [GridLen, NP, ND, DroppedAnswers]).
