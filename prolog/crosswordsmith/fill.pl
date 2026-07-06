% fill.pl - Flavour-B grid-first, open-dictionary auto-fill (design-spec §8.4;
% OD-1..4 resolved by DP-1/DP-2). Take a pre-validated legal blocked grid (a §8.3
% stock-grid mask), derive its slots, and fill every slot with a dictionary word
% subject to crossing constraints, with the user's words pinned as seeds (the
% §6.6 fragment primitive). Output is the canonical layout, so it composes with
% lint/export.
%
% The grid model: each white cell is one logical VARIABLE, shared between its
% across and down slot. Assigning a word to a slot unifies the slot's variable
% list with the word's letters; crossings are then consistent by construction
% (the shared variable), and dead-ends backtrack. MRV ordering (fewest matching
% candidates first) + a node/inference budget make it deterministic and bounded.
%
% Exports: fill_solve/4 (the CLI seam) plus the prebuilt-index seams
% (fill_solve_index/5, fill_save_index/2,3); every other predicate is
% internal (tests reach them as crosswordsmith_fill:Pred(...)). All
% project dependencies are explicit imports below (stockgrid, metrics,
% arrange, core).

:- module(crosswordsmith_fill,
          [ fill_solve/4,
            fill_solve_index/5,
            fill_save_index/2,
            fill_save_index/3
          ]).

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5). (No JSON
% import: this module's emit delegates to core's emit_json/3 / arrange's
% emit_arrange/4 — the old library(http/json) line here was vestigial.)
:- use_module(library(apply), [convlist/3, exclude/3, foldl/4, maplist/3]).
:- use_module(library(lists),
              [append/3, member/2, nextto/3, nth0/3, select/3, subtract/3]).
:- use_module(library(ordsets),
              [list_to_ord_set/2, ord_intersect/2, ord_intersection/3]).
:- use_module(library(assoc),
              [ assoc_to_list/2, assoc_to_values/2, gen_assoc/3,
                get_assoc/3, list_to_assoc/2, ord_list_to_assoc/2 ]).
:- use_module(library(pairs),
              [group_pairs_by_key/2, map_list_to_pairs/3, pairs_keys_values/3]).
% F-L2: dict-file SHA-256 (artifact integrity). NB library(sha)'s file defines
% module `crypto_hash` on 10.1.10; both predicates are on its export list.
:- use_module(library(sha), [sha_hash/3, hash_atom/2]).
% F-L2: fast_write/fast_read artifact I/O.
:- use_module(library(fastrw), [fast_read/2, fast_write/2]).
% option/2,3: the artifact Meta is a keyed list of unary Key(Value) terms -
% literally an SWI option list - and fill_save_index's Options is one too.
:- use_module(library(option), [option/2, option/3]).
% read_file_to_string/3: whole-file reads (dict lines + the SHA-256
% fingerprint) without the hand-rolled open/read_string/close dance (C39).
:- use_module(library(readutil), [read_file_to_string/3]).

% Slot derivation from a stock-grid mask.
:- use_module(crosswordsmith(stockgrid),
              [stockgrid_load/2, mask_white_cells/3, grid_run/4]).

% word_letters/3: the separator-stripped placement footprint of a seed answer.
:- use_module(crosswordsmith(metrics), [word_letters/3]).

% load_fragment/3: seeds arrive in the §6.6 fragment format (canonical form
% only - /3 rejects the thin arrange-only list form with a shaped error,
% since fill matches seeds by raw cell numbers on its own grid). emit_arrange/4:
% emit_fill's `max` mode delegates the cropped emit to arrange.
:- use_module(crosswordsmith(arrange), [load_fragment/3, emit_arrange/4]).

% Numbering + the canonical JSON emit for the filled layout, plus the placed-word
% record (pw/8) accessor fill uses to recover each answer for the emit metadata.
% check_unique_answers/1: the emit-boundary unique-answers guard the solve path
% runs in crossword/4; fill bypasses crossword/4, so emit_fill/4 runs it itself.
:- use_module(crosswordsmith(core),
              [ assign_clue_numbers/2, check_unique_answers/1, emit_json/3,
                verbose_report/2, pw_answer/2 ]).

fill_budget(800_000_000).   % inference budget (determinism via INV-2, bounded)


% --- build the slots from a stock-grid mask ----------------------------------
% Each slot is slot(Start, Dir, Cells, Vars): Vars are the shared per-cell
% logical variables (so an across and a down slot crossing at a cell share that
% cell's variable). Returns the slots + the cell->var assoc.
fill_grid(GridFile, Size, Slots, CellVar) :-
    stockgrid_load(GridFile, grid(_Name, Size, Mask)),
    mask_white_cells(Mask, Size, WhiteSet),
    init_cell_vars(WhiteSet, CellVar),
    % findall COPIES its template, which would rename the shared cell variables
    % apart and break crossings - so collect ground slot specs first, then wire
    % in the shared variables (from the one CellVar assoc) outside findall.
    findall(spec(Start, Dir, Cells),
            ( member(Dir, [across, down]),
              grid_run(Size, WhiteSet, Dir, Cells),
              Cells = [Start|_] ),
            Specs),
    maplist(spec_slot(CellVar), Specs, Slots).

spec_slot(CellVar, spec(Start, Dir, Cells), slot(Start, Dir, Cells, Vars)) :-
    slot_vars(Cells, CellVar, Vars).

% WhiteCells is a strictly-ascending ordset (mask_white_cells/3 ends in
% list_to_ord_set/2), so ord_list_to_assoc/2 constructs the balanced assoc
% directly - no per-element put_assoc rebalance (C31).
init_cell_vars(WhiteCells, Assoc) :-
    pairs_keys_values(Pairs, WhiteCells, _FreshVars),
    ord_list_to_assoc(Pairs, Assoc).

% NB: a yall lambda `[C,V]>>get_assoc(C, CellVar, V)` would copy_term its free
% variable CellVar on every call - cloning the assoc and the unbound cell
% variables inside it, so nothing would be shared. A named helper passes the
% assoc as a plain argument (no copy), keeping the crossing cells' variables
% genuinely shared between the across and down slots.
slot_vars(Cells, CellVar, Vars) :-
    maplist(cell_var(CellVar), Cells, Vars).
cell_var(CellVar, Cell, Var) :- get_assoc(Cell, CellVar, Var).


% --- seeds: pin fragment words at their slots (hard pins, OD-3) ---------------
% A seed answer is unified into the slot whose cells/direction match; a seed that
% matches no slot, or clashes with another seed at a crossing, is an error.
% Pin every seed. SeededKeys collects each pinned slot's Start-Dir so the search
% and the unfillable report can EXEMPT them: a seed is a hard pin (a setter's
% theme/own word) and need NOT be a dictionary word (R2/OD-3). Its letters still
% constrain crossing slots via the shared cell variables, and a slot completed
% by crossings (not a seed) is still validated against the dictionary.
apply_seeds(SeedFile, Slots, SeededKeys) :-
    load_fragment(SeedFile, _GridLen, Frags),
    check_unique_seed_answers(Frags),
    foldl(apply_seed(Slots), Frags, [], SeededKeys).

% Two seeds pinning the SAME answer would collide at emit (crosswords do not
% repeat answers; the emit metadata join is keyed on the answer), and neither
% slot is searched, so the search-side dedup (seed_used/3) cannot catch that
% shape. Reject up front with a clean hooked error, consistent with
% fill_seed_no_slot / fill_seed_clash. Detection mirrors core's
% check_unique_answers/1: msort keeps duplicates, equal neighbours = a dup.
check_unique_seed_answers(Frags) :-
    findall(A, member(frag(A, _, _, _), Frags), Answers),
    msort(Answers, Sorted),
    (   nextto(Dup, Dup, Sorted)
    ->  throw(error(fill_seed_duplicate(Dup), _))
    ;   true
    ).

apply_seed(Slots, frag(Answer, Dir, _Start, CellNums), SeededIn, SeededOut) :-
    ( member(slot(SStart, Dir, CellNums, Vars), Slots)
    ->  word_letters([Answer], Letters, _),
        ( Vars = Letters
        ->  SeededOut = [SStart-Dir|SeededIn]
        ;   throw(error(fill_seed_clash(Answer), _)) )
    ;   throw(error(fill_seed_no_slot(Answer), _)) ).

% A slot is a seed pin (exempt from search + empty_slots) iff its Start-Dir was
% recorded by apply_seed (Start+Dir uniquely identify a slot).
seeded_slot(SeededKeys, slot(Start, Dir, _, _)) :- memberchk(Start-Dir, SeededKeys).


% --- dictionary: in-memory pattern index (OD-2) ------------------------------
% DictByLen: assoc Len -> list of words (each a list of upper-case char atoms).
% Index: assoc k(Len,Pos,Char) -> ordset of indices into DictByLen[Len]. A slot's
% candidates = the words of its length matching its already-bound positions
% (intersection of the index sets), or all words of that length if nothing bound.
load_dict(File, DictByLen, Index) :-
    read_file_lines(File, Lines),
    % convlist, not findall+member+guard (C32): the same word list in the same
    % order with no findall record/copy per word. NB named helper, not a yall
    % lambda (the F-L1 rationale, alpha_chars/2 below) - this is the gated
    % load_inf path.
    convlist(line_word, Lines, Ws0),
    sort(Ws0, Words),                       % dedupe + a stable, deterministic order
    map_list_to_pairs(length, Words, LPairs),
    keysort(LPairs, LSorted),
    group_pairs_by_key(LSorted, LGroups),
    list_to_assoc(LGroups, DictByLen),
    build_index(DictByLen, Index).

read_file_lines(File, Lines) :-
    read_file_to_string(File, Str, []),
    split_string(Str, "\n", "\r \t", Parts),
    exclude(==(""), Parts, Lines).

% One dictionary line -> its normalized letters; fails (and convlist skips the
% line) when nothing alphabetic remains.
line_word(Line, Letters) :-
    normalize_word(Line, Letters),
    Letters \== [].

% NB first-order on purpose: an include([C]>>char_type(C,alpha), ...) filter
% pays the yall meta-call + lambda copy PER CHARACTER (this module never
% imports library(yall), so the lambda is not compile-expanded): ~12 inf/char
% = 19.4M of full ENABLE's 26.6M dict-load inferences (P-F1 attribution).
% alpha_chars/2 makes the identical char_type(C, alpha) decision per char at
% ~2.2 inf/char (F-L1: load_inf 26.60M -> 10.72M at 172k, output-identical).
normalize_word(Line, Letters) :-
    string_upper(Line, U), string_chars(U, Cs),
    alpha_chars(Cs, Letters).

alpha_chars([], []).
alpha_chars([C|Cs], Letters) :-
    (   char_type(C, alpha)
    ->  Letters = [C|Rest]
    ;   Letters = Rest
    ),
    alpha_chars(Cs, Rest).

build_index(DictByLen, Index) :-
    findall(k(Len, P, Ch)-Idx,
            ( gen_assoc(Len, DictByLen, Words),
              nth0(Idx, Words, W), nth0(P, W, Ch) ),
            Triples),
    keysort(Triples, Sorted),
    group_pairs_by_key(Sorted, Grouped),
    % NB named helper, not a yall lambda (the F-L1 rationale, alpha_chars/2
    % above): un-imported `>>` is interpreted - a meta-call + lambda copy per
    % index key group, on the load_inf-gated dict-load path.
    maplist(pair_ordset, Grouped, GroupedSets),
    list_to_assoc(GroupedSets, Index).

pair_ordset(K-Is, K-Set) :- list_to_ord_set(Is, Set).

% The slot's bound positions as ascending P-V pairs. NB first-order on purpose
% (the F-L1 rationale, alpha_chars/2 above): the previous
% findall(P-V, (nth0(P, Vars, V), nonvar(V)), Bound) paid the findall
% record/collect/copy machinery plus nth0's enumerate-by-backtracking
% choicepoint per cell, at EVERY recount at EVERY search node - the hottest
% kernel in the engine. This positional walk makes the identical nonvar/1
% decision per cell; Bound's content and order (ascending P, V shared not
% copied - always a ground char atom here) are identical by construction.
% Shared by slot_bucket/5 and the masks kernel of candidate_count/5.
bound_positions(Vars, Bound) :-
    bound_positions(Vars, 0, Bound).

bound_positions([], _, []).
bound_positions([V|Vs], P, Bound) :-
    (   nonvar(V)
    ->  Bound = [P-V|Rest]
    ;   Bound = Rest
    ),
    P1 is P + 1,
    bound_positions(Vs, P1, Rest).

% Resolve a slot to its length-bucket Words plus a selector Sel over that bucket:
% `all` (no cell bound yet - every word of the length is a candidate) or
% idx(Indices) (the ordset of bucket indices matching the bound cells). Shared by
% candidates/4 (which materializes the words) and candidate_count/4 (which needs
% only the size), so the Bound + index-intersection work is written once.
slot_bucket(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    bound_positions(Vars, Bound),
    ( Bound == []
    ->  Sel = all
    ;   index_intersection(Bound, Len, Index, Indices), Sel = idx(Indices)
    ).

% Candidate words for a slot, given its currently-bound positions.
candidates(Vars, DictByLen, Index, Cands) :-
    slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> Cands = Words
    ; Sel = idx(Indices), maplist(nth0_of(Words), Indices, Cands)
    ).

nth0_of(Words, I, W) :- nth0(I, Words, W).

% Number of candidate words for a slot WITHOUT materializing the word list (P3).
% The MRV metric only needs the size, so count the matching index ordset (or the
% whole length bucket) directly and skip the per-slot maplist(nth0_of...) word
% build that select_mrv/6 otherwise ran for every unfilled slot at every search
% node. |Indices| == |maplist(nth0_of(Words), Indices)| by construction (every
% index into Words is valid), so the count is identical to length(Cands).
%
% F-H2 (bitset counting via artifact v2). candidate_count/5 dispatches on a
% threaded Masks context (see the search below): `none` keeps the ordset kernel
% (raw text mode and every non-artifact / v1-shaped path - byte-identical); a
% masks(MaskAssoc) context (present ONLY in a v2 artifact) counts the SAME slot
% with bignum `/\` chains + popcount over per-(len,pos,char) bit masks. The masks
% are derived FROM these very ordsets at --save-index time (bit i == bucket index
% i), so popcount(chain) == |ord_intersection(chain)| for every pattern by
% construction. Enumeration (candidates/4) stays on ordsets everywhere; only the
% count seam changes. candidate_count/4 is the ordset reference retained for the
% white-box tests (P3) and select_mrv/6.
candidate_count(Vars, DictByLen, Index, Count) :-
    candidate_count(none, Vars, DictByLen, Index, Count).

% The Masks context is ARGUMENT 1 so SWI's first-argument indexing dispatches
% `none` vs masks(_) deterministically - no choicepoint, no cut (on this SWI a
% non-first-argument scan does NOT prune the choicepoint, so the old arg-4
% shape needed a leading cut and still paid the CP-push per call). The clauses
% are disjoint by construction, and the reference selector select_mrv/6 stays
% choicepoint-free (P13).
%
% ORDSET kernel (Masks == none): identical work to the pre-F-H2 candidate_count/4.
candidate_count(none, Vars, DictByLen, Index, Count) :-
    slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> length(Words, Count)
    ; Sel = idx(Indices), length(Indices, Count)
    ).
% BIGNUM kernel (Masks == masks(MaskAssoc), v2 artifact mode). Mirrors
% slot_bucket's branch structure EXACTLY - the `all` branch (no cell bound) still
% counts the whole length bucket; the idx branch AND-folds the bound cells' masks
% and popcounts. No ord_intersection is computed on this path (that is the win).
candidate_count(masks(MaskAssoc), Vars, DictByLen, _Index, Count) :-
    length(Vars, Len),
    bound_positions(Vars, Bound),
    ( Bound == []
    ->  ( get_assoc(Len, DictByLen, Words) -> length(Words, Count) ; Count = 0 )
    ;   mask_count(Bound, Len, MaskAssoc, Count)
    ).

% popcount(AND over the bound cells' masks). Mirrors index_intersection/4: the
% first bound cell seeds the accumulator, the rest AND into it; an absent key is
% mask 0 (mirrors index_set's `S = []`), so a dead cell zeroes the chain and the
% count is 0 - exactly as ord_intersection with [] yields the empty set.
mask_count([P0-Ch0|Rest], Len, MaskAssoc, Count) :-
    mask_set(Len, P0, Ch0, MaskAssoc, M0),
    foldl(mask_and(Len, MaskAssoc), Rest, M0, MAnd),
    Count is popcount(MAnd).
mask_and(Len, MaskAssoc, P-Ch, Acc, Acc1) :-
    mask_set(Len, P, Ch, MaskAssoc, M), Acc1 is Acc /\ M.
mask_set(Len, P, Ch, MaskAssoc, M) :-
    ( get_assoc(k(Len, P, Ch), MaskAssoc, M0) -> M = M0 ; M = 0 ).

index_intersection([P-Ch|Rest], Len, Index, Indices) :-
    index_set(Len, P, Ch, Index, S0),
    foldl(index_intersect(Len, Index), Rest, S0, Indices).
index_intersect(Len, Index, P-Ch, Acc, Acc1) :-
    index_set(Len, P, Ch, Index, S), ord_intersection(Acc, S, Acc1).
index_set(Len, P, Ch, Index, S) :-
    ( get_assoc(k(Len, P, Ch), Index, S0) -> S = S0 ; S = [] ).


% --- the MRV backtracking search ---------------------------------------------
% Fill all slots. At each step take the unfilled slot with the FEWEST matching
% candidates (most-constrained-first; ties broken by lowest start cell, for
% determinism), try each candidate (in dictionary order), unify it into the grid
% (binding the crossing cells), and recurse. Used answers are not repeated.
% F-H1 (incremental candidate counts). A slot's candidate count is a pure
% function of its cells' binding state, and placing a word binds ONLY the
% previously-free cells of the placed slot. So instead of recounting EVERY
% unfilled slot at EVERY node (select_mrv's per-node full recount, measured at
% 79.5-85.6% of search_inf in P-F1), carry the per-slot counts as backtrack-
% restored threaded state: after a placement recount EXACTLY the slots that
% cross a NEWLY-BOUND cell of the placed word; every other count carries over
% unchanged. Backtracking restores the previous counts for free (pure threaded
% state - the caller's structure is untouched, and the carried terms are built
% by unification, never findall/copy, so the shared crossing variables survive).
%
% EQUIVALENCE (the tree is preserved node-for-node):
%  (1) COUNTS EXACT BY CONSTRUCTION. The initial counts are a full candidate_count
%      per slot; after each placement, a slot whose cells did not change keeps its
%      exact count (candidate_count is a pure function of the cells' binding state)
%      and a slot crossing a newly-bound cell is recounted with the same
%      candidate_count/4. Counts NEVER read Used (the \+ memberchk(Word, Used)
%      filter is try-time only, below), so they are exact regardless of which
%      words are already placed.
%  (2) SAME SELECTION. select_min_count picks the minimum c(Count, Start, Dir) in
%      standard order - byte-identical to select_mrv's sort(0, @=<, ...) + head.
%      Start+Dir uniquely identify a slot, so this is a total order: no genuine
%      ties, and the winner is independent of the carried list's order.
%  (3) SAME WINNER MATERIALIZATION. The winner's candidate list is built by the
%      same candidates/4.
%  (4) COMPLETED SLOTS STAY IN THE SET. A slot fully bound by crossings (count 0
%      or 1) is NOT dropped, so a 0-count dead slot is still selected first (it
%      has the lowest count) and fails the branch exactly as today.
% Masks is the F-H2 counting context threaded down the whole search: `none` (raw
% text / v1 artifact - ordset kernel, byte-identical to pre-F-H2) or
% masks(MaskAssoc) (a v2 artifact - bignum popcount kernel). Threading a ground
% atom/compound through the recursion adds ZERO inferences (SWI clause indexing
% dispatches candidate_count/5 with no penalty), so the `none` path is
% inference-identical to the pre-F-H2 engine.
fill_search(Slots, DictByLen, Index, Masks, Used) :-
    maplist(slot_with_count(DictByLen, Index, Masks), Slots, Counted),
    fill_search_inc(Counted, DictByLen, Index, Masks, Used).

% Pair each slot with its current candidate count, keeping the ORIGINAL slot
% term (shared cell variables intact - no reconstruction, no copy). This is the
% one full per-slot recount, paid once at the root; thereafter it is incremental.
slot_with_count(DictByLen, Index, Masks, Slot, cnt(Count, Slot)) :-
    Slot = slot(_, _, _, Vars),
    candidate_count(Masks, Vars, DictByLen, Index, Count).

% Counted is a list of cnt(Count, Slot). Select the min-count slot, place a word,
% recount only the crossings the placement disturbed, recurse on the rest.
% The []/[_|_] head pair is SWI's two-clause special case (jitindex §2.17):
% clause selection is deterministic with no choicepoint and no cut, and the
% recursive clause never pays a failed []-head unification per node.
fill_search_inc([], _DictByLen, _Index, _Masks, _Used).
fill_search_inc([C0|Cs], DictByLen, Index, Masks, Used) :-
    select_min_count([C0|Cs], cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    newly_bound_cells(BestCells, BestVars, NewCells),   % free cells, PRE-placement
    candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    BestVars = Word,                      % unify into the shared cell variables
    recount_crossing(Rest, NewCells, DictByLen, Index, Masks, Rest1),
    fill_search_inc(Rest1, DictByLen, Index, Masks, [Word|Used]).

% Winner = the cnt/2 with the minimum c(Count, Start, Dir) in standard order
% (Count, then Start, then Dir; all ground). This reproduces select_mrv's
% sort(0, @=<, ...)+head+once(select) exactly. Rest is the other cnt/2s with
% their order + term-sharing preserved (mirrors once(select(Best, Slots, Rest))).
select_min_count([H|T], Min, Rest) :-
    min_count_walk(T, H, Min),
    Min = cnt(_, slot(MStart, MDir, _, _)),
    remove_slot(MStart, MDir, [H|T], Rest).

min_count_walk([], Best, Best).
min_count_walk([X|Xs], Best0, Best) :-
    ( count_le(Best0, X) -> Best1 = Best0 ; Best1 = X ),
    min_count_walk(Xs, Best1, Best).

% Best0 stays iff its (Count, Start, Dir) key is @=< X's - the same key order
% select_mrv's sort/4 imposes on the c/3 terms.
count_le(cnt(C0, slot(S0, D0, _, _)), cnt(C1, slot(S1, D1, _, _))) :-
    c(C0, S0, D0) @=< c(C1, S1, D1).

% Drop the winning slot by its ground (Start, Dir) key, compared with == so no
% variable is bound. Preserves the order and term-sharing of the remaining cnt/2s.
remove_slot(MStart, MDir, [cnt(_, slot(S, D, _, _))|T], Rest) :-
    S == MStart, D == MDir, !, Rest = T.
remove_slot(MStart, MDir, [X|T], [X|Rest]) :-
    remove_slot(MStart, MDir, T, Rest).

% The placed slot's cells that were still FREE before placement - exactly the
% cells this placement newly binds. Its already-bound cells are unchanged, so
% crossings through them cannot change count. Cells and Vars are positionally
% aligned (slot_vars/3). NewCells is a plain list of ground cell NUMBERS: the
% crossing test is by cell identity, so no cell variable is ever copied.
newly_bound_cells([], [], []).
newly_bound_cells([Cell|Cs], [Var|Vs], New) :-
    ( var(Var) -> New = [Cell|New1] ; New = New1 ),
    newly_bound_cells(Cs, Vs, New1).

% Recount EXACTLY the carried slots that cross a newly-bound cell; carry every
% other count unchanged. cnt(Count, Slot) keeps the ORIGINAL Slot term (shared
% cell variables survive - no findall/copy). A crossing is a shared cell NUMBER
% (ground), so detection needs no variable comparison. Both Cells and NewCells
% are strictly ascending cell numbers by construction - Cells is a contiguous
% grid_run row/column run (stockgrid: numlist / col_cells + grab_run), and
% NewCells is an order-preserving filter of the placed slot's Cells
% (newly_bound_cells/3) - i.e. both are ordsets, so the shared-cell test is
% ord_intersect/2: one O(|Cells|+|NewCells|) merge walk, choicepoint-free,
% instead of a memberchk scan of NewCells per cell.
recount_crossing([], _NewCells, _DictByLen, _Index, _Masks, []).
recount_crossing([cnt(C0, Slot)|T], NewCells, DictByLen, Index, Masks, [cnt(C1, Slot)|T1]) :-
    Slot = slot(_, _, Cells, Vars),
    ( ord_intersect(Cells, NewCells)
    ->  candidate_count(Masks, Vars, DictByLen, Index, C1)
    ;   C1 = C0 ),
    recount_crossing(T, NewCells, DictByLen, Index, Masks, T1).

% select_mrv/6 is RETAINED as the reference selector: the white-box tests
% (tests/fill.plt R6/P13) call it directly, and it documents the exact selection
% rule fill_search_inc reproduces incrementally. It is no longer on the engine's
% hot path (fill_search_inc drives the search).
% Pick the slot with the fewest current candidates (>=0); deterministic
% tie-break: lowest start cell, then direction. Recovering the slot by BOTH
% start AND direction is load-bearing when a cell begins an across and a down
% slot (e.g. cell 1): the slot whose count was the minimum must be the one
% expanded, not whichever shares the start and appears first in the list (R6).
select_mrv(Slots, DictByLen, Index, Best, Rest, BestCands) :-
    maplist(slot_candidate_count(DictByLen, Index), Slots, Counted),
    sort(0, @=<, Counted, [c(_, BestStart, BestDir)|_]),
    Best = slot(BestStart, BestDir, _, _),
    once(select(Best, Slots, Rest)),   % Start+Dir is unique, so prune the CP (P13)
    Best = slot(_, _, _, Vars),
    candidates(Vars, DictByLen, Index, BestCands).

% c(Count, Start, Dir) - sorted by Count, then Start, then Dir gives
% most-constrained-first deterministically AND uniquely identifies the slot
% (start+dir), so select_mrv recovers exactly the slot the count was computed for.
slot_candidate_count(DictByLen, Index, slot(Start, Dir, _, Vars), c(Count, Start, Dir)) :-
    candidate_count(Vars, DictByLen, Index, Count).


% --- emit the filled layout --------------------------------------------------
slots_to_layout(Slots, Numbered, InputWords) :-
    maplist(slot_to_word, Slots, Placed),
    % No once/1: assign_clue_numbers/2 is deterministic at source (X6.B4;
    % probe re-verified on this call shape) - the wrap was defensive (C23).
    assign_clue_numbers(Placed, Numbered),
    % maplist, not findall+member (C35): a deterministic 1:1 projection in the
    % same order, with no findall copy per entry.
    maplist(pw_answer_entry, Placed, InputWords).

pw_answer_entry(PW, [A]) :- pw_answer(PW, A).

slot_to_word(slot(Start, Dir, Cells, Vars),
             pw(A, Vars, Cells, Dir, Len, Start, _End, _Num)) :-
    length(Vars, Len), atom_chars(A, Vars).


% --- entry point -------------------------------------------------------------
% seed_used(+AllSlots, +SearchSlots, -Used0): the search's INITIAL no-duplicate
% set = the seed answers. The seeded slots are exactly AllSlots minus
% SearchSlots (fill_prepare's exclude/3 keeps the same slot terms, so
% subtract/3's memberchk recovers them; Start+Dir are ground and unique, so no
% free cell variable can be bound by the scan). A seeded slot's Vars are fully
% bound by apply_seed and ARE the placed word in the exact shape Used holds (a
% list of char atoms), so fill_search_inc's `\+ memberchk(Word, Used)` now
% dedups searched slots against the seed pins too (regression: a seed answer
% the dictionary could also place was re-placed in a searched slot and the
% duplicate blew up the emit; see docs/plans/fill-seed-pin-crash-fix.md).
% Do NOT collapse Vars to an atom (atom_chars): Used holds char LISTS and an
% atom never matches, silently disabling the dedup. With no seeds,
% SearchSlots == AllSlots: the guard below yields Used0 == [] without the
% subtract scan, so the benched fill_attempt/8 window stays neutral on unseeded
% rungs (the seed_used cost lands only on seeded rungs). subtract(X, X, R) would
% also give R == [], so the fast path is byte-identical, only cheaper.
seed_used(AllSlots, SearchSlots, Used0) :-
    (   AllSlots == SearchSlots
    ->  Used0 = []
    ;   subtract(AllSlots, SearchSlots, SeededSlots),
        maplist(slot_word, SeededSlots, Used0)
    ).

slot_word(slot(_, _, _, Vars), Vars).

% Outcome: filled | infeasible | not_proven. seeds/dict applied before search;
% a seed clash / no-slot / duplicate-answer seed throws (reported before
% searching). SearchSlots are the slots the engine fills (all slots minus seed
% pins); AllSlots are emitted (so seed pins appear in the layout too).
fill_attempt(SearchSlots, AllSlots, DictByLen, Index, Outcome, Numbered, InputWords) :-
    fill_budget(B),
    fill_attempt(SearchSlots, AllSlots, DictByLen, Index, B, Outcome, Numbered, InputWords).

% Budget-explicit form: Budget is the inference budget. The default /7 reads
% fill_budget/1; passing a tiny Budget drives the AC-FILL-1 "not proven within
% budget" path (tests/fill.plt).
fill_attempt(SearchSlots, AllSlots, DictByLen, Index, Budget, Outcome, Numbered, InputWords) :-
    % No catch/3: see construct_one/7 in arrange.pl - call_with_inference_limit/3
    % binds Limit = inference_limit_exceeded on the budget path and only re-throws
    % genuine errors. Infeasibility is a search FAILURE (R == exhausted), so a
    % thrown error is a real bug and must surface, not be masked as infeasible.
    % Ordset entry (Masks == none): fill_search/5 with `none` does byte-identical
    % work to the pre-F-H2 fill_search/4, so this bench/test seam is inference-
    % identical. The product artifact path uses fill_attempt_masked/9 below.
    % Used0 is computed BEFORE call_with_inference_limit/3 so it is never
    % charged against the search Budget.
    seed_used(AllSlots, SearchSlots, Used0),
    call_with_inference_limit(
        ( fill_search(SearchSlots, DictByLen, Index, none, Used0) -> R = ok ; R = exhausted ),
        Budget, Limit),
    (   Limit == inference_limit_exceeded
    ->  Outcome = not_proven, Numbered = [], InputWords = []
    ;   R == ok
    ->  Outcome = filled, slots_to_layout(AllSlots, Numbered, InputWords)
    ;   Outcome = infeasible, Numbered = [], InputWords = []
    ).

% Masks-carrying twin of fill_attempt/8 (product path only - NOT the gated bench
% seam). Identical control flow; the only difference is the threaded Masks (none
% for the raw CLI, masks(_) for a v2 artifact). Kept separate so fill_attempt/8
% above stays byte-identical for the ratchet.
fill_attempt_masked(SearchSlots, AllSlots, DictByLen, Index, Budget, Masks,
                    Outcome, Numbered, InputWords) :-
    seed_used(AllSlots, SearchSlots, Used0),
    call_with_inference_limit(
        ( fill_search(SearchSlots, DictByLen, Index, Masks, Used0) -> R = ok ; R = exhausted ),
        Budget, Limit),
    (   Limit == inference_limit_exceeded
    ->  Outcome = not_proven, Numbered = [], InputWords = []
    ;   R == ok
    ->  Outcome = filled, slots_to_layout(AllSlots, Numbered, InputWords)
    ;   Outcome = infeasible, Numbered = [], InputWords = []
    ).

%!  fill_solve(+GridFile:atom, +SeedFileOrNone, +DictFile:atom,
%!             +SizeMode:oneof([fixed,max])) is semidet.
%
%   The `fill` CLI seam: derive the slots of the stock grid GridFile, pin the
%   seeds of SeedFileOrNone (a §6.6 fragment file, or the atom `none`), fill
%   every remaining slot from the text dictionary DictFile, and emit the
%   filled canonical layout on stdout framed by SizeMode. FAILS (no stdout)
%   on any non-filled outcome - not_proven (budget) or infeasible - after
%   reporting the unfillable slot(s) on stderr (INV-3, never silent). Throws
%   on a malformed grid/seed file, an unmatchable seed, or two seeds pinning
%   the same answer (fill_seed_duplicate).
fill_solve(GridFile, SeedFileOrNone, DictFile, SizeMode) :-
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    load_dict(DictFile, DictByLen, Index),
    fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, none, SizeMode).

%!  fill_solve_index(+GridFile:atom, +SeedFileOrNone, +IndexFile:atom,
%!                   +DictFileOrNone, +SizeMode:oneof([fixed,max])) is semidet.
%
%   Artifact-consuming twin of fill_solve/4 (F-L2): load a prebuilt,
%   verified index artifact (fill_save_index/2,3) instead of parsing a text
%   dictionary, then fill IDENTICALLY - the filled layout is byte-for-byte
%   identical to the raw path (proven by the identity oracle in both modes).
%   DictFileOrNone: an atom path checks the artifact's embedded SHA-256;
%   `none` skips that check (version + SWI are always checked); any mismatch
%   throws a clear rebuild error. Same failure contract as fill_solve/4.

% Grid/seed derivation, search, and emit are the shared body below.
% A v2 artifact carries precomputed bitset Masks; F-H2's bignum counting kernel
% runs iff Masks == masks(_). A `none` here (no masks in the artifact) transparently
% falls back to the ordset kernel, so this entry is correct for any artifact shape
% the loader accepts.
fill_solve_index(GridFile, SeedFileOrNone, IndexFile, DictFileOrNone, SizeMode) :-
    fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, Masks),
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, Masks, SizeMode).

% Grid + seed derivation, shared by both entry points (identical to the raw
% path's front matter). Slots are ALL slots (emitted); SearchSlots exclude seed
% pins (what the engine fills).
fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots) :-
    fill_grid(GridFile, Size, Slots, _CellVar),
    ( SeedFileOrNone == none
    ->  SeededKeys = []
    ;   apply_seeds(SeedFileOrNone, Slots, SeededKeys) ),
    exclude(seeded_slot(SeededKeys), Slots, SearchSlots).

% Search + emit, shared by both entry points. Fails (no stdout) on any
% non-filled outcome (INV-3): report the unfillable slot(s), never silent.
% Masks threads F-H2's counting context: `none` runs the ordset kernel (raw CLI
% path - byte-identical to pre-F-H2, exercised by the identity oracle), masks(_)
% runs the bignum kernel (v2 artifact). The raw branch calls the unchanged
% fill_attempt/7 so the raw CLI path is untouched; only the artifact branch takes
% the masked core.
fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, Masks, SizeMode) :-
    (   Masks == none
    ->  fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, InputWords)
    ;   fill_budget(B),
        fill_attempt_masked(SearchSlots, Slots, DictByLen, Index, B, Masks,
                            Outcome, Numbered, InputWords)
    ),
    (   Outcome == filled
    ->  emit_fill(Numbered, InputWords, Size, SizeMode),
        length(Numbered, NS),
        verbose_report("fill: grid ~wx~w, filled ~w slots~n", [Size, Size, NS])
    ;   fill_report_failure(Outcome, SearchSlots, DictByLen, Index, Size),
        fail
    ).

% The fill emit boundary. First re-assert the unique-answers invariant the
% solve path checks in crossword/4 (fill bypasses crossword/4, so without this
% it never ran here): any duplicate that still reaches emit reports as the
% clean hooked duplicate_answer error (exit 1), never the raw
% domain_error(unique_key_pairs) from core's answer_meta_assoc/2. Defense in
% depth only - seed_used/3 + check_unique_seed_answers/1 are the real fixes,
% and this guard alone could NOT be (it would fail solvable grids the search
% must instead fill with a non-duplicate word). The check sits HERE, outside
% fill_attempt/8, so the gated bench window (call_time over fill_attempt/8)
% is untouched on every rung.
emit_fill(Numbered, InputWords, Size, SizeMode) :-
    check_unique_answers(InputWords),
    emit_fill_mode(SizeMode, Numbered, InputWords, Size).

% fixed: the exact Size x Size canonical layout (blocks as null). (`max` would
% crop, but a stock grid is already its own frame, so fixed is the norm.)
emit_fill_mode(fixed, Numbered, InputWords, Size) :-
    emit_json(Numbered, InputWords, Size).
emit_fill_mode(max, Numbered, InputWords, Size) :-
    emit_arrange(Numbered, InputWords, Size, max).

fill_report_failure(not_proven, _Slots, _D, _I, Size) :-
    format(user_error,
           "fill: not proven within budget on ~wx~w grid (search did not complete)~n",
           [Size, Size]).
fill_report_failure(infeasible, Slots, DictByLen, Index, Size) :-
    ( empty_slots(Slots, DictByLen, Index, Bad), Bad = [_|_]
    ->  format(user_error,
               "fill: no dictionary word fits slot(s) ~w on the ~wx~w grid~n",
               [Bad, Size, Size])
    ;   format(user_error,
               "fill: no complete fill exists for this grid + dictionary (+ seeds)~n",
               [])
    ).

% Slots (by start cell) that already have ZERO candidate words - genuinely
% unfillable, reportable up front. "Count is zero" is asked of the purpose-
% built counter (C36), not by materializing candidates/4's word list and
% matching [].
empty_slots(Slots, DictByLen, Index, Bad) :-
    findall(Start,
            ( member(slot(Start, _, _, Vars), Slots),
              candidate_count(Vars, DictByLen, Index, 0) ),
            Bad).


% --- precomputed index artifact (F-L2) ---------------------------------------
% The dictionary parse + index build is a PURE FUNCTION of the frozen dict file
% and dominates end-to-end latency (P-F1: dict_load is 58-84% of CLI wall on
% 10/11 rungs; F-L1 cut the parse, leaving the build_index keysort + GC as the
% inference-blind residue). This serializes the EXACT runtime structures (the
% DictByLen length buckets + the assoc Index) once, so an interactive fill loads
% them back instead of recomputing them every invocation. The raw text-dict path
% (load_dict/build_index) is untouched: this builder CALLS load_dict, and the
% loader reconstructs terms that are ==-identical to a fresh build.
%
% ARTIFACT TERM (versioned + extensible):
%   fill_index(Version, Meta, DictByLen, Index)
%     Version - integer artifact-SCHEMA version (fill_index_format_version/1). A
%               schema change bumps it; the loader refuses an unknown version.
%     Meta    - a keyed list [Key(Value), ...] carrying integrity + provenance
%               AND, since schema v2 (F-H2), OPTIONALLY the precomputed bitset
%               masks under a masks(...) key (the realized extension point the
%               F-L2 design reserved; the Version bump 1->2 makes it a clean
%               break - v1 artifacts are refused, never misread).
%       dict_sha256(Hex) - SHA-256 of the source dict file's bytes (staleness)
%       swi_version(V)   - the SWI-Prolog that built it (binary-format guard)
%       words(N)         - dictionary word count (informational)
%       source(Path)     - the dict path as given at build (informational)
%       built_epoch(E)   - build time, integer seconds (informational)
%       masks(MaskAssoc) - [v2, OPTIONAL - `--save-index --masks` only] assoc
%                          k(Len,Pos,Char) -> bignum; bit i set iff bucket index
%                          i is in Index's ordset for that key. The F-H2 counting
%                          kernel AND-folds + popcounts these instead of
%                          ord_intersection'ing the ordsets (wall win, count
%                          identical). Built ONLY here (amortized), never at
%                          load; absent -> the loader hands back Masks = none and
%                          counting stays on the ordset kernel (no size/load tax).
%
% ON-DISK FORMAT: fast_write/fast_read (library(fastrw)) - MEASURED the fastest
% candidate (F-L2 results doc): ~20x the post-F-L1 warm raw load at 172k, smaller
% on disk than the .qlf of the equivalent fact, ~0 inferences. fast_read's binary
% format is SWI-version-bound (documented), so the artifact embeds swi_version
% and the loader REFUSES a mismatch (rebuild explicitly, never silently).

% Schema version 2 (F-H2): the Meta list MAY carry a masks(MaskAssoc) key
% (precomputed bignum bitsets for the counting kernel). Masks are OPTIONAL within
% v2: the default build omits them (F-H2 follow-up finding - masks inflate the
% artifact ~32% and tax EVERY load ~+30% fast_read wall, a net end-to-end LOSS on
% load-dominated fills), and `--save-index --masks` opts in for search-heavy /
% WASM deployments where the 23-29% fill-phase win outweighs the load tax. The
% loader accepts both shapes (absent key -> ordset kernel). The version bump is a
% hard break - a v1 artifact (or any Version != 2) is refused by the loader with
% the existing rebuild message. Masks are built ONLY here (the amortized
% --save-index step - the gate probe's binding condition), never at load or fill.
fill_index_format_version(2).

%!  fill_save_index(+DictFile:atom, +OutFile:atom) is det.
%!  fill_save_index(+DictFile:atom, +OutFile:atom, +Options:list) is det.
%
%   BUILD the precomputed index artifact: load_dict the frozen DictFile and
%   serialize the exact runtime structures + integrity/provenance meta to
%   OutFile (fast_write binary). This is the one-off cost (~ current load +
%   a fast_write); the CLI's `fill --save-index FILE` seam calls it. /2 = the
%   default build (no masks). /3: Options is a keyed list; masks(true)
%   additionally derives the parallel bitset masks FROM the ordset Index (so
%   popcount agreement is by construction) and embeds them as the masks(...)
%   Meta key. Any other Options content is ignored (same permissive
%   convention as Meta itself). Throws on an unreadable dictionary file.
fill_save_index(DictFile, OutFile) :-
    fill_save_index(DictFile, OutFile, []).

fill_save_index(DictFile, OutFile, Options) :-
    load_dict(DictFile, DictByLen, Index),
    dict_word_count(DictByLen, NWords),
    fill_file_sha256(DictFile, Sha),
    fill_swi_version(Swi),
    get_time(TF), Epoch is round(TF),
    fill_index_format_version(Version),
    Meta0 = [ dict_sha256(Sha), swi_version(Swi), words(NWords),
              source(DictFile), built_epoch(Epoch) ],
    (   option(masks(true), Options)
    ->  build_masks(Index, MaskAssoc),
        append(Meta0, [masks(MaskAssoc)], Meta)
    ;   Meta = Meta0
    ),
    Artifact = fill_index(Version, Meta, DictByLen, Index),
    setup_call_cleanup(open(OutFile, write, S, [type(binary)]),
                       fast_write(S, Artifact),
                       close(S)).

% Derive the mask assoc from the ordset Index: for each k(Len,Pos,Char) key whose
% value is an ordset of bucket indices, the mask is the bignum with bit i set for
% every index i in that ordset. Because masks and buckets share index_set's
% ordering (the masks ARE the ordsets, re-encoded), popcount(mask chain) equals
% |ord_intersection(set chain)| for every pattern - the F-H2 equivalence, by
% construction. Same key set as Index, so an absent key means mask 0 (dead cell).
build_masks(Index, MaskAssoc) :-
    assoc_to_list(Index, Pairs),
    maplist(ordset_mask_pair, Pairs, MaskPairs),
    list_to_assoc(MaskPairs, MaskAssoc).
ordset_mask_pair(K-Set, K-Mask) :- ordset_to_mask(Set, 0, Mask).
ordset_to_mask([], M, M).
ordset_to_mask([I|Is], M0, M) :- M1 is M0 \/ (1 << I), ordset_to_mask(Is, M1, M).

% LOAD: read the artifact, verify schema version + SWI version (+ the dict hash
% iff a dict path is supplied), and hand back the reconstructed DictByLen +
% Index. Any integrity failure throws a clear fill_index_* error - no silent
% rebuild, no fall-through to the raw path.
% /4: the pre-F-H2 signature, retained for callers/tests that only need the
% dictionary structures (e.g. the roundtrip identity test). Discards the masks.
fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index) :-
    fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, _Masks).

% /5: also returns the F-H2 counting context - masks(MaskAssoc) if the (v2)
% artifact carries a masks key, else none (the ordset kernel still works). Every
% integrity gate (version, SWI, hash) runs BEFORE masks are extracted, so a bad
% artifact throws exactly as before; masks come off the happy path only.
fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, Masks) :-
    ( exists_file(IndexFile) -> true
    ; throw(error(fill_index_missing(IndexFile), _)) ),
    % Keep the root cause (C33): the context slot carries the original error
    % term (fastrw version mismatch vs permissions vs truncation), so it shows
    % under verbose error printing instead of being swallowed.
    catch(fill_read_artifact(IndexFile, Artifact),
          ReadErr,
          throw(error(fill_index_unreadable(IndexFile), context(_, ReadErr)))),
    ( Artifact = fill_index(Version, Meta, DictByLen0, Index0)
    -> true
    ;  throw(error(fill_index_malformed(IndexFile), _)) ),
    fill_index_format_version(Want),
    ( Version == Want -> true
    ; throw(error(fill_index_version(Version, Want), _)) ),
    fill_swi_version(CurSwi),
    ( option(swi_version(ArtSwi), Meta), ArtSwi == CurSwi -> true
    ; option(swi_version(ArtSwi2), Meta, unknown),
      throw(error(fill_index_swi(ArtSwi2, CurSwi), _)) ),
    ( DictFileOrNone == none
    -> true
    ;  fill_file_sha256(DictFileOrNone, CurSha),
       ( option(dict_sha256(ArtSha), Meta), ArtSha == CurSha -> true
       ; option(dict_sha256(ArtSha2), Meta, unknown),
         throw(error(fill_index_hash(DictFileOrNone, ArtSha2, CurSha), _)) ) ),
    DictByLen = DictByLen0, Index = Index0,
    ( option(masks(MaskAssoc), Meta) -> Masks = masks(MaskAssoc) ; Masks = none ).

fill_read_artifact(IndexFile, Artifact) :-
    setup_call_cleanup(open(IndexFile, read, S, [type(binary)]),
                       fast_read(S, Artifact),
                       close(S)).

% Word count over the length buckets (same fold run_fill's metadata uses).
dict_word_count(DictByLen, N) :-
    assoc_to_values(DictByLen, Buckets),
    foldl(bucket_len, Buckets, 0, N).
bucket_len(B, A0, A1) :- length(B, L), A1 is A0 + L.

% SHA-256 of the file's raw bytes - a self-consistent build/verify fingerprint
% that also matches coreutils sha256sum (encoding(octet) hashes each byte as-is).
fill_file_sha256(File, Hex) :-
    read_file_to_string(File, Bytes, [encoding(octet)]),
    sha_hash(Bytes, Digest, [algorithm(sha256), encoding(octet)]),
    hash_atom(Digest, Hex).

% "Ma.Mi.Pa" - the same version string run_fill.pl records, so the artifact's
% guard and the bench provenance agree.
fill_swi_version(Ver) :-
    current_prolog_flag(version_data, V),
    ( V = swi(Ma, Mi, Pa, _) -> format(atom(Ver), '~d.~d.~d', [Ma, Mi, Pa]) ; Ver = V ).


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(fill_seed_no_slot(A)) -->
    [ 'fill: seed ~q does not match any slot of the grid (check its cells/direction)'-[A] ].
prolog:error_message(fill_seed_clash(A)) -->
    [ 'fill: seed ~q clashes with another seed at a shared cell'-[A] ].
prolog:error_message(fill_seed_duplicate(A)) -->
    [ 'fill: seed ~q is pinned more than once; answers must be unique'-[A] ].
% F-L2 index-artifact integrity failures (refuse, never silently rebuild).
prolog:error_message(fill_index_missing(F)) -->
    [ 'fill: index artifact ~w not found (build one with `fill --dict DICT --save-index ~w`)'-[F, F] ].
prolog:error_message(fill_index_unreadable(F)) -->
    [ 'fill: index artifact ~w could not be read (corrupt, or built by a different SWI-Prolog; rebuild with --save-index)'-[F] ].
prolog:error_message(fill_index_malformed(F)) -->
    [ 'fill: index artifact ~w is not a fill_index/4 artifact (rebuild with --save-index)'-[F] ].
prolog:error_message(fill_index_version(Got, Want)) -->
    [ 'fill: index artifact schema version ~w is not supported (this build reads version ~w); rebuild with --save-index'-[Got, Want] ].
prolog:error_message(fill_index_swi(Art, Cur)) -->
    [ 'fill: index artifact was built by SWI-Prolog ~w but this is ~w (the binary format is version-bound); rebuild with --save-index'-[Art, Cur] ].
prolog:error_message(fill_index_hash(File, Art, Cur)) -->
    [ 'fill: index artifact does not match --dict ~w (artifact dict SHA-256 ~w, file ~w); the dictionary changed - rebuild with --save-index'-[File, Art, Cur] ].
