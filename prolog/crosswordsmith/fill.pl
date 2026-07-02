% fill.pl - Flavour-B grid-first, open-dictionary auto-fill (design-spec §8.4;
% OD-1..4 resolved by DP-1/DP-2). Take a pre-validated legal blocked grid (a §8.3
% stock-grid mask), derive its slots, and fill every slot with a dictionary word
% subject to crossing constraints, with the user's words pinned as seeds (the
% §6.6 fragment primitive). Output is the canonical layout, so it composes with
% lint/export. Consult AFTER core.pl, metrics.pl, lint.pl, stockgrid.pl and
% arrange.pl (it reuses load_fragment/3).
%
% The grid model: each white cell is one logical VARIABLE, shared between its
% across and down slot. Assigning a word to a slot unifies the slot's variable
% list with the word's letters; crossings are then consistent by construction
% (the shared variable), and dead-ends backtrack. MRV ordering (fewest matching
% candidates first) + a node/inference budget make it deterministic and bounded.

:- use_module(library(http/json)).
:- use_module(library(apply)).
:- use_module(library(lists)).
:- use_module(library(ordsets)).
:- use_module(library(assoc)).
:- use_module(library(pairs)).

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

init_cell_vars(WhiteCells, Assoc) :-
    empty_assoc(A0), foldl(add_cell_var, WhiteCells, A0, Assoc).
add_cell_var(C, AIn, AOut) :- put_assoc(C, AIn, _FreshVar, AOut).

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
    foldl(apply_seed(Slots), Frags, [], SeededKeys).

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
    findall(W, ( member(L, Lines), normalize_word(L, W), W \== [] ), Ws0),
    sort(Ws0, Words),                       % dedupe + a stable, deterministic order
    map_list_to_pairs(length, Words, LPairs),
    keysort(LPairs, LSorted),
    group_pairs_by_key(LSorted, LGroups),
    list_to_assoc(LGroups, DictByLen),
    build_index(DictByLen, Index).

read_file_lines(File, Lines) :-
    setup_call_cleanup(open(File, read, S), read_string(S, _, Str), close(S)),
    split_string(Str, "\n", "\r \t", Parts),
    exclude(==(""), Parts, Lines).

normalize_word(Line, Letters) :-
    string_upper(Line, U), string_chars(U, Cs),
    include([C]>>char_type(C, alpha), Cs, Letters).

build_index(DictByLen, Index) :-
    findall(k(Len, P, Ch)-Idx,
            ( gen_assoc(Len, DictByLen, Words),
              nth0(Idx, Words, W), nth0(P, W, Ch) ),
            Triples),
    keysort(Triples, Sorted),
    group_pairs_by_key(Sorted, Grouped),
    maplist([K-Is, K-Set]>>list_to_ord_set(Is, Set), Grouped, GroupedSets),
    list_to_assoc(GroupedSets, Index).

% Resolve a slot to its length-bucket Words plus a selector Sel over that bucket:
% `all` (no cell bound yet - every word of the length is a candidate) or
% idx(Indices) (the ordset of bucket indices matching the bound cells). Shared by
% candidates/4 (which materializes the words) and candidate_count/4 (which needs
% only the size), so the Bound + index-intersection work is written once.
slot_bucket(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    findall(P-V, ( nth0(P, Vars, V), nonvar(V) ), Bound),
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
candidate_count(Vars, DictByLen, Index, Count) :-
    slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> length(Words, Count)
    ; Sel = idx(Indices), length(Indices, Count)
    ).

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
fill_search([], _DictByLen, _Index, _Used) :- !.
fill_search(Slots, DictByLen, Index, Used) :-
    select_mrv(Slots, DictByLen, Index, slot(_, _, _, Vars), Rest, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    Vars = Word,                          % unify into the shared cell variables
    fill_search(Rest, DictByLen, Index, [Word|Used]).

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
    once(assign_clue_numbers(Placed, Numbered)),
    findall([A], member(word{answer:A, letters:_, cells:_, dir:_, len:_, start:_}, Placed),
            InputWords).

slot_to_word(slot(Start, Dir, Cells, Vars),
             word{answer:A, letters:Vars, cells:Cells, dir:Dir, len:Len, start:Start}) :-
    length(Vars, Len), atom_chars(A, Vars).


% --- entry point -------------------------------------------------------------
% Outcome: filled | infeasible | not_proven. seeds/dict applied before search;
% a seed clash / no-slot throws (reported before searching).
% SearchSlots are the slots the engine fills (all slots minus seed pins);
% AllSlots are emitted (so seed pins appear in the layout too).
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
    call_with_inference_limit(
        ( once(fill_search(SearchSlots, DictByLen, Index, [])) -> R = ok ; R = exhausted ),
        Budget, Limit),
    (   Limit == inference_limit_exceeded
    ->  Outcome = not_proven, Numbered = [], InputWords = []
    ;   R == ok
    ->  Outcome = filled, slots_to_layout(AllSlots, Numbered, InputWords)
    ;   Outcome = infeasible, Numbered = [], InputWords = []
    ).

% Construct + emit a fill on stdout; report on stderr. Fails (no stdout) on any
% non-filled outcome (INV-3: report the unfillable slot(s), never silent).
fill_solve(GridFile, SeedFileOrNone, DictFile, SizeMode) :-
    fill_grid(GridFile, Size, Slots, _CellVar),
    ( SeedFileOrNone == none
    ->  SeededKeys = []
    ;   apply_seeds(SeedFileOrNone, Slots, SeededKeys) ),
    exclude(seeded_slot(SeededKeys), Slots, SearchSlots),
    load_dict(DictFile, DictByLen, Index),
    fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, InputWords),
    (   Outcome == filled
    ->  emit_fill(Numbered, InputWords, Size, SizeMode),
        length(Numbered, NS),
        format(user_error, "fill: grid ~wx~w, filled ~w slots~n", [Size, Size, NS])
    ;   fill_report_failure(Outcome, SearchSlots, DictByLen, Index, Size),
        fail
    ).

% fixed: the exact Size x Size canonical layout (blocks as null). (`max` would
% crop, but a stock grid is already its own frame, so fixed is the norm.)
emit_fill(Numbered, InputWords, Size, fixed) :-
    emit_json(Numbered, InputWords, Size).
emit_fill(Numbered, InputWords, Size, max) :-
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
% unfillable, reportable up front.
empty_slots(Slots, DictByLen, Index, Bad) :-
    findall(Start,
            ( member(slot(Start, _, _, Vars), Slots),
              candidates(Vars, DictByLen, Index, []) ),
            Bad).


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(fill_seed_no_slot(A)) -->
    [ 'fill: seed ~q does not match any slot of the grid (check its cells/direction)'-[A] ].
prolog:error_message(fill_seed_clash(A)) -->
    [ 'fill: seed ~q clashes with another seed at a shared cell'-[A] ].
