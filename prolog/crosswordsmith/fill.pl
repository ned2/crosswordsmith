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
% Exports: fill_solve/4 (the CLI seam) — nothing else; every other predicate
% is internal (tests reach them as crosswordsmith_fill:Pred(...)). All
% project dependencies are explicit imports below (stockgrid, metrics,
% arrange, core).

:- module(crosswordsmith_fill,
          [ fill_solve/4,
            fill_solve_index/5,
            fill_save_index/2
          ]).

:- use_module(library(http/json)).
:- use_module(library(apply)).
:- use_module(library(lists)).
:- use_module(library(ordsets)).
:- use_module(library(assoc)).
:- use_module(library(pairs)).
:- use_module(library(sha)).        % F-L2: dict-file SHA-256 (artifact integrity)
:- use_module(library(fastrw)).     % F-L2: fast_write/fast_read artifact I/O

% Slot derivation from a stock-grid mask.
:- use_module(crosswordsmith(stockgrid),
              [stockgrid_load/2, mask_white_cells/3, grid_run/4]).

% word_letters/3: the separator-stripped placement footprint of a seed answer.
:- use_module(crosswordsmith(metrics), [word_letters/3]).

% load_fragment/3: seeds arrive in the §6.6 fragment format. emit_arrange/4:
% emit_fill's `max` mode delegates the cropped emit to arrange.
:- use_module(crosswordsmith(arrange), [load_fragment/3, emit_arrange/4]).

% Numbering + the canonical JSON emit for the filled layout, plus the placed-word
% record (pw/8) accessor fill uses to recover each answer for the emit metadata.
:- use_module(crosswordsmith(core),
              [assign_clue_numbers/2, emit_json/3, verbose_report/2, pw_answer/2]).

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
fill_search(Slots, DictByLen, Index, Used) :-
    maplist(slot_with_count(DictByLen, Index), Slots, Counted),
    fill_search_inc(Counted, DictByLen, Index, Used).

% Pair each slot with its current candidate count, keeping the ORIGINAL slot
% term (shared cell variables intact - no reconstruction, no copy). This is the
% one full per-slot recount, paid once at the root; thereafter it is incremental.
slot_with_count(DictByLen, Index, Slot, cnt(Count, Slot)) :-
    Slot = slot(_, _, _, Vars),
    candidate_count(Vars, DictByLen, Index, Count).

% Counted is a list of cnt(Count, Slot). Select the min-count slot, place a word,
% recount only the crossings the placement disturbed, recurse on the rest.
fill_search_inc([], _DictByLen, _Index, _Used) :- !.
fill_search_inc(Counted, DictByLen, Index, Used) :-
    select_min_count(Counted, cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    newly_bound_cells(BestCells, BestVars, NewCells),   % free cells, PRE-placement
    candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    BestVars = Word,                      % unify into the shared cell variables
    recount_crossing(Rest, NewCells, DictByLen, Index, Rest1),
    fill_search_inc(Rest1, DictByLen, Index, [Word|Used]).

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
% (ground), so detection needs no variable comparison.
recount_crossing([], _NewCells, _DictByLen, _Index, []).
recount_crossing([cnt(C0, Slot)|T], NewCells, DictByLen, Index, [cnt(C1, Slot)|T1]) :-
    Slot = slot(_, _, Cells, Vars),
    ( shares_cell(Cells, NewCells)
    ->  candidate_count(Vars, DictByLen, Index, C1)
    ;   C1 = C0 ),
    recount_crossing(T, NewCells, DictByLen, Index, T1).

shares_cell([C|_], NewCells) :- memberchk(C, NewCells), !.
shares_cell([_|Cs], NewCells) :- shares_cell(Cs, NewCells).

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
    once(assign_clue_numbers(Placed, Numbered)),
    findall([A], ( member(PW, Placed), pw_answer(PW, A) ), InputWords).

slot_to_word(slot(Start, Dir, Cells, Vars),
             pw(A, Vars, Cells, Dir, Len, Start, _End, _Num)) :-
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
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    load_dict(DictFile, DictByLen, Index),
    fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, SizeMode).

% Artifact-consuming twin of fill_solve/4 (F-L2): load a prebuilt, verified
% index artifact instead of parsing a text dictionary, then fill IDENTICALLY.
% Grid/seed derivation, search, and emit are the shared body below, so the
% filled layout is byte-for-byte identical to the raw path (proven by the
% identity oracle in both modes). DictFileOrNone: an atom path checks the
% artifact's embedded SHA-256; `none` skips that check (version + SWI are
% always checked). fill_load_index throws a clear error on any mismatch.
fill_solve_index(GridFile, SeedFileOrNone, IndexFile, DictFileOrNone, SizeMode) :-
    fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index),
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, SizeMode).

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
fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, SizeMode) :-
    fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, InputWords),
    (   Outcome == filled
    ->  emit_fill(Numbered, InputWords, Size, SizeMode),
        length(Numbered, NS),
        verbose_report("fill: grid ~wx~w, filled ~w slots~n", [Size, Size, NS])
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
%               AND the F-H2 EXTENSION POINT: precomputed bitset masks will slot
%               in here as an added masks(...) key under a Version bump, so
%               adding them is not a format break (old readers never look for
%               the key; new readers gate on Version before using it).
%       dict_sha256(Hex) - SHA-256 of the source dict file's bytes (staleness)
%       swi_version(V)   - the SWI-Prolog that built it (binary-format guard)
%       words(N)         - dictionary word count (informational)
%       source(Path)     - the dict path as given at build (informational)
%       built_epoch(E)   - build time, integer seconds (informational)
%
% ON-DISK FORMAT: fast_write/fast_read (library(fastrw)) - MEASURED the fastest
% candidate (F-L2 results doc): ~20x the post-F-L1 warm raw load at 172k, smaller
% on disk than the .qlf of the equivalent fact, ~0 inferences. fast_read's binary
% format is SWI-version-bound (documented), so the artifact embeds swi_version
% and the loader REFUSES a mismatch (rebuild explicitly, never silently).

fill_index_format_version(1).

% BUILD: load_dict the frozen file, then serialize the exact structures + meta.
% This is the one-off cost (~ current load + a fast_write); the CLI's
% `fill --save-index FILE` seam calls it.
fill_save_index(DictFile, OutFile) :-
    load_dict(DictFile, DictByLen, Index),
    dict_word_count(DictByLen, NWords),
    fill_file_sha256(DictFile, Sha),
    fill_swi_version(Swi),
    get_time(TF), Epoch is round(TF),
    fill_index_format_version(Version),
    Meta = [ dict_sha256(Sha), swi_version(Swi), words(NWords),
             source(DictFile), built_epoch(Epoch) ],
    Artifact = fill_index(Version, Meta, DictByLen, Index),
    setup_call_cleanup(open(OutFile, write, S, [type(binary)]),
                       fast_write(S, Artifact),
                       close(S)).

% LOAD: read the artifact, verify schema version + SWI version (+ the dict hash
% iff a dict path is supplied), and hand back the reconstructed DictByLen +
% Index. Any integrity failure throws a clear fill_index_* error - no silent
% rebuild, no fall-through to the raw path.
fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index) :-
    ( exists_file(IndexFile) -> true
    ; throw(error(fill_index_missing(IndexFile), _)) ),
    catch(fill_read_artifact(IndexFile, Artifact),
          _ReadErr,
          throw(error(fill_index_unreadable(IndexFile), _))),
    ( Artifact = fill_index(Version, Meta, DictByLen0, Index0)
    -> true
    ;  throw(error(fill_index_malformed(IndexFile), _)) ),
    fill_index_format_version(Want),
    ( Version == Want -> true
    ; throw(error(fill_index_version(Version, Want), _)) ),
    fill_swi_version(CurSwi),
    ( meta_value(Meta, swi_version, ArtSwi), ArtSwi == CurSwi -> true
    ; meta_value_or(Meta, swi_version, unknown, ArtSwi2),
      throw(error(fill_index_swi(ArtSwi2, CurSwi), _)) ),
    ( DictFileOrNone == none
    -> true
    ;  fill_file_sha256(DictFileOrNone, CurSha),
       ( meta_value(Meta, dict_sha256, ArtSha), ArtSha == CurSha -> true
       ; meta_value_or(Meta, dict_sha256, unknown, ArtSha2),
         throw(error(fill_index_hash(DictFileOrNone, ArtSha2, CurSha), _)) ) ),
    DictByLen = DictByLen0, Index = Index0.

fill_read_artifact(IndexFile, Artifact) :-
    setup_call_cleanup(open(IndexFile, read, S, [type(binary)]),
                       fast_read(S, Artifact),
                       close(S)).

% Meta is a keyed list of unary Key(Value) terms (see the artifact term above).
meta_value(Meta, Key, Value) :- Probe =.. [Key, Value], memberchk(Probe, Meta).
meta_value_or(Meta, Key, _Default, Value) :- meta_value(Meta, Key, Value), !.
meta_value_or(_Meta, _Key, Default, Default).

% Word count over the length buckets (same fold run_fill's metadata uses).
dict_word_count(DictByLen, N) :-
    assoc_to_values(DictByLen, Buckets),
    foldl(bucket_len, Buckets, 0, N).
bucket_len(B, A0, A1) :- length(B, L), A1 is A0 + L.

% SHA-256 of the file's raw bytes - a self-consistent build/verify fingerprint
% that also matches coreutils sha256sum (encoding(octet) hashes each byte as-is).
fill_file_sha256(File, Hex) :-
    setup_call_cleanup(open(File, read, S, [encoding(octet)]),
                       read_string(S, _, Bytes),
                       close(S)),
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
