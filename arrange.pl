% arrange.pl - Flavour-A layout engine: the scoring oracle + the
% construct/rescore/emit driver, fragment seeding, and diverse candidates.
% See docs/arrange-implementation-plan.md and docs/design-spec.md §7.
%
% This file deliberately does NOT define `main`/initialization, so it can be
% consulted by a harness without side effects. It consults crossword.pl (which
% in turn loads quality.pl), reusing the legality core, the MRV-inc branch step,
% and the existing checking metrics verbatim. (The Phase-1.5 `gate_*`
% measurement harness that lived here was removed once the gate's descope
% decision was recorded in the implementation plan.)

:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, 'crossword.pl', CrosswordFile),
   ensure_loaded(CrosswordFile).

:- use_module(library(apply)).
:- use_module(library(aggregate)).


% ---------------------------------------------------------------------------
% Phase 1 - the scoring oracle (from-scratch capped reward over a placement).
% ---------------------------------------------------------------------------

% Default integer weights: WCap:WTail = 5:1  => epsilon = WTail/WCap = 0.2.
arrange_weights(5, 1).

% Per-word checking target: ceil(L/2) by default (== word_meets_half/2's
% threshold in quality.pl; L is the placed dict's `len`). The §7.2 reachability
% escape hatch `--check-target N` lowers the ceiling to min(ceil(L/2), N) via
% check_target_override/1 (absent => the default). min/2 means N only ever
% LOWERS the target - an N above ceil(L/2) is a no-op - matching the spec intent
% ("lowering target below ceil(L/2) where half is unreachable").
% Invariant: at most ONE check_target_override/1 fact. The CLI's set_check_target/1
% always retractall's before assertz, so check_target/2's `->` (which commits to
% the first solution) reads the single current override. Do not assertz a second
% without retracting (P17).
:- dynamic check_target_override/1.
check_target(L, T) :-
    T0 is (L + 1) // 2,
    ( check_target_override(K) -> T is min(T0, K) ; T = T0 ).

% Per-word integer reward contribution:
%   WCap * min(checked, target) + WTail * checked
% checked(w) is reused from quality.pl's word_checked_count/3 (W's cells that
% are also covered by a perpendicular placed word).
word_reward(WCap, WTail, PW, Placed, R) :-
    get_dict(len, PW, L),
    check_target(L, T),
    word_checked_count(PW, Placed, C),
    Capped is min(C, T),
    R is WCap * Capped + WTail * C.

% Total objective over a complete (or partial) placement - the from-scratch
% oracle Phase 2's incremental delta must match.
layout_reward(WCap, WTail, Placed, Reward) :-
    foldl([PW, A, A1]>>( word_reward(WCap, WTail, PW, Placed, R), A1 is A + R ),
          Placed, 0, Reward).

% How many placed words reach their cap (checked >= ceil(L/2)) - the
% reachability probe: if this is ~0, the cap is inert and the objective has
% degenerated to plain total-crossings.
cap_binding_count(Placed, Count) :-
    aggregate_all(count,
                  ( member(PW, Placed),
                    get_dict(len, PW, L), check_target(L, T),
                    word_checked_count(PW, Placed, C),
                    C >= T ),
                  Count).

% §7.2 (INV-3): when the cap binds on NO placed word, the objective has
% degenerated to plain total-crossings - report it (on stderr), pointing the
% user at --check-target as the knob to lower an unreachable target.
cap_status_note(Placed, Note) :-
    cap_binding_count(Placed, CB),
    ( CB =:= 0
    ->  Note = ' (cap inert: objective = total-crossings; tune with --check-target)'
    ;   Note = '' ).


% ---------------------------------------------------------------------------
% Phase 2-3 - strict layout (construct + rescore + emit) with size framing.
%
% Per the Phase-1.5 gate the B&B search is descoped: strict mode constructs a
% complete placement via the reused MRV-inc path over the four canonical start
% corners (cheap construction diversity), rescores each with layout_reward/4,
% and emits the best. place-all-or-fail with budget-aware 3-outcome semantics.
% ---------------------------------------------------------------------------

% Construction budget (inferences) for finding ONE complete placement per
% corner. First solutions on puzzle-shaped inputs are 5-462 ms (well under).
arrange_budget(500_000_000).

% Best-scoring complete placement over the four start corners, rescored and
% picked by reward (deterministic: reward desc, ties by start-corner order).
% Outcome: placed | not_proven (budget hit) | infeasible (search completed,
% no full placement).
arrange_best_layout(Words, GridLen, Numbered, Reward, Outcome) :-
    arrange_budget(Budget),
    arrange_best_layout(Words, GridLen, Budget, Numbered, Reward, Outcome).

% Budget-parameterized form: Budget is the per-corner inference budget. The
% default /5 reads arrange_budget/1; passing a tiny Budget drives the
% AC-ARR-1c / AC-ARR-10 "not proven within budget" path (tests/arrange.plt).
arrange_best_layout(Words, GridLen, Budget, Numbered, Reward, Outcome) :-
    arrange_weights(WCap, WTail),
    start_locs(Locs),
    construct_corners(Locs, Words, GridLen, WCap, WTail, Budget, Results),
    findall(R-P, member(ok(R, P), Results), OKs),
    (   OKs = [_|_]
    ->  sort(1, @>=, OKs, [Reward-BestPlaced|_]),
        % once/1: assign_clue_numbers leaves a spurious choicepoint (second-level
        % list indexing in add_clue_nums); arrange_best_layout is single-valued.
        once(assign_clue_numbers(BestPlaced, Numbered)),
        Outcome = placed
    ;   Numbered = [], Reward = -1,
        ( memberchk(budget, Results) -> Outcome = not_proven ; Outcome = infeasible )
    ).

% Try each start corner under ONE shared inference budget (§7.3 reads the budget
% as a single cap on the operation, not per-corner): the running budget is
% decremented by what each corner actually consumes, so a hard/infeasible input
% costs at most ~Budget total, not Budget x |corners|. On realistic inputs every
% corner completes in << Budget, so all corners still run and the
% best-of-corners selection is unchanged - only the worst case is bounded. Once
% the running budget is spent, the remaining corners are tagged `budget` without
% running. Deterministic: inference counts are reproducible (INV-2).
construct_corners([], _Words, _GL, _WC, _WT, _Budget, []).
construct_corners([Loc|Locs], Words, GL, WC, WT, Budget, [Res|Rest]) :-
    (   Budget < 1
    ->  Res = budget,
        construct_corners(Locs, Words, GL, WC, WT, Budget, Rest)
    ;   statistics(inferences, I0),
        construct_one(Loc, Words, GL, WC, WT, Budget, Res),
        statistics(inferences, I1),
        Used is I1 - I0,
        Budget1 is Budget - Used,
        construct_corners(Locs, Words, GL, WC, WT, Budget1, Rest)
    ).

% First MRV-inc complete placement from one corner, rescored, under a budget.
% Always succeeds with a Res tag: ok(Reward,Placed) | budget | exhausted.
construct_one(Loc, Words, GridLen, WCap, WTail, Budget, Res) :-
    % No catch/3 here: call_with_inference_limit/3 handles the budget itself - on
    % the limit it binds Outcome = inference_limit_exceeded and SUCCEEDS, and it
    % only re-throws GENUINE exceptions (metacall.md). Infeasibility is a search
    % FAILURE (Found == none), never an exception, so a thrown error is a real
    % bug: let it surface (main/0 reports it, exit 1) rather than be silently
    % reclassified as `exhausted`.
    call_with_inference_limit(
        ( once(find_crossword(mrv_inc, GridLen, Words, Loc, _G, P)) -> Found = P
        ; Found = none ),
        Budget, Outcome),
    (   Outcome == inference_limit_exceeded -> Res = budget
    ;   Found == none                       -> Res = exhausted
    ;   layout_reward(WCap, WTail, Found, R), Res = ok(R, Found)
    ).

% Construct the best strict layout, emit it (framed by SizeMode) on stdout,
% report on stderr. On any non-placed outcome: report + fail with no stdout
% (the json-output-spec §6.3 "no output / non-zero exit" contract).
arrange_strict_solve(Words, GridLen, SizeMode) :-
    check_unique_answers(Words),
    arrange_best_layout(Words, GridLen, Numbered, Reward, Outcome),
    (   Outcome == placed
    ->  emit_arrange(Numbered, Words, GridLen, SizeMode),
        length(Numbered, NP),
        cap_status_note(Numbered, Note),
        format(user_error, "arrange: grid ~w, mode ~w, placed ~w, reward ~w~w~n",
               [GridLen, SizeMode, NP, Reward, Note])
    ;   arrange_report_failure(Outcome, Words, GridLen),
        fail
    ).

arrange_report_failure(not_proven, _Words, GridLen) :-
    format(user_error,
           "arrange: not proven within budget on ~wx~w grid \c
(search did not complete; feasibility unresolved)~n",
           [GridLen, GridLen]).
arrange_report_failure(infeasible, Words, GridLen) :-
    (   unplaceable_words(Words, Bad), Bad = [_|_]
    ->  format(user_error,
               "arrange: cannot place all words on ~wx~w grid; \c
words with no possible crossing: ~w~n",
               [GridLen, GridLen, Bad])
    ;   format(user_error,
               "arrange: no complete placement on ~wx~w grid \c
(grid too small or no full interlock)~n",
               [GridLen, GridLen])
    ).

% A word sharing no letter with any other word can never be crossed - the
% spec's "genuinely unplaceable (no legal intersection anywhere)".
unplaceable_words(Words, Bad) :-
    findall(A,
            ( member(E, Words), E = [A|_], \+ word_shares_letter(E, Words) ),
            Bad).

word_shares_letter(Entry, Words) :-
    Entry = [A|_], word_letters(Entry, Ls, _),
    member(Other, Words), Other = [B|_], B \== A,
    word_letters(Other, OLs, _),
    shares_letter(Ls, OLs).   % short-circuits at the first shared letter (P6)


% --- Phase 4: best-effort (drop) via the greedy constructor ----------------
% Per the Phase-1.5 result, best-effort is served by quality.pl's greedy
% constructor (which drops words it cannot place), NOT a drop-branch on the
% strict DFS. Construct over the seed x start-corner sweep on the given grid,
% rescore, and pick lexicographically: most words placed, then highest reward.

% Best greedy layout on GridLen. The key score(NumPlaced, Reward) compared in
% standard term order ranks most-placed first, then highest reward.
arrange_best_effort(Words, GridLen, Numbered, Reward, NumPlaced, Dropped) :-
    arrange_weights(WCap, WTail),
    seed_candidates(Words, Seeds),
    start_locs(Locs),
    findall(score(NP, R)-pd(Placed, DroppedAnswers),
            ( member(Loc, Locs),
              member(Seed, Seeds),
              greedy_construct(Words, GridLen, Loc, Seed, Placed, DroppedEntries),
              length(Placed, NP),
              layout_reward(WCap, WTail, Placed, R),
              findall(A, member([A|_], DroppedEntries), DroppedAnswers) ),
            Results),
    Results = [_|_],
    sort(1, @>=, Results, [score(NumPlaced, Reward)-pd(BestPlaced, Dropped)|_]),
    once(assign_clue_numbers(BestPlaced, Numbered)).

% Construct the best-effort layout, emit it, report placed/dropped on stderr.
% Fails only when nothing at all is placeable (AC-ARR-2).
arrange_best_effort_solve(Words, GridLen, SizeMode) :-
    check_unique_answers(Words),
    (   arrange_best_effort(Words, GridLen, Numbered, Reward, NP, Dropped)
    ->  emit_arrange(Numbered, Words, GridLen, SizeMode),
        length(Dropped, ND),
        format(user_error,
               "arrange: grid ~w, mode ~w, placed ~w, dropped ~w ~w, reward ~w~n",
               [GridLen, SizeMode, NP, ND, Dropped, Reward])
    ;   format(user_error, "arrange: nothing placeable on ~wx~w grid~n",
               [GridLen, GridLen]),
        fail
    ).


% --- emit framing (Phase 3) ------------------------------------------------
% fixed : exactly GridLen x GridLen (the canonical emit_json/3).
% max   : crop to the tight SQUARE enclosing the content (side = max(H,W)),
%         anchored at the bbox top-left. Square keeps the single-`gridLength`
%         output schema (rectangular grids are out of scope, design-spec §3).

emit_arrange(Numbered, Words, GridLen, SizeMode) :-
    arrange_layout_dict(Numbered, Words, GridLen, SizeMode, Payload),
    current_output(Out),
    json_write_dict(Out, Payload),
    nl(Out).

% Build (without writing) the canonical layout dict for one placement, framed by
% SizeMode. Splitting dict-construction from writing lets --candidates emit an
% array of these (emit_candidates/4) while single emit writes exactly one.
%   fixed : the canonical N x N dict (identical to emit_json/3's payload).
%   max   : the tight enclosing-square crop dict.
arrange_layout_dict(Numbered, Words, GridLen, fixed,
                    _{gridLength: GridLen, grid: Rows, words: WordObjs}) :-
    build_grid_rows(Numbered, GridLen, Rows),
    build_words(Numbered, Words, GridLen, WordObjs).
arrange_layout_dict(Numbered, Words, GridLen, max, Dict) :-
    placed_bbox(Numbered, GridLen, bbox(MinR, MaxR, MinC, MaxC), _Area),
    H is MaxR - MinR + 1, W is MaxC - MinC + 1, S is max(H, W),
    cropped_layout_dict(Numbered, Words, GridLen, MinR, MinC, S, Dict).

% The S x S crop dict of a GridLen grid, origin at (MinR,MinC); same JSON shape
% as emit_json/3, reusing add_word_cells/3, cell_coord/3 and answer_meta_assoc/2.
cropped_layout_dict(PlacedWords, Words, GridLen, MinR, MinC, S,
                    _{gridLength: S, grid: Rows, words: WordObjs}) :-
    empty_assoc(A0),
    foldl(add_word_cells, PlacedWords, A0, CellMap),
    Smax is S - 1, numlist(0, Smax, Idxs),
    maplist(cropped_row(CellMap, GridLen, MinR, MinC, Idxs), Idxs, Rows),
    answer_meta_assoc(Words, MetaAssoc),
    maplist(cropped_word(MetaAssoc, GridLen, MinR, MinC), PlacedWords, WordObjs).

cropped_row(CellMap, GridLen, MinR, MinC, CIdxs, R, Row) :-
    maplist(cropped_cell(CellMap, GridLen, MinR, MinC, R), CIdxs, Row).

cropped_cell(CellMap, GridLen, MinR, MinC, R, C, Json) :-
    OrigR is MinR + R, OrigC is MinC + C,
    (   OrigR < GridLen, OrigC < GridLen,
        OrigCell is OrigR * GridLen + OrigC + 1,
        get_assoc(OrigCell, CellMap, cell(Letter, Ac, Dn, N))
    ->  Json = _{letter: Letter, number: N, across: Ac, down: Dn}
    ;   Json = null
    ).

cropped_word(MetaAssoc, GridLen, MinR, MinC, PW, WordObj) :-
    get_dict(answer, PW, Answer),
    get_dict(dir, PW, Dir),
    get_dict(num, PW, Num),
    get_dict(cells, PW, Cells),
    maplist(cropped_coord(GridLen, MinR, MinC), Cells, Coords),
    get_assoc(Answer, MetaAssoc, Meta),
    WordObj = _{number: Num, direction: Dir, answer: Answer,
                cells: Coords, meta: Meta}.

cropped_coord(GridLen, MinR, MinC, Cell, [R, C]) :-
    cell_coord(GridLen, Cell, [R0, C0]),
    R is R0 - MinR, C is C0 - MinC.

% Drop-contract dispatch (strict | best_effort) -> the matching solver. The
% crosswordsmith CLI is the entry point; it calls arrange_solve/4 (and the
% fragment/candidates/enumerate solvers) directly on a loaded Words list.
arrange_solve(Words, GridLen, strict, SizeMode) :-
    arrange_strict_solve(Words, GridLen, SizeMode).
arrange_solve(Words, GridLen, best_effort, SizeMode) :-
    arrange_best_effort_solve(Words, GridLen, SizeMode).


% ===========================================================================
% Phase 5 - fragment-grid seeding (the anchor mechanism). design-spec §6.6.
%
% A fragment grid is the emit JSON made partial (words-only v1): presence =
% fixed. The engine pins every fragment word at exactly its cells, then solves
% the unmatched remainder with the SAME construct path as the unseeded engine.
% This generalises the single seed: instead of seeding one word at a start
% corner, we seed the whole fragment into (Placed, Grid) and shrink Words.
%
% Load-bearing properties:
%   AC-EMIT-2 emit -> re-ingest as fragment -> identical layout (an emitted
%             layout IS already a valid fragment).
%   AC-FRAG-3 pinned words sit at exactly their fragment positions.
% Everything here reuses crossword.pl's legality core (assign_word and the
% adjacency/merge/prev-cell checks it calls), geometry (cell_coord, word_cells,
% fits_on_grid) and the MRV-inc remainder search; nothing is re-derived.
% ===========================================================================

:- multifile prolog:error_message//1.

% --- fragment schema parsing (the emit format, made partial) ----------------
% A fragment dict -> GridLen + [frag(Answer, Dir, Start, CellNums)]. Validates
% shape only; reconciliation against --input and legality are seed_from_fragment's
% job. Throws shaped error/2 terms (rendered by the hooks at the foot of file).

fragment_dict_words(Dict, GridLen, Frags) :-
    (   is_dict(Dict), get_dict(gridLength, Dict, GridLen),
        integer(GridLen), GridLen > 0
    ->  true
    ;   throw(error(fragment_no_grid_length, _))
    ),
    (   get_dict(words, Dict, WordEntries), is_list(WordEntries)
    ->  true
    ;   throw(error(fragment_no_words_array, _))
    ),
    maplist(fragment_word(GridLen), WordEntries, Frags).

% One fragment word entry -> frag(Answer, Dir, Start, CellNums). Only
% answer/direction/cells are read; emit's `number`/`meta` are ignored (number
% is reassigned by clue numbering, meta rejoined from --input at emit time).
% Start is the 1-based cell number of cells[0]; CellNums is the full declared
% run as cell numbers (validated against the answer in seed_from_fragment).
fragment_word(GridLen, Entry, frag(Answer, Dir, Start, CellNums)) :-
    (   is_dict(Entry), get_dict(answer, Entry, RawAns), fragment_atom(RawAns, Answer)
    ->  true
    ;   throw(error(fragment_invalid_answer(Entry), _))
    ),
    (   get_dict(direction, Entry, RawDir), fragment_dir(RawDir, Dir)
    ->  true
    ;   throw(error(fragment_invalid_direction(Answer), _))
    ),
    (   get_dict(cells, Entry, RawCells), is_list(RawCells), RawCells = [_|_]
    ->  true
    ;   throw(error(fragment_no_cells(Answer), _))
    ),
    maplist(cell_pair_to_num(GridLen, Answer), RawCells, CellNums),
    CellNums = [Start|_].

% Normalise a JSON string-or-atom to an atom (independent of the double_quotes
% flag, which differs between consulted source and json_read_dict output).
fragment_atom(X, A) :- atom(X), !, A = X.
fragment_atom(X, A) :- string(X), !, atom_string(A, X).

fragment_dir(Raw, Dir) :- fragment_atom(Raw, A), fragment_dir_atom(A, Dir).
fragment_dir_atom(across, across).
fragment_dir_atom(down,   down).

% A [Row, Col] pair (0-based, on-grid) -> 1-based row-major cell number.
cell_pair_to_num(GridLen, Answer, Pair, Num) :-
    (   Pair = [R, C], integer(R), integer(C),
        R >= 0, R < GridLen, C >= 0, C < GridLen
    ->  Num is R * GridLen + C + 1
    ;   throw(error(fragment_invalid_cell(Answer, Pair), _))
    ).

% Read a fragment file (JSON - the emit format) into GridLen + frags.
load_fragment(File, GridLen, Frags) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Dict), close(S)),
    fragment_dict_words(Dict, GridLen, Frags).

% The fragment's gridLength sets N; an explicit --size is redundant and an
% error if it disagrees (design-spec §6.6). Ready for the Phase-7 CLI.
reconcile_fragment_size(FragGridLen, none, FragGridLen) :- !.
reconcile_fragment_size(FragGridLen, FragGridLen, FragGridLen) :- !.
reconcile_fragment_size(FragGridLen, OptSize, _) :-
    throw(error(fragment_size_mismatch(FragGridLen, OptSize), _)).


% --- seeding: pin the fragment, validate up front, shrink Words -------------
% Pin every fragment word into (Placed, Grid) via the legality core, returning
% the seeded state and the remaining (unpinned) input words. ALL validation
% happens here, before any search: AC-FRAG-1 (answer in --input), the declared
% cells being a legal straight run, and AC-FRAG-2 (the pin is legal against the
% already-pinned words). Any violation throws a shaped error.
seed_from_fragment(Frags, Words, GridLen, SeededPlaced, SeededGrid, Remaining) :-
    check_fragment_unique(Frags),
    sort_frags_by_start(Frags, Ordered),
    init_grid(GridLen, G0),
    foldl(pin_fragment_word(Words, GridLen), Ordered, []-G0, SeededPlaced-SeededGrid),
    fragment_answers(Ordered, Pinned),
    remaining_words(Words, Pinned, Remaining).

% Pin one fragment word. Order of checks is deliberate: input-membership first
% (the clearest authoring error), then declared-cells consistency, then the
% on-grid legality pin (whose failure is diagnosed by fragment_conflict/7).
pin_fragment_word(Words, GridLen, frag(Answer, Dir, Start, FragCellNums),
                  PlacedIn-GridIn, [PW|PlacedIn]-GridOut) :-
    (   reconcile_answer(Answer, Words, Entry)
    ->  true
    ;   throw(error(fragment_word_not_in_input(Answer), _))
    ),
    word_letters(Entry, Letters, WLen),
    (   expected_run(Start, Dir, WLen, GridLen, FragCellNums)
    ->  true
    ;   throw(error(fragment_cells_inconsistent(Answer), _))
    ),
    (   assign_word(Answer, Letters, WLen, Start, Dir, GridLen,
                    PlacedIn, GridIn, PW, GridOut)
    ->  true
    ;   fragment_conflict(Answer, Letters, Start, Dir, WLen, GridLen, GridIn)
    ).

% The declared cells must be exactly the straight, on-grid run of the answer's
% length from Start in Dir - i.e. the fragment's geometry agrees with its own
% answer/direction. (For a fragment produced by emit this always holds.)
expected_run(Start, Dir, WLen, GridLen, FragCellNums) :-
    fits_on_grid(Dir, Start, WLen, GridLen),
    word_cells(Start, Dir, WLen, GridLen, Run),
    Run == FragCellNums.

% Diagnose (and throw on) an illegal pin. A clashing fixed letter is the common,
% precisely-reportable case; anything else (adjacency/merge/abutment) is reported
% generically. fragment_conflict/7 always throws.
fragment_conflict(Answer, Letters, Start, Dir, WLen, GridLen, Grid) :-
    word_cells(Start, Dir, WLen, GridLen, Cells),
    (   clashing_cell(Cells, Letters, Grid, Cell, Existing, Wanted)
    ->  cell_coord(GridLen, Cell, RC),
        throw(error(fragment_letter_clash(Answer, RC, Existing, Wanted), _))
    ;   throw(error(fragment_illegal_pin(Answer), _))
    ).

% First cell whose grid letter is already fixed to something the word cannot
% supply (an overlap conflict). Cells and Letters are the parallel run.
clashing_cell([Cell|_], [L|_], Grid, Cell, Existing, L) :-
    get_assoc(Cell, Grid, Existing),
    Existing \== empty,
    Existing \== L,
    !.
clashing_cell([_|Cs], [_|Ls], Grid, Cell, Existing, Wanted) :-
    clashing_cell(Cs, Ls, Grid, Cell, Existing, Wanted).

% Reconcile a fragment answer to its --input entry (==-keyed on the answer).
reconcile_answer(Answer, Words, Entry) :-
    member(Entry, Words), Entry = [A|_], A == Answer, !.

fragment_answers(Frags, Answers) :-
    findall(A, member(frag(A, _, _, _), Frags), Answers).

remaining_words(Words, Pinned, Remaining) :-
    findall(E,
            ( member(E, Words), E = [A|_], \+ memberchk(A, Pinned) ),
            Remaining).

% Stable order by start cell so pin order, error reporting and numbering are
% deterministic (an across+down sharing a start keep fragment order).
sort_frags_by_start(Frags, Ordered) :-
    map_list_to_pairs(frag_start, Frags, Pairs),
    keysort(Pairs, Sorted),
    pairs_values(Sorted, Ordered).
frag_start(frag(_, _, Start, _), Start).

% A fragment must not pin the same answer twice (silent double-place guard;
% mirrors check_unique_answers/1 for the input set).
check_fragment_unique(Frags) :-
    fragment_answers(Frags, As),
    msort(As, Sorted),
    (   append(_, [D, D|_], Sorted)
    ->  throw(error(fragment_duplicate_answer(D), _))
    ;   true
    ).


% --- solving from the seed --------------------------------------------------
% Generalises the single-seed construct: start from the fragment-seeded
% (Placed, Grid) and place Remaining via the reused MRV-inc path (select_inc's
% non-empty-Placed clause fires immediately). One deterministic construction -
% the fragment fixes the canvas and anchors, so there is no corner/seed sweep.
construct_from_seed(Remaining, SeededPlaced, SeededGrid, GridLen, AllPlaced) :-
    assign_words_inc(Remaining, SeededPlaced, none, GridLen, _S, _D,
                     SeededGrid, _GOut, AllPlaced).

% Strict fragment solve: pin (validating up front), construct the remainder
% under a budget, rescore, number. Outcome: placed | not_proven | infeasible.
% seed_from_fragment runs OUTSIDE the budget so conflicts are reported before
% any search (the §6.6 "validate up front" contract).
arrange_fragment_strict(Words, Frags, GridLen, Numbered, Reward, Outcome) :-
    arrange_weights(WCap, WTail),
    seed_from_fragment(Frags, Words, GridLen, SeededPlaced, SeededGrid, Remaining),
    arrange_budget(Budget),
    construct_fragment_one(Remaining, SeededPlaced, SeededGrid, GridLen, Budget, Res),
    (   Res = ok(AllPlaced)
    ->  layout_reward(WCap, WTail, AllPlaced, Reward),
        once(assign_clue_numbers(AllPlaced, Numbered)),
        Outcome = placed
    ;   Numbered = [], Reward = -1,
        ( Res == budget -> Outcome = not_proven ; Outcome = infeasible )
    ).

% First MRV-inc completion from the seed, under a budget. Always succeeds with a
% tag: ok(Placed) | budget | exhausted (mirrors construct_one/7).
construct_fragment_one(Remaining, SeededPlaced, SeededGrid, GridLen, Budget, Res) :-
    % No catch/3: see construct_one/7 - the budget is handled inside
    % call_with_inference_limit/3, and a thrown error is a real bug to surface.
    call_with_inference_limit(
        ( once(construct_from_seed(Remaining, SeededPlaced, SeededGrid, GridLen, P))
        ->  Found = P
        ;   Found = none ),
        Budget, Outcome),
    (   Outcome == inference_limit_exceeded -> Res = budget
    ;   Found == none                       -> Res = exhausted
    ;   Res = ok(Found)
    ).

% Best-effort fragment solve: pin, then place the remainder with the greedy
% constructor (which drops what it cannot place). greedy_loop/6 already accepts a
% starting (Placed, Grid), and only ever drops from Remaining - the pins are
% never dropped (AC-FRAG-3 holds under drop too).
arrange_fragment_best_effort(Words, Frags, GridLen, Numbered, Reward, NumPlaced, Dropped) :-
    arrange_weights(WCap, WTail),
    seed_from_fragment(Frags, Words, GridLen, SeededPlaced, SeededGrid, Remaining),
    greedy_loop(Remaining, SeededPlaced, GridLen, SeededGrid, FinalPlaced, DroppedEntries),
    layout_reward(WCap, WTail, FinalPlaced, Reward),
    length(FinalPlaced, NumPlaced),
    findall(A, member([A|_], DroppedEntries), Dropped),
    once(assign_clue_numbers(FinalPlaced, Numbered)).


% --- fragment solve entry points (emit + report) ---------------------------
% Mirror arrange_strict_solve/arrange_best_effort_solve: emit on stdout, report
% on stderr; on any non-placed outcome report and fail with no stdout layout.

arrange_fragment_solve(Words, Frags, GridLen, strict, SizeMode) :-
    check_unique_answers(Words),
    arrange_fragment_strict(Words, Frags, GridLen, Numbered, Reward, Outcome),
    (   Outcome == placed
    ->  emit_arrange(Numbered, Words, GridLen, SizeMode),
        length(Numbered, NP),
        format(user_error,
               "arrange: fragment-seeded, grid ~w, mode ~w, placed ~w, reward ~w~n",
               [GridLen, SizeMode, NP, Reward])
    ;   arrange_fragment_report_failure(Outcome, Words, GridLen),
        fail
    ).
arrange_fragment_solve(Words, Frags, GridLen, best_effort, SizeMode) :-
    check_unique_answers(Words),
    (   arrange_fragment_best_effort(Words, Frags, GridLen, Numbered, Reward, NP, Dropped)
    ->  emit_arrange(Numbered, Words, GridLen, SizeMode),
        length(Dropped, ND),
        format(user_error,
               "arrange: fragment-seeded, grid ~w, mode ~w, placed ~w, dropped ~w ~w, reward ~w~n",
               [GridLen, SizeMode, NP, ND, Dropped, Reward])
    ;   format(user_error,
               "arrange: nothing placeable around the fragment on ~wx~w grid~n",
               [GridLen, GridLen]),
        fail
    ).

arrange_fragment_report_failure(not_proven, _Words, GridLen) :-
    format(user_error,
           "arrange: fragment search not proven within budget on ~wx~w grid \c
(search did not complete; feasibility unresolved)~n",
           [GridLen, GridLen]).
arrange_fragment_report_failure(infeasible, Words, GridLen) :-
    (   unplaceable_words(Words, Bad), Bad = [_|_]
    ->  format(user_error,
               "arrange: cannot place all words around the fragment on ~wx~w grid; \c
words with no possible crossing: ~w~n",
               [GridLen, GridLen, Bad])
    ;   format(user_error,
               "arrange: no complete placement around the fragment on ~wx~w grid \c
(remaining words cannot all interlock with the pinned seed)~n",
               [GridLen, GridLen])
    ).


% --- fragment error messages ------------------------------------------------
prolog:error_message(fragment_no_grid_length) -->
    [ 'fragment: expected a JSON object with a positive integer "gridLength"' ].
prolog:error_message(fragment_no_words_array) -->
    [ 'fragment: expected a "words" array' ].
prolog:error_message(fragment_invalid_answer(Entry)) -->
    [ 'fragment: every word needs a string "answer" (offending entry: ~q)'-[Entry] ].
prolog:error_message(fragment_invalid_direction(Answer)) -->
    [ 'fragment: word ~q needs a "direction" of "across" or "down"'-[Answer] ].
prolog:error_message(fragment_no_cells(Answer)) -->
    [ 'fragment: word ~q needs a non-empty "cells" array'-[Answer] ].
prolog:error_message(fragment_invalid_cell(Answer, Pair)) -->
    [ 'fragment: word ~q has an invalid cell ~q (need [row,col] within the grid)'-[Answer, Pair] ].
prolog:error_message(fragment_duplicate_answer(Answer)) -->
    [ 'fragment: answer ~q appears more than once; a fragment pins each answer at most once'-[Answer] ].
prolog:error_message(fragment_word_not_in_input(Answer)) -->
    [ 'fragment: word ~q is not one of the --input answers (a fragment word must come from the input set)'-[Answer] ].
prolog:error_message(fragment_cells_inconsistent(Answer)) -->
    [ 'fragment: the cells of word ~q are not a straight on-grid run matching its direction and length'-[Answer] ].
prolog:error_message(fragment_letter_clash(Answer, RC, Existing, Wanted)) -->
    [ 'fragment: word ~q conflicts at cell ~w - already-pinned letter ~w but this word needs ~w'-[Answer, RC, Existing, Wanted] ].
prolog:error_message(fragment_illegal_pin(Answer)) -->
    [ 'fragment: word ~q cannot be legally pinned (it abuts or merges with another pinned word)'-[Answer] ].
prolog:error_message(fragment_size_mismatch(FragGridLen, OptSize)) -->
    [ 'fragment: gridLength ~w disagrees with the requested --size ~w'-[FragGridLen, OptSize] ].


% ===========================================================================
% Phase 6 - candidates (diverse alternative layouts). design-spec §7.4, AC-ARR-7.
%
% Default output is a single deterministic best. `--candidates K` opts into up
% to K *meaningfully distinct* layouts. Distinctness comes from CONSTRUCTOR
% BREADTH (the greedy constructor over seed x start-corner) + GREEDY DIVERSITY:
% rank the pool best-first, take the best, then each next-best whose
% placement-distance >= tau from ALL already-picked. This is why candidates ride
% the greedy path, not a single deterministic search's near-duplicate leaves.
%
% Placement distance is TRANSLATION-INVARIANT: each word's position is taken
% relative to its layout's bounding-box origin, so two layouts that differ only
% by a global shift count as the same candidate (the same crossword), matching
% AC-ARR-7's "meaningfully distinct". Distance = fraction of answers whose
% (rel-row, rel-col, direction) differs; integer-compared (no floats) for
% deterministic output (INV-2).
% ===========================================================================

% tau = 0.30 placement-distance threshold, as integer percent (OD-9: tunable /
% to be calibrated against the fixtures alongside epsilon/target).
candidate_tau_pct(30).

% The candidate pool: every greedy construction over seed x start-corner,
% ranked best-first by score(NumPlaced, Reward) in standard term order. Under
% --strict only full placements are eligible. Each entry is tagged with its
% translation-invariant placement assoc (answer -> RelRow-RelCol-Dir).
arrange_candidate_pool(Words, GridLen, DropContract, Pool) :-
    arrange_weights(WCap, WTail),
    length(Words, Total),
    seed_candidates(Words, Seeds),
    start_locs(Locs),
    findall(score(NP, R)-Placed,
            ( member(Loc, Locs),
              member(Seed, Seeds),
              greedy_construct(Words, GridLen, Loc, Seed, Placed, _Dropped),
              length(Placed, NP),
              ( DropContract == strict -> NP =:= Total ; true ),
              layout_reward(WCap, WTail, Placed, R) ),
            Raw),
    sort(1, @>=, Raw, Sorted),                 % best-first; stable for ties
    pairs_values(Sorted, Placeds),
    maplist(tag_with_assoc(GridLen), Placeds, Pool).

tag_with_assoc(GridLen, Placed, c(Placed, Assoc)) :-
    placement_assoc(Placed, GridLen, Assoc).

% answer -> (RelRow - RelCol - Dir), positions relative to the layout's bbox
% origin (translation-invariant). Answers are unique within a layout.
placement_assoc(Placed, GridLen, Assoc) :-
    placed_bbox(Placed, GridLen, bbox(MinR, _MaxR, MinC, _MaxC), _Area),
    findall(A-(RR-RC-D),
            ( member(PW, Placed),
              get_dict(answer, PW, A), get_dict(start, PW, S), get_dict(dir, PW, D),
              cell_coord(GridLen, S, [SR, SC]),
              RR is SR - MinR, RC is SC - MinC ),
            Pairs),
    list_to_assoc(Pairs, Assoc).

% Greedy diversity selection: best-first, keep each layout that is >= tau from
% every already-kept one, until K kept or the pool is exhausted.
pick_diverse(Pool, TauPct, Total, K, Picked) :-
    pick_diverse_(Pool, TauPct, Total, K, [], Rev),
    reverse(Rev, Picked).

pick_diverse_([], _TauPct, _Total, _K, Acc, Acc).
pick_diverse_([c(P, A)|Cs], TauPct, Total, K, Acc, Out) :-
    length(Acc, L),
    (   L >= K
    ->  Out = Acc
    ;   Acc == []
    ->  pick_diverse_(Cs, TauPct, Total, K, [c(P, A)], Out)        % best: always kept
    ;   far_from_all(A, Acc, TauPct, Total)
    ->  pick_diverse_(Cs, TauPct, Total, K, [c(P, A)|Acc], Out)
    ;   pick_diverse_(Cs, TauPct, Total, K, Acc, Out)
    ).

far_from_all(Assoc, Acc, TauPct, Total) :-
    forall(member(c(_, A2), Acc),
           ( pos_diff_count(Assoc, A2, Diff), Diff * 100 >= TauPct * Total )).

% Number of answers placed differently between two layouts (present-vs-absent or
% a differing relative position/direction both count as a difference).
pos_diff_count(A1, A2, Diff) :-
    assoc_to_keys(A1, K1), assoc_to_keys(A2, K2),
    ord_union(K1, K2, Keys),
    aggregate_all(count, ( member(K, Keys), \+ same_pos(K, A1, A2) ), Diff).

same_pos(K, A1, A2) :- get_assoc(K, A1, P), get_assoc(K, A2, P).

% Up to K diverse numbered layouts for the word set. Returned =< K; fewer only
% when fewer >= tau-distinct layouts exist (reported by the solve wrapper).
arrange_candidates(Words, GridLen, DropContract, K, Numbered, Returned) :-
    arrange_candidate_pool(Words, GridLen, DropContract, Pool),
    Pool = [_|_],
    length(Words, Total),
    candidate_tau_pct(TauPct),
    pick_diverse(Pool, TauPct, Total, K, Picked),
    findall(N, ( member(c(P, _), Picked), once(assign_clue_numbers(P, N)) ), Numbered),
    length(Numbered, Returned).

% Emit the candidates as a JSON ARRAY of canonical layout dicts (each element is
% itself a valid fragment for re-ingestion). Single-layout emit stays a bare
% object; the array shape is what `--candidates` opts into.
emit_candidates(NumberedLayouts, Words, GridLen, SizeMode) :-
    maplist(candidate_dict(Words, GridLen, SizeMode), NumberedLayouts, Dicts),
    current_output(Out),
    json_write_dict(Out, Dicts),
    nl(Out).
candidate_dict(Words, GridLen, SizeMode, Numbered, Dict) :-
    arrange_layout_dict(Numbered, Words, GridLen, SizeMode, Dict).

% Solve + emit the candidates array on stdout, report requested/returned/tau on
% stderr (INV-3: a short return is reported, never silent). Fails (no stdout) if
% nothing is placeable.
arrange_candidates_solve(Words, GridLen, DropContract, SizeMode, K) :-
    check_unique_answers(Words),
    (   arrange_candidates(Words, GridLen, DropContract, K, Layouts, Returned),
        Returned > 0
    ->  emit_candidates(Layouts, Words, GridLen, SizeMode),
        candidate_tau_pct(TauPct),
        (   Returned < K
        ->  format(user_error,
                   "arrange: ~w candidate(s) requested, ~w returned (tau ~w%); \c
fewer >=tau-distinct layouts exist~n",
                   [K, Returned, TauPct])
        ;   format(user_error,
                   "arrange: ~w candidate(s) requested, ~w returned (tau ~w%)~n",
                   [K, Returned, TauPct])
        )
    ;   format(user_error, "arrange: no candidate layout on ~wx~w grid~n",
               [GridLen, GridLen]),
        fail
    ).


% --- enumerate (count all feasible full placements) ------------------------
% The `--enumerate` engine seam (design-spec §7.1, AC-ARR-8): the exhaustive
% count of every feasible full placement on an N x N grid, over all four start
% corners. Reuses crossword.pl's all_crossword/5 with the production default
% strategy, so the count matches the old `--all` exactly (strategies only
% reorder the same search tree). The branch step / legality core are untouched.
arrange_enumerate(Words, GridLen, Num) :-
    default_strategy(Strat),
    all_crossword(Strat, GridLen, Words, _StartLoc, Num).

% Count + print the enumeration on stdout (a bare integer), as the old CLI did.
arrange_enumerate_solve(Words, GridLen) :-
    check_unique_answers(Words),
    arrange_enumerate(Words, GridLen, Num),
    current_output(Out),
    writeln(Out, Num).
