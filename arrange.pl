% arrange.pl - Flavour-A layout engine: scoring oracle (Phase 1) + the
% search-value gate (Phase 1.5). See docs/arrange-implementation-plan.md.
%
% This file deliberately does NOT define `main`/initialization, so it can be
% consulted by a harness (or `swipl -g gate_run`) without side effects. It
% consults crossword.pl (which in turn loads quality.pl), reusing the legality
% core, the MRV-inc branch step, and the existing checking metrics verbatim.

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

% ceil(L/2): the per-word checking cap (== word_meets_half/2's threshold in
% quality.pl). L is the word length in cells (the placed dict's `len`).
check_target(L, T) :- T is (L + 1) // 2.

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


% ---------------------------------------------------------------------------
% Phase 1.5 - the search-value gate.
%
% A minimal branch-and-bound over crossword.pl's MRV-inc branch step
% (select_inc -> find_intersecting_word -> assign_word, all reused verbatim):
% record the FIRST complete placement (the MRV incumbent), then continue the
% DFS recording any STRICTLY BETTER complete placement, pruned by a cheap,
% admissible upper bound. Node and prune counters + an inference budget make
% it a measurement: does any search beat the first incumbent, and does the
% bound ever fire?
% ---------------------------------------------------------------------------

% Inference budget per fixture (call_with_inference_limit/3). Generous: the
% point is to give the search its best shot at finding a better layout.
gate_budget(800_000_000).

incr(Key) :- nb_getval(Key, V), V1 is V + 1, nb_setval(Key, V1).

record_first(R) :-
    nb_getval(arrange_first, F),
    ( F == none -> nb_setval(arrange_first, R) ; true ).

record_best(Placed, R) :-
    nb_getval(arrange_best, _-BestR),
    ( R > BestR -> nb_setval(arrange_best, Placed-R) ; true ).

% Cheap, loose, ADMISSIBLE upper bound on the reward of any completion from
% this node: current partial reward + optimistic room on each placed word
% (fill all its remaining cells as checked, up to the cap) + optimistic max for
% each unplaced word (fully capped and fully checked). Never under-estimates,
% so pruning a branch with UB =< incumbent cannot discard a better layout.
upper_bound(Unplaced, Placed, WCap, WTail, UB) :-
    layout_reward(WCap, WTail, Placed, PartialR),
    foldl(placed_room(WCap, WTail, Placed), Placed, 0, PRoom),
    foldl(unplaced_max(WCap, WTail), Unplaced, 0, UMax),
    UB is PartialR + PRoom + UMax.

placed_room(WCap, WTail, Placed, PW, A, A1) :-
    get_dict(len, PW, L),
    check_target(L, T),
    word_checked_count(PW, Placed, C),
    Room is WCap * (T - min(C, T)) + WTail * (L - C),
    A1 is A + Room.

unplaced_max(WCap, WTail, Entry, A, A1) :-
    word_letters(Entry, _, L),
    check_target(L, T),
    A1 is A + WCap * T + WTail * L.

% Complete placement: record incumbents, then FAIL to keep searching.
gate_search([], Placed, _GridLen, WCap, WTail, _G, _St, _Start, _Dir) :-
    incr(arrange_nodes),
    layout_reward(WCap, WTail, Placed, R),
    record_first(R),
    record_best(Placed, R),
    fail.
% Partial: count the node, prune if the bound says this subtree cannot beat the
% incumbent, else branch via the reused MRV-inc step and recurse.
gate_search([W|Ws], Placed, GridLen, WCap, WTail, G, St, Start, Dir) :-
    incr(arrange_nodes),
    nb_getval(arrange_best, _-BestR),
    (   BestR >= 0,
        upper_bound([W|Ws], Placed, WCap, WTail, UB),
        UB =< BestR
    ->  incr(arrange_prunes),
        fail
    ;   select_inc([W|Ws], Placed, St, GridLen, Start, Dir, G, Entry, RemWords, St1),
        Entry = [Word|_],
        word_letters(Entry, Letters, WLen),
        find_intersecting_word(Letters, WLen, Placed, GridLen, Start, Dir),
        assign_word(Word, Letters, WLen, Start, Dir, GridLen, Placed, G, PW, G1),
        gate_search(RemWords, [PW|Placed], GridLen, WCap, WTail, G1, St1, _S, _D)
    ).

% Run the exhaustive-with-budget search for one word set on an N x N grid,
% seeded at the canonical topleft_across corner (mirrors find_crossword/6's
% mrv_inc path). Records first/best/nodes/prunes/outcome in globals.
gate_search_run(Words, GridLen, Outcome) :-
    arrange_weights(WCap, WTail),
    nb_setval(arrange_first, none),
    nb_setval(arrange_best, none-(-1)),
    nb_setval(arrange_nodes, 0),
    nb_setval(arrange_prunes, 0),
    init_grid(GridLen, G0),
    start_loc(topleft_across, GridLen, StartNum, StartDir),
    gate_budget(Budget),
    catch(
        (   call_with_inference_limit(
                ( gate_search(Words, [], GridLen, WCap, WTail, G0, none,
                              StartNum, StartDir), fail ; true ),
                Budget, Outcome0)
        ->  Outcome = Outcome0
        ;   Outcome = failed
        ),
        E, Outcome = error(E)).

% A greedy-construct layout's reward, for comparison (it may drop words).
greedy_reward(Words, GridLen, R, NP, ND) :-
    arrange_weights(WCap, WTail),
    (   seed_candidates(Words, [Seed|_]),
        greedy_construct(Words, GridLen, topleft_across, Seed, Placed, Dropped)
    ->  layout_reward(WCap, WTail, Placed, R),
        length(Placed, NP), length(Dropped, ND)
    ;   R = -1, NP = 0, ND = -1
    ).


% ---------------------------------------------------------------------------
% Gate driver + report.
% ---------------------------------------------------------------------------

% (fixture file, canonical grid size) - the puzzle-shaped inputs where the DFS
% seeds quickly (probe: 5-462 ms). Grids per benchmarks/fixtures + the probe.
gate_fixtures([
    fix('fixtures/benchmark_08_words.pl', 13),
    fix('fixtures/benchmark_14_words.pl', 17),
    fix('fixtures/benchmark_20_words.pl', 37),
    fix('fixtures/benchmark_26_words.pl', 49),
    fix('fixtures/toc_demo.pl',           25)
]).

gate_one(File, GridLen) :-
    load_clues(File, Words),
    length(Words, NW),
    gate_search_run(Words, GridLen, Outcome),
    nb_getval(arrange_first, FirstR),
    nb_getval(arrange_best, BestPlaced-BestR),
    nb_getval(arrange_nodes, Nodes),
    nb_getval(arrange_prunes, Prunes),
    ( BestPlaced == none -> CapBind = 'n/a (no full placement)'
    ; cap_binding_count(BestPlaced, CapBind) ),
    greedy_reward(Words, GridLen, GR, GNP, GND),
    format("~n=== ~w  (~w words, grid ~w) ===~n", [File, NW, GridLen]),
    format("  search : first=~w  best=~w  delta=~w  nodes=~w  prunes=~w  outcome=~w~n",
           [FirstR, BestR, '?', Nodes, Prunes, Outcome]),
    report_delta(FirstR, BestR),
    format("  capbind: ~w of ~w placed words meet ceil(L/2) in best layout~n",
           [CapBind, NW]),
    format("  greedy : reward=~w  placed=~w  dropped=~w~n", [GR, GNP, GND]).

report_delta(none, _) :- !,
    format("  verdict: NO complete placement found within budget~n").
report_delta(F, B) :-
    D is B - F,
    ( D =:= 0
    ->  format("  verdict: best == first (search added NOTHING over the first incumbent)~n")
    ;   format("  verdict: best > first by ~w (search improved the layout)~n", [D])
    ).

gate_run :-
    gate_fixtures(Fixtures),
    format("arrange Phase 1.5 search-value gate  (weights 5:1, budget per gate_budget/1)~n"),
    forall(member(fix(File, Grid), Fixtures), gate_one(File, Grid)),
    nl.


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
    arrange_weights(WCap, WTail),
    arrange_budget(Budget),
    start_locs(Locs),
    findall(Res,
            ( member(Loc, Locs),
              construct_one(Loc, Words, GridLen, WCap, WTail, Budget, Res) ),
            Results),
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

% First MRV-inc complete placement from one corner, rescored, under a budget.
% Always succeeds with a Res tag: ok(Reward,Placed) | budget | exhausted.
construct_one(Loc, Words, GridLen, WCap, WTail, Budget, Res) :-
    catch(
        call_with_inference_limit(
            ( once(find_crossword(mrv_inc, GridLen, Words, Loc, _G, P)) -> Found = P
            ; Found = none ),
            Budget, Outcome),
        _Err, (Outcome = error, Found = none)),
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
        format(user_error, "arrange: grid ~w, mode ~w, placed ~w, reward ~w~n",
               [GridLen, SizeMode, NP, Reward])
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
    intersection(Ls, OLs, [_|_]).


% --- emit framing (Phase 3) ------------------------------------------------
% fixed : exactly GridLen x GridLen (the canonical emit_json/3).
% max   : crop to the tight SQUARE enclosing the content (side = max(H,W)),
%         anchored at the bbox top-left. Square keeps the single-`gridLength`
%         output schema (rectangular grids are out of scope, design-spec §3).

emit_arrange(Numbered, Words, GridLen, fixed) :-
    emit_json(Numbered, Words, GridLen).
emit_arrange(Numbered, Words, GridLen, max) :-
    placed_bbox(Numbered, GridLen, bbox(MinR, MaxR, MinC, MaxC), _Area),
    H is MaxR - MinR + 1, W is MaxC - MinC + 1, S is max(H, W),
    emit_cropped(Numbered, Words, GridLen, MinR, MinC, S).

% Emit an S x S crop of a GridLen grid, origin at (MinR,MinC); same JSON shape
% as emit_json/3, reusing add_word_cells/3, cell_coord/3 and answer_meta/3.
emit_cropped(PlacedWords, Words, GridLen, MinR, MinC, S) :-
    empty_assoc(A0),
    foldl(add_word_cells, PlacedWords, A0, CellMap),
    Smax is S - 1, numlist(0, Smax, Idxs),
    maplist(cropped_row(CellMap, GridLen, MinR, MinC, Idxs), Idxs, Rows),
    maplist(cropped_word(Words, GridLen, MinR, MinC), PlacedWords, WordObjs),
    Payload = _{gridLength: S, grid: Rows, words: WordObjs},
    current_output(Out),
    json_write_dict(Out, Payload),
    nl(Out).

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

cropped_word(Words, GridLen, MinR, MinC, PW, WordObj) :-
    get_dict(answer, PW, Answer),
    get_dict(dir, PW, Dir),
    get_dict(num, PW, Num),
    get_dict(cells, PW, Cells),
    maplist(cropped_coord(GridLen, MinR, MinC), Cells, Coords),
    answer_meta(Answer, Words, Meta),
    WordObj = _{number: Num, direction: Dir, answer: Answer,
                cells: Coords, meta: Meta}.

cropped_coord(GridLen, MinR, MinC, Cell, [R, C]) :-
    cell_coord(GridLen, Cell, [R0, C0]),
    R is R0 - MinR, C is C0 - MinC.

% Convenience runner: load a clue file and emit a strict arrange layout.
arrange_run(File, GridLen, SizeMode) :-
    load_clues(File, Words),
    arrange_strict_solve(Words, GridLen, SizeMode).
