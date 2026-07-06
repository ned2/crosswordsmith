% arrange.pl - Flavour-A layout engine: the scoring oracle + the
% construct/rescore/emit driver, fragment seeding, diverse candidates, and the
% greedy density constructor (moved here from metrics.pl - see the section at
% the foot of this file). See docs/arrange-implementation-plan.md and
% docs/design-spec.md §7.
%
% This file deliberately does NOT define `main`/initialization, so it can be
% loaded by a harness without side effects. It reuses core.pl's legality core
% and MRV-inc branch step and the checking metrics — both explicit imports
% below (crosswordsmith_core, crosswordsmith_metrics). (The Phase-1.5
% `gate_*` measurement harness that lived here was removed once the gate's
% descope decision was recorded in the implementation plan.)
%
% Exports: the four *_solve CLI seams, the fragment API the driver and fill
% share (load_fragment/3, reconcile_fragment_size/3), and set_check_target/1
% (the only sanctioned writer of the module-private check_target_override/1
% dynamic). White-box test helpers are NOT exported — arrange.plt reaches
% them as crosswordsmith_arrange:Pred(...).

:- module(crosswordsmith_arrange,
          [ arrange_solve/4,
            arrange_fragment_solve/5,
            arrange_candidates_solve/5,
            arrange_enumerate_solve/2,
            % the browser/WASM envelope seam (no emit, tagged outcome) —
            % crosswordsmith_browser is its consumer
            arrange_outcome/5,
            load_fragment/3,
            reconcile_fragment_size/3,
            set_check_target/1,
            % fill's emit_fill(max) delegates its cropped emit here
            emit_arrange/4
          ]).

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5).
% library(json), NOT the legacy library(http/json) alias — the alias does not
% resolve in the WASM image (C6): fragment load + the JSON emits.
:- use_module(library(json), [json_read_dict/2, json_write_dict/2]).
:- use_module(library(apply), [foldl/4, maplist/3]).
:- use_module(library(aggregate), [aggregate_all/3]).
:- use_module(library(assoc),
              [assoc_to_keys/2, empty_assoc/1, get_assoc/3, list_to_assoc/2]).
:- use_module(library(lists),
              [append/3, member/2, numlist/3, reverse/2, selectchk/3]).
:- use_module(library(ordsets), [ord_union/3]).
% map_list_to_pairs/3, pairs_values/2 (fragment ordering, candidate rescore,
% and the greedy constructor's seed selection).
:- use_module(library(pairs), [map_list_to_pairs/3, pairs_values/2]).

% The shared metric layer: optimizer signals (word_checked_count/3,
% placed_bbox/4, word_cells/5), the answer footprint (word_letters/3), and
% cell_rc/4 for the constructor's bbox-growth scoring.
:- use_module(crosswordsmith(metrics),
              [ word_checked_count/3,
                word_letters/3,
                placed_bbox/4,
                word_cells/5,
                cell_rc/4
              ]).

% The shared substrate: legality core + MRV-inc branch step, layout build +
% numbering, geometry, and the strategy registry.
:- use_module(crosswordsmith(core),
              [ verbose_report/2,
                assign_clue_numbers/2,
                build_grid_rows/3,
                build_words/4,
                answer_meta_assoc/2,
                add_word_cells/3,
                cell_coord/3,
                init_gs/2,
                start_loc/4,
                start_locs/1,
                next_cell/4,
                fits_on_grid/4,
                assign_word/9,
                check_word_fits/5,
                find_intersecting_word/6,
                assign_words_inc/9,
                all_crossword/5,
                % the C1/C48 memo-hygiene seam: run ONCE at each of this
                % module's top-level search entries (see core.pl)
                reset_search_memos/0,
                default_strategy/1,
                pw_answer/2,
                pw_cells/2,
                pw_dir/2,
                pw_len/2,
                pw_start/2,
                pw_num/2,
                remove_x/3,
                shares_letter/2,
                check_unique_answers/1
              ]).


% ---------------------------------------------------------------------------
% Phase 1 - the scoring oracle (from-scratch capped reward over a placement).
% ---------------------------------------------------------------------------

% Default integer weights: WCap:WTail = 5:1  => epsilon = WTail/WCap = 0.2.
arrange_weights(5, 1).

% Per-word checking target: ceil(L/2) by default (== word_meets_half/2's
% threshold in metrics.pl; L is the placed dict's `len`). The §7.2 reachability
% escape hatch `--check-target N` lowers the ceiling to min(ceil(L/2), N) via
% check_target_override/1 (absent => the default). min/2 means N only ever
% LOWERS the target - an N above ceil(L/2) is a no-op - matching the spec intent
% ("lowering target below ceil(L/2) where half is unreachable").
% Invariant: at most ONE check_target_override/1 fact. set_check_target/1 (the
% exported setter - the ONLY sanctioned writer; the dynamic itself is
% module-private) always retractall's before assertz, so check_target/2's `->`
% (which commits to the first solution) reads the single current override. Do
% not assertz a second without retracting (P17).
:- dynamic check_target_override/1.
check_target(L, T) :-
    T0 is (L + 1) // 2,
    ( check_target_override(K) -> T is min(T0, K) ; T = T0 ).

% set_check_target(N): install the global --check-target ceiling (spec §7.2
% reachability escape hatch). -1 (the CLI's "not given" sentinel) clears any
% override; N>=1 lowers each word's target to min(ceil(L/2), N). Argument
% validation (and the CLI error message) stays with the driver; this setter
% fails on anything else rather than half-updating the state.
set_check_target(-1) :- !, retractall(check_target_override(_)).
set_check_target(N) :-
    integer(N), N >= 1,
    retractall(check_target_override(_)),
    assertz(check_target_override(N)).

% Per-word integer reward contribution:
%   WCap * min(checked, target) + WTail * checked
% checked(w) is reused from metrics.pl's word_checked_count/3 (W's cells that
% are also covered by a perpendicular placed word).
word_reward(WCap, WTail, PW, Placed, R) :-
    pw_len(PW, L),
    check_target(L, T),
    word_checked_count(PW, Placed, C),
    Capped is min(C, T),
    R is WCap * Capped + WTail * C.

% Total objective over a complete (or partial) placement - the from-scratch
% oracle Phase 2's incremental delta must match.
% NB named helper, not a yall lambda: this module never imports library(yall),
% so a `>>` here would be interpreted (meta-call + lambda copy per placed word
% per rescore) - the cost mode fill.pl's alpha_chars/2 comment documents
% avoiding. First-order closure instead.
layout_reward(WCap, WTail, Placed, Reward) :-
    foldl(add_word_reward(WCap, WTail, Placed), Placed, 0, Reward).

add_word_reward(WCap, WTail, Placed, PW, A0, A) :-
    word_reward(WCap, WTail, PW, Placed, R),
    A is A0 + R.

% How many placed words reach their cap (checked >= ceil(L/2)) - the
% reachability probe: if this is ~0, the cap is inert and the objective has
% degenerated to plain total-crossings.
cap_binding_count(Placed, Count) :-
    aggregate_all(count,
                  ( member(PW, Placed),
                    pw_len(PW, L), check_target(L, T),
                    word_checked_count(PW, Placed, C),
                    C >= T ),
                  Count).

% §7.2 (INV-3): when the cap binds on NO placed word, the objective has
% degenerated to plain total-crossings. That compromise - like a best-effort
% drop - is reported in the emitted payload's diagnostics property
% (json-output-spec §6.4) rather than on stderr: arrange output is best-effort
% by nature, so quality caveats are data for the consumer, not terminal noise.
cap_inert(Placed, CapInert) :-
    cap_binding_count(Placed, CB),
    (   CB =:= 0 -> CapInert = true ; CapInert = false ).

% The diagnostics.arrange sub-dict for an emitted layout: what the solver
% alone knows about how the result compromised. Recomputed from the final
% placement so every emit path (strict, best-effort, fragment, candidates)
% reports identically. dropped preserves input order. `seed` is added ONLY when
% --seed/--shuffle perturbed the search (provenance to reproduce the layout via
% --seed N); a deterministic run omits it, so its output is byte-unchanged.
arrange_diag_dict(Numbered, Words, Diag) :-
    arrange_weights(WCap, WTail),
    layout_reward(WCap, WTail, Numbered, Reward),
    cap_inert(Numbered, CapInert),
    dropped_answers(Words, Numbered, Dropped),
    Base = _{capInert: CapInert, dropped: Dropped, reward: Reward},
    (   current_search_seed(Seed)
    ->  put_dict(seed, Base, Seed, Diag)
    ;   Diag = Base
    ).

dropped_answers(Words, Placed, Dropped) :-
    maplist(pw_answer, Placed, PlacedAnswers),
    findall(A, ( member([A|_], Words), \+ memberchk(A, PlacedAnswers) ),
            Dropped).

% Human-readable cap-status suffix for the --verbose summary line (the
% payload's capInert, surfaced for interactive runs that opted in).
cap_status_note(Placed, Note) :-
    (   cap_inert(Placed, true)
    ->  Note = " (cap inert: objective = total-crossings; tune with --check-target)"
    ;   Note = ""
    ).


% ---------------------------------------------------------------------------
% Phase 2-3 - strict layout (construct + rescore + emit) with size framing.
%
% Per the Phase-1.5 gate the B&B search is descoped: strict mode constructs a
% complete placement via the reused MRV-inc path over the strict start corners
% (arrange_corners/1 - one per transpose-pair; see the note there and E-H1),
% rescores each with layout_reward/4, and emits the best. place-all-or-fail with
% budget-aware 3-outcome semantics.
% ---------------------------------------------------------------------------

% Construction budget (inferences) for finding ONE complete placement per
% corner. First solutions on puzzle-shaped inputs are 5-462 ms (well under).
arrange_budget(500_000_000).

% The start corners the STRICT construct path sweeps. This is deliberately a
% 2-corner subset of core.pl's start_locs/1 ([topleft_across, topleft_down,
% topright, bottomleft]): the four corners form two TRANSPOSE-PAIRS,
% {topleft_across, topleft_down} and {topright, bottomleft}. Transposing the
% grid (swap rows<->cols) maps across<->down while preserving letter order, so
% the mrv_inc search trees of a pair are exact transposes of each other - same
% solutions modulo transposition, and layout_reward is transpose-invariant (it
% counts checked cells), so a pair's first-solution rewards are identical.
% Sweeping all four therefore does ~2x redundant work in strict mode. We keep
% ONE corner per pair (topleft_across, topright) - the same two the reward-tie
% break already prefers (sort(1,@>=) is stable, so ties resolve to the earlier
% corner). Verified reward-equal AND literal-transpose across every ladder rung
% in benchmarks/workloads.pl (grids 9/15/21, 8..80 words). See experiment E-H1
% in docs/experiments.md.
%
% SCOPE: strict construct only. start_locs/1 stays 4-corner for the
% greedy/best-effort/candidates paths, where transposed layouts are legitimately
% DISTINCT candidates (candidate diversity), and for the fragment/enumerate
% paths. Do not route those through arrange_corners/1.
arrange_corners([topleft_across, topright]).

% Best-scoring complete placement over the strict start corners, rescored and
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
    % Top-level search entry: one memo reset per request (core.pl
    % reset_search_memos/0, C1/C48). The corner sweep below SHARES the
    % pair_crossings/answer_letters memos - the reset must never fire
    % per corner (construct_one bypasses find_crossword/6's seam).
    reset_search_memos,
    arrange_weights(WCap, WTail),
    arrange_corners(Locs),
    construct_corners(Locs, Words, GridLen, WCap, WTail, Budget, Results),
    findall(R-P, member(ok(R, P), Results), OKs),
    (   OKs = [_|_]
    ->  sort(1, @>=, OKs, [Reward-BestPlaced|_]),
        assign_clue_numbers(BestPlaced, Numbered),
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
        ( corner_search(Loc, Words, GridLen, P) -> Found = P
        ; Found = none ),
        Budget, Outcome),
    (   Outcome == inference_limit_exceeded -> Res = budget
    ;   Found == none                       -> Res = exhausted
    ;   layout_reward(WCap, WTail, Found, R), Res = ok(R, Found)
    ).

% One strict corner: core's mrv_inc driver composed DIRECTLY (init_gs +
% start_loc + assign_words_inc, the same primitives construct_from_seed
% composes for the fragment path), NOT via find_crossword/6 - that entry runs
% reset_search_memos/0 on every call (it is its own top-level seam), which
% would flush the memo tables between corners. arrange_best_layout/6 owns
% this request's single reset; the corners share the memos (C1 consolidation:
% exactly one reset per external call).
corner_search(Loc, Words, GridLen, PlacedWords) :-
    init_gs(GridLen, G1),
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words_inc(Words, [], none, GridLen, StartNum, StartDir,
                     G1, _Grid, PlacedWords).

% Construct the best strict layout, emit it (framed by SizeMode) on stdout,
% report on stderr. On any non-placed outcome: report + fail with no stdout
% (the json-output-spec §6.3 "no output / non-zero exit" contract).
arrange_strict_solve(Words, GridLen, SizeMode) :-
    check_unique_answers(Words),
    arrange_best_layout(Words, GridLen, Numbered, Reward, Outcome),
    (   Outcome == placed
    ->  emit_arrange_diag(Numbered, Words, GridLen, SizeMode),
        length(Numbered, NP),
        cap_status_note(Numbered, Note),
        verbose_report("arrange: grid ~w, mode ~w, placed ~w, reward ~w~w~n",
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
% Per the Phase-1.5 result, best-effort is served by metrics.pl's greedy
% constructor (which drops words it cannot place), NOT a drop-branch on the
% strict DFS. Construct over the seed x start-corner sweep on the given grid,
% rescore, and pick lexicographically: most words placed, then highest reward.

% Best greedy layout on GridLen. The key score(NumPlaced, Reward) compared in
% standard term order ranks most-placed first, then highest reward.
arrange_best_effort(Words, GridLen, Numbered, Reward, NumPlaced, Dropped) :-
    % Top-level search entry: one memo reset per request (core.pl
    % reset_search_memos/0, C1/C48 - this greedy path never reset before, so
    % its tables grew without bound in persistent processes and its inference
    % counts depended on call history). The seed x corner sweep below SHARES
    % the memos - do not reset per greedy_construct.
    reset_search_memos,
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
    assign_clue_numbers(BestPlaced, Numbered).

% Construct the best-effort layout, emit it, report placed/dropped on stderr.
% Fails only when nothing at all is placeable (AC-ARR-2).
arrange_best_effort_solve(Words, GridLen, SizeMode) :-
    check_unique_answers(Words),
    (   arrange_best_effort(Words, GridLen, Numbered, Reward, NP, Dropped)
    ->  % the dropped set rides the payload's diagnostics (INV-3, AC-ARR-2);
        % the stderr summary is routine -> --verbose only
        emit_arrange_diag(Numbered, Words, GridLen, SizeMode),
        length(Dropped, ND),
        verbose_report(
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

% Plain emit: layout only, no diagnostics. This is the variant fill's
% emit_fill(max) delegates to - fill has no arrange-solver compromises to
% report, so its payload must not grow a diagnostics.arrange sub-object.
emit_arrange(Numbered, Words, GridLen, SizeMode) :-
    arrange_layout_dict(Numbered, Words, GridLen, SizeMode, Payload),
    current_output(Out),
    json_write_dict(Out, Payload),
    nl(Out).

% arrange's own emit: the layout dict plus the diagnostics property
% (json-output-spec §6.4) carrying the solver's quality caveats.
emit_arrange_diag(Numbered, Words, GridLen, SizeMode) :-
    arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, Payload),
    current_output(Out),
    json_write_dict(Out, Payload),
    nl(Out).

arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, Payload) :-
    arrange_layout_dict(Numbered, Words, GridLen, SizeMode, Layout),
    arrange_diag_dict(Numbered, Words, Diag),
    put_dict(diagnostics, Layout, _{arrange: Diag}, Payload).

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
    pw_answer(PW, Answer),
    pw_dir(PW, Dir),
    pw_num(PW, Num),
    pw_cells(PW, Cells),
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

% arrange_outcome(+Words, +GridLen, +Drop, +SizeMode, -Outcome)
% The browser/WASM envelope seam (wasm-sdk-strategy §3, DEC-7): same search as
% arrange_solve/4 but EMITS NOTHING and NEVER FAILS — it returns a tagged
% outcome carrying the diagnostics-bearing layout dict, so the caller
% (crosswordsmith_browser) can wrap every result in a discriminated envelope
% instead of collapsing the non-happy paths into a bare fail. Tags:
%   placed(Dict)       strict search placed a full interlock
%   not_proven         the inference budget elapsed first (feasibility open)
%   infeasible(Bad)    strict search completed with no full placement; Bad =
%                      the genuinely-uncrossable answers ([] = interlock/size)
%   best_effort(Dict)  best-effort placement (drops ride
%                      Dict.diagnostics.arrange.dropped)
%   nothing_placeable  best-effort could not place a single word
% Throws duplicate_answer(_) like the solve seams; Dict is built by
% arrange_diag_layout_dict/5 (identical to the emitted payload, framed by
% SizeMode) but never written, so the envelope nests it as a live dict (DEC-8).
arrange_outcome(Words, GridLen, strict, SizeMode, Outcome) :-
    !,
    check_unique_answers(Words),
    arrange_best_layout(Words, GridLen, Numbered, _Reward, Outcome0),
    (   Outcome0 == placed
    ->  % once/1: arrange_layout_dict/5's fixed|max clauses leave a spurious
        % choicepoint (no arg-4 indexing); the dict is single-valued.
        once(arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, Dict)),
        Outcome = placed(Dict)
    ;   Outcome0 == not_proven
    ->  Outcome = not_proven
    ;   unplaceable_words(Words, Bad),
        Outcome = infeasible(Bad)
    ).
arrange_outcome(Words, GridLen, best_effort, SizeMode, Outcome) :-
    !,
    check_unique_answers(Words),
    (   arrange_best_effort(Words, GridLen, Numbered, _Reward, _NP, _Dropped)
    ->  once(arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, Dict)),
        Outcome = best_effort(Dict)
    ;   Outcome = nothing_placeable
    ).


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
% Everything here reuses core.pl's legality core (assign_word and the
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
    init_gs(GridLen, G0),
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
                    GridIn, PW, GridOut)
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
% generically. fragment_conflict/7 always throws. It receives the gs/2 bundle and
% reads the LETTER grid (a boundary-cell violation has no letter to report and
% falls through to the generic branch).
fragment_conflict(Answer, Letters, Start, Dir, WLen, GridLen, gs(Grid, _BGrid)) :-
    word_cells(Start, Dir, WLen, GridLen, Cells),
    (   clashing_cell(Cells, Letters, Grid, Cell, Existing, Wanted)
    ->  cell_coord(GridLen, Cell, RC),
        throw(error(fragment_letter_clash(Answer, RC, Existing, Wanted), _))
    ;   throw(error(fragment_illegal_pin(Answer), _))
    ).

% First cell whose grid letter is already fixed to something the word cannot
% supply (an overlap conflict). Cells and Letters are the parallel run.
clashing_cell([Cell|_], [L|_], Grid, Cell, Existing, L) :-
    arg(Cell, Grid, Existing),
    nonvar(Existing),               % a filled cell (an empty cell is an unbound var)
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
    % Top-level search entry: one memo reset per request (core.pl
    % reset_search_memos/0, C1/C48 - this path reaches assign_words_inc
    % directly and never reset before).
    reset_search_memos,
    arrange_weights(WCap, WTail),
    seed_from_fragment(Frags, Words, GridLen, SeededPlaced, SeededGrid, Remaining),
    arrange_budget(Budget),
    construct_fragment_one(Remaining, SeededPlaced, SeededGrid, GridLen, Budget, Res),
    (   Res = ok(AllPlaced)
    ->  layout_reward(WCap, WTail, AllPlaced, Reward),
        assign_clue_numbers(AllPlaced, Numbered),
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
        ( construct_from_seed(Remaining, SeededPlaced, SeededGrid, GridLen, P)
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
    % Top-level search entry: one memo reset per request (core.pl
    % reset_search_memos/0, C1/C48).
    reset_search_memos,
    arrange_weights(WCap, WTail),
    seed_from_fragment(Frags, Words, GridLen, SeededPlaced, SeededGrid, Remaining),
    greedy_loop(Remaining, SeededPlaced, GridLen, SeededGrid, FinalPlaced, DroppedEntries),
    layout_reward(WCap, WTail, FinalPlaced, Reward),
    length(FinalPlaced, NumPlaced),
    findall(A, member([A|_], DroppedEntries), Dropped),
    assign_clue_numbers(FinalPlaced, Numbered).


% --- fragment solve entry points (emit + report) ---------------------------
% Mirror arrange_strict_solve/arrange_best_effort_solve: emit on stdout, report
% on stderr; on any non-placed outcome report and fail with no stdout layout.

arrange_fragment_solve(Words, Frags, GridLen, strict, SizeMode) :-
    check_unique_answers(Words),
    arrange_fragment_strict(Words, Frags, GridLen, Numbered, Reward, Outcome),
    (   Outcome == placed
    ->  emit_arrange_diag(Numbered, Words, GridLen, SizeMode),
        length(Numbered, NP),
        cap_status_note(Numbered, Note),
        verbose_report(
               "arrange: fragment-seeded, grid ~w, mode ~w, placed ~w, reward ~w~w~n",
               [GridLen, SizeMode, NP, Reward, Note])
    ;   arrange_fragment_report_failure(Outcome, Words, GridLen),
        fail
    ).
arrange_fragment_solve(Words, Frags, GridLen, best_effort, SizeMode) :-
    check_unique_answers(Words),
    (   arrange_fragment_best_effort(Words, Frags, GridLen, Numbered, Reward, NP, Dropped)
    ->  % same contract as arrange_best_effort_solve: drops ride the payload's
        % diagnostics, the stderr summary is --verbose only
        emit_arrange_diag(Numbered, Words, GridLen, SizeMode),
        length(Dropped, ND),
        verbose_report(
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
              pw_answer(PW, A), pw_start(PW, S), pw_dir(PW, D),
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
    % Top-level search entry: one memo reset per request (core.pl
    % reset_search_memos/0, C1/C48). The pool's K greedy constructions
    % SHARE the memos - the reset lives here, never inside the pool sweep.
    reset_search_memos,
    arrange_candidate_pool(Words, GridLen, DropContract, Pool),
    Pool = [_|_],
    length(Words, Total),
    candidate_tau_pct(TauPct),
    pick_diverse(Pool, TauPct, Total, K, Picked),
    findall(N, ( member(c(P, _), Picked), assign_clue_numbers(P, N) ), Numbered),
    length(Numbered, Returned).

% Emit the candidates as a JSON ARRAY of canonical layout dicts (each element is
% itself a valid fragment for re-ingestion). Single-layout emit stays a bare
% object; the array shape is what `--candidates` opts into.
emit_candidates(NumberedLayouts, Words, GridLen, SizeMode) :-
    maplist(candidate_dict(Words, GridLen, SizeMode), NumberedLayouts, Dicts),
    current_output(Out),
    json_write_dict(Out, Dicts),
    nl(Out).
% Each candidate carries its own diagnostics: quality caveats are per-layout
% (rewards differ; under a drops contract the dropped sets can too).
candidate_dict(Words, GridLen, SizeMode, Numbered, Dict) :-
    arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, Dict).

% Solve + emit the candidates array on stdout, report requested/returned/tau on
% stderr (INV-3: a short return is reported, never silent). Fails (no stdout) if
% nothing is placeable.
arrange_candidates_solve(Words, GridLen, DropContract, SizeMode, K) :-
    check_unique_answers(Words),
    (   arrange_candidates(Words, GridLen, DropContract, K, Layouts, Returned),
        Returned > 0
    ->  emit_candidates(Layouts, Words, GridLen, SizeMode),
        candidate_tau_pct(TauPct),
        % fewer-than-K is a reportable compromise (AC-ARR-7): unconditional;
        % the got-what-you-asked-for summary is routine -> --verbose only
        (   Returned < K
        ->  format(user_error,
                   "arrange: ~w candidate(s) requested, ~w returned (tau ~w%); \c
fewer >=tau-distinct layouts exist~n",
                   [K, Returned, TauPct])
        ;   verbose_report(
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
% corners. Reuses core.pl's all_crossword/5 with the production default
% strategy, so the count matches the old `--all` exactly (strategies only
% reorder the same search tree). The branch step / legality core are untouched.
arrange_enumerate(Words, GridLen, Num) :-
    default_strategy(Strat),
    % No reset_search_memos/0 here: all_crossword/5 runs the count through
    % find_crossword/6, whose own top-level seam resets once per call -
    % resetting here too would double-fire (C1: one reset per external call).
    all_crossword(Strat, GridLen, Words, _StartLoc, Num).

% Count + print the enumeration on stdout (a bare integer), as the old CLI did.
arrange_enumerate_solve(Words, GridLen) :-
    check_unique_answers(Words),
    arrange_enumerate(Words, GridLen, Num),
    current_output(Out),
    writeln(Out, Num).


% ---------------------------------------------------------------------------
% The greedy density constructor (moved here from metrics.pl in Phase 3 of the
% source-structure migration - arrange is its only consumer, and its project
% dependencies are core's search primitives, not the metrics).
%
% arrange's best-effort and candidates path: greedily place the
% highest-(new-crossings, then lowest-bbox-growth) legal placement each step,
% dropping words that cannot be placed. One pass, no backtracking, cut-free.
% ---------------------------------------------------------------------------

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

neg_answer_len(Entry, NL) :- word_letters(Entry, _, WLen), NL is -WLen.

% --- greedy construction from one start location --------------------------

greedy_construct(Words, GridLen, Loc, Seed, Placed, Dropped) :-
    init_gs(GridLen, G0),
    start_loc(Loc, GridLen, StartNum, StartDir),
    remove_x(Seed, Words, Rest),
    seed_word(Seed, StartNum, StartDir, GridLen, G0, SeedPW, G1),
    greedy_loop(Rest, [SeedPW], GridLen, G1, Placed, Dropped).

% Seed: place the chosen word at the fixed start cell/direction (no crossings).
seed_word(Entry, Start, Dir, GridLen, GIn, PW, GOut) :-
    word_letters(Entry, Letters, WLen),
    fits_on_grid(Dir, Start, WLen, GridLen),
    Entry = [Word|_],
    assign_word(Word, Letters, WLen, Start, Dir, GridLen, GIn, PW, GOut).

% Place the globally best-scoring placeable word, repeat; drop the rest. The
% construction is cut-free: instead of an `( Best -> place ; stop )` if-then-else,
% the next move (or `none`) is reified as a term and dispatched on the mutually
% exclusive clause heads of best_move/5 ([] vs [_|_]). No !, ->, or \+ here or in
% word_best_placement below; only findall/sort/arithmetic and the winner's one
% assign_word remain. (This replaced an equivalent cut-based version; it is
% identical-output and faster - see docs/cryptic-layout-spec.md v1b.1.)
greedy_loop(Remaining, Placed, GridLen, Grid, FinalPlaced, Dropped) :-
    next_move(Remaining, Placed, GridLen, Grid, Move),
    apply_move(Move, Remaining, Placed, GridLen, FinalPlaced, Dropped).

apply_move(none, Remaining, Placed, _GridLen, Placed, Remaining).
apply_move(move(Answer, NewPW, NewGrid), Remaining, Placed, GridLen, FinalPlaced, Dropped) :-
    % Answer is the winner's ground answer ATOM (atoms survive next_move's
    % findall copy identically); the full entry in Remaining can be non-ground
    % (a .pl fixture's [Answer, _{...}] has an unbound dict tag), so an
    % ==-keyed remove_x on a copied Entry would drop nothing and greedy_loop
    % would re-offer the just-placed word forever. Key the removal on the
    % answer atom instead (answers are unique, check_unique_answers/1);
    % selectchk keeps the original, uncopied tail. (Same term-copy footgun the
    % MRV path avoids via map_list_to_pairs - see core.pl select_inc.)
    selectchk([Answer|_], Remaining, Remaining1),
    greedy_loop(Remaining1, [NewPW|Placed], GridLen, NewGrid, FinalPlaced, Dropped).

% The best placeable word as move(Answer,PW,Grid), or `none` when nothing fits.
% Among remaining words with a legal placement, the one whose best placement
% scores highest (most new crossings, then least bbox growth). The [] vs [_|_]
% dispatch in best_move/5 is what lets greedy_loop avoid an if-then-else.
%
% The findall collects only a lightweight GROUND descriptor per word
% (Score-best(Answer,Letters,WLen,Start,Dir)) - legality inside
% word_best_placement is the pure probe check_word_fits/5, so no mutated grid
% is materialized or snapshotted per candidate; assign_word/9 runs ONCE, for
% the winner, in best_move. (This replaced a version whose findall templates
% snapshotted the whole gs/2 bundle - 2 x N^2 args - per legal placement AND
% per word; identical output, the probe/assign accept sets and the stable
% first-tie-wins sorts are unchanged.) sort(1,@>=) compares KEYS only and is
% stable, so the winner is the FIRST-generated candidate with the maximal
% score, exactly as before.
next_move(Remaining, Placed, GridLen, Grid, Move) :-
    placed_bbox(Placed, GridLen, BBox, _),
    findall(Score-Best,
            ( member(Entry, Remaining),
              word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, Best) ),
            Cands),
    best_move(Cands, Placed, GridLen, Grid, Move).

best_move([], _Placed, _GridLen, _Grid, none).
best_move([C|Cs], _Placed, GridLen, Grid, move(Answer, PW, G1)) :-
    sort(1, @>=, [C|Cs], [_-best(Answer, Letters, WLen, Start, Dir)|_]),
    % Realize the winning placement. check_word_fits/5 accepted exactly this
    % (Letters, Start, Dir) against exactly this Grid - the probe binds
    % nothing, so the grid is untouched since - and it accepts iff
    % assign_word/9 accepts (core.pl check_word_fits doc), so this call
    % cannot fail. G1 is Grid itself, mutated in place (letter cells bound,
    % boundary cells marked, undone by the trail on backtracking), not a copy.
    assign_word(Answer, Letters, WLen, Start, Dir, GridLen, Grid, PW, G1).

% Best legal placement of one word on the current grid, as the lightweight
% ground descriptor best(Answer,Letters,WLen,Start,Dir) - cut-free: legality
% (the pure probe check_word_fits/5, which binds no cell) is folded INTO the
% findall generator so only legal placements are collected (keyed by the
% density score), and the head of the @>=-sorted list is the best - no
% first-solution `!` (see spec v1b.1).
%
% This used to run assign_word/9 per legal candidate and let findall snapshot
% the mutated gs/2 bundle (2 x N^2 args) per placement - the largest copy-term
% cost in the codebase. check_word_fits/5 answers "would assign_word accept
% this candidate?" with NO side effects and the SAME accept/reject set
% (core.pl:731-749), so the collected Key list - and therefore the stable-sort
% winner - is unchanged, while nothing grid-sized is copied. The winner's
% assign_word runs once, in best_move above.
%
% Order matters: placement_key MUST run BEFORE check_word_fits. Not for grid
% state (unlike assign_word, the probe never binds a cell), but because
% find_intersecting_word/6 can hand us an off-grid Start < 1 on which a raw
% arg/3 read would THROW; both placement_key (crossing_count's Start >= 1
% guard) and check_word_fits reject it by failure, and keeping the scorer
% first preserves the old generator's exact evaluation order.
word_best_placement(Entry, Placed, GridLen, Grid, BBox, Score, Best) :-
    word_letters(Entry, Letters, WLen),
    Entry = [Answer|_],
    % Grid is the gs/2 bundle; the scorer reads letters only, so it gets the
    % letter grid, while the legality probe reads both grids of the bundle.
    Grid = gs(LGrid, _BGrid),
    findall(Key-(Start-Dir),
            ( find_intersecting_word(Letters, WLen, Placed, GridLen, Start, Dir),
              placement_key(Letters, Start, Dir, WLen, GridLen, LGrid, BBox, Key),
              check_word_fits(Letters, Start, Dir, GridLen, Grid) ),
            Keyed),
    Keyed = [_|_],
    sort(1, @>=, Keyed, [Score-(BestStart-BestDir)|_]),
    Best = best(Answer, Letters, WLen, BestStart, BestDir).

% Density score for a candidate placement: crossings dominate, bbox-growth
% breaks ties (smaller is better). 10000 > any plausible bbox area.
placement_key(Letters, Start, Dir, WLen, GridLen, Grid, BBox, Key) :-
    crossing_count(Letters, Start, Dir, GridLen, Grid, Crossings),
    word_cells(Start, Dir, WLen, GridLen, Cells),
    bbox_growth(BBox, Cells, GridLen, Growth),
    Key is Crossings * 10000 - Growth.

crossing_count(Letters, Start, Dir, GridLen, Grid, Count) :-
    % Reject an off-grid (underflow) start before indexing the grid: since
    % word_best_placement now scores BEFORE assign_word, find_intersecting_word/6
    % can hand us a Start < 1 that assign_word would go on to reject. With the
    % var-cell grid, arg/3 THROWS on a negative index (the old get_assoc/3 failed
    % silently), so fail here instead - the candidate is dropped either way. With
    % Start >= 1 and a fitted run, every cell index stays >= 1 (greedy-path only).
    Start >= 1,
    cc_(Letters, Start, Dir, GridLen, Grid, 0, Count).
cc_([], _, _, _, _, A, A).
cc_([L|Ls], Num, Dir, GridLen, Grid, A0, Count) :-
    % A crossing = the cell already holds this exact letter. `==` is a pure test
    % (never binds), so scoring a candidate leaves the grid untouched; an empty
    % (unbound) cell is never == an atom, so it counts as no crossing.
    arg(Num, Grid, Cell),
    ( Cell == L -> A1 is A0 + 1 ; A1 = A0 ),
    next_cell(Dir, Num, GridLen, Num2),
    cc_(Ls, Num2, Dir, GridLen, Grid, A1, Count).

% Bounding-box growth a candidate placement would cause (0 when it stays
% inside the current bbox). cell_rc/4 stays in metrics.pl (lint uses it too).
bbox_growth(bbox(MinR, MaxR, MinC, MaxC), NewCells, GridLen, Growth) :-
    OldArea is (MaxR - MinR + 1) * (MaxC - MinC + 1),
    foldl(extend_cell(GridLen), NewCells, b(MinR, MaxR, MinC, MaxC), b(R0, R1, C0, C1)),
    NewArea is (R1 - R0 + 1) * (C1 - C0 + 1),
    Growth is NewArea - OldArea.

extend_cell(GridLen, Cell, b(MinR, MaxR, MinC, MaxC), b(MinR2, MaxR2, MinC2, MaxC2)) :-
    cell_rc(Cell, GridLen, R, C),
    MinR2 is min(MinR, R), MaxR2 is max(MaxR, R),
    MinC2 is min(MinC, C), MaxC2 is max(MaxC, C).
