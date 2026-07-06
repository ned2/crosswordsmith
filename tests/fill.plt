% tests/fill.plt - plunit suite for fill.pl (grid-first auto-fill, §8.4).
%
% Assumes the project is loaded (via load.pl) by the runner before this file.
% fill exports only fill_solve/4, so these white-box tests reach its internals
% as crosswordsmith_fill:Pred(...) (migration plan, Resolved Decision 3).
% Covers slot derivation + the SHARED cell-variable invariant (the crossing
% constraint), the fill search (AC-FILL-1), seeds (AC-FILL-2), determinism
% (AC-FILL-3), and the no-fill report. The byte-exact fill is pinned by
% tests/golden/fill_3.json.

:- use_module(library(plunit)).
:- use_module(library(fastrw)).   % F-L2 tests hand-craft artifacts to prove refusals
:- use_module(library(assoc)).
% lists/apply: test-body helpers; explicit so the suite also runs under
% autoload(false) (P11/C5).
:- use_module(library(lists)).
:- use_module(library(apply)).

% Run a fill and return the across/down answer lists (deterministic).
do_fill(GridFile, SeedFile, DictFile, Across, Down) :-
    crosswordsmith_fill:fill_grid(GridFile, _Size, Slots, _CellVar),
    ( SeedFile == none -> SeededKeys = [] ; crosswordsmith_fill:apply_seeds(SeedFile, Slots, SeededKeys) ),
    exclude(crosswordsmith_fill:seeded_slot(SeededKeys), Slots, SearchSlots),
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    crosswordsmith_fill:fill_attempt(SearchSlots, Slots, DictByLen, Index, filled, Numbered, _),
    findall(A, ( member(W, Numbered), pw_dir(W, across), pw_answer(W, A) ), Across),
    findall(A, ( member(W, Numbered), pw_dir(W, down),   pw_answer(W, A) ), Down).


:- begin_tests(fill).

% --- slots + the shared-variable invariant (the crossing constraint) ---------
% The 3x3 has 3 across + 3 down slots.
test(fill_derives_slots) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', 3, Slots, _),
    length(Slots, 6).

% REGRESSION: a crossing cell's variable must be SHARED between its across and
% down slot (cell 1 starts both row-0-across and col-0-down). A findall-copied
% or yall-lambda-copied template breaks this and the fill is inconsistent.
test(fill_crossing_cells_are_shared) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', 3, Slots, _),
    once(( member(slot(1, across, [1,2,3], [A1|_]), Slots) )),
    once(( member(slot(1, down,   [1,4,7], [D1|_]), Slots) )),
    A1 == D1.

% --- the fill (AC-FILL-1) ----------------------------------------------------
% A complete, consistent fill: every across is a row, every down is a column
% (deterministic - pinned to the exact result).
test(fill_produces_consistent_square) :-
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', Across, Down),
    Across == ['CAT', 'ORE', 'WED'],
    Down   == ['COW', 'ARE', 'TED'].

% No dictionary word of the slot length -> infeasible (reported, not silent).
test(fill_infeasible_when_no_matching_words) :-
    tmp_file_stream(text, F, S), write(S, "FOUR\nFIVE\n"), close(S),   % 4-letter only
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict(F, DictByLen, Index),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DictByLen, Index, Outcome, _, _),
    delete_file(F),
    Outcome == infeasible.

% Regression for P2: fill_attempt/8 must NOT swallow a genuine error from
% fill_search and report it as `infeasible`. call_with_inference_limit/3 handles
% the budget itself and re-throws real exceptions; infeasibility is a search
% FAILURE (R == exhausted), never a throw. Feeding a malformed pattern index (a
% non-assoc) makes get_assoc/3 raise a type_error, which must propagate rather
% than be masked as infeasibility. The old broad catch/3 returned `infeasible`.
test(fill_propagates_genuine_error,
     [throws(error(type_error(btree, _), _))]) :-
    crosswordsmith_fill:fill_attempt([slot(a, across, 3, [_,_,_])], [slot(a, across, 3, [_,_,_])],
                 [], not_an_assoc, 1_000_000, _Outcome, _Numbered, _InputWords).

% --- seeds (AC-FILL-2) -------------------------------------------------------
% Pinning COW across row 0 forces the transpose square; the pin appears at its
% slot, and the unseeded default (CAT across) is overridden.
test(fill_respects_seed) :-
    do_fill('fixtures/fill_grid_3.json', 'fixtures/fill_seed_3.json',
            'fixtures/wordlist_sample.txt', Across, _Down),
    Across == ['COW', 'ARE', 'TED'],
    Across \== ['CAT', 'ORE', 'WED'].

% R2 (revamp audit): a seed is a HARD PIN, not required to be a dictionary word
% (a setter's theme/own word). With CAT pinned across row 0 and a dictionary
% that OMITS CAT, the grid still fills around the pin (CAT/ORE/WED across) - the
% old code reported the user's own pinned slot as unfillable and failed.
test(fill_seed_need_not_be_in_dict) :-
    tmp_file_stream(text, F, S),
    write(S, "ORE\nWED\nCOW\nARE\nTED\n"), close(S),   % the sample wordlist minus CAT
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    foldl(crosswordsmith_fill:apply_seed(Slots), [frag('CAT', across, 1, [1, 2, 3])], [], SeededKeys),
    exclude(crosswordsmith_fill:seeded_slot(SeededKeys), Slots, SearchSlots),
    crosswordsmith_fill:load_dict(F, DictByLen, Index),
    crosswordsmith_fill:fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, _),
    delete_file(F),
    Outcome == filled,
    findall(A, ( member(W, Numbered), pw_dir(W, across), pw_answer(W, A) ), Across),
    Across == ['CAT', 'ORE', 'WED'].

% A seed whose answer matches no slot of the grid is rejected.
test(fill_seed_no_slot_throws, [throws(error(fill_seed_no_slot('CAT'), _))]) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    % CAT down at col 0 cells [1,4,7] is a real slot; ask for it as a 2-cell
    % run that no slot has -> no match.
    apply_seeds_frags(Slots, [frag('CAT', across, 2, [2, 3])]).

% --- seed answers vs the search's no-duplicate rule (fill-seed-pin-crash-fix) --
% fill_grid_split3 (mask ["...","###","..."]) has two INDEPENDENT 3-cell across
% slots (the blocked middle row leaves no down run >= 2), so the searched row is
% filled independently of the row-0 seed. Before the fix the search's Used set
% started empty (seed answers were exempt from the no-duplicate rule), so it
% could re-place the seed word and the duplicate blew up the emit with an
% uncaught domain_error(unique_key_pairs).

% Trigger A on a SOLVABLE grid: with dict [COW,PIG] and COW pinned, the searched
% slot must take the non-duplicate alternative PIG (the search used to pick COW
% - alphabetically first - and crash). Asserting the exact answer set guards the
% CORRECTNESS of the fix, incl. the char-list-vs-atom trap (an atom-shaped Used
% silently no-ops the dedup and this test would see COW twice), not just
% no-throw.
test(fill_seed_answer_reused_takes_alternative) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', _Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_top.json', Slots, SeededKeys),
    exclude(crosswordsmith_fill:seeded_slot(SeededKeys), Slots, SearchSlots),
    crosswordsmith_fill:load_dict('fixtures/dict_cow_pig.txt', DictByLen, Index),
    crosswordsmith_fill:fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, Numbered, _),
    Outcome == filled,
    findall(A, ( member(W, Numbered), pw_answer(W, A) ), Answers0),
    msort(Answers0, Answers),
    Answers == ['COW', 'PIG'].

% Trigger A when NO non-duplicate fill exists (the dictionary holds only the
% seed word): a clean infeasible outcome. The body must not throw - plunit
% fails an uncaught error, so this is the direct crash-regression guard (it
% threw domain_error(unique_key_pairs) before the fix).
test(fill_seed_answer_only_word_is_infeasible) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', _Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_top.json', Slots, SeededKeys),
    exclude(crosswordsmith_fill:seeded_slot(SeededKeys), Slots, SearchSlots),
    crosswordsmith_fill:load_dict('fixtures/dict_cow.txt', DictByLen, Index),
    crosswordsmith_fill:fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, _, _),
    Outcome == infeasible.

% Trigger B: two seeds pinning the SAME answer are rejected up front with the
% clean hooked error (neither slot is searched, so the Used-set fix cannot
% catch this shape; reported before searching, like fill_seed_no_slot/clash).
test(fill_duplicate_seed_answers_rejected,
     [throws(error(fill_seed_duplicate('COW'), _))]) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', _Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_both.json', Slots, _).

% Defense in depth: the fill emit boundary re-runs core's unique-answers check
% (the solve path already does, in crossword/4), so a duplicate that ever
% reaches emit again reports the clean hooked duplicate_answer error - never
% the raw domain_error(unique_key_pairs) from answer_meta_assoc/2.
test(fill_emit_boundary_rejects_duplicate_answers,
     [throws(error(duplicate_answer('COW'), _))]) :-
    crosswordsmith_fill:emit_fill([], [['COW'], ['COW']], 3, fixed).

% R11 (revamp audit): AC-FILL-1 "not proven within budget". A budget too small
% to complete the search yields Outcome==not_proven - distinct from infeasible
% (a genuinely 0-candidate slot) and from filled.
test(fill_budget_exhausted_not_proven) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt', DictByLen, Index),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DictByLen, Index, 100, Outcome, Numbered, _),
    Outcome == not_proven,
    Numbered == [].

% R11 (revamp audit): AC-FILL-4 - a produced fill re-validates as a legal
% layout. Linting the fill output gives a clean structural verdict: the 3x3
% full grid is every-cell-checked, connected, and symmetric, so it PASSes the
% strictest structural profile (american). Letters do not affect structure.
test(fill_revalidates_under_lint) :-
    lint_load('tests/golden/fill_3.json', GridLen, Placed),
    lint_run(Placed, GridLen, american, false, Report),
    get_dict(verdict, Report, V),
    V == 'PASS'.

% R6 (revamp audit): when a start cell begins both an across and a down slot,
% select_mrv must expand the slot that actually has the fewest candidates, not
% whichever shares the start and was generated first. Here start 1's down slot
% is fully constrained to CAT (1 candidate) while its across slot has 3 (C__);
% the down slot (the minimum) must be chosen.
test(select_mrv_recovers_correct_direction_on_tie) :-
    tmp_file_stream(text, F, S), write(S, "CAT\nCOW\nCUB\nDOG\n"), close(S),
    crosswordsmith_fill:load_dict(F, DictByLen, Index), delete_file(F),
    Across = slot(1, across, [1, 2, 3], ['C', _, _]),       % C__ -> CAT/COW/CUB (3)
    Down   = slot(1, down,   [1, 4, 7], ['C', 'A', 'T']),   % CAT (1, the minimum)
    crosswordsmith_fill:select_mrv([Across, Down], DictByLen, Index, Best, _Rest, Cands),
    Best = slot(1, down, _, _),
    Cands == [['C', 'A', 'T']].

% P13: select_mrv/6 must succeed DETERMINISTICALLY - the winning slot's Start+Dir
% is unique, so once(select/3) prunes the spurious choicepoint plain select/3
% would otherwise leave. The winner (Down, 1 candidate) is placed FIRST in the
% slot list so a stray choicepoint from select/3 (scanning the rest of the list)
% is possible - that is exactly what once/1 removes.
test(select_mrv_leaves_no_choicepoint) :-
    tmp_file_stream(text, F, S), write(S, "CAT\nCOW\nCUB\nDOG\n"), close(S),
    crosswordsmith_fill:load_dict(F, DictByLen, Index), delete_file(F),
    Across = slot(1, across, [1, 2, 3], ['C', _, _]),
    Down   = slot(1, down,   [1, 4, 7], ['C', 'A', 'T']),
    crosswordsmith_fill:select_mrv([Down, Across], DictByLen, Index, Best, _Rest, _Cands),
    Best = slot(1, down, _, _),
    deterministic(Det),
    Det == true.

% P3: candidate_count/4 (the count-only MRV metric, no word materialization) must
% agree with length(candidates/4) on BOTH branches - `all` (no cell bound) and
% idx (some cells bound). select_mrv orders slots by these counts, so if the two
% ever diverged the MRV choice would silently differ from the true candidate set.
test(candidate_count_matches_candidates) :-
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt', DictByLen, Index),
    % all-unbound (`all` branch): count == number of length-3 words materialized
    length(Free, 3),
    crosswordsmith_fill:candidates(Free, DictByLen, Index, CAll), length(CAll, NAll),
    crosswordsmith_fill:candidate_count(Free, DictByLen, Index, NAll),
    NAll > 0,
    % one cell bound to 'C' (idx branch): only CAT, COW match
    Bound = ['C', _, _],
    crosswordsmith_fill:candidates(Bound, DictByLen, Index, CB), length(CB, NB),
    crosswordsmith_fill:candidate_count(Bound, DictByLen, Index, NB),
    NB =:= 2.

% --- determinism (AC-FILL-3) -------------------------------------------------
test(fill_deterministic) :-
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A1, D1),
    do_fill('fixtures/fill_grid_3.json', none, 'fixtures/wordlist_sample.txt', A2, D2),
    A1 == A2, D1 == D2.

% --- F-L2: precomputed index artifact ----------------------------------------
% Write a temp dict, return its path (the caller deletes it).
tmp_dict(Text, Path) :- tmp_file_stream(text, Path, S), write(S, Text), close(S).
tmp_path(Path) :- tmp_file_stream(binary, Path, S), close(S).

% EQUIVALENCE LEMMA (the core F-L2 correctness claim): the DictByLen buckets and
% the assoc Index reconstructed from a saved artifact are ==-IDENTICAL to a fresh
% load_dict of the same file. The whole speed-up rests on this, so pin it
% directly (scaled up to 10k/50k/172k out-of-band by the results-doc harness).
test(fill_index_roundtrip_identical) :-
    tmp_dict("CAT\nCOW\nCUB\nDOG\nORE\nWED\nCUBIT\nTOMMY\n", DF), tmp_path(AF),
    crosswordsmith_fill:fill_save_index(DF, AF),
    crosswordsmith_fill:load_dict(DF, DBL0, Idx0),
    crosswordsmith_fill:fill_load_index(AF, none, DBL1, Idx1),
    delete_file(DF), delete_file(AF),
    DBL1 == DBL0,
    Idx1 == Idx0.

% An artifact-mode fill produces the SAME layout as the raw path (the CLI
% identity oracle proves this on all 11 rungs; this is the unit-level witness).
test(fill_index_fill_matches_raw) :-
    tmp_path(AF),
    crosswordsmith_fill:fill_save_index('fixtures/wordlist_sample.txt', AF),
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    crosswordsmith_fill:fill_load_index(AF, none, DBL, Idx),
    delete_file(AF),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    findall(A, ( member(W, Numbered), pw_dir(W, across), pw_answer(W, A) ), Across),
    Across == ['CAT', 'ORE', 'WED'].

% Passing --dict that MATCHES the artifact's source verifies the SHA-256 and
% loads (integrity check on the happy path).
test(fill_index_hash_ok_when_dict_matches) :-
    tmp_dict("CAT\nCOW\nCUB\nDOG\n", DF), tmp_path(AF),
    crosswordsmith_fill:fill_save_index(DF, AF),
    crosswordsmith_fill:fill_load_index(AF, DF, _DBL, _Idx),
    delete_file(DF), delete_file(AF).

% A dict whose bytes differ from the artifact's source is REFUSED (no silent
% rebuild, no stale fill) - the staleness gate.
test(fill_index_hash_mismatch_refused,
     [throws(error(fill_index_hash(_, _, _), _))]) :-
    tmp_dict("CAT\nCOW\nCUB\nDOG\n", DF1), tmp_dict("CAT\nCOW\nCUB\nDOGS\n", DF2),
    tmp_path(AF),
    crosswordsmith_fill:fill_save_index(DF1, AF),
    setup_call_cleanup(true,
                       crosswordsmith_fill:fill_load_index(AF, DF2, _, _),
                       ( delete_file(DF1), delete_file(DF2), delete_file(AF) )).

% An unrecognised SCHEMA version is refused (forward/backward-compat guard; this
% is also where F-H2's masks bump the version). Hand-craft the artifact so the
% test does not depend on ever shipping a second version.
test(fill_index_version_mismatch_refused,
     [throws(error(fill_index_version(999, _), _))]) :-
    empty_assoc(E), tmp_path(AF),
    setup_call_cleanup(open(AF, write, S, [type(binary)]),
                       fast_write(S, fill_index(999, [swi_version('x')], E, E)),
                       close(S)),
    setup_call_cleanup(true,
                       crosswordsmith_fill:fill_load_index(AF, none, _, _),
                       delete_file(AF)).

% An artifact built by a different SWI-Prolog (binary format is version-bound)
% is refused even without --dict.
test(fill_index_swi_mismatch_refused,
     [throws(error(fill_index_swi('0.0.0', _), _))]) :-
    empty_assoc(E), tmp_path(AF),
    % version 2 (current schema) so the SWI gate is reached, not the version gate.
    setup_call_cleanup(open(AF, write, S, [type(binary)]),
                       fast_write(S, fill_index(2, [swi_version('0.0.0'), dict_sha256(x)], E, E)),
                       close(S)),
    setup_call_cleanup(true,
                       crosswordsmith_fill:fill_load_index(AF, none, _, _),
                       delete_file(AF)).

% A file that is not a fill_index/4 artifact is refused, not misinterpreted.
test(fill_index_malformed_refused,
     [throws(error(fill_index_malformed(_), _))]) :-
    tmp_path(AF),
    setup_call_cleanup(open(AF, write, S, [type(binary)]),
                       fast_write(S, not_an_artifact(42)),
                       close(S)),
    setup_call_cleanup(true,
                       crosswordsmith_fill:fill_load_index(AF, none, _, _),
                       delete_file(AF)).

% A missing artifact file is refused with a clear error.
test(fill_index_missing_refused,
     [throws(error(fill_index_missing(_), _))]) :-
    crosswordsmith_fill:fill_load_index('/nonexistent/f-l2/no.idx', none, _, _).

% --- F-H2: v2 bitset masks (bignum counting; masks OPT-IN via --masks) --------
% The whole F-H2 correctness claim: for EVERY bound pattern, popcount of the mask
% AND-chain equals the length of the ordset intersection. A masks-mode v2
% artifact (fill_save_index/3 with masks(true)) must carry masks
% (Masks = masks(_)), and mask_count/4 must agree with index_intersection/4 +
% length on a broad deterministic sample at 10k. If these ever diverged the MRV
% counts would differ and the filled grid could change.
test(fill_index_v2_mask_count_matches_ordset) :-
    tmp_path(AF),
    crosswordsmith_fill:fill_save_index('fixtures/dict/enable_10k.txt', AF, [masks(true)]),
    crosswordsmith_fill:fill_load_index(AF, none, _DBL, Idx, Masks),
    delete_file(AF),
    Masks = masks(MA),                       % a masks-mode artifact MUST carry masks
    findall(L-B, sample_pattern(L, B), Patterns),
    Patterns = [_|_],                        % non-empty sample
    forall( member(Len-Bound, Patterns),
            ( crosswordsmith_fill:index_intersection(Bound, Len, Idx, OrdI),
              length(OrdI, OrdN),
              crosswordsmith_fill:mask_count(Bound, Len, MA, MaskN),
              OrdN =:= MaskN ) ).

% The DEFAULT v2 build (fill_save_index/2, no --masks) carries NO masks: the
% loader hands back Masks == none and counting runs the ordset kernel - the
% candidate_count/5 `none` clause agrees with the /4 reference on both branches
% (all-unbound and bound), and an artifact-mode fill completes. This pins the
% no-size-tax default path (the F-H2 follow-up's shipping shape).
test(fill_index_v2_default_no_masks_counts_via_ordsets) :-
    tmp_path(AF),
    crosswordsmith_fill:fill_save_index('fixtures/wordlist_sample.txt', AF),
    crosswordsmith_fill:fill_load_index(AF, none, DBL, Idx, Masks),
    delete_file(AF),
    Masks == none,                           % default v2 MUST NOT carry masks
    % `all` branch and bound branch agree with the /4 reference kernel
    length(Free, 3),
    crosswordsmith_fill:candidate_count(Free, DBL, Idx, NAll),
    crosswordsmith_fill:candidate_count(Masks, Free, DBL, Idx, NAll),
    Bound = ['C', _, _],
    crosswordsmith_fill:candidate_count(Bound, DBL, Idx, NB),
    crosswordsmith_fill:candidate_count(Masks, Bound, DBL, Idx, NB),
    NB =:= 2,
    % and a no-masks artifact-mode fill completes on the ordset kernel
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _Size, Slots, _),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    findall(A, ( member(W, Numbered), pw_dir(W, across), pw_answer(W, A) ), Across),
    Across == ['CAT', 'ORE', 'WED'].

% Deterministic bound patterns, lengths 3-6: single-bound cells over a letter
% spread (incl. rare Q/Z -> often empty intersections) and 0-anchored two-cell
% patterns. Covers the all-idx branch, dead cells (mask 0), and multi-cell chains.
sample_pattern(Len, [P-C]) :-
    member(Len, [3,4,5,6]), between(0, 5, P), P < Len,
    member(C, ['A','E','I','O','S','T','R','N','Q','Z']).
sample_pattern(Len, [0-C0, P1-C1]) :-
    member(Len, [4,5,6]), between(1, 5, P1), P1 < Len,
    member(C0, ['S','C','T','B']), member(C1, ['A','E','O']).

:- end_tests(fill).

% Helper: apply a list of frag/4 seeds directly (mirrors apply_seeds/3 minus the
% file read), for the no-slot test; the seeded-key accumulator is discarded.
apply_seeds_frags(Slots, Frags) :- foldl(crosswordsmith_fill:apply_seed(Slots), Frags, [], _).
