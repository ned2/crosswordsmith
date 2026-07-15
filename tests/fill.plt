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
    crosswordsmith_fill:fill_grid(GridFile, Size, Slots, _CellVar),
    ( SeedFile == none -> SeededKeys = [] ; crosswordsmith_fill:apply_seeds(SeedFile, Size, Slots, SeededKeys) ),
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
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_top.json', Size, Slots, SeededKeys),
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
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_top.json', Size, Slots, SeededKeys),
    exclude(crosswordsmith_fill:seeded_slot(SeededKeys), Slots, SearchSlots),
    crosswordsmith_fill:load_dict('fixtures/dict_cow.txt', DictByLen, Index),
    crosswordsmith_fill:fill_attempt(SearchSlots, Slots, DictByLen, Index, Outcome, _, _),
    Outcome == infeasible.

% Trigger B: two seeds pinning the SAME answer are rejected up front with the
% clean hooked error (neither slot is searched, so the Used-set fix cannot
% catch this shape; reported before searching, like fill_seed_no_slot/clash).
test(fill_duplicate_seed_answers_rejected,
     [throws(error(fill_seed_duplicate('COW'), _))]) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_split3.json', Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_both.json', Size, Slots, _).

% --- thin seed form (§6.6 convenience form at the fill boundary) --------------

% The AC-FRAG-4 analogue for fill: the thin spelling of the COW seed (desugared
% on the fill grid's own side, the SizeCtx apply_seeds passes) produces the
% identical fill as the canonical spelling. Asserts the expected answers too,
% so identical-but-wrong cannot pass.
test(fill_thin_seed_identical_to_canonical) :-
    do_fill('fixtures/fill_grid_split3.json', 'fixtures/fill_seed_cow_top.json',
            'fixtures/dict_cow_pig.txt', AC, DC),
    do_fill('fixtures/fill_grid_split3.json', 'fixtures/fill_seed_cow_top_thin.json',
            'fixtures/dict_cow_pig.txt', AT, DT),
    AC-DC == AT-DT,
    AC == ['COW', 'PIG'],
    DC == [].

% A thin seed is framed by THIS grid's side (3 here): COW from [0,1] across
% runs off it and is rejected at parse with the grid's own dimensions in the
% error. (The other thin parse errors - bad dir, bad row/col, bad top level -
% are consumer-independent and covered in arrange.plt.)
test(fill_thin_seed_off_grid_throws,
     [throws(error(fragment_thin_off_grid('COW', 0, 1, across, 3), _))]) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', Size, Slots, _),
    crosswordsmith_fill:apply_seeds('fixtures/fill_seed_cow_off_thin.json', Size, Slots, _).

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
    % the CURRENT schema version so the SWI gate is reached, not the version
    % gate (hardcoding it broke once, at the v2->v3 Unicode-hardening bump).
    crosswordsmith_fill:fill_index_format_version(V),
    setup_call_cleanup(open(AF, write, S, [type(binary)]),
                       fast_write(S, fill_index(V, [swi_version('0.0.0'), dict_sha256(x)], E, E)),
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


% --- Unicode/locale hardening (docs/plans/fill-dict-unicode-normalization.md)
% Policy: fold-then-strict. NFC and NFD forms must agree for EVERY class;
% pure-ASCII dictionaries must normalize exactly as before the hardening;
% dropped words warn unconditionally on stderr (INV-3), silent when zero.
% NB every non-ASCII char below is a \uXXXX escape, never a literal: the
% locale-matrix run consults this file under LC_ALL=C, where literal UTF-8
% source bytes would themselves mojibake - and escapes make the NFC-vs-NFD
% pairs visibly different in source.

% Write a temp dict with EXPLICIT utf8 encoding (tmp_dict/2 writes in the
% locale encoding - wrong vehicle for non-ASCII fixtures under LC_ALL=C).
tmp_dict_utf8(Text, Path) :-
    tmp_file_stream(binary, Path, S0), close(S0),
    setup_call_cleanup(open(Path, write, S, [encoding(utf8)]),
                       write(S, Text), close(S)).

% Capture user_error during Goal (the drop warning is unconditional stderr).
with_stderr_string(Goal, Str) :-
    stream_property(OldErr, alias(user_error)),
    with_output_to(string(Str),
        ( current_output(Cap),
          setup_call_cleanup(set_stream(Cap, alias(user_error)),
                             Goal,
                             set_stream(OldErr, alias(user_error))) )).

% The pre-hardening pipeline, verbatim (string_upper + char_type(alpha)
% squeeze): the ASCII-freeze witness below proves the new pipeline is
% output-identical to it on every shipped dictionary line.
old_line_word(Line, Letters) :-
    string_upper(Line, U), string_chars(U, Cs),
    old_alpha_chars(Cs, Letters),
    Letters \== [].
old_alpha_chars([], []).
old_alpha_chars([C|Cs], Letters) :-
    (   char_type(C, alpha)
    ->  Letters = [C|Rest]
    ;   Letters = Rest
    ),
    old_alpha_chars(Cs, Rest).

% cafe with precomposed e-acute (NFC) vs e + combining acute (NFD): same
% text, historically different letter lists by byte encoding.
test(norm_nfc_nfd_agree_cafe) :-
    crosswordsmith_fill:line_word("caf\u00e9", L1),
    crosswordsmith_fill:line_word("cafe\u0301", L2),
    L1 == ['C','A','F','E'],
    L2 == L1.

test(norm_sharp_s_folds_to_ss) :-
    crosswordsmith_fill:line_word("Stra\u00dfe", L),
    L == ['S','T','R','A','S','S','E'].

test(norm_multichar_fold_mid_word) :-
    crosswordsmith_fill:line_word("\u0132sselmeer", L),   % IJ ligature
    L == ['I','J','S','S','E','L','M','E','E','R'].

test(norm_typographic_apostrophe_squeezes) :-
    crosswordsmith_fill:line_word("don't", L1),
    crosswordsmith_fill:line_word("don\u2019t", L2),
    L1 == ['D','O','N','T'],
    L2 == L1.

% Unfoldable letters drop the WORD - in both encodings (the review's
% consistency requirement): Romanian s-comma is Latin Ext-B (not in the
% table) and its NFD mark U+0326 is not an allowed pair.
test(norm_unfoldable_drops_nfc, [fail]) :-
    crosswordsmith_fill:line_word("\u0219ir", _).
test(norm_unfoldable_drops_nfd, [fail]) :-
    crosswordsmith_fill:line_word("s\u0326ir", _).

% Pair-keyed marks: w-acute is a real letter (U+1E83) but NOT in the table,
% so BOTH its precomposed and decomposed forms drop - squeezing the mark just
% because "acute is a known mark" would keep NFD while dropping NFC.
test(norm_untabled_pair_drops_nfc, [fail]) :-
    crosswordsmith_fill:line_word("\u1e83ord", _).
test(norm_untabled_pair_drops_nfd, [fail]) :-
    crosswordsmith_fill:line_word("w\u0301ord", _).

% Mark chains never squeeze (e-acute + a second acute has no table home in
% either form; marks reset the previous-kept-char state).
test(norm_mark_chain_drops, [fail]) :-
    crosswordsmith_fill:line_word("caf\u00e9\u0301", _).

test(norm_cyrillic_word_drops, [fail]) :-
    crosswordsmith_fill:line_word("\u0434\u043e\u043c", _).   % "dom"

% Digit/punctuation squeeze is UNCHANGED (pre-hardening behavior).
test(norm_ascii_squeeze_unchanged) :-
    crosswordsmith_fill:line_word("3d", L),
    L == ['D'].

% A line with nothing representable is a silent SKIP (uncounted), while an
% unfoldable word is a counted DROP - load_dict distinguishes them: the mixed
% dict below has one skip (123), one drop (Cyrillic), two keeps.
test(load_dict_drop_count_and_skip,
     [cleanup(delete_file(DF))]) :-
    tmp_dict_utf8("CAT\n123\n\u0434\u043e\u043c\ncaf\u00e9\n", DF),
    with_stderr_string(crosswordsmith_fill:load_dict(DF, DBL, _), Err),
    assoc_to_list(DBL, Buckets),
    Buckets == [3-[['C','A','T']], 4-[['C','A','F','E']]],
    once(sub_string(Err, _, _, _, "dropped 1 dictionary word")),
    nb_getval(fill_dict_dropped_words, 1).

% Pure-ASCII dictionaries load in SILENCE (fill's quiet-by-default stderr
% contract rests on the zero-drop path printing nothing).
test(load_dict_ascii_silent) :-
    with_stderr_string(
        crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt', _, _),
        Err),
    Err == "".

% The encoding(utf8) pin: raw UTF-8 bytes (0xC3 0xA9 = e-acute) written with
% NO Prolog encoding layer load identically regardless of the process locale
% (without the pin, LC_ALL=C decodes these bytes as two mojibake chars and
% the word is dropped instead of folded).
test(load_dict_utf8_pin,
     [cleanup(delete_file(DF))]) :-
    tmp_file_stream(binary, DF, S),
    format(S, "caf~s~n", [[0xC3, 0xA9]]),
    close(S),
    crosswordsmith_fill:load_dict(DF, DBL, _),
    assoc_to_list(DBL, [4-[['C','A','F','E']]]).

% ASCII-freeze witness: on every line of the CLI default dictionary AND the
% full frozen ENABLE list, the new pipeline is output-identical to the
% pre-hardening pipeline (the hard no-output-change acceptance).
test(norm_ascii_freeze_shipped_dicts) :-
    forall(member(F, ['fixtures/wordlist_sample.txt', 'fixtures/dict/enable1.txt']),
           ( crosswordsmith_fill:read_file_lines(F, Lines),
             convlist(old_line_word, Lines, Old),
             convlist(crosswordsmith_fill:line_word, Lines, New),
             New == Old )).

% Oracle: the generated fold table and mark pairs agree with the Unicode
% database (library(unicode)/utf8proc - native-only, NEVER imported by
% fill.pl: the WASM build does not ship it). Guards against a hand-edited or
% regeneration-drifted table.
unicode_oracle_available :-
    catch(use_module(library(unicode)), _, fail).

test(fold_table_matches_unicode_oracle, [condition(unicode_oracle_available)]) :-
    forall(crosswordsmith_fill:fold_char(C, F), fold_oracle_ok(C, F)),
    forall(crosswordsmith_fill:allowed_mark(B, M), mark_oracle_ok(B, M)).

% The explicitly anchored specials (casefold / NFKC / orthographic
% convention): sharp-s, AE + OE ligatures, IJ ligature, O-slash.
special_fold('\u00df', ['S','S']).   % sharp s
special_fold('\u00c6', ['A','E']).   % AE ligature
special_fold('\u00e6', ['A','E']).
special_fold('\u0152', ['O','E']).   % OE ligature
special_fold('\u0153', ['O','E']).
special_fold('\u0132', ['I','J']).   % IJ ligature
special_fold('\u0133', ['I','J']).
special_fold('\u00d8', ['O']).       % O with stroke
special_fold('\u00f8', ['O']).

% Every non-special table entry must be exactly base + one combining mark
% under unicode_nfd/2, folding to the upcased base.
fold_oracle_ok(C, F) :-
    (   special_fold(C, SF)
    ->  F == SF
    ;   atom_chars(A, [C]),
        unicode_nfd(A, D),
        atom_codes(D, [Base, Mark]),
        Mark >= 0x0300, Mark =< 0x036F,
        (   Base >= 0'a, Base =< 0'z
        ->  Up is Base - 32
        ;   Up = Base
        ),
        char_code(UC, Up),
        F == [UC]
    ).

% Every allowed (base, mark) pair must recompose to a fold-table letter that
% folds back to exactly that base (the pair-keyed consistency invariant).
mark_oracle_ok(B, M) :-
    atom_chars(Decomposed, [B, M]),
    unicode_nfc(Decomposed, Composed),
    atom_chars(Composed, [P]),
    crosswordsmith_fill:fold_char(P, [B]).

% --- §8.4a scored ingestion + --min-score prune + ordering (AC-FILL-5/6/8) ---
% The scored fixture (fixtures/dict_scored_sample.txt) is ORIGINAL (INV-4):
% a score-0 blocklist entry (AAA), unscored lines (RAW/RED/ROW -> uniform 1),
% one malformed line (BAD;lots -> dropped + reported), and equal-score ties
% (CAT/COW at 95, DOG/WED at 80) - the mainline tiebreak case.

% Bucket order is the §8.4a total order: score-descending, then lex (ties
% collapse to dictionary order). Pins the whole 3-letter bucket, incl. the
% default `score >= 1` prune (DP-5): AAA (score 0) is OUT, unscored words
% (score 1) are IN and sort last.
test(scored_bucket_order_score_desc_then_lex) :-
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [], DBL, _, scores(_)),
    assoc_to_list(DBL, [3-Bucket]),
    findall(A, ( member(W, Bucket), atom_chars(A, W) ), Atoms),
    Atoms == ['CAT','COW','ORE','DOG','WED','ERA','DEW','ARE','GAD',
              'AWE','OWE','TED','CAR','COT','RAW','RED','ROW'].

% The score assoc keeps PRE-prune scores (a pruned word is still scorable by
% the report); unscored lines carry the uniform score 1 (DP-5, AC-FILL-8).
test(scored_assoc_uniform_and_blocklist_scores) :-
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [], _, _, scores(SA)),
    get_assoc(['R','E','D'], SA, 1),
    get_assoc(['A','A','A'], SA, 0),
    get_assoc(['C','A','T'], SA, 95).

% An explicit --min-score 0 re-admits the score-0 floor (the DP-5 default is
% a default, not a hard rule); AAA sorts last (lowest score).
test(scored_min_score_0_includes_blocklist) :-
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [min_score(0)], DBL, _, _),
    assoc_to_list(DBL, [3-Bucket]),
    last(Bucket, ['A','A','A']).

% The malformed scored line is dropped + counted + reported (INV-3,
% AC-FILL-8) - never squeezed into a word (the old plain reading would have
% kept BADLOTS).
test(scored_malformed_line_reported) :-
    with_stderr_string(
        crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                      [], DBL, _, _),
        Err),
    once(sub_string(Err, _, _, _, "dropped 1 malformed scored line")),
    nb_getval(fill_dict_malformed_lines, 1),
    assoc_to_list(DBL, [3-_]).             % nothing but the 3-letter bucket

% A duplicate word keeps its MAXIMUM score, deterministically and
% input-order-independently (both orders -> the same single entry).
test(scored_duplicate_word_keeps_max,
     [cleanup(( delete_file(D1), delete_file(D2) ))]) :-
    tmp_dict("FOO;30\nFOO;50\nFOO\n", D1),
    tmp_dict("FOO;50\nFOO\nFOO;30\n", D2),
    crosswordsmith_fill:load_dict(D1, [], DBL1, _, scores(SA1)),
    crosswordsmith_fill:load_dict(D2, [], DBL2, _, scores(SA2)),
    get_assoc(['F','O','O'], SA1, 50),
    get_assoc(['F','O','O'], SA2, 50),
    assoc_to_list(DBL1, [3-[['F','O','O']]]),
    DBL1 == DBL2.

% The hard prune (AC-FILL-5): --min-score 50 removes every word scoring < 50
% from the buckets (so no slot domain can see one), and the prune is
% reported unconditionally on stderr (INV-3).
test(scored_min_score_prunes_and_reports) :-
    with_stderr_string(
        crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                      [min_score(50)], DBL, _, scores(SA)),
        Err),
    once(sub_string(Err, _, _, _, "--min-score 50 pruned 9 of 18")),
    assoc_to_list(DBL, [3-Bucket]),
    length(Bucket, 9),
    forall(member(W, Bucket), ( get_assoc(W, SA, S), S >= 50 )).

% A prune that empties the dictionary reports the targeted hint (max score;
% the `--min-score 2` vs uniform-score-1 degenerate case), and the empty
% domains downstream are the ordinary AC-FILL-1 infeasible outcome.
test(scored_prune_all_hints_max_score) :-
    with_stderr_string(
        crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt',
                                      [min_score(2)], DBL, Idx, _),
        Err),
    once(sub_string(Err, _, _, _, "exceeds the dictionary's maximum score 1")),
    once(sub_string(Err, _, _, _, "unscored words score 1")),
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, Outcome, _, _),
    Outcome == infeasible.

% AC-FILL-6, structural half: a plain wordlist loads IDENTICALLY through the
% options path (Scores == uniform, same buckets, same index), so the
% pre-scoring engine's behavior is reproduced by construction.
test(scored_plain_path_identity) :-
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt', DBL3, Idx3),
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt',
                                  [], DBL5, Idx5, Scores),
    Scores == uniform,
    DBL3 == DBL5,
    Idx3 == Idx5.

% AC-FILL-6, uniform-score half: the same words all at one explicit score
% produce the exact plain-path buckets + index (score-desc collapses to lex).
test(scored_uniform_score_collapses_to_lex,
     [cleanup(( delete_file(DP), delete_file(DS) ))]) :-
    tmp_dict("CAT\nCOW\nARE\nTED\nORE\nWED\n", DP),
    tmp_dict("CAT;7\nCOW;7\nARE;7\nTED;7\nORE;7\nWED;7\n", DS),
    crosswordsmith_fill:load_dict(DP, [], DBLP, IdxP, uniform),
    crosswordsmith_fill:load_dict(DS, [], DBLS, IdxS, scores(_)),
    DBLP == DBLS,
    IdxP == IdxS.

% Score-descending ordering drives the SEARCH: with two disjoint word
% squares available, the fill lands on the high-scored one even though the
% low-scored square wins lexicographically (CAT < DOG).
test(scored_ordering_prefers_high_score_square,
     [cleanup(delete_file(DF))]) :-
    tmp_dict("CAT;10\nCOW;10\nARE;10\nTED;10\nORE;90\nWED;90\nDOG;90\nERA;90\nDEW;90\nGAD;90\n", DF),
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict(DF, [], DBL, Idx, _),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    findall(A, ( member(PW, Numbered), pw_answer(PW, A) ), As),
    msort(As, Sorted),
    Sorted == ['DEW','DOG','ERA','GAD','ORE','WED'].

% AC-FILL-5 end to end on the fixture: a --min-score 50 fill contains no
% word scoring < 50 (and it exists - the DOG/ERA/WED square is all >= 50).
test(scored_min50_fill_all_clean) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [min_score(50)], DBL, Idx, scores(SA)),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    forall(( member(PW, Numbered), pw_answer(PW, A) ),
           ( atom_chars(A, W), get_assoc(W, SA, S), S >= 50 )).

% --- the fill-quality sidecar report (AC-FILL-7) ------------------------------
% Stats over the default fixture fill (CAT/COW/ARE/TED/ORE/WED at
% 95/95/60/30/90/80): n=6, mean 75.0, min 30 (TED), one entry below the 50
% clean floor. The same numbers score_fill.py computes post hoc.
test(quality_stats_scored_default) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [], DBL, Idx, Scores),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    crosswordsmith_fill:fill_quality_stats(Numbered, Scores, DBL, St),
    St == stats(6, 75.0, 30, 1).

% A uniform (plain) dict scores 1 per in-dict word: degenerate but defined
% (n=6, mean 1.0, min 1, all 6 below the 50 floor).
test(quality_stats_uniform_dict) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt',
                                  [], DBL, Idx, Scores),
    Scores == uniform,
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    crosswordsmith_fill:fill_quality_stats(Numbered, Scores, DBL, St),
    St == stats(6, 1.0, 1, 6).

% An answer absent from the dictionary (a non-dict seed pin) scores 0 -
% score_fill.py's junk convention - on both the uniform and the scored path.
test(quality_absent_word_scores_zero) :-
    crosswordsmith_fill:load_dict('fixtures/wordlist_sample.txt',
                                  [], DBLU, _, uniform),
    crosswordsmith_fill:placed_word_score(uniform, DBLU,
                                          pw('ZZZ', _, _, _, _, _, _, _), SU),
    SU == 0,
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt',
                                  [], DBLS, _, Scores),
    crosswordsmith_fill:placed_word_score(Scores, DBLS,
                                          pw('ZZZ', _, _, _, _, _, _, _), SS),
    SS == 0.

% The report JSON is one sorted-key object, byte-pinned (AC-X-2 spirit; the
% CLI-level file lands in the goldens).
test(quality_report_json_bytes, [cleanup(delete_file(RF))]) :-
    tmp_path(RF),
    crosswordsmith_fill:write_quality_json(RF, 6, 75.0, 30, 1, 50),
    read_file_to_string(RF, S, []),
    S == "{\"belowThreshold\":1,\"mean\":75.0,\"min\":30,\"n\":6,\"threshold\":50}\n".

% --- §8.4b search levers (DP-6): --budget (AC-FILL-9) -------------------------

% The -1 sentinel (and an absent option) selects the engine default;
% fill_budget/1 stays the single source of truth.
test(budget_option_default_and_sentinel) :-
    crosswordsmith_fill:fill_effective_budget([], B0),
    crosswordsmith_fill:fill_budget(Default),
    B0 == Default,
    crosswordsmith_fill:fill_effective_budget([budget(-1)], B1),
    B1 == Default.

test(budget_option_override) :-
    crosswordsmith_fill:fill_effective_budget([budget(12345)], B),
    B == 12345.

% White-box callers get the same clean rejection the CLI gives (AC-FILL-9).
test(budget_option_rejects_zero, [throws(error(type_error(positive_integer, 0), _))]) :-
    crosswordsmith_fill:fill_effective_budget([budget(0)], _).

test(budget_option_rejects_nonint, [throws(error(type_error(positive_integer, lots), _))]) :-
    crosswordsmith_fill:fill_effective_budget([budget(lots)], _).

% A too-small budget yields the ordinary AC-FILL-1 not-proven failure (no
% stdout), through the full fill_solve/5 option path. (Budget 10: the 3x3
% fixture fill completes within a few hundred inferences, so a four-digit
% "tiny" budget is not tiny enough - found the hard way.)
test(budget_tiny_fill_solve_fails) :-
    with_output_to(string(S),
        \+ crosswordsmith_fill:fill_solve('fixtures/fill_grid_3.json', none,
                                          'fixtures/dict_scored_sample.txt',
                                          fixed, [budget(10)])),
    S == "".

% AC-FILL-9's core clause: a budget change never alters the CONTENT of a
% produced fill - default and a vastly larger budget emit byte-identical
% layouts.
test(budget_never_changes_content) :-
    with_output_to(string(S1),
        crosswordsmith_fill:fill_solve('fixtures/fill_grid_3.json', none,
                                       'fixtures/dict_scored_sample.txt',
                                       fixed, [])),
    with_output_to(string(S2),
        crosswordsmith_fill:fill_solve('fixtures/fill_grid_3.json', none,
                                       'fixtures/dict_scored_sample.txt',
                                       fixed, [budget(100_000_000_000)])),
    S1 == S2,
    S1 \== "".

% DP-5's default in action: the score-0 blocklist entry never appears in a
% default-flags fill.
test(scored_default_fill_excludes_score0) :-
    crosswordsmith_fill:fill_grid('fixtures/fill_grid_3.json', _, Slots, _),
    crosswordsmith_fill:load_dict('fixtures/dict_scored_sample.txt', DBL, Idx),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DBL, Idx, filled, Numbered, _),
    findall(A, ( member(PW, Numbered), pw_answer(PW, A) ), As),
    \+ memberchk('AAA', As),
    msort(As, Sorted),
    Sorted == ['ARE','CAT','COW','ORE','TED','WED'].

:- end_tests(fill).

% Helper: apply a list of frag/4 seeds directly (mirrors apply_seeds/4 minus the
% file read), for the no-slot test; the seeded-key accumulator is discarded.
apply_seeds_frags(Slots, Frags) :- foldl(crosswordsmith_fill:apply_seed(Slots), Frags, [], _).
