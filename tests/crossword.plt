% tests/crossword.plt - plunit test suite for the shared substrate (core.pl, formerly crossword.pl)
%
% These tests assume core.pl has already been consulted into the
% `user` module (the runner, tests/run_tests.pl, does this before loading
% this file). Run them via:
%
%     ./run_tests.sh        (or)        make test
%
% The tests are grouped into units: pure grid geometry, generic list
% utilities, the solver, and clue numbering.

:- use_module(library(plunit)).
:- use_module(library(json)).
% lists/apply/assoc/aggregate: test-body helpers; explicit so the suite also
% runs under autoload(false) (P11/C5).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(assoc)).
:- use_module(library(aggregate)).
:- use_module(library(solution_sequences)).   % limit/2 in the count_upto2 unit

bundled_words(Words) :-
    load_clues('fixtures/bundled_17_clues.pl', Words).


% Grid geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cells are numbered 1..GridLen*GridLen in row-major order. `across` steps
% by +1 along a row; `down` steps by +GridLen down a column.

:- begin_tests(geometry).

test(swap_dir_across) :- crosswordsmith_core:swap_dir(across, down).
test(swap_dir_down)   :- crosswordsmith_core:swap_dir(down, across).

% A word fits across only if there is room to its right in the row.
test(fits_across_full_row)      :- fits_on_grid(across, 1, 17, 17).
test(fits_across_too_long, [fail]) :- fits_on_grid(across, 1, 18, 17).
% Cell 10 is column 10 of 17, so only 8 cells remain (10..17).
test(fits_across_mid_row)       :- fits_on_grid(across, 10, 8, 17).
test(fits_across_overflow, [fail]) :- fits_on_grid(across, 10, 9, 17).
% A cell at the end of a row (mod == 0) cannot start an across word.
test(fits_across_row_end, [fail]) :- fits_on_grid(across, 17, 2, 17).

% A word fits down only if it does not run off the bottom of the grid.
test(fits_down_full_col)        :- fits_on_grid(down, 1, 17, 17).
test(fits_down_too_long, [fail]) :- fits_on_grid(down, 1, 18, 17).

% calc_num and calc_start are inverses: given a word starting at WStart,
% the letter at position P has cell number N; calc_start must recover
% WStart from N and P.
test(calc_num_across, [true(N =:= 7)])  :- crosswordsmith_core:calc_num(across, 17, 3, 5, N).
test(calc_num_down,   [true(N =:= 35)]) :- crosswordsmith_core:calc_num(down, 17, 3, 1, N).
test(calc_roundtrip_across, [true(S =:= 5)]) :-
    crosswordsmith_core:calc_num(across, 17, 3, 5, N), crosswordsmith_core:calc_start(across, 17, 3, N, S).
test(calc_roundtrip_down, [true(S =:= 1)]) :-
    crosswordsmith_core:calc_num(down, 17, 4, 1, N), crosswordsmith_core:calc_start(down, 17, 4, N, S).

test(next_cell_across, [true(C =:= 6)])  :- next_cell(across, 5, 17, C).
test(next_cell_down,   [true(C =:= 22)]) :- next_cell(down, 5, 17, C).
test(prev_cell_across, [true(C =:= 4)])  :- crosswordsmith_core:prev_cell(across, 5, 17, C).
test(prev_cell_down,   [true(C =:= 5)])  :- crosswordsmith_core:prev_cell(down, 22, 17, C).

% First column starts across rows; first row starts down columns.
test(is_start_across)            :- crosswordsmith_core:is_start_cell(across, 18, 17).   % col 1, row 2
test(is_start_across_no, [fail]) :- crosswordsmith_core:is_start_cell(across, 19, 17).
test(is_start_down)              :- crosswordsmith_core:is_start_cell(down, 5, 17).      % row 1
test(is_start_down_no, [fail])   :- crosswordsmith_core:is_start_cell(down, 18, 17).
test(is_end_across)              :- crosswordsmith_core:is_end_cell(across, 17, 17).     % last col
test(is_end_down)                :- crosswordsmith_core:is_end_cell(down, 289, 17).      % last row
% P1 regression (crosswordsmith_core:is_end_cell(down) off-by-one): on a 5x5 the bottom row is
% cells 21..25 and cell 20 is the LAST cell of row 4 - cell 25 lies directly
% below it - so cell 20 is NOT a down end cell. The old `>=` (Num >= 20)
% misclassified it as one; the fix is `>` (Num > 20).
test(is_end_down_penultimate_row_not_end, [fail]) :- crosswordsmith_core:is_end_cell(down, 20, 5).
test(is_end_down_bottom_row_start) :- crosswordsmith_core:is_end_cell(down, 21, 5).
test(is_end_down_bottom_row_end)   :- crosswordsmith_core:is_end_cell(down, 25, 5).

% The four named start locations resolve to the expected cell + direction.
test(start_topleft_across, [true(N-D == 1-across)])  :- start_loc(topleft_across, 17, N, D).
test(start_topleft_down,   [true(N-D == 1-down)])    :- start_loc(topleft_down, 17, N, D).
test(start_topright,       [true(N-D == 17-down)])   :- start_loc(topright, 17, N, D).
test(start_bottomleft,     [true(N-D == 273-across)]):- start_loc(bottomleft, 17, N, D).

:- end_tests(geometry).


% Generic list utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(utilities).

% remove_x removes only the first occurrence. (nondet: it leaves a harmless
% choicepoint; we only care about the first answer.)
test(remove_x_first, [true(R == [a,c,b]), nondet]) :- remove_x(b, [a,b,c,b], R).
test(remove_x_absent, [true(R == [a,b,c]), nondet]) :- remove_x(z, [a,b,c], R).

% init_grid builds a grid(...) term of GridLen*GridLen unbound (empty) cells;
% cell Num is arg Num.
test(init_grid_size, [true(Arity =:= 9)]) :-
    init_grid(3, G), functor(G, grid, Arity).
test(init_grid_empty) :-
    init_grid(3, G), arg(1, G, C1), var(C1), arg(9, G, C9), var(C9).
test(init_grid_no_extra, [fail]) :-
    init_grid(3, G), arg(10, G, _).

% init_gs builds the gs(LetterGrid, BoundaryGrid) bundle the search threads:
% two independent all-var grids of the same shape.
test(init_gs_shape) :-
    init_gs(3, gs(L, B)),
    functor(L, grid, 9), functor(B, grid, 9),
    arg(5, L, C), var(C), arg(5, B, M), var(M).

% answer_meta_assoc/2 (P5): the answer->meta assoc that replaces the O(n^2)
% member/2 rescan in the emit join. An entry with metadata maps to its dict; an
% entry with none ([Answer]) maps to an empty dict.
test(answer_meta_assoc_join) :-
    Words = [['CAT', _{clue:"feline"}], ['DOG']],
    answer_meta_assoc(Words, A),
    get_assoc('CAT', A, M1), get_dict(clue, M1, "feline"),
    get_assoc('DOG', A, M2), dict_pairs(M2, _, []).

:- end_tests(utilities).


% Saturating solution counter (mrv_count's hot Cap=2 path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count_upto2/2 must reproduce aggregate_all(count, limit(2, Goal), _)
% EXACTLY: 0 / 1 / 2-saturated, and - like findall - leave NO residual
% bindings from Goal on exit (the arrange search relies on the counted
% Start/Dir staying untouched).

:- begin_tests(count_upto2).

test(zero, [true(N =:= 0)]) :-
    crosswordsmith_core:count_upto2(member(_, []), N).
test(one, [true(N =:= 1)]) :-
    crosswordsmith_core:count_upto2(member(_, [a]), N).
test(two_exact, [true(N =:= 2)]) :-
    crosswordsmith_core:count_upto2(member(_, [a,b]), N).
test(saturates_at_two, [true(N =:= 2)]) :-
    crosswordsmith_core:count_upto2(member(_, [a,b,c,d,e]), N).
test(filtered_goal_one, [true(N =:= 1)]) :-
    crosswordsmith_core:count_upto2((member(X, [1,2,3]), X =:= 2), N).
test(filtered_goal_zero, [true(N =:= 0)]) :-
    crosswordsmith_core:count_upto2((member(X, [1,2,3]), X > 9), N).

% No residual bindings: the counted variable is unbound after counting
% (exercises both the >=2 early-exit path and the exhausted path).
test(no_residual_binding_saturated, [true(var(X))]) :-
    crosswordsmith_core:count_upto2(member(X, [a,b,c]), _).
test(no_residual_binding_single, [true(var(X))]) :-
    crosswordsmith_core:count_upto2(member(X, [only]), _).
test(no_residual_binding_empty, [true(var(X))]) :-
    crosswordsmith_core:count_upto2(member(X, []), _).

% Semidet: yields exactly one solution, no choicepoint left behind.
test(deterministic) :-
    findall(N, crosswordsmith_core:count_upto2(member(_, [a,b,c]), N), Ns),
    Ns == [2].

% Matches the reference (aggregate_all + limit/2) across a range of arities.
test(matches_reference, [forall(member(L, [[], [x], [x,y], [x,y,z], [x,y,z,w]]))]) :-
    aggregate_all(count, limit(2, member(_, L)), Ref),
    crosswordsmith_core:count_upto2(member(_, L), Got),
    Got =:= Ref.

:- end_tests(count_upto2).


% The solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(solver).

% Two words that cross on a shared letter are both placed. OMEGA POINT
% (across, starting cell 1) and GNOSTIC GOSPELS (down) cross on G.
% (nondet: find_crossword backtracks over many solutions; we assert one exists.)
test(places_two_crossing_words, [nondet]) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    crosswordsmith_core:find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 2),
    % the first word is laid across starting at cell 1
    member(PW, Placed), pw_dir(PW, across), pw_start(PW, 1).

% The full bundled clue set has a solution on a 17x17 grid.
test(bundled_clues_solve_at_17, [nondet]) :-
    bundled_words(Words),
    crosswordsmith_core:find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 6).

% A grid that is too small for the words has no solution.
test(too_small_grid_fails, [fail]) :-
    bundled_words(Words),
    crosswordsmith_core:find_crossword(3, Words, topleft_across, _Grid, _Placed).

% End-to-end: the top-level crossword/3 emits one JSON object describing the
% solution. Parse it and assert the expected shape. Output is captured both to
% keep logs clean and to inspect it.
test(full_pipeline_emits_json) :-
    bundled_words(Words),
    with_output_to(string(S), crosswordsmith_core:crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Grid), length(Grid, 17),
    forall(member(Row, Grid), length(Row, 17)),
    get_dict(words, Dict, WordObjs), length(WordObjs, 6).

% A crossing cell carries BOTH directions, so both links remain reachable (G2).
% Row 0, col 3 is where OMEGA POINT (across 1) and GNOSTIC GOSPELS (down 2) cross.
test(crossing_cell_has_both_directions) :-
    bundled_words(Words),
    with_output_to(string(S), crosswordsmith_core:crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(grid, Dict, Grid),
    nth0(0, Grid, Row0),
    nth0(3, Row0, Cell),
    get_dict(across, Cell, 1),
    get_dict(down, Cell, 2).

% Per-word metadata is rejoined verbatim under `meta`.
test(word_metadata_under_meta) :-
    bundled_words(Words),
    with_output_to(string(S), crosswordsmith_core:crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(words, Dict, WordObjs),
    once(( member(W, WordObjs), get_dict(answer, W, "OMEGA POINT") )),
    get_dict(meta, W, Meta),
    get_dict(clue, Meta, "Transcending entropy"),
    get_dict(link, Meta, "http://en.wikipedia.org/wiki/Omega_Point").

% The emit-time join requires unique answers; a duplicate is rejected up front.
test(duplicate_answer_rejected, [throws(error(duplicate_answer('CAT'), _))]) :-
    crosswordsmith_core:crossword(5, [['CAT', _{}], ['CAT', _{}]], topleft_across).

% Regression for I5: a word must stay a MAXIMAL run. Placing DFEEE across over
% an already-placed DFEE (same start cell) would make DFEE a "word inside a
% word" - two same-direction answers at one start, which assign_clue_numbers/2
% cannot number. assign_word/9 must reject it (no_word_merge_bg/2 via the
% boundary grid marks maintained by mark_boundary/5).
test(rejects_collinear_prefix_overlap, [fail]) :-
    init_gs(10, G0),
    assign_word('DFEE',  [d,f,e,e],   4, 1, across, 10, G0, _P1, G1),
    assign_word('DFEEE', [d,f,e,e,e], 5, 1, across, 10, G1, _, _).

% ...but a collinear word separated by the required empty boundary cell is fine.
test(allows_separated_collinear_word, [nondet]) :-
    init_gs(10, G0),
    assign_word('DFEE', [d,f,e,e], 4, 1, across, 10, G0, _P1, G1),
    assign_word('DOG',  [d,o,g],   3, 6, across, 10, G1, _, _).

% The boundary grid after a placement: exactly the word's on-grid before-start /
% after-end cells are marked, edge-abutting sides are omitted, and the LETTER
% grid's boundary cells stay unbound (they must still read as empty to the
% adjacency checks - only BoundaryGrid carries the mark).
test(boundary_marks_after_placement, [nondet]) :-
    init_gs(10, G0), G0 = gs(L, B),
    % starts at column 1 (grid edge: no before-start mark), ends at cell 4
    assign_word('DFEE', [d,f,e,e], 4, 1, across, 10, G0, _P1, _G1),
    arg(5, B, M5), nonvar(M5),          % after-end boundary marked
    arg(5, L, C5), var(C5),             % ...but still EMPTY in the letter grid
    forall(between(1, 100, N), (N =:= 5 -> true ; arg(N, B, M), var(M))).

% check_word_fits/5 (the counting path's non-binding legality probe) must agree
% with assign_word/9 on accept AND reject, and must leave the grid untouched.
test(check_word_fits_agrees_with_assign_word, [nondet]) :-
    init_gs(10, GS), GS = gs(L, B),
    assign_word('DFEE', [d,f,e,e], 4, 1, across, 10, GS, _P1, _),
    % DOG down crossing the D at cell 1: legal for both
    crosswordsmith_core:check_word_fits([d,o,g], 1, down, 10, GS),
    % DFEEE across at 1 fills DFEE's boundary cell 5: illegal for both
    \+ crosswordsmith_core:check_word_fits([d,f,e,e,e], 1, across, 10, GS),
    % ...and the probe bound nothing: cells 11/21 (DOG's tail) still empty,
    % boundary grid unchanged beyond DFEE's own mark at 5
    arg(11, L, C11), var(C11), arg(21, L, C21), var(C21),
    arg(6, B, M6), var(M6).

% Regression for P1: a down word ending at cell (L-1)*L must still honour the
% "cell directly below must be empty" guard (check_next_cell/4). On a 5x5, a
% down word ending at cell 20 sits directly above cell 25; if cell 25 already
% holds a letter, placing that down word would splice a 4th cell onto its run -
% the forbidden collinear "word inside a word" merge (core.pl, no_word_merge).
% The old crosswordsmith_core:is_end_cell(down) `>=` classified cell 20 as a bottom-row end and so
% SKIPPED the below-must-be-empty check, wrongly allowing the placement.
test(rejects_down_word_ending_above_filled_cell, [fail]) :-
    init_gs(5, G0),
    % fill the whole bottom row across, so cell 25 (below cell 20) is occupied
    assign_word('VWXYZ', [v,w,x,y,z], 5, 21, across, 5, G0, _PA, G1),
    % a down word ending exactly at cell 20 must be rejected: cell 25 is filled
    assign_word('PQR', [p,q,r], 3, 10, down, 5, G1, _, _).

% all_crossword/5 (the --all/count path) totals every solution for a start
% position. It must agree with independently enumerating find_crossword's
% solutions, and do so deterministically (aggregate_all runs the search once).
test(all_crossword_counts_all_solutions, [true(Num =:= Expected)]) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    findall(G, crosswordsmith_core:find_crossword(baseline, 17, Words, topleft_across, G, _), Sols),
    length(Sols, Expected),
    all_crossword(baseline, 17, Words, topleft_across, Num).

% Exactly one answer, no leftover choicepoint (guards against a relapse to the
% length/findall idiom, which re-ran the search per candidate length).
test(all_crossword_is_deterministic) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    findall(N, all_crossword(baseline, 17, Words, topleft_across, N), Counts),
    Counts = [_].

:- end_tests(solver).


% Clue numbering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A placed word is the record pw(Answer,Letters,Cells,Dir,Len,Start,End,Num).
% assign_clue_numbers/2 sorts by start and fills in num.

:- begin_tests(clue_numbering).

% Helper: clue number assigned to the word with the given text.
clue_num_of(Word, Placed, Num) :-
    member(PW, Placed), pw_answer(PW, Word), pw_num(PW, Num).

% Two words with distinct start cells get sequential numbers, in start order.
test(distinct_starts_numbered_in_order, [nondet]) :-
    W1 = pw('CAT', [c,a,t], [1,2,3],   across, 3, 1, 3,  _),
    W2 = pw('DOG', [d,o,g], [5,12,19], down,   3, 5, 19, _),
    assign_clue_numbers([W2,W1], Placed),   % deliberately unsorted input
    clue_num_of('CAT', Placed, 1),
    clue_num_of('DOG', Placed, 2).

% Regression test for the (now-fixed) add_clue_nums/3 arity bug: when one
% cell is the start of BOTH an across and a down word, the two clues must
% share a clue number. Before the fix the buggy clause had the wrong arity
% and this scenario failed clue numbering entirely.
test(shared_start_cell_shares_number, [nondet]) :-
    Wa = pw('CAT', [c,a,t], [1,2,3],   across, 3, 1, 3,  _),
    Wd = pw('COW', [c,o,w], [1,18,35], down,   3, 1, 35, _),
    assign_clue_numbers([Wa,Wd], Placed),
    clue_num_of('CAT', Placed, N),
    clue_num_of('COW', Placed, N),
    integer(N).

:- end_tests(clue_numbering).


% JSON clue input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JSON input reads into the same internal Words list as Prolog fixtures, so
% the whole pipeline downstream is unchanged.
% Paths are relative to the repo root (the runner consults from there).

:- begin_tests(json_input).

% Reduce one emitted words[] entry back to an input clue {answer, meta}.
reduce_word(W, _{answer:A, meta:M}) :-
    get_dict(answer, W, A),
    get_dict(meta, W, M).

% The fixture loads to the internal Words form: atom answers, dict metadata.
test(json_loads_to_words) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    length(Words, 6),
    once((member([Answer, Meta], Words), Answer == 'OMEGA POINT')),
    atom(Answer),
    % json_read_dict yields string meta *values* (vs the Prolog path's atoms);
    % json_write_dict renders both to the same JSON token (input spec D3).
    get_dict(clue, Meta, "Transcending entropy"),
    get_dict(link, Meta, "http://en.wikipedia.org/wiki/Omega_Point").

% An omitted `meta` key and an explicit empty `{}` both yield an empty dict.
test(json_absent_meta_is_empty_dict) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    once((member([A1, M1], Words), A1 == 'GNOSTIC GOSPELS')),  % no meta key
    once((member([A2, M2], Words), A2 == 'BIAS')),             % meta: {}
    dict_pairs(M1, _, []),
    dict_pairs(M2, _, []).

% The JSON-loaded words drive the existing pipeline to a valid JSON solution.
test(json_pipeline_emits_json) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    with_output_to(string(S), crosswordsmith_core:crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Grid), length(Grid, 17),
    forall(member(Row, Grid), length(Row, 17)),
    get_dict(words, Dict, WordObjs), length(WordObjs, 6).

% Schema violations throw error/2 terms (rendered by error_message//1).
test(json_no_clues_array_throws, [throws(error(json_no_clues_array, _))]) :-
    crosswordsmith_core:doc_to_words(_{version: 1}, _).            % no `clues` key

test(json_clues_not_list_throws, [throws(error(json_no_clues_array, _))]) :-
    crosswordsmith_core:doc_to_words(_{clues: "nope"}, _).         % `clues` not an array

test(json_missing_answer_throws, [throws(error(json_invalid_answer(_), _))]) :-
    crosswordsmith_core:doc_to_words(_{clues: [_{meta: _{}}]}, _). % entry with no `answer`

test(json_non_string_answer_throws, [throws(error(json_invalid_answer(_), _))]) :-
    crosswordsmith_core:doc_to_words(_{clues: [_{answer: 42}]}, _).% numeric answer (would coerce)

test(json_non_object_meta_throws, [throws(error(json_invalid_meta('FLOW'), _))]) :-
    crosswordsmith_core:doc_to_words(_{clues: [_{answer: "FLOW", meta: "no"}]}, _).

% Malformed JSON and a missing file surface as standard ISO errors.
test(json_malformed_throws, [throws(_)]) :-
    setup_call_cleanup(open_string("{ not json", S),
                       json_read_dict(S, _),
                       close(S)).

test(json_missing_file_throws, [throws(error(existence_error(source_sink, _), _))]) :-
    crosswordsmith_core:read_clues_json('tests/no_such_file.json', _).

% Symmetry (G5): reduce the emitted output to {answer, meta} entries and
% confirm they round-trip as valid input — they load and re-solve, the answer
% set is preserved, and metadata survives. (Asserted as validity, not a byte
% diff: solver layout depends on input order, which the reduction reshuffles.)
test(json_output_reduces_to_valid_input) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    with_output_to(string(S1), crosswordsmith_core:crossword(17, Words, topleft_across)),
    atom_json_dict(S1, Dict, []),
    get_dict(words, Dict, WordObjs),
    maplist(reduce_word, WordObjs, ReducedClues),
    crosswordsmith_core:doc_to_words(_{clues: ReducedClues}, Words2),
    % the reduced output is itself valid input: it re-solves and re-emits JSON
    with_output_to(string(S2), crosswordsmith_core:crossword(17, Words2, topleft_across)),
    atom_json_dict(S2, Dict2, []),
    get_dict(words, Dict2, WordObjs2), length(WordObjs2, 6),
    % the answer set survives the round trip (atoms on both sides)
    findall(A, member([A, _], Words),  AsOrig0),  msort(AsOrig0,  AsOrig),
    findall(A, member([A, _], Words2), AsRound0), msort(AsRound0, AsRound),
    AsOrig == AsRound,
    % and metadata survives for a representative entry
    once((member([RA, RM], Words2), RA == 'OMEGA POINT')),
    get_dict(clue, RM, "Transcending entropy"),
    get_dict(link, RM, "http://en.wikipedia.org/wiki/Omega_Point").

:- end_tests(json_input).


% CLI option plumbing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main/0 parses argv with library(optparse) and dispatches via run/2; here we
% test the supporting predicates directly (load_clues/2, with_output/2,
% valid_loc/1). End-to-end argv handling is covered by the golden regression.

:- begin_tests(cli).

% load_clues/2 reads either Prolog fixtures or JSON input files.
test(load_clues_from_prolog_fixture) :-
    load_clues('fixtures/bundled_17_clues.pl', W),
    length(W, 6),
    W = [[A|_]|_], A == 'OMEGA POINT'.

test(load_clues_from_json_file) :-
    load_clues('tests/clues.json', W),
    length(W, 6),
    once((member([Ans, _], W), Ans == 'BIAS')).

test(load_clues_rejects_unsupported_extension,
     [throws(error(unsupported_clue_file('tests/clues.txt', txt), _))]) :-
    load_clues('tests/clues.txt', _).

% valid_loc/1 accepts the four named start locations and nothing else.
test(valid_loc_accepts_known) :- crosswordsmith_core:valid_loc(topleft_across).
test(valid_loc_rejects_unknown, [fail]) :- crosswordsmith_core:valid_loc(nowhere).

% with_output('', Goal) writes straight to the current output (stdout path).
test(with_output_stdout_passthrough) :-
    with_output_to(string(S), with_output('', write(hello))),
    S == "hello".

% with_output(File, Goal) sends Goal's output to File, byte-for-byte the same
% as the stdout path would produce.
test(with_output_file_matches_stdout) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    with_output_to(string(StdoutText), crosswordsmith_core:crossword(17, Words, topleft_across)),
    tmp_file_stream(text, Tmp, S0), close(S0),
    with_output(Tmp, crosswordsmith_core:crossword(17, Words, topleft_across)),
    setup_call_cleanup(open(Tmp, read, S), read_string(S, _, FileText), close(S)),
    delete_file(Tmp),
    FileText == StdoutText,
    atom_json_dict(FileText, Dict, []),         % and it is valid JSON
    get_dict(words, Dict, Ws), length(Ws, 6).

% On a goal that fails (grid too small for a layout), no file is written, so a
% failed run never leaves an empty output file behind.
test(with_output_no_file_on_failure) :-
    crosswordsmith_core:read_clues_json('tests/clues.json', Words),
    tmp_file_stream(text, Tmp, S0), close(S0),
    delete_file(Tmp),                           % ensure it is absent
    \+ with_output(Tmp, crosswordsmith_core:crossword(3, Words, topleft_across)),
    \+ exists_file(Tmp).

:- end_tests(cli).


% Shared layout metrics (metrics.pl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metrics.pl is loaded alongside core.pl (via its ensure_loaded). It holds the
% shared metric predicates (consumed by arrange.pl as optimizer signals and by
% lint as validators); the greedy density constructor moved to arrange.pl in
% Phase 3 of the source-structure migration, so its tests (seed selection) now
% live in arrange.plt. These tests pin the metrics.

:- begin_tests(quality).

% --- fixture: a hand-built layout with cell_rc-consistent geometry -----------
% A placed word is a pw/8 record; the metric predicates read only its answer, dir
% and cells fields, so the rest are left unbound. Cells are numbered as in the
% live engine.

% An across word crossed by two down words at its END cells (17-grid). Across A
% spans row 0, cols 0..2; each down spans a column over rows 0..2, sharing one
% cell with A. So A is checked at 2 of 3 cells; each down at only 1 of 3.
sample_cross([Wa, Wd1, Wd2]) :-
    Wa  = pw('ABC', _, [1,2,3],    across, _, _, _, _),
    Wd1 = pw('ADE', _, [1,18,35],  down,   _, _, _, _),
    Wd2 = pw('CFG', _, [3,20,37],  down,   _, _, _, _).

% --- scoring (checked_cells / dir_cells) -------------------------------------

% dir_cells/3 collects, per direction, the sorted ordset of covered cells.
test(dir_cells_partition_by_direction) :-
    sample_cross(Placed),
    crosswordsmith_metrics:dir_cells(Placed, across, AcrossCells), AcrossCells == [1,2,3],
    crosswordsmith_metrics:dir_cells(Placed, down,   DownCells),   DownCells   == [1,3,18,20,35,37].

% checked_cells/2 = cells in BOTH an across and a down word (the two crossings).
test(checked_cells_counts_crossings, [true(N =:= 2)]) :-
    sample_cross(Placed),
    crosswordsmith_metrics:checked_cells(Placed, N).

% --- per-word metrics (the lint-rule primitives) -----------------------------

% word_checked_count/3 = how many of a word's cells a perpendicular word covers.
test(word_checked_count_per_word) :-
    sample_cross(Placed),
    Placed = [Wa, Wd1, _],
    word_checked_count(Wa,  Placed, 2),   % across A is crossed at 2 of its 3 cells
    word_checked_count(Wd1, Placed, 1).   % down D1 only at its shared cell

% word_meets_half/2: at least ceil(len/2) of a word's cells must be checked.
test(word_meets_half_needs_half_its_cells_checked) :-
    sample_cross(Placed),
    Placed = [Wa, Wd1, _],
    crosswordsmith_metrics:word_meets_half(Wa, Placed),          % 2 >= ceil(3/2)=2
    \+ crosswordsmith_metrics:word_meets_half(Wd1, Placed).      % 1  < 2

% word_half_threshold/2 (P10): the single ceil(L/2) definition that both
% word_meets_half and lint's checked_half rule reuse.
test(word_half_threshold_is_ceil_half) :-
    word_half_threshold(1, 1),
    word_half_threshold(3, 2),
    word_half_threshold(4, 2),
    word_half_threshold(5, 3).

% word_max_unch_run/3: longest run of consecutive UNchecked cells in the word.
test(word_max_unch_run_longest_gap) :-
    sample_cross(Placed),
    Placed = [Wa, Wd1, _],
    crosswordsmith_metrics:word_max_unch_run(Wa,  Placed, 1),    % A: checked,unchecked,checked
    crosswordsmith_metrics:word_max_unch_run(Wd1, Placed, 2).    % D1: checked,unchecked,unchecked

% P4: the hoisted canonical bitmap. layout_dir_cells/2 computes both directions
% ONCE; word_checked_bitmap/3 turns it into the word's 1/0 checked flags, and the
% derived count/max-run (bits_*) are the single source of truth that the
% convenience (W, Placed) forms also use - so they must agree.
test(word_checked_bitmap_canonical_and_derived) :-
    sample_cross(Placed),
    Placed = [Wa, Wd1, _],
    layout_dir_cells(Placed, DirCells),
    word_checked_bitmap(Wa,  DirCells, BitsA), BitsA == [1, 0, 1],   % crossed at 1 & 3
    word_checked_bitmap(Wd1, DirCells, BitsD), BitsD == [1, 0, 0],   % crossed at start only
    bits_checked_count(BitsA, 2), bits_max_unch_run(BitsA, 1),
    bits_checked_count(BitsD, 1), bits_max_unch_run(BitsD, 2),
    % identical to the single-direction (W, Placed) forms used by arrange:
    word_checked_count(Wa, Placed, 2), crosswordsmith_metrics:word_max_unch_run(Wd1, Placed, 2).

% --- helpers (word_letters) ---------------------------------------------------

% word_letters/3 yields the space-stripped char list and its length.
test(word_letters_strips_spaces) :-
    word_letters(['NEW YORK', _{}], Letters, WLen),
    Letters == ['N','E','W','Y','O','R','K'],
    WLen =:= 7.

:- end_tests(quality).
