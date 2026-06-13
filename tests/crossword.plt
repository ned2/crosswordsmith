% tests/crossword.plt - plunit test suite for crossword.pl
%
% These tests assume crossword.pl has already been consulted into the
% `user` module (the runner, tests/run_tests.pl, does this before loading
% this file). Run them via:
%
%     ./run_tests.sh        (or)        make test
%
% The tests are grouped into units: pure grid geometry, generic list
% utilities, the solver, and clue numbering.

:- use_module(library(plunit)).
:- use_module(library(http/json)).


% Grid geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cells are numbered 1..GridLen*GridLen in row-major order. `across` steps
% by +1 along a row; `down` steps by +GridLen down a column.

:- begin_tests(geometry).

test(swap_dir_across) :- swap_dir(across, down).
test(swap_dir_down)   :- swap_dir(down, across).

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
test(calc_num_across, [true(N =:= 7)])  :- calc_num(across, 17, 3, 5, N).
test(calc_num_down,   [true(N =:= 35)]) :- calc_num(down, 17, 3, 1, N).
test(calc_roundtrip_across, [true(S =:= 5)]) :-
    calc_num(across, 17, 3, 5, N), calc_start(across, 17, 3, N, S).
test(calc_roundtrip_down, [true(S =:= 1)]) :-
    calc_num(down, 17, 4, 1, N), calc_start(down, 17, 4, N, S).

test(next_cell_across, [true(C =:= 6)])  :- next_cell(across, 5, 17, C).
test(next_cell_down,   [true(C =:= 22)]) :- next_cell(down, 5, 17, C).
test(prev_cell_across, [true(C =:= 4)])  :- prev_cell(across, 5, 17, C).
test(prev_cell_down,   [true(C =:= 5)])  :- prev_cell(down, 22, 17, C).

% First column starts across rows; first row starts down columns.
test(is_start_across)            :- is_start_cell(across, 18, 17).   % col 1, row 2
test(is_start_across_no, [fail]) :- is_start_cell(across, 19, 17).
test(is_start_down)              :- is_start_cell(down, 5, 17).      % row 1
test(is_start_down_no, [fail])   :- is_start_cell(down, 18, 17).
test(is_end_across)              :- is_end_cell(across, 17, 17).     % last col
test(is_end_down)                :- is_end_cell(down, 289, 17).      % last row

% The four named start locations resolve to the expected cell + direction.
test(start_topleft_across, [true(N-D == 1-across)])  :- start_loc(topleft_across, 17, N, D).
test(start_topleft_down,   [true(N-D == 1-down)])    :- start_loc(topleft_down, 17, N, D).
test(start_topright,       [true(N-D == 17-down)])   :- start_loc(topright, 17, N, D).
test(start_bottomleft,     [true(N-D == 273-across)]):- start_loc(bottomleft, 17, N, D).

:- end_tests(geometry).


% Generic list utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(utilities).

% position/3 yields every 1-based index of X, on backtracking.
test(position_first, [true(P == 1)])  :- once(position(a, [a,b,c], P)).
test(position_all,   [all(P == [1,3])]) :- position(a, [a,b,a], P).
test(position_none,  [fail])          :- position(z, [a,b,c], _).

% remove_x removes only the first occurrence. (nondet: it leaves a harmless
% choicepoint; we only care about the first answer.)
test(remove_x_first, [true(R == [a,c,b]), nondet]) :- remove_x(b, [a,b,c,b], R).
test(remove_x_absent, [true(R == [a,b,c]), nondet]) :- remove_x(z, [a,b,c], R).

% shuffle is a permutation of its input (we cannot assert order: it is
% random, so just check the multiset of elements is preserved).
test(shuffle_is_permutation, [true(Sorted == [1,2,3,4,5])]) :-
    shuffle([1,2,3,4,5], Out), msort(Out, Sorted).

% init_grid builds a GridLen*GridLen assoc of `empty` cells.
test(init_grid_size, [true(Len =:= 9)]) :-
    init_grid(3, G), assoc_to_list(G, L), length(L, Len).
test(init_grid_empty) :-
    init_grid(3, G), get_assoc(1, G, empty), get_assoc(9, G, empty).
test(init_grid_no_extra, [fail]) :-
    init_grid(3, G), get_assoc(10, G, _).

:- end_tests(utilities).


% The solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- begin_tests(solver).

% Two words that cross on a shared letter are both placed. OMEGA POINT
% (across, starting cell 1) and GNOSTIC GOSPELS (down) cross on G.
% (nondet: find_crossword backtracks over many solutions; we assert one exists.)
test(places_two_crossing_words, [nondet]) :-
    Words = [['OMEGA POINT', _{}], ['GNOSTIC GOSPELS', _{}]],
    find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 2),
    % the first word is laid across starting at cell 1
    member(PW, Placed), get_dict(dir, PW, across), get_dict(start, PW, 1).

% The full bundled clue set has a solution on a 17x17 grid.
test(bundled_clues_solve_at_17, [nondet]) :-
    clues(Words),
    find_crossword(17, Words, topleft_across, _Grid, Placed),
    length(Placed, 6).

% A grid that is too small for the words has no solution.
test(too_small_grid_fails, [fail]) :-
    clues(Words),
    find_crossword(3, Words, topleft_across, _Grid, _Placed).

% End-to-end: the top-level crossword/3 emits one JSON object describing the
% solution. Parse it and assert the expected shape. Output is captured both to
% keep logs clean and to inspect it.
test(full_pipeline_emits_json) :-
    clues(Words),
    with_output_to(string(S), crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Grid), length(Grid, 17),
    forall(member(Row, Grid), length(Row, 17)),
    get_dict(words, Dict, WordObjs), length(WordObjs, 6).

% A crossing cell carries BOTH directions, so both links remain reachable (G2).
% Row 0, col 3 is where OMEGA POINT (across 1) and GNOSTIC GOSPELS (down 2) cross.
test(crossing_cell_has_both_directions) :-
    clues(Words),
    with_output_to(string(S), crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(grid, Dict, Grid),
    nth0(0, Grid, Row0),
    nth0(3, Row0, Cell),
    get_dict(across, Cell, 1),
    get_dict(down, Cell, 2).

% Per-word metadata is rejoined verbatim under `meta`.
test(word_metadata_under_meta) :-
    clues(Words),
    with_output_to(string(S), crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(words, Dict, WordObjs),
    once(( member(W, WordObjs), get_dict(answer, W, "OMEGA POINT") )),
    get_dict(meta, W, Meta),
    get_dict(clue, Meta, "Transcending entropy"),
    get_dict(link, Meta, "http://en.wikipedia.org/wiki/Omega_Point").

% The emit-time join requires unique answers; a duplicate is rejected up front.
test(duplicate_answer_rejected, [throws(error(duplicate_answer('CAT'), _))]) :-
    crossword(5, [['CAT', _{}], ['CAT', _{}]], topleft_across).

:- end_tests(solver).


% Clue numbering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A placed word is word{answer,letters,cells,dir,len,start,num}.
% assign_clue_numbers/2 sorts by start and fills in num.

:- begin_tests(clue_numbering).

% Helper: clue number assigned to the word with the given text.
clue_num_of(Word, Placed, Num) :-
    member(PW, Placed), get_dict(answer, PW, Word), get_dict(num, PW, Num).

% Two words with distinct start cells get sequential numbers, in start order.
test(distinct_starts_numbered_in_order, [nondet]) :-
    W1 = word{answer:'CAT', letters:[c,a,t], cells:[1,2,3],   dir:across, len:3, start:1},
    W2 = word{answer:'DOG', letters:[d,o,g], cells:[5,12,19], dir:down,   len:3, start:5},
    assign_clue_numbers([W2,W1], Placed),   % deliberately unsorted input
    clue_num_of('CAT', Placed, 1),
    clue_num_of('DOG', Placed, 2).

% Regression test for the (now-fixed) add_clue_nums/3 arity bug: when one
% cell is the start of BOTH an across and a down word, the two clues must
% share a clue number. Before the fix the buggy clause had the wrong arity
% and this scenario failed clue numbering entirely.
test(shared_start_cell_shares_number, [nondet]) :-
    Wa = word{answer:'CAT', letters:[c,a,t], cells:[1,2,3],   dir:across, len:3, start:1},
    Wd = word{answer:'COW', letters:[c,o,w], cells:[1,18,35], dir:down,   len:3, start:1},
    assign_clue_numbers([Wa,Wd], Placed),
    clue_num_of('CAT', Placed, N),
    clue_num_of('COW', Placed, N),
    integer(N).

:- end_tests(clue_numbering).


% JSON clue input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The --clues path reads a JSON file into the same internal Words list the
% bundled clues/1 produces, so the whole pipeline downstream is unchanged.
% Paths are relative to the repo root (the runner consults from there).

:- begin_tests(json_input).

% Reduce one emitted words[] entry back to an input clue {answer, meta}.
reduce_word(W, _{answer:A, meta:M}) :-
    get_dict(answer, W, A),
    get_dict(meta, W, M).

% The fixture loads to the internal Words form: atom answers, dict metadata.
test(json_loads_to_words) :-
    read_clues_json('tests/clues.json', Words),
    length(Words, 6),
    once((member([Answer, Meta], Words), Answer == 'OMEGA POINT')),
    atom(Answer),
    % json_read_dict yields string meta *values* (vs the Prolog path's atoms);
    % json_write_dict renders both to the same JSON token (input spec D3).
    get_dict(clue, Meta, "Transcending entropy"),
    get_dict(link, Meta, "http://en.wikipedia.org/wiki/Omega_Point").

% An omitted `meta` key and an explicit empty `{}` both yield an empty dict.
test(json_absent_meta_is_empty_dict) :-
    read_clues_json('tests/clues.json', Words),
    once((member([A1, M1], Words), A1 == 'GNOSTIC GOSPELS')),  % no meta key
    once((member([A2, M2], Words), A2 == 'BIAS')),             % meta: {}
    dict_pairs(M1, _, []),
    dict_pairs(M2, _, []).

% The JSON-loaded words drive the existing pipeline to a valid JSON solution.
test(json_pipeline_emits_json) :-
    read_clues_json('tests/clues.json', Words),
    with_output_to(string(S), crossword(17, Words, topleft_across)),
    atom_json_dict(S, Dict, []),
    get_dict(gridLength, Dict, 17),
    get_dict(grid, Dict, Grid), length(Grid, 17),
    forall(member(Row, Grid), length(Row, 17)),
    get_dict(words, Dict, WordObjs), length(WordObjs, 6).

% Schema violations throw error/2 terms (rendered by error_message//1).
test(json_no_clues_array_throws, [throws(error(json_no_clues_array, _))]) :-
    doc_to_words(_{version: 1}, _).            % no `clues` key

test(json_clues_not_list_throws, [throws(error(json_no_clues_array, _))]) :-
    doc_to_words(_{clues: "nope"}, _).         % `clues` not an array

test(json_missing_answer_throws, [throws(error(json_invalid_answer(_), _))]) :-
    doc_to_words(_{clues: [_{meta: _{}}]}, _). % entry with no `answer`

test(json_non_string_answer_throws, [throws(error(json_invalid_answer(_), _))]) :-
    doc_to_words(_{clues: [_{answer: 42}]}, _).% numeric answer (would coerce)

test(json_non_object_meta_throws, [throws(error(json_invalid_meta('FLOW'), _))]) :-
    doc_to_words(_{clues: [_{answer: "FLOW", meta: "no"}]}, _).

% Malformed JSON and a missing file surface as standard ISO errors.
test(json_malformed_throws, [throws(_)]) :-
    setup_call_cleanup(open_string("{ not json", S),
                       json_read_dict(S, _),
                       close(S)).

test(json_missing_file_throws, [throws(error(existence_error(source_sink, _), _))]) :-
    read_clues_json('tests/no_such_file.json', _).

% Symmetry (G5): reduce the emitted output to {answer, meta} entries and
% confirm they round-trip as valid input — they load and re-solve, the answer
% set is preserved, and metadata survives. (Asserted as validity, not a byte
% diff: solver layout depends on input order, which the reduction reshuffles.)
test(json_output_reduces_to_valid_input) :-
    read_clues_json('tests/clues.json', Words),
    with_output_to(string(S1), crossword(17, Words, topleft_across)),
    atom_json_dict(S1, Dict, []),
    get_dict(words, Dict, WordObjs),
    maplist(reduce_word, WordObjs, ReducedClues),
    doc_to_words(_{clues: ReducedClues}, Words2),
    % the reduced output is itself valid input: it re-solves and re-emits JSON
    with_output_to(string(S2), crossword(17, Words2, topleft_across)),
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

% load_clues/2: '' selects the bundled clues/1; a path reads the JSON file.
test(load_clues_default_is_bundled) :-
    load_clues('', W),
    length(W, 6),
    W = [[A|_]|_], A == 'OMEGA POINT'.

test(load_clues_from_file) :-
    load_clues('tests/clues.json', W),
    length(W, 6),
    once((member([Ans, _], W), Ans == 'BIAS')).

% valid_loc/1 accepts the four named start locations and nothing else.
test(valid_loc_accepts_known) :- valid_loc(topleft_across).
test(valid_loc_rejects_unknown, [fail]) :- valid_loc(nowhere).

% with_output('', Goal) writes straight to the current output (stdout path).
test(with_output_stdout_passthrough) :-
    with_output_to(string(S), with_output('', write(hello))),
    S == "hello".

% with_output(File, Goal) sends Goal's output to File, byte-for-byte the same
% as the stdout path would produce.
test(with_output_file_matches_stdout) :-
    read_clues_json('tests/clues.json', Words),
    with_output_to(string(StdoutText), crossword(17, Words, topleft_across)),
    tmp_file_stream(text, Tmp, S0), close(S0),
    with_output(Tmp, crossword(17, Words, topleft_across)),
    setup_call_cleanup(open(Tmp, read, S), read_string(S, _, FileText), close(S)),
    delete_file(Tmp),
    FileText == StdoutText,
    atom_json_dict(FileText, Dict, []),         % and it is valid JSON
    get_dict(words, Dict, Ws), length(Ws, 6).

% On a goal that fails (grid too small for a layout), no file is written, so a
% failed run never leaves an empty output file behind.
test(with_output_no_file_on_failure) :-
    read_clues_json('tests/clues.json', Words),
    tmp_file_stream(text, Tmp, S0), close(S0),
    delete_file(Tmp),                           % ensure it is absent
    \+ with_output(Tmp, crossword(3, Words, topleft_across)),
    \+ exists_file(Tmp).

:- end_tests(cli).
