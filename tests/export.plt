% tests/export.plt - plunit suite for export.pl (Flavour-B interchange, §8.2).
%
% Assumes core.pl and export.pl are consulted by the runner before this
% file. Covers the enumeration deriver (AC-ENUM-1), the ipuz v2 transform
% (AC-EXP-1 structure) and the Exolve transform (AC-EXP-2 structure), plus the
% preservation of clue text + enumerations (AC-EXP-3). True third-party
% ingestion (kotwords / Exet) is a manual step; the byte-exact output is pinned
% by the goldens (tests/golden/export_bundled_17.{ipuz,exolve}).

:- use_module(library(plunit)).
% lists: test-body helpers (memberchk); explicit so the suite also runs under
% autoload(false) (P11/C5).
:- use_module(library(lists)).

% A tiny canonical layout: CAT across row0 + COW down col0 of a 3x3, crossing at
% C. Blocks are the bare atom `null` (a non-dict cell), as in a parsed layout.
mini_canonical(
    _{ gridLength:3,
       grid:[
         [ _{letter:"C",number:1,   across:1,   down:1},
           _{letter:"A",number:null,across:1,   down:null},
           _{letter:"T",number:null,across:1,   down:null} ],
         [ _{letter:"O",number:null,across:null,down:1}, null, null ],
         [ _{letter:"W",number:null,across:null,down:1}, null, null ] ],
       words:[
         _{number:1, direction:"across", answer:"CAT", cells:[[0,0],[0,1],[0,2]], meta:_{clue:"Feline"}},
         _{number:1, direction:"down",   answer:"COW", cells:[[0,0],[1,0],[2,0]], meta:_{}} ] }).


:- begin_tests(export).

% --- enumeration deriver (AC-ENUM-1 / AC-EXP-3) ------------------------------
test(enum_single_word)   :- answer_enumeration('CAT', '(3)').
test(enum_two_words)     :- answer_enumeration('OMEGA POINT', '(5,5)').
test(enum_hyphenated)    :- answer_enumeration('WELL-BEING', '(4-5)').
test(enum_mixed_hyphens) :- answer_enumeration('MOTHER-IN-LAW', '(6-2-3)').
test(enum_from_string)   :- answer_enumeration("NEW YORK CITY", '(3,4,4)').

% R5 (revamp audit): leading/trailing/adjacent separators (authoring typos)
% must NOT emit 0-length segments like (2,0).
test(enum_trailing_space)  :- answer_enumeration('TT ', '(2)').
test(enum_leading_space)   :- answer_enumeration(' AB', '(2)').
test(enum_double_space)    :- answer_enumeration('ROCK  ROLL', '(4,4)').
test(enum_trailing_hyphen) :- answer_enumeration('RE-', '(2)').

% --- ipuz v2 (AC-EXP-1 structure) --------------------------------------------

% version / kind / dimensions are the ipuz v2 crossword identifiers.
test(ipuz_header) :-
    mini_canonical(D), layout_to_ipuz(D, I),
    get_dict(version, I, "http://ipuz.org/v2"),
    get_dict(kind, I, ["http://ipuz.org/crossword#1"]),
    get_dict(title, I, "Untitled"),   % non-null title (Exet imports ipuz->exolve; V1)
    get_dict(dimensions, I, Dim), get_dict(width, Dim, 3), get_dict(height, Dim, 3).

% puzzle: clue-start cell -> its number, white-but-unnumbered -> 0, block -> "#".
test(ipuz_puzzle_cell_encoding) :-
    mini_canonical(D), layout_to_ipuz(D, I),
    get_dict(puzzle, I, P),
    P == [[1, 0, 0], [0, '#', '#'], [0, '#', '#']].

% solution: letters for white, "#" for block.
test(ipuz_solution_letters) :-
    mini_canonical(D), layout_to_ipuz(D, I),
    get_dict(solution, I, S),
    S == [["C","A","T"], ["O",'#','#'], ["W",'#','#']].

% clues grouped Across/Down, each carrying number, clue text and enumeration.
test(ipuz_clues_carry_text_and_enumeration) :-
    mini_canonical(D), layout_to_ipuz(D, I),
    get_dict(clues, I, C),
    get_dict('Across', C, [AC]), get_dict('Down', C, [DN]),
    get_dict(number, AC, 1), get_dict(clue, AC, "Feline"), get_dict(enumeration, AC, '(3)'),
    get_dict(number, DN, 1), get_dict(clue, DN, ""),       get_dict(enumeration, DN, '(3)').

% --- Exolve (AC-EXP-2 structure / AC-EXP-3) ----------------------------------

% The grid is letters for white cells and '.' for blocks; the framing directives
% and the across/down clue lines (number, clue text, enumeration) are present.
test(exolve_structure) :-
    mini_canonical(D), layout_to_exolve(D, Text),
    atomic_list_concat(Lines, '\n', Text),
    memberchk('exolve-begin', Lines),
    memberchk('  exolve-title: Untitled', Lines),  % non-null title (Exet Save needs it; AC-EXP-2)
    memberchk('  exolve-width: 3', Lines),
    memberchk('  exolve-height: 3', Lines),
    memberchk('    CAT', Lines),                 % row 0: all white
    memberchk('    O..', Lines),                 % row 1: white then two blocks
    memberchk('    W..', Lines),
    memberchk('  exolve-across:', Lines),
    memberchk('    1 Feline (3)', Lines),         % number + clue + enumeration
    memberchk('  exolve-down:', Lines),
    memberchk('    1  (3)', Lines),               % empty clue (no data invented)
    memberchk('exolve-end', Lines).

% --- validation --------------------------------------------------------------
% export_load/2 is module-internal (not exported); white-box tests reach it
% via module qualification (migration plan, Resolved Decision 3).
test(export_rejects_non_layout, [throws(error(export_invalid_layout, _))]) :-
    tmp_file_stream(text, F, S), write(S, '{"foo": 1}'), close(S),
    call_cleanup(crosswordsmith_export:export_load(F, _), delete_file(F)).

% The C50 gate depth: export_layout_dict/2 is a SHAPE gate over everything the
% transforms dereference structurally — each grid row a list, each words[]
% entry an object carrying answer/direction/cells. Pre-fix, the first payload
% below made both transforms plain-fail (surfacing as browser.pl's mislabelled
% internal "engine bug" envelope) and the second was blessed as a degenerate
% empty-ipuz success; both (and the same-class shapes) must throw the same
% typed formal the presence checks already threw.
test(export_gate_rejects_nonlist_grid_row,          % C50 probe payload 1
     [throws(error(export_invalid_layout, _))]) :-
    export_layout_dict(_{gridLength: 1, grid: [5], words: []}, _).
test(export_gate_rejects_schemaless_word,           % C50 probe payload 2
     [throws(error(export_invalid_layout, _))]) :-
    export_layout_dict(_{gridLength: 3, grid: [], words: [_{bogus: 1}]}, _).
test(export_gate_rejects_word_not_object,
     [throws(error(export_invalid_layout, _))]) :-
    export_layout_dict(_{gridLength: 3, grid: [], words: [7]}, _).
test(export_gate_rejects_word_cells_not_list,
     [throws(error(export_invalid_layout, _))]) :-
    export_layout_dict(_{gridLength: 3, grid: [],
                         words: [_{answer: "CAT", direction: "across",
                                   cells: "garbage"}]}, _).
% ...and the canonical shape still passes, unchanged (the happy-path lock; the
% golden suite pins the full-size layouts).
test(export_gate_accepts_canonical_layout) :-
    mini_canonical(D),
    export_layout_dict(D, D2),
    D2 == D.

:- end_tests(export).
