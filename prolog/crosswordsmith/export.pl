% export.pl - Flavour-B standard-format interchange (design-spec §8.2). These
% are TRANSFORMATIONS of the canonical layout JSON (what `arrange` emits), not
% new emitters:
%
%   --to ipuz   : ipuz v2 crossword (a JSON object; reachable onward to
%                 .puz/.jpz/PDF via off-the-shelf kotwords, §3).
%   --to exolve : Exolve plain text (git-diffable; round-trips to Exet).
%
% Enumerations are derived from the answer's spaces/hyphens (§6.3); clue text
% rides in meta.clue. Nothing is invented (AC-EXP-3): a missing clue is empty.
% Zero project dependencies (it only needs JSON + list/pairs utilities).
%
% Exports: export_solve/2 is the CLI seam; export_layout_dict/2 is the
% file-free shape gate (browser.pl's params path); the layout_to_*/2
% transforms and answer_enumeration/2 are deliberate library API. Everything
% else is internal — tests reach internals as crosswordsmith_export:Pred(...).

:- module(crosswordsmith_export,
          [ export_solve/2,
            export_layout_dict/2,
            layout_to_ipuz/2,
            layout_to_exolve/2,
            answer_enumeration/2
          ]).

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5).
% library(json), NOT the legacy library(http/json) alias — the alias does not
% resolve in the WASM image (C6).
:- use_module(library(json), [json_read_dict/2, json_write_dict/2]).
:- use_module(library(apply), [maplist/2, maplist/3]).
:- use_module(library(lists), [append/2, member/2]).
:- use_module(library(pairs), [pairs_values/2]).


% --- load + validate a canonical layout --------------------------------------
export_load(File, Dict) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Dict0), close(S)),
    export_layout_dict(Dict0, Dict).

% The shape gate, factored file-free so the browser path validates through the
% SAME predicate the CLI does (no duplicated validation to drift). Depth (C50):
% a cheap SHAPE check on everything the transforms dereference structurally —
% each grid row must be a list and each words[] entry an object carrying
% answer/direction/cells — so gate-passing garbage can no longer surface as a
% mislabelled internal "engine bug" (a non-list grid row made both transforms
% plain-fail into browser.pl's never-plain-fail backstop) or as a degenerate
% success (a schema-less words entry was silently skipped by the clue
% collectors, blessing garbage with an empty ipuz). CONTENT stays the
% transforms' domain: on a malformed cell they still fail by design with no
% partial output (P17) — this is a shape gate, not schema validation.
export_layout_dict(Dict0, Dict) :-
    (   is_dict(Dict0),
        get_dict(gridLength, Dict0, N), integer(N), N > 0,
        get_dict(grid, Dict0, Rows), is_list(Rows),
        maplist(is_list, Rows),
        get_dict(words, Dict0, Ws), is_list(Ws),
        maplist(export_word_shape, Ws)
    ->  Dict = Dict0
    ;   throw(error(export_invalid_layout, _))
    ).

% One canonical words[] entry: an object with the keys the transforms (and the
% enumeration deriver) read. Mirrors lint_word/3's required key set — answer,
% direction, cells — `number` is optional there and likewise unchecked here.
export_word_shape(W) :-
    is_dict(W),
    get_dict(answer, W, _),
    get_dict(direction, W, _),
    get_dict(cells, W, Cells),
    is_list(Cells).

% A grid cell is white iff it is a dict (a filled cell object); anything else
% (JSON null, however the reader represents it) is a black square. This keeps
% the transforms independent of the null representation.
white_cell(Cell) :- is_dict(Cell).

% Match a word's (string) direction against an atom.
word_dir(W, Dir) :-
    get_dict(direction, W, Raw),
    ( atom(Raw) -> A = Raw ; atom_string(A, Raw) ),
    A == Dir.

% Clue text from meta.clue, or "" when absent (no data invented).
word_clue_text(W, Text) :-
    ( get_dict(meta, W, M), is_dict(M), get_dict(clue, M, C) -> Text = C ; Text = "" ).


% --- enumeration: derive "(5,5)" / "(4-5)" from spaces/hyphens (§6.3) ---------
answer_enumeration(Answer, Enum) :-
    ( atom(Answer) -> atom_chars(Answer, Chars) ; string_chars(Answer, Chars) ),
    enum_tokens(Chars, 0, Tokens),
    enum_clean(Tokens, Segs),
    ( Segs == [] -> Body = 0 ; atomic_list_concat(Segs, Body) ),
    atomic_list_concat(['(', Body, ')'], Enum).

% Tokenize into letter-run lengths and the separators between them (',' for a
% space, '-' for a hyphen): the result is R0, S0, R1, S1, ..., Rk - always
% run-delimited (runs may be 0 for adjacent/leading/trailing separators).
enum_tokens([], N, [N]).
enum_tokens([C|Cs], N, Toks) :-
    ( sep_char(C, S) -> Toks = [N, S | Rest], enum_tokens(Cs, 0, Rest)
    ; N1 is N + 1, enum_tokens(Cs, N1, Toks)
    ).

sep_char(' ', ',').
sep_char('-', '-').

% Keep only positive runs; a separator is emitted ONLY between two kept runs
% (the one immediately preceding the later run). Leading/trailing/adjacent
% separators and 0-length runs (authoring typos) thus produce no empty segment
% like (2,0) (R5); an all-separator answer collapses to body 0.
enum_clean([R0|Pairs], Clean) :-
    ( R0 > 0 -> Clean = [R0|More], enum_clean_pairs(Pairs, true, More)
    ;          enum_clean_pairs(Pairs, false, Clean) ).

enum_clean_pairs([], _Emitted, []).
enum_clean_pairs([S, R | Rest], Emitted, Out) :-
    ( R > 0
    ->  ( Emitted == true -> Out = [S, R | More] ; Out = [R | More] ),
        enum_clean_pairs(Rest, true, More)
    ;   enum_clean_pairs(Rest, Emitted, Out)
    ).


% --- ipuz v2 -----------------------------------------------------------------
% Invert null->"#", split the merged cell into parallel puzzle/solution arrays,
% group words by direction into a clues dict, add version/kind/dimensions.
% ipuz carries no symmetry field - symmetry stays a crosswordsmith lint concept.
layout_to_ipuz(Dict,
    _{ version: "http://ipuz.org/v2",
       kind: ["http://ipuz.org/crossword#1"],
       % A default title (optional in ipuz, so kotwords is unaffected). Exet
       % imports ipuz by converting it to Exolve and its format-agnostic
       % fileTitle() crashes on a null title, so a non-null one keeps the ipuz
       % round-trippable through Exet too - mirrors the Exolve default (V1).
       title: "Untitled",
       dimensions: _{width:N, height:N},
       empty: 0,
       puzzle: Puzzle,
       solution: Solution,
       clues: _{'Across':Across, 'Down':Down} }) :-
    get_dict(gridLength, Dict, N),
    get_dict(grid, Dict, Rows),
    maplist(ipuz_puzzle_row, Rows, Puzzle),
    maplist(ipuz_solution_row, Rows, Solution),
    get_dict(words, Dict, Words),
    ipuz_clues(Words, across, Across),
    ipuz_clues(Words, down, Down).

ipuz_puzzle_row(Row, Out) :- maplist(ipuz_puzzle_cell, Row, Out).
ipuz_puzzle_cell(Cell, Out) :-
    ( white_cell(Cell)
    -> ( get_dict(number, Cell, Num), integer(Num) -> Out = Num ; Out = 0 )
    ;  Out = '#' ).

% Invariant: a white cell always carries `letter` and each word carries `answer`
% (the canonical layout contract, structure-checked by export_load/2 up front).
% So the get_dict/3 reads here and in exolve_cell_char / ipuz_clue /
% exolve_clue_line need no fallback; on a malformed cell the transform simply
% fails (no partial ipuz/exolve is emitted) rather than throwing (P17).
ipuz_solution_row(Row, Out) :- maplist(ipuz_solution_cell, Row, Out).
ipuz_solution_cell(Cell, Out) :-
    ( white_cell(Cell) -> get_dict(letter, Cell, Out) ; Out = '#' ).

% Collect one item per word in direction Dir, ordered by clue number. MkItem is
% called as MkItem(W, Num, Item); shared by the ipuz and Exolve clue collectors
% (same findall -> keysort -> pairs_values skeleton, only the per-item goal differs).
:- meta_predicate collect_by_number(+, +, 3, -).
collect_by_number(Words, Dir, MkItem, Items) :-
    findall(Num-Item,
            ( member(W, Words), word_dir(W, Dir), call(MkItem, W, Num, Item) ),
            Pairs),
    keysort(Pairs, Sorted),
    pairs_values(Sorted, Items).

% Clue objects {number, clue, enumeration}, sorted by number.
ipuz_clues(Words, Dir, Clues) :- collect_by_number(Words, Dir, ipuz_clue, Clues).
ipuz_clue(W, Num, _{number:Num, clue:Text, enumeration:Enum}) :-
    get_dict(number, W, Num),
    word_clue_text(W, Text),
    get_dict(answer, W, Answer),
    answer_enumeration(Answer, Enum).


% --- Exolve ------------------------------------------------------------------
% Plain text: a grid of letters (white) and '.' (block), then the across/down
% clue lines (number, clue, enumeration). Exolve numbers the grid itself; our
% standard numbering matches, so the clue-line numbers align.
layout_to_exolve(Dict, Text) :-
    get_dict(gridLength, Dict, N),
    get_dict(grid, Dict, Rows),
    get_dict(words, Dict, Words),
    maplist(exolve_grid_row, Rows, GridLines),
    exolve_clue_lines(Words, across, AcrossLines),
    exolve_clue_lines(Words, down, DownLines),
    format(atom(WLine), "  exolve-width: ~w", [N]),
    format(atom(HLine), "  exolve-height: ~w", [N]),
    append([ ['exolve-begin',
              '  exolve-id: crosswordsmith-export',
              % A default title: Exet requires a non-null title to save (its
              % fileTitle() does title.replace(...) and crashes on a null/absent
              % title), so emit one even though the canonical layout has none.
              % It is an editable placeholder in any downstream tool. (AC-EXP-2.)
              '  exolve-title: Untitled',
              WLine, HLine,
              '  exolve-grid:' ],
             GridLines,
             ['  exolve-across:'], AcrossLines,
             ['  exolve-down:'],   DownLines,
             ['exolve-end'] ],
           AllLines),
    atomic_list_concat(AllLines, '\n', Body),
    atom_concat(Body, '\n', Text).

exolve_grid_row(Row, Line) :-
    maplist(exolve_cell_char, Row, Chars),
    atomic_list_concat(Chars, Body),
    atom_concat('    ', Body, Line).
exolve_cell_char(Cell, Ch) :-
    ( white_cell(Cell) -> get_dict(letter, Cell, L), atom_string(Ch, L) ; Ch = '.' ).

exolve_clue_lines(Words, Dir, Lines) :- collect_by_number(Words, Dir, exolve_clue_line, Lines).
exolve_clue_line(W, Num, Line) :-
    get_dict(number, W, Num),
    word_clue_text(W, Clue),
    get_dict(answer, W, Answer),
    answer_enumeration(Answer, Enum),
    format(atom(Line), "    ~w ~w ~w", [Num, Clue, Enum]).


% --- entry point -------------------------------------------------------------
export_solve(File, ipuz) :-
    export_load(File, Dict),
    layout_to_ipuz(Dict, Ipuz),
    current_output(Out),
    json_write_dict(Out, Ipuz),
    nl(Out).
export_solve(File, exolve) :-
    export_load(File, Dict),
    layout_to_exolve(Dict, Text),
    current_output(Out),
    write(Out, Text).


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(export_invalid_layout) -->
    [ 'export: input is not a canonical layout (needs gridLength, grid rows, and words each carrying answer/direction/cells)' ].
