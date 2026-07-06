% stockgrid.pl - the bundled stock-grid library (design-spec §8.3; schema fixed
% by DP-1 / OD-5). A stock grid is a curated, pre-validated legal blocked
% template, stored as a black-square MASK (the single source of truth):
%
%   { "name": "...", "size": N, "mask": ["#.....", ...] }
%   % one string per row; '#' = block, else light
%
% Required fields: name, size, mask. An optional "symmetry" annotation MAY
% appear in the file (the shipped grids carry "rot180") but is NOT trusted or
% required - stockgrid_load/2 ignores it and the validator RE-DERIVES the actual
% symmetry from the mask, so a declared value is descriptive only.
%
% Slots (lights) are DERIVED on load - not stored - by run-scanning the mask,
% then validated by `lint --profile blocked-uk` (the metric predicates as
% template validators). This is a library (loader + validator); stock grids
% feed `lint` profiles and `fill`. There is no CLI verb.
%
% Exports: fill uses stockgrid_load/2, mask_white_cells/3 and grid_run/4; the
% validate/report trio is deliberate library API. Internals are reached by
% tests as crosswordsmith_stockgrid:Pred(...).

:- module(crosswordsmith_stockgrid,
          [ stockgrid_load/2,
            mask_white_cells/3,
            grid_run/4,
            stockgrid_validate/3,
            stockgrid_validate_file/3,
            stockgrid_report/1
          ]).

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5).
% library(json), NOT the legacy library(http/json) alias — the alias does not
% resolve in the WASM image (C6).
:- use_module(library(json), [json_read_dict/2]).
:- use_module(library(apply), [maplist/2, maplist/3]).
:- use_module(library(lists),
              [append/2, append/3, member/2, nth0/3, numlist/3]).
:- use_module(library(ordsets),
              [list_to_ord_set/2, ord_memberchk/2, ord_subtract/3]).

% lint_run/5: the blocked-uk template validation.
:- use_module(crosswordsmith(lint), [lint_run/5]).

% assign_clue_numbers/2: number the derived lights for validation/reporting.
% pw_cells/2: the placed-word record (pw/8) cell accessor the lights are read via.
:- use_module(crosswordsmith(core), [assign_clue_numbers/2, pw_cells/2]).


% --- load + parse the mask ---------------------------------------------------

%!  stockgrid_load(+File:atom, -Grid) is det.
%
%   Read a stock-grid JSON file into grid(Name, Size, Mask), Mask a list of
%   Size row char-lists. Requires {name, size:N, mask:[N strings of length
%   N]}; throws error(stockgrid_invalid(File), _) otherwise. An optional
%   "symmetry" annotation in the file is ignored (the validator RE-DERIVES
%   the actual symmetry from the mask).
stockgrid_load(File, grid(Name, Size, Mask)) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Dict), close(S)),
    (   is_dict(Dict),
        get_dict(name, Dict, RawName), atom_string(Name, RawName),
        get_dict(size, Dict, Size), integer(Size), Size > 0,
        get_dict(mask, Dict, Rows), is_list(Rows), length(Rows, Size),
        forall(member(R, Rows), ( string(R), string_length(R, Size) ))
    ->  maplist(row_chars, Rows, Mask)          % each row -> list of chars
    ;   throw(error(stockgrid_invalid(File), _))
    ).

row_chars(RowStr, Chars) :- string_chars(RowStr, Chars).

%!  mask_white_cells(+Mask:list, +Size:integer, -WhiteSet:ordset) is det.
%
%   The white (light) cell numbers of a mask (1-based, row-major): any mask
%   char that is not '#'. An ordset, for the run/validation set tests.
mask_white_cells(Mask, Size, WhiteSet) :-
    findall(Cell,
            ( nth0(R, Mask, Chars), nth0(C, Chars, Ch), Ch \== '#',
              Cell is R * Size + C + 1 ),
            Ws),
    list_to_ord_set(Ws, WhiteSet).


% --- derive the lights (maximal white runs of length >= 2) -------------------
% A light is a maximal run of >= 2 consecutive white cells in a direction. A
% length-1 white run is not a light in that direction; the cell must still belong
% to a perpendicular light (else it is an isolated cell - reported separately).
grid_lights(Size, WhiteSet, Numbered) :-
    findall(W,
            ( grid_run(Size, WhiteSet, across, Cells), light_word(across, Cells, W) ),
            Across),
    findall(W,
            ( grid_run(Size, WhiteSet, down, Cells), light_word(down, Cells, W) ),
            Down),
    append(Across, Down, Lights),
    % No once/1: assign_clue_numbers/2 is deterministic at source (X6.B4;
    % probe re-verified on the shipped grids) - the wrap was defensive (C23).
    assign_clue_numbers(Lights, Numbered).

%!  grid_run(+Size:integer, +WhiteSet:ordset, +Dir:oneof([across,down]),
%!           -Cells:list(integer)) is nondet.
%
%   Backtrack over the light cell-lists (maximal white runs of length >= 2)
%   of one direction. Exposed so `fill` can build its slots from the same
%   derivation the validator uses.
grid_run(Size, WhiteSet, across, Cells) :-
    Smax is Size - 1, between(0, Smax, R),
    Lo is R * Size + 1, Hi is R * Size + Size, numlist(Lo, Hi, RowCells),
    line_lights(RowCells, WhiteSet, Cells).
grid_run(Size, WhiteSet, down, Cells) :-
    Smax is Size - 1, between(0, Smax, C),
    col_cells(C, Size, ColCells),
    line_lights(ColCells, WhiteSet, Cells).

col_cells(C, Size, Cells) :-
    Smax is Size - 1,
    findall(Cell, ( between(0, Smax, R), Cell is R * Size + C + 1 ), Cells).

% Backtrack over the maximal white runs (length >= 2) of one line of cell numbers.
line_lights(Cells, WhiteSet, Run) :-
    split_runs(Cells, WhiteSet, Runs),
    member(Run, Runs),
    Run = [_,_|_].   % length >= 2, in O(1) (P7)

split_runs([], _WS, []).
split_runs([Cell|Cs], WS, Runs) :-
    ( ord_memberchk(Cell, WS)
    -> grab_run([Cell|Cs], WS, Run, Rest), Runs = [Run|More], split_runs(Rest, WS, More)
    ;  split_runs(Cs, WS, Runs)
    ).
grab_run([Cell|Cs], WS, [Cell|Run], Rest) :-
    ord_memberchk(Cell, WS), !,
    grab_run(Cs, WS, Run, Rest).
grab_run(Cs, _WS, [], Cs).

% A light as a placed-word record pw/8 (synthetic answer of the right length; the
% metric rules read only cells/dir/len/start, so Letters/End are left unbound).
light_word(Dir, Cells, pw(Synth, _Letters, Cells, Dir, Len, Start, _End, _Num)) :-
    length(Cells, Len), Cells = [Start|_],
    length(As, Len), maplist(=('A'), As), atom_chars(Synth, As).

% White cells covered by at least one light (the union of all light cells).
light_cell_set(Lights, Set) :-
    findall(C, ( member(W, Lights), pw_cells(W, Cs), member(C, Cs) ), All),
    list_to_ord_set(All, Set).


% --- validate (the design-time legality check; reuses lint blocked-uk) -------

%!  stockgrid_validate(+Grid, -Verdict:atom, -Fails:list) is det.
%
%   Validate a grid(Name, Size, Mask) term: derive its lights and run lint's
%   blocked-uk profile over them. Verdict is 'PASS' / 'WARN' / 'FAIL' (an
%   isolated white cell forces 'FAIL'). A grid is legal iff Verdict \==
%   'FAIL'. Fails lists the blocking reasons as fail(Rule, NumOrGrid,
%   CellsOrNone) terms.
stockgrid_validate(grid(_Name, Size, Mask), Verdict, Fails) :-
    mask_white_cells(Mask, Size, WhiteSet),
    grid_lights(Size, WhiteSet, Words),
    light_cell_set(Words, LightCells),
    ord_subtract(WhiteSet, LightCells, Isolated),       % white cells in no light
    lint_run(Words, Size, 'blocked-uk', false, Report),
    get_dict(verdict, Report, V),
    ( Isolated == [] -> Verdict = V ; Verdict = 'FAIL' ),
    % lint_run reports words in input order, so pair each report entry back to
    % its source pw/8 POSITIONALLY (nth0/3 as an in-order index/element
    % generator) - no re-search of Ws, and no reliance on report dicts being
    % pairwise distinct.
    findall(fail(Rule, Num, Cells),
            ( get_dict(words, Report, Ws),
              nth0(I, Ws, WR),
              nth0(I, Words, SrcW),
              get_dict(results, WR, Rs), member(Res, Rs),
              get_dict(severity, Res, 'FAIL'), get_dict(rule, Res, Rule),
              get_dict(number, WR, Num),
              pw_cells(SrcW, Cells) ),
            WordFails),
    findall(fail(Rule, grid, none),
            ( get_dict(grid, Report, GRs), member(GR, GRs),
              get_dict(severity, GR, 'FAIL'), get_dict(rule, GR, Rule) ),
            GridFails),
    ( Isolated == [] -> IsoFails = [] ; IsoFails = [fail(isolated_cells, grid, Isolated)] ),
    append([WordFails, GridFails, IsoFails], Fails).

%!  stockgrid_validate_file(+File:atom, -Verdict:atom, -Fails:list) is det.
%
%   stockgrid_load/2 + stockgrid_validate/3 on one grid file. Throws
%   stockgrid_invalid on a malformed file.
stockgrid_validate_file(File, Verdict, Fails) :-
    stockgrid_load(File, Grid),
    stockgrid_validate(Grid, Verdict, Fails).

%!  stockgrid_report(+File:atom) is semidet.
%
%   Design-time report on stderr for one grid file (the "draft -> lint ->
%   iterate" loop); succeeds iff the grid is legal (Verdict \== 'FAIL').
stockgrid_report(File) :-
    stockgrid_load(File, Grid),
    Grid = grid(Name, Size, _Mask),
    stockgrid_validate(Grid, Verdict, Fails),
    format(user_error, "~w (~wx~w): ~w~n", [Name, Size, Size, Verdict]),
    forall(member(F, Fails), format(user_error, "    ~q~n", [F])),
    Verdict \== 'FAIL'.


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(stockgrid_invalid(File)) -->
    [ 'stock grid ~q: needs {name, size:N, mask:[N strings of length N]}'-[File] ].
