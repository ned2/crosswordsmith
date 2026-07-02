% lint.pl - Flavour-B grid validator (design-spec §8.1). Consumes a canonical
% layout (what `arrange` emits) and reports PASS / WARN / FAIL per rule, per
% word, plus a summary verdict, under a named profile. It needs no engine: it
% reuses core.pl's geometry/JSON and metrics.pl's shared metric predicates
% (word_checked_bitmap, layout_dir_cells, word_checked_count, dir_cells, ...).
% Consult AFTER core.pl (and metrics.pl).

:- use_module(library(http/json)).
:- use_module(library(apply)).
:- use_module(library(aggregate)).
:- use_module(library(ordsets)).


% --- profiles: EXACTLY the rule/severity set each one applies (AC-LINT-4) ----
% Per-word rules: min_length, checked_half, checked_full, max_unch_run,
% double_unch_ends, odd_even. Grid rules: connectivity, symmetry. A rule absent
% from a profile is simply not evaluated; severity is fail | warn. Advisory
% checks (odd_even) never sit at fail.
lint_profile(toc,
    [ min_length-warn, checked_half-warn, max_unch_run-warn,
      double_unch_ends-warn, odd_even-warn, connectivity-warn, symmetry-warn ]).
lint_profile('blocked-uk',
    [ min_length-fail, checked_half-fail, max_unch_run-fail,
      double_unch_ends-warn, odd_even-warn, connectivity-fail, symmetry-fail ]).
lint_profile(american,
    [ min_length-fail, checked_full-fail, connectivity-fail, symmetry-fail ]).
% barred-ximenean (OD-7, resolved): the Ximenean per-length unchecked-letter band
% (checked_band) with RELAXED symmetry (advisory). The band table is primary-
% sourced (Ximenes 1966 + Azed); per-publication symmetry codes are not modelled
% (a documented v1 simplification - symmetry is advisory here).
lint_profile('barred-ximenean',
    [ min_length-fail, checked_band-fail, connectivity-fail, symmetry-warn ]).

lint_known_profile(P) :- lint_profile(P, _).

per_word_rule(min_length).
per_word_rule(checked_half).
per_word_rule(checked_full).
per_word_rule(checked_band).
per_word_rule(max_unch_run).
per_word_rule(double_unch_ends).
per_word_rule(odd_even).
grid_rule(connectivity).
grid_rule(symmetry).


% --- load a canonical layout (the emit format) into placed-word dicts --------
lint_load(File, GridLen, PlacedWords) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Dict), close(S)),
    lint_dict_layout(Dict, GridLen, PlacedWords).

lint_dict_layout(Dict, GridLen, PlacedWords) :-
    (   is_dict(Dict), get_dict(gridLength, Dict, GridLen),
        integer(GridLen), GridLen > 0
    ->  true
    ;   throw(error(lint_no_grid_length, _))
    ),
    (   get_dict(words, Dict, Ws), is_list(Ws)
    ->  true
    ;   throw(error(lint_no_words_array, _))
    ),
    maplist(lint_word(GridLen), Ws, PlacedWords).

% One canonical words[] entry -> word{answer, dir, cells (1-based nums), len, num}.
lint_word(GridLen, Entry, word{answer:A, dir:Dir, cells:CellNums, len:Len, num:Num}) :-
    (   is_dict(Entry), get_dict(answer, Entry, RawA), lint_atom(RawA, A)
    ->  true
    ;   throw(error(lint_invalid_word(Entry), _))
    ),
    (   get_dict(direction, Entry, RawD), lint_dir(RawD, Dir)
    ->  true
    ;   throw(error(lint_invalid_direction(A), _))
    ),
    (   get_dict(cells, Entry, RawCells), is_list(RawCells), RawCells = [_|_]
    ->  true
    ;   throw(error(lint_no_cells(A), _))
    ),
    maplist(lint_cell_num(GridLen, A), RawCells, CellNums),
    length(CellNums, Len),
    ( get_dict(number, Entry, Num0), integer(Num0) -> Num = Num0 ; Num = 0 ).

lint_atom(X, A) :- atom(X), !, A = X.
lint_atom(X, A) :- string(X), !, atom_string(A, X).
lint_dir(X, D) :- lint_atom(X, A), ( A == across -> D = across ; A == down, D = down ).
lint_cell_num(GridLen, A, Pair, Num) :-
    (   Pair = [R, C], integer(R), integer(C),
        R >= 0, R < GridLen, C >= 0, C < GridLen
    ->  Num is R * GridLen + C + 1
    ;   throw(error(lint_invalid_cell(A, Pair), _))
    ).

% White (filled) cells = the union of all word cells (the black-square pattern
% is the complement). An ordset, for the connectivity/symmetry set tests.
filled_cells(PlacedWords, Filled) :-
    findall(C, ( member(W, PlacedWords), get_dict(cells, W, Cs), member(C, Cs) ), All),
    sort(All, Filled).


% --- the run: build the report dict ------------------------------------------
lint_run(PlacedWords, GridLen, Profile, AllowAsym, Report) :-
    lint_profile(Profile, Rules),
    filled_cells(PlacedWords, Filled),
    % compute both directions' covered-cell sets ONCE (P4): every checked-based
    % per-word rule derives from the word's bitmap, so dir_cells/3 (a findall+sort
    % over ALL placed words) is not re-run per rule per word.
    layout_dir_cells(PlacedWords, DirCells),
    findall(_{number:Num, direction:Dir, answer:A, results:Results},
            ( member(W, PlacedWords),
              word_rule_results(Rules, W, DirCells, Results),
              get_dict(answer, W, A), get_dict(num, W, Num), get_dict(dir, W, Dir) ),
            WordReports),
    findall(GR,
            ( member(Rule-Sev0, Rules), grid_rule(Rule),
              symmetry_override(Rule, Sev0, AllowAsym, Sev),
              eval_grid_rule(Rule, Sev, Filled, GridLen, Res),
              result_to_dict(Res, GR) ),
            GridReports),
    all_severities(WordReports, GridReports, Sevs),
    worst_severity(Sevs, Verdict),
    tally(Sevs, Pass, Warn, Fail),
    Report = _{ profile:Profile, allowAsymmetry:AllowAsym, verdict:Verdict,
                summary:_{pass:Pass, warn:Warn, fail:Fail},
                grid:GridReports, words:WordReports }.

word_rule_results(Rules, W, DirCells, Results) :-
    % one bitmap per word; every per-word checked rule reads it (double_unch_ends
    % / odd_even directly, checked_* / max_unch_run via the derived count/run).
    word_checked_bitmap(W, DirCells, Bits),
    findall(RD,
            ( member(Rule-Sev, Rules), per_word_rule(Rule),
              eval_word_rule(Rule, Sev, W, Bits, Res),
              result_to_dict(Res, RD) ),
            Results).

% --allow-asymmetry never lets symmetry hard-FAIL (AC-LINT-3): downgrade to warn.
symmetry_override(symmetry, fail, true, warn) :- !.
symmetry_override(_Rule, Sev, _Allow, Sev).

result_to_dict(result(Rule, Sev, Detail), D) :-
    upcase_sev(Sev, SU),
    ( Detail == null -> D = _{rule:Rule, severity:SU}
    ; D = _{rule:Rule, severity:SU, detail:Detail} ).

upcase_sev(pass, 'PASS').
upcase_sev(warn, 'WARN').
upcase_sev(fail, 'FAIL').


% --- per-word rule evaluators (Sev is the configured violation severity) -----
% Each yields result(Rule, pass|Sev, Detail) - Detail null on PASS.
% Sev is the configured violation severity; Bits is the word's precomputed
% checked bitmap (word_rule_results/4). checked-count and max-unchecked-run come
% from the shared bitmap primitives (metrics.pl); the two bit-pattern rules read
% Bits directly.
eval_word_rule(min_length, CS, W, _Bits, result(min_length, Sev, Detail)) :-
    get_dict(len, W, L),
    ( L >= 3 -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail), "length ~w is below 3", [L]) ).
eval_word_rule(checked_half, CS, W, Bits, result(checked_half, Sev, Detail)) :-
    get_dict(len, W, L), bits_checked_count(Bits, C), word_half_threshold(L, T),
    ( C >= T -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail), "checked ~w of ~w, need ~w", [C, L, T]) ).
eval_word_rule(checked_full, CS, W, Bits, result(checked_full, Sev, Detail)) :-
    get_dict(len, W, L), bits_checked_count(Bits, C),
    ( C >= L -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail), "checked ~w of ~w (every cell must be checked)", [C, L]) ).
% The Ximenean barred-grid band: unchecked count must not exceed barred_max_unch/2.
eval_word_rule(checked_band, CS, W, Bits, result(checked_band, Sev, Detail)) :-
    get_dict(len, W, L), bits_checked_count(Bits, C),
    Unch is L - C, barred_max_unch(L, Max),
    ( Unch =< Max -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail),
                       "~w unchecked of ~w (Ximenean max ~w at this length)", [Unch, L, Max]) ).
eval_word_rule(max_unch_run, CS, _W, Bits, result(max_unch_run, Sev, Detail)) :-
    bits_max_unch_run(Bits, R),
    ( R =< 2 -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail), "unchecked run of ~w (max 2)", [R]) ).
eval_word_rule(double_unch_ends, CS, _W, Bits, result(double_unch_ends, Sev, Detail)) :-
    ( double_unch_end(Bits, Which)
    -> Sev = CS, format(atom(Detail), "double-unchecked at the ~w", [Which])
    ; Sev = pass, Detail = null ).
eval_word_rule(odd_even, CS, _W, Bits, result(odd_even, Sev, Detail)) :-
    ( lopsided_parity(Bits)
    -> Sev = CS, Detail = "checked cells all share one parity (odd/even imbalance)"
    ; Sev = pass, Detail = null ).

% Maximum unchecked letters per entry length for the barred-ximenean band
% (primary-sourced: Ximenes' 1966 "Ximenes on the Art of the Crossword" + Azed
% slip conventions): none in a 3, one in 4-5, two in 6-7, three in 8 (Azed
% self-limits to 2), and no more than a third (floor L/3) in entries of 9+.
barred_max_unch(L, 0) :- L =< 3, !.
barred_max_unch(4, 1) :- !.
barred_max_unch(5, 1) :- !.
barred_max_unch(6, 2) :- !.
barred_max_unch(7, 2) :- !.
barred_max_unch(8, 3) :- !.
barred_max_unch(L, M) :- L >= 9, M is L // 3.

% Bit-pattern rule helpers over a word's checked bitmap (metrics.pl builds it).
double_unch_end(Bits, start) :- Bits = [0, 0|_], !.
double_unch_end(Bits, end)   :- reverse(Bits, [0, 0|_]).

lopsided_parity(Bits) :-
    findall(I, nth1(I, Bits, 1), Checked),
    Checked = [_, _|_],
    ( forall(member(I, Checked), 1 is I mod 2)
    ; forall(member(I, Checked), 0 is I mod 2) ).


% --- grid rule evaluators ----------------------------------------------------
eval_grid_rule(connectivity, CS, Filled, GridLen, result(connectivity, Sev, Detail)) :-
    ( connected(Filled, GridLen) -> Sev = pass, Detail = null
    ; Sev = CS, Detail = "filled cells are not all connected" ).
eval_grid_rule(symmetry, CS, Filled, GridLen, result(symmetry, Sev, Detail)) :-
    symmetry_deficit(Filled, GridLen, Deficit),
    ( Deficit =:= 0 -> Sev = pass, Detail = null
    ; Sev = CS, format(atom(Detail), "180-rotational asymmetry at ~w cell(s)", [Deficit]) ).

% All filled cells reachable from the first via 4-adjacency through filled cells.
connected([], _GridLen) :- !.
connected(Filled, GridLen) :-
    Filled = [Start|_],
    % Filled is already an ordset (filled_cells/2 sorts it), so use it directly (P16).
    reach([Start], Filled, GridLen, [Start], Reached),
    length(Filled, N), length(Reached, N).

reach([], _FilledSet, _GridLen, Acc, Acc).
reach([Cell|Q], FilledSet, GridLen, Acc, Reached) :-
    findall(Nb,
            ( cell_neighbour(Cell, GridLen, Nb),
              ord_memberchk(Nb, FilledSet), \+ ord_memberchk(Nb, Acc) ),
            NewNbs0),
    sort(NewNbs0, NewNbs),
    ord_union(Acc, NewNbs, Acc1),
    append(Q, NewNbs, Q1),
    reach(Q1, FilledSet, GridLen, Acc1, Reached).

cell_neighbour(Cell, GridLen, Nb) :-
    cell_rc(Cell, GridLen, R, C),
    Last is GridLen - 1,
    ( R > 0,    Nb is Cell - GridLen
    ; R < Last, Nb is Cell + GridLen
    ; C > 0,    Nb is Cell - 1
    ; C < Last, Nb is Cell + 1 ).

% Count of filled cells whose 180-rotation is NOT filled (0 iff symmetric).
symmetry_deficit(Filled, GridLen, Deficit) :-
    % Filled is already an ordset (filled_cells/2 sorts it), so use it directly (P16).
    aggregate_all(count,
                  ( member(Cell, Filled), rot180(Cell, GridLen, Cell2),
                    \+ ord_memberchk(Cell2, Filled) ),
                  Deficit).

rot180(Cell, GridLen, Cell2) :-
    cell_rc(Cell, GridLen, R, C),
    R2 is GridLen - 1 - R, C2 is GridLen - 1 - C,
    Cell2 is R2 * GridLen + C2 + 1.


% --- verdict / summary -------------------------------------------------------
all_severities(WordReports, GridReports, Sevs) :-
    findall(S,
            ( member(WR, WordReports), get_dict(results, WR, Rs),
              member(R, Rs), get_dict(severity, R, S) ),
            WS),
    findall(S, ( member(GR, GridReports), get_dict(severity, GR, S) ), GS),
    append(WS, GS, Sevs).

worst_severity(Sevs, Verdict) :-
    ( memberchk('FAIL', Sevs) -> Verdict = 'FAIL'
    ; memberchk('WARN', Sevs) -> Verdict = 'WARN'
    ; Verdict = 'PASS' ).

tally(Sevs, Pass, Warn, Fail) :-
    include(==('PASS'), Sevs, Ps), length(Ps, Pass),
    include(==('WARN'), Sevs, Ws), length(Ws, Warn),
    include(==('FAIL'), Sevs, Fs), length(Fs, Fail).


% --- entry point -------------------------------------------------------------
% Emit the report on stdout (always, regardless of verdict - AC-LINT-1) and bind
% Verdict for the caller's exit code (AC-LINT-2). Malformed input throws.
lint_solve(File, Profile, AllowAsym, Verdict) :-
    lint_load(File, GridLen, PlacedWords),
    lint_run(PlacedWords, GridLen, Profile, AllowAsym, Report),
    current_output(Out),
    json_write_dict(Out, Report),
    nl(Out),
    get_dict(verdict, Report, Verdict).


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(lint_no_grid_length) -->
    [ 'lint: layout needs an integer "gridLength"' ].
prolog:error_message(lint_no_words_array) -->
    [ 'lint: layout needs a "words" array' ].
prolog:error_message(lint_invalid_word(Entry)) -->
    [ 'lint: every word needs a string "answer" (offending entry: ~q)'-[Entry] ].
prolog:error_message(lint_invalid_direction(A)) -->
    [ 'lint: word ~q needs a "direction" of "across" or "down"'-[A] ].
prolog:error_message(lint_no_cells(A)) -->
    [ 'lint: word ~q needs a non-empty "cells" array'-[A] ].
prolog:error_message(lint_invalid_cell(A, Pair)) -->
    [ 'lint: word ~q has an invalid cell ~q (need [row,col] within the grid)'-[A, Pair] ].
