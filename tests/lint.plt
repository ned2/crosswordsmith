% tests/lint.plt - plunit suite for lint.pl (Flavour-B grid validator, §8.1).
%
% Assumes core.pl, metrics.pl (via core.pl) and lint.pl are consulted
% by the runner (tests/run_tests.pl) before this file. Run via `make test`.
%
% Covers the loader, the rule evaluators, the per-profile rule/severity sets
% (AC-LINT-4), the symmetry override (AC-LINT-3), and the report shape +
% verdict (AC-LINT-1/2). lint_run/5 is exercised directly on small hand-built
% layouts so each severity is checkable; the CLI exit-code mapping and the
% byte-exact report live in the golden (tests/golden/lint_bundled_17_toc.json).

:- use_module(library(plunit)).
% lists: test-body helpers; explicit so the suite also runs under
% autoload(false) (P11/C5).
:- use_module(library(lists)).

% A tiny controlled layout on a 5x5 grid: across CAT (cells 1,2,3) crossing down
% COT (cells 1,6,11) at cell 1. Each word is checked at only its shared cell, so
% neither is half-checked; each has an unchecked run of 2 ending the word.
mini_layout([ pw('CAT', _, [1,2,3],  across, 3, _, _, 1),
              pw('COT', _, [1,6,11], down,   3, _, _, 1) ]).

% A fully-interlocked 2x2 mesh (every cell checked); also symmetric + connected.
dense_mesh([ pw('AB', _, [1,2], across, 2, _, _, 1),
             pw('CD', _, [3,4], across, 2, _, _, 2),
             pw('AC', _, [1,3], down,   2, _, _, 1),
             pw('BD', _, [2,4], down,   2, _, _, 2) ]).

grid_rule_sev(Report, Rule, Sev) :-
    get_dict(grid, Report, Grid),
    once(( member(GR, Grid), get_dict(rule, GR, Rule), get_dict(severity, GR, Sev) )).


:- begin_tests(lint).

% --- loader ------------------------------------------------------------------

% A canonical words[] entry parses to a placed word: [r,c] cells -> 1-based
% numbers, direction -> atom, len from the cell count.
test(lint_loads_canonical_word) :-
    crosswordsmith_lint:lint_dict_layout(
        _{gridLength:5, words:[
            _{number:1, direction:"across", answer:"CAT", cells:[[0,0],[0,1],[0,2]]}]},
        5, [W]),
    pw_answer(W, 'CAT'), pw_dir(W, across),
    pw_cells(W, [1,2,3]), pw_len(W, 3), pw_num(W, 1).

test(lint_layout_no_grid_length_throws, [throws(error(lint_no_grid_length, _))]) :-
    crosswordsmith_lint:lint_dict_layout(_{words:[]}, _, _).

% --- report shape (AC-LINT-1) ------------------------------------------------

test(lint_report_shape) :-
    mini_layout(L),
    lint_run(L, 5, toc, false, R),
    get_dict(verdict, R, V), memberchk(V, ['PASS','WARN','FAIL']),
    get_dict(summary, R, Sum),
    get_dict(pass, Sum, _), get_dict(warn, Sum, _), get_dict(fail, Sum, _),
    get_dict(words, R, Ws), length(Ws, 2),
    once(( member(W, Ws), get_dict(results, W, Rs), Rs = [_|_] )).

% --- verdict / severity (AC-LINT-2) ------------------------------------------

% toc is advisory-only: never FAIL, so the mini layout verdicts WARN.
test(lint_verdict_warn_under_toc) :-
    mini_layout(L),
    lint_run(L, 5, toc, false, R),
    get_dict(verdict, R, 'WARN').

% Under blocked-uk the same layout FAILs (sub-half checking + asymmetry).
test(lint_verdict_fail_under_blocked_uk) :-
    mini_layout(L),
    lint_run(L, 5, 'blocked-uk', false, R),
    get_dict(verdict, R, 'FAIL'),
    % the checked-fraction rule is the FAIL for both words
    get_dict(words, R, Ws),
    forall( member(W, Ws),
            once(( get_dict(results, W, Rs), member(Res, Rs),
                   get_dict(rule, Res, checked_half), get_dict(severity, Res, 'FAIL') )) ).

% A fully-checked mesh passes checked_full under american (every cell checked).
test(lint_checked_full_passes_on_dense_mesh) :-
    dense_mesh(M),
    lint_run(M, 2, american, false, R),
    get_dict(words, R, Ws),
    forall( member(W, Ws),
            once(( get_dict(results, W, Rs), member(Res, Rs),
                   get_dict(rule, Res, checked_full), get_dict(severity, Res, 'PASS') )) ).

% --- per-profile rule sets (AC-LINT-4) ---------------------------------------

% american applies EXACTLY min_length + checked_full per word (its documented
% per-word set), nothing else.
test(lint_american_applies_exact_per_word_ruleset) :-
    mini_layout(L),
    lint_run(L, 5, american, false, R),
    get_dict(words, R, Ws),
    once(( member(W, Ws), get_dict(results, W, Rs),
           findall(Rule, ( member(Res, Rs), get_dict(rule, Res, Rule) ), Rules),
           Rules == [min_length, checked_full] )).

% --- connectivity ------------------------------------------------------------

% Two non-touching words on one grid form two components -> connectivity trips.
test(lint_connectivity_detects_split) :-
    Split = [ pw('CAT', _, [1,2,3], across, 3, _, _, 1),
              pw('DOG', _, [13,14,15], across, 3, _, _, 2) ],
    lint_run(Split, 5, toc, false, R),
    grid_rule_sev(R, connectivity, 'WARN').

% A single connected word passes connectivity.
test(lint_connectivity_passes_connected) :-
    One = [ pw('CAT', _, [4,5,6], across, 3, _, _, 1) ],
    lint_run(One, 3, toc, false, R),
    grid_rule_sev(R, connectivity, 'PASS').

% --- symmetry + the --allow-asymmetry override (AC-LINT-3) -------------------

% A centred 3x3 row {4,5,6} is 180-symmetric -> PASS even under blocked-uk.
test(lint_symmetry_passes_on_symmetric_pattern) :-
    Sym = [ pw('CAT', _, [4,5,6], across, 3, _, _, 1) ],
    lint_run(Sym, 3, 'blocked-uk', false, R),
    grid_rule_sev(R, symmetry, 'PASS').

% The top row {1,2,3} is asymmetric -> FAIL under blocked-uk...
test(lint_symmetry_fails_asymmetric_under_blocked_uk) :-
    Asym = [ pw('CAT', _, [1,2,3], across, 3, _, _, 1) ],
    lint_run(Asym, 3, 'blocked-uk', false, R),
    grid_rule_sev(R, symmetry, 'FAIL').

% ...never hard-FAILs under --allow-asymmetry (downgraded to WARN)...
test(lint_symmetry_allow_asymmetry_downgrades_to_warn) :-
    Asym = [ pw('CAT', _, [1,2,3], across, 3, _, _, 1) ],
    lint_run(Asym, 3, 'blocked-uk', true, R),
    grid_rule_sev(R, symmetry, 'WARN').

% ...and never FAILs at all under a profile that does not enforce it (toc).
test(lint_symmetry_advisory_under_toc) :-
    Asym = [ pw('CAT', _, [1,2,3], across, 3, _, _, 1) ],
    lint_run(Asym, 3, toc, false, R),
    grid_rule_sev(R, symmetry, 'WARN').

% --- profile registry --------------------------------------------------------

test(lint_known_profiles) :-
    forall(member(P, [toc, 'blocked-uk', american, 'barred-ximenean']),
           lint_known_profile(P)),
    \+ lint_known_profile(nope).

% --- barred-ximenean: the Ximenean per-length unchecked-letter band (OD-7) ----

% The primary-sourced table: none in a 3, one in 4-5, two in 6-7, three in 8,
% a third (floor L/3) in 9+.
test(barred_band_table) :-
    crosswordsmith_lint:barred_max_unch(3, 0), crosswordsmith_lint:barred_max_unch(4, 1), crosswordsmith_lint:barred_max_unch(5, 1),
    crosswordsmith_lint:barred_max_unch(6, 2), crosswordsmith_lint:barred_max_unch(7, 2), crosswordsmith_lint:barred_max_unch(8, 3),
    crosswordsmith_lint:barred_max_unch(9, 3), crosswordsmith_lint:barred_max_unch(12, 4).

% A 5-letter entry checked at only 1 cell has 4 unches > the max of 1 -> FAIL.
% (P4: eval_word_rule/5 now takes the word's precomputed checked bitmap; build it
% via the same layout_dir_cells/word_checked_bitmap path lint_run uses.)
test(barred_band_fails_underchecked_five) :-
    W = pw('ABCDE', _, [1,2,3,4,5], across, 5, _, _, 1),
    D = pw('XYZ', _, [3,20,37], down, 3, _, _, 2),
    layout_dir_cells([W, D], DirCells),
    word_checked_bitmap(W, DirCells, Bits),
    crosswordsmith_lint:eval_word_rule(checked_band, fail, W, Bits, result(checked_band, Sev, _)),
    Sev == fail.

% A fully-checked entry has 0 unches, within any band -> PASS.
test(barred_band_passes_fully_checked) :-
    dense_mesh(M),
    layout_dir_cells(M, DirCells),
    once(( member(W, M),
           word_checked_bitmap(W, DirCells, Bits),
           crosswordsmith_lint:eval_word_rule(checked_band, fail, W, Bits, result(checked_band, pass, _)) )).

% The profile applies the band per word and RELAXES symmetry to advisory (WARN).
test(barred_profile_applies_band_relaxed_symmetry) :-
    mini_layout(L),
    lint_run(L, 5, 'barred-ximenean', false, R),
    get_dict(words, R, Ws),
    once(( member(WR, Ws), get_dict(results, WR, Rs),
           member(Res, Rs), get_dict(rule, Res, checked_band) )),
    grid_rule_sev(R, symmetry, 'WARN').

:- end_tests(lint).
