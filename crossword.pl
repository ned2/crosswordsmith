#!/usr/bin/swipl

:- set_prolog_flag(verbose, silent).
% crossword.pl - migration-message shim ONLY. The implementation moved to
% prolog/crosswordsmith/ (core.pl et al.; see load.pl) and the CLI is the
% `crosswordsmith` script (design-spec §5.3). Running this file just points
% the user at the new CLI and exits non-zero; it loads nothing.
:- initialization(legacy_main, main).

legacy_main :-
    format(user_error,
"crossword.pl is now an internal library; the CLI is the `crosswordsmith` script.~n\c
  old: ./crossword.pl --input F <N> <loc>   ->  crosswordsmith arrange --strict --size-mode fixed --size <N> --input F~n\c
  old: ./crossword.pl --input F --quality   ->  crosswordsmith arrange --best-effort --size-mode max --input F~n\c
  old: ./crossword.pl --input F --all <N>   ->  crosswordsmith arrange --enumerate --size <N> --input F~n\c
Run `crosswordsmith` with no arguments for usage. (--shuffle is removed: output is deterministic.)~n",
           []),
    halt(1).
