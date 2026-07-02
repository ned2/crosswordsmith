% load.pl - single owner of the project load order (design-spec §4).
%
% Consulting this file loads the whole crosswordsmith implementation. The CLI
% driver, tests/run_tests.pl, and the benchmark harnesses all load the project
% through this file, so the load order lives in exactly one place. It is
% deliberately side-effect free beyond loading: no initialization directives
% and no flag changes, so consulting it from any harness is safe.
%
% The `crosswordsmith` file-search alias resolves relative to this file's own
% directory (prolog_load_context/2), so loads are cwd-independent. It currently
% points at the repository root; the source-structure migration
% (docs/source-structure-migration-plan.md Phase 1) repoints it at
% prolog/crosswordsmith/ when the implementation files move there.
:- prolog_load_context(directory, Dir),
   (   user:file_search_path(crosswordsmith, Dir)
   ->  true
   ;   assertz(user:file_search_path(crosswordsmith, Dir))
   ).

% Known-good load order (the driver's historical order): arrange chain-loads
% crossword.pl (which chain-loads quality.pl); lint must precede stockgrid
% (stockgrid calls lint_run/5); fill last (uses stockgrid + arrange + metrics).
% lint/export/stockgrid/fill perform no project loads of their own.
:- ensure_loaded(crosswordsmith(arrange)).
:- ensure_loaded(crosswordsmith(lint)).
:- ensure_loaded(crosswordsmith(export)).
:- ensure_loaded(crosswordsmith(stockgrid)).
:- ensure_loaded(crosswordsmith(fill)).
