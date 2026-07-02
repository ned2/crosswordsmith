% load.pl - single owner of the project load order (design-spec §4).
%
% Consulting this file loads the whole crosswordsmith implementation. The CLI
% driver, tests/run_tests.pl, and the benchmark harnesses all load the project
% through this file, so the load order lives in exactly one place. It is
% deliberately side-effect free beyond loading: no initialization directives
% and no flag changes, so consulting it from any harness is safe.
%
% The `crosswordsmith` file-search alias points at prolog/crosswordsmith/,
% resolved relative to this file's own directory (prolog_load_context/2), so
% loads are cwd-independent.
:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, 'prolog/crosswordsmith', LibDir),
   (   user:file_search_path(crosswordsmith, LibDir)
   ->  true
   ;   assertz(user:file_search_path(crosswordsmith, LibDir))
   ).

% Known-good load order (the driver's historical order): arrange chain-loads
% core.pl (which chain-loads metrics.pl); lint must precede stockgrid
% (stockgrid calls lint_run/5); fill last (uses stockgrid + arrange + metrics).
% lint/export/stockgrid/fill perform no project loads of their own.
% Module-ized files are use_module'd from this (user) context so their exports
% land in `user`, where the not-yet-module-ized plain files still resolve them
% (Phase-4 bridge; see the migration plan).
:- ensure_loaded(crosswordsmith(arrange)).
:- ensure_loaded(crosswordsmith(lint)).
:- use_module(crosswordsmith(export)).
:- ensure_loaded(crosswordsmith(stockgrid)).
:- ensure_loaded(crosswordsmith(fill)).
