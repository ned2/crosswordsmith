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

% Known-good load order. Module files are use_module'd from this (user)
% context so their exports land in `user`; the still-plain core.pl resolves
% them there, and the modules' own calls to core predicates resolve via
% inheritance from `user` (Phase-4 bridge; see the migration plan).
% Constraints: metrics before core (core's transitional chain-load of
% metrics.pl then no-ops); core (plain, into `user`) before the arrange
% module, which relies on `user` inheritance for core's predicates until
% Phase 4.7; lint before stockgrid (stockgrid imports lint_run/5); fill last
% (imports stockgrid, metrics).
:- use_module(crosswordsmith(metrics)).
:- ensure_loaded(crosswordsmith(core)).
:- use_module(crosswordsmith(arrange)).
:- use_module(crosswordsmith(lint)).
:- use_module(crosswordsmith(export)).
:- use_module(crosswordsmith(stockgrid)).
:- use_module(crosswordsmith(fill)).
