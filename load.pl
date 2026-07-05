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

% Every implementation file is a module with explicit imports, so inter-file
% dependencies no longer ride on this order — use_module'ing the top-level
% modules here (a) pulls in the whole dependency closure and (b) lands every
% export in `user`, where the CLI driver, the .plt suites, and the benchmark
% harnesses (all compiled into `user`) resolve them unqualified.
:- use_module(crosswordsmith(core)).
:- use_module(crosswordsmith(metrics)).
:- use_module(crosswordsmith(arrange)).
:- use_module(crosswordsmith(lint)).
:- use_module(crosswordsmith(export)).
:- use_module(crosswordsmith(stockgrid)).
:- use_module(crosswordsmith(fill)).
:- use_module(crosswordsmith(browser)).
