:- module(bench_paths,
          [ benchmark_path/2,
            repo_path/2
          ]).

:- use_module(library(filesex), [directory_file_path/3]).

%!  benchmark_path(+Relative:atom, -Path:atom) is det.
%
%   Resolve Relative against the directory containing this module. Resolution
%   is independent of the process working directory and stores no mutable root.
benchmark_path(Relative, Path) :-
    benchmark_directory(Directory),
    directory_file_path(Directory, Relative, Path).

%!  repo_path(+Relative:atom, -Path:atom) is det.
%
%   Resolve Relative against the repository root containing benchmarks/.
repo_path(Relative, Path) :-
    benchmark_directory(BenchmarkDirectory),
    file_directory_name(BenchmarkDirectory, RepositoryDirectory),
    directory_file_path(RepositoryDirectory, Relative, Path).

benchmark_directory(Directory) :-
    once(module_property(bench_paths, file(Source))),
    file_directory_name(Source, Directory).
