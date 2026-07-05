% benchmarks/probe_fh2/build_probe.pl - F-H2 v2 artifact build cost + size.
% In-process (excludes SWI startup): time load_dict (unchanged), build_masks (the
% F-H2 construction tax), and the whole fill_save_index (v2). Report file size.
%   swipl -q benchmarks/probe_fh2/build_probe.pl -g 'run, halt' </dev/null
:- set_prolog_flag(verbose, silent).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', RepoRoot, [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(RepoRoot, 'load.pl', LoadFile), consult(LoadFile),
   nb_setval(repo, RepoRoot).
:- use_module(library(apply)).

scale(enable_10k, 'fixtures/dict/enable_10k.txt', '/tmp/claude-1000/f-h2/enable_10k_v2.idx').
scale(enable_50k, 'fixtures/dict/enable_50k.txt', '/tmp/claude-1000/f-h2/enable_50k_v2.idx').
scale(enable1,    'fixtures/dict/enable1.txt',    '/tmp/claude-1000/f-h2/enable1_v2.idx').

t(Goal, Wall) :- call_time(Goal, T), get_dict(wall, T, Wall).

one(Name, Rel, Out) :-
    nb_getval(repo, Root), directory_file_path(Root, Rel, Dict),
    t(crosswordsmith_fill:load_dict(Dict, _DBL0, Idx0), Wload),
    t(crosswordsmith_fill:build_masks(Idx0, _MA), Wmask),
    % warm, then time full save
    crosswordsmith_fill:fill_save_index(Dict, Out),
    t(crosswordsmith_fill:fill_save_index(Dict, Out), Wsave),
    size_file(Out, Bytes),
    format('~w: load=~3fs  build_masks(TAX)=~3fs  full_save_v2=~3fs  size=~D bytes~n',
           [Name, Wload, Wmask, Wsave, Bytes]).

run :-
    format('# F-H2 v2 artifact build cost + size (in-process, warm)~n'),
    forall(scale(N,R,O), one(N,R,O)).
