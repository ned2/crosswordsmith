% benchmarks/probe_fh2/build_probe.pl - F-H2 v2 artifact build cost + size.
% In-process (excludes SWI startup): time load_dict (unchanged), build_masks (the
% F-H2 construction tax), and fill_save_index in BOTH v2 flavours (default =
% no masks; masks(true) = the --masks opt-in). Report both file sizes.
%   swipl -q -g "consult('benchmarks/probe_fh2/build_probe.pl'), run, halt" </dev/null
:- set_prolog_flag(verbose, silent).
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', RepoRoot, [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(RepoRoot, 'load.pl', LoadFile), consult(LoadFile),
   nb_setval(repo, RepoRoot).
:- use_module(library(apply)).

scale(enable_10k, 'fixtures/dict/enable_10k.txt', '/tmp/claude-1000/f-h2/enable_10k_v2def.idx', '/tmp/claude-1000/f-h2/enable_10k_v2.idx').
scale(enable_50k, 'fixtures/dict/enable_50k.txt', '/tmp/claude-1000/f-h2/enable_50k_v2def.idx', '/tmp/claude-1000/f-h2/enable_50k_v2.idx').
scale(enable1,    'fixtures/dict/enable1.txt',    '/tmp/claude-1000/f-h2/enable1_v2def.idx',    '/tmp/claude-1000/f-h2/enable1_v2.idx').

t(Goal, Wall) :- call_time(Goal, T), get_dict(wall, T, Wall).

one(Name, Rel, OutDef, OutMask) :-
    nb_getval(repo, Root), directory_file_path(Root, Rel, Dict),
    t(crosswordsmith_fill:load_dict(Dict, _DBL0, Idx0), Wload),
    t(crosswordsmith_fill:build_masks(Idx0, _MA), Wmask),
    % warm each flavour once, then time
    crosswordsmith_fill:fill_save_index(Dict, OutDef),
    t(crosswordsmith_fill:fill_save_index(Dict, OutDef), Wdef),
    crosswordsmith_fill:fill_save_index(Dict, OutMask, [masks(true)]),
    t(crosswordsmith_fill:fill_save_index(Dict, OutMask, [masks(true)]), Wm),
    size_file(OutDef, BDef), size_file(OutMask, BMask),
    format('~w: load=~3fs  save_v2_default=~3fs (~D bytes)  build_masks(TAX)=~3fs  save_v2_masks=~3fs (~D bytes)~n',
           [Name, Wload, Wdef, BDef, Wmask, Wm, BMask]).

run :-
    format('# F-H2 v2 artifact build cost + size, default vs masks (in-process, warm)~n'),
    forall(scale(N,R,OD,OM), one(N,R,OD,OM)).
