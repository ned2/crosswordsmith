% benchmarks/probe_fh2/phase_b.pl - F-H2 Phase B controlled kernel ablation.
%
% Load a v2 artifact once (DictByLen + Index + masks). Then time the fill(search)
% phase with the ORDSET kernel (fill_search .. none - byte-identical to the
% pre-F-H2 / v1 engine) vs the BIGNUM kernel (fill_search .. masks(MA)). Same
% load, same dict structures, same fresh slots per sample: the ONLY difference is
% the count kernel, so the wall delta is F-H2's fill-phase win, cleanly isolated
% from load (the honest artifact-mode wall metric). Warm, median of reps.
%
%   swipl -q benchmarks/probe_fh2/phase_b.pl -g 'fh2b_run, halt' </dev/null

:- set_prolog_flag(verbose, silent).
:- dynamic repo_root/1.
:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(Dir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', LoadFile),
   consult(LoadFile),
   directory_file_path(RepoRoot, 'benchmarks/fill_subjects.pl', Subjects),
   use_module(Subjects).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(pairs)).

abspath(Rel, Abs) :- repo_root(Root), directory_file_path(Root, Rel, Abs).

median(Xs, M) :-
    msort(Xs, S), length(S, L), L > 0, Mid is L // 2,
    ( 1 is L mod 2 -> nth0(Mid, S, M)
    ; I1 is Mid-1, nth0(I1,S,A), nth0(Mid,S,B), M is (A+B)/2 ).

reps(9).

rung(g09_full, 'fixtures/fill_grid_09a.json', '/tmp/claude-1000/f-h2/enable1_v2.idx',    none).
rung(g17_50k,  'fixtures/fill_grid_17a.json', '/tmp/claude-1000/f-h2/enable_50k_v2.idx', none).
rung(g21_full, 'fixtures/fill_grid_21a.json', '/tmp/claude-1000/f-h2/enable1_v2.idx',    none).

sol_words(Slots, Ws) :- findall(A, ( member(slot(_,_,_,V), Slots), atom_chars(A,V) ), Ws).

% One timed fill-phase run with the given Masks context, on FRESH slots.
run(GridFile, Seeds, DBL, Idx, Masks, Wall, Inf, Words) :-
    build_search_slots(GridFile, Seeds, slots(Slots, _)),
    call_time(once(crosswordsmith_fill:fill_search(Slots, DBL, Idx, Masks, [])), T),
    get_dict(wall, T, Wall), get_dict(inferences, T, Inf),
    sol_words(Slots, Words).

measure(Rung) :-
    rung(Rung, GridRel, Idx2, Seeds), abspath(GridRel, GridFile),
    crosswordsmith_fill:fill_load_index(Idx2, none, DBL, Idx, Masks),
    ( Masks = masks(_) -> true ; format('~w: NO MASKS in artifact!~n',[Rung]), fail ),
    reps(R),
    % warm both kernels
    run(GridFile, Seeds, DBL, Idx, none,  _, _, W0),
    run(GridFile, Seeds, DBL, Idx, Masks, _, _, W1),
    ( W0 == W1 -> Eq = 'IDENTICAL' ; Eq = 'DIVERGED' ),
    % ordset (v1) reps
    findall(Wc-Ic, ( between(1,R,_), run(GridFile,Seeds,DBL,Idx,none,Wc,Ic,_) ), OrdPairs),
    pairs_keys_values(OrdPairs, OrdW, OrdI), median(OrdW, OrdMed), OrdI=[OrdInf|_],
    % bignum (v2) reps
    findall(Wm-Im, ( between(1,R,_), run(GridFile,Seeds,DBL,Idx,Masks,Wm,Im,_) ), BigPairs),
    pairs_keys_values(BigPairs, BigW, BigI), median(BigW, BigMed), BigI=[BigInf|_],
    Delta is 100*(OrdMed-BigMed)/OrdMed,
    format('~n=== ~w (median of ~w warm reps) ===~n', [Rung, R]),
    format('  equivalence (ordset vs bignum filled grid): ~w~n', [Eq]),
    format('  ordset kernel (v1/none): fill wall=~4f s   search_inf=~D~n', [OrdMed, OrdInf]),
    format('  bignum kernel (v2/masks): fill wall=~4f s   search_inf=~D~n', [BigMed, BigInf]),
    format('  fill-phase WALL improvement: ~2f%%   (inf collapse ~2fx - NOT the metric)~n',
           [Delta, OrdInf/max(1,BigInf)]).

fh2b_run :-
    format('# F-H2 Phase B kernel ablation (artifact mode, load excluded)~n'),
    forall(rung(R,_,_,_), measure(R)).
