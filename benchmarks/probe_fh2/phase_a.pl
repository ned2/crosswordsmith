% benchmarks/probe_fh2/phase_a.pl - F-H2 Phase A re-attribution gate (MEASUREMENT).
%
% On the CURRENT post-F-H1 engine, in artifact-consuming mode, attribute the
% fill(search)-phase WALL and inferences into three buckets:
%   (1) counting kernel  = candidate_count/4 (index_set + ord_intersection + count)
%   (2) materialization  = candidates/4 (winner-path intersection + nth0_of walk)
%   (3) other            = the remainder (unify, select_min_count, recount scan, ...)
% F-H2 can only cheapen bucket (1) (masks replace candidate_count's ordset kernel;
% candidates/4 stays on ordsets because materialization needs the ordset indices).
%
% METHOD (P-F1 precedent, sharpened for the LCO-buried ord_intersection wall):
% a verbatim copy of the search hot path (fill_search/fill_search_inc + helpers)
% that calls the REAL crosswordsmith_fill:candidate_count/4 and candidates/4
% through call_time/2 accumulators. Counts/tree are the engine's own (the copy
% only reorders wrappers), verified by a term-identical filled-grid check against
% the clean engine. The clean fill-phase wall (W_real) is the denominator; the
% call_time accumulators (outside the wrapped goal) give each bucket's wall+inf.
%
% Instruments stay branch-only (P-F1/gate-probe convention). Run:
%   swipl -q benchmarks/probe_fh2/phase_a.pl -g 'fh2_run, halt' -- </dev/null

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

% Resolve a repo-relative fixture path to an absolute one (cwd-independent).
abspath(Rel, Abs) :- repo_root(Root), directory_file_path(Root, Rel, Abs).

% ---- accumulators -----------------------------------------------------------
reset_acc :-
    forall(member(K, [cc_wall, cc_inf, cc_n, cd_wall, cd_inf, cd_n]),
           nb_setval(K, 0)).

accrue(WKey, IKey, NKey, T) :-
    get_dict(wall, T, W), get_dict(inferences, T, I),
    nb_getval(WKey, W0), W1 is W0 + W, nb_setval(WKey, W1),
    nb_getval(IKey, I0), I1 is I0 + I, nb_setval(IKey, I1),
    nb_getval(NKey, N0), N1 is N0 + 1, nb_setval(NKey, N1).

% ---- instrumented wrappers over the REAL engine predicates ------------------
i_candidate_count(Vars, DictByLen, Index, Count) :-
    call_time(crosswordsmith_fill:candidate_count(Vars, DictByLen, Index, Count), T),
    accrue(cc_wall, cc_inf, cc_n, T).

i_candidates(Vars, DictByLen, Index, Cands) :-
    call_time(crosswordsmith_fill:candidates(Vars, DictByLen, Index, Cands), T),
    accrue(cd_wall, cd_inf, cd_n, T).

% ---- verbatim copy of the search hot path (fill.pl:238-312), with the two
%      count/materialize seams routed through the instrumented wrappers --------
i_fill_search(Slots, DictByLen, Index, Used) :-
    maplist(i_slot_with_count(DictByLen, Index), Slots, Counted),
    i_fill_search_inc(Counted, DictByLen, Index, Used).

i_slot_with_count(DictByLen, Index, Slot, cnt(Count, Slot)) :-
    Slot = slot(_, _, _, Vars),
    i_candidate_count(Vars, DictByLen, Index, Count).

i_fill_search_inc([], _DictByLen, _Index, _Used) :- !.
i_fill_search_inc(Counted, DictByLen, Index, Used) :-
    crosswordsmith_fill:select_min_count(Counted, cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    crosswordsmith_fill:newly_bound_cells(BestCells, BestVars, NewCells),
    i_candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    BestVars = Word,
    i_recount_crossing(Rest, NewCells, DictByLen, Index, Rest1),
    i_fill_search_inc(Rest1, DictByLen, Index, [Word|Used]).

i_recount_crossing([], _NewCells, _DictByLen, _Index, []).
i_recount_crossing([cnt(C0, Slot)|T], NewCells, DictByLen, Index, [cnt(C1, Slot)|T1]) :-
    Slot = slot(_, _, Cells, Vars),
    ( crosswordsmith_fill:shares_cell(Cells, NewCells)
    ->  i_candidate_count(Vars, DictByLen, Index, C1)
    ;   C1 = C0 ),
    i_recount_crossing(T, NewCells, DictByLen, Index, T1).

% ---- DECOMPOSITION of candidate_count: which sub-part does F-H2 replace? -----
% F-H2 (masks) replaces index_intersection + the final length; it does NOT
% replace slot_bucket's Bound findall (the mask still needs the bound positions)
% nor length(Vars)/get_assoc. So the ACHIEVABLE ceiling is the isx+len sub-share
% of candidate_count, not the whole bucket. accrue into isx / bnd / len.
reset_acc2 :-
    forall(member(K, [isx_wall,isx_inf,isx_n, bnd_wall,bnd_inf,bnd_n,
                      len_wall,len_inf,len_n]),
           nb_setval(K, 0)).

i2_slot_bucket(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    call_time(findall(P-V, ( nth0(P, Vars, V), nonvar(V) ), Bound), TB),
    accrue(bnd_wall, bnd_inf, bnd_n, TB),
    ( Bound == []
    ->  Sel = all
    ;   call_time(crosswordsmith_fill:index_intersection(Bound, Len, Index, Ix), TX),
        accrue(isx_wall, isx_inf, isx_n, TX),
        Sel = idx(Ix) ).

i2_candidate_count(Vars, DictByLen, Index, Count) :-
    i2_slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all
    ->  call_time(length(Words, Count), TL), accrue(len_wall, len_inf, len_n, TL)
    ;   Sel = idx(Ix),
        call_time(length(Ix, Count), TL), accrue(len_wall, len_inf, len_n, TL) ).

i2_fill_search(Slots, DictByLen, Index, Used) :-
    maplist(i2_slot_with_count(DictByLen, Index), Slots, Counted),
    i2_fill_search_inc(Counted, DictByLen, Index, Used).
i2_slot_with_count(DictByLen, Index, Slot, cnt(Count, Slot)) :-
    Slot = slot(_, _, _, Vars),
    i2_candidate_count(Vars, DictByLen, Index, Count).
i2_fill_search_inc([], _, _, _) :- !.
i2_fill_search_inc(Counted, DictByLen, Index, Used) :-
    crosswordsmith_fill:select_min_count(Counted, cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    crosswordsmith_fill:newly_bound_cells(BestCells, BestVars, NewCells),
    crosswordsmith_fill:candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands), \+ memberchk(Word, Used), BestVars = Word,
    i2_recount_crossing(Rest, NewCells, DictByLen, Index, Rest1),
    i2_fill_search_inc(Rest1, DictByLen, Index, [Word|Used]).
i2_recount_crossing([], _, _, _, []).
i2_recount_crossing([cnt(C0, Slot)|T], NewCells, DictByLen, Index, [cnt(C1, Slot)|T1]) :-
    Slot = slot(_, _, Cells, Vars),
    ( crosswordsmith_fill:shares_cell(Cells, NewCells)
    ->  i2_candidate_count(Vars, DictByLen, Index, C1)
    ;   C1 = C0 ),
    i2_recount_crossing(T, NewCells, DictByLen, Index, T1).

decomp_rung(Rung) :-
    rung(Rung, GridRel, IdxFile, Seeds), abspath(GridRel, GridFile),
    crosswordsmith_fill:fill_load_index(IdxFile, none, DictByLen, Index),
    reps(R),
    reset_acc2,
    build_search_slots(GridFile, Seeds, slots(WS, _)),      % warm
    once(i2_fill_search(WS, DictByLen, Index, [])),
    findall(Isx-Bnd-Len, ( between(1, R, _),
              reset_acc2,
              build_search_slots(GridFile, Seeds, slots(S, _)),
              once(i2_fill_search(S, DictByLen, Index, [])),
              nb_getval(isx_wall, Isx), nb_getval(bnd_wall, Bnd),
              nb_getval(len_wall, Len) ), Trips),
    findall(X, member(X-_-_, Trips), Isxs), median(Isxs, IsxW),
    findall(Y, member(_-Y-_, Trips), Bnds), median(Bnds, BndW),
    findall(Z, member(_-_-Z, Trips), Lens), median(Lens, LenW),
    nb_getval(isx_inf, IsxI), nb_getval(bnd_inf, BndI), nb_getval(len_inf, LenI),
    nb_getval(isx_n, IsxN),
    CcW is IsxW + BndW + LenW,       % ~ candidate_count wall (decomposed)
    ( CcW > 0 -> ReplPct is 100*(IsxW+LenW)/CcW, BndPct is 100*BndW/CcW ; ReplPct = 0, BndPct = 0),
    format('~n--- ~w candidate_count decomposition (median ~w reps) ---~n', [Rung, R]),
    format('  intersection (index_intersection): wall=~4f s  inf=~D  (n=~D)~n', [IsxW, IsxI, IsxN]),
    format('  count length (idx|bucket)        : wall=~4f s  inf=~D~n', [LenW, LenI]),
    format('  Bound findall (SURVIVES F-H2)    : wall=~4f s  inf=~D~n', [BndW, BndI]),
    format('  => F-H2-replaceable (isx+len) = ~2f%% of candidate_count wall; Bound survives = ~2f%%~n',
           [ReplPct, BndPct]).

decomp_run :-
    format('~n# candidate_count sub-decomposition (F-H2 ceiling)~n'),
    forall(rung(R,_,_,_), decomp_rung(R)).

% ---- rung table -------------------------------------------------------------
rung(g09_full, 'fixtures/fill_grid_09a.json', '/tmp/claude-1000/f-h2/enable1_v1.idx',    none).
rung(g17_50k,  'fixtures/fill_grid_17a.json', '/tmp/claude-1000/f-h2/enable_50k_v1.idx', none).
rung(g21_full, 'fixtures/fill_grid_21a.json', '/tmp/claude-1000/f-h2/enable1_v1.idx',    none).

% Solve words of a filled grid (for the term-identity equivalence check).
solution_words(Slots, Words) :-
    findall(A, ( member(slot(_,_,_,Vars), Slots), atom_chars(A, Vars) ), Words).

median(Xs, M) :-
    msort(Xs, Sorted), length(Sorted, L), L > 0,
    Mid is L // 2,
    ( 1 is L mod 2
    -> nth0(Mid, Sorted, M)
    ;  I1 is Mid - 1, nth0(I1, Sorted, A), nth0(Mid, Sorted, B), M is (A + B) / 2 ).

reps(5).

% One clean timed fill-phase run on FRESH slots -> wall.
clean_run(GridFile, Seeds, DictByLen, Index, Wall, Inf, Words) :-
    build_search_slots(GridFile, Seeds, slots(Slots, _All)),
    % none = the ordset kernel (Phase A measures the CURRENT counting seam); on the
    % post-F-H2 engine this is byte-identical to the pre-F-H2 fill_search/4.
    call_time(once(crosswordsmith_fill:fill_search(Slots, DictByLen, Index, none, [])), T),
    get_dict(wall, T, Wall), get_dict(inferences, T, Inf),
    solution_words(Slots, Words).

% One instrumented timed run on FRESH slots -> per-run bucket wall accruals.
instr_run(GridFile, Seeds, DictByLen, Index, CcW, CdW, Winstr, Words) :-
    reset_acc,
    build_search_slots(GridFile, Seeds, slots(Slots, _All)),
    call_time(once(i_fill_search(Slots, DictByLen, Index, [])), T),
    get_dict(wall, T, Winstr),
    nb_getval(cc_wall, CcW), nb_getval(cd_wall, CdW),
    solution_words(Slots, Words).

% ---- one rung: warm, median-of-reps, thermally consistent -------------------
measure_rung(Rung) :-
    rung(Rung, GridRel, IdxFile, Seeds),
    abspath(GridRel, GridFile),
    crosswordsmith_fill:fill_load_index(IdxFile, none, DictByLen, Index),
    reps(R),

    % WARM UP (2 clean + 1 instrumented, untimed use of result).
    clean_run(GridFile, Seeds, DictByLen, Index, _, _, CleanWords),
    clean_run(GridFile, Seeds, DictByLen, Index, _, _, _),
    instr_run(GridFile, Seeds, DictByLen, Index, _, _, _, InstrWords),
    ( CleanWords == InstrWords -> Equiv = 'IDENTICAL' ; Equiv = 'DIVERGED' ),

    % CLEAN reps -> median fill-phase wall (the denominator). Inf is deterministic.
    findall(W-I, ( between(1, R, _),
                   clean_run(GridFile, Seeds, DictByLen, Index, W, I, _) ),
            CleanPairs),
    pairs_keys_values(CleanPairs, CleanWalls, CleanInfs),
    median(CleanWalls, WrealMed),
    CleanInfs = [InfReal|_],

    % INSTRUMENTED reps -> median bucket walls. Inf/n deterministic (read last run).
    findall(cc(CcW,CdW), ( between(1, R, _),
              instr_run(GridFile, Seeds, DictByLen, Index, CcW, CdW, _, _) ), Accs),
    findall(C, member(cc(C,_), Accs), CcWalls),
    findall(D, member(cc(_,D), Accs), CdWalls),
    median(CcWalls, CcWMed), median(CdWalls, CdWMed),
    nb_getval(cc_inf, CcI0), nb_getval(cc_n, CcN),
    nb_getval(cd_inf, CdI0), nb_getval(cd_n, CdN),
    CcI is CcI0 + 2*CcN, CdI is CdI0 + 2*CdN,   % undo call_time's ~2-inf/call bias

    OtherW is WrealMed - CcWMed - CdWMed,
    CcWpct  is 100 * CcWMed / WrealMed,
    CdWpct  is 100 * CdWMed / WrealMed,
    OthWpct is 100 * OtherW / WrealMed,
    CcIpct  is 100 * CcI / InfReal,
    CdIpct  is 100 * CdI / InfReal,
    OthI    is InfReal - CcI - CdI,
    OthIpct is 100 * OthI / InfReal,
    ( CcWpct >= 20.0 -> Gate = 'PASS(>=20%)' ; Gate = 'below-20%' ),

    format('~n=== ~w  (median of ~w warm reps) ===~n', [Rung, R]),
    format('  equivalence: ~w~n', [Equiv]),
    format('  clean fill-phase wall (median): ~4f s   search_inf=~D~n', [WrealMed, InfReal]),
    format('  counting kernel candidate_count/4: n=~D  wall=~4f s (~2f%% wall)  inf=~D (~2f%% inf)~n',
           [CcN, CcWMed, CcWpct, CcI, CcIpct]),
    format('  materializ.     candidates/4     : n=~D  wall=~4f s (~2f%% wall)  inf=~D (~2f%% inf)~n',
           [CdN, CdWMed, CdWpct, CdI, CdIpct]),
    format('  other  unify/select/recount-scan :        wall=~4f s (~2f%% wall)  inf=~D (~2f%% inf)~n',
           [OtherW, OthWpct, OthI, OthIpct]),
    format('  GATE (20%% wall on counting kernel): ~w~n', [Gate]).

fh2_run :-
    format('# F-H2 Phase A re-attribution (post-F-H1, artifact mode)~n'),
    forall(rung(R,_,_,_), measure_rung(R)).
