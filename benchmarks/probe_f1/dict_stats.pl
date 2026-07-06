#!/usr/bin/env swipl
% benchmarks/probe_f1/dict_stats.pl - P-F1 probe: per-length dictionary stats +
% per-grid wasted-length fractions (question A, slot-length-filtered index
% sizing). PROBE-ONLY.
%
% For a dictionary: words per length, index keys per length, triples per length
% (triples = words*len - the unit of build_index findall/keysort/group work).
% For each bench grid: the slot-length set (from the mask), and the fraction of
% triples/keys/words on lengths the grid CANNOT use.
%
% Usage: swipl -q benchmarks/probe_f1/dict_stats.pl -- fixtures/dict/enable1.txt

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(assoc)).
:- use_module(library(pairs)).
:- use_module(library(ordsets)).

:- dynamic repo_root/1.
:- prolog_load_context(directory, ProbeDir),
   absolute_file_name('../..', RepoRoot,
                      [ relative_to(ProbeDir), file_type(directory), access(read) ]),
   asserta(repo_root(RepoRoot)),
   directory_file_path(RepoRoot, 'load.pl', Load),
   consult(Load).

:- initialization(main, main).

grid_file('fixtures/fill_grid_04a.json').
grid_file('fixtures/fill_grid_05a.json').
grid_file('fixtures/fill_grid_09a.json').
grid_file('fixtures/fill_grid_11a.json').
grid_file('fixtures/fill_grid_13a.json').
grid_file('fixtures/fill_grid_15a.json').
grid_file('fixtures/fill_grid_17a.json').
grid_file('fixtures/fill_grid_21a.json').

main :-
    current_prolog_flag(argv, [DictRel|_]),
    repo_root(Root),
    directory_file_path(Root, DictRel, File),
    crosswordsmith_fill:load_dict(File, DictByLen, Index),
    % per-length: words, triples(=words*len), index keys
    assoc_to_list(DictByLen, LenBuckets),
    findall(Len-stat(NW, TR), ( member(Len-Ws, LenBuckets), length(Ws, NW),
                                TR is NW*Len ), LenStats),
    assoc_to_keys(Index, IKeys),
    findall(L, member(k(L,_,_), IKeys), ILens),
    msort(ILens, ILensSorted),
    clumped(ILensSorted, KeysPerLen),
    format("len,words,triples,index_keys~n"),
    forall(member(Len-stat(NW, TR), LenStats),
           ( ( memberchk(Len-NK, KeysPerLen) -> true ; NK = 0 ),
             format("~d,~d,~d,~d~n", [Len, NW, TR, NK]) )),
    findall(TR, member(_-stat(_, TR), LenStats), TRs), sum_list(TRs, TotTR),
    findall(NW, member(_-stat(NW, _), LenStats), NWs), sum_list(NWs, TotNW),
    length(IKeys, TotNK),
    format("TOTAL words=~d triples=~d keys=~d~n", [TotNW, TotTR, TotNK]),
    % per-grid waste
    format("~ngrid,slot_lens,pct_words_unusable,pct_triples_unusable,pct_keys_unusable~n"),
    forall(grid_file(GRel),
           ( directory_file_path(Root, GRel, GFile),
             crosswordsmith_fill:fill_grid(GFile, _Size, Slots, _CV),
             findall(SL, ( member(slot(_,_,_,Vars), Slots), length(Vars, SL) ), SLs),
             sort(SLs, SlotLens),
             waste(LenStats, KeysPerLen, SlotLens, TotNW, TotTR, TotNK, PW, PT, PK),
             format("~w,~w,~2f,~2f,~2f~n", [GRel, SlotLens, PW, PT, PK]) )).

waste(LenStats, KeysPerLen, SlotLens, TotNW, TotTR, TotNK, PW, PT, PK) :-
    findall(NW-TR, ( member(Len-stat(NW, TR), LenStats),
                     \+ memberchk(Len, SlotLens) ), Wasted),
    findall(NW, member(NW-_, Wasted), WNW), sum_list(WNW, WastedW),
    findall(TR, member(_-TR, Wasted), WTR), sum_list(WTR, WastedT),
    findall(NK, ( member(Len-NK, KeysPerLen), \+ memberchk(Len, SlotLens) ), WNK),
    sum_list(WNK, WastedK),
    PW is 100.0*WastedW/TotNW,
    PT is 100.0*WastedT/TotTR,
    PK is 100.0*WastedK/TotNK.
