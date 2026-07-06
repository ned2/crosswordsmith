% F-H2 gate probe — real index-set distribution + construction tax.
% Loads the actual dict via the engine's OWN load_dict/3 (crosswordsmith:fill),
% then (a) reports the index-set size distribution + bucket sizes (grounds the
% kernel operand sizes), (b) times building a PARALLEL bignum bitset index from
% the existing ordsets (the load-time tax F-H2 would add).
% MEASUREMENT ONLY. Engine untouched. Requires load.pl consulted first
% (sets the crosswordsmith alias + loads crosswordsmith_fill).
:- use_module(library(assoc)).
:- use_module(library(apply)).

timed(Goal, WallSec, Inf) :-
    garbage_collect,
    statistics(inferences, I0), get_time(T0),
    call(Goal),
    get_time(T1), statistics(inferences, I1),
    WallSec is T1 - T0, Inf is I1 - I0.

% Build one bignum mask from a sorted ordset (the natural F-H2 construction).
ordset_to_mask(Set, Mask) :- foldl(set_bit, Set, 0, Mask).
set_bit(I, A, A1) :- A1 is A \/ (1<<I).

% Build the parallel bitset index: one bignum per index key.
build_bitset_index(Index, BIndex) :-
    assoc_to_list(Index, Pairs),
    maplist(key_ordset_to_mask, Pairs, BPairs),
    list_to_assoc(BPairs, BIndex).
key_ordset_to_mask(K-Set, K-Mask) :- ordset_to_mask(Set, Mask).

quantile(Sorted, Len, Frac, V) :-
    Idx0 is truncate(Frac*(Len-1)), nth0(Idx0, Sorted, V).

run(Runtime, File) :-
    timed(crosswordsmith_fill:load_dict(File, DictByLen, Index), WLoad, InfLoad),
    % bucket sizes (per length)
    assoc_to_list(DictByLen, LB),
    maplist([_-Ws, N]>>length(Ws, N), LB, BucketSizes0),
    max_list(BucketSizes0, MaxBucket), sum_list(BucketSizes0, TotWords),
    % index-set size distribution
    assoc_to_values(Index, Sets),
    length(Sets, NKeys),
    maplist(length, Sets, Sizes0),
    msort(Sizes0, SS), max_list(SS, MaxSet), sum_list(SS, SumSet),
    MeanSet is SumSet / NKeys,
    quantile(SS, NKeys, 0.25, Q1), quantile(SS, NKeys, 0.50, Med),
    quantile(SS, NKeys, 0.75, Q3), quantile(SS, NKeys, 0.90, P90),
    quantile(SS, NKeys, 0.99, P99),
    format("DIST ~w file=~w | load_inf=~w load_wall=~3fs | nwords=~w nkeys=~w maxbucket=~w | set q1=~w med=~w q3=~w p90=~w p99=~w max=~w mean=~2f~n",
           [Runtime, File, InfLoad, WLoad, TotWords, NKeys, MaxBucket, Q1, Med, Q3, P90, P99, MaxSet, MeanSet]),
    % construction tax: build parallel bignum index from existing ordsets
    timed(build_bitset_index(Index, BIndex), WBuild, InfBuild),
    ( InfLoad > 0 -> PctInf is 100.0*InfBuild/InfLoad ; PctInf = 0.0 ),
    ( WLoad  > 0 -> PctWall is 100.0*WBuild/WLoad ; PctWall = 0.0 ),
    % sanity: popcount of a mask == size of its ordset (byte-identity of the count)
    assoc_to_list(Index, [K1-S1|_]), get_assoc(K1, BIndex, M1),
    length(S1, LS1), PC1 is popcount(M1),
    ( LS1 =:= PC1 -> Ok = count_identical ; Ok = count_mismatch ),
    format("BUILD ~w file=~w | build_inf=~w (~2f%% of load_inf) build_wall=~3fs (~2f%% of load_wall) | check ~w (len=~w popcount=~w)~n",
           [Runtime, File, InfBuild, PctInf, WBuild, PctWall, Ok, LS1, PC1]).
