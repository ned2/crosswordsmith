% F-H2 gate probe — kernel micro-bench (runs identically native + swipl-wasm/node).
% Models the fill counting kernel (fill.pl:148-184):
%   (i)  ordset path : ord_intersection/3 chain + length/2   -> per-slot count
%   (ii) bignum path : X is A /\ B ; C is popcount(X)         -> proposed F-H2
% MEASUREMENT ONLY. Not loaded by the engine.
:- use_module(library(ordsets)).

% ---- operand builders -------------------------------------------------------
% N distinct ints spread across universe U (density N/U). Sorted ascending.
make_list(N, U, L) :-
    ( N =< 0 -> L = []
    ; Stride is max(1, U // N), make_list_(0, N, Stride, U, L) ).
make_list_(_, 0, _, _, []) :- !.
make_list_(V, K, S, U, [V|T]) :-
    ( V < U -> true ; throw(overflow(V,U)) ),
    K1 is K-1, V1 is V+S, make_list_(V1, K1, S, U, T).

make_ordset(N, U, Set) :- make_list(N, U, L), list_to_ord_set(L, Set).
% second operand: same size, offset by half a stride -> ~50% overlap (cost is
% O(|A|+|B|) regardless of overlap; result size is realistic-ish).
make_ordset2(N, U, Set) :-
    ( N =< 0 -> Set = []
    ; Stride is max(1, U // N), Off is Stride // 2,
      make_list_(Off, N, Stride, U, L), list_to_ord_set(L, Set) ).

make_bignum(N, U, Mask) :-
    ( N =< 0 -> Mask = 0
    ; Stride is max(1, U // N), make_bignum_(0, N, Stride, U, 0, Mask) ).
make_bignum_(_, 0, _, _, M, M) :- !.
make_bignum_(V, K, S, U, M0, M) :-
    ( V < U -> true ; throw(overflow(V,U)) ),
    M1 is M0 \/ (1<<V), K1 is K-1, V1 is V+S, make_bignum_(V1, K1, S, U, M1, M).
make_bignum2(N, U, Mask) :-
    ( N =< 0 -> Mask = 0
    ; Stride is max(1, U // N), Off is Stride // 2,
      make_bignum_(Off, N, Stride, U, 0, Mask) ).

% ---- timed loops ------------------------------------------------------------
timed(Goal, WallSec, Inf) :-
    garbage_collect,
    statistics(inferences, I0), get_time(T0),
    call(Goal),
    get_time(T1), statistics(inferences, I1),
    WallSec is T1 - T0, Inf is I1 - I0.

loop_ord(0, _, _, Acc, Acc) :- !.
loop_ord(K, A, B, Acc0, Acc) :-
    ord_intersection(A, B, C), length(C, Len),
    Acc1 is Acc0 + Len, K1 is K-1, loop_ord(K1, A, B, Acc1, Acc).

loop_big(0, _, _, Acc, Acc) :- !.
loop_big(K, MA, MB, Acc0, Acc) :-
    X is MA /\ MB, P is popcount(X),
    Acc1 is Acc0 + P, K1 is K-1, loop_big(K1, MA, MB, Acc1, Acc).

% empty loop to subtract per-iteration bookkeeping inferences.
loop_noop(0, Acc, Acc) :- !.
loop_noop(K, Acc0, Acc) :- Acc1 is Acc0 + 0, K1 is K-1, loop_noop(K1, Acc1, Acc).

% ---- chain (foldl over M sets), models index_intersection over M bound cells -
loop_ord_chain(0, _, _, Acc, Acc) :- !.
loop_ord_chain(K, [S0|Ss], _, Acc0, Acc) :-
    foldl(ord_intersection_r, Ss, S0, C), length(C, Len),
    Acc1 is Acc0 + Len, K1 is K-1, loop_ord_chain(K1, [S0|Ss], _, Acc1, Acc).
ord_intersection_r(S, Acc, Acc1) :- ord_intersection(Acc, S, Acc1).

loop_big_chain(0, _, _, Acc, Acc) :- !.
loop_big_chain(K, [M0|Ms], _, Acc0, Acc) :-
    foldl(band_r, Ms, M0, X), P is popcount(X),
    Acc1 is Acc0 + P, K1 is K-1, loop_big_chain(K1, [M0|Ms], _, Acc1, Acc).
band_r(M, Acc, Acc1) :- Acc1 is Acc /\ M.

% ---- reps table (target ~0.5-2s wall per cell on native; wasm ~3x) ----------
reps(100,   1500000).
reps(1000,   400000).
reps(10000,   20000).
reps(24000,   10000).

per_op_us(WallSec, Reps, Us) :- Us is (WallSec / Reps) * 1.0e6.

run_one(Runtime, Size) :- reps(Size, Reps), run_one_reps(Runtime, Size, Reps).

run_one_reps(Runtime, Size, Reps) :-
    U = Size,                         % headline "same-N" view: universe = set size
    make_ordset(Size, U, A), make_ordset2(Size, U, B),
    make_bignum(Size, U, MA), make_bignum2(Size, U, MB),
    timed(loop_noop(Reps, 0, _), _, InfNoop),
    timed(loop_ord(Reps, A, B, 0, OA), WOrd, InfOrd),
    timed(loop_big(Reps, MA, MB, 0, OB), WBig, InfBig),
    per_op_us(WOrd, Reps, UsOrd), per_op_us(WBig, Reps, UsBig),
    InfOrdOp is (InfOrd - InfNoop) / Reps,
    InfBigOp is (InfBig - InfNoop) / Reps,
    ( UsBig > 0 -> Ratio is UsOrd / UsBig ; Ratio = inf ),
    format("KERNEL ~w size=~w reps=~w | ord_us/op=~4f big_us/op=~6f ratio=~2fx | ord_inf/op=~1f big_inf/op=~2f | ord_wall=~3fs big_wall=~3fs | chk ord=~w big=~w~n",
           [Runtime, Size, Reps, UsOrd, UsBig, Ratio, InfOrdOp, InfBigOp, WOrd, WBig, OA, OB]).

% realistic pairing: ordset of size S within a bucket of universe U;
% bignum spans the whole bucket U (its cost is U-bound, not S-bound).
run_pair(Runtime, Label, S, U, Reps) :-
    make_ordset(S, U, A), make_ordset2(S, U, B),
    make_bignum(S, U, MA), make_bignum2(S, U, MB),
    timed(loop_noop(Reps, 0, _), _, InfNoop),
    timed(loop_ord(Reps, A, B, 0, OA), WOrd, InfOrd),
    timed(loop_big(Reps, MA, MB, 0, OB), WBig, _),
    per_op_us(WOrd, Reps, UsOrd), per_op_us(WBig, Reps, UsBig),
    InfOrdOp is (InfOrd - InfNoop) / Reps,
    ( UsBig > 0 -> Ratio is UsOrd / UsBig ; Ratio = inf ),
    format("PAIR ~w ~w setS=~w bucketU=~w reps=~w | ord_us/op=~4f big_us/op=~6f ratio=~2fx | ord_inf/op=~1f | chk ~w/~w~n",
           [Runtime, Label, S, U, Reps, UsOrd, UsBig, Ratio, InfOrdOp, OA, OB]).

run_all(Runtime) :-
    forall(member(Sz, [100,1000,10000,24000]), run_one(Runtime, Sz)).
