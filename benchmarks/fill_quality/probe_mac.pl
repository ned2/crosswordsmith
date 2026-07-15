% benchmarks/fill_quality/probe_mac.pl - DP-7 evidence spike: deterministic MAC fill.
%
% NOT SHIPPABLE ENGINEERING. This is a throwaway probe measuring which parts
% of ingrid_core's search stack close the DP-6 reference row (blocked_13a,
% STW, --min-score <= 30) that budget (x20) and ordering (0/8 seeds)
% measurably cannot: probe_mac/4 is plain maintained arc consistency (AC-3
% to fixpoint over the crossing graph), probe_mac_restarts/5,6 adds seeded
% restarts and an optional fillability value order, and probe_mac_dwd/6 adds
% dom/wdeg conflict-weight slot ordering (the DP-7 amendment: the ingredient
% that actually closes the row - see ./README.md). Pure bignum-mask
% bookkeeping over the F-H2 build_masks/2 bit-space (bit i = bucket index i;
% buckets are score-desc, so ascending-bit materialization preserves the
% 8.4a value order). Search does no cell-variable unification and no emit -
% completion, node count, and wall time are the outputs; probe_mac_dwd
% verifies a completed fill afterwards by binding it onto the grid's shared
% cell variables.
%
% Run from the repo root (dict = a scored wordlist, e.g. STW - not bundled):
%   swipl -q -l load.pl -l benchmarks/fill_quality/probe_mac.pl \
%     -g "probe_mac:probe_mac('grids/blocked_13a.json', 'stw.txt', 30, 600)" -g halt
% DP-7 results against this probe are recorded in ./README.md (MAC probes
% section); the SWI \ int64-boundary arithmetic quirk that bit this probe
% is documented there too.
:- module(probe_mac, [probe_mac/4, probe_mac_restarts/5, probe_mac_restarts/6,
                      probe_mac_dwd/6, probe_mac_dwd/7, replay_engine_fill/2]).

:- use_module(library(assoc)).
:- use_module(library(apply)).
:- use_module(library(lists)).
:- use_module(library(pairs)).
:- use_module(library(time)).
:- use_module(library(yall)).

probe_mac(GridFile, DictFile, MinScore, TimeoutSecs) :-
    crosswordsmith_fill:fill_grid(GridFile, Size, Slots, _CellVar),
    length(Slots, NSlots),
    format(user_error, 'mac: grid ~w, ~w slots~n', [Size, NSlots]),
    crosswordsmith_fill:load_dict(DictFile, [min_score(MinScore)], DictByLen, Index, _Scores),
    crosswordsmith_fill:build_masks(Index, MaskAssoc0),
    flat_masks(MaskAssoc0, LmA),
    format(user_error, 'mac: dict + masks built~n', []),
    % Slot ids 0..N-1; per-slot: length, bucket compound (arg/3 O(1) access).
    numlist(0, NSlots, _), NMax is NSlots - 1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, Slots),
    maplist(slot_len, Slots, Lens),
    pairs_keys_values(IdLens, Ids, Lens),
    list_to_assoc(IdLens, LenA),
    buckets_compound(DictByLen, BucketA),
    % Crossing edges: cell shared by slots A@i and B@j -> two directed edges.
    crossing_edges(IdSlots, EdgeA),
    % Initial domains: the full length bucket (no seeds in this probe).
    maplist(full_domain(DictByLen), Lens, Doms0),
    pairs_keys_values(IdDoms, Ids, Doms0),
    list_to_assoc(IdDoms, DomA0),
    nb_setval(mac_nodes, 0),
    nb_setval(mac_props, 0),
    nb_setval(mac_timed_out, false),
    get_time(T0),
    (   catch(call_with_time_limit(TimeoutSecs,
                  mac_run(Ids, DomA0, LenA, BucketA, EdgeA, LmA, Fill)),
              time_limit_exceeded,
              ( report(T0, timeout), nb_setval(mac_timed_out, true), fail ))
    ->  report(T0, completed),
        length(Fill, NF),
        format(user_error, 'mac: FILLED ~w slots~n', [NF]),
        msort(Fill, SortedFill), length(SortedFill, NF0), sort(Fill, Dedup), length(Dedup, NF1),
        ( NF0 == NF1 -> true ; format(user_error, 'mac: WARNING duplicate answers~n', []) )
    ;   % A timeout also lands here (the catch handler fails); only an
        % un-timed-out failure is a genuine exhaustion proof.
        (   nb_getval(mac_timed_out, true)
        ->  format(user_error, 'mac: TIMED OUT - completion unknown (NOT an infeasibility proof)~n', [])
        ;   report(T0, exhausted),
            format(user_error, 'mac: search space exhausted - INFEASIBLE (proven)~n', [])
        ),
        fail
    ).

slot_len(slot(_, _, _, Vars), Len) :- length(Vars, Len).

% Debug: solve GridFile with the ENGINE, then replay its solution word-by-word
% through the probe's domain machinery (place_word + full AC propagate),
% reporting the first step where the known-good fill is rejected.
replay_engine_fill(GridFile, DictFile) :-
    crosswordsmith_fill:fill_grid(GridFile, _, Slots, _),
    crosswordsmith_fill:load_dict(DictFile, DictByLen, Index),
    copy_term(Slots, FreshSlots),
    crosswordsmith_fill:fill_attempt(Slots, Slots, DictByLen, Index, filled, _, _),
    crosswordsmith_fill:build_masks(Index, MaskAssoc0),
    flat_masks(MaskAssoc0, LmA),
    length(FreshSlots, N), NMax is N - 1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, FreshSlots),
    maplist(slot_len, FreshSlots, Lens),
    pairs_keys_values(IdLens, Ids, Lens), list_to_assoc(IdLens, LenA),
    crossing_edges(IdSlots, EdgeA),
    maplist(full_domain(DictByLen), Lens, Doms0),
    pairs_keys_values(IdDoms, Ids, Doms0), list_to_assoc(IdDoms, DomA0),
    nb_setval(mac_props, 0),
    pairs_keys_values(IdBound, Ids, Slots),   % bound twins, same ids
    replay_steps(IdBound, Ids, DomA0, LenA, DictByLen, EdgeA, LmA),
    format(user_error, 'replay done~n', []).

replay_steps([], _, _, _, _, _, _).
replay_steps([Id-slot(St, Dir, _, BVars)|Steps], Unfilled, DomA0, LenA, DictByLen, EdgeA, LmA) :-
    get_assoc(Id, DomA0, Dom),
    get_assoc(Id, LenA, Len),
    get_assoc(Len, DictByLen, Bucket),
    (   nth0(BitI, Bucket, BVars)
    ->  (   Dom /\ (1 << BitI) =\= 0
        ->  true
        ;   format(user_error, 'REPLAY FAIL: slot ~w ~w: solution word ~w pruned from domain~n',
                   [St, Dir, BVars])
        )
    ;   format(user_error, 'REPLAY ODD: slot ~w ~w word ~w not in bucket~n', [St, Dir, BVars])
    ),
    selectchk(Id, Unfilled, Rest),
    edges_of(EdgeA, Id, Es),
    (   place_word(Es, BVars, Rest, LenA, LmA, DomA0, DomA1, [], Dirty)
    ->  (   propagate(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2)
        ->  true
        ;   format(user_error, 'REPLAY FAIL: propagate wiped a domain after ~w ~w (~w)~n',
                   [St, Dir, BVars]),
            DomA2 = DomA1
        )
    ;   format(user_error, 'REPLAY FAIL: place_word zeroed a crosser at ~w ~w (~w)~n',
               [St, Dir, BVars]),
        DomA2 = DomA0
    ),
    replay_steps(Steps, Rest, DomA2, LenA, DictByLen, EdgeA, LmA).

report(T0, What) :-
    get_time(T1), W is T1 - T0,
    nb_getval(mac_nodes, Nodes),
    nb_getval(mac_props, Props),
    format(user_error, 'mac: ~w after ~2f s, ~w nodes, ~w propagation domain-updates~n',
           [What, W, Nodes, Props]).

% root propagation, then search
mac_run(Ids, DomA0, LenA, BucketA, EdgeA, LmA, Fill) :-
    propagate(Ids, DomA0, LenA, EdgeA, LmA, [], DomA1),
    format(user_error, 'mac: root AC fixpoint reached~n', []),
    mac_search(Ids, DomA1, LenA, BucketA, EdgeA, LmA, [], Fill).

% --- MAC + restarts (the ingrid recipe minus dom/wdeg) -------------------------
% Weighted-random pick among the top 3 candidates (ingrid's
% RANDOM_WORD_WEIGHTS = [4,2,1]) via core's splitmix64, restart on a growing
% node cap (500 x1.5 per attempt, ingrid's 500 x1.1 tightened for spike wall
% economics), fresh derived seed per attempt. Slot order stays deterministic
% MRV; no conflict-weight learning. probe_mac_restarts/5 shares all setup.
probe_mac_restarts(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed) :-
    probe_mac_restarts(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed, score).

% ValueOrder: `score` (the §8.4a bucket order - quality-first) or `fillability`
% (ingrid's dominant value signal, grid_config.rs:209-231: order candidates by
% how much crossing support their letters keep, x900 vs score x5 - i.e. the
% OPPOSITE priority to the §8.4a quality contract; here: sum over positions of
% log10(letter-mask popcount + 1), desc, precomputed per bucket).
probe_mac_restarts(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed, ValueOrder) :-
    crosswordsmith_fill:fill_grid(GridFile, Size, Slots, _CellVar),
    length(Slots, NSlots),
    format(user_error, 'macr: grid ~w, ~w slots~n', [Size, NSlots]),
    crosswordsmith_fill:load_dict(DictFile, [min_score(MinScore)], DictByLen, Index, _Scores),
    crosswordsmith_fill:build_masks(Index, MaskAssoc0),
    flat_masks(MaskAssoc0, LmA),
    format(user_error, 'macr: dict + masks built~n', []),
    NMax is NSlots - 1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, Slots),
    maplist(slot_len, Slots, Lens),
    pairs_keys_values(IdLens, Ids, Lens), list_to_assoc(IdLens, LenA),
    buckets_compound(DictByLen, BucketA),
    crossing_edges(IdSlots, EdgeA),
    maplist(full_domain(DictByLen), Lens, Doms0),
    pairs_keys_values(IdDoms, Ids, Doms0), list_to_assoc(IdDoms, DomA0),
    nb_setval(mac_nodes, 0),
    nb_setval(mac_props, 0),
    nb_setval(mac_timed_out, false),
    get_time(T0),
    (   ValueOrder == fillability
    ->  fill_orders(BucketA, LmA, LenA, FOrdA), VO = ford(FOrdA),
        format(user_error, 'macr: fillability value order precomputed~n', [])
    ;   VO = score
    ),
    (   catch(call_with_time_limit(TimeoutSecs,
                  ( propagate(Ids, DomA0, LenA, EdgeA, LmA, [], DomA1),
                    format(user_error, 'macr: root AC fixpoint reached~n', []),
                    restart_loop(1, 500, BaseSeed, T0,
                                 Ids, DomA1, LenA, BucketA, EdgeA, LmA, VO, Fill) )),
              time_limit_exceeded,
              ( report(T0, 'timeout (restarts)'), nb_setval(mac_timed_out, true), fail ))
    ->  report(T0, 'completed (restarts)'),
        length(Fill, NF),
        format(user_error, 'macr: FILLED ~w slots~n', [NF])
    ;   (   nb_getval(mac_timed_out, true)
        ->  format(user_error, 'macr: TIMED OUT - completion unknown~n', [])
        ;   report(T0, 'attempts exhausted (restarts)')
        ),
        fail
    ).

restart_loop(Attempt, Cap, BaseSeed, T0, Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Fill) :-
    Seed is BaseSeed + Attempt,
    crosswordsmith_core:set_search_seed(Seed),
    nb_setval(mac_attempt_nodes, 0),
    get_time(T), El is T - T0,
    format(user_error, 'macr: attempt ~w cap ~w seed ~w (t=~1f s)~n', [Attempt, Cap, Seed, El]),
    (   catch(mac_search_r(Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Cap, [], Fill),
              mac_cap, fail)
    ->  true
    ;   Attempt < 10000,
        Cap1 is ceiling(Cap * 1.5),
        A1 is Attempt + 1,
        restart_loop(A1, Cap1, BaseSeed, T0, Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Fill)
    ).

mac_search_r([], _DomA, _LenA, _BucketA, _EdgeA, _LmA, _VO, _Cap, _Used, []).
mac_search_r(Unfilled, DomA, LenA, BucketA, EdgeA, LmA, VO, Cap, Used, [Word|Fill]) :-
    Unfilled = [_|_],
    select_mrv_id(Unfilled, DomA, Best),
    selectchk(Best, Unfilled, Rest),
    get_assoc(Best, DomA, Dom),
    get_assoc(Best, LenA, Len),
    get_assoc(Len, BucketA, Bucket),
    materialize(VO, Len, Dom, Bucket, Cands0),
    pick_top3(Cands0, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    nb_getval(mac_nodes, N0), N1 is N0 + 1, nb_setval(mac_nodes, N1),
    nb_getval(mac_attempt_nodes, A0), A1 is A0 + 1, nb_setval(mac_attempt_nodes, A1),
    ( A1 > Cap -> throw(mac_cap) ; true ),
    edges_of(EdgeA, Best, Es),
    place_word(Es, Word, Rest, LenA, LmA, DomA, DomA1, [], Dirty),
    propagate(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2),
    mac_search_r(Rest, DomA2, LenA, BucketA, EdgeA, LmA, VO, Cap, [Word|Used], Fill).

materialize(score, _Len, Dom, Bucket, Cands) :-
    mask_words(Dom, Bucket, Cands).
materialize(ford(FOrdA), Len, Dom, Bucket, Cands) :-
    get_assoc(Len, FOrdA, Order),
    ford_words(Order, Dom, Bucket, Cands).

% walk the precomputed fillability-desc index order, keeping live bits
ford_words([], _, _, []).
ford_words([I|Is], Dom, Bucket, Cands) :-
    (   Dom /\ (1 << I) =\= 0
    ->  J is I + 1, arg(J, Bucket, W), Cands = [W|Ws]
    ;   Cands = Ws
    ),
    ford_words(Is, Dom, Bucket, Ws).

% Per-length fillability order: word score = sum over positions of
% log10(popcount(letter mask at that position) + 1), descending. Popcounts
% come straight off the flat letter-mask table.
fill_orders(BucketA, LmA, _LenA, FOrdA) :-
    assoc_to_list(BucketA, Pairs),
    maplist(fill_order_len(LmA), Pairs, FPairs),
    list_to_assoc(FPairs, FOrdA).
fill_order_len(LmA, Len-Bucket, Len-Order) :-
    functor(Bucket, b, N),
    numlist(1, N, Js),
    maplist(word_fillability(LmA, Len, Bucket), Js, Keyed),
    msort(Keyed, Sorted),
    findall(I, member(_-I, Sorted), Order).
word_fillability(LmA, Len, Bucket, J, NegF-I) :-
    I is J - 1,
    arg(J, Bucket, W),
    word_fill_score(W, 0, Len, LmA, 0.0, F),
    NegF is -F.
word_fill_score([], _, _, _, F, F).
word_fill_score([Ch|Cs], Pos, Len, LmA, F0, F) :-
    lm_of(LmA, Len, Pos, Lm),
    letter_arg(Ch, ChI),
    arg(ChI, Lm, M),
    P is popcount(M),
    F1 is F0 + log10(P + 1),
    Pos1 is Pos + 1,
    word_fill_score(Cs, Pos1, Len, LmA, F1, F).

% weighted-random promotion of one of the top 3 candidates ([4,2,1]); the
% rest keep their score-desc order, so the attempt remains complete.
pick_top3([], []).
pick_top3([A], [A]).
pick_top3([A, B], Out) :-
    crosswordsmith_core:prng_draw(V), W is V mod 6,
    ( W < 4 -> Out = [A, B] ; Out = [B, A] ).
pick_top3([A, B, C|Rest], [P|Ps]) :-
    crosswordsmith_core:prng_draw(V), W is V mod 7,
    (   W < 4 -> P = A, Ps = [B, C|Rest]
    ;   W < 6 -> P = B, Ps = [A, C|Rest]
    ;   P = C, Ps = [A, B|Rest]
    ).

% --- MAC + restarts + dom/wdeg (the full ingrid recipe at spike speeds) --------
% The DP-7 residual: conflict-weight VARIABLE ordering (Boussemart dom/wdeg)
% was never probed. Classic formulation: every crossing (binary constraint)
% carries a weight, init 1.0, +1 whenever ITS revision wipes a domain (in
% place_word or the AC worklist); variable order picks the unfilled slot
% minimizing popcount(dom) / sum(weights of crossings to unfilled neighbors).
% Weights are GLOBAL LEARNING STATE: they survive backtracking (flat compound
% + nb_setarg) and persist across restarts with x0.99 aging per attempt
% (ingrid's recipe per the 1.3.1 source read). Edges gain an id
% (e(Pos,T,TPos,EId) shared by both directions), so this variant has its own
% place/revise/select twins - the three DP-7 configurations are untouched.
probe_mac_dwd(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed, ValueOrder) :-
    probe_mac_dwd(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed, ValueOrder, top3).

% Pick: `top3` (ingrid's weighted-random [4,2,1] promotion - needs the PRNG)
% or `greedy` (strict best-first; the PRNG is NEVER drawn, so the whole run
% is RNG-free - restart diversity comes only from the learned weights
% reordering slot selection between attempts).
probe_mac_dwd(GridFile, DictFile, MinScore, TimeoutSecs, BaseSeed, ValueOrder, Pick) :-
    nb_setval(mac_pick, Pick),
    crosswordsmith_fill:fill_grid(GridFile, Size, Slots, _CellVar),
    length(Slots, NSlots),
    format(user_error, 'macw: grid ~w, ~w slots~n', [Size, NSlots]),
    crosswordsmith_fill:load_dict(DictFile, [min_score(MinScore)], DictByLen, Index, Scores),
    crosswordsmith_fill:build_masks(Index, MaskAssoc0),
    flat_masks(MaskAssoc0, LmA),
    format(user_error, 'macw: dict + masks built~n', []),
    NMax is NSlots - 1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, Slots),
    maplist(slot_len, Slots, Lens),
    pairs_keys_values(IdLens, Ids, Lens), list_to_assoc(IdLens, LenA),
    buckets_compound(DictByLen, BucketA),
    crossing_edges_w(IdSlots, EdgeA, NEdges),
    format(user_error, 'macw: ~w crossings~n', [NEdges]),
    init_weights(NEdges),
    maplist(full_domain(DictByLen), Lens, Doms0),
    pairs_keys_values(IdDoms, Ids, Doms0), list_to_assoc(IdDoms, DomA0),
    nb_setval(mac_nodes, 0),
    nb_setval(mac_props, 0),
    nb_setval(mac_timed_out, false),
    get_time(T0),
    (   ValueOrder == fillability
    ->  fill_orders(BucketA, LmA, LenA, FOrdA), VO = ford(FOrdA),
        format(user_error, 'macw: fillability value order precomputed~n', [])
    ;   VO = score
    ),
    (   catch(call_with_time_limit(TimeoutSecs,
                  ( propagate_w(Ids, DomA0, LenA, EdgeA, LmA, [], DomA1),
                    format(user_error, 'macw: root AC fixpoint reached~n', []),
                    restart_loop_w(1, 500, BaseSeed, T0,
                                   Ids, DomA1, LenA, BucketA, EdgeA, LmA, VO, Fill) )),
              time_limit_exceeded,
              ( report(T0, 'timeout (dom/wdeg)'), nb_setval(mac_timed_out, true), fail ))
    ->  report(T0, 'completed (dom/wdeg)'),
        length(Fill, NF),
        format(user_error, 'macw: FILLED ~w slots~n', [NF]),
        maplist([_-W]>>(atom_chars(A, W), format(user_error, '  ~w~n', [A])), Fill),
        % hard verification: bind every placed word onto the grid's SHARED
        % cell variables - crossing inconsistency makes unification fail.
        (   maplist(bind_fill(IdSlots), Fill)
        ->  format(user_error, 'macw: crossing-consistency VERIFIED (all shared cells unify)~n', [])
        ;   format(user_error, 'macw: CONSISTENCY FAILURE - fill is bogus~n', [])
        ),
        fill_score_stats(Fill, Scores)
    ;   (   nb_getval(mac_timed_out, true)
        ->  format(user_error, 'macw: TIMED OUT - completion unknown~n', [])
        ;   report(T0, 'attempts exhausted (dom/wdeg)')
        ),
        fail
    ).

restart_loop_w(Attempt, Cap, BaseSeed, T0, Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Fill) :-
    Seed is BaseSeed + Attempt,
    crosswordsmith_core:set_search_seed(Seed),
    nb_setval(mac_attempt_nodes, 0),
    get_time(T), El is T - T0,
    format(user_error, 'macw: attempt ~w cap ~w seed ~w (t=~1f s)~n', [Attempt, Cap, Seed, El]),
    (   catch(mac_search_w(Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Cap, [], Fill),
              mac_cap, fail)
    ->  true
    ;   Attempt < 10000,
        age_weights(0.99),                       % weights persist, aged
        Cap1 is ceiling(Cap * 1.5),
        A1 is Attempt + 1,
        restart_loop_w(A1, Cap1, BaseSeed, T0, Ids, DomA, LenA, BucketA, EdgeA, LmA, VO, Fill)
    ).

bind_fill(IdSlots, Id-Word) :-
    memberchk(Id-slot(_, _, _, Vars), IdSlots),
    Vars = Word.

fill_score_stats(Fill, Scores) :-
    (   Scores = scores(SA),
        findall(S, ( member(_-W, Fill), get_assoc(W, SA, S) ), Ss),
        Ss \== []
    ->  sum_list(Ss, Sum), length(Ss, N), min_list(Ss, Min), length(Fill, NF),
        Mean is Sum / N,
        format(user_error, 'macw: fill quality mean ~2f min ~w (~w of ~w scored)~n',
               [Mean, Min, N, NF])
    ;   true
    ).

mac_search_w([], _DomA, _LenA, _BucketA, _EdgeA, _LmA, _VO, _Cap, _Used, []).
mac_search_w(Unfilled, DomA, LenA, BucketA, EdgeA, LmA, VO, Cap, Used, [Best-Word|Fill]) :-
    Unfilled = [_|_],
    select_dwd_id(Unfilled, DomA, EdgeA, Best),
    selectchk(Best, Unfilled, Rest),
    get_assoc(Best, DomA, Dom),
    get_assoc(Best, LenA, Len),
    get_assoc(Len, BucketA, Bucket),
    materialize(VO, Len, Dom, Bucket, Cands0),
    (   nb_current(mac_pick, greedy)
    ->  Cands = Cands0
    ;   pick_top3(Cands0, Cands)
    ),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    nb_getval(mac_nodes, N0), N1 is N0 + 1, nb_setval(mac_nodes, N1),
    ( N1 mod 10000 =:= 0 -> report(0.0, progress) ; true ),
    nb_getval(mac_attempt_nodes, A0), A1 is A0 + 1, nb_setval(mac_attempt_nodes, A1),
    ( A1 > Cap -> throw(mac_cap) ; true ),
    edges_of(EdgeA, Best, Es),
    place_word_w(Es, Word, Rest, LenA, LmA, DomA, DomA1, [], Dirty),
    propagate_w(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2),
    mac_search_w(Rest, DomA2, LenA, BucketA, EdgeA, LmA, VO, Cap, [Word|Used], Fill).

% dom/wdeg selection: minimize C/W. Compared cross-multiplied to stay in
% arithmetic (C1*W2 < C2*W1); W >= its crossing count >= 0, +0.001 floor
% guards fully-isolated slots.
select_dwd_id([Id|Ids], DomA, EdgeA, Best) :-
    dwd_key(Id, DomA, EdgeA, [Id|Ids], C, W),
    select_dwd_walk(Ids, DomA, EdgeA, [Id|Ids], Id, C, W, Best).
select_dwd_walk([], _, _, _, Best, _, _, Best).
select_dwd_walk([Id|Ids], DomA, EdgeA, Unfilled, B0, C0, W0, Best) :-
    dwd_key(Id, DomA, EdgeA, Unfilled, C, W),
    (   C * W0 < C0 * W
    ->  select_dwd_walk(Ids, DomA, EdgeA, Unfilled, Id, C, W, Best)
    ;   select_dwd_walk(Ids, DomA, EdgeA, Unfilled, B0, C0, W0, Best)
    ).
dwd_key(Id, DomA, EdgeA, Unfilled, C, W) :-
    get_assoc(Id, DomA, D), C is popcount(D),
    edges_of(EdgeA, Id, Es),
    nb_getval(mac_wts, WT),
    dwd_wdeg(Es, Unfilled, WT, 0.001, W).
dwd_wdeg([], _, _, W, W).
dwd_wdeg([e(_, T, _, EId)|Es], Unfilled, WT, W0, W) :-
    (   memberchk(T, Unfilled)
    ->  arg(EId, WT, V), W1 is W0 + V
    ;   W1 = W0
    ),
    dwd_wdeg(Es, Unfilled, WT, W1, W).

% global constraint weights: flat w/NEdges compound in a global variable,
% destructively bumped (nb_setarg) so learning survives backtracking.
init_weights(N) :-
    length(Ones, N), maplist(=(1.0), Ones),
    WT =.. [w|Ones],
    nb_setval(mac_wts, WT).
bump_edge(EId) :-
    nb_getval(mac_wts, WT),
    arg(EId, WT, V), V1 is V + 1.0,
    nb_setarg(EId, WT, V1).
age_weights(Factor) :-
    nb_getval(mac_wts, WT),
    functor(WT, w, N),
    forall(between(1, N, I),
           ( arg(I, WT, V), V1 is V * Factor, nb_setarg(I, WT, V1) )).

% place/revise twins over e/4 edges: identical pruning, but a wipeout bumps
% the responsible crossing's weight BEFORE failing.
place_word_w([], _Word, _Unfilled, _LenA, _LmA, DomA, DomA, Dirty, Dirty).
place_word_w([e(Pos, T, TPos, EId)|Es], Word, Unfilled, LenA, LmA, DomA0, DomA, Dirty0, Dirty) :-
    (   memberchk(T, Unfilled)
    ->  nth0(Pos, Word, Ch),
        get_assoc(T, LenA, TLen),
        lm_of(LmA, TLen, TPos, LmT),
        letter_arg(Ch, ChI),
        arg(ChI, LmT, ChMask),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ ChMask,
        (   TDom1 =:= 0
        ->  bump_edge(EId), fail
        ;   true
        ),
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Dirty1 = Dirty0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            Dirty1 = [T|Dirty0],
            nb_getval(mac_props, P0), P1 is P0 + 1, nb_setval(mac_props, P1)
        )
    ;   DomA1 = DomA0, Dirty1 = Dirty0
    ),
    place_word_w(Es, Word, Unfilled, LenA, LmA, DomA1, DomA, Dirty1, Dirty).

propagate_w([], DomA, _LenA, _EdgeA, _LmA, _Unfilled, DomA).
propagate_w([S|Queue], DomA0, LenA, EdgeA, LmA, Unfilled, DomA) :-
    edges_of(EdgeA, S, Es),
    get_assoc(S, DomA0, SDom),
    get_assoc(S, LenA, SLen),
    revise_edges_w(Es, S, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA1, Queue, Queue1),
    propagate_w(Queue1, DomA1, LenA, EdgeA, LmA, Unfilled, DomA).

revise_edges_w([], _S, _SDom, _SLen, _LenA, _LmA, _Unfilled, DomA, DomA, Q, Q).
revise_edges_w([e(Pos, T, TPos, EId)|Es], S, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA, Q0, Q) :-
    (   memberchk(T, Unfilled)
    ->  get_assoc(T, LenA, TLen),
        lm_of(LmA, SLen, Pos, LmS),
        lm_of(LmA, TLen, TPos, LmT),
        support_mask(LmS, SDom, LmT, Supp),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ Supp,
        (   TDom1 =:= 0
        ->  bump_edge(EId), fail
        ;   true
        ),
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Q1 = Q0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            ( memberchk(T, Q0) -> Q1 = Q0 ; Q1 = [T|Q0] ),
            nb_getval(mac_props, P0), P1 is P0 + 1, nb_setval(mac_props, P1)
        )
    ;   DomA1 = DomA0, Q1 = Q0
    ),
    revise_edges_w(Es, S, SDom, SLen, LenA, LmA, Unfilled, DomA1, DomA, Q1, Q).

% crossing edges with a shared undirected edge id (1-based, for arg/nb_setarg)
crossing_edges_w(IdSlots, EdgeA, NEdges) :-
    findall(Cell-(Id-Pos),
            ( member(Id-slot(_, _, Cells, _), IdSlots),
              nth0(Pos, Cells, Cell) ),
            CellOccs),
    msort(CellOccs, Sorted),
    group_pairs_by_key(Sorted, Groups),
    findall(x(A, PA, B, PB), member(_-[A-PA, B-PB], Groups), Crossings),
    length(Crossings, NEdges),
    findall(Edge,
            ( nth1(EId, Crossings, x(A, PA, B, PB)),
              ( Edge = A-e(PA, B, PB, EId) ; Edge = B-e(PB, A, PA, EId) ) ),
            AllEdges),
    msort(AllEdges, SortedEdges),
    group_pairs_by_key(SortedEdges, EdgeGroups),
    list_to_assoc(EdgeGroups, EdgeA).

% buckets: assoc Len -> b(W1, ..., Wn) compound (or absent -> no words)
buckets_compound(DictByLen, BucketA) :-
    assoc_to_list(DictByLen, Pairs),
    maplist([L-Ws, L-T]>>(T =.. [b|Ws]), Pairs, BPairs),
    list_to_assoc(BPairs, BucketA).

full_domain(DictByLen, Len, Dom) :-
    ( get_assoc(Len, DictByLen, Ws) -> length(Ws, N), Dom is (1 << N) - 1 ; Dom = 0 ).

% Crossing edges from shared cell ids: EdgeA maps Id -> list of e(Pos, OtherId, OtherPos).
crossing_edges(IdSlots, EdgeA) :-
    findall(Cell-(Id-Pos),
            ( member(Id-slot(_, _, Cells, _), IdSlots),
              nth0(Pos, Cells, Cell) ),
            CellOccs),
    msort(CellOccs, Sorted),
    group_pairs_by_key(Sorted, Groups),
    findall(A-e(PA, B, PB),
            ( member(_-[A-PA, B-PB], Groups) ),
            Fwd),
    findall(B-e(PB, A, PA),
            ( member(_-[A-PA, B-PB], Groups) ),
            Bwd),
    append(Fwd, Bwd, AllEdges),
    msort(AllEdges, SortedEdges),
    group_pairs_by_key(SortedEdges, EdgeGroups),
    list_to_assoc(EdgeGroups, EdgeA).

edges_of(EdgeA, Id, Es) :- ( get_assoc(Id, EdgeA, Es) -> true ; Es = [] ).

% --- the MAC search -----------------------------------------------------------
mac_search([], _DomA, _LenA, _BucketA, _EdgeA, _LmA, _Used, []).
mac_search(Unfilled, DomA, LenA, BucketA, EdgeA, LmA, Used, [Word|Fill]) :-
    Unfilled = [_|_],
    select_mrv_id(Unfilled, DomA, Best),
    selectchk(Best, Unfilled, Rest),
    get_assoc(Best, DomA, Dom),
    get_assoc(Best, LenA, Len),
    get_assoc(Len, BucketA, Bucket),
    mask_words(Dom, Bucket, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    nb_getval(mac_nodes, N0), N1 is N0 + 1, nb_setval(mac_nodes, N1),
    ( N1 mod 10000 =:= 0 -> report(0.0, progress) ; true ),
    (   N1 < 60
    ->  length(Rest, Depth), atom_chars(WA, Word), length(Cands, NC),
        format(user_error, '  try n~w rest~w slot~w ~w (of ~w cands)~n',
               [N1, Depth, Best, WA, NC])
    ;   true
    ),
    % place: AND each unfilled crosser's domain with the placed letter's mask
    edges_of(EdgeA, Best, Es),
    place_word(Es, Word, Rest, LenA, LmA, DomA, DomA1, [], Dirty),
    propagate(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2),
    mac_search(Rest, DomA2, LenA, BucketA, EdgeA, LmA, [Word|Used], Fill).

select_mrv_id([Id|Ids], DomA, Best) :-
    get_assoc(Id, DomA, D), C is popcount(D),
    select_mrv_walk(Ids, DomA, Id-C, Best).
select_mrv_walk([], _, Best-_, Best).
select_mrv_walk([Id|Ids], DomA, B0-C0, Best) :-
    get_assoc(Id, DomA, D), C is popcount(D),
    ( C < C0 -> select_mrv_walk(Ids, DomA, Id-C, Best)
    ; select_mrv_walk(Ids, DomA, B0-C0, Best)
    ).

% materialize words from mask bits, ascending bit index (= bucket order =
% score-desc-then-lex: the §8.4a value order, preserved).
mask_words(0, _, []) :- !.
mask_words(M, Bucket, [W|Ws]) :-
    B is lsb(M),
    I is B + 1, arg(I, Bucket, W),
    % xor, NOT `M /\ \(1<<B)`: SWI's `\` misevaluates at the int64 boundary
    % (\(1<<63) yields 2^63-1, wiping bits >= 63 - found the hard way)
    M1 is M xor (1 << B),
    mask_words(M1, Bucket, Ws).

% placement: direct letter AND on unfilled crossers; collect dirtied ids
place_word([], _Word, _Unfilled, _LenA, _LmA, DomA, DomA, Dirty, Dirty).
place_word([e(Pos, T, TPos)|Es], Word, Unfilled, LenA, LmA, DomA0, DomA, Dirty0, Dirty) :-
    (   memberchk(T, Unfilled)
    ->  nth0(Pos, Word, Ch),
        get_assoc(T, LenA, TLen),
        lm_of(LmA, TLen, TPos, LmT),
        letter_arg(Ch, ChI),
        arg(ChI, LmT, ChMask),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ ChMask,
        TDom1 =\= 0,
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Dirty1 = Dirty0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            Dirty1 = [T|Dirty0],
            nb_getval(mac_props, P0), P1 is P0 + 1, nb_setval(mac_props, P1)
        )
    ;   DomA1 = DomA0, Dirty1 = Dirty0
    ),
    place_word(Es, Word, Unfilled, LenA, LmA, DomA1, DomA, Dirty1, Dirty).

% --- AC-3 fixpoint over the crossing graph -------------------------------------
% worklist of slot ids whose domain changed; re-filter each unfilled crosser's
% domain to the letters the changed slot still supports at the shared cell.
propagate(_, DomA, _LenA, _EdgeA, _LmA, _Unfilled, DomA) :-
    nb_current(mac_no_ac, true), !.   % FC-1-only debug mode
propagate([], DomA, _LenA, _EdgeA, _LmA, _Unfilled, DomA).
propagate([S|Queue], DomA0, LenA, EdgeA, LmA, Unfilled, DomA) :-
    edges_of(EdgeA, S, Es),
    get_assoc(S, DomA0, SDom),
    get_assoc(S, LenA, SLen),
    revise_edges(Es, S, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA1, Queue, Queue1),
    propagate(Queue1, DomA1, LenA, EdgeA, LmA, Unfilled, DomA).

revise_edges([], _S, _SDom, _SLen, _LenA, _LmA, _Unfilled, DomA, DomA, Q, Q).
revise_edges([e(Pos, T, TPos)|Es], S, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA, Q0, Q) :-
    (   memberchk(T, Unfilled)
    ->  get_assoc(T, LenA, TLen),
        lm_of(LmA, SLen, Pos, LmS),
        lm_of(LmA, TLen, TPos, LmT),
        support_mask(LmS, SDom, LmT, Supp),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ Supp,
        TDom1 =\= 0,
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Q1 = Q0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            ( memberchk(T, Q0) -> Q1 = Q0 ; Q1 = [T|Q0] ),
            nb_getval(mac_props, P0), P1 is P0 + 1, nb_setval(mac_props, P1)
        )
    ;   DomA1 = DomA0, Q1 = Q0
    ),
    revise_edges(Es, S, SDom, SLen, LenA, LmA, Unfilled, DomA1, DomA, Q1, Q).

% letters S still supports at Pos, projected into T's (TPos) mask space.
% Flat form: LmS/LmT are lm/26 compounds (arg I = letter I's mask, 0 absent),
% so the 26-letter sweep costs 2 arg/3 calls + one AND test per letter
% instead of two logarithmic assoc lookups (the naive form measured ~14ms/node
% on the 13x13/STW row - assoc + bignum overhead dominated).
support_mask(LmS, SDom, LmT, Supp) :-
    support_fold(1, LmS, SDom, LmT, 0, Supp).
support_fold(27, _, _, _, Supp, Supp) :- !.
support_fold(I, LmS, SDom, LmT, Acc, Supp) :-
    arg(I, LmS, SM),
    (   SDom /\ SM =\= 0
    ->  arg(I, LmT, TM),
        Acc1 is Acc \/ TM
    ;   Acc1 = Acc
    ),
    I1 is I + 1,
    support_fold(I1, LmS, SDom, LmT, Acc1, Supp).

letters(['A','B','C','D','E','F','G','H','I','J','K','L','M',
         'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']).

% Flat letter-mask table: assoc (Len-Pos) -> lm(MA, ..., MZ). Built once from
% build_masks/2's k(Len,Pos,Ch) assoc.
flat_masks(MaskAssoc, LmA) :-
    assoc_to_list(MaskAssoc, Pairs),
    maplist([k(L, P, Ch)-M, (L-P)-(Ch-M)]>>true, Pairs, LPPairs),
    msort(LPPairs, Sorted),
    group_pairs_by_key(Sorted, Groups),
    maplist(lm_compound, Groups, LmPairs),
    list_to_assoc(LmPairs, LmA).
lm_compound(LP-ChMs, LP-Lm) :-
    letters(Chs),
    maplist(lm_arg(ChMs), Chs, Ms),
    Lm =.. [lm|Ms].
lm_arg(ChMs, Ch, M) :- ( memberchk(Ch-M0, ChMs) -> M = M0 ; M = 0 ).

lm_of(LmA, Len, Pos, Lm) :-
    ( get_assoc(Len-Pos, LmA, Lm0) -> Lm = Lm0 ; Lm = lm(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) ).

letter_arg(Ch, I) :- char_code(Ch, C), I is C - 64.   % 'A' -> 1
