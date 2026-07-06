#!/usr/bin/env swipl
% MEASUREMENT PROBE P1 driver (throwaway; never merged).
%
% Runs the production mrv_inc first-solution arrange search (find_crossword/6)
% from the two production corners [topleft_across, topright] over the dense
% ladder rungs + one easy control, with the core.pl P1 instrumentation enabled,
% and post-processes the event accumulators into:
%   (a) per-rung places/unplaces, retreat-distance histogram + max + %dist-1
%   (b) per-rung top-10 wipeout culprits, top-10 most-unplaced words, concentration
% Writes the full tables to a scratch file and echoes them to stdout.
%
% Usage: swipl -q benchmarks/probe_backtrack.pl -- <out-file>

:- set_prolog_flag(verbose, silent).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(pairs)).

:- prolog_load_context(directory, D),
   absolute_file_name('..', Root, [relative_to(D), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load),
   consult(Load),
   nb_setval(repo_root, Root).

:- initialization(main, main).

% The four subjects: fixture, grid size, tag.
subject('fixtures/ladder_15x15_12w.pl', 15, control_12w).
subject('fixtures/ladder_15x15_34w.pl', 15, dense_34w).
subject('fixtures/ladder_15x15_36w.pl', 15, dense_36w).
subject('fixtures/ladder_21x21_80w.pl', 21, dense_80w).

corner(topleft_across).
corner(topright).

read_clues(File, Words) :-
    setup_call_cleanup(open(File, read, S), read_loop(S, Words), close(S)).
read_loop(S, Words) :-
    read_term(S, T, []),
    ( T == end_of_file -> throw(no_clues)
    ; T = clues(Words)  -> true
    ; read_loop(S, Words) ).

% Run one corner with a fresh accumulator; return the stats dict.
run_corner(Words, Size, Loc, Placed, Stats) :-
    crosswordsmith_core:probe_reset,
    crosswordsmith_core:probe_enable,
    ( once(crosswordsmith_core:find_crossword(mrv_inc, Size, Words, Loc, _G, P))
      -> Placed = P ; Placed = none ),
    crosswordsmith_core:probe_disable,
    crosswordsmith_core:probe_stats(Stats).

main :-
    ( current_prolog_flag(argv, [Out|_]) -> true ; Out = 'scratch_probe_p1.txt' ),
    findall(Line, subject_report(_, Line), Nested),
    append(Nested, Lines),
    atomic_list_concat(Lines, '\n', Blob),
    open(Out, write, S), write(S, Blob), nl(S), close(S),
    format("~w~n", [Blob]),
    format("~n[written to ~w]~n", [Out]).

% Emit the full report block (list of atom lines) for one subject.
subject_report(Tag, Lines) :-
    subject(Fixture, Size, Tag),
    nb_getval(repo_root, Root),
    directory_file_path(Root, Fixture, File),
    read_clues(File, Words),
    length(Words, NW),
    findall(Loc-St-NP,
            ( corner(Loc),
              run_corner(Words, Size, Loc, Pl, St),
              ( Pl == none -> NP = 0 ; length(Pl, NP) ) ),
            Corners),
    % Aggregate the two corners.
    foldl(merge_corner, Corners, agg(0,0,0,0,[],[],[]), Agg),
    Agg = agg(Places, Unplaces, Wipeouts, MaxRun, Runs, WipeW, UnplW),
    % Per-corner one-liners.
    findall(CL,
            ( member(Loc-St-NP, Corners),
              St = stats{places:P,unplaces:U,wipeouts:W,max_run:MR,runs:_,wipe_words:_,unplace_words:_},
              format(atom(CL), "    ~w~t~20|placed=~d  places=~d  unplaces=~d  wipeouts=~d  max_retreat=~d",
                     [Loc, NP, P, U, W, MR]) ),
            CLines),
    % Histogram over merged run lengths.
    histo(Runs, H1, H2, H35, H610, Hgt, TotalRuns),
    ( Unplaces > 0 -> Dist1Pct is 100.0 * H1 / Unplaces ; Dist1Pct = 0.0 ),
    % Concentration of wipeout culprits.
    sort_desc_by_count(WipeW, WipeSorted),
    sort_desc_by_count(UnplW, UnplSorted),
    top_n(10, WipeSorted, WipeTop),
    top_n(10, UnplSorted, UnplTop),
    top5_share(WipeSorted, Wipeouts, Wipe5Pct),
    top5_share(UnplSorted, Unplaces, Unpl5Pct),
    format_pairs(WipeTop, WipeStr),
    format_pairs(UnplTop, UnplStr),
    format(atom(Title), "~w   (~w, grid ~w, ~d words)", [Tag, Fixture, Size, NW]),
    format(atom(A1), "(a) TREE CHURN + RETREAT DISTANCE (aggregate of both corners)", []),
    format(atom(A2), "    total places   = ~d", [Places]),
    format(atom(A3), "    total unplaces = ~d   (search-tree churn)", [Unplaces]),
    format(atom(A4), "    total wipeouts = ~d   (Placeable==[] deadends)", [Wipeouts]),
    format(atom(A5), "    retreat runs   = ~d", [TotalRuns]),
    format(atom(A6), "    retreat-distance histogram (words retracted per forward move):", []),
    format(atom(A7), "        dist 1    : ~d", [H1]),
    format(atom(A8), "        dist 2    : ~d", [H2]),
    format(atom(A9), "        dist 3-5  : ~d", [H35]),
    format(atom(A10),"        dist 6-10 : ~d", [H610]),
    format(atom(A11),"        dist >10  : ~d", [Hgt]),
    format(atom(A12),"        max retreat distance = ~d", [MaxRun]),
    format(atom(A13),"    %% of unplace events in distance-1 retreats = ~2f%", [Dist1Pct]),
    format(atom(B1), "(b) FAILURE CLUSTERING", []),
    format(atom(B2), "    top wipeout culprits (last-placed word at a deadend):", []),
    format(atom(B3), "        ~w", [WipeStr]),
    format(atom(B4), "        top-5 culprits = ~2f%% of all wipeouts", [Wipe5Pct]),
    format(atom(B5), "    top most-unplaced words (retraction churn):", []),
    format(atom(B6), "        ~w", [UnplStr]),
    format(atom(B7), "        top-5 unplaced = ~2f%% of all unplaces", [Unpl5Pct]),
    append([[ '',
              '============================================================',
              Title,
              '------------------------------------------------------------',
              '  per-corner:' ],
            CLines,
            [ '', A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13,
              '', B1, B2, B3, B4, B5, B6, B7 ]],
           Lines).

merge_corner(_Loc-St-_NP, agg(P0,U0,W0,MR0,R0,WW0,UW0), agg(P,U,W,MR,R,WW,UW)) :-
    St = stats{places:P1,unplaces:U1,wipeouts:W1,max_run:MR1,
               runs:R1,wipe_words:WW1,unplace_words:UW1},
    P is P0+P1, U is U0+U1, W is W0+W1, MR is max(MR0,MR1),
    add_count_lists(R0, R1, R),
    add_count_lists(WW0, WW1, WW),
    add_count_lists(UW0, UW1, UW).

% Merge two Key-Count lists by summing counts on shared keys.
add_count_lists(A, B, Merged) :-
    append(A, B, All),
    keysort(All, Sorted),
    group_pairs_by_key(Sorted, Groups),
    findall(K-Sum, (member(K-Vs, Groups), sum_list(Vs, Sum)), Merged).

% Histogram buckets over RunLen-Count pairs.
histo(Runs, H1, H2, H35, H610, Hgt, Total) :-
    foldl(histo_one, Runs, h(0,0,0,0,0,0), h(H1,H2,H35,H610,Hgt,Total)).
histo_one(Len-Cnt, h(A,B,C,D,E,T), h(A2,B2,C2,D2,E2,T2)) :-
    T2 is T + Cnt,
    ( Len =:= 1        -> A2 is A+Cnt, B2=B, C2=C, D2=D, E2=E
    ; Len =:= 2        -> B2 is B+Cnt, A2=A, C2=C, D2=D, E2=E
    ; Len =< 5         -> C2 is C+Cnt, A2=A, B2=B, D2=D, E2=E
    ; Len =< 10        -> D2 is D+Cnt, A2=A, B2=B, C2=C, E2=E
    ; E2 is E+Cnt, A2=A, B2=B, C2=C, D2=D ).

sort_desc_by_count(Pairs, Sorted) :-
    map_list_to_pairs(neg_count, Pairs, Keyed),
    keysort(Keyed, K1),
    pairs_values(K1, Sorted).
neg_count(_Word-Count, Neg) :- Neg is -Count.

top_n(N, List, Top) :- length(Pre, N), append(Pre, _, List), !, Top = Pre.
top_n(_, List, List).

top5_share(Sorted, Total, Pct) :-
    top_n(5, Sorted, Top5),
    pairs_values(Top5, Vs), sum_list(Vs, S),
    ( Total > 0 -> Pct is 100.0 * S / Total ; Pct = 0.0 ).

format_pairs([], '(none)').
format_pairs([P|Ps], Str) :-
    findall(A, (member(W-C, [P|Ps]), format(atom(A), "~w:~d", [W,C])), As),
    atomic_list_concat(As, '  ', Str).
