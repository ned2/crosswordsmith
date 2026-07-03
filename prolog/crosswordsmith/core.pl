% core.pl - the shared substrate: grid model, free-canvas legality core, clue
% numbering, JSON emit, input loading (design-spec §4/§6). Formerly the root
% crossword.pl; the CLI lives in the `crosswordsmith` script. This file carries
% NO initialization directive and no shebang, so loading it (via load.pl or
% directly) is side-effect free and its predicates are unit-testable in
% isolation.
%
% The export list is the verified union of every consumer's needs (metrics,
% arrange, stockgrid, fill, the CLI driver, the benchmarks — see the
% migration plan §4.7). Internals (incl. the legacy crossword/3,4 top-level,
% a benchmark-only research surface) are reached by tests as
% crosswordsmith_core:Pred(...).
%
% Copyright (C) 2011  Ned Letcher - nedned.net

:- module(crosswordsmith_core,
          [ % I/O + emit
            with_output/2,
            load_clues/2,
            emit_json/3,
            set_verbose/1,
            verbose_report/2,
            % clue numbering + layout build
            assign_clue_numbers/2,
            build_grid_rows/3,
            build_words/4,
            answer_meta_assoc/2,
            add_word_cells/3,
            % grid geometry
            cell_coord/3,
            init_grid/2,
            start_loc/4,
            start_locs/1,
            valid_loc/1,
            next_cell/4,
            fits_on_grid/4,
            % search primitives
            assign_word/10,
            find_intersecting_word/6,
            assign_words_inc/9,
            find_crossword/6,
            all_crossword/5,
            % strategy registry
            strategies/1,
            default_strategy/1,
            valid_strategy/1,
            require_strategy/1,
            % utilities
            remove_x/3,
            shares_letter/2,
            check_unique_answers/1,
            % search perturbation (the arrange --seed / --shuffle knobs)
            set_search_seed/1,
            set_shuffle_seed/1,
            current_search_seed/1
          ]).


% The program uses two simple data structures. The first is a list of
% words that have been placed on the crossword grid with each element
% being a dict of attributes of the word. See the placed words utility
% predicate section below for the structure of this dict.
%
% The second data structure holds the contents of the grid. It is a single
% compound term grid(C1, ..., C(GridLen*GridLen)) whose Num-th argument is cell
% Num (1 through GridLen*GridLen, row-major): an UNBOUND variable when the cell
% is empty, bound to a letter atom when filled. Reads are arg/3 + var/nonvar
% (O(1), allocation-free); a placement unifies the cell variable with its
% letter and is undone automatically by the trail on backtracking, so a grid's
% GIn and GOut are one and the same term (the old immutable-version threading
% collapses onto the trail). See docs/experiments.md (E-H2).

% Used by the JSON output emitter (canonically library(json) on SWI 10+).
:- use_module(library(http/json)).

% random_permutation/2 for the opt-in search perturbation (set_search_seed/1).
% NOT reached on the deterministic path (search_seed/1 unset -> identity), so a
% default run never draws from — nor even seeds — the RNG.
:- use_module(library(random)).

% aggregate_all/3, used to count solutions in all_crossword/5.
:- use_module(library(aggregate)).

% program predicates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The variable-ordering strategies the solver can run. `baseline` is the
% original input-order search; `mrv`, `mrv_capped` and `mrv_inc` are the
% fail-first variants (see select_word/9 and assign_words_inc/9). Adding a
% strategy is a one-line entry here plus a select_word/9 clause.
strategies([baseline, mrv, mrv_capped, mrv_inc]).

% The production default when no --strategy is given. mrv_inc (incremental
% capped MRV) is the best general strategy: it tames the pathological dense
% search (~1300x fewer inferences than baseline) while keeping the large-grid
% per-node cost low. See docs/experiments.md (entry E5) for the evidence.
default_strategy(mrv_inc).

valid_strategy(S) :-
    strategies(Ss),
    memberchk(S, Ss).

require_strategy(S) :-
    valid_strategy(S),
    !.
require_strategy(S) :-
    throw(error(unknown_strategy(S), _)).

valid_loc(Loc) :-
    start_locs(Locs),
    memberchk(Loc, Locs).

% Run Goal with its output sent to a file, or straight to stdout when '' (no
% --out). For a file the output is captured first and written only if Goal
% succeeds, so a no-solution run leaves no empty file behind; the stdout path
% is the unchanged direct write.
:- meta_predicate with_output(+, 0).   % Goal is called (P12)
with_output('', Goal) :-
    !,
    call(Goal).
with_output(File, Goal) :-
    with_output_to(string(Text), Goal),
    setup_call_cleanup(open(File, write, S), write(S, Text), close(S)).

% Success-summary verbosity (design-spec §5.1). Routine success summaries
% ("arrange: grid ..., reward ...", "fill: ... filled N slots") print on
% stderr only under --verbose; verbose_report/2 is their one gate. Warnings
% and compromises (dropped words, fewer-than-K candidates, cap-inert
% degeneration) and all failure reports stay unconditional per INV-3 - they
% must NOT go through verbose_report.
:- dynamic verbose_mode/0.

set_verbose(true)  :- ( verbose_mode -> true ; assertz(verbose_mode) ).
set_verbose(false) :- retractall(verbose_mode).

verbose_report(Fmt, Args) :-
    ( verbose_mode -> format(user_error, Fmt, Args) ; true ).


% Clue input loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input files converge on the same internal Words list
% ([[Answer, MetaDict], ...]) before the pipeline runs. JSON is the
% interchange format; Prolog fixture files define a clues/1 term in the same
% shape. Everything after this loader is unchanged.

load_clues(File, Words) :-
    file_name_extension(_, Ext0, File),
    downcase_atom(Ext0, Ext),
    load_clues_by_extension(Ext, File, Words).

load_clues_by_extension(json, File, Words) :-
    !,
    read_clues_json(File, Words).
load_clues_by_extension(pl, File, Words) :-
    !,
    read_clues_prolog(File, Words).
load_clues_by_extension(Ext, File, _Words) :-
    throw(error(unsupported_clue_file(File, Ext), _)).

% Read a Prolog fixture file containing a clues/1 term. This reads terms rather
% than consulting the file, so benchmark/main input fixtures do not define or
% redefine global predicates.
read_clues_prolog(File, Words) :-
    setup_call_cleanup(open(File, read, S),
                       read_clues_prolog_term(S, File, Words),
                       close(S)).

read_clues_prolog_term(S, File, Words) :-
    read_term(S, Term, []),
    (   Term == end_of_file
    ->  throw(error(prolog_no_clues_term(File), _))
    ;   Term = clues(Words)
    ->  true
    ;   read_clues_prolog_term(S, File, Words)
    ).

% Read and validate a JSON clue file into the internal Words list. A missing
% file or malformed JSON throws standard ISO errors (existence_error /
% syntax_error) that SWI's default handler renders clearly, so they need no
% handling here; schema violations are caught by doc_to_words/2.
read_clues_json(File, Words) :-
    setup_call_cleanup(open(File, read, S), json_read_dict(S, Doc), close(S)),
    doc_to_words(Doc, Words).

% Map the parsed document to [[AtomAnswer, MetaDict], ...], the form clues/1
% already returns. json_read_dict does not check shape and atom_string/2
% silently coerces (a number or JSON null would become a bogus atom), so each
% type is guarded *before* conversion. Schema violations throw error/2 terms
% rendered by the prolog:error_message//1 clauses below. Answers are
% normalised to atoms so the emit-time answer_meta_assoc/2 join and
% check_unique_answers/1 (both ==-based) behave identically across sources.
doc_to_words(Doc, Words) :-
    (   is_dict(Doc), get_dict(clues, Doc, Clues), is_list(Clues)
    ->  true
    ;   throw(error(json_no_clues_array, _))
    ),
    maplist(entry_to_word, Clues, Words).

% One JSON clue entry -> [AtomAnswer, MetaDict]. `answer` is required and must
% be a string; `meta` is an optional object (default _{}), copied verbatim.
entry_to_word(Entry, [Answer, Meta]) :-
    (   is_dict(Entry), get_dict(answer, Entry, RawAnswer), string(RawAnswer)
    ->  atom_string(Answer, RawAnswer)
    ;   throw(error(json_invalid_answer(Entry), _))
    ),
    (   get_dict(meta, Entry, RawMeta)
    ->  (   is_dict(RawMeta)
        ->  Meta = RawMeta
        ;   throw(error(json_invalid_meta(Answer), _))
        )
    ;   Meta = _{}
    ).

% prolog:error_message//1 renders the formal term of an error(Formal, _). It is
% verified working on the pinned SWI 10.0.2 but is undocumented in the manual
% (which documents prolog:message//1); do not migrate without re-checking.
:- multifile prolog:error_message//1.
prolog:error_message(json_no_clues_array) -->
    [ 'clues file: expected a JSON object with a "clues" array' ].
prolog:error_message(json_invalid_answer(Entry)) -->
    [ 'clues file: every entry needs a string "answer" (offending entry: ~q)'-[Entry] ].
prolog:error_message(json_invalid_meta(Answer)) -->
    [ 'clues file: "meta" for answer ~q must be a JSON object'-[Answer] ].
prolog:error_message(unknown_strategy(S)) -->
    [ 'unknown --strategy ~q; expected one of baseline, mrv, mrv_capped, mrv_inc'-[S] ].
prolog:error_message(unsupported_clue_file(File, Ext)) -->
    [ 'clues file ~q has unsupported extension ~q; expected .json or .pl'-
      [File, Ext] ].
prolog:error_message(prolog_no_clues_term(File)) -->
    [ 'Prolog clues file ~q does not contain a clues/1 term'-[File] ].


% Top level predicate for solving the crossword with a specified
% starting position. Emits the solution as a single JSON object.
% crossword/3 uses the production default strategy (default_strategy/1);
% crossword/4 takes an explicit strategy.
crossword(GridLen, Words, StartLoc) :-
    default_strategy(Strategy),
    crossword(Strategy, GridLen, Words, StartLoc).

crossword(Strategy, GridLen, Words, StartLoc) :-
    check_unique_answers(Words),
    find_crossword(Strategy, GridLen, Words, StartLoc, _Grid, PlacedWords),
    assign_clue_numbers(PlacedWords, NumberedPlacedWords),
    emit_json(NumberedPlacedWords, Words, GridLen).


% Top level predicate for finding the number of solutions for the crossword
% for a specific starting position (or all start positions when StartLoc is
% unbound). aggregate_all/3 runs the search once and counts deterministically;
% the previous length(Sols,Num)/findall idiom re-ran the entire search once per
% candidate length and left a (failing) choicepoint behind.
all_crossword(Strategy, GridLen, Words, StartLoc, Num) :-
    aggregate_all(count,
                  find_crossword(Strategy, GridLen, Words, StartLoc, _Grid, _Placed),
                  Num).


% The driver predicate used to solve the crossword. find_crossword/5 uses the
% production default strategy (default_strategy/1); find_crossword/6 takes an
% explicit strategy. The benchmark harness calls /6 directly per strategy.
find_crossword(GridLen, Words, Loc, Grid, PlacedWords) :-
    default_strategy(Strategy),
    find_crossword(Strategy, GridLen, Words, Loc, Grid, PlacedWords).

% mrv_inc threads an incremental count cache, so it has its own driver
% (assign_words_inc/9) rather than sharing the stateless assign_words/9 path.
find_crossword(mrv_inc, GridLen, Words, Loc, Grid, PlacedWords) :-
    !,
    % Reset the per-search crossing memo (see pair_crossings/3): the tabled
    % (Val,PPos,Pos) sequences are a pure function of two answers' letters and
    % never change, so they amortize across the millions of nodes in THIS
    % search - but abolishing at the search boundary keeps each search's
    % inference count self-contained and order-independent (a warm table
    % carried across searches/workloads would make counts depend on run order,
    % breaking the benchmark's determinism invariant). The rebuild is ~|Words|^2
    % pairs, negligible against the search it serves.
    abolish_table_subgoals(pair_crossings(_,_,_)),
    abolish_table_subgoals(answer_letters(_,_)),
    init_grid(GridLen, G1),
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words_inc(Words, [], none, GridLen, StartNum, StartDir, G1, Grid, PlacedWords).
find_crossword(Strategy, GridLen, Words, Loc, Grid, PlacedWords) :-
    init_grid(GridLen, G1),
    % Get the cell number and direction for start loc
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words(Strategy, Words, [], GridLen, StartNum, StartDir, G1, Grid, PlacedWords).


% Assign all words. The starting location is selected by locating an
% intersecting word from the words already placed.
assign_words(_Strategy, [], P, _, _, _, G, G, P).
assign_words(Strategy, Words, PlacedWords, GridLen, Start, Dir, GIn, GOut, PlacedWordsOut) :-
    % select_word/9 chooses the next word to place (and the rest, RemWords)
    % according to the strategy; it is the only thing the strategies vary.
    select_word(Strategy, Words, PlacedWords, GridLen, Start, Dir, GIn, Entry, RemWords),
    Entry = [Word|_],   % the solver uses only the answer; metadata is ignored
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    % on first pass, Start and Dir will be grounded with the start values
    % then afterwards will be unground, with find_intersecting_word grounding them
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Start, Dir, GridLen, PlacedWords, GIn, Placed, G1),
    assign_words(Strategy, RemWords, [Placed|PlacedWords], GridLen, _Start, _Dir, G1, GOut, PlacedWordsOut).


% Variable ordering, pluggable per strategy. Each clause picks the next Entry
% to place and returns the remaining words; selection stays backtrackable, so
% strategies only REORDER the same search tree (completeness is unchanged).
%
% baseline: take words in input order. member/2 is the backtrack point and
% find_intersecting_word/6 (in assign_words) rejects words that cannot connect
% to the current grid - exactly the original behaviour.
select_word(baseline, Words, _Placed, _GridLen, _Start, _Dir, _GIn, Entry, RemWords) :-
    member(Entry, Words),
    remove_x(Entry, Words, RemWords).

% mrv / mrv_capped, first word (grid empty): Start/Dir are the ground seed and
% every word has the same single placement there, so ranking is moot - branch
% over all words, keeping the seed-word choice a backtrack point.
select_word(Strategy, Words, [], _GridLen, _Start, _Dir, _GIn, Entry, RemWords) :-
    mrv_strategy(Strategy),
    !,
    member(Entry, Words),
    remove_x(Entry, Words, RemWords).

% mrv / mrv_capped, otherwise: fail-first ordering. Count each word's viable
% placements on the current grid, keep only those that can connect now
% (count > 0), and offer them MOST-CONSTRAINED FIRST via member/2 (still
% backtrackable). A word's placements must cross an already-placed word, so a
% count of 0 means "not connectable yet", NOT "dead" - hence we filter rather
% than fail. The branch is dead only when NOTHING is placeable, which leaves
% Ordered empty and fails here (a sound forward check).
select_word(Strategy, Words, PlacedWords, GridLen, Start, Dir, GIn, Entry, RemWords) :-
    mrv_strategy(Strategy),
    mrv_cap(Strategy, Cap),
    map_list_to_pairs(mrv_count(Cap, PlacedWords, GridLen, Start, Dir, GIn),
                      Words, Pairs),
    include(positive_key, Pairs, Placeable),
    keysort(Placeable, Sorted),
    pairs_values(Sorted, Ordered),
    member(Entry, Ordered),
    remove_x(Entry, Words, RemWords).

mrv_strategy(mrv).
mrv_strategy(mrv_capped).

% Placement-count cap. `mrv` counts every viable placement (exact MRV);
% `mrv_capped` saturates at 2 - the ordering/forward-check only needs the
% buckets 0 / 1 / >=2, so it stops enumerating after the 2nd hit and avoids
% full MRV's per-node enumeration cost on grids where words have many slots.
mrv_cap(mrv, unbounded).
mrv_cap(mrv_capped, 2).

positive_key(Count-_) :- Count > 0.

% Count Entry's viable placements right now, bounded by Cap. Each candidate is a
% (Start,Dir) from find_intersecting_word that also survives the assign_word
% adjacency/bounds checks. Like findall, counting UNDOES the goal's bindings, so
% the caller's Start/Dir (ground on the first word, unbound after) are left
% untouched.
%
% Cap is `unbounded` (exact MRV, `mrv` strategy) or the integer 2 (mrv_capped /
% mrv_inc, which only ever need the 0 / 1 / >=2 buckets - see select_word/9 and
% assign_words_inc/9). The unbounded path keeps aggregate_all/3; the hot Cap=2
% path uses count_upto2/2, a hand-rolled saturating counter that avoids
% aggregate_all's spec-checking (error:has_type/2) and limit/2's per-solution
% nb-state machinery. The cut makes the two clauses mutually exclusive and keeps
% each call deterministic (a bare variable-cap head would leave a choicepoint).
mrv_count(unbounded, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count) :-
    !,
    mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal),
    aggregate_all(count, Goal, Count).
mrv_count(2, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count) :-
    mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal),
    count_upto2(Goal, Count).

% The placement-enumeration goal shared by both count paths: a candidate
% (Start,Dir) from find_intersecting_word that also survives assign_word's
% adjacency/bounds checks. On backtracking it re-binds Start/Dir to the next
% placement (or leaves them unbound when it fails), so the counter around it must
% undo those bindings.
mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal) :-
    Entry = [Word|_],
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    Goal = ( find_intersecting_word(Letters2, WLen, PlacedWords, GridLen,
                                    Start, Dir),
             assign_word(Word, Letters2, WLen, Start, Dir, GridLen,
                         PlacedWords, GIn, _Placed, _G1) ).

% count_upto2(:Goal, -Count): Count is the number of solutions of Goal SATURATED
% at 2 (so Count is 0, 1, or 2, where 2 means ">=2"). Like findall/aggregate_all
% it leaves NO residual bindings from Goal on exit - the search relies on the
% caller's Start/Dir being untouched.
%
% This is the SWI manual's non-backtrackable solution counter (a fresh mutable
% holder incremented per solution via nb_setarg/3 - see manipterm, succeeds_n_times)
% with an early exit: the failure-driven loop stops the moment the 2nd solution
% lands (N >= 2), so it never enumerates a word's remaining placements. Wrapping
% the loop in \+ discards the loop's bindings whether it stopped early (inner goal
% succeeded -> \+ fails) or ran dry (inner goal failed -> \+ succeeds); either way
% the nb_setarg side effect on the holder survives (it is non-backtrackable) and
% carries the count out. A fresh holder per call keeps this reentrant - no global
% key to collide if counts ever nest.
:- meta_predicate count_upto2(0, -).
count_upto2(Goal, Count) :-
    Counter = counter(0),
    (   \+ ( Goal,
             arg(1, Counter, N0),
             N is N0 + 1,
             nb_setarg(1, Counter, N),
             N >= 2 )
    ->  true
    ;   true
    ),
    arg(1, Counter, Count).


% mrv_inc: capped MRV with an INCREMENTAL count cache (idea I1 in
% docs/experiments.md). mrv_capped recomputes every remaining word's placement
% count at every node; that per-node recount is its residual cost, which grows
% with word count and grid size. mrv_inc instead recomputes a count only for
% words that can have CHANGED since the last placement, and carries the rest
% forward.
%
% Correctness rests on this invariant: a placement must cross an already-placed
% word, so placing word W can only ADD options to words that SHARE A LETTER
% with W (a new crossing target); for every other word, placing W can only
% remove options (blocking), never add. Hence a carried-forward (not recounted)
% count is always >= the true current count - never an under-count. An
% over-count is safe: such a word is merely tried when it has no real
% placement, and find_intersecting_word/6 then fails it, so member/2 moves on.
% An under-count to 0 would be unsafe (it could wrongly prune a placeable word
% and break completeness) - and the invariant guarantees that never happens.
%
% The cache is threaded as a plain argument, so backtracking restores prior
% caches automatically. Counts are capped at 2 (as in mrv_capped).

% State is `none` (no cache yet) or state(CountAssoc, LastPlacedLetters).
assign_words_inc([], P, _State, _, _, _, G, G, P).
assign_words_inc(Words, PlacedWords, StateIn, GridLen, Start, Dir, GIn, GOut, Out) :-
    select_inc(Words, PlacedWords, StateIn, GridLen, Start, Dir, GIn,
               Entry, RemWords, StateOut),
    Entry = [Word|_],
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Start, Dir, GridLen, PlacedWords, GIn, Placed, G1),
    assign_words_inc(RemWords, [Placed|PlacedWords], StateOut, GridLen,
                     _Start, _Dir, G1, GOut, Out).

% --- search perturbation (the arrange --seed knob) --------------------------
% The mrv_inc search is DETERMINISTIC by default. With no seed installed,
% seed_word_order/2 and order_candidates/2 are pure identities and no RNG
% builtin is reachable, so a default run's output is byte-identical (the golden
% invariant holds) and the RNG is never even seeded. A seed opts into a
% REPRODUCIBLE pseudo-random perturbation of the SAME search tree: it reorders
% branches, never prunes or adds them, so completeness is unchanged. The seed
% lives only on the perturbed path — it can never influence the deterministic
% one.
:- dynamic search_seed/1.

% set_search_seed(N): install the search seed. -1 (the CLI "not given"
% sentinel) clears it -> the fully deterministic path. N>=0 records the seed
% and seeds SWI's global RNG ONCE, so a given N always reproduces the same
% layout. Mirrors set_check_target/1: at most one search_seed/1 fact, retract
% before assert, the sole sanctioned writer of the module-private dynamic.
set_search_seed(-1) :- !, retractall(search_seed(_)).
set_search_seed(N) :-
    integer(N), N >= 0,
    retractall(search_seed(_)),
    set_random(seed(N)),
    assertz(search_seed(N)).

% current_search_seed(-N): the installed search seed, or FAILS if none (the
% deterministic default). Lets the emitter record the seed as provenance in the
% output without threading it through every emit signature. Semidet.
current_search_seed(N) :- search_seed(N).

% set_shuffle_seed(-N): the --shuffle path. Draw an UNPREDICTABLE seed from the
% OS entropy source (a different N on each process run -> a different layout),
% then install it as a normal search seed. Returning N keeps shuffle
% RECOVERABLE: the run is still reproducible via --seed N, so a liked layout is
% never lost. This is the only place that reseeds from entropy; --seed N and the
% deterministic default never do.
set_shuffle_seed(N) :-
    set_random(seed(random)),
    random_between(0, 1_000_000_000, N),
    set_search_seed(N).

% Seed-word order: which word anchors the layout. Deterministic = input order;
% seeded = a shuffled copy (the single biggest source of layout variety).
seed_word_order(Words, Order) :-
    ( search_seed(_) -> random_permutation(Words, Order) ; Order = Words ).

% Interior value ordering over the keysorted Count-Entry pairs (most-constrained
% first). Deterministic = the plain values; seeded = shuffled WITHIN each
% equal-count bucket only, so the fail-first heuristic (and thus solution
% quality) is preserved while equally-ranked branches are explored in a
% seed-varied order.
order_candidates(Sorted, Ordered) :-
    (   search_seed(_)
    ->  group_pairs_by_key(Sorted, Groups),
        maplist(shuffle_bucket, Groups, Buckets),
        append(Buckets, Ordered)
    ;   pairs_values(Sorted, Ordered)
    ).

shuffle_bucket(_Count-Entries, Shuffled) :- random_permutation(Entries, Shuffled).

% Seed word (grid empty): branch over all words, as the other MRV strategies
% do; the next node builds the full cache (StateOut = none).
select_inc(Words, [], _StateIn, _GridLen, _Start, _Dir, _GIn, Entry, RemWords, none) :-
    !,
    seed_word_order(Words, SeedOrder),
    member(Entry, SeedOrder),
    remove_x(Entry, Words, RemWords).
% Otherwise: get current counts (full or incremental), order most-constrained
% first over the connectable words, pick backtrackably. The cache passed
% forward is this node's count map tagged with the chosen word's letters, so
% the next node only recounts words sharing a letter with it.
select_inc(Words, PlacedWords, StateIn, GridLen, Start, Dir, GIn, Entry, RemWords, StateOut) :-
    inc_counts(StateIn, Words, PlacedWords, GridLen, Start, Dir, GIn, CountMap),
    % map_list_to_pairs keeps the ORIGINAL word terms (findall would copy them,
    % breaking the ==-based remove_x below and looping forever).
    map_list_to_pairs(count_of(CountMap), Words, Pairs),
    include(positive_key, Pairs, Placeable),
    keysort(Placeable, Sorted),
    order_candidates(Sorted, Ordered),
    member(Entry, Ordered),
    remove_x(Entry, Words, RemWords),
    entry_letters(Entry, EntryLetters),
    StateOut = state(CountMap, EntryLetters).

count_of(CountMap, Word, Count) :-
    Word = [A|_],
    get_assoc(A, CountMap, Count).

% Build the count map for the current words. With no cache yet, count them all;
% otherwise recount only words sharing a letter with the last-placed word and
% carry the rest forward from the previous map.
inc_counts(none, Words, PlacedWords, GridLen, Start, Dir, GIn, CountMap) :-
    !,
    empty_assoc(A0),
    foldl(full_count(PlacedWords, GridLen, Start, Dir, GIn), Words, A0, CountMap).
inc_counts(state(PrevMap, LastLetters), Words, PlacedWords, GridLen, Start, Dir, GIn, CountMap) :-
    empty_assoc(A0),
    foldl(inc_count(PrevMap, LastLetters, PlacedWords, GridLen, Start, Dir, GIn),
          Words, A0, CountMap).

full_count(PlacedWords, GridLen, Start, Dir, GIn, Entry, AIn, AOut) :-
    Entry = [A|_],
    mrv_count(2, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count),
    put_assoc(A, AIn, Count, AOut).

inc_count(PrevMap, LastLetters, PlacedWords, GridLen, Start, Dir, GIn, Entry, AIn, AOut) :-
    Entry = [A|_],
    entry_letters(Entry, ELetters),
    (   shares_letter(ELetters, LastLetters)
    ->  mrv_count(2, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count)
    ;   get_assoc(A, PrevMap, Count)            % carry forward (>= true count)
    ->  true
    ;   mrv_count(2, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count)
    ),
    put_assoc(A, AIn, Count, AOut).

entry_letters([Word|_], Letters) :-
    answer_letters(Word, Letters).

% answer_letters(+Word, -Letters): the placement footprint (letters with space
% and hyphen separators dropped) of a single answer atom. Word is ground, so
% this is a pure function of the atom and is tabled: the mrv_inc counting path
% (mrv_count/select_inc/inc_count) re-derives a word's letters per node, and the
% atom_chars + two delete/3 passes are memoized to one lookup after the first.
% Reset per search alongside pair_crossings/3 (find_crossword), so the inference
% count stays self-contained across benchmark iterations.
:- table answer_letters/2.
answer_letters(Word, Letters) :-
    atom_chars(Word, L0),
    delete(L0, ' ', L1),
    delete(L1, '-', Letters).

shares_letter(Letters, OtherLetters) :-
    member(L, Letters),
    memberchk(L, OtherLetters),
    !.


% Given a Word and a set of Placed words, locates a candidate
% start cell and direction that intersects with an existing word.

% No placed words; just use grounded Start and Dir values
find_intersecting_word(_Letters, _WLen, [], _GridLen, _Start, _Dir).

find_intersecting_word(Letters, WLen, PlacedWords, GridLen, Start, Dir) :-
    member(PW, PlacedWords),
    get_dict(letters, PW, PLetters),
    get_dict(dir, PW, PDir),
    get_dict(start, PW, PStart),
    % The (PPos,Pos) crossing sequence for this (new letters, placed letters)
    % pair is grid-independent and unchanging, so it is memoized once by
    % pair_crossings/3 (tabled) and replayed here. member/2 walks the precomputed
    % list in the SAME order the former inline intersection/list_to_set/nth1/nth1
    % conjunction produced (HARD requirement: order is golden-visible).
    pair_crossings(Letters, PLetters, Crossings),
    member(x(PPos, Pos), Crossings),
    calc_num(PDir, GridLen, PPos, PStart, PNum),
    swap_dir(PDir, Dir),
    calc_start(Dir, GridLen, Pos, PNum, Start),
    fits_on_grid(Dir, Start, WLen, GridLen).

% pair_crossings(+Letters, +PLetters, -Crossings): the ordered list of x(PPos,Pos)
% cells at which a NEW word (letters Letters) can cross an already-placed word
% (letters PLetters) - PPos the crossing position in the placed word, Pos the
% position in the new word. Both letter-lists are GROUND (derived from the ground
% answer atom; the non-ground part of a [Answer,_{...}] entry never reaches here),
% so this is a pure function of the two letter sequences and is memoized with
% tabling: a search reuses each pair across all its nodes instead of re-deriving
% intersection/3 + list_to_set/2 + two nth1/3 backtracking passes per call.
%
% ORDER IS LOAD-BEARING. The findall reproduces the former inline conjunction
% verbatim: intersection/3 keeps Letters-order-with-duplicates, list_to_set/2
% dedupes to the new word's distinct crossing letters in first-occurrence order,
% then for each such Val, PPos ascends over the placed word and Pos ascends over
% the new word. Replaying the collected list preserves first (and all) solutions
% byte-for-byte. Val itself is not retained - the caller only needs the two
% positions. Tabling is supported under the SWI WASM build.
:- table pair_crossings/3.
pair_crossings(Letters, PLetters, Crossings) :-
    findall(x(PPos, Pos),
            ( intersection(Letters, PLetters, Vals),
              list_to_set(Vals, Vals2),
              member(Val, Vals2),
              nth1(PPos, PLetters, Val),
              nth1(Pos, Letters, Val) ),
            Crossings).


assign_word(Word, Letters, WLen, Start, Dir, GridLen, PlacedWords, GIn, Placed, GOut) :-
    % Reject an off-grid (underflow) start cell. find_intersecting_word/6 can
    % compute a Start < 1 that still passes fits_on_grid (the end lands on-grid),
    % which the old assoc grid rejected implicitly: get_assoc/3 FAILED for the
    % absent negative key. arg/3 instead THROWS for a negative index (0 and
    % >arity still fail, as get_assoc did), so we reject Start < 1 here, once per
    % word. With Start >= 1 and a valid run, every downstream cell index is >= 1,
    % so no read can go negative (an over-grid index fails, never throws).
    Start >= 1,
    % make sure previous cell does not have a letter
    check_prev_cell(Dir, Start, GridLen, GIn),
    assign_letters(Letters, Start, Dir, GridLen, Cells, GIn, GOut),
    % keep every word a maximal run: this word must not fill the boundary cell
    % of an already-placed word (see no_word_merge/3 and docs/experiments.md I5)
    no_word_merge(Cells, GridLen, PlacedWords),
    % `num` is added later by assign_clue_numbers/2
    Placed = word{answer:Word, letters:Letters, cells:Cells,
                  dir:Dir, len:WLen, start:Start}.


% A valid crossword needs every word to be a MAXIMAL run of cells - bounded on
% both ends, in its own direction, by an empty cell or the grid edge.
% check_prev_cell/4 and check_next_cell/4 enforce that only at the moment a word
% is placed; a later, longer collinear word can pass straight through a shorter
% word's cells (the X==L branch of assign_letters/7 skips adjacency checks) and
% fill the cell just past its end, retroactively making the shorter word a
% "word inside a word". That produces two same-direction answers sharing a start
% cell, which assign_clue_numbers/2 cannot number. We forbid it here: the word
% being placed (cells Cells) must not occupy the boundary cell of any
% already-placed word. (The reverse - an existing word occupying THIS word's
% boundary - is already caught by check_prev_cell/check_next_cell, so together
% the invariant is maintained inductively.)
no_word_merge(Cells, GridLen, PlacedWords) :-
    \+ ( member(PW, PlacedWords),
         placed_boundary_cell(PW, GridLen, Boundary),
         memberchk(Boundary, Cells) ).

% Each on-grid boundary cell of a placed word: the cell just before its start,
% and the cell just after its end, in the word's own direction (omitting either
% when the word abuts the grid edge there).
placed_boundary_cell(PW, GridLen, Before) :-
    get_dict(dir, PW, Dir),
    get_dict(start, PW, Start),
    \+ is_start_cell(Dir, Start, GridLen),
    prev_cell(Dir, Start, GridLen, Before).
placed_boundary_cell(PW, GridLen, After) :-
    get_dict(dir, PW, Dir),
    get_dict(cells, PW, Cells),
    last(Cells, End),
    \+ is_end_cell(Dir, End, GridLen),
    next_cell(Dir, End, GridLen, After).


% Previous cell before start of word. Make sure it doesn't contain
% anything.
check_prev_cell(Dir, Num, GridLen, Grid) :-
    (
     % don't check if start letter is start of a row/col
     is_start_cell(Dir, Num, GridLen)
    ;
     % otherwise prev cell must be empty (an unbound var)
     prev_cell(Dir, Num, GridLen, Prev),
     arg(Prev, Grid, Cell), var(Cell)
    ), !.


% Next cell after end of word. Make sure it doesn't contain anything.
check_next_cell(Dir, Num, GridLen, Grid) :-
    prev_cell(Dir, Num, GridLen, Prev),
    (
     % no need to check if prev was end of row/col
     is_end_cell(Dir, Prev, GridLen)
    ;
     % then this cell must be empty (an unbound var)
     arg(Num, Grid, Cell), var(Cell)
    ), !.


% Assign each letter of the word, checking that adjacent cells are empty to
% prevent words being placed next to each other. The grid is mutated in place by
% binding empty cell vars, so GIn and GOut are the SAME term (backtracking
% unbinds via the trail); they thread as one.

% Last letter of word, make sure next cell is free
assign_letters([], Num, Dir, GridLen, [], Grid, Grid) :-
    check_next_cell(Dir, Num, GridLen, Grid).


assign_letters([L|Ls], Num, Dir, GridLen, [Num|RestCells], Grid, GridOut) :-
    arg(Num, Grid, Cell),
    (
     % existing letter in this cell matches letter being placed,
     % nothing needs doing, we can continue to next letter
     nonvar(Cell),
     Cell == L
    ;
     % empty cell (an unbound var), so check adjacent cells are free
     % and then bind this cell to the letter
     var(Cell),
     adj_is_free(Dir, Num, GridLen, Grid),
     Cell = L
    ), !,
    next_cell(Dir, Num, GridLen, Num2),
    assign_letters(Ls, Num2, Dir, GridLen, RestCells, Grid, GridOut).


% check that adjacent cells are empty (unbound vars). arg/3 + var/1 is a pure
% test - it never binds a cell (aliasing a fresh local to an empty cell's var
% leaves both unbound), so a CHECK can never accidentally fill a cell.
adj_is_free(down, Num, GridLen, Grid) :-
    N1 is Num - 1,
    N2 is Num + 1,
    M is Num mod GridLen,
    (
     M == 0 -> % last cell in row
     arg(N1, Grid, C1), var(C1)
    ;
     M == 1 -> % first cell in row
     arg(N2, Grid, C2), var(C2)
    ;
     arg(N1, Grid, C1), var(C1),
     arg(N2, Grid, C2), var(C2)
    ).

adj_is_free(across, Num, GridLen, Grid) :-
    N1 is Num - GridLen,
    N2 is Num + GridLen,
    LastCell is (GridLen * GridLen),
    (
     N1 =< 0 -> % before beginning of grid
     arg(N2, Grid, C2), var(C2)
    ;
     N2 > LastCell -> % after end of grid
     arg(N1, Grid, C1), var(C1)
    ;
     arg(N1, Grid, C1), var(C1),
     arg(N2, Grid, C2), var(C2)
    ).


% Takes the placed words and works out the clue numbers of each
% clue/word. This works by interating through the list of words which
% are sorted by their start numbers.  Clue numbers are then assigned
% in the order that words occur in this list.
assign_clue_numbers(PlacedWords, WordsClues) :-
    % sort placed words by their start number
    map_list_to_pairs(start_is, PlacedWords, Pairs),
    keysort(Pairs, SortedPairs),

    % group clues that have the same start num together (ie a cell
    % that is start of both across and down clue) as these will
    % recieve the same clue number
    group_pairs_by_key(SortedPairs, GroupedPairs),

    % assign the clue numbers, starting at 1
    add_clue_nums(GroupedPairs, 1, WordsClues).


% Updates the placed words list by appending the clue number to the end
% of each word. The list now looks like this:
% placed words -- [word, letters, cells, dir, len, start, clue_num]
%
% Invariant: each key-group (words sharing a start cell) holds 1 or 2 words - a
% cell begins at most one across AND one down word, never 3+. So there are
% deliberately only [W] and [W1,W2] clauses; a >2 group is structurally
% impossible and (by design) has no clause (P17).
add_clue_nums([], _, []).

% a word whose start cell only belongs to a down or an across word
add_clue_nums([_-[W]|Rest], ClueNum, [WClue|RestClues]) :-    
    add_clue_word(W, ClueNum, WClue),
    ClueNum2 is ClueNum + 1,
    add_clue_nums(Rest, ClueNum2, RestClues).    
    
% a word whose start cell belongs to both a down and an across word
add_clue_nums([_-[W1,W2]|Rest], ClueNum, [WClue1,WClue2|RestClues]) :-
    add_clue_word(W1, ClueNum, WClue1),
    add_clue_word(W2, ClueNum, WClue2),
    ClueNum2 is ClueNum + 1,
    add_clue_nums(Rest, ClueNum2, RestClues).    


% JSON output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The emit-time metadata join is keyed on the answer string, so fail loudly if
% two input entries share an answer (crosswords do not repeat answers). A
% standard error/2 term plus the error_message//1 hook below renders cleanly;
% a bare custom throw would print as "Unknown message".
check_unique_answers(Words) :-
    findall(A, member([A|_], Words), Answers),
    msort(Answers, Sorted),
    (
     append(_, [Dup, Dup|_], Sorted)
    ->
     throw(error(duplicate_answer(Dup), _))
    ;
     true
    ).

prolog:error_message(duplicate_answer(Answer)) -->
    [ 'duplicate answer ~q in clues; answers must be unique'-[Answer] ].


% Emit the solved crossword as a single JSON object on the current output.
% Words is the original input list, used to rejoin per-word metadata.
emit_json(NumberedPlacedWords, Words, GridLen) :-
    build_grid_rows(NumberedPlacedWords, GridLen, GridRows),
    build_words(NumberedPlacedWords, Words, GridLen, WordObjs),
    Payload = _{gridLength: GridLen, grid: GridRows, words: WordObjs},
    current_output(Out),
    json_write_dict(Out, Payload),
    nl(Out).


% Build the dense GridLen x GridLen array, row-major. Empty cells are the atom
% `null` (written as JSON null); filled cells are dicts.
build_grid_rows(PlacedWords, GridLen, Rows) :-
    empty_assoc(A0),
    foldl(add_word_cells, PlacedWords, A0, CellMap),
    NumCells is GridLen * GridLen,
    numlist(1, NumCells, AllCells),
    maplist(cell_to_json(CellMap), AllCells, FlatCells),
    rows_of(FlatCells, GridLen, Rows).


% Fold one placed word's cells into the cell map.
add_word_cells(PW, AIn, AOut) :-
    get_dict(cells, PW, Cells),
    get_dict(letters, PW, Letters),
    get_dict(dir, PW, Dir),
    get_dict(num, PW, Num),
    get_dict(start, PW, Start),
    foldl(add_cell(Dir, Num, Start), Cells, Letters, AIn, AOut).


% Record one cell's letter, its across/down clue numbers, and (if it is the
% start cell of a word) its corner-label number.
add_cell(Dir, Num, Start, Cell, Letter, AIn, AOut) :-
    (
     get_assoc(Cell, AIn, cell(_, Ac0, Dn0, N0))
    ->
     true
    ;
     Ac0 = null, Dn0 = null, N0 = null
    ),
    ( Dir == across -> Ac = Num, Dn = Dn0 ; Ac = Ac0, Dn = Num ),
    ( Cell =:= Start -> N = Num ; N = N0 ),
    put_assoc(Cell, AIn, cell(Letter, Ac, Dn, N), AOut).


cell_to_json(CellMap, Cell, Json) :-
    (
     get_assoc(Cell, CellMap, cell(Letter, Ac, Dn, N))
    ->
     Json = _{letter: Letter, number: N, across: Ac, down: Dn}
    ;
     Json = null
    ).


% Group a flat list of cells into rows of length GridLen.
rows_of([], _, []) :- !.
rows_of(Flat, GridLen, [Row|Rows]) :-
    length(Row, GridLen),
    append(Row, Rest, Flat),
    rows_of(Rest, GridLen, Rows).


% Build the `words` array. Metadata is rejoined from the input list by answer -
% via a one-shot answer->meta assoc so the join is O(n log n), not a member/2
% rescan of the whole input per placed word (P5).
build_words(PlacedWords, Words, GridLen, WordObjs) :-
    canonical_word_order(PlacedWords, Ordered),
    answer_meta_assoc(Words, MetaAssoc),
    maplist(placed_to_word(MetaAssoc, GridLen), Ordered, WordObjs).

% answer -> metadata dict for every input entry. An entry may omit metadata
% ([Answer] -> {}); answers are unique (check_unique_answers/1), so no key clash.
answer_meta_assoc(Words, Assoc) :-
    findall(A-Meta,
            ( member(Entry, Words), Entry = [A|_],
              ( Entry = [_, M] -> Meta = M ; Meta = _{} ) ),
            Pairs),
    list_to_assoc(Pairs, Assoc).

% Canonical emit order for the `words` array: by clue number, then across before
% down (atom order 'across' @< 'down'). Two words sharing a start cell (same
% number, one across + one down) are otherwise emitted in placement order, which
% differs between the strict DFS and a fragment re-solve - breaking the
% emit -> re-ingest-as-fragment byte-identity (AC-EMIT-2, R3). num+dir is unique
% per word, so this is a total, stable order. The grid rows are built from cell
% occupancy and are unaffected. map_list_to_pairs (not findall) keeps the
% original word dicts (no copy).
canonical_word_order(PlacedWords, Ordered) :-
    map_list_to_pairs(word_canon_key, PlacedWords, Keyed),
    keysort(Keyed, Sorted),
    pairs_values(Sorted, Ordered).

word_canon_key(PW, Num-Dir) :- get_dict(num, PW, Num), get_dict(dir, PW, Dir).


placed_to_word(MetaAssoc, GridLen, PW, WordObj) :-
    get_dict(answer, PW, Answer),
    get_dict(dir, PW, Dir),
    get_dict(num, PW, Num),
    get_dict(cells, PW, Cells),
    maplist(cell_coord(GridLen), Cells, Coords),
    get_assoc(Answer, MetaAssoc, Meta),
    WordObj = _{number: Num, direction: Dir, answer: Answer,
                cells: Coords, meta: Meta}.


% Cell number (1-based, row-major) to 0-based [Row, Col].
cell_coord(GridLen, Cell, [Row, Col]) :-
    Row is (Cell - 1) // GridLen,
    Col is (Cell - 1) mod GridLen.



% Placed word utility predicates:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A placed word is a dict:
%   word{answer:Answer, letters:Letters, cells:Cells, dir:Dir,
%        len:Len, start:Start, num:ClueNum}
% `num` is absent until assign_clue_numbers/2 fills it in.


start_is(PW, Start) :- get_dict(start, PW, Start).

% Add the assigned clue number to a placed word.
add_clue_word(PW, ClueNum, NumberedPW) :-
    put_dict(num, PW, ClueNum, NumberedPW).



% grid utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

swap_dir(down, across).
swap_dir(across, down).


% The list of start locations.
start_locs([topleft_across, topleft_down, topright, bottomleft]).

% Get the cell number of each possible start location. 
start_loc(topleft_across, _GridLen, 1, across).
start_loc(topleft_down, _GridLen, 1, down).
start_loc(topright, GridLen, GridLen, down).
start_loc(bottomleft, GridLen, StartNum, across) :- 
    StartNum is (GridLen * GridLen) - (GridLen - 1).


% make sure word fits in the row...
fits_on_grid(across, Start, WLen, GridLen) :- 
    M is Start mod GridLen,
    M \== 0,
    Space is GridLen - (M - 1),
    WLen =< Space.

% Make sure word fits in the column...
fits_on_grid(down, Start, WLen, GridLen) :- 
     EndNum is Start + (GridLen * (WLen - 1)),
     EndNum =< GridLen * GridLen.


is_start_cell(across, Num, Length) :- 1 is Num mod Length.
is_start_cell(down, Num, Length) :- Num =< Length.
is_end_cell(across, Num, Length) :- 0 is Num mod Length.
is_end_cell(down, Num, Length) :- Num > (Length - 1) * Length.


prev_cell(across, Num, _Length, Prev) :- Prev is Num - 1.
prev_cell(down, Num, Length, Prev) :- Prev is Num - Length.


next_cell(across, Num, _Length, Next) :- Next is Num + 1.
next_cell(down, Num, Length, Next) :- Next is Num + Length.


% Calculates the start of a word given the number
% of a letter that occurs at a certain position
calc_start(across, _GridLen, PPos, WNum, Start) :-
    Start is WNum - (PPos - 1).

calc_start(down, GridLen, PPos, WNum, Start) :-
    Start is WNum - (GridLen * (PPos - 1)).


% Calulates the number of a letter in a word given
% the position in the word and the number it starts at
calc_num(across, _GridLen, WPos, WStart, WNum) :-
    WNum is  WStart + (WPos - 1).

calc_num(down, GridLen, WPos, WStart, WNum) :-
    WNum is  WStart + (GridLen * (WPos - 1)).


% The grid is the compound term grid(C1,...,C(N*N)); cell Num is arg Num, an
% unbound var when empty. functor/3 builds it with N*N fresh variables in one
% allocation (441 args at 21x21 is fine). See the data-structure note above.
init_grid(GridLen, Grid) :-
    NumTiles is GridLen * GridLen,
    functor(Grid, grid, NumTiles).



% (the solved crossword is emitted as JSON by emit_json/3; see the
% "JSON output" section above)



% generic utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove_x(X,L,R) :- R is L with first occurrence of X removed from it.
remove_x(Y,[X|Xs],[X|Tail]) :-
	Y \== X,
	remove_x(Y,Xs,Tail).
remove_x(X,[X|Xs],Xs) :- !.
remove_x(_,[],[]).
