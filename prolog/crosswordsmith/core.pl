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
            % the file-free JSON-document input seam (browser.pl's params path)
            doc_to_words/2,
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
            init_gs/2,
            start_loc/4,
            start_locs/1,
            valid_loc/1,
            next_cell/4,
            fits_on_grid/4,
            % search primitives
            assign_word/9,
            check_word_fits/5,
            find_intersecting_word/6,
            assign_words_inc/9,
            find_crossword/6,
            all_crossword/5,
            % the memo-hygiene seam every top-level search entry runs
            % (arrange.pl's entries are the out-of-module consumers)
            reset_search_memos/0,
            % strategy registry
            strategies/1,
            default_strategy/1,
            valid_strategy/1,
            require_strategy/1,
            % placed-word record (pw/8) accessors
            pw_answer/2,
            pw_letters/2,
            pw_cells/2,
            pw_dir/2,
            pw_len/2,
            pw_start/2,
            pw_end/2,
            pw_num/2,
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
% being a fixed-arity compound record pw/8 of the word's attributes. See
% the placed words utility predicate section below for its structure.
%
% The second data structure holds the contents of the grid. It is a single
% compound term grid(C1, ..., C(GridLen*GridLen)) whose Num-th argument is cell
% Num (1 through GridLen*GridLen, row-major): an UNBOUND variable when the cell
% is empty, bound to a letter atom when filled. Reads are arg/3 + var/nonvar
% (O(1), allocation-free); a placement unifies the cell variable with its
% letter and is undone automatically by the trail on backtracking, so a grid's
% GIn and GOut are one and the same term (the old immutable-version threading
% collapses onto the trail). See docs/experiments.md (E-H2).

% All library imports below carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5).

% Used by the JSON output emitter — library(json), NOT the legacy
% library(http/json) alias: the alias does not resolve in the WASM image
% (C6; wasm/README.md "Browser gotchas"). NB the swap moves the CLI fill's
% one-off index build across a global-stack growth threshold (one extra stack
% SHIFT, a transient ~30 MB peak-RSS blip on the ENABLE-scale build); live
% usage (statistics globalused) and all gated inference counts are identical.
:- use_module(library(json), [json_read_dict/2, json_write_dict/2]).

% library(random) serves only set_shuffle_seed/1's entropy draw (which WANTS
% unpredictability). The seeded search deliberately does NOT use the VM RNG:
% set_random/
% random_permutation delegate to whatever backend the build linked (GMP's
% generator natively, SWI's builtin under the USE_GMP=OFF wasm build), so the
% same seed would produce different layouts on the CLI vs the browser. The
% seeded path draws from a module-owned portable PRNG instead (splitmix64,
% below). Neither is reached on the deterministic path (search_seed/1 unset ->
% identity), so a default run never draws from — nor even seeds — any RNG.
% (set_random/1 itself is a system builtin, not a library import.)
:- use_module(library(random), [random_between/3]).

% aggregate_all/3, used to count solutions in all_crossword/5.
:- use_module(library(aggregate), [aggregate_all/3]).

% The list/assoc/pairs workhorses (previously autoload-supplied).
:- use_module(library(apply),
              [exclude/3, foldl/4, foldl/5, include/3, maplist/3]).
:- use_module(library(assoc),
              [empty_assoc/1, get_assoc/3, list_to_assoc/2, put_assoc/4]).
:- use_module(library(lists),
              [ append/2, append/3, intersection/3, last/2, list_to_set/2,
                member/2, nth1/3, numlist/3 ]).
:- use_module(library(pairs),
              [group_pairs_by_key/2, map_list_to_pairs/3, pairs_values/2]).

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
% The expected-strategies list is rendered from strategies/1 at message time
% (via {}/1) so the registry stays the single source of truth: adding a
% strategy keeps this diagnostic current with no edit here.
prolog:error_message(unknown_strategy(S)) -->
    { strategies(Ss),
      atomic_list_concat(Ss, ', ', List)
    },
    [ 'unknown --strategy ~q; expected one of ~w'-[S, List] ].
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


% reset_search_memos: abolish the two tabled search memos (pair_crossings/3,
% answer_letters/2) so the next search starts from empty tables.
%
% THE SEAM INVARIANT (audit C1/C48): every TOP-LEVEL search entry runs this
% exactly ONCE on entry, all strategies and paths alike - find_crossword/6
% below (which crossword/3,4 and all_crossword/5 route through; the benchmark
% samplers and tests call it directly) and arrange.pl's entry predicates
% (arrange_best_layout/6, arrange_best_effort/6, arrange_candidates/6,
% arrange_fragment_strict/6, arrange_fragment_best_effort/7). That restores
% two properties at once: (a) in-process inference counts are
% HISTORY-INDEPENDENT - a search's count never depends on which searches ran
% before it in the same process (the bench ratchet's determinism invariant;
% previously only the mrv_inc clause abolished, so baseline/greedy counts
% were call-order-dependent), and (b) table growth in a persistent process
% (WASM worker, native SDK embedding, this test suite) is BOUNDED to one
% request's residue, flushed by the next entry's reset, instead of one memo
% variant per (Letters, PLetters) pair ever seen. WITHIN one entry the memos
% are deliberately shared (corner sweeps, candidate pools, greedy per-move
% scoring): never call this from a per-corner/per-candidate loop.
%
% abolish_module_tables/1 "removes all tables that belong to predicates in
% Module" (tabling-preds.md); this module's tabled predicates are exactly the
% two memos, so the scope is precise - if a table with a DIFFERENT lifecycle
% ever lands in this module, revisit this call. It is chosen over two
% abolish_table_subgoals/1 calls on probe evidence (SWI-Prolog 10.1.10): the
% subgoal form's flush cost VARIES with the residue being flushed (70/74/16
% inferences after a strict/greedy/no residue), leaking a +-4-inference
% history signature into every entry's count, while abolish_module_tables/1
% costs a CONSTANT 67 inferences for any non-empty residue (9 when already
% empty - only ever a process's very first search, whose count carries
% one-time JIT/autoload noise anyway), making steady-state entry counts
% exactly history-independent. Both forms were probe-verified to flush all
% variants; scoped abolition is preferred over abolish_all_tables/0 so any
% future table outside this module keeps its own lifecycle. (Test-side trap:
% current_table/2 VARIANT-matches its goal argument, so an open pattern
% silently enumerates nothing - census code must enumerate unbound and
% filter, never pattern-match.) Cost of the flush: the rebuild is ~|Words|^2
% pairs, negligible against the search it serves; the tabled values are pure
% functions of the answers' letters, so abolition can never change results -
% only memory and inference counts.
% Known measurement floor (probe-verified, irreducible from Prolog): the
% first tabled call after a reset costs +-2 inferences depending on the SHAPE
% of the residue the reset destroyed (e.g. whether the previous request ever
% populated answer_letters/2 - the greedy path does not). So an entry's total
% count is exact-deterministic given (its input, the predecessor entry's
% type), and cross-history counts agree to within 4 inferences (~0.002%);
% the memo WORK itself is exactly history-independent (decomposition probe:
% identical search counts after a manual pre-reset). Canonicalizing the
% post-abolish state (dummy seed variants, double abolish) was probed and
% does NOT remove this - the variance just moves into the seeding calls.
reset_search_memos :-
    abolish_module_tables(crosswordsmith_core).

% The driver predicate used to solve the crossword. find_crossword/5 uses the
% production default strategy (default_strategy/1); find_crossword/6 takes an
% explicit strategy. The benchmark harness calls /6 directly per strategy.
find_crossword(GridLen, Words, Loc, Grid, PlacedWords) :-
    default_strategy(Strategy),
    find_crossword(Strategy, GridLen, Words, Loc, Grid, PlacedWords).

% find_crossword/6 is a TOP-LEVEL search entry, so it owns a
% reset_search_memos/0 seam: exactly one memo reset per external call, ALL
% strategies alike (see reset_search_memos/0 above). arrange.pl's strict
% corner sweep deliberately does NOT come through here - it calls the mrv_inc
% driver (init_gs/start_loc/assign_words_inc) directly under its own
% entry-level reset, so one arrange request's corners share the memos.
find_crossword(Strategy, GridLen, Words, Loc, Grid, PlacedWords) :-
    reset_search_memos,
    find_crossword_strategy(Strategy, GridLen, Words, Loc, Grid, PlacedWords).

% mrv_inc threads an incremental count cache, so it has its own driver
% (assign_words_inc/9) rather than sharing the stateless assign_words/9 path.
% The strategies are ENUMERATED clauses (one constant head per strategy, no
% var catch-all, no cut): first-argument indexing selects deterministically,
% and an unknown strategy FAILS here rather than falling through to the
% stateless path - require_strategy/1 is the seam that turns an unknown
% strategy atom into an error, and every CLI/API entry runs it first.
find_crossword_strategy(mrv_inc, GridLen, Words, Loc, Grid, PlacedWords) :-
    init_gs(GridLen, G1),
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words_inc(Words, [], none, GridLen, StartNum, StartDir, G1, Grid, PlacedWords).
find_crossword_strategy(baseline, GridLen, Words, Loc, Grid, PlacedWords) :-
    find_crossword_stateless(baseline, GridLen, Words, Loc, Grid, PlacedWords).
find_crossword_strategy(mrv, GridLen, Words, Loc, Grid, PlacedWords) :-
    find_crossword_stateless(mrv, GridLen, Words, Loc, Grid, PlacedWords).
find_crossword_strategy(mrv_capped, GridLen, Words, Loc, Grid, PlacedWords) :-
    find_crossword_stateless(mrv_capped, GridLen, Words, Loc, Grid, PlacedWords).

% Shared driver for the stateless (non-incremental) strategies.
find_crossword_stateless(Strategy, GridLen, Words, Loc, Grid, PlacedWords) :-
    init_gs(GridLen, G1),
    % Get the cell number and direction for start loc
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words(Strategy, Words, [], GridLen, StartNum, StartDir, G1, Grid, PlacedWords).


% Assign all words. The starting location is selected by locating an
% intersecting word from the words already placed.
% The [_|_] head keeps the clauses index-disjoint on arg 2 ([] vs cons): a
% terminal-[] call exits without leaving a choicepoint into this clause.
assign_words(_Strategy, [], P, _, _, _, G, G, P).
assign_words(Strategy, [W|Ws], PlacedWords, GridLen, Start, Dir, GIn, GOut, PlacedWordsOut) :-
    Words = [W|Ws],
    % select_word/9 chooses the next word to place (and the rest, RemWords)
    % according to the strategy; it is the only thing the strategies vary.
    select_word(Strategy, Words, PlacedWords, GridLen, Start, Dir, GIn, Entry, RemWords),
    Entry = [Word|_],   % the solver uses only the answer; metadata is ignored
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    % on first pass, Start and Dir will be grounded with the start values
    % then afterwards will be unground, with find_intersecting_word grounding them
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Start, Dir, GridLen, GIn, Placed, G1),
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
% nb-state machinery. The two clauses dispatch on distinct first-arg constants
% (`unbounded` vs `2`), so clause indexing keeps each call deterministic - no
% cut needed.
mrv_count(unbounded, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count) :-
    mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal),
    aggregate_all(count, Goal, Count).
mrv_count(2, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count) :-
    mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal),
    count_upto2(Goal, Count).

% The placement-enumeration goal shared by both count paths: a candidate
% (Start,Dir) from find_intersecting_word that also survives the legality
% checks. On backtracking it re-binds Start/Dir to the next placement (or
% leaves them unbound when it fails), so the counter around it must undo those
% bindings.
%
% Counting only needs to know whether a candidate placement is LEGAL, not to
% realize it: check_word_fits/5 runs the SAME checks as assign_word/9 but
% binds no grid cell and builds neither the Cells list nor the pw/8 record
% (all of which assign_word materialized only for count_upto2 to discard).
% Same accept/reject set in the same enumeration order -> every count is
% identical and the search tree/output/goldens stay byte-identical.
mrv_count_goal(PlacedWords, GridLen, Start, Dir, GIn, Entry, Goal) :-
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    Goal = ( find_intersecting_word(Letters2, WLen, PlacedWords, GridLen,
                                    Start, Dir),
             check_word_fits(Letters2, Start, Dir, GridLen, GIn) ).

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
% The [_|_] head keeps the clauses index-disjoint on arg 1 ([] vs cons): a
% terminal-[] call exits without leaving a choicepoint into this clause.
assign_words_inc([], P, _State, _, _, _, G, G, P).
assign_words_inc([W|Ws], PlacedWords, StateIn, GridLen, Start, Dir, GIn, GOut, Out) :-
    Words = [W|Ws],
    select_inc(Words, PlacedWords, StateIn, GridLen, Start, Dir, GIn,
               Entry, RemWords, StateOut),
    Entry = [Word|_],
    entry_letters(Entry, Letters2),
    length(Letters2, WLen),
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Start, Dir, GridLen, GIn, Placed, G1),
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

% prng_state/1: the advancing state of the module-owned PRNG; exists iff
% search_seed/1 does. Like the global RNG it replaced, it is NON-backtrackable
% (retract/assert survives backtracking into the search) and shared with
% {engine:true} search engines (dynamics live in the global database).
:- dynamic prng_state/1.

% set_search_seed(N): install the search seed. -1 (the CLI "not given"
% sentinel) clears it -> the fully deterministic path. N>=0 records the seed
% and seeds the module's PORTABLE PRNG (splitmix64) ONCE, so a given N always
% reproduces the same layout — on every build of every platform. We own the
% algorithm precisely so reproducibility is NOT engine-build-scoped: SWI's
% set_random/1 seeds GMP's generator natively but SWI's builtin under the
% USE_GMP=OFF wasm build, which made the same seed diverge CLI vs browser
% (finding 2026-07-06, wasm-sdk-strategy §10). Mirrors set_check_target/1: at
% most one search_seed/1 fact, retract before assert, the sole sanctioned
% writer of the two module-private dynamics.
set_search_seed(-1) :- !,
    retractall(search_seed(_)),
    retractall(prng_state(_)).
set_search_seed(N) :-
    integer(N), N >= 0,
    retractall(search_seed(_)),
    retractall(prng_state(_)),
    S0 is N /\ 0xFFFFFFFFFFFFFFFF,
    assertz(prng_state(S0)),
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

% The portable PRNG: splitmix64 (public-domain algorithm; reference vectors
% locked in tests/arrange.plt). Pure unbounded-integer arithmetic masked to 64
% bits, so it is bit-identical under GMP and LibBF — the same guarantee the
% inference-parity certification rests on.
splitmix64(S0, V, S1) :-
    S1 is (S0 + 0x9E3779B97F4A7C15) /\ 0xFFFFFFFFFFFFFFFF,
    Z1 is ((S1 xor (S1 >> 30)) * 0xBF58476D1CE4E5B9) /\ 0xFFFFFFFFFFFFFFFF,
    Z2 is ((Z1 xor (Z1 >> 27)) * 0x94D049BB133111EB) /\ 0xFFFFFFFFFFFFFFFF,
    V is Z2 xor (Z2 >> 31).

% Draw one value and advance the state destructively (mirrors the global RNG's
% non-backtracking stream). Throws if no seed is installed — every caller is
% guarded by search_seed(_), so reaching this seedless is a program error.
prng_draw(V) :-
    (   retract(prng_state(S0))
    ->  splitmix64(S0, V, S1),
        assertz(prng_state(S1))
    ;   throw(error(existence_error(prng_state, unseeded), _))
    ).

% seeded_permutation(+List, -Perm): selection shuffle driven by prng_draw/1 —
% one draw per element, index by V mod remaining-length. The tiny modulo bias
% (~len/2^64) is irrelevant here: the contract is reproducibility, not
% statistical uniformity. O(n^2), fine at word-list scale.
seeded_permutation([], []) :- !.
seeded_permutation(List, [X|Perm]) :-
    length(List, N),
    prng_draw(V),
    I is V mod N,
    list_nth0_rest(I, List, X, Rest),
    seeded_permutation(Rest, Perm).

list_nth0_rest(0, [X|Rest], X, Rest) :- !.
list_nth0_rest(I, [H|T], X, [H|Rest]) :-
    I > 0,
    I1 is I - 1,
    list_nth0_rest(I1, T, X, Rest).

% Seed-word order: which word anchors the layout. Deterministic = input order;
% seeded = a shuffled copy (the single biggest source of layout variety).
seed_word_order(Words, Order) :-
    ( search_seed(_) -> seeded_permutation(Words, Order) ; Order = Words ).

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

shuffle_bucket(_Count-Entries, Shuffled) :- seeded_permutation(Entries, Shuffled).

% Seed word (grid empty): branch over all words, as the other MRV strategies
% do; the next node builds the full cache (StateOut = none).
select_inc(Words, [], _StateIn, _GridLen, _Start, _Dir, _GIn, Entry, RemWords, none) :-
    seed_word_order(Words, SeedOrder),
    member(Entry, SeedOrder),
    remove_x(Entry, Words, RemWords).
% Otherwise: get current counts (full or incremental), order most-constrained
% first over the connectable words, pick backtrackably. The cache passed
% forward is this node's count map tagged with the chosen word's letters, so
% the next node only recounts words sharing a letter with it. The [_|_] head
% keeps the two clauses index-disjoint on arg 2 ([] vs cons), so no cut is
% needed to stop a seed-node call falling into this clause on backtracking.
select_inc(Words, [P|Ps], StateIn, GridLen, Start, Dir, GIn, Entry, RemWords, StateOut) :-
    PlacedWords = [P|Ps],
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
% atom_chars + separator strip are memoized to one lookup after the first.
% Reset alongside pair_crossings/3 by reset_search_memos/0 at every top-level
% search entry, so inference counts stay self-contained per entry.
:- table answer_letters/2.
answer_letters(Word, Letters) :-
    atom_chars(Word, L0),
    exclude(separator_char, L0, Letters).

% Word-separator characters: enumeration markers, not grid cells. (delete/3,
% the previous strip, is deprecated in SWI 10 and took two passes.)
separator_char(' ').
separator_char('-').

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
    pw_letters(PW, PLetters),
    pw_dir(PW, PDir),
    pw_start(PW, PStart),
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


% The working grid is the bundle gs(LetterGrid, BoundaryGrid) (see init_gs/2):
% two grid-sized compound terms threaded as one. LetterGrid holds placed letters
% (empty = unbound var); BoundaryGrid marks the permanently-unfillable boundary
% cells of placed words (marked = atom `b`, empty = unbound var). Bundling keeps
% the deep counting chain (mrv_count / select_inc / inc_count) a pure pass-through
% - it never inspects the grid, only carries it - so only the sites that actually
% read cells (assign_word, check_word_fits, the greedy scorer, the fragment pin)
% destructure the bundle.
% NB no placed-words argument: placement legality depends only on the grid
% state - the boundary grid subsumed the old placed-words scan, and the
% vestigial _PlacedWords parameter the arity-10 signature still advertised
% was dropped with it (C15).
assign_word(Word, Letters, WLen, Start, Dir, GridLen,
            gs(LGrid, BGrid), Placed, gs(LGrid, BGrid)) :-
    % Reject an off-grid (underflow) start cell. find_intersecting_word/6 can
    % compute a Start < 1 that still passes fits_on_grid (the end lands on-grid),
    % which the old assoc grid rejected implicitly: get_assoc/3 FAILED for the
    % absent negative key. arg/3 instead THROWS for a negative index (0 and
    % >arity still fail, as get_assoc did), so we reject Start < 1 here, once per
    % word. With Start >= 1 and a valid run, every downstream cell index is >= 1,
    % so no read can go negative (an over-grid index fails, never throws).
    Start >= 1,
    % make sure previous cell does not have a letter
    check_prev_cell(Dir, Start, GridLen, LGrid),
    assign_letters(Letters, Start, Dir, GridLen, Cells, LGrid, LGrid),
    % keep every word a maximal run: this word must not fill the boundary cell
    % of an already-placed word (docs/experiments.md I5). O(1)/cell via the
    % boundary grid: no cell of this word may be a marked boundary cell.
    no_word_merge_bg(Cells, BGrid),
    % Precompute the word's END cell (last of the cells run) once here so
    % mark_boundary/5 needs no last/2, and record it in the pw/8 record.
    last(Cells, End),
    % Incrementally maintain the boundary structure: mark THIS word's 0-2 on-grid
    % boundary cells (before-start / after-end). The trail undoes these marks on
    % backtrack, exactly as it undoes the letter bindings.
    mark_boundary(Start, End, Dir, GridLen, BGrid),
    Placed = pw(Word, Letters, Cells, Dir, WLen, Start, End, _Num).


% A valid crossword needs every word to be a MAXIMAL run of cells - bounded on
% both ends, in its own direction, by an empty cell or the grid edge.
% check_prev_cell/4 and check_next_cell/4 enforce that only at the moment a word
% is placed; a later, longer collinear word can pass straight through a shorter
% word's cells (the X==L branch of assign_letters/7 skips adjacency checks) and
% fill the cell just past its end, retroactively making the shorter word a
% "word inside a word". That produces two same-direction answers sharing a start
% cell, which assign_clue_numbers/2 cannot number. We forbid it here: the word
% being placed (cells Cells) must not occupy a marked boundary cell. (The reverse
% - an existing word occupying THIS word's boundary - is already caught by
% check_prev_cell/check_next_cell, so together the invariant holds inductively.)
%
% BoundaryGrid carries, for the current placed set, EXACTLY the union of every
% placed word's on-grid boundary cells - the same set the former placed-words
% scan (member/placed_boundary_cell/memberchk) recomputed per candidate - so the
% accept/reject is identical, at O(1) per cell instead of O(placed words).
no_word_merge_bg(Cells, BGrid) :-
    \+ ( member(Cell, Cells),
         arg(Cell, BGrid, Mark),
         nonvar(Mark) ).

% Mark a placed word's on-grid boundary cells in BoundaryGrid: the cell just
% before its start, and the cell just after its end, in the word's own direction
% (omitting either when the word abuts the grid edge there - the exact edge
% conditions the former placed_boundary_cell/3 used). arg/3 with `b` binds an
% empty boundary cell, or unifies harmlessly when two words' boundaries coincide.
mark_boundary(Start, End, Dir, GridLen, BGrid) :-
    (   is_start_cell(Dir, Start, GridLen)
    ->  true
    ;   prev_cell(Dir, Start, GridLen, Before),
        arg(Before, BGrid, b)
    ),
    (   is_end_cell(Dir, End, GridLen)
    ->  true
    ;   next_cell(Dir, End, GridLen, After),
        arg(After, BGrid, b)
    ).


% check_word_fits/5: a pure legality PROBE for the mrv_count Cap=2 counting path
% (mrv_count_goal) and the greedy constructor's candidate scan (arrange.pl
% word_best_placement/7). It answers "would assign_word/9 accept this candidate?" with
% NO side effects - it binds no letter cell, marks no boundary cell, and builds
% neither the Cells run nor the pw/8 record (all of which assign_word materialized
% only for count_upto2 to discard) - so it is cheaper per counted candidate and
% leaves nothing for count_upto2's \+ to unwind (no trail churn in the loop).
%
% CORRECTNESS (must match assign_word exactly): Start >= 1 ; check_prev_cell ;
% per cell (not a marked boundary cell, then a matching existing letter OR an
% empty cell with free perpendicular neighbours) ; check_next_cell. Not binding is
% sound for a one-shot check: a straight run's cells are all distinct (no
% self-collision within one word) and adj_is_free reads only the two PERPENDICULAR
% off-line neighbours - never a cell this word occupies - so nothing a check would
% bind is ever re-read by the same check. Each counted candidate is thus evaluated
% against the same grid, exactly as assign_word is (whose bindings the counting
% loop backtracks away between candidates anyway). The boundary test is folded
% into the per-cell walk (an O(1) BoundaryGrid read) rather than run as a separate
% pass; the accept/reject SET is unchanged (a pure conjunction of the same
% conditions), so every count stays byte-identical to the placed-words scan.
check_word_fits(Letters, Start, Dir, GridLen, gs(LGrid, BGrid)) :-
    Start >= 1,
    check_prev_cell(Dir, Start, GridLen, LGrid),
    check_letters(Letters, Start, Dir, GridLen, BGrid, LGrid).

% Non-binding twin of assign_letters/7 with the O(1) boundary test folded in:
% same per-cell letter test (matching existing letter, or empty cell with free
% perpendicular neighbours) and same trailing check_next_cell/4, but it binds no
% cell, marks no boundary, and builds no Cells run or pw/8 record.
check_letters([], Num, Dir, GridLen, _BGrid, LGrid) :-
    check_next_cell(Dir, Num, GridLen, LGrid).
check_letters([L|Ls], Num, Dir, GridLen, BGrid, LGrid) :-
    arg(Num, BGrid, Mark), var(Mark),      % not a placed word's boundary cell
    arg(Num, LGrid, Cell),
    (   nonvar(Cell)
    ->  Cell == L
    ;   adj_is_free(Dir, Num, GridLen, LGrid)
    ),
    next_cell(Dir, Num, GridLen, Num2),
    check_letters(Ls, Num2, Dir, GridLen, BGrid, LGrid).


% Previous cell before start of word. Make sure it doesn't contain
% anything.
check_prev_cell(Dir, Num, GridLen, Grid) :-
    (   % don't check if start letter is start of a row/col
        is_start_cell(Dir, Num, GridLen)
    ->  true
    ;   % otherwise prev cell must be empty (an unbound var)
        prev_cell(Dir, Num, GridLen, Prev),
        arg(Prev, Grid, Cell), var(Cell)
    ).


% Next cell after end of word. Make sure it doesn't contain anything.
check_next_cell(Dir, Num, GridLen, Grid) :-
    prev_cell(Dir, Num, GridLen, Prev),
    (   % no need to check if prev was end of row/col
        is_end_cell(Dir, Prev, GridLen)
    ->  true
    ;   % then this cell must be empty (an unbound var)
        arg(Num, Grid, Cell), var(Cell)
    ).


% Assign each letter of the word, checking that adjacent cells are empty to
% prevent words being placed next to each other. The grid is mutated in place by
% binding empty cell vars, so GIn and GOut are the SAME term (backtracking
% unbinds via the trail); they thread as one.

% Last letter of word, make sure next cell is free
assign_letters([], Num, Dir, GridLen, [], Grid, Grid) :-
    check_next_cell(Dir, Num, GridLen, Grid).


assign_letters([L|Ls], Num, Dir, GridLen, [Num|RestCells], Grid, GridOut) :-
    arg(Num, Grid, Cell),
    (   nonvar(Cell)
    ->  % existing letter in this cell matches letter being placed,
        % nothing needs doing, we can continue to next letter
        Cell == L
    ;   % empty cell (an unbound var), so check adjacent cells are free
        % and then bind this cell to the letter
        adj_is_free(Dir, Num, GridLen, Grid),
        Cell = L
    ),
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
% cell begins at most one across AND one down word, never 3+. The group's first
% word is consumed here; the tail goes through add_group_tail/4, which has
% deliberately only [] and [W2] clauses - a >2 group is structurally impossible
% and (by design) has no clause (P17).
%
% Head shapes are []/[_|_] on both predicates, so SWI's two-clause special-case
% indexing selects deterministically - assign_clue_numbers/2 succeeds WITHOUT a
% choicepoint, and callers need no defensive once/1 (determinism-audit L2).
add_clue_nums([], _, []).
add_clue_nums([_-[W|More]|Rest], ClueNum, [WClue|Clues]) :-
    add_clue_word(W, ClueNum, WClue),
    add_group_tail(More, ClueNum, Clues, RestClues),
    ClueNum2 is ClueNum + 1,
    add_clue_nums(Rest, ClueNum2, RestClues).

% The group's second word shares the first word's clue number (a cell that
% starts both an across and a down word). Difference-list tail: RestClues is
% where the remaining groups' clues continue.
add_group_tail([], _, Clues, Clues).
add_group_tail([W2], ClueNum, [WClue2|Clues], Clues) :-
    add_clue_word(W2, ClueNum, WClue2).


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
    pw_cells(PW, Cells),
    pw_letters(PW, Letters),
    pw_dir(PW, Dir),
    pw_num(PW, Num),
    pw_start(PW, Start),
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

word_canon_key(PW, Num-Dir) :- pw_num(PW, Num), pw_dir(PW, Dir).


placed_to_word(MetaAssoc, GridLen, PW, WordObj) :-
    pw_answer(PW, Answer),
    pw_dir(PW, Dir),
    pw_num(PW, Num),
    pw_cells(PW, Cells),
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

% A placed word is a fixed-arity compound record (formerly a word{...} dict):
%   pw(Answer, Letters, Cells, Dir, Len, Start, End, Num)
% End is the word's last cell (precomputed in assign_word). `Num` (the clue
% number) is a fresh, unbound var until assign_clue_numbers/2 rebuilds the record
% with it bound (a pure, setarg-free rebuild). The accessors below are arg/3-style
% one-liners the compiler indexes cheaply; call sites read like the old dict gets.
pw_answer( pw(A,_,_,_,_,_,_,_), A).
pw_letters(pw(_,L,_,_,_,_,_,_), L).
pw_cells(  pw(_,_,C,_,_,_,_,_), C).
pw_dir(    pw(_,_,_,D,_,_,_,_), D).
pw_len(    pw(_,_,_,_,L,_,_,_), L).
pw_start(  pw(_,_,_,_,_,S,_,_), S).
pw_end(    pw(_,_,_,_,_,_,E,_), E).
pw_num(    pw(_,_,_,_,_,_,_,N), N).


start_is(PW, Start) :- pw_start(PW, Start).

% Bind the assigned clue number into a placed word by rebuilding the record with
% Num filled (pure; no setarg, no dict mutation). The unnumbered record's Num var
% is left untouched, so this is safe to backtrack over.
add_clue_word(pw(A,L,C,D,Len,S,E,_), ClueNum, pw(A,L,C,D,Len,S,E,ClueNum)).



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

% The working-grid bundle threaded through the search: gs(LetterGrid,
% BoundaryGrid). LetterGrid is the letter grid above; BoundaryGrid is a second
% grid-sized all-var term in which assign_word marks the placed words' boundary
% cells (atom `b`), giving no_word_merge_bg/2 and check_word_fits/5 an O(1)
% arg/3+nonvar boundary test instead of a scan over all placed words. Boundary
% cells stay UNBOUND in LetterGrid (they must still read as empty to
% check_prev/next_cell and adj_is_free - a nonvar sentinel there would change
% the search); the mark lives only in BoundaryGrid. Both structures backtrack
% via the trail, so the bundle threads exactly as the single grid did.
init_gs(GridLen, gs(LGrid, BGrid)) :-
    init_grid(GridLen, LGrid),
    init_grid(GridLen, BGrid).



% (the solved crossword is emitted as JSON by emit_json/3; see the
% "JSON output" section above)



% generic utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove_x(X,L,R) :- R is L with the first ==-identical occurrence of X
% removed from it (identity, not unification: entries are selected from L by
% member/2 upstream, so the match must be the very term, never a unifiable
% lookalike). The []/[_|_] heads are index-disjoint and the if-then-else
% leaves no choicepoint, so the traversal is deterministic - the search
% backtracks straight to member/2's real alternative instead of grinding
% through one dead remove_x choicepoint per traversed element.
remove_x(_, [], []).
remove_x(Y, [X|Xs], R) :-
	(   Y == X
	->  R = Xs
	;   R = [X|Tail],
	    remove_x(Y, Xs, Tail)
	).
