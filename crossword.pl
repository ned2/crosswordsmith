#!/usr/bin/swipl

:- set_prolog_flag(verbose, silent).
% initialization(main, main) (vs plain `initialization main`) means main is
% only run when this file is executed as a script. When the file is loaded
% (consulted) by the test suite, main does not fire, so the predicates below
% can be unit tested in isolation. It also gives a proper process exit code.
:- initialization(main, main).

% crossword.pl - A crossword layout generator in Prolog
% Copyright (C) 2011  Ned Letcher - nedned.net


% See the file README.mk for background and a general overview..
% 
% This program can be run as an SWI PrologScript like this:
%
% $ ./crossword.pl --input <input_file> <grid_length> <start_loc>
%
% Or you can find a vaguely random solution by shuffling the input words                                
% and the order of the start locations:
%
% $ ./crossword.pl --input <input_file> --shuffle <grid_length>
%
% Where grid_length is integer specifying the dimensions of the crossword and
% start_loc specifies where the first word of the crossword is placed and can be
% one of {topleft_down, topleft_across, topright, bottomleft}.


% The program uses two simple data structures. The first is a list of
% words that have been placed on the crossword grid with each element
% being a dict of attributes of the word. See the placed words utility
% predicate section below for the structure of this dict.
%
% The second data structure is as association list used to to store
% the contents of each cell in the grid. The keys are the numbers of
% the cells (1 through to GridLen*GridLen) and the values are the
% contents of the cell.

% Used by the JSON output emitter (canonically library(json) on SWI 10+).
:- use_module(library(http/json)).

% Command-line option parsing (--input/--shuffle/--all/--out + positionals).
:- use_module(library(optparse)).

% limit/2, used by the capped placement count in the mrv_capped strategy.
:- use_module(library(solution_sequences)).

% The greedy quality layout engine (docs/cryptic-layout-spec.md), reached via
% --quality. Loaded from the same directory as this script so it resolves
% regardless of the working directory. ensure_loaded avoids a double-load when
% a harness consults both files.
:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, 'quality.pl', QualityFile),
   ensure_loaded(QualityFile).



% program predicates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Command-line entry point. Options are parsed with library(optparse), which
% strips the flags from argv and returns the leftover positional arguments
% (grid_length and an optional start_loc, as atoms). So the flags compose
% freely and in any order, and adding another is a one-line spec entry rather
% than another arg-matching clause. The dispatch lives in run/2.
main :-
    current_prolog_flag(argv, Argv),
    opts_spec(Spec),
    opt_parse(Spec, Argv, Opts, Positional),
    run(Opts, Positional).

% The option set. `input`/`out` default to '' meaning "not given". The
% `strategy` default is the production default (default_strategy/1), so the
% CLI and the find_crossword/5 + crossword/3 convenience wrappers all agree.
opts_spec(Spec) :-
    default_strategy(DefStrategy),
    Spec =
    [ [opt(help),    type(boolean), default(false),
       shortflags([h]), longflags([help]),
       help('show this help and exit')],
      [opt(input),   type(atom),    default(''), meta('FILE'),
       longflags([input]),
       help('load clues from FILE (.json or .pl)')],
      [opt(shuffle), type(boolean), default(false),
       longflags([shuffle]),
       help('shuffle the words and start positions for a random layout')],
      [opt(all),     type(boolean), default(false),
       longflags([all]),
       help('count every solution instead of emitting one (see README)')],
      [opt(out),     type(atom),    default(''), meta('FILE'),
       longflags([out]),
       help('write output to FILE instead of stdout')],
      [opt(strategy), type(atom),   default(DefStrategy), meta('STRAT'),
       longflags([strategy]),
       help('variable-ordering strategy: baseline, mrv, mrv_capped, or mrv_inc')],
      [opt(quality), type(boolean), default(false),
       longflags([quality]),
       help('cryptic-style quality layout: greedy density construction, the engine picks the grid, words may be dropped (no grid/start args needed)')],
      [opt(min_half), type(boolean), default(false),
       longflags(['min-half']),
       help('quality floor: every word at least half-checked (drops words that cannot be)')],
      [opt(max_unch), type(integer), default(-1), meta('K'),
       longflags(['max-unch']),
       help('quality floor: no word has more than K unchecked cells in a row')],
      [opt(all_words), type(boolean), default(false),
       longflags(['all-words']),
       help('quality floor: every input word must be placed (fail rather than drop)')]
    ].

% Dispatch on the parsed options. The positional grammar is:
%   --input FILE <grid_length> <start_loc>   solve and emit one JSON solution
%   --input FILE --shuffle <grid_length>     ditto, words/start shuffled
%   --input FILE --all <grid_length> [loc]   count solutions (all start_locs if loc omitted)
% --out FILE is orthogonal and composes with any of these.
run(Opts, _) :-
    memberchk(help(true), Opts),
    !,
    print_usage.
run(Opts, Positional) :-
    memberchk(input(InputFile), Opts),
    require_input_file(InputFile),
    memberchk(out(OutFile), Opts),
    memberchk(strategy(Strategy), Opts),
    require_strategy(Strategy),
    load_clues(InputFile, Words),
    (   memberchk(quality(true), Opts)
    ->  build_floors(Opts, Floors),                       % no grid/start needed
        with_output(OutFile, quality_solve(Words, Floors))
    ;   memberchk(all(true), Opts)
    ->  positional_all(Positional, GridLen, StartLoc),   % StartLoc unbound if absent
        with_output(OutFile, count_solutions(Strategy, GridLen, Words, StartLoc))
    ;   memberchk(shuffle(true), Opts)
    ->  Positional = [GridLenArg],
        atom_number(GridLenArg, GridLen),
        solve_shuffled(Strategy, GridLen, Words, OutFile)
    ;   Positional = [GridLenArg, StartLocArg],
        atom_number(GridLenArg, GridLen),
        valid_loc(StartLocArg),
        with_output(OutFile, crossword(Strategy, GridLen, Words, StartLocArg))
    ).

print_usage :-
    opts_spec(Spec),
    opt_help(Spec, Help),
    format("Usage: crossword.pl --input <file> [options] <grid_length> [<start_loc>]~n~n~w", [Help]).

require_input_file('') :-
    !,
    throw(error(missing_input_file, _)).
require_input_file(_).

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

% Build the quality-engine Floors dict from the CLI flags (groundedness ->
% mode: no floors = auto, some = manual/partial, all = strict). See quality.pl.
build_floors(Opts, floors{min_half:MH, max_unch_run:MU, all_words:AW}) :-
    ( memberchk(min_half(true), Opts)              -> MH = on ; MH = off ),
    ( memberchk(max_unch(K), Opts), integer(K), K >= 0 -> MU = K ; MU = off ),
    ( memberchk(all_words(true), Opts)             -> AW = on ; AW = off ).

% Positional args for the --all path: grid_length is required; start_loc is
% optional and left unbound when absent, so all four start positions are
% enumerated by the solver.
positional_all([GridLenArg], GridLen, _StartLoc) :-
    atom_number(GridLenArg, GridLen).
positional_all([GridLenArg, StartLoc], GridLen, StartLoc) :-
    atom_number(GridLenArg, GridLen),
    valid_loc(StartLoc).

% Solve with shuffled words and start order; member/2 is the backtrack point
% that tries successive shuffled start positions until one yields a layout.
solve_shuffled(Strategy, GridLen, Words, OutFile) :-
    start_locs(Locs),
    shuffle(Words, UseWords),
    shuffle(Locs, ShuffledLocs),
    member(StartLoc, ShuffledLocs),
    with_output(OutFile, crossword(Strategy, GridLen, UseWords, StartLoc)).

% Count solutions for one start position (or all, with StartLoc unbound).
count_solutions(Strategy, GridLen, Words, StartLoc) :-
    all_crossword(Strategy, GridLen, Words, StartLoc, Num),
    writeln(Num).

valid_loc(Loc) :-
    start_locs(Locs),
    memberchk(Loc, Locs).

% Run Goal with its output sent to a file, or straight to stdout when '' (no
% --out). For a file the output is captured first and written only if Goal
% succeeds, so a no-solution run leaves no empty file behind; the stdout path
% is the unchanged direct write.
with_output('', Goal) :-
    !,
    call(Goal).
with_output(File, Goal) :-
    with_output_to(string(Text), Goal),
    setup_call_cleanup(open(File, write, S), write(S, Text), close(S)).


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
% normalised to atoms so the emit-time answer_meta/3 join and
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

:- multifile prolog:error_message//1.
prolog:error_message(json_no_clues_array) -->
    [ 'clues file: expected a JSON object with a "clues" array' ].
prolog:error_message(json_invalid_answer(Entry)) -->
    [ 'clues file: every entry needs a string "answer" (offending entry: ~q)'-[Entry] ].
prolog:error_message(json_invalid_meta(Answer)) -->
    [ 'clues file: "meta" for answer ~q must be a JSON object'-[Answer] ].
prolog:error_message(missing_input_file) -->
    [ 'missing required --input FILE option' ].
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


% Top level predicate for finding the number of solutions for the
% crossword for a specific starting position.
all_crossword(Strategy, GridLen, Words, StartLoc, Num) :-
    length(Sols, Num),
    findall(Grid, find_crossword(Strategy, GridLen, Words, StartLoc, Grid, _), Sols).


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
    atom_chars(Word, Letters),
    delete(Letters, ' ', Letters2),
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

% Count Entry's viable placements right now, bounded by Cap. findall enumerates
% each candidate (Start,Dir) from find_intersecting_word that also survives the
% assign_word adjacency/bounds checks; findall undoes the goal's bindings, so
% the caller's Start/Dir (ground on the first word, unbound after) are left
% untouched.
mrv_count(Cap, PlacedWords, GridLen, Start, Dir, GIn, Entry, Count) :-
    Entry = [Word|_],
    atom_chars(Word, Letters),
    delete(Letters, ' ', Letters2),
    length(Letters2, WLen),
    findall(t,
            capped(Cap,
                   ( find_intersecting_word(Letters2, WLen, PlacedWords, GridLen,
                                            Start, Dir),
                     assign_word(Word, Letters2, WLen, Start, Dir, GridLen,
                                 PlacedWords, GIn, _Placed, _G1) )),
            Ts),
    length(Ts, Count).

capped(unbounded, Goal) :- call(Goal).
capped(N, Goal) :- integer(N), limit(N, Goal).


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
    atom_chars(Word, Letters),
    delete(Letters, ' ', Letters2),
    length(Letters2, WLen),
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Start, Dir, GridLen, PlacedWords, GIn, Placed, G1),
    assign_words_inc(RemWords, [Placed|PlacedWords], StateOut, GridLen,
                     _Start, _Dir, G1, GOut, Out).

% Seed word (grid empty): branch over all words, as the other MRV strategies
% do; the next node builds the full cache (StateOut = none).
select_inc(Words, [], _StateIn, _GridLen, _Start, _Dir, _GIn, Entry, RemWords, none) :-
    !,
    member(Entry, Words),
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
    pairs_values(Sorted, Ordered),
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
    atom_chars(Word, L0),
    delete(L0, ' ', Letters).

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
    intersection(Letters, PLetters, Vals),
    list_to_set(Vals, Vals2),
    member(Val, Vals2),
    position(Val, PLetters, PPos),
    position(Val, Letters, Pos),
    calc_num(PDir, GridLen, PPos, PStart, PNum),
    swap_dir(PDir, Dir),
    calc_start(Dir, GridLen, Pos, PNum, Start),
    fits_on_grid(Dir, Start, WLen, GridLen).


assign_word(Word, Letters, WLen, Start, Dir, GridLen, PlacedWords, GIn, Placed, GOut) :-
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
check_prev_cell(Dir, Num, GridLen, G) :-
    (
     % don't check if start letter is start of a row/col
     is_start_cell(Dir, Num, GridLen)
    ;
     % otherwise prev cell must be empty
     prev_cell(Dir, Num, GridLen, Prev),
     get_assoc(Prev, G, empty)
    ), !.


% Next cell after end of word. Make sure it doesn't contain anything.
check_next_cell(Dir, Num, GridLen, G) :-
    prev_cell(Dir, Num, GridLen, Prev),
    (
     % no need to check if prev was end of row/col
     is_end_cell(Dir, Prev, GridLen)
    ;
     % then this cell must be empty
     get_assoc(Num, G, empty)
    ), !.


% Assign each letter of the word, checking that adjacent cells
% are empty to prevent words being placed next to each other.

% Last letter of word, make sure next cell is free
assign_letters([], Num, Dir, GridLen, [], G, G) :- 
    check_next_cell(Dir, Num, GridLen, G).


assign_letters([L|Ls], Num, Dir, GridLen, [Num|RestCells], GIn, GOut) :-
    get_assoc(Num, GIn, X),
    (
     % existing letter in this cell matches letter being placed,
     % nothing needs doing, we can continue to next letter
     X == L,
     G1 = GIn
    ;
     % no letter in this cell, so check adjacent cells are free
     % and then add letter to this cell
     X == empty,
     adj_is_free(Dir, Num, GridLen, GIn),
     put_assoc(Num, GIn, L, G1)
    ), !,
    next_cell(Dir, Num, GridLen, Num2),
    assign_letters(Ls, Num2, Dir, GridLen, RestCells, G1, GOut).


% check that adjacent cells are empty
adj_is_free(down, Num, GridLen, G) :-
    N1 is Num - 1,
    N2 is Num + 1,
    M is Num mod GridLen,
    (
     M == 0 -> % last cell in row
     get_assoc(N1, G, empty)
    ;
     M == 1 -> % first cell in row
     get_assoc(N2, G, empty)
    ;
     get_assoc(N1, G, empty),
     get_assoc(N2, G, empty)
    ), !.

adj_is_free(across, Num, GridLen, G) :-
    N1 is Num - GridLen,
    N2 is Num + GridLen,
    LastCell is (GridLen * GridLen),
    (
     N1 =< 0 -> % before beginning of grid
     get_assoc(N2, G, empty)
    ;
     N2 > LastCell -> % after end of grid
     get_assoc(N1, G, empty)
    ;
     get_assoc(N1, G, empty),
     get_assoc(N2, G, empty)
    ), !.


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

:- multifile prolog:error_message//1.
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


% Build the `words` array. Metadata is rejoined from the input list by answer.
build_words(PlacedWords, Words, GridLen, WordObjs) :-
    maplist(placed_to_word(Words, GridLen), PlacedWords, WordObjs).


placed_to_word(Words, GridLen, PW, WordObj) :-
    get_dict(answer, PW, Answer),
    get_dict(dir, PW, Dir),
    get_dict(num, PW, Num),
    get_dict(cells, PW, Cells),
    maplist(cell_coord(GridLen), Cells, Coords),
    answer_meta(Answer, Words, Meta),
    WordObj = _{number: Num, direction: Dir, answer: Answer,
                cells: Coords, meta: Meta}.


% Look up the (opaque) metadata for an answer in the input list. An entry may
% omit metadata ([Answer]); answers are unique (check_unique_answers/1).
answer_meta(Answer, Words, Meta) :-
    member(Entry, Words),
    Entry = [A|_],
    A == Answer,
    !,
    ( Entry = [_, M] -> Meta = M ; Meta = _{} ).


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

direction(down).
direction(across).


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
is_end_cell(down, Num, Length) :- Num >= (Length - 1) * Length.


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


new_tile(Num, Num-empty).

init_grid(GridLen, Grid) :-
    NumTiles is GridLen * GridLen,
    numlist(1, NumTiles, Tiles),
    maplist(new_tile, Tiles, TupleList),
    list_to_assoc(TupleList, Grid).



% (the solved crossword is emitted as JSON by emit_json/3; see the
% "JSON output" section above)



% generic utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

head([X|_], X).

length_sort(Words, OrderedWords) :-
    map_list_to_pairs(atom_length, Words, Pairs),
    keysort(Pairs, SortedPairs),
    pairs_values(SortedPairs, SortedWords),
    reverse(SortedWords, OrderedWords).


% finds all positions X occurs in List (over backtracking)
position(X, List, Pos) :- x_position(List, X, 1, Pos).

x_position([], _, _, _) :- false.
x_position([X|_], X, Pos, Pos).
x_position([_|Ys], X, N, Pos) :-
    N2 is N + 1,
    x_position(Ys, X, N2, Pos).


% remove_x(X,L,R) :- R is L with first occurrence of X removed from it.
remove_x(Y,[X|Xs],[X|Tail]) :-
	Y \== X,
	remove_x(Y,Xs,Tail).
remove_x(X,[X|Xs],Xs) :- !.
remove_x(_,[],[]).


%% shuffle(ListIn, ListOut) - randomly shuffles
%% ListIn and unifies it with ListOut
shuffle([], []) :- !.
shuffle(List, [Element|Rest]) :-
    choose(List, Element),
    delete(List, Element, NewList),
    shuffle(NewList, Rest).


%% choose(List, Elt) - chooses a random element
%% in List and unifies it with Elt.
choose(List, Elt) :-
    length(List, Length),
    Index is random(Length),
    nth0(Index, List, Elt).
