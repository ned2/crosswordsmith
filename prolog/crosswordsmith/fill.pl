% fill.pl - Flavour-B grid-first, open-dictionary auto-fill (design-spec §8.4;
% OD-1..4 resolved by DP-1/DP-2). Take a pre-validated legal blocked grid (a §8.3
% stock-grid mask), derive its slots, and fill every slot with a dictionary word
% subject to crossing constraints, with the user's words pinned as seeds (the
% §6.6 fragment primitive). Output is the canonical layout, so it composes with
% lint/export.
%
% The grid model: each white cell is one logical VARIABLE, shared between its
% across and down slot. Assigning a word to a slot unifies the slot's variable
% list with the word's letters; crossings are then consistent by construction
% (the shared variable), and dead-ends backtrack. MRV ordering (fewest matching
% candidates first) + a node/inference budget make it deterministic and bounded.
%
% Exports: fill_solve/4 (the CLI seam), the prebuilt-index seams
% (fill_solve_index/5, fill_save_index/2,3), and the five benchmark
% measurement seams (load_dict/3, fill_grid/4, fill_attempt/8, apply_seeds/4,
% seeded_slot/2) the fill product bench times directly — see the export list
% below. Every other predicate is internal (tests reach them as
% crosswordsmith_fill:Pred(...)). All project dependencies are explicit
% imports below (stockgrid, metrics, arrange, core).

:- module(crosswordsmith_fill,
          [ fill_solve/4,
            fill_solve/5,
            fill_solve_index/5,
            fill_solve_index/6,
            fill_save_index/2,
            fill_save_index/3,
            % the benchmark measurement seams (benchmarks/fill_subjects.pl,
            % run_fill.pl): the dict_load/grid/search buckets time load_dict,
            % fill_grid and fill_attempt directly, and the search sampler
            % rebuilds fresh seeded slots via apply_seeds + seeded_slot.
            % Every recorded fill baseline count (fill_baseline.json,
            % fill_identity.sha256) is defined against THESE predicates —
            % benchmarks must never re-point the seams elsewhere.
            load_dict/3,
            fill_grid/4,
            fill_attempt/8,
            apply_seeds/4,
            seeded_slot/2
          ]).

% All library imports carry explicit import lists so a
% qsave_program(..., [autoload(false)]) build resolves them (P11/C5). (No JSON
% import: this module's emit delegates to core's emit_json/3 / arrange's
% emit_arrange/4 — the old library(http/json) line here was vestigial.)
:- use_module(library(apply),
              [convlist/3, exclude/3, foldl/4, include/3, maplist/2, maplist/3]).
:- use_module(library(lists),
              [ append/2, append/3, member/2, min_list/2, nextto/3, nth0/3,
                nth1/3, numlist/3, select/3, selectchk/3, subtract/3,
                sum_list/2 ]).
:- use_module(library(ordsets),
              [list_to_ord_set/2, ord_intersect/2, ord_intersection/3]).
:- use_module(library(assoc),
              [ assoc_to_list/2, assoc_to_values/2, gen_assoc/3,
                get_assoc/3, is_assoc/1, list_to_assoc/2,
                ord_list_to_assoc/2, put_assoc/4 ]).
:- use_module(library(pairs),
              [group_pairs_by_key/2, map_list_to_pairs/3, pairs_keys_values/3,
               pairs_values/2]).
% F-L2: dict-file SHA-256 (artifact integrity). NB library(sha)'s file defines
% module `crypto_hash` on 10.1.10; both predicates are on its export list.
:- use_module(library(sha), [sha_hash/3, hash_atom/2]).
% F-L2: fast_write/fast_read artifact I/O.
:- use_module(library(fastrw), [fast_read/2, fast_write/2]).
% option/2,3: the artifact Meta is a keyed list of unary Key(Value) terms -
% literally an SWI option list - and fill_save_index's Options is one too.
:- use_module(library(option), [option/2, option/3]).
:- use_module(library(error), [must_be/2]).
% read_file_to_string/3: whole-file reads (dict lines + the SHA-256
% fingerprint) without the hand-rolled open/read_string/close dance (C39).
:- use_module(library(readutil), [read_file_to_string/3]).

% Slot derivation from a stock-grid mask.
:- use_module(crosswordsmith(stockgrid),
              [stockgrid_load/2, mask_white_cells/3, grid_run/4]).

% word_letters/3: the separator-stripped placement footprint of a seed answer.
:- use_module(crosswordsmith(metrics), [word_letters/3]).

% load_fragment/4: seeds arrive in either §6.6 fragment form - canonical, or
% the thin [{answer,row,col,dir}] list desugared on the fill grid's own side
% (the SizeCtx apply_seeds passes), so thin seeds land on exactly the grid
% they pin. emit_arrange/4: emit_fill's `max` mode delegates the cropped
% emit to arrange.
:- use_module(crosswordsmith(arrange), [load_fragment/4, emit_arrange/4]).

% Numbering + the canonical JSON emit for the filled layout, plus the placed-word
% record (pw/8) accessor fill uses to recover each answer for the emit metadata.
% check_unique_answers/1: the unique-answers guard emit_json/3's metadata join
% requires; nothing upstream of fill's emit runs it, so emit_fill/4 runs it
% itself.
:- use_module(crosswordsmith(core),
              [ assign_clue_numbers/2, check_unique_answers/1, emit_json/3,
                verbose_report/2, pw_answer/2,
                % §8.4b (DP-6) seed perturbation: fill reads the SAME
                % module-owned PRNG as arrange's §7.6 knobs (wasm parity).
                % prng_draw/1: the §8.4c top-3 pick (DP-8) draws the same
                % stream directly; prng_seed_state/1 reseeds it from the
                % pinned constant for default-path restarts (C3).
                current_search_seed/1, seeded_permutation/2, prng_draw/1,
                prng_seed_state/1 ]).

fill_budget(800_000_000).   % inference budget (determinism via INV-2, bounded)


% --- build the slots from a stock-grid mask ----------------------------------

%!  fill_grid(+GridFile:atom, -Size:integer, -Slots:list, -CellVar:assoc) is det.
%
%   Derive the slot set of the stock-grid mask GridFile (§8.3). Each slot is
%   slot(Start, Dir, Cells, Vars): Vars are the shared per-cell logical
%   variables (so an across and a down slot crossing at a cell share that
%   cell's variable). Returns the slots + the cell->var assoc. Throws on a
%   malformed grid file (stockgrid_load/2). Benchmark seam: the fill bench's
%   `grid` bucket (grid_inf) times exactly this goal, and its search sampler
%   rebuilds fresh slots through it outside the timed window.
fill_grid(GridFile, Size, Slots, CellVar) :-
    stockgrid_load(GridFile, grid(_Name, Size, Mask)),
    mask_white_cells(Mask, Size, WhiteSet),
    init_cell_vars(WhiteSet, CellVar),
    % findall COPIES its template, which would rename the shared cell variables
    % apart and break crossings - so collect ground slot specs first, then wire
    % in the shared variables (from the one CellVar assoc) outside findall.
    findall(spec(Start, Dir, Cells),
            ( member(Dir, [across, down]),
              grid_run(Size, WhiteSet, Dir, Cells),
              Cells = [Start|_] ),
            Specs),
    maplist(spec_slot(CellVar), Specs, Slots).

spec_slot(CellVar, spec(Start, Dir, Cells), slot(Start, Dir, Cells, Vars)) :-
    slot_vars(Cells, CellVar, Vars).

% WhiteCells is a strictly-ascending ordset (mask_white_cells/3 ends in
% list_to_ord_set/2), so ord_list_to_assoc/2 constructs the balanced assoc
% directly - no per-element put_assoc rebalance (C31).
init_cell_vars(WhiteCells, Assoc) :-
    pairs_keys_values(Pairs, WhiteCells, _FreshVars),
    ord_list_to_assoc(Pairs, Assoc).

% NB: a yall lambda `[C,V]>>get_assoc(C, CellVar, V)` would copy_term its free
% variable CellVar on every call - cloning the assoc and the unbound cell
% variables inside it, so nothing would be shared. A named helper passes the
% assoc as a plain argument (no copy), keeping the crossing cells' variables
% genuinely shared between the across and down slots.
slot_vars(Cells, CellVar, Vars) :-
    maplist(cell_var(CellVar), Cells, Vars).
cell_var(CellVar, Cell, Var) :- get_assoc(Cell, CellVar, Var).


% --- seeds: pin fragment words at their slots (hard pins, OD-3) ---------------

%!  apply_seeds(+SeedFile:atom, +Size:integer, +Slots:list,
%!              -SeededKeys:list) is det.
%
%   Pin every seed of the §6.6 fragment file SeedFile into Slots. A seed
%   answer is unified into the slot whose cells/direction match; a seed that
%   matches no slot, clashes with another seed at a crossing, or duplicates
%   another seed's answer THROWS (hooked errors, reported before searching).
%   SeededKeys collects each pinned slot's Start-Dir so the search and the
%   unfillable report can EXEMPT them: a seed is a hard pin (a setter's
%   theme/own word) and need NOT be a dictionary word (R2/OD-3). Its letters
%   still constrain crossing slots via the shared cell variables, and a slot
%   completed by crossings (not a seed) is still validated against the
%   dictionary. Size is the fill grid's own side, passed as load_fragment's
%   SizeCtx: it frames a thin seed file's desugar so its cell numbers land on
%   THIS grid (a canonical file ignores it - matching stays purely by
%   cells/direction). Benchmark seam: the search sampler's seeded-rung slot
%   rebuild (benchmarks/fill_subjects.pl build_search_slots/3).
apply_seeds(SeedFile, Size, Slots, SeededKeys) :-
    load_fragment(SeedFile, Size, _FragGridLen, Frags),
    check_unique_seed_answers(Frags),
    foldl(apply_seed(Slots), Frags, [], SeededKeys).

% Two seeds pinning the SAME answer would collide at emit (crosswords do not
% repeat answers; the emit metadata join is keyed on the answer), and neither
% slot is searched, so the search-side dedup (seed_used/3) cannot catch that
% shape. Reject up front with a clean hooked error, consistent with
% fill_seed_no_slot / fill_seed_clash. Detection mirrors core's
% check_unique_answers/1: msort keeps duplicates, equal neighbours = a dup.
check_unique_seed_answers(Frags) :-
    findall(A, member(frag(A, _, _, _), Frags), Answers),
    msort(Answers, Sorted),
    (   nextto(Dup, Dup, Sorted)
    ->  throw(error(fill_seed_duplicate(Dup), _))
    ;   true
    ).

apply_seed(Slots, frag(Answer, Dir, _Start, CellNums), SeededIn, SeededOut) :-
    ( member(slot(SStart, Dir, CellNums, Vars), Slots)
    ->  word_letters([Answer], Letters, _),
        ( Vars = Letters
        ->  SeededOut = [SStart-Dir|SeededIn]
        ;   throw(error(fill_seed_clash(Answer), _)) )
    ;   throw(error(fill_seed_no_slot(Answer), _)) ).

%!  seeded_slot(+SeededKeys:list, +Slot) is semidet.
%
%   A slot is a seed pin (exempt from search + empty_slots) iff its Start-Dir
%   was recorded by apply_seed (Start+Dir uniquely identify a slot).
%   Benchmark seam: the search sampler's exclude/3 filter that derives
%   SearchSlots from AllSlots, mirroring fill_prepare's own use.
seeded_slot(SeededKeys, slot(Start, Dir, _, _)) :- memberchk(Start-Dir, SeededKeys).


% --- dictionary: in-memory pattern index (OD-2) ------------------------------

%!  load_dict(+File:atom, -DictByLen:assoc, -Index:assoc) is det.
%
%   Read the text dictionary File into the in-memory pattern index.
%   DictByLen: assoc Len -> list of words (each a list of upper-case char
%   atoms). Index: assoc k(Len,Pos,Char) -> ordset of indices into
%   DictByLen[Len]. A slot's candidates = the words of its length matching
%   its already-bound positions (intersection of the index sets), or all
%   words of that length if nothing bound. Throws on an unreadable File.
%   Delegate of load_dict/5 with default options: scored `word;score` lines
%   are ingested per §8.4a (default prune `score >= 1`, DP-5), and a plain
%   wordlist (no `;` anywhere in the file) takes the pre-scoring path
%   verbatim, so this seam's behavior and counts are unchanged for every
%   existing plain-dict rung.
%
%   The File is DEFINED as UTF-8 (decoded with an explicit encoding(utf8)
%   pin - never the process locale). Each line is normalized fold-then-strict
%   into A-Z letters: ASCII case-mapped, Latin-1/Extended-A diacritics folded
%   through the explicit fold_char/2 table (NFC and NFD forms normalize
%   identically), ASCII punctuation/digits and a small typographic set
%   squeezed; a word containing any letter that cannot be faithfully folded
%   (Cyrillic, Greek, un-tabled Latin, orphan combining marks) is DROPPED
%   whole, and the drop count is reported unconditionally on stderr (INV-3:
%   dropped words are a compromise, never verbose-gated). Pure-ASCII
%   dictionaries normalize exactly as before this hardening. Policy + review
%   record: docs/plans/fill-dict-unicode-normalization.md.
%   Benchmark seam: the fill bench's `dict_load` bucket (load_inf) times
%   exactly this goal - it is the gated load path (F-L1).
load_dict(File, DictByLen, Index) :-
    load_dict(File, [], DictByLen, Index, _Scores).

% load_dict(+File, +Options, -DictByLen, -Index, -Scores) is det (internal).
% The §8.4a scored-ingestion entry the product path uses. Options:
%   min_score(N) - the hard candidate prune (default 1, DP-5: exclude only
%                  the score-0 blocklist floor; N is in the dict's NATIVE
%                  units, no normalisation).
% Scores is `uniform` (plain wordlist: every word scores 1, DP-5) or
% scores(Assoc) (assoc: normalized letter-list -> score, PRE-prune, so the
% fill-quality report can score any placed word - incl. a pruned seed - the
% way score_fill.py would). Dispatch is per FILE, not per line: a file with
% no `;` at all is a plain wordlist and takes the pre-scoring path verbatim
% (one C-level sub_string scan decides), keeping AC-FILL-6's byte-identity
% and the gated load_inf counts on every plain-dict rung. A file with any
% `;` is scored: `word;score` lines carry an integer score >= 0 in native
% units; `;`-less lines in it are unscored (uniform score 1); a line whose
% score part is not a non-negative integer is dropped + counted + reported
% (INV-3, AC-FILL-8); a duplicate word keeps its maximum score (the same
% silent, deterministic dedupe class as the plain path's sort/2).
load_dict(File, Options, DictByLen, Index, Scores) :-
    option(min_score(MinScore), Options, 1),
    read_file_lines(File, Lines, Str),
    (   MinScore =< 1,
        \+ sub_string(Str, _, _, _, ";")
    ->  plain_dict_words(Lines, File, Words),
        Scores = uniform
    ;   scored_dict_words(Lines, File, MinScore, Words, ScoreAssoc),
        Scores = scores(ScoreAssoc)
    ),
    index_words(Words, DictByLen, Index).

% The pre-scoring plain-wordlist path, byte-for-byte (the gated F-L1 body).
plain_dict_words(Lines, File, Words) :-
    % convlist, not findall+member+guard (C32): the same word list in the same
    % order with no findall record/copy per word. NB named helper, not a yall
    % lambda (the F-L1 rationale, fold_chars/3 below) - this is the gated
    % load_inf path. The drop counter is nb_ state (zeroed here, bumped only
    % on the drop path, read once after): zero inference cost per kept word.
    nb_setval(fill_dict_dropped_words, 0),
    convlist(line_word, Lines, Ws0),
    nb_getval(fill_dict_dropped_words, Dropped),
    warn_dropped_words(Dropped, File),
    sort(Ws0, Sorted),                      % dedupe + a stable, deterministic order
    seed_perturb_plain(Sorted, Words).

% §8.4b (AC-FILL-10): permute the whole plain word list - a uniform
% (unscored) dict's "equal-score tie" is the entire length bucket, and
% bucket order inherits Words' order (keysort/2 is stable). Seeded: the
% user's PRNG stream (the DP-6 semantics). Default path (§8.4c/C3): the
% SAME permutation on the PINNED engine constant - still a pure function of
% the input, so determinism (AC-FILL-13 as amended) holds. This permutation
% is load-bearing for search power, not cosmetics: a lexicographic bucket's
% prefix is an alphabetical clump (AAH/AAL/AAS...) with near-zero letter
% diversity, and the C3 campaign showed clumped value order defeats the
% restart portfolio on open grids over big uniform dicts (g17_50k) where
% permuted order converges in a handful of attempts. Scored dicts are NOT
% touched on the default path - §8.4a's score-desc-then-lex order is the
% quality contract. NB the default path uses an O(n log n) draw-key shuffle,
% not seeded_permutation/2's O(n^2) selection walk (minutes at ENABLE
% scale); the seeded path keeps the O(n^2) form byte-for-byte (its draws
% are contract - DP-6 goldens/fuzz reproduce per seed).
seed_perturb_plain(Ws, Out) :-
    (   current_search_seed(_)
    ->  seeded_permutation(Ws, Out)
    ;   prng_seed_state(0x9E3779B9),     % pinned load stream (Weyl constant)
        prng_shuffle(Ws, Out)
    ).

% O(n log n) pinned shuffle: one prng_draw key per word, keysort, strip.
% keysort/2 is stable and 64-bit key collisions are ~n^2/2^64, so the result
% is deterministic on every build (GMP and LibBF alike).
prng_shuffle(Ws, Out) :-
    maplist(prng_key_pair, Ws, KWs),
    keysort(KWs, Sorted),
    pairs_values(Sorted, Out).
prng_key_pair(W, K-W) :- prng_draw(K).

% Scored ingestion (§8.4a): parse `word;score` lines, dedupe to the max
% score, prune below MinScore (reported, INV-3), and order every survivor
% score-DESCENDING then lexicographic - the §8.4a total order. The ordering
% lands here (not in the search): candidate enumeration follows the length
% buckets' internal order, so a reordered bucket IS the new candidate order
% and the MRV search itself needs zero edits. For a uniform-score file the
% order collapses to the plain path's lexicographic sort (AC-FILL-6).
scored_dict_words(Lines, File, MinScore, Words, ScoreAssoc) :-
    nb_setval(fill_dict_dropped_words, 0),
    nb_setval(fill_dict_malformed_lines, 0),
    convlist(line_word_score, Lines, Pairs0),
    nb_getval(fill_dict_dropped_words, Dropped),
    warn_dropped_words(Dropped, File),
    nb_getval(fill_dict_malformed_lines, Malformed),
    warn_malformed_lines(Malformed, File),
    msort(Pairs0, Sorted),                  % word, then score ASCENDING
    dedupe_max_score(Sorted, Pairs),        % duplicate word -> its max score
    list_to_assoc(Pairs, ScoreAssoc),
    prune_min_score(Pairs, MinScore, File, Kept),
    order_score_desc(Kept, Words).

% Shared tail: length buckets + the pattern index. keysort/2 is STABLE, so
% each length bucket inherits Words' relative order - lexicographic on the
% plain path, score-descending-then-lex on the scored path.
index_words(Words, DictByLen, Index) :-
    map_list_to_pairs(length, Words, LPairs),
    keysort(LPairs, LSorted),
    group_pairs_by_key(LSorted, LGroups),
    list_to_assoc(LGroups, DictByLen),
    build_index(DictByLen, Index).

% encoding(utf8), NOT the locale default: under LC_ALL=C the default decode
% mangles multibyte input into U+FFFD mojibake, making the loaded dictionary
% depend on the invoking environment. Dictionary files are DEFINED as UTF-8
% (fixtures/dict/README.md); pure-ASCII files decode identically either way.
% /3 also hands back the raw file string so load_dict/5's plain-vs-scored
% dispatch is one whole-string scan, not a second file read.
read_file_lines(File, Lines) :-
    read_file_lines(File, Lines, _Str).
read_file_lines(File, Lines, Str) :-
    read_file_to_string(File, Str, [encoding(utf8)]),
    split_string(Str, "\n", "\r \t", Parts),
    exclude(==(""), Parts, Lines).

% One dictionary line -> its normalized letters. Two distinct failure kinds,
% both making convlist skip the line:
%   - nothing left after squeezing (punctuation/digit-only line): plain skip,
%     uncounted - exactly the pre-hardening behavior;
%   - normalize_word/2 FAILED (a letter that cannot be faithfully folded to
%     A-Z): the whole word is dropped and counted; load_dict reports the
%     count unconditionally on stderr per INV-3 (see core.pl's verbose_report
%     contract note: dropped words must NOT go through verbose_report).
line_word(Line, Letters) :-
    (   normalize_word(Line, Letters0)
    ->  Letters0 \== [],
        Letters = Letters0
    ;   dict_drop_note,
        fail
    ).

% Bump the dropped-word counter. nb_ globals, not assert/retract: runs ONLY
% on the (rare) drop path, so the ASCII fast path pays nothing; WASM-safe.
% The load-time init below makes the global always exist (load_dict re-zeroes
% it per call); without it a white-box line_word call would error on the
% first-ever drop.
:- nb_setval(fill_dict_dropped_words, 0).
dict_drop_note :-
    nb_getval(fill_dict_dropped_words, N0),
    N is N0 + 1,
    nb_setval(fill_dict_dropped_words, N).

% Unconditional stderr per INV-3 (a compromise report, like arrange's dropped
% words) - quiet exactly when nothing was dropped, so pure-ASCII dictionaries
% keep fill's quiet-by-default stderr contract.
warn_dropped_words(Dropped, File) :-
    (   Dropped =:= 0
    ->  true
    ;   format(user_error,
               "fill: dropped ~w dictionary word(s) from ~w (unrepresentable in A-Z after folding)~n",
               [Dropped, File])
    ).

% One `word;score` line -> Letters-Score. A `;`-less line inside a scored
% file is unscored: uniform score 1 (DP-5, "unrated" != "blocklisted"). A
% scored line splits at its LAST `;` (mirroring score_fill.py's rsplit); any
% earlier `;` stays in the word part, where normalize_word squeezes it like
% all ASCII punctuation. A score part that is not a non-negative integer
% drops the line, counted + reported (INV-3, AC-FILL-8) - unlike the word
% part, malformed metadata cannot be "folded" into anything faithful. The
% word part then goes through the exact plain-path line_word/2 (normalize,
% empty-skip, unfoldable-drop accounting).
line_word_score(Line, Letters-Score) :-
    split_string(Line, ";", " \t", Parts),
    (   Parts = [WordStr]
    ->  Score = 1
    ;   % first solution = the split at the LAST `;` is the single solution
        % with a 1-element suffix; once/1 prunes append/3's dead choicepoint.
        once(append(WordParts, [ScoreStr], Parts)),
        (   number_string(Score0, ScoreStr), integer(Score0), Score0 >= 0
        ->  Score = Score0,
            atomic_list_concat(WordParts, ';', WordStr)
        ;   dict_malformed_note,
            fail
        )
    ),
    line_word(WordStr, Letters).

% Malformed scored-line counter: same nb_ shape + load-time init as the
% dropped-words counter above (scored_dict_words re-zeroes it per load).
:- nb_setval(fill_dict_malformed_lines, 0).
dict_malformed_note :-
    nb_getval(fill_dict_malformed_lines, N0),
    N is N0 + 1,
    nb_setval(fill_dict_malformed_lines, N).

warn_malformed_lines(Malformed, File) :-
    (   Malformed =:= 0
    ->  true
    ;   format(user_error,
               "fill: dropped ~w malformed scored line(s) from ~w (expected word;integer-score, score >= 0)~n",
               [Malformed, File])
    ).

% msort/2 sorts Letters-Score pairs by word then score ASCENDING, so within
% a duplicate-word group the LAST pair carries the maximum score - keep it.
% Deterministic and input-order-independent.
dedupe_max_score([], []).
dedupe_max_score([W-S|T], Out) :-
    (   T = [W2-_|_], W == W2
    ->  dedupe_max_score(T, Out)
    ;   Out = [W-S|Rest],
        dedupe_max_score(T, Rest)
    ).

% The §8.4a hard prune: drop every word scoring < MinScore, BEFORE the
% bucket/index build (so slot domains never see a pruned word, AC-FILL-5).
% Any nonzero prune is reported unconditionally on stderr (INV-3) - the
% DP-5 default (1) prunes nothing from a dict without a score-0 floor, so
% the quiet-success stderr contract holds wherever it held before. A prune
% that empties the whole dictionary additionally reports the observed
% maximum score - the targeted hint for `--min-score N` above a uniform
% (unscored, all-1) dictionary's ceiling.
prune_min_score(Pairs, MinScore, File, Kept) :-
    partition_min_score(Pairs, MinScore, Kept, NPruned),
    (   NPruned =:= 0
    ->  true
    ;   length(Pairs, Total),
        format(user_error,
               "fill: --min-score ~w pruned ~w of ~w dictionary word(s) from ~w~n",
               [MinScore, NPruned, Total, File]),
        (   Kept == []
        ->  max_pair_score(Pairs, 0, Max),
            format(user_error,
                   "fill: --min-score ~w exceeds the dictionary's maximum score ~w (unscored words score 1); no words remain~n",
                   [MinScore, Max])
        ;   true
        )
    ).

partition_min_score([], _, [], 0).
partition_min_score([W-S|T], Min, Kept, N) :-
    (   S >= Min
    ->  Kept = [W-S|K1],
        partition_min_score(T, Min, K1, N)
    ;   partition_min_score(T, Min, Kept, N1),
        N is N1 + 1
    ).

max_pair_score([], M, M).
max_pair_score([_-S|T], M0, M) :-
    M1 is max(M0, S),
    max_pair_score(T, M1, M).

% Score-descending, then a within-band shuffle (§8.4a's order as amended at
% §8.4c/C3): msort over k(-Score, Letters) keys - ascending -Score =
% descending score (the quality contract, untouched) - then every
% equal-score band is reordered by seed_tie_shuffle/2 below. Lexicographic
% order survives only as the pre-shuffle canonical form (it makes the
% shuffle input - and so the output - deterministic).
order_score_desc(Pairs, Words) :-
    maplist(neg_score_key, Pairs, Keyed),
    msort(Keyed, SortedKeys),
    seed_tie_shuffle(SortedKeys, Perturbed),
    maplist(key_word, Perturbed, Words).

neg_score_key(W-S, k(NS, W)) :- NS is -S.
key_word(k(_, W), W).

% §8.4b (AC-FILL-10): reorder WITHIN each equal-score group only -
% score-descending stays the primary order (the §8.4a quality contract),
% only the lexicographic tiebreak varies. Seeded: the user's PRNG stream
% (DP-6 semantics, byte-for-byte). Default path (§8.4c/C3): pinned shuffle
% IFF the dict is a single band (all scores equal) - i.e. exactly when it
% carries no ordering information and is semantically a plain list, which
% keeps AC-FILL-6's uniform-dict identity (same input order, same reseed,
% same draw-per-element walk as seed_perturb_plain/2's whole-list shuffle).
% A dict with REAL bands keeps the strict §8.4a score-desc-then-lex order:
% the C3 campaign measured both directions - lex order defeats the restart
% portfolio on big uniform dicts (g17_50k), yet on the scored reference row
% (blocked_13a @30) it is the lex bands that converge (the DP-7 probe's 3/3
% fast seeds ran top3 over lex bands) while in-band permutation starves the
% budget (seed-7 and pinned-shuffle runs both died). Shuffle where order is
% meaningless, preserve it where it is contract.
seed_tie_shuffle(Keyed, Perturbed) :-
    (   current_search_seed(_)
    ->  maplist(ns_group_key, Keyed, NSPairs),
        group_pairs_by_key(NSPairs, Groups),
        maplist(shuffle_score_group, Groups, Buckets),
        append(Buckets, Perturbed)
    ;   maplist(ns_group_key, Keyed, NSPairs),
        group_pairs_by_key(NSPairs, Groups),
        (   Groups = [_-OnlyBand]
        ->  prng_seed_state(0x9E3779B9), % pinned load stream (Weyl constant)
            prng_shuffle(OnlyBand, Perturbed)
        ;   Perturbed = Keyed
        )
    ).

ns_group_key(k(NS, W), NS-k(NS, W)).
shuffle_score_group(_NS-Group, Shuffled) :- seeded_permutation(Group, Shuffled).

% NB first-order on purpose: an include([C]>>char_type(C,alpha), ...) filter
% pays the yall meta-call + lambda copy PER CHARACTER (this module never
% imports library(yall), so the lambda is not compile-expanded): ~12 inf/char
% = 19.4M of full ENABLE's 26.6M dict-load inferences (P-F1 attribution).
% The az/1 first-arg-indexed fact lookup makes the keep/other decision per
% char at the same measured ~2.1 inf/char as the char_type(C, alpha) scan it
% replaced (F-L1: load_inf 26.60M -> 10.72M at 172k, output-identical), and
% unlike char_type it is LOCALE-INDEPENDENT: both string_upper/2 and
% char_type/2 classify non-ASCII differently under LC_ALL=C vs UTF-8 locales,
% so no non-ASCII decision below consults either (string_upper is retained
% only for its locale-stable ASCII a-z -> A-Z mapping; the fold table is
% double-cased to cover non-ASCII letters whether or not the locale's
% string_upper case-mapped them).
normalize_word(Line, Letters) :-
    string_upper(Line, U), string_chars(U, Cs),
    fold_chars(Cs, none, Letters).

% fold_chars(+Chars, +PrevKept, -Letters): fold-then-strict walk. PrevKept is
% the immediately preceding kept A-Z char (the only legal combining-mark
% base), or none - reset by every char that is not a kept A-Z letter, so mark
% chains never squeeze. Together with the PAIR-keyed allowed_mark/2 test this
% makes NFC and NFD forms agree for EVERY input class: a decomposed sequence
% squeezes its mark only when the recomposed letter is in the fold table
% (allowed_mark holds the (base, mark) pairs of the table's own NFD
% decompositions), and anything else - orphan marks, mark chains, marks on
% the wrong base - drops the word exactly as its precomposed form would.
% Fails (word dropped) on any char that is neither kept, folded, a legal
% mark, nor squeezable punctuation.
fold_chars([], _, []).
fold_chars([C|Cs], Prev, Letters) :-
    (   az(C)
    ->  Letters = [C|Rest],
        fold_chars(Cs, C, Rest)
    ;   fold_other(C, Prev, Cs, Letters)
    ).

% Non-A-Z char (rare path; never entered for pure-ASCII dictionaries).
fold_other(C, Prev, Cs, Letters) :-
    (   fold_char(C, F)
    ->  append(F, Rest, Letters),
        fold_chars(Cs, none, Rest)
    ;   combining_mark(C)
    ->  Prev \== none,
        allowed_mark(Prev, C),
        fold_chars(Cs, none, Letters)
    ;   squeeze_char(C),
        fold_chars(Cs, none, Letters)
    ).

combining_mark(C) :-
    char_code(C, Code),
    Code >= 0x0300, Code =< 0x036F.

% ASCII non-letters (digits, punctuation - a-z cannot appear post-upper) plus
% an explicit typographic set squeeze silently, matching the pre-hardening
% punctuation behavior ("don't" with either the ASCII or the U+2019 apostrophe -> DONT).
squeeze_char(C) :-
    char_code(C, Code),
    (   Code < 0x80
    ->  true
    ;   typographic_squeeze(Code)
    ).

typographic_squeeze(0x00AD).    % soft hyphen
typographic_squeeze(0x2013).    % en dash
typographic_squeeze(0x2014).    % em dash
typographic_squeeze(0x2018).    % left single quote
typographic_squeeze(0x2019).    % right single quote (typographic apostrophe)

% The A-Z fast path: 26 first-arg-indexed facts, measured inference-identical
% to the char_type(C, alpha) call they replace (see normalize_word note).
az('A'). az('B'). az('C'). az('D'). az('E'). az('F'). az('G'). az('H').
az('I'). az('J'). az('K'). az('L'). az('M'). az('N'). az('O'). az('P').
az('Q'). az('R'). az('S'). az('T'). az('U'). az('V'). az('W'). az('X').
az('Y'). az('Z').

% The fold table: every NFD-decomposable Latin-1 / Latin Extended-A letter ->
% its A-Z base, DOUBLE-CASED (both the upper- and lower-case precomposed char
% appear, so the outcome does not depend on whether the locale's string_upper
% case-mapped it), plus five explicitly anchored specials: sharp-s->SS (Unicode
% casefold), AE/OE ligatures (established English orthographic equivalence),
% IJ ligature (NFKC), O-slash->O. Thorn/Eth are EXCLUDED by
% review: no normative crossword-convention source for a transliteration, so
% words containing them drop. GENERATED from library(unicode) (utf8proc,
% Unicode 16) by the plan's generator script and verified against the same
% library by the fill.plt oracle test - never hand-edit an entry; regenerate.
% \uXXXX escapes keep every FACT byte-ASCII: a table entry can never be
% corrupted by source-encoding mishandling, and diffs stay greppable.
fold_char('\u00c0', ['A']).
fold_char('\u00c1', ['A']).
fold_char('\u00c2', ['A']).
fold_char('\u00c3', ['A']).
fold_char('\u00c4', ['A']).
fold_char('\u00c5', ['A']).
fold_char('\u00c6', ['A','E']).
fold_char('\u00c7', ['C']).
fold_char('\u00c8', ['E']).
fold_char('\u00c9', ['E']).
fold_char('\u00ca', ['E']).
fold_char('\u00cb', ['E']).
fold_char('\u00cc', ['I']).
fold_char('\u00cd', ['I']).
fold_char('\u00ce', ['I']).
fold_char('\u00cf', ['I']).
fold_char('\u00d1', ['N']).
fold_char('\u00d2', ['O']).
fold_char('\u00d3', ['O']).
fold_char('\u00d4', ['O']).
fold_char('\u00d5', ['O']).
fold_char('\u00d6', ['O']).
fold_char('\u00d8', ['O']).
fold_char('\u00d9', ['U']).
fold_char('\u00da', ['U']).
fold_char('\u00db', ['U']).
fold_char('\u00dc', ['U']).
fold_char('\u00dd', ['Y']).
fold_char('\u00df', ['S','S']).
fold_char('\u00e0', ['A']).
fold_char('\u00e1', ['A']).
fold_char('\u00e2', ['A']).
fold_char('\u00e3', ['A']).
fold_char('\u00e4', ['A']).
fold_char('\u00e5', ['A']).
fold_char('\u00e6', ['A','E']).
fold_char('\u00e7', ['C']).
fold_char('\u00e8', ['E']).
fold_char('\u00e9', ['E']).
fold_char('\u00ea', ['E']).
fold_char('\u00eb', ['E']).
fold_char('\u00ec', ['I']).
fold_char('\u00ed', ['I']).
fold_char('\u00ee', ['I']).
fold_char('\u00ef', ['I']).
fold_char('\u00f1', ['N']).
fold_char('\u00f2', ['O']).
fold_char('\u00f3', ['O']).
fold_char('\u00f4', ['O']).
fold_char('\u00f5', ['O']).
fold_char('\u00f6', ['O']).
fold_char('\u00f8', ['O']).
fold_char('\u00f9', ['U']).
fold_char('\u00fa', ['U']).
fold_char('\u00fb', ['U']).
fold_char('\u00fc', ['U']).
fold_char('\u00fd', ['Y']).
fold_char('\u00ff', ['Y']).
fold_char('\u0100', ['A']).
fold_char('\u0101', ['A']).
fold_char('\u0102', ['A']).
fold_char('\u0103', ['A']).
fold_char('\u0104', ['A']).
fold_char('\u0105', ['A']).
fold_char('\u0106', ['C']).
fold_char('\u0107', ['C']).
fold_char('\u0108', ['C']).
fold_char('\u0109', ['C']).
fold_char('\u010a', ['C']).
fold_char('\u010b', ['C']).
fold_char('\u010c', ['C']).
fold_char('\u010d', ['C']).
fold_char('\u010e', ['D']).
fold_char('\u010f', ['D']).
fold_char('\u0112', ['E']).
fold_char('\u0113', ['E']).
fold_char('\u0114', ['E']).
fold_char('\u0115', ['E']).
fold_char('\u0116', ['E']).
fold_char('\u0117', ['E']).
fold_char('\u0118', ['E']).
fold_char('\u0119', ['E']).
fold_char('\u011a', ['E']).
fold_char('\u011b', ['E']).
fold_char('\u011c', ['G']).
fold_char('\u011d', ['G']).
fold_char('\u011e', ['G']).
fold_char('\u011f', ['G']).
fold_char('\u0120', ['G']).
fold_char('\u0121', ['G']).
fold_char('\u0122', ['G']).
fold_char('\u0123', ['G']).
fold_char('\u0124', ['H']).
fold_char('\u0125', ['H']).
fold_char('\u0128', ['I']).
fold_char('\u0129', ['I']).
fold_char('\u012a', ['I']).
fold_char('\u012b', ['I']).
fold_char('\u012c', ['I']).
fold_char('\u012d', ['I']).
fold_char('\u012e', ['I']).
fold_char('\u012f', ['I']).
fold_char('\u0130', ['I']).
fold_char('\u0132', ['I','J']).
fold_char('\u0133', ['I','J']).
fold_char('\u0134', ['J']).
fold_char('\u0135', ['J']).
fold_char('\u0136', ['K']).
fold_char('\u0137', ['K']).
fold_char('\u0139', ['L']).
fold_char('\u013a', ['L']).
fold_char('\u013b', ['L']).
fold_char('\u013c', ['L']).
fold_char('\u013d', ['L']).
fold_char('\u013e', ['L']).
fold_char('\u0143', ['N']).
fold_char('\u0144', ['N']).
fold_char('\u0145', ['N']).
fold_char('\u0146', ['N']).
fold_char('\u0147', ['N']).
fold_char('\u0148', ['N']).
fold_char('\u014c', ['O']).
fold_char('\u014d', ['O']).
fold_char('\u014e', ['O']).
fold_char('\u014f', ['O']).
fold_char('\u0150', ['O']).
fold_char('\u0151', ['O']).
fold_char('\u0152', ['O','E']).
fold_char('\u0153', ['O','E']).
fold_char('\u0154', ['R']).
fold_char('\u0155', ['R']).
fold_char('\u0156', ['R']).
fold_char('\u0157', ['R']).
fold_char('\u0158', ['R']).
fold_char('\u0159', ['R']).
fold_char('\u015a', ['S']).
fold_char('\u015b', ['S']).
fold_char('\u015c', ['S']).
fold_char('\u015d', ['S']).
fold_char('\u015e', ['S']).
fold_char('\u015f', ['S']).
fold_char('\u0160', ['S']).
fold_char('\u0161', ['S']).
fold_char('\u0162', ['T']).
fold_char('\u0163', ['T']).
fold_char('\u0164', ['T']).
fold_char('\u0165', ['T']).
fold_char('\u0168', ['U']).
fold_char('\u0169', ['U']).
fold_char('\u016a', ['U']).
fold_char('\u016b', ['U']).
fold_char('\u016c', ['U']).
fold_char('\u016d', ['U']).
fold_char('\u016e', ['U']).
fold_char('\u016f', ['U']).
fold_char('\u0170', ['U']).
fold_char('\u0171', ['U']).
fold_char('\u0172', ['U']).
fold_char('\u0173', ['U']).
fold_char('\u0174', ['W']).
fold_char('\u0175', ['W']).
fold_char('\u0176', ['Y']).
fold_char('\u0177', ['Y']).
fold_char('\u0178', ['Y']).
fold_char('\u0179', ['Z']).
fold_char('\u017a', ['Z']).
fold_char('\u017b', ['Z']).
fold_char('\u017c', ['Z']).
fold_char('\u017d', ['Z']).
fold_char('\u017e', ['Z']).

% allowed_mark(BaseUpper, Mark): the (base, mark) NFD-decomposition pairs of
% the fold table's single-base entries - the ONLY combining marks fold_chars
% squeezes, and only immediately after that base. Generated + oracle-checked
% alongside fold_char/2.
allowed_mark('A', '\u0300').
allowed_mark('A', '\u0301').
allowed_mark('A', '\u0302').
allowed_mark('A', '\u0303').
allowed_mark('A', '\u0304').
allowed_mark('A', '\u0306').
allowed_mark('A', '\u0308').
allowed_mark('A', '\u030a').
allowed_mark('A', '\u0328').
allowed_mark('C', '\u0301').
allowed_mark('C', '\u0302').
allowed_mark('C', '\u0307').
allowed_mark('C', '\u030c').
allowed_mark('C', '\u0327').
allowed_mark('D', '\u030c').
allowed_mark('E', '\u0300').
allowed_mark('E', '\u0301').
allowed_mark('E', '\u0302').
allowed_mark('E', '\u0304').
allowed_mark('E', '\u0306').
allowed_mark('E', '\u0307').
allowed_mark('E', '\u0308').
allowed_mark('E', '\u030c').
allowed_mark('E', '\u0328').
allowed_mark('G', '\u0302').
allowed_mark('G', '\u0306').
allowed_mark('G', '\u0307').
allowed_mark('G', '\u0327').
allowed_mark('H', '\u0302').
allowed_mark('I', '\u0300').
allowed_mark('I', '\u0301').
allowed_mark('I', '\u0302').
allowed_mark('I', '\u0303').
allowed_mark('I', '\u0304').
allowed_mark('I', '\u0306').
allowed_mark('I', '\u0307').
allowed_mark('I', '\u0308').
allowed_mark('I', '\u0328').
allowed_mark('J', '\u0302').
allowed_mark('K', '\u0327').
allowed_mark('L', '\u0301').
allowed_mark('L', '\u030c').
allowed_mark('L', '\u0327').
allowed_mark('N', '\u0301').
allowed_mark('N', '\u0303').
allowed_mark('N', '\u030c').
allowed_mark('N', '\u0327').
allowed_mark('O', '\u0300').
allowed_mark('O', '\u0301').
allowed_mark('O', '\u0302').
allowed_mark('O', '\u0303').
allowed_mark('O', '\u0304').
allowed_mark('O', '\u0306').
allowed_mark('O', '\u0308').
allowed_mark('O', '\u030b').
allowed_mark('R', '\u0301').
allowed_mark('R', '\u030c').
allowed_mark('R', '\u0327').
allowed_mark('S', '\u0301').
allowed_mark('S', '\u0302').
allowed_mark('S', '\u030c').
allowed_mark('S', '\u0327').
allowed_mark('T', '\u030c').
allowed_mark('T', '\u0327').
allowed_mark('U', '\u0300').
allowed_mark('U', '\u0301').
allowed_mark('U', '\u0302').
allowed_mark('U', '\u0303').
allowed_mark('U', '\u0304').
allowed_mark('U', '\u0306').
allowed_mark('U', '\u0308').
allowed_mark('U', '\u030a').
allowed_mark('U', '\u030b').
allowed_mark('U', '\u0328').
allowed_mark('W', '\u0302').
allowed_mark('Y', '\u0301').
allowed_mark('Y', '\u0302').
allowed_mark('Y', '\u0308').
allowed_mark('Z', '\u0301').
allowed_mark('Z', '\u0307').
allowed_mark('Z', '\u030c').

build_index(DictByLen, Index) :-
    findall(k(Len, P, Ch)-Idx,
            ( gen_assoc(Len, DictByLen, Words),
              nth0(Idx, Words, W), nth0(P, W, Ch) ),
            Triples),
    keysort(Triples, Sorted),
    group_pairs_by_key(Sorted, Grouped),
    % NB named helper, not a yall lambda (the F-L1 rationale, alpha_chars/2
    % above): un-imported `>>` is interpreted - a meta-call + lambda copy per
    % index key group, on the load_inf-gated dict-load path.
    maplist(pair_ordset, Grouped, GroupedSets),
    list_to_assoc(GroupedSets, Index).

pair_ordset(K-Is, K-Set) :- list_to_ord_set(Is, Set).

% The slot's bound positions as ascending P-V pairs. NB first-order on purpose
% (the F-L1 rationale, alpha_chars/2 above): the previous
% findall(P-V, (nth0(P, Vars, V), nonvar(V)), Bound) paid the findall
% record/collect/copy machinery plus nth0's enumerate-by-backtracking
% choicepoint per cell, at EVERY recount at EVERY search node - the hottest
% kernel in the engine. This positional walk makes the identical nonvar/1
% decision per cell; Bound's content and order (ascending P, V shared not
% copied - always a ground char atom here) are identical by construction.
% Shared by slot_bucket/5 and the masks kernel of candidate_count/5.
bound_positions(Vars, Bound) :-
    bound_positions(Vars, 0, Bound).

bound_positions([], _, []).
bound_positions([V|Vs], P, Bound) :-
    (   nonvar(V)
    ->  Bound = [P-V|Rest]
    ;   Bound = Rest
    ),
    P1 is P + 1,
    bound_positions(Vs, P1, Rest).

% Resolve a slot to its length-bucket Words plus a selector Sel over that bucket:
% `all` (no cell bound yet - every word of the length is a candidate) or
% idx(Indices) (the ordset of bucket indices matching the bound cells). Shared by
% candidates/4 (which materializes the words) and candidate_count/4 (which needs
% only the size), so the Bound + index-intersection work is written once.
slot_bucket(Vars, DictByLen, Index, Words, Sel) :-
    length(Vars, Len),
    ( get_assoc(Len, DictByLen, Words) -> true ; Words = [] ),
    bound_positions(Vars, Bound),
    ( Bound == []
    ->  Sel = all
    ;   index_intersection(Bound, Len, Index, Indices), Sel = idx(Indices)
    ).

% Candidate words for a slot, given its currently-bound positions.
candidates(Vars, DictByLen, Index, Cands) :-
    slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> Cands = Words
    ; Sel = idx(Indices), maplist(nth0_of(Words), Indices, Cands)
    ).

nth0_of(Words, I, W) :- nth0(I, Words, W).

% Number of candidate words for a slot WITHOUT materializing the word list (P3).
% The MRV metric only needs the size, so count the matching index ordset (or the
% whole length bucket) directly and skip the per-slot maplist(nth0_of...) word
% build that select_mrv/6 otherwise ran for every unfilled slot at every search
% node. |Indices| == |maplist(nth0_of(Words), Indices)| by construction (every
% index into Words is valid), so the count is identical to length(Cands).
%
% F-H2 (bitset counting via artifact v2). candidate_count/5 dispatches on a
% threaded Masks context (see the search below): `none` keeps the ordset kernel
% (raw text mode and every non-artifact / v1-shaped path - byte-identical); a
% masks(MaskAssoc) context (present ONLY in a v2 artifact) counts the SAME slot
% with bignum `/\` chains + popcount over per-(len,pos,char) bit masks. The masks
% are derived FROM these very ordsets at --save-index time (bit i == bucket index
% i), so popcount(chain) == |ord_intersection(chain)| for every pattern by
% construction. Enumeration (candidates/4) stays on ordsets everywhere; only the
% count seam changes. candidate_count/4 is the ordset reference retained for the
% white-box tests (P3) and select_mrv/6.
candidate_count(Vars, DictByLen, Index, Count) :-
    candidate_count(none, Vars, DictByLen, Index, Count).

% The Masks context is ARGUMENT 1 so SWI's first-argument indexing dispatches
% `none` vs masks(_) deterministically - no choicepoint, no cut (on this SWI a
% non-first-argument scan does NOT prune the choicepoint, so the old arg-4
% shape needed a leading cut and still paid the CP-push per call). The clauses
% are disjoint by construction, and the reference selector select_mrv/6 stays
% choicepoint-free (P13).
%
% ORDSET kernel (Masks == none): identical work to the pre-F-H2 candidate_count/4.
candidate_count(none, Vars, DictByLen, Index, Count) :-
    slot_bucket(Vars, DictByLen, Index, Words, Sel),
    ( Sel == all -> length(Words, Count)
    ; Sel = idx(Indices), length(Indices, Count)
    ).
% BIGNUM kernel (Masks == masks(MaskAssoc), v2 artifact mode). Mirrors
% slot_bucket's branch structure EXACTLY - the `all` branch (no cell bound) still
% counts the whole length bucket; the idx branch AND-folds the bound cells' masks
% and popcounts. No ord_intersection is computed on this path (that is the win).
candidate_count(masks(MaskAssoc), Vars, DictByLen, _Index, Count) :-
    length(Vars, Len),
    bound_positions(Vars, Bound),
    ( Bound == []
    ->  ( get_assoc(Len, DictByLen, Words) -> length(Words, Count) ; Count = 0 )
    ;   mask_count(Bound, Len, MaskAssoc, Count)
    ).

% popcount(AND over the bound cells' masks). Mirrors index_intersection/4: the
% first bound cell seeds the accumulator, the rest AND into it; an absent key is
% mask 0 (mirrors index_set's `S = []`), so a dead cell zeroes the chain and the
% count is 0 - exactly as ord_intersection with [] yields the empty set.
mask_count([P0-Ch0|Rest], Len, MaskAssoc, Count) :-
    mask_set(Len, P0, Ch0, MaskAssoc, M0),
    foldl(mask_and(Len, MaskAssoc), Rest, M0, MAnd),
    Count is popcount(MAnd).
mask_and(Len, MaskAssoc, P-Ch, Acc, Acc1) :-
    mask_set(Len, P, Ch, MaskAssoc, M), Acc1 is Acc /\ M.
mask_set(Len, P, Ch, MaskAssoc, M) :-
    ( get_assoc(k(Len, P, Ch), MaskAssoc, M0) -> M = M0 ; M = 0 ).

index_intersection([P-Ch|Rest], Len, Index, Indices) :-
    index_set(Len, P, Ch, Index, S0),
    foldl(index_intersect(Len, Index), Rest, S0, Indices).
index_intersect(Len, Index, P-Ch, Acc, Acc1) :-
    index_set(Len, P, Ch, Index, S), ord_intersection(Acc, S, Acc1).
index_set(Len, P, Ch, Index, S) :-
    ( get_assoc(k(Len, P, Ch), Index, S0) -> S = S0 ; S = [] ).


% --- the MRV backtracking search ---------------------------------------------
% Fill all slots. At each step take the unfilled slot with the FEWEST matching
% candidates (most-constrained-first; ties broken by lowest start cell, for
% determinism), try each candidate (in dictionary order), unify it into the grid
% (binding the crossing cells), and recurse. Used answers are not repeated.
% F-H1 (incremental candidate counts). A slot's candidate count is a pure
% function of its cells' binding state, and placing a word binds ONLY the
% previously-free cells of the placed slot. So instead of recounting EVERY
% unfilled slot at EVERY node (select_mrv's per-node full recount, measured at
% 79.5-85.6% of search_inf in P-F1), carry the per-slot counts as backtrack-
% restored threaded state: after a placement recount EXACTLY the slots that
% cross a NEWLY-BOUND cell of the placed word; every other count carries over
% unchanged. Backtracking restores the previous counts for free (pure threaded
% state - the caller's structure is untouched, and the carried terms are built
% by unification, never findall/copy, so the shared crossing variables survive).
%
% EQUIVALENCE (the tree is preserved node-for-node):
%  (1) COUNTS EXACT BY CONSTRUCTION. The initial counts are a full candidate_count
%      per slot; after each placement, a slot whose cells did not change keeps its
%      exact count (candidate_count is a pure function of the cells' binding state)
%      and a slot crossing a newly-bound cell is recounted with the same
%      candidate_count/4. Counts NEVER read Used (the \+ memberchk(Word, Used)
%      filter is try-time only, below), so they are exact regardless of which
%      words are already placed.
%  (2) SAME SELECTION. select_min_count picks the minimum c(Count, Start, Dir) in
%      standard order - byte-identical to select_mrv's sort(0, @=<, ...) + head.
%      Start+Dir uniquely identify a slot, so this is a total order: no genuine
%      ties, and the winner is independent of the carried list's order.
%  (3) SAME WINNER MATERIALIZATION. The winner's candidate list is built by the
%      same candidates/4.
%  (4) COMPLETED SLOTS STAY IN THE SET. A slot fully bound by crossings (count 0
%      or 1) is NOT dropped, so a 0-count dead slot is still selected first (it
%      has the lowest count) and fails the branch exactly as today.
% Masks is the F-H2 counting context threaded down the whole search: `none` (raw
% text / v1 artifact - ordset kernel, byte-identical to pre-F-H2) or
% masks(MaskAssoc) (a v2 artifact - bignum popcount kernel). Threading a ground
% atom/compound through the recursion adds ZERO inferences (SWI clause indexing
% dispatches candidate_count/5 with no penalty), so the `none` path is
% inference-identical to the pre-F-H2 engine.
% §8.4c (DP-8): fill_search/5 IS the MAC + dom/wdeg search core. The former
% MRV/backtracking loops (fill_search_inc/_seeded, below) are RETAINED as
% reference predicates (select_mrv precedent): white-box tests exercise them
% and they document the §8.4 tree this core replaced at the DP-8 version
% bump. Masks: a v2 artifact's masks(MA) is reused directly; `none` (raw
% text path) builds the same structure from the Index here — the core NEEDS
% the letter masks (they are its domains), so the F-H2 "masks optional"
% counting split no longer applies to the search proper.
fill_search(Slots, DictByLen, Index, Masks, Used) :-
    (   Slots == []
    ->  true
    ;   mac_search_top(Slots, DictByLen, Index, Masks, Used)
    ).

% --- the §8.4c search core (MAC + dom/wdeg + restarts; DP-8) -------------------
% Domains are bignum candidate masks over the §8.4a score-ordered per-length
% buckets (bit i = bucket index i, so ascending-bit materialization IS the
% score-desc-then-dictionary candidate order — the quality contract is
% structural). Placement ANDs each crossing slot's domain with the placed
% letter's mask and propagates letter-support to an AC-3 fixpoint; a wiped
% domain fails the branch immediately. Slot selection is dom/wdeg: every
% crossing carries a conflict weight (+1 when its revision wipes a domain),
% held in a per-search nb_setarg'd global so learning survives backtracking,
% and the next slot minimizes domain-size / Σ(weights to unfilled
% neighbours). Attempts run under a growing node cap (500 ×1.5, weights aged
% ×0.99 between attempts); a cap-abandoned attempt restarts with better
% weights, an attempt that exhausts WITHOUT hitting its cap is a completed
% search — genuine infeasibility (AC-FILL-14: the cap throw and plain
% failure are distinct by construction, so a timeout can never masquerade as
% a proof). Candidate picks (mac_attempt_pick/3): default-path attempt 1 is
% strict best-first (no PRNG), later default attempts and every seeded
% attempt use weighted top-3 promotion ([4,2,1] on the engine-internal
% PRNG) — default attempts on a PINNED constant stream, so the path stays a
% pure function of the input (AC-FILL-13 as amended); seed presence is read
% once per search.
% Search assigns words WITHOUT binding cell variables (pure mask
% bookkeeping); the winning assignment is bound onto the shared cell
% variables at the end (mac_bind_fill/2), which re-checks crossing
% consistency by unification as a hard invariant.
mac_search_top(Slots, DictByLen, Index, Masks, Used) :-
    % P2 guard (tests/fill.plt fill_propagates_genuine_error): a malformed
    % dict/index must THROW, never fail into `infeasible`. The pre-§8.4c
    % engine got this for free (get_assoc's C builtin type-checks); this
    % core's first structural touch is assoc_to_list/2, which fails SILENTLY
    % on a non-assoc — so the type check is explicit, with the same error
    % class the old engine raised.
    (   is_assoc(DictByLen) -> true
    ;   throw(error(type_error(btree, DictByLen), _))
    ),
    (   Masks = masks(MA0)
    ->  MA = MA0
    ;   is_assoc(Index) -> build_masks(Index, MA)
    ;   throw(error(type_error(btree, Index), _))
    ),
    mac_flat_masks(MA, LmA),
    mac_buckets(DictByLen, BucketA),
    length(Slots, N), NMax is N - 1, numlist(0, NMax, Ids),
    pairs_keys_values(IdSlots, Ids, Slots),
    maplist(mac_slot_len, Slots, Lens),
    pairs_keys_values(IdLens, Ids, Lens),
    list_to_assoc(IdLens, LenA),
    mac_edges(IdSlots, EdgeA, NEdges),
    mac_init_weights(NEdges),
    maplist(mac_init_domain(DictByLen, LmA), IdSlots, IdDoms),
    % An initially-empty domain is a proof up front (no word matches the
    % slot's seed-bound letters): fail -> infeasible, named by empty_slots.
    \+ ( member(_-D0, IdDoms), D0 =:= 0 ),
    list_to_assoc(IdDoms, DomA0),
    % §8.4b dispatch, read ONCE per search: `seeded` = the user gave a seed
    % (--seed/--shuffle; every attempt is top3 on their PRNG stream);
    % `default` = attempt 1 is strict greedy, attempts >= 2 diversify on the
    % pinned engine constant (mac_attempt_pick below — the C3 discovery).
    (   current_search_seed(_) -> Mode = seeded ; Mode = default ),
    % Root AC fixpoint: failure here is likewise a proof (a wiped domain
    % with nothing placed).
    mac_propagate(Ids, DomA0, LenA, EdgeA, LmA, Ids, DomA1),
    mac_restart_loop(1, 500, Ids, DomA1, LenA, BucketA, EdgeA, LmA, Mode,
                     Used, Assign),
    mac_bind_fill(Assign, IdSlots).

mac_slot_len(slot(_, _, _, Vars), Len) :- length(Vars, Len).

% Bind the winning assignment onto the grid's shared cell variables. Words
% and Vars are both char lists, and crossing slots share cell variables, so
% this unification is also a full crossing-consistency check of the mask
% bookkeeping — it cannot fail unless the core is buggy, in which case
% failing (-> infeasible, never a bad emit) is the safe outcome.
mac_bind_fill([], _).
mac_bind_fill([Id-Word|As], IdSlots) :-
    memberchk(Id-slot(_, _, _, Vars), IdSlots),
    Vars = Word,
    mac_bind_fill(As, IdSlots).

% Per-length bucket compounds: arg/3 gives O(1) word lookup from a mask's
% bit index.
mac_buckets(DictByLen, BucketA) :-
    assoc_to_list(DictByLen, Pairs),
    maplist(mac_bucket_pair, Pairs, BPairs),
    list_to_assoc(BPairs, BucketA).
mac_bucket_pair(L-Ws, L-T) :- T =.. [b|Ws].

% Flat letter-mask table: assoc (Len-Pos) -> lm(MA,...,MZ), built once per
% search from build_masks/2's k(Len,Pos,Ch) assoc. arg/3 access replaces two
% logarithmic assoc lookups in the propagation hot loop (the naive form
% measured ~14ms/node on the 13x13/STW row in the DP-7 spike).
mac_flat_masks(MaskAssoc, LmA) :-
    assoc_to_list(MaskAssoc, Pairs),
    maplist(mac_lp_pair, Pairs, LPPairs),
    msort(LPPairs, Sorted),
    group_pairs_by_key(Sorted, Groups),
    maplist(mac_lm_compound, Groups, LmPairs),
    list_to_assoc(LmPairs, LmA).
mac_lp_pair(k(L, P, Ch)-M, (L-P)-(Ch-M)).
mac_lm_compound(LP-ChMs, LP-Lm) :-
    mac_letters(Chs),
    maplist(mac_lm_arg(ChMs), Chs, Ms),
    Lm =.. [lm|Ms].
mac_lm_arg(ChMs, Ch, M) :- ( memberchk(Ch-M0, ChMs) -> M = M0 ; M = 0 ).

mac_letters(['A','B','C','D','E','F','G','H','I','J','K','L','M',
             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']).

mac_lm_of(LmA, Len, Pos, Lm) :-
    (   get_assoc(Len-Pos, LmA, Lm0) -> Lm = Lm0
    ;   Lm = lm(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    ).

mac_letter_arg(Ch, I) :- char_code(Ch, C), I is C - 64.   % 'A' -> 1

% Initial domain: the full length bucket, intersected with the letter mask of
% every cell already bound (seed pins / seed crossings). Positions with free
% variables constrain nothing yet.
mac_init_domain(DictByLen, LmA, Id-slot(_, _, _, Vars), Id-Dom) :-
    length(Vars, Len),
    (   get_assoc(Len, DictByLen, Ws)
    ->  length(Ws, NW), Full is (1 << NW) - 1
    ;   Full = 0
    ),
    mac_constrain_bound(Vars, 0, Len, LmA, Full, Dom).
mac_constrain_bound([], _, _, _, D, D).
mac_constrain_bound([V|Vs], Pos, Len, LmA, D0, D) :-
    (   var(V)
    ->  D1 = D0
    ;   mac_lm_of(LmA, Len, Pos, Lm),
        mac_letter_arg(V, ChI),
        arg(ChI, Lm, M),
        D1 is D0 /\ M
    ),
    Pos1 is Pos + 1,
    mac_constrain_bound(Vs, Pos1, Len, LmA, D1, D).

% Crossing edges among the SEARCH slots, from shared ground cell numbers:
% EdgeA maps Id -> [e(Pos, OtherId, OtherPos, EdgeId), ...] with EdgeId
% (1-based) shared by both directions of a crossing — the index into the
% conflict-weight compound. Cells shared with SEEDED slots are already bound
% and constrain the initial domains instead (seeded slots are not searched).
mac_edges(IdSlots, EdgeA, NEdges) :-
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

mac_edges_of(EdgeA, Id, Es) :- ( get_assoc(Id, EdgeA, Es0) -> Es = Es0 ; Es = [] ).

% Conflict weights: a flat w/NEdges float compound in a per-search global,
% destructively bumped (nb_setarg) so the learning survives backtracking.
% Re-initialized on every fill_search/5 entry, so consecutive fills in one
% process are independent (INV-2). nb_setval/nb_getval are thread-local.
mac_init_weights(N) :-
    length(Ones, N), maplist(=(1.0), Ones),
    WT =.. [w|Ones],
    nb_setval(cw_mac_weights, WT).
mac_bump_edge(EId) :-
    nb_getval(cw_mac_weights, WT),
    arg(EId, WT, V), V1 is V + 1.0,
    nb_setarg(EId, WT, V1).
mac_age_weights(Factor) :-
    nb_getval(cw_mac_weights, WT),
    functor(WT, w, N),
    forall(between(1, N, I),
           ( arg(I, WT, V), V1 is V * Factor, nb_setarg(I, WT, V1) )).

% Restart loop: growing node cap, weights persist (aged) across attempts.
% The three attempt outcomes are structurally distinct (AC-FILL-14):
%   - success: mac_search/10 succeeds (Recover stays unbound);
%   - cap abandoned: throw(mac_cap) unwinds the attempt, the catch recovers
%     binding Recover = restart -> retry with a x1.5 cap and aged weights;
%   - exhausted: mac_search/10 FAILS -> the whole loop fails. That failure
%     is a completed search over the full tree (propagation and weights
%     reorder branches, never prune legal ones), i.e. genuine infeasibility
%     - it must NEVER be retried (the pre-fix conflation looped forever /
%     float-overflowed the cap on provably-infeasible inputs).
% No attempt-count bound: the inference budget (fill_attempt/8's
% call_with_inference_limit) is the outer stop and maps to not_proven.
mac_restart_loop(Attempt, Cap, Ids, DomA, LenA, BucketA, EdgeA, LmA, Mode,
                 Used, Assign) :-
    mac_attempt_pick(Mode, Attempt, Pick),
    nb_setval(cw_mac_attempt_nodes, 0),
    verbose_report("fill: search attempt ~w (node cap ~w)~n", [Attempt, Cap]),
    catch(mac_search(Ids, DomA, LenA, BucketA, EdgeA, LmA, Pick, Cap,
                     Used, Assign),
          mac_cap, Recover = restart),
    (   Recover == restart
    ->  mac_age_weights(0.99),
        Cap1 is (Cap * 3 + 1) // 2,   % ceiling(x1.5), integer-only arithmetic
        A1 is Attempt + 1,
        mac_restart_loop(A1, Cap1, Ids, DomA, LenA, BucketA, EdgeA, LmA, Mode,
                         Used, Assign)
    ;   true
    ).

% Per-attempt pick policy. Seeded (--seed/--shuffle): every attempt is top3
% on the user's PRNG stream (set_search_seed/1 already seeded it). Default:
% attempt 1 is STRICT GREEDY — the pure candidate order, no draws — so any
% grid the first attempt completes gets the order-maximal fill; attempts
% >= 2 are top3 on a stream reseeded ONCE (at attempt 2) from a pinned
% engine constant (distinct from the load-shuffle constant, so the two
% seams never replay each other's draw sequence). Diversified restarts are
% load-bearing: restarts only pay off when successive attempts explore
% DIFFERENT prefixes, and greedy-only restarts re-tread the same doomed
% subtree (C3: open 17x17 x full ENABLE never converged in 23 greedy
% attempts / ~1M nodes; top3 filled it in a handful). Value-order
% diversity, in turn, comes from the load shuffle (seed_perturb_plain) —
% top3 over a lexicographic clump diversifies nothing (the g17_50k
% failure). Everything is a pure function of the input: byte-identical
% across runs and CLI/WASM (AC-FILL-13 as amended — deterministic, not
% draw-free).
mac_attempt_pick(seeded, _, top3).
mac_attempt_pick(default, 1, greedy) :- !.
mac_attempt_pick(default, 2, top3) :- !,
    prng_seed_state(0xCC9E2D51).         % pinned search stream (murmur3 c1)
mac_attempt_pick(default, _, top3).

mac_search([], _DomA, _LenA, _BucketA, _EdgeA, _LmA, _Pick, _Cap, _Used, []).
mac_search(Unfilled, DomA, LenA, BucketA, EdgeA, LmA, Pick, Cap, Used,
           [Best-Word|Assign]) :-
    Unfilled = [_|_],
    mac_select_dwd(Unfilled, DomA, EdgeA, Best),
    selectchk(Best, Unfilled, Rest),
    get_assoc(Best, DomA, Dom),
    get_assoc(Best, LenA, Len),
    get_assoc(Len, BucketA, Bucket),
    mac_candidate(Pick, Dom, Bucket, Word),
    \+ memberchk(Word, Used),
    nb_getval(cw_mac_attempt_nodes, A0), A1 is A0 + 1,
    nb_setval(cw_mac_attempt_nodes, A1),
    ( A1 > Cap -> throw(mac_cap) ; true ),
    mac_edges_of(EdgeA, Best, Es),
    mac_place(Es, Word, Rest, LenA, LmA, DomA, DomA1, [], Dirty),
    mac_propagate(Dirty, DomA1, LenA, EdgeA, LmA, Rest, DomA2),
    mac_search(Rest, DomA2, LenA, BucketA, EdgeA, LmA, Pick, Cap,
               [Word|Used], Assign).

% LAZY candidate enumeration in the §8.4a order (C3 fix): the eager form
% (materialize the whole domain, then member/2) cost O(popcount x mask-width)
% bignum work per node — on open grids over the full dictionary that burned
% the inference budget before the search could converge (g17_full: 800M gone
% in ~4k nodes). Enumerating set bits on demand makes the common case (an
% early candidate sticks) one lsb/arg lookup; backtracking walks exactly as
% far as the search actually retries. mac_words/3 below is the eager twin,
% retained as the order oracle for the white-box tests.
%
% greedy: strict best-first — no PRNG.
% top3: weighted promotion of one of the first three candidates ([4,2,1] via
% prng_draw); the other two keep their order, then the residue mask continues
% lazily — branch ORDER changes, the branch SET does not, so the attempt
% remains complete (reordered, never dropped).
mac_candidate(greedy, Dom, Bucket, Word) :-
    mac_next_word(Dom, Bucket, Word).
mac_candidate(top3, Dom, Bucket, Word) :-
    mac_firstn(3, Dom, Bucket, First, Residue),
    mac_top3(First, Ordered),
    (   member(Word, Ordered)
    ;   mac_next_word(Residue, Bucket, Word)
    ).

% Nondeterministic ascending-bit walk of a candidate mask: semantically
% member/2 over mac_words/3, materializing nothing. Bits clear with xor, not
% `/\ \Bit` (the int64 `\` quirk — see mac_words).
mac_next_word(M, Bucket, Word) :-
    M =\= 0,
    B is lsb(M),
    (   I is B + 1, arg(I, Bucket, Word)
    ;   M1 is M xor (1 << B),
        mac_next_word(M1, Bucket, Word)
    ).

% First (up to) N set bits' words, plus the mask with those bits cleared.
mac_firstn(0, M, _, [], M) :- !.
mac_firstn(_, 0, _, [], 0) :- !.
mac_firstn(N, M, Bucket, [W|Ws], Residue) :-
    B is lsb(M),
    I is B + 1, arg(I, Bucket, W),
    M1 is M xor (1 << B),
    N1 is N - 1,
    mac_firstn(N1, M1, Bucket, Ws, Residue).

mac_top3([], []).
mac_top3([A], [A]).
mac_top3([A, B], Out) :-
    prng_draw(V), W is V mod 6,
    ( W < 4 -> Out = [A, B] ; Out = [B, A] ).
mac_top3([A, B, C|Rest], [P|Ps]) :-
    prng_draw(V), W is V mod 7,
    (   W < 4 -> P = A, Ps = [B, C|Rest]
    ;   W < 6 -> P = B, Ps = [A, C|Rest]
    ;   P = C, Ps = [A, B|Rest]
    ).

% dom/wdeg selection: minimize popcount(dom)/wdeg, compared cross-multiplied.
% wdeg = Σ weights of crossings to still-unfilled neighbours (+0.001 floor
% for a slot with no unfilled neighbour, e.g. the last slot standing).
mac_select_dwd([Id|Ids], DomA, EdgeA, Best) :-
    nb_getval(cw_mac_weights, WT),
    mac_dwd_key(Id, DomA, EdgeA, [Id|Ids], WT, C, W),
    mac_dwd_walk(Ids, DomA, EdgeA, [Id|Ids], WT, Id, C, W, Best).
mac_dwd_walk([], _, _, _, _, Best, _, _, Best).
mac_dwd_walk([Id|Ids], DomA, EdgeA, Unfilled, WT, B0, C0, W0, Best) :-
    mac_dwd_key(Id, DomA, EdgeA, Unfilled, WT, C, W),
    (   C * W0 < C0 * W
    ->  mac_dwd_walk(Ids, DomA, EdgeA, Unfilled, WT, Id, C, W, Best)
    ;   mac_dwd_walk(Ids, DomA, EdgeA, Unfilled, WT, B0, C0, W0, Best)
    ).
mac_dwd_key(Id, DomA, EdgeA, Unfilled, WT, C, W) :-
    get_assoc(Id, DomA, D), C is popcount(D),
    mac_edges_of(EdgeA, Id, Es),
    mac_dwd_wdeg(Es, Unfilled, WT, 0.001, W).
mac_dwd_wdeg([], _, _, W, W).
mac_dwd_wdeg([e(_, T, _, EId)|Es], Unfilled, WT, W0, W) :-
    (   memberchk(T, Unfilled)
    ->  arg(EId, WT, V), W1 is W0 + V
    ;   W1 = W0
    ),
    mac_dwd_wdeg(Es, Unfilled, WT, W1, W).

% Materialize words from mask bits, ascending bit index (= bucket order =
% score-desc-then-dictionary: the §8.4a value order, preserved structurally).
% NOT on the hot path since the C3 lazy fix (mac_candidate/4 above): retained
% as the order oracle for the white-box tests — it must stay bit-for-bit
% equivalent to a full mac_next_word/3 walk.
% Bits are cleared with xor, NOT `M /\ \(1<<B)`: SWI's `\` misevaluates at
% the int64 boundary (\(1<<63) yields 2^63-1, silently wiping bits >= 63 of
% a bignum AND-chain — found the hard way in the DP-7 spike).
mac_words(0, _, []) :- !.
mac_words(M, Bucket, [W|Ws]) :-
    B is lsb(M),
    I is B + 1, arg(I, Bucket, W),
    M1 is M xor (1 << B),
    mac_words(M1, Bucket, Ws).

% Placement: AND each unfilled crosser's domain with the placed letter's
% mask; collect dirtied ids for the propagation worklist. A wipeout bumps
% the responsible crossing's weight (dom/wdeg's learning signal) and fails.
mac_place([], _Word, _Unfilled, _LenA, _LmA, DomA, DomA, Dirty, Dirty).
mac_place([e(Pos, T, TPos, EId)|Es], Word, Unfilled, LenA, LmA, DomA0, DomA,
          Dirty0, Dirty) :-
    (   memberchk(T, Unfilled)
    ->  nth0(Pos, Word, Ch),
        get_assoc(T, LenA, TLen),
        mac_lm_of(LmA, TLen, TPos, LmT),
        mac_letter_arg(Ch, ChI),
        arg(ChI, LmT, ChMask),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ ChMask,
        (   TDom1 =:= 0
        ->  mac_bump_edge(EId), fail
        ;   true
        ),
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Dirty1 = Dirty0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            Dirty1 = [T|Dirty0]
        )
    ;   DomA1 = DomA0, Dirty1 = Dirty0
    ),
    mac_place(Es, Word, Unfilled, LenA, LmA, DomA1, DomA, Dirty1, Dirty).

% AC-3 to fixpoint: worklist of slot ids whose domain changed; re-filter each
% unfilled crosser's domain to the letters the changed slot still supports at
% the shared cell. A wipeout bumps the crossing's weight and fails the branch.
mac_propagate([], DomA, _LenA, _EdgeA, _LmA, _Unfilled, DomA).
mac_propagate([S|Queue], DomA0, LenA, EdgeA, LmA, Unfilled, DomA) :-
    mac_edges_of(EdgeA, S, Es),
    get_assoc(S, DomA0, SDom),
    get_assoc(S, LenA, SLen),
    mac_revise(Es, SDom, SLen, LenA, LmA, Unfilled, DomA0, DomA1, Queue, Queue1),
    mac_propagate(Queue1, DomA1, LenA, EdgeA, LmA, Unfilled, DomA).

mac_revise([], _SDom, _SLen, _LenA, _LmA, _Unfilled, DomA, DomA, Q, Q).
mac_revise([e(Pos, T, TPos, EId)|Es], SDom, SLen, LenA, LmA, Unfilled,
           DomA0, DomA, Q0, Q) :-
    (   memberchk(T, Unfilled)
    ->  get_assoc(T, LenA, TLen),
        mac_lm_of(LmA, SLen, Pos, LmS),
        mac_lm_of(LmA, TLen, TPos, LmT),
        mac_support(1, LmS, SDom, LmT, 0, Supp),
        get_assoc(T, DomA0, TDom),
        TDom1 is TDom /\ Supp,
        (   TDom1 =:= 0
        ->  mac_bump_edge(EId), fail
        ;   true
        ),
        (   TDom1 =:= TDom
        ->  DomA1 = DomA0, Q1 = Q0
        ;   put_assoc(T, DomA0, TDom1, DomA1),
            ( memberchk(T, Q0) -> Q1 = Q0 ; Q1 = [T|Q0] )
        )
    ;   DomA1 = DomA0, Q1 = Q0
    ),
    mac_revise(Es, SDom, SLen, LenA, LmA, Unfilled, DomA1, DomA, Q1, Q).

% Letters S still supports at the shared cell, projected into T's mask space:
% a 26-way sweep of two lm/26 compounds (arg/3 + one AND test per letter).
mac_support(27, _, _, _, Supp, Supp) :- !.
mac_support(I, LmS, SDom, LmT, Acc, Supp) :-
    arg(I, LmS, SM),
    (   SDom /\ SM =\= 0
    ->  arg(I, LmT, TM),
        Acc1 is Acc \/ TM
    ;   Acc1 = Acc
    ),
    I1 is I + 1,
    mac_support(I1, LmS, SDom, LmT, Acc1, Supp).

% --- the RETAINED pre-§8.4c MRV search (reference; not on the hot path) --------
fill_search_mrv(Slots, DictByLen, Index, Masks, Used) :-
    maplist(slot_with_count(DictByLen, Index, Masks), Slots, Counted),
    % §8.4b dispatch, read ONCE per search (not per node).
    (   current_search_seed(_)
    ->  fill_search_inc_seeded(Counted, DictByLen, Index, Masks, Used)
    ;   fill_search_inc(Counted, DictByLen, Index, Masks, Used)
    ).

% Pair each slot with its current candidate count, keeping the ORIGINAL slot
% term (shared cell variables intact - no reconstruction, no copy). This is the
% one full per-slot recount, paid once at the root; thereafter it is incremental.
slot_with_count(DictByLen, Index, Masks, Slot, cnt(Count, Slot)) :-
    Slot = slot(_, _, _, Vars),
    candidate_count(Masks, Vars, DictByLen, Index, Count).

% Counted is a list of cnt(Count, Slot). Select the min-count slot, place a word,
% recount only the crossings the placement disturbed, recurse on the rest.
% The []/[_|_] head pair is SWI's two-clause special case (jitindex §2.17):
% clause selection is deterministic with no choicepoint and no cut, and the
% recursive clause never pays a failed []-head unification per node.
fill_search_inc([], _DictByLen, _Index, _Masks, _Used).
fill_search_inc([C0|Cs], DictByLen, Index, Masks, Used) :-
    select_min_count([C0|Cs], cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    newly_bound_cells(BestCells, BestVars, NewCells),   % free cells, PRE-placement
    candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    BestVars = Word,                      % unify into the shared cell variables
    recount_crossing(Rest, NewCells, DictByLen, Index, Masks, Rest1),
    fill_search_inc(Rest1, DictByLen, Index, Masks, [Word|Used]).

% §8.4b (AC-FILL-10) seeded twin of fill_search_inc - the F-H2 twin idiom
% (fill_attempt_masked precedent): kept separate so the deterministic loop
% above stays byte-identical for the ratchet. Two differences only: slot
% selection picks among the equal-min-count TIES via the core PRNG (the
% deterministic loop's tie-break is (Start, Dir) order - exactly the tie the
% DP-6 contract lets the seed vary), and candidate order arrives already
% seed-perturbed from the load seams. Reorders branches, never prunes or adds
% them: completeness and the budget contract are unchanged.
fill_search_inc_seeded([], _DictByLen, _Index, _Masks, _Used).
fill_search_inc_seeded([C0|Cs], DictByLen, Index, Masks, Used) :-
    select_min_count_seeded([C0|Cs], cnt(_, Best), Rest),
    Best = slot(_, _, BestCells, BestVars),
    newly_bound_cells(BestCells, BestVars, NewCells),
    candidates(BestVars, DictByLen, Index, Cands),
    member(Word, Cands),
    \+ memberchk(Word, Used),
    BestVars = Word,
    recount_crossing(Rest, NewCells, DictByLen, Index, Masks, Rest1),
    fill_search_inc_seeded(Rest1, DictByLen, Index, Masks, [Word|Used]).

% The min COUNT, then a PRNG pick among every cnt/2 sharing it (a one-element
% seeded_permutation head-pick: deterministic per seed, advances the shared
% stream like arrange's bucket shuffles). Rest preserves order + term-sharing,
% same as the deterministic remove_slot path.
select_min_count_seeded(Counted, Min, Rest) :-
    Counted = [H|T],
    min_count_walk(T, H, cnt(MinC, _)),
    include(tie_count(MinC), Counted, Ties),
    seeded_permutation(Ties, [Min|_]),
    Min = cnt(_, slot(MStart, MDir, _, _)),
    remove_slot(MStart, MDir, Counted, Rest).

tie_count(C, cnt(C0, _)) :- C0 == C.

% Winner = the cnt/2 with the minimum c(Count, Start, Dir) in standard order
% (Count, then Start, then Dir; all ground). This reproduces select_mrv's
% sort(0, @=<, ...)+head+once(select) exactly. Rest is the other cnt/2s with
% their order + term-sharing preserved (mirrors once(select(Best, Slots, Rest))).
select_min_count([H|T], Min, Rest) :-
    min_count_walk(T, H, Min),
    Min = cnt(_, slot(MStart, MDir, _, _)),
    remove_slot(MStart, MDir, [H|T], Rest).

min_count_walk([], Best, Best).
min_count_walk([X|Xs], Best0, Best) :-
    ( count_le(Best0, X) -> Best1 = Best0 ; Best1 = X ),
    min_count_walk(Xs, Best1, Best).

% Best0 stays iff its (Count, Start, Dir) key is @=< X's - the same key order
% select_mrv's sort/4 imposes on the c/3 terms.
count_le(cnt(C0, slot(S0, D0, _, _)), cnt(C1, slot(S1, D1, _, _))) :-
    c(C0, S0, D0) @=< c(C1, S1, D1).

% Drop the winning slot by its ground (Start, Dir) key, compared with == so no
% variable is bound. Preserves the order and term-sharing of the remaining cnt/2s.
remove_slot(MStart, MDir, [cnt(_, slot(S, D, _, _))|T], Rest) :-
    S == MStart, D == MDir, !, Rest = T.
remove_slot(MStart, MDir, [X|T], [X|Rest]) :-
    remove_slot(MStart, MDir, T, Rest).

% The placed slot's cells that were still FREE before placement - exactly the
% cells this placement newly binds. Its already-bound cells are unchanged, so
% crossings through them cannot change count. Cells and Vars are positionally
% aligned (slot_vars/3). NewCells is a plain list of ground cell NUMBERS: the
% crossing test is by cell identity, so no cell variable is ever copied.
newly_bound_cells([], [], []).
newly_bound_cells([Cell|Cs], [Var|Vs], New) :-
    ( var(Var) -> New = [Cell|New1] ; New = New1 ),
    newly_bound_cells(Cs, Vs, New1).

% Recount EXACTLY the carried slots that cross a newly-bound cell; carry every
% other count unchanged. cnt(Count, Slot) keeps the ORIGINAL Slot term (shared
% cell variables survive - no findall/copy). A crossing is a shared cell NUMBER
% (ground), so detection needs no variable comparison. Both Cells and NewCells
% are strictly ascending cell numbers by construction - Cells is a contiguous
% grid_run row/column run (stockgrid: numlist / col_cells + grab_run), and
% NewCells is an order-preserving filter of the placed slot's Cells
% (newly_bound_cells/3) - i.e. both are ordsets, so the shared-cell test is
% ord_intersect/2: one O(|Cells|+|NewCells|) merge walk, choicepoint-free,
% instead of a memberchk scan of NewCells per cell.
recount_crossing([], _NewCells, _DictByLen, _Index, _Masks, []).
recount_crossing([cnt(C0, Slot)|T], NewCells, DictByLen, Index, Masks, [cnt(C1, Slot)|T1]) :-
    Slot = slot(_, _, Cells, Vars),
    ( ord_intersect(Cells, NewCells)
    ->  candidate_count(Masks, Vars, DictByLen, Index, C1)
    ;   C1 = C0 ),
    recount_crossing(T, NewCells, DictByLen, Index, Masks, T1).

% select_mrv/6 is RETAINED as the reference selector: the white-box tests
% (tests/fill.plt R6/P13) call it directly, and it documents the exact selection
% rule fill_search_inc reproduces incrementally. It is no longer on the engine's
% hot path (fill_search_inc drives the search).
% Pick the slot with the fewest current candidates (>=0); deterministic
% tie-break: lowest start cell, then direction. Recovering the slot by BOTH
% start AND direction is load-bearing when a cell begins an across and a down
% slot (e.g. cell 1): the slot whose count was the minimum must be the one
% expanded, not whichever shares the start and appears first in the list (R6).
select_mrv(Slots, DictByLen, Index, Best, Rest, BestCands) :-
    maplist(slot_candidate_count(DictByLen, Index), Slots, Counted),
    sort(0, @=<, Counted, [c(_, BestStart, BestDir)|_]),
    Best = slot(BestStart, BestDir, _, _),
    once(select(Best, Slots, Rest)),   % Start+Dir is unique, so prune the CP (P13)
    Best = slot(_, _, _, Vars),
    candidates(Vars, DictByLen, Index, BestCands).

% c(Count, Start, Dir) - sorted by Count, then Start, then Dir gives
% most-constrained-first deterministically AND uniquely identifies the slot
% (start+dir), so select_mrv recovers exactly the slot the count was computed for.
slot_candidate_count(DictByLen, Index, slot(Start, Dir, _, Vars), c(Count, Start, Dir)) :-
    candidate_count(Vars, DictByLen, Index, Count).


% --- emit the filled layout --------------------------------------------------
slots_to_layout(Slots, Numbered, InputWords) :-
    maplist(slot_to_word, Slots, Placed),
    % No once/1: assign_clue_numbers/2 is deterministic at source (X6.B4;
    % probe re-verified on this call shape) - the wrap was defensive (C23).
    assign_clue_numbers(Placed, Numbered),
    % maplist, not findall+member (C35): a deterministic 1:1 projection in the
    % same order, with no findall copy per entry.
    maplist(pw_answer_entry, Placed, InputWords).

pw_answer_entry(PW, [A]) :- pw_answer(PW, A).

slot_to_word(slot(Start, Dir, Cells, Vars),
             pw(A, Vars, Cells, Dir, Len, Start, _End, _Num)) :-
    length(Vars, Len), atom_chars(A, Vars).


% --- entry point -------------------------------------------------------------
% seed_used(+AllSlots, +SearchSlots, -Used0): the search's INITIAL no-duplicate
% set = the seed answers. The seeded slots are exactly AllSlots minus
% SearchSlots (fill_prepare's exclude/3 keeps the same slot terms, so
% subtract/3's memberchk recovers them; Start+Dir are ground and unique, so no
% free cell variable can be bound by the scan). A seeded slot's Vars are fully
% bound by apply_seed and ARE the placed word in the exact shape Used holds (a
% list of char atoms), so fill_search_inc's `\+ memberchk(Word, Used)` now
% dedups searched slots against the seed pins too (regression: a seed answer
% the dictionary could also place was re-placed in a searched slot and the
% duplicate blew up the emit; see docs/plans/fill-seed-pin-crash-fix.md).
% Do NOT collapse Vars to an atom (atom_chars): Used holds char LISTS and an
% atom never matches, silently disabling the dedup. With no seeds,
% SearchSlots == AllSlots: the guard below yields Used0 == [] without the
% subtract scan, so the benched fill_attempt/8 window stays neutral on unseeded
% rungs (the seed_used cost lands only on seeded rungs). subtract(X, X, R) would
% also give R == [], so the fast path is byte-identical, only cheaper.
seed_used(AllSlots, SearchSlots, Used0) :-
    (   AllSlots == SearchSlots
    ->  Used0 = []
    ;   subtract(AllSlots, SearchSlots, SeededSlots),
        maplist(slot_word, SeededSlots, Used0)
    ).

slot_word(slot(_, _, _, Vars), Vars).

% Outcome: filled | infeasible | not_proven. seeds/dict applied before search;
% a seed clash / no-slot / duplicate-answer seed throws (reported before
% searching). SearchSlots are the slots the engine fills (all slots minus seed
% pins); AllSlots are emitted (so seed pins appear in the layout too).
fill_attempt(SearchSlots, AllSlots, DictByLen, Index, Outcome, Numbered, InputWords) :-
    fill_budget(B),
    fill_attempt(SearchSlots, AllSlots, DictByLen, Index, B, Outcome, Numbered, InputWords).

%!  fill_attempt(+SearchSlots:list, +AllSlots:list, +DictByLen:assoc,
%!               +Index:assoc, +Budget:integer,
%!               -Outcome:oneof([filled,infeasible,not_proven]),
%!               -Numbered:list, -InputWords:list) is det.
%
%   Budget-explicit form: Budget is the inference budget. The default /7
%   above reads fill_budget/1 and stays internal (the default-budget
%   delegate, mirroring arrange_best_layout/5); passing a tiny Budget drives
%   the AC-FILL-1 "not proven within budget" path (tests/fill.plt). Always
%   succeeds with exactly one Outcome - infeasible is an Outcome, never a
%   failure. Benchmark seam: the fill bench's GATED `search` bucket
%   (search_inf, the metric of record) times exactly this goal with a
%   pre-loaded dict and fresh slots rebuilt outside the timed window.
fill_attempt(SearchSlots, AllSlots, DictByLen, Index, Budget, Outcome, Numbered, InputWords) :-
    % No catch/3: see construct_one/7 in arrange.pl - call_with_inference_limit/3
    % binds Limit = inference_limit_exceeded on the budget path and only re-throws
    % genuine errors. Infeasibility is a search FAILURE (R == exhausted), so a
    % thrown error is a real bug and must surface, not be masked as infeasible.
    % Ordset entry (Masks == none): fill_search/5 with `none` does byte-identical
    % work to the pre-F-H2 fill_search/4, so this bench/test seam is inference-
    % identical. The product artifact path uses fill_attempt_masked/9 below.
    % Used0 is computed BEFORE call_with_inference_limit/3 so it is never
    % charged against the search Budget.
    seed_used(AllSlots, SearchSlots, Used0),
    call_with_inference_limit(
        ( fill_search(SearchSlots, DictByLen, Index, none, Used0) -> R = ok ; R = exhausted ),
        Budget, Limit),
    (   Limit == inference_limit_exceeded
    ->  Outcome = not_proven, Numbered = [], InputWords = []
    ;   R == ok
    ->  Outcome = filled, slots_to_layout(AllSlots, Numbered, InputWords)
    ;   Outcome = infeasible, Numbered = [], InputWords = []
    ).

% Masks-carrying twin of fill_attempt/8 (product path only - NOT the gated bench
% seam). Identical control flow; the only difference is the threaded Masks (none
% for the raw CLI, masks(_) for a v2 artifact). Kept separate so fill_attempt/8
% above stays byte-identical for the ratchet.
fill_attempt_masked(SearchSlots, AllSlots, DictByLen, Index, Budget, Masks,
                    Outcome, Numbered, InputWords) :-
    seed_used(AllSlots, SearchSlots, Used0),
    call_with_inference_limit(
        ( fill_search(SearchSlots, DictByLen, Index, Masks, Used0) -> R = ok ; R = exhausted ),
        Budget, Limit),
    (   Limit == inference_limit_exceeded
    ->  Outcome = not_proven, Numbered = [], InputWords = []
    ;   R == ok
    ->  Outcome = filled, slots_to_layout(AllSlots, Numbered, InputWords)
    ;   Outcome = infeasible, Numbered = [], InputWords = []
    ).

%!  fill_solve(+GridFile:atom, +SeedFileOrNone, +DictFile:atom,
%!             +SizeMode:oneof([fixed,max])) is semidet.
%
%   The `fill` CLI seam: derive the slots of the stock grid GridFile, pin the
%   seeds of SeedFileOrNone (a §6.6 fragment file, or the atom `none`), fill
%   every remaining slot from the text dictionary DictFile, and emit the
%   filled canonical layout on stdout framed by SizeMode. FAILS (no stdout)
%   on any non-filled outcome - not_proven (budget) or infeasible - after
%   reporting the unfillable slot(s) on stderr (INV-3, never silent). Throws
%   on a malformed grid/seed file, an unmatchable seed, or two seeds pinning
%   the same answer (fill_seed_duplicate). Default-options delegate of
%   fill_solve/5.
fill_solve(GridFile, SeedFileOrNone, DictFile, SizeMode) :-
    fill_solve(GridFile, SeedFileOrNone, DictFile, SizeMode, []).

%!  fill_solve(+GridFile:atom, +SeedFileOrNone, +DictFile:atom,
%!             +SizeMode:oneof([fixed,max]), +Options:list) is semidet.
%
%   §8.4a/§8.4b options-carrying form of fill_solve/4 (same contract). Options:
%     * min_score(N): the hard candidate prune, in the dictionary's native
%       units (default 1, DP-5 - excludes only the score-0 blocklist floor).
%     * report_json(FileOrNone): write the fill-quality report as a one-line
%       sorted-key JSON object to File on success (default `none`). The
%       stderr human summary rides `--verbose` (§5.1's quiet-success stderr
%       contract); the canonical stdout layout is untouched (AC-FILL-7).
%     * budget(N): the §8.4b inference budget (AC-FILL-9). -1 (the default
%       sentinel) means the engine default, fill_budget/1 - which stays the
%       single source of truth. The budget never changes WHICH fill is
%       found, only fill-vs-not-proven; it is an escape hatch, NOT a
%       completion fix (DP-6, measured).
fill_solve(GridFile, SeedFileOrNone, DictFile, SizeMode, Options) :-
    option(min_score(MinScore), Options, 1),
    option(report_json(ReportFile), Options, none),
    fill_effective_budget(Options, Budget),
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    fill_capacity_guard(
        load_dict(DictFile, [min_score(MinScore)], DictByLen, Index, Scores),
        load(DictFile)),
    fill_capacity_guard(
        fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, none,
                            SizeMode, report_ctx(Scores, ReportFile), Budget),
        search(Size)).

% Resolve the effective inference budget from a budget(N) option. The -1
% sentinel (or no option) selects the engine default fill_budget/1; anything
% else must be a positive integer (the CLI validates too, but white-box
% callers get the same clean error).
fill_effective_budget(Options, Budget) :-
    option(budget(Budget0), Options, -1),
    (   Budget0 == -1
    ->  fill_budget(Budget)
    ;   must_be(positive_integer, Budget0),
        Budget = Budget0
    ).

%!  fill_solve_index(+GridFile:atom, +SeedFileOrNone, +IndexFile:atom,
%!                   +DictFileOrNone, +SizeMode:oneof([fixed,max])) is semidet.
%
%   Artifact-consuming twin of fill_solve/4 (F-L2): load a prebuilt,
%   verified index artifact (fill_save_index/2,3) instead of parsing a text
%   dictionary, then fill IDENTICALLY - the filled layout is byte-for-byte
%   identical to the raw path (proven by the identity oracle in both modes).
%   DictFileOrNone: an atom path checks the artifact's embedded SHA-256;
%   `none` skips that check (version + SWI are always checked); any mismatch
%   throws a clear rebuild error. Same failure contract as fill_solve/4.

% Grid/seed derivation, search, and emit are the shared body below.
% A v2 artifact carries precomputed bitset Masks; F-H2's bignum counting kernel
% runs iff Masks == masks(_). A `none` here (no masks in the artifact) transparently
% falls back to the ordset kernel, so this entry is correct for any artifact shape
% the loader accepts.
% The artifact path carries no scores (v1: --min-score/--report-json are
% text-dict-only, rejected upstream at the CLI), so its report context is
% report_ctx(none, none) - no quality report on this path.
fill_solve_index(GridFile, SeedFileOrNone, IndexFile, DictFileOrNone, SizeMode) :-
    fill_solve_index(GridFile, SeedFileOrNone, IndexFile, DictFileOrNone,
                     SizeMode, []).

%!  fill_solve_index(+GridFile:atom, +SeedFileOrNone, +IndexFile:atom,
%!                   +DictFileOrNone, +SizeMode:oneof([fixed,max]),
%!                   +Options:list) is semidet.
%
%   §8.4b options-carrying form of fill_solve_index/5 (same contract).
%   Options: budget(N) only - `--budget` is a pure search parameter and
%   composes with the artifact path (AC-FILL-9/AC-FILL-11), unlike the
%   score/seed flags (artifacts carry no scores and pin bucket order).
fill_solve_index(GridFile, SeedFileOrNone, IndexFile, DictFileOrNone,
                 SizeMode, Options) :-
    fill_effective_budget(Options, Budget),
    fill_capacity_guard(
        fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, Masks),
        load(IndexFile)),
    fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots),
    fill_capacity_guard(
        fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, Masks,
                            SizeMode, report_ctx(none, none), Budget),
        search(Size)).

% Grid + seed derivation, shared by both entry points (identical to the raw
% path's front matter). Slots are ALL slots (emitted); SearchSlots exclude seed
% pins (what the engine fills).
fill_prepare(GridFile, SeedFileOrNone, Size, Slots, SearchSlots) :-
    fill_grid(GridFile, Size, Slots, _CellVar),
    ( SeedFileOrNone == none
    ->  SeededKeys = []
    ;   apply_seeds(SeedFileOrNone, Size, Slots, SeededKeys) ),
    exclude(seeded_slot(SeededKeys), Slots, SearchSlots).

% Search + emit, shared by both entry points. Fails (no stdout) on any
% non-filled outcome (INV-3): report the unfillable slot(s), never silent.
% Masks threads F-H2's counting context: `none` runs the ordset kernel (raw CLI
% path - byte-identical to pre-F-H2, exercised by the identity oracle), masks(_)
% runs the bignum kernel (v2 artifact). Budget is the §8.4b effective budget
% (fill_effective_budget/2: the fill_budget/1 default unless a budget(N)
% option overrode it) - the raw branch feeds it to the budget-explicit
% fill_attempt/8 (the gated bench seam, arity untouched); only the artifact
% branch takes the masked core.
fill_place_and_emit(Size, Slots, SearchSlots, DictByLen, Index, Masks, SizeMode,
                    ReportCtx, Budget) :-
    (   Masks == none
    ->  fill_attempt(SearchSlots, Slots, DictByLen, Index, Budget, Outcome,
                     Numbered, InputWords)
    ;   fill_attempt_masked(SearchSlots, Slots, DictByLen, Index, Budget, Masks,
                            Outcome, Numbered, InputWords)
    ),
    (   Outcome == filled
    ->  emit_fill(Numbered, InputWords, Size, SizeMode),
        length(Numbered, NS),
        verbose_report("fill: grid ~wx~w, filled ~w slots~n", [Size, Size, NS]),
        fill_quality_report(ReportCtx, Numbered, DictByLen)
    ;   fill_report_failure(Outcome, SearchSlots, DictByLen, Index, Size),
        fail
    ).

% --- the §8.4a fill-quality sidecar report (AC-FILL-7, INV-3) -----------------
% A SIDECAR: the canonical stdout layout is byte-untouched. The human summary
% rides the --verbose stderr channel (§5.1: quiet on clean success; the report
% is informational, not a compromise); report_json(File) is the unconditional
% machine channel when requested. The artifact path carries no scores
% (Scores == none) and emits no report. Written only on a FILLED outcome, so
% --report-json never leaves a partial/failed report file (the §5.2 --out
% discipline). The threshold is the documented clean-floor convention (50, in
% the dict's native units - DP-5's recommended `--min-score 50`), matching
% benchmarks/fill_quality/score_fill.py's below-50 column by construction.
fill_quality_threshold(50).

fill_quality_report(report_ctx(Scores, ReportFile), Numbered, DictByLen) :-
    (   Scores == none
    ->  true
    ;   fill_quality_stats(Numbered, Scores, DictByLen, stats(N, Mean, Min, Below)),
        fill_quality_threshold(T),
        verbose_report("fill: quality n=~w mean=~w min=~w below~w=~w~n",
                       [N, Mean, Min, T, Below]),
        (   ReportFile == none
        ->  true
        ;   write_quality_json(ReportFile, N, Mean, Min, Below, T)
        )
    ).

% n / mean (1 d.p.) / min / below-threshold over every placed answer,
% mirroring score_fill.py: an answer absent from the dictionary (a non-dict
% seed pin) scores 0. Answers are normalized through the same
% word_letters/3 footprint the seed path uses, so multi-word/hyphenated
% seeds score by their placed letters, exactly as the benchmark scores them.
fill_quality_stats(Numbered, Scores, DictByLen, stats(N, Mean, Min, Below)) :-
    maplist(placed_word_score(Scores, DictByLen), Numbered, Vals),
    length(Vals, N),
    sum_list(Vals, Sum),
    Mean is round(Sum * 10 / N) / 10.0,
    min_list(Vals, Min),
    fill_quality_threshold(T),
    foldl(count_below(T), Vals, 0, Below).

count_below(T, V, A0, A) :- ( V < T -> A is A0 + 1 ; A = A0 ).

placed_word_score(Scores, DictByLen, PW, Score) :-
    pw_answer(PW, A),
    word_letters([A], Letters, _),
    (   Scores = scores(SA)
    ->  ( get_assoc(Letters, SA, S) -> Score = S ; Score = 0 )
    ;   % uniform: every in-dict word scores 1 (DP-5); absent -> 0
        length(Letters, L),
        (   get_assoc(L, DictByLen, Bucket),
            memberchk(Letters, Bucket)
        ->  Score = 1
        ;   Score = 0
        )
    ).

% One sorted-key JSON object, hand-formatted: four integers and a 1-d.p.
% float need no JSON library, and the shape is pinned by golden + plunit.
write_quality_json(File, N, Mean, Min, Below, T) :-
    setup_call_cleanup(
        open(File, write, S, [encoding(utf8)]),
        format(S, '{"belowThreshold":~w,"mean":~w,"min":~w,"n":~w,"threshold":~w}~n',
               [Below, Mean, Min, N, T]),
        close(S)).

% The fill emit boundary. First re-assert the unique-answers invariant that
% emit_json/3's metadata join requires (nothing upstream on fill's pipeline
% runs check_unique_answers/1, so without this it never ran here): any
% duplicate that still reaches emit reports as the
% clean hooked duplicate_answer error (exit 1), never the raw
% domain_error(unique_key_pairs) from core's answer_meta_assoc/2. Defense in
% depth only - seed_used/3 + check_unique_seed_answers/1 are the real fixes,
% and this guard alone could NOT be (it would fail solvable grids the search
% must instead fill with a non-duplicate word). The check sits HERE, outside
% fill_attempt/8, so the gated bench window (call_time over fill_attempt/8)
% is untouched on every rung.
emit_fill(Numbered, InputWords, Size, SizeMode) :-
    check_unique_answers(InputWords),
    emit_fill_mode(SizeMode, Numbered, InputWords, Size).

% fixed: the exact Size x Size canonical layout (blocks as null). (`max` would
% crop, but a stock grid is already its own frame, so fixed is the norm.)
emit_fill_mode(fixed, Numbered, InputWords, Size) :-
    emit_json(Numbered, InputWords, Size).
emit_fill_mode(max, Numbered, InputWords, Size) :-
    emit_arrange(Numbered, InputWords, Size, max).

% --- DP-10: clean capacity failure (AC-FILL-15) -------------------------------
% SWI renders resource_error(stack) as a multi-line raw dump (stack stats + a
% backtrace that includes kilobyte bignum mask spew from the §8.4c kernels) -
% it reads as an engine bug and names no remedy. Both entry points guard
% their two phases with this catch, so a capacity overflow becomes the
% ordinary report-and-fail failure shape (the budget-exhaustion style above);
% every other exception rethrows untouched. The capacity ENVELOPE itself is
% unchanged and documented (benchmarks/fill_quality/README.md, DP-9 findings).
% The relaunch hint is verified to work: the shebang launcher runs normally
% under an explicit `swipl --stack-limit=<size> ./crosswordsmith ...`.
fill_capacity_guard(Goal, Phase) :-
    catch(Goal, error(resource_error(stack), _), fill_capacity_fail(Phase)).

fill_capacity_fail(load(File)) :-
    fill_stack_limit_desc(Lim),
    format(user_error,
           "fill: capacity exceeded loading/indexing ~w (SWI stack limit ~w)~n",
           [File, Lim]),
    fill_capacity_hint('shrink the dictionary (e.g. --min-score 50)'),
    fail.
fill_capacity_fail(search(Size)) :-
    fill_stack_limit_desc(Lim),
    format(user_error,
           "fill: capacity exceeded during search on ~wx~w grid (SWI stack limit ~w)~n",
           [Size, Size, Lim]),
    fill_capacity_hint('raise --min-score or use a smaller dictionary'),
    fail.

fill_capacity_hint(Remedy) :-
    format(user_error,
           "fill: hint: ~w, or relaunch with a larger stack: swipl --stack-limit=4g ./crosswordsmith ...~n",
           [Remedy]).

fill_stack_limit_desc(Desc) :-
    current_prolog_flag(stack_limit, Limit),
    (   Limit >= 1073741824
    ->  Gb is Limit / 1073741824.0,
        format(atom(Desc), '~1fGb', [Gb])
    ;   Mb is Limit / 1048576.0,
        format(atom(Desc), '~0fMb', [Mb])
    ).

fill_report_failure(not_proven, _Slots, _D, _I, Size) :-
    format(user_error,
           "fill: not proven within budget on ~wx~w grid (search did not complete)~n",
           [Size, Size]).
fill_report_failure(infeasible, Slots, DictByLen, Index, Size) :-
    ( empty_slots(Slots, DictByLen, Index, Bad), Bad = [_|_]
    ->  format(user_error,
               "fill: no dictionary word fits slot(s) ~w on the ~wx~w grid~n",
               [Bad, Size, Size])
    ;   format(user_error,
               "fill: no complete fill exists for this grid + dictionary (+ seeds)~n",
               [])
    ).

% Slots (by start cell) that already have ZERO candidate words - genuinely
% unfillable, reportable up front. "Count is zero" is asked of the purpose-
% built counter (C36), not by materializing candidates/4's word list and
% matching [].
empty_slots(Slots, DictByLen, Index, Bad) :-
    findall(Start,
            ( member(slot(Start, _, _, Vars), Slots),
              candidate_count(Vars, DictByLen, Index, 0) ),
            Bad).


% --- precomputed index artifact (F-L2) ---------------------------------------
% The dictionary parse + index build is a PURE FUNCTION of the frozen dict file
% and dominates end-to-end latency (P-F1: dict_load is 58-84% of CLI wall on
% 10/11 rungs; F-L1 cut the parse, leaving the build_index keysort + GC as the
% inference-blind residue). This serializes the EXACT runtime structures (the
% DictByLen length buckets + the assoc Index) once, so an interactive fill loads
% them back instead of recomputing them every invocation. The raw text-dict path
% (load_dict/build_index) is untouched: this builder CALLS load_dict, and the
% loader reconstructs terms that are ==-identical to a fresh build.
%
% ARTIFACT TERM (versioned + extensible):
%   fill_index(Version, Meta, DictByLen, Index)
%     Version - integer artifact-SCHEMA version (fill_index_format_version/1). A
%               schema change bumps it; the loader refuses an unknown version.
%     Meta    - a keyed list [Key(Value), ...] carrying integrity + provenance
%               AND, since schema v2 (F-H2), OPTIONALLY the precomputed bitset
%               masks under a masks(...) key (the realized extension point the
%               F-L2 design reserved; the Version bump 1->2 makes it a clean
%               break - v1 artifacts are refused, never misread).
%       dict_sha256(Hex) - SHA-256 of the source dict file's bytes (staleness)
%       swi_version(V)   - the SWI-Prolog that built it (binary-format guard)
%       words(N)         - dictionary word count (informational)
%       source(Path)     - the dict path as given at build (informational)
%       built_epoch(E)   - build time, integer seconds (informational)
%       masks(MaskAssoc) - [v2, OPTIONAL - `--save-index --masks` only] assoc
%                          k(Len,Pos,Char) -> bignum; bit i set iff bucket index
%                          i is in Index's ordset for that key. The F-H2 counting
%                          kernel AND-folds + popcounts these instead of
%                          ord_intersection'ing the ordsets (wall win, count
%                          identical). Built ONLY here (amortized), never at
%                          load; absent -> the loader hands back Masks = none and
%                          counting stays on the ordset kernel (no size/load tax).
%
% ON-DISK FORMAT: fast_write/fast_read (library(fastrw)) - MEASURED the fastest
% candidate (F-L2 results doc): ~20x the post-F-L1 warm raw load at 172k, smaller
% on disk than the .qlf of the equivalent fact, ~0 inferences. fast_read's binary
% format is SWI-version-bound (documented), so the artifact embeds swi_version
% and the loader REFUSES a mismatch (rebuild explicitly, never silently).

% Schema version 2 (F-H2): the Meta list MAY carry a masks(MaskAssoc) key
% (precomputed bignum bitsets for the counting kernel). Masks are OPTIONAL within
% v2: the default build omits them (F-H2 follow-up finding - masks inflate the
% artifact ~32% and tax EVERY load ~+30% fast_read wall, a net end-to-end LOSS on
% load-dominated fills), and `--save-index --masks` opts in for search-heavy /
% WASM deployments where the 23-29% fill-phase win outweighs the load tax. The
% loader accepts both shapes (absent key -> ordset kernel). The version bump is a
% hard break - a v1 artifact (or any Version != 2) is refused by the loader with
% the existing rebuild message. Masks are built ONLY here (the amortized
% --save-index step - the gate probe's binding condition), never at load or fill.
%
% Schema version 3 (Unicode hardening): the TERM SHAPE is identical to v2; the
% bump is SEMANTIC. The artifact bakes in normalize_word's output, and the
% dict_sha256 staleness check hashes the SOURCE file, so it cannot see a
% normalization-policy change - v2 artifacts built from a non-ASCII dictionary
% embed the old locale-dependent normalization. The bump routes them into the
% existing fill_index_version refusal (explicit rebuild, never silent misread).
% For pure-ASCII dictionaries a rebuilt artifact is ==-identical to its v2
% predecessor (the hardening changes nothing for ASCII input).
fill_index_format_version(3).

%!  fill_save_index(+DictFile:atom, +OutFile:atom) is semidet.
%!  fill_save_index(+DictFile:atom, +OutFile:atom, +Options:list) is semidet.
%
%   BUILD the precomputed index artifact: load_dict the frozen DictFile and
%   serialize the exact runtime structures + integrity/provenance meta to
%   OutFile (fast_write binary). This is the one-off cost (~ current load +
%   a fast_write); the CLI's `fill --save-index FILE` seam calls it. /2 = the
%   default build (no masks). /3: Options is a keyed list; masks(true)
%   additionally derives the parallel bitset masks FROM the ordset Index (so
%   popcount agreement is by construction) and embeds them as the masks(...)
%   Meta key. Any other Options content is ignored (same permissive
%   convention as Meta itself). Throws on an unreadable dictionary file;
%   fails (after the AC-FILL-15 capacity report, DP-10) when loading/indexing
%   the dictionary overflows the stack limit - det on every in-envelope input.
fill_save_index(DictFile, OutFile) :-
    fill_save_index(DictFile, OutFile, []).

fill_save_index(DictFile, OutFile, Options) :-
    fill_capacity_guard(load_dict(DictFile, DictByLen, Index), load(DictFile)),
    dict_word_count(DictByLen, NWords),
    fill_file_sha256(DictFile, Sha),
    fill_swi_version(Swi),
    get_time(TF), Epoch is round(TF),
    fill_index_format_version(Version),
    Meta0 = [ dict_sha256(Sha), swi_version(Swi), words(NWords),
              source(DictFile), built_epoch(Epoch) ],
    (   option(masks(true), Options)
    ->  build_masks(Index, MaskAssoc),
        append(Meta0, [masks(MaskAssoc)], Meta)
    ;   Meta = Meta0
    ),
    Artifact = fill_index(Version, Meta, DictByLen, Index),
    setup_call_cleanup(open(OutFile, write, S, [type(binary)]),
                       fast_write(S, Artifact),
                       close(S)).

% Derive the mask assoc from the ordset Index: for each k(Len,Pos,Char) key whose
% value is an ordset of bucket indices, the mask is the bignum with bit i set for
% every index i in that ordset. Because masks and buckets share index_set's
% ordering (the masks ARE the ordsets, re-encoded), popcount(mask chain) equals
% |ord_intersection(set chain)| for every pattern - the F-H2 equivalence, by
% construction. Same key set as Index, so an absent key means mask 0 (dead cell).
build_masks(Index, MaskAssoc) :-
    assoc_to_list(Index, Pairs),
    maplist(ordset_mask_pair, Pairs, MaskPairs),
    list_to_assoc(MaskPairs, MaskAssoc).
ordset_mask_pair(K-Set, K-Mask) :- ordset_to_mask(Set, 0, Mask).
ordset_to_mask([], M, M).
ordset_to_mask([I|Is], M0, M) :- M1 is M0 \/ (1 << I), ordset_to_mask(Is, M1, M).

% LOAD: read the artifact, verify schema version + SWI version (+ the dict hash
% iff a dict path is supplied), and hand back the reconstructed DictByLen +
% Index. Any integrity failure throws a clear fill_index_* error - no silent
% rebuild, no fall-through to the raw path.
% /4: the pre-F-H2 signature, retained for callers/tests that only need the
% dictionary structures (e.g. the roundtrip identity test). Discards the masks.
fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index) :-
    fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, _Masks).

% /5: also returns the F-H2 counting context - masks(MaskAssoc) if the (v2)
% artifact carries a masks key, else none (the ordset kernel still works). Every
% integrity gate (version, SWI, hash) runs BEFORE masks are extracted, so a bad
% artifact throws exactly as before; masks come off the happy path only.
fill_load_index(IndexFile, DictFileOrNone, DictByLen, Index, Masks) :-
    ( exists_file(IndexFile) -> true
    ; throw(error(fill_index_missing(IndexFile), _)) ),
    % Keep the root cause (C33): the context slot carries the original error
    % term (fastrw version mismatch vs permissions vs truncation), so it shows
    % under verbose error printing instead of being swallowed.
    catch(fill_read_artifact(IndexFile, Artifact),
          ReadErr,
          throw(error(fill_index_unreadable(IndexFile), context(_, ReadErr)))),
    ( Artifact = fill_index(Version, Meta, DictByLen0, Index0)
    -> true
    ;  throw(error(fill_index_malformed(IndexFile), _)) ),
    fill_index_format_version(Want),
    ( Version == Want -> true
    ; throw(error(fill_index_version(Version, Want), _)) ),
    fill_swi_version(CurSwi),
    ( option(swi_version(ArtSwi), Meta), ArtSwi == CurSwi -> true
    ; option(swi_version(ArtSwi2), Meta, unknown),
      throw(error(fill_index_swi(ArtSwi2, CurSwi), _)) ),
    ( DictFileOrNone == none
    -> true
    ;  fill_file_sha256(DictFileOrNone, CurSha),
       ( option(dict_sha256(ArtSha), Meta), ArtSha == CurSha -> true
       ; option(dict_sha256(ArtSha2), Meta, unknown),
         throw(error(fill_index_hash(DictFileOrNone, ArtSha2, CurSha), _)) ) ),
    DictByLen = DictByLen0, Index = Index0,
    ( option(masks(MaskAssoc), Meta) -> Masks = masks(MaskAssoc) ; Masks = none ).

fill_read_artifact(IndexFile, Artifact) :-
    setup_call_cleanup(open(IndexFile, read, S, [type(binary)]),
                       fast_read(S, Artifact),
                       close(S)).

% Word count over the length buckets (same fold run_fill's metadata uses).
dict_word_count(DictByLen, N) :-
    assoc_to_values(DictByLen, Buckets),
    foldl(bucket_len, Buckets, 0, N).
bucket_len(B, A0, A1) :- length(B, L), A1 is A0 + L.

% SHA-256 of the file's raw bytes - a self-consistent build/verify fingerprint
% that also matches coreutils sha256sum (encoding(octet) hashes each byte as-is).
fill_file_sha256(File, Hex) :-
    read_file_to_string(File, Bytes, [encoding(octet)]),
    sha_hash(Bytes, Digest, [algorithm(sha256), encoding(octet)]),
    hash_atom(Digest, Hex).

% "Ma.Mi.Pa" - the same version string run_fill.pl records, so the artifact's
% guard and the bench provenance agree.
fill_swi_version(Ver) :-
    current_prolog_flag(version_data, V),
    ( V = swi(Ma, Mi, Pa, _) -> format(atom(Ver), '~d.~d.~d', [Ma, Mi, Pa]) ; Ver = V ).


% --- error messages ----------------------------------------------------------
:- multifile prolog:error_message//1.
prolog:error_message(fill_seed_no_slot(A)) -->
    [ 'fill: seed ~q does not match any slot of the grid (check its cells/direction)'-[A] ].
prolog:error_message(fill_seed_clash(A)) -->
    [ 'fill: seed ~q clashes with another seed at a shared cell'-[A] ].
prolog:error_message(fill_seed_duplicate(A)) -->
    [ 'fill: seed ~q is pinned more than once; answers must be unique'-[A] ].
% F-L2 index-artifact integrity failures (refuse, never silently rebuild).
prolog:error_message(fill_index_missing(F)) -->
    [ 'fill: index artifact ~w not found (build one with `fill --dict DICT --save-index ~w`)'-[F, F] ].
prolog:error_message(fill_index_unreadable(F)) -->
    [ 'fill: index artifact ~w could not be read (corrupt, or built by a different SWI-Prolog; rebuild with --save-index)'-[F] ].
prolog:error_message(fill_index_malformed(F)) -->
    [ 'fill: index artifact ~w is not a fill_index/4 artifact (rebuild with --save-index)'-[F] ].
prolog:error_message(fill_index_version(Got, Want)) -->
    [ 'fill: index artifact schema version ~w is not supported (this build reads version ~w); rebuild with --save-index'-[Got, Want] ].
prolog:error_message(fill_index_swi(Art, Cur)) -->
    [ 'fill: index artifact was built by SWI-Prolog ~w but this is ~w (the binary format is version-bound); rebuild with --save-index'-[Art, Cur] ].
prolog:error_message(fill_index_hash(File, Art, Cur)) -->
    [ 'fill: index artifact does not match --dict ~w (artifact dict SHA-256 ~w, file ~w); the dictionary changed - rebuild with --save-index'-[File, Art, Cur] ].
