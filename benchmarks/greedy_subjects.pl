% benchmarks/greedy_subjects.pl - product samplers and exact-replay observation
% twin for the greedy arrange benchmark. Product code is never instrumented.

:- module(greedy_subjects,
          [ load_words/2,
            construction_sampler/6,
            sweep_sampler/4,
            build_raw_pool/4,
            postprocess_sampler/6,
            command_sampler/7,
            semantic_counters/4,
            replay_equivalent/4,
            identity_row/10
          ]).

:- use_module(library(apply), [foldl/4, maplist/3]).
:- use_module(library(assoc), [assoc_to_list/2]).
:- use_module(library(sha), [sha_hash/3, hash_atom/2]).
:- use_module(library(lists), [member/2, selectchk/3]).
:- use_module(library(pairs), [pairs_values/2]).
:- use_module(library(process), [process_create/3, process_wait/2]).

:- use_module('bench_core.pl').

load_words(File, Words) :-
    load_clues(File, Words).

command_k(candidates(_, K), K).
command_k(best_effort, 1).

entry_answer([Answer|_], Answer).

entry_by_answer(Words, Answer, Entry) :-
    member(Entry, Words),
    Entry = [A|_],
    A == Answer,
    !.

% One pinned construction. The reset is part of the top-level operation and
% occurs exactly once; greedy_construct/6 itself deliberately does not reset.
construction_sampler(Words, GridLen, SeedAnswer, Corner, Expected, Sample) :-
    inproc_sampler(
        greedy_subjects:construction_operation(
            Words, GridLen, SeedAnswer, Corner, Expected, _Outcome),
        Sample).

construction_operation(Words, GridLen, SeedAnswer, Corner, Expected, Outcome) :-
    crosswordsmith_core:reset_search_memos,
    entry_by_answer(Words, SeedAnswer, Seed),
    (   crosswordsmith_arrange:greedy_construct(
            Words, GridLen, Corner, Seed, Placed, Dropped)
    ->  length(Placed, NP), length(Dropped, ND), Outcome = completed(NP, ND)
    ;   Outcome = setup_failed
    ),
    ( Outcome == Expected -> true
    ; throw(error(greedy_construction_outcome(Expected, Outcome), _))
    ).

% Primary gated subject: unchanged product greedy_construct/layout_reward calls
% over the direct member(Loc),member(Seed) product. No semantic counter executes
% here. The single reset is outside the corner/seed loops, so all constructions
% share product memos exactly as arrange_candidate_pool/4 does.
sweep_sampler(Words, GridLen, Command, Sample) :-
    inproc_sampler(greedy_subjects:sweep_operation(Words, GridLen, Command, _), Sample).

sweep_operation(Words, GridLen, Command, Raw) :-
    crosswordsmith_core:reset_search_memos,
    raw_pool_no_reset(Command, Words, GridLen, Raw).

build_raw_pool(Words, GridLen, Command, Raw) :-
    crosswordsmith_core:reset_search_memos,
    raw_pool_no_reset(Command, Words, GridLen, Raw).

% Exact pre-sort twin of arrange_candidate_pool/4, lines 1053-1064. Keep the
% conjunction and template in lockstep with product code: candidates Raw carries
% Placed only and strict eligibility is tested before layout_reward/4.
raw_pool_no_reset(candidates(DropContract, _K), Words, GridLen, Raw) :-
    crosswordsmith_arrange:arrange_weights(WCap, WTail),
    length(Words, Total),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    start_locs(Locs),
    findall(score(NP, Reward)-Placed,
            ( member(Loc, Locs),
              member(Seed, Seeds),
              crosswordsmith_arrange:greedy_construct(
                  Words, GridLen, Loc, Seed, Placed, _Dropped),
              length(Placed, NP),
              ( DropContract == strict -> NP =:= Total ; true ),
              crosswordsmith_arrange:layout_reward(
                  WCap, WTail, Placed, Reward) ),
            Raw).

% Exact pre-sort twin of arrange_best_effort/6, lines 418-428. Best-effort Raw
% carries ordered dropped ANSWERS, with that map performed inside the sweep.
raw_pool_no_reset(best_effort, Words, GridLen, Raw) :-
    crosswordsmith_arrange:arrange_weights(WCap, WTail),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    start_locs(Locs),
    findall(score(NP, Reward)-pd(Placed, DroppedAnswers),
            ( member(Loc, Locs),
              member(Seed, Seeds),
              crosswordsmith_arrange:greedy_construct(
                  Words, GridLen, Loc, Seed, Placed, DroppedEntries),
              length(Placed, NP),
              crosswordsmith_arrange:layout_reward(
                  WCap, WTail, Placed, Reward),
              maplist(entry_answer, DroppedEntries, DroppedAnswers) ),
            Raw).

% Raw is prebuilt by the caller outside call_time/2. Clause dispatch preserves
% each mode's actual product work rather than charging candidate diversity to
% best-effort.
postprocess_sampler(Raw, GridLen, Total, Command, K, Sample) :-
    inproc_sampler(
        greedy_subjects:postprocess(Command, Raw, GridLen, Total, K, _),
        Sample).

% Exact post-sweep work from arrange_candidate_pool/4 + arrange_candidates/6:
% stable rank, placement assocs, non-empty gate, diversity, and numbering.
postprocess(candidates(_Drop, _CommandK), Raw, GridLen, Total, K,
            candidates(Numbered, Returned)) :-
    sort(1, @>=, Raw, Sorted),
    pairs_values(Sorted, Placeds),
    maplist(crosswordsmith_arrange:tag_with_assoc(GridLen), Placeds, Pool),
    Pool = [_|_],
    crosswordsmith_arrange:candidate_tau_pct(TauPct),
    crosswordsmith_arrange:pick_diverse(Pool, TauPct, Total, K, Picked),
    maplist(crosswordsmith_arrange:numbered_candidate, Picked, Numbered),
    length(Numbered, Returned).

% Exact post-sweep work from arrange_best_effort/6: non-empty gate, stable
% best-first selection, and numbering. No assoc or diversity work occurs here.
postprocess(best_effort, Raw, _GridLen, _Total, _K,
            best_effort(Numbered, Reward, NumPlaced, Dropped)) :-
    Raw = [_|_],
    sort(1, @>=, Raw,
         [score(NumPlaced, Reward)-pd(BestPlaced, Dropped)|_]),
    crosswordsmith_core:assign_clue_numbers(BestPlaced, Numbered).

command_sampler(Exe, File, GridLen, Framing, Command, ExpectedExit,
                _{wall:Wall, rss:Rss}) :-
    command_args(File, GridLen, Framing, Command, '/dev/null', Args),
    bench_core:process_sampler(Exe, Args, Wall, Rss, Exit),
    ( Exit =:= ExpectedExit -> true
    ; throw(error(greedy_command_exit(ExpectedExit, Exit), _))
    ).

command_args(File, GridLen, Framing, Command, Out, Args) :-
    framing_args(Framing, GridLen, FrameArgs),
    mode_args(Command, ModeArgs),
    append([['arrange', '--input', File], FrameArgs, ModeArgs,
            ['--out', Out]], Args).

framing_args(size, N, ['--size', A]) :- atom_number(A, N).
framing_args(max_size, N, ['--max-size', A]) :- atom_number(A, N).

mode_args(candidates(strict, K), ['--strict', '--candidates', A]) :-
    atom_number(A, K).
mode_args(best_effort, ['--best-effort']).

% ---------------------------------------------------------------------------
% Benchmark-only exact replay twin and semantic counters.
% ---------------------------------------------------------------------------

new_counter(counter(0, 0, 0, 0, 0, 0)).

counter_inc(I, Counter) :-
    arg(I, Counter, N0),
    N is N0 + 1,
    nb_setarg(I, Counter, N).

counter_add(I, Add, Counter) :-
    arg(I, Counter, N0),
    N is N0 + Add,
    nb_setarg(I, Counter, N).

counter_dict(C, _{generated_crossing_descriptors:Generated,
                  legality_probes:Probes,
                  scored_candidates:Scored,
                  score_cell_visits:Visits,
                  greedy_steps:Steps,
                  completed_constructions:Completed}) :-
    arg(1, C, Generated), arg(2, C, Probes), arg(3, C, Scored),
    arg(4, C, Visits), arg(5, C, Steps), arg(6, C, Completed).

% Counter semantics: generated increments for every descriptor yielded by
% find_intersecting_word/6; scored after placement_key/8 succeeds; score-cell
% visits counts both full word walks performed by placement_key (crossing count
% and bbox extension); legality increments only after scoring, preserving the
% current score-before-legality order; a greedy step is one realized non-seed
% winner; completed is one construction whose seed setup and loop complete.
semantic_counters(Words, GridLen, Command, Counters) :-
    crosswordsmith_core:reset_search_memos,
    new_counter(C),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    start_locs(Locs),
    forall((member(Loc, Locs), member(Seed, Seeds)),
           replay_pair_checked(Words, GridLen, Loc, Seed, C)),
    % Command is intentionally accepted to pin the same manifest subject. The
    % semantic sweep observes every attempted direct construction, including
    % setup failures, before strict eligibility filtering.
    nonvar(Command),
    counter_dict(C, Counters).

replay_equivalent(Words, GridLen, Corner, SeedAnswer) :-
    entry_by_answer(Words, SeedAnswer, Seed),
    crosswordsmith_core:reset_search_memos,
    construction_result(product, Words, GridLen, Corner, Seed, Product),
    crosswordsmith_core:reset_search_memos,
    new_counter(C),
    construction_result(replay(C), Words, GridLen, Corner, Seed, Replay),
    equivalent_result(Product, Replay).

replay_pair_checked(Words, GridLen, Corner, Seed, Counter) :-
    construction_result(product, Words, GridLen, Corner, Seed, Product),
    construction_result(replay(Counter), Words, GridLen, Corner, Seed, Replay),
    ( equivalent_result(Product, Replay) -> true
    ; Seed = [Answer|_],
      throw(error(greedy_replay_diverged(Answer, Corner, Product, Replay), _))
    ).

construction_result(product, Words, GridLen, Corner, Seed, Result) :-
    ( crosswordsmith_arrange:greedy_construct(
          Words, GridLen, Corner, Seed, Placed, Dropped)
    -> Result = ok(Placed, Dropped)
    ;  Result = setup_failed
    ).
construction_result(replay(C), Words, GridLen, Corner, Seed, Result) :-
    ( replay_construct(Words, GridLen, Corner, Seed, C, Placed, Dropped)
    -> Result = ok(Placed, Dropped)
    ;  Result = setup_failed
    ).

equivalent_result(setup_failed, setup_failed).
equivalent_result(ok(P1, D1), ok(P2, D2)) :-
    P1 =@= P2,
    maplist(entry_answer, D1, A1),
    maplist(entry_answer, D2, A2),
    A1 == A2,
    crosswordsmith_arrange:arrange_weights(WC, WT),
    crosswordsmith_arrange:layout_reward(WC, WT, P1, R1),
    crosswordsmith_arrange:layout_reward(WC, WT, P2, R2),
    R1 =:= R2.

replay_construct(Words, GridLen, Corner, Seed, Counter, Placed, Dropped) :-
    crosswordsmith_core:init_gs(GridLen, G0),
    crosswordsmith_core:start_loc(Corner, GridLen, Start, Dir),
    crosswordsmith_core:remove_x(Seed, Words, Rest),
    crosswordsmith_arrange:seed_word(Seed, Start, Dir, GridLen, G0, SeedPW, G1),
    replay_loop(Rest, [SeedPW], GridLen, G1, Counter, Placed, Dropped),
    counter_inc(6, Counter).

replay_loop(Remaining, Placed, GridLen, Grid, Counter, Final, Dropped) :-
    replay_next_move(Remaining, Placed, GridLen, Grid, Counter, Move),
    replay_apply_move(Move, Remaining, Placed, GridLen, Counter, Final, Dropped).

replay_apply_move(none, Remaining, Placed, _GridLen, _Counter, Placed, Remaining).
replay_apply_move(move(Answer, NewPW, NewGrid), Remaining, Placed, GridLen,
                  Counter, Final, Dropped) :-
    selectchk([Answer|_], Remaining, Remaining1),
    counter_inc(5, Counter),
    replay_loop(Remaining1, [NewPW|Placed], GridLen, NewGrid,
                Counter, Final, Dropped).

replay_next_move(Remaining, Placed, GridLen, Grid, Counter, Move) :-
    crosswordsmith_metrics:placed_bbox(Placed, GridLen, BBox, _),
    findall(Score-Best,
            ( member(Entry, Remaining),
              replay_word_best(Entry, Placed, GridLen, Grid, BBox,
                               Counter, Score, Best) ),
            Candidates),
    replay_best_move(Candidates, GridLen, Grid, Move).

replay_best_move([], _GridLen, _Grid, none).
replay_best_move([C|Cs], GridLen, Grid, move(Answer, PW, G1)) :-
    foldl(crosswordsmith_arrange:max_key_first, Cs, C,
          _-best(Answer, Letters, WLen, Start, Dir)),
    crosswordsmith_core:assign_word(
        Answer, Letters, WLen, Start, Dir, GridLen, Grid, PW, G1).

replay_word_best(Entry, Placed, GridLen, Grid, BBox, Counter, Score, Best) :-
    crosswordsmith_metrics:word_letters(Entry, Letters, WLen),
    Entry = [Answer|_],
    Grid = gs(LGrid, _),
    findall(Key-(Start-Dir),
            ( crosswordsmith_core:find_intersecting_word(
                  Letters, WLen, Placed, GridLen, Start, Dir),
              counter_inc(1, Counter),
              crosswordsmith_arrange:placement_key(
                  Letters, Start, Dir, WLen, GridLen, LGrid, BBox, Key),
              counter_inc(3, Counter),
              Visits is 2 * WLen,
              counter_add(4, Visits, Counter),
              counter_inc(2, Counter),
              crosswordsmith_core:check_word_fits(
                  Letters, Start, Dir, GridLen, Grid) ),
            Keyed),
    Keyed = [K0|Ks],
    foldl(crosswordsmith_arrange:max_key_first, Ks, K0,
          Score-(BestStart-BestDir)),
    Best = best(Answer, Letters, WLen, BestStart, BestDir).

% ---------------------------------------------------------------------------
% Semantic identity document.
% ---------------------------------------------------------------------------

identity_row(Id, File, GridLen, Framing, Command, SeedAnswer, Corner,
             Expected, Exe, Row) :-
    load_words(File, Words),
    format(user_error, "heartbeat: identity ~w direct-attempts~n", [Id]),
    direct_attempts(Words, GridLen, Command, Attempts),
    format(user_error, "heartbeat: identity ~w raw-pool~n", [Id]),
    identity_raw_pool(Words, GridLen, Command, Raw),
    format(user_error, "heartbeat: identity ~w selected~n", [Id]),
    selected_layouts(Words, GridLen, Command, Selected),
    selected_signatures(Words, GridLen, Framing, Selected,
                        SelectedSigs, Distances),
    format(user_error, "heartbeat: identity ~w cli~n", [Id]),
    cli_identity(Exe, File, GridLen, Framing, Command, Cli),
    expected_json(Expected, ExpectedJson),
    command_json(Command, CommandJson),
    length(SelectedSigs, CandidateCount),
    Row = _{rung:Id, fixture:File, size:GridLen, framing:Framing,
            mode:CommandJson, construction:_{seed_answer:SeedAnswer,
                corner:Corner, expected:ExpectedJson},
            direct_attempts:Attempts, raw_pool:Raw, selected:SelectedSigs,
            candidate_count:CandidateCount,
            pairwise_distances:Distances, cli:Cli}.

% Identity-only observer: one explicit row for every direct corner x seed slot,
% including failed setup and strict-ineligible completed constructions. It runs
% under its own reset and never participates in gated sweep measurement.
direct_attempts(Words, GridLen, Command, Attempts) :-
    crosswordsmith_core:reset_search_memos,
    crosswordsmith_arrange:arrange_weights(WC, WT),
    crosswordsmith_arrange:seed_candidates(Words, Seeds),
    start_locs(Locs),
    length(Words, Total),
    findall(Attempt,
            ( member(Corner, Locs), member(Seed, Seeds),
              Seed = [SeedAnswer|_],
              direct_attempt(Words, GridLen, Command, Total, WC, WT,
                             Corner, Seed, SeedAnswer, Attempt) ),
            Attempts).

direct_attempt(Words, GridLen, Command, Total, WC, WT,
               Corner, Seed, SeedAnswer, Attempt) :-
    (   crosswordsmith_arrange:greedy_construct(
            Words, GridLen, Corner, Seed, Placed, Dropped)
    ->  length(Placed, NP), length(Dropped, ND),
        attempt_eligible(Command, NP, Total, Eligible),
        crosswordsmith_arrange:layout_reward(WC, WT, Placed, Reward),
        placed_signature(Placed, PlacedSig),
        maplist(entry_answer, Dropped, DroppedSig),
        Attempt = _{corner:Corner,seed_answer:SeedAnswer,outcome:completed,
                    eligibility:Eligible,
                    score:_{placed:NP,dropped:ND,reward:Reward},
                    placed_signature:PlacedSig,dropped_signature:DroppedSig}
    ;   Attempt = _{corner:Corner,seed_answer:SeedAnswer,outcome:setup_failed,
                    eligibility:false,score:null,
                    placed_signature:null,dropped_signature:null}
    ).

attempt_eligible(candidates(Drop, _), NP, Total, Eligible) :-
    ( Drop == strict
    -> ( NP =:= Total -> Eligible=true ; Eligible=false )
    ;  Eligible=true
    ).
attempt_eligible(best_effort, _NP, _Total, true).

% The exact entries that reach product's pre-sort pool, observed under a fresh
% reset so the richer direct-attempt traversal cannot perturb memo lifecycle.
identity_raw_pool(Words, GridLen, Command, Rows) :-
    build_raw_pool(Words, GridLen, Command, Raw),
    maplist(identity_raw_entry(Command), Raw, Rows).

identity_raw_entry(candidates(_, _), score(NP, Reward)-Placed,
                   _{score:_{placed:NP,reward:Reward},
                     placed_signature:PlacedSig}) :-
    !,
    placed_signature(Placed, PlacedSig).
identity_raw_entry(best_effort,
                   score(NP, Reward)-pd(Placed, DroppedAnswers),
                   _{score:_{placed:NP,reward:Reward},
                     placed_signature:PlacedSig,
                     dropped_answers:DroppedAnswers}) :-
    placed_signature(Placed, PlacedSig).

selected_layouts(Words, GridLen, candidates(Drop, K), Layouts) :-
    crosswordsmith_arrange:arrange_candidates(
        Words, GridLen, Drop, K, Layouts, _).
selected_layouts(Words, GridLen, best_effort, [Layout]) :-
    crosswordsmith_arrange:arrange_best_effort(
        Words, GridLen, Layout, _Reward, _NP, _Dropped).

selected_signatures(Words, GridLen, Framing, Layouts, Signatures, Distances) :-
    maplist(selected_signature(Words, GridLen, Framing), Layouts, Signatures),
    maplist(layout_assoc(GridLen), Layouts, Assocs),
    findall(_{left:I,right:J,distance:Diff},
            ( nth0(I, Assocs, A1), nth0(J, Assocs, A2), I < J,
              crosswordsmith_arrange:pos_diff_count(A1, A2, Diff) ),
            Distances).

layout_assoc(GridLen, Layout, Assoc) :-
    crosswordsmith_arrange:placement_assoc(Layout, GridLen, Assoc).

selected_signature(Words, GridLen, Framing, Layout, Sig) :-
    crosswordsmith_arrange:arrange_weights(WC, WT),
    crosswordsmith_arrange:layout_reward(WC, WT, Layout, Reward),
    crosswordsmith_arrange:dropped_answers(Words, Layout, Dropped),
    placed_signature(Layout, Placed),
    crosswordsmith_arrange:placement_assoc(Layout, GridLen, Assoc),
    assoc_signature(Assoc, AssocSig),
    emit_framing(Framing, SizeMode),
    with_output_to(string(Bytes),
        crosswordsmith_arrange:emit_arrange(Layout, Words, GridLen, SizeMode)),
    sha256(Bytes, OutputSha),
    Sig = _{reward:Reward, dropped:Dropped, placed_signature:Placed,
            normalized_assoc:AssocSig, output_sha256:OutputSha}.

emit_framing(size, fixed).
emit_framing(max_size, max).

placed_signature(Placed, Signature) :-
    maplist(placed_word_signature, Placed, Signature).

placed_word_signature(PW, Sig) :-
    pw_answer(PW, Answer), pw_start(PW, Start), pw_dir(PW, Dir),
    pw_cells(PW, Cells), pw_num(PW, Num0),
    ( integer(Num0) -> Num = Num0 ; Num = null ),
    Sig = _{answer:Answer,start:Start,direction:Dir,cells:Cells,number:Num}.

assoc_signature(Assoc, Signature) :-
    assoc_to_list(Assoc, Pairs),
    maplist(assoc_pair_signature, Pairs, Signature).

assoc_pair_signature(Answer-(Row-Col-Dir),
                     _{answer:Answer,row:Row,col:Col,direction:Dir}).

cli_identity(Exe, File, GridLen, Framing, Command, Sig) :-
    command_args(File, GridLen, Framing, Command, '/dev/stdout', Args),
    process_create(Exe, Args,
                   [stdout(pipe(Out)), stderr(pipe(Err)), process(PID)]),
    read_string(Out, _, Stdout), close(Out),
    read_string(Err, _, Stderr), close(Err),
    process_wait(PID, Status),
    ( Status = exit(Exit) -> true ; Exit = -1 ),
    sha256(Stdout, OutSha), sha256(Stderr, ErrSha),
    Sig = _{exit:Exit,stdout_sha256:OutSha,stderr_sha256:ErrSha}.

sha256(Data, Hash) :-
    sha_hash(Data, Bytes, [algorithm(sha256)]),
    hash_atom(Bytes, Hash).

expected_json(completed(P, D), _{outcome:completed,placed:P,dropped:D}).
expected_json(setup_failed, _{outcome:setup_failed}).

command_json(candidates(strict, K), _{kind:candidates,drop:strict,k:K}).
command_json(best_effort, _{kind:best_effort}).

:- multifile prolog:error_message//1.
prolog:error_message(greedy_construction_outcome(Expected, Got)) -->
    [ 'greedy benchmark construction expected ~q, got ~q'-[Expected, Got] ].
prolog:error_message(greedy_command_exit(Expected, Got)) -->
    [ 'greedy benchmark command expected exit ~w, got ~w'-[Expected, Got] ].
prolog:error_message(greedy_replay_diverged(Answer, Corner, _, _)) -->
    [ 'greedy benchmark replay diverged for seed ~q at ~w'-[Answer, Corner] ].
