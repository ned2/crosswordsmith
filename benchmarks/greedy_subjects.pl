% benchmarks/greedy_subjects.pl - product samplers and semantic identity for
% the greedy arrange benchmark. Product code is never instrumented.

:- module(greedy_subjects,
           [ construction_sampler/6,
             sweep_sampler/4,
             build_raw_pool/4,
             postprocess_sampler/6,
             command_sampler/7,
             identity_row/11
           ]).

:- use_module(library(apply), [maplist/3]).
:- use_module(library(assoc), [assoc_to_list/2]).
:- use_module(library(sha), [sha_hash/3, hash_atom/2]).
:- use_module(library(lists), [append/2, member/2, nth0/3]).
:- use_module(library(pairs), [pairs_values/2]).
:- use_module('bench_core.pl', [inproc_sampler/2, process_sampler/5]).
:- use_module('bench_fixture.pl', [load_arrange_fixture/2]).
:- use_module('bench_process.pl', [capture_process/6]).

entry_answer([Answer|_], Answer).

entry_by_answer(Words, Answer, Entry) :-
    member(Entry, Words),
    Entry = [A|_],
    A == Answer,
    !.

%!  construction_sampler(+Words:list, +GridLen:integer, +SeedAnswer:atom,
%!                        +Corner:atom, +Expected, -Sample:dict) is det.
%
% One pinned construction. The reset is part of the top-level operation and
% occurs exactly once; greedy_construct/6 itself deliberately does not reset.
% SeedAnswer must name an entry in Words, and Corner/Expected must come from the
% workload manifest.
construction_sampler(Words, GridLen, SeedAnswer, Corner, Expected, Sample) :-
    require_seed_answer(Words, SeedAnswer),
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

%!  sweep_sampler(+Words:list, +GridLen:integer, +Command,
%!                -Sample:dict) is det.
%
% Primary gated subject: the product's two directly searched blocks plus their
% derived visible partners. No semantic counter executes here. The single reset
% is outside the construction helper, so both direct blocks share product memos
% exactly as arrange_candidate_pool/4 does. Command is candidates(strict,K),
% with positive K, or best_effort.
sweep_sampler(Words, GridLen, Command, Sample) :-
    require_command(Command),
    inproc_sampler(greedy_subjects:sweep_operation(Words, GridLen, Command, _), Sample).

sweep_operation(Words, GridLen, Command, Raw) :-
    crosswordsmith_core:reset_search_memos,
    raw_pool_no_reset(Command, Words, GridLen, Raw).

%!  build_raw_pool(+Words:list, +GridLen:integer, +Command, -Raw:list) is det.
%
build_raw_pool(Words, GridLen, Command, Raw) :-
    require_command(Command),
    crosswordsmith_core:reset_search_memos,
    raw_pool_no_reset(Command, Words, GridLen, Raw).

% Exact pre-sort twin of arrange_candidate_pool/4. Candidates Raw carries Placed
% only and strict eligibility is tested before layout_reward/4.
raw_pool_no_reset(candidates(DropContract, _K), Words, GridLen, Raw) :-
    crosswordsmith_arrange:arrange_weights(WCap, WTail),
    length(Words, Total),
    crosswordsmith_arrange:greedy_constructions(Words, GridLen, Constructions),
    findall(score(NP, Reward)-Placed,
            ( member(gc(Placed, _Dropped), Constructions),
              length(Placed, NP),
              ( DropContract == strict -> NP =:= Total ; true ),
              crosswordsmith_arrange:layout_reward(
                  WCap, WTail, Placed, Reward) ),
            Raw).

% Exact pre-sort twin of arrange_best_effort/6. Best-effort Raw carries ordered
% dropped ANSWERS, with that map performed inside the sweep.
raw_pool_no_reset(best_effort, Words, GridLen, Raw) :-
    crosswordsmith_arrange:arrange_weights(WCap, WTail),
    crosswordsmith_arrange:greedy_constructions(Words, GridLen, Constructions),
    findall(score(NP, Reward)-pd(Placed, DroppedAnswers),
            ( member(gc(Placed, DroppedEntries), Constructions),
              length(Placed, NP),
              crosswordsmith_arrange:layout_reward(
                  WCap, WTail, Placed, Reward),
              maplist(entry_answer, DroppedEntries, DroppedAnswers) ),
            Raw).

% Raw is prebuilt by the caller outside call_time/2. Clause dispatch preserves
% each mode's actual product work rather than charging candidate diversity to
% best-effort.
%!  postprocess_sampler(+Raw:list, +GridLen:integer, +Total:integer, +Command,
%!                      +K:integer, -Sample:dict) is det.
%
%   Raw must be the nonempty, mode-compatible pool from build_raw_pool/4, and
%   Command must satisfy sweep_sampler/4's command contract.
postprocess_sampler(Raw, GridLen, Total, Command, K, Sample) :-
    require_command(Command),
    require_nonempty_raw(Raw),
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

%!  command_sampler(+Exe, +File:atom, +GridLen:integer, +Framing:atom,
%!                   +Command, +ExpectedExit:integer, -Sample:dict) is det.
%
%   Framing is size or max_size; Command satisfies sweep_sampler/4's contract.
command_sampler(Exe, File, GridLen, Framing, Command, ExpectedExit,
                _{wall:Wall, rss:Rss}) :-
    require_framing(Framing),
    require_command(Command),
    command_args(File, GridLen, Framing, Command, '/dev/null', Args),
    process_sampler(Exe, Args, Wall, Rss, Exit),
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
% Semantic identity document.
% ---------------------------------------------------------------------------

%!  identity_row(+Id, +Fixture, +File, +GridLen, +Framing, +Command,
%!               +SeedAnswer, +Corner, +Expected, +Exe, -Row) is det.
%
%   Build one semantic identity row. Fixture is the stable repository-relative
%   label recorded in the document; File is its absolute operational path. The
%   framing, command, seed, corner, and expected outcome come from the manifest.
identity_row(Id, Fixture, File, GridLen, Framing, Command, SeedAnswer, Corner,
             Expected, Exe, Row) :-
    require_framing(Framing),
    require_command(Command),
    require_expected(Expected),
    load_arrange_fixture(File, Words),
    require_seed_answer(Words, SeedAnswer),
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
    Row = _{rung:Id, fixture:Fixture, size:GridLen, framing:Framing,
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
    !,
    ( crosswordsmith_arrange:arrange_candidates(
          Words, GridLen, Drop, K, Selected, _)
    -> Layouts = Selected
    ;  Layouts = []
    ).
selected_layouts(Words, GridLen, best_effort, Layouts) :-
    ( crosswordsmith_arrange:arrange_best_effort(
          Words, GridLen, Layout, _Reward, _NP, _Dropped)
    -> Layouts = [Layout]
    ;  Layouts = []
    ).

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
    capture_process(Exe, Args, capture, Stdout, Stderr, Status),
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

require_seed_answer(Words, SeedAnswer) :-
    ( entry_by_answer(Words, SeedAnswer, _) -> true
    ; throw(error(greedy_seed_answer_missing(SeedAnswer), _))
    ).

require_command(candidates(strict, K)) :- integer(K), K >= 1, !.
require_command(best_effort) :- !.
require_command(Command) :-
    throw(error(greedy_benchmark_command(Command), _)).

require_framing(size) :- !.
require_framing(max_size) :- !.
require_framing(Framing) :-
    throw(error(greedy_benchmark_framing(Framing), _)).

require_expected(completed(Placed, Dropped)) :-
    integer(Placed), Placed >= 0,
    integer(Dropped), Dropped >= 0,
    !.
require_expected(setup_failed) :- !.
require_expected(Expected) :-
    throw(error(greedy_benchmark_expected(Expected), _)).

require_nonempty_raw([_|_]) :- !.
require_nonempty_raw([]) :-
    throw(error(greedy_empty_raw_pool, _)).

:- multifile prolog:error_message//1.
prolog:error_message(greedy_construction_outcome(Expected, Got)) -->
    [ 'greedy benchmark construction expected ~q, got ~q'-[Expected, Got] ].
prolog:error_message(greedy_command_exit(Expected, Got)) -->
    [ 'greedy benchmark command expected exit ~w, got ~w'-[Expected, Got] ].
