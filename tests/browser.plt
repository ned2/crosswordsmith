% tests/browser.plt - plunit suite for browser.pl (the WASM dispatch spine).
%
% Assumes the project is loaded (via load.pl) by the runner
% (tests/run_tests.pl) before this file. browser.pl's internals are reached as
% crosswordsmith_browser:Pred(...) — never add an export for a test. Run via:
%
%     ./run_tests.sh        (or)        make test
%
% This suite is the regression lock for wasm-sdk-strategy.md §7's headline
% invariant: every dispatch path yields a TYPED envelope — a tagged
% success/failure or a mapped error — and never the bare-fail path that
% produced the spike's generic "no layout (search failed)" masquerade.
% Covers: the malformed-params fuzz, protocol errors, the edge-case matrix,
% same-instance determinism (the reset lock), native value-parity against the
% committed CLI goldens (DEC-8: value-equality of the live `result` dict, not
% bytes), resolver parity with the CLI binary, and classifier totality.

:- use_module(library(plunit)).
% library(json), not the legacy library(http/json) alias (mirrors browser.pl,
% the module under test — the alias does not resolve in the WASM image).
:- use_module(library(json),
              [ atom_json_dict/3,
                json_read_dict/3,
                json_write_dict/2
              ]).
:- use_module(library(process)).
% lists/apply/aggregate: test-body helpers; explicit so the suite also runs
% under autoload(false) (P11/C5).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(aggregate)).
:- use_module(library(yall)).    % test-body lambdas
:- use_module(library(debug)).   % assertion/1
:- use_module(library(readutil)). % read_file_to_string/3 (golden parity)

% --- helpers (bt_ prefix: .plt helpers share the `user` namespace) -----------

bt_toy_clues([ _{answer: "CAT"}, _{answer: "CAR"}, _{answer: "ARC"},
               _{answer: "RAT"}, _{answer: "TAR"} ]).

bt_toy_params(P) :-
    bt_toy_clues(Clues),
    P = _{clues: Clues, size: 5}.

% Round-trip one request through browser_dispatch/3 exactly as the Worker
% would: envelope dict -> JSON payload atom -> dispatch -> JSON.parse'd reply.
bt_dispatch(Verb, ReqDict, Env) :-
    atom_json_dict(Payload, ReqDict, []),
    bt_dispatch_payload(Verb, Payload, Env).

bt_dispatch_payload(Verb, Payload, Env) :-
    browser_dispatch(Verb, Payload, Json),
    atom_json_dict(Json, Env, [default_tag(json)]).

bt_arrange(Params, Env) :-
    bt_dispatch(arrange, _{v: 1, id: "t-1", verb: "arrange", params: Params}, Env).

% The masquerade lock: an error envelope must be typed and must never carry
% the spike's generic wording.
bt_error(Env, Type) :-
    get_dict(status, Env, "error"),
    get_dict(error, Env, E),
    get_dict(type, E, Type),
    get_dict(message, E, Msg),
    string(Msg), Msg \== "",
    \+ sub_string(Msg, _, _, _, "no layout").

bt_failure(Env, Reason) :-
    get_dict(status, Env, "failure"),
    get_dict(detail, Env, D),
    get_dict(reason, D, Reason).

bt_success_result(Env, R) :-
    get_dict(status, Env, "success"),
    get_dict(result, Env, R).

% Fixture -> browser params (the .pl clue fixtures feed the goldens; the
% browser path speaks JSON, so translate [[Answer, Meta], ...] verbatim).
bt_word_clue([A, M], _{answer: S, meta: M}) :- !, atom_string(A, S).
bt_word_clue([A], _{answer: S}) :- atom_string(A, S).

bt_fixture_params(Fixture, Extra, Params) :-
    load_clues(Fixture, Words),
    maplist(bt_word_clue, Words, Clues),
    put_dict(Extra, _{clues: Clues}, Params).

bt_golden(File, Golden) :-
    setup_call_cleanup(open(File, read, S),
                       json_read_dict(S, Golden, [default_tag(json)]),
                       close(S)).

% Run the real CLI binary and capture its exit code (stdout/stderr discarded)
% — the resolver-parity reference. Assumes cwd = repo root (run_tests.sh).
bt_cli_exit(Args, Code) :-
    process_create('./crosswordsmith', Args,
                   [stdout(null), stderr(null), process(PID)]),
    process_wait(PID, exit(Code)).

% lint/export both take a canonical layout; the committed arrange golden is
% the shared input — exactly what `make regen-goldens` feeds the CLI to
% produce the lint/export goldens, so deep-equality against those locks the
% browser path to the CLI byte-for-value.
bt_layout(Layout) :-
    bt_golden('tests/golden/arrange_bundled_17_fixed.json', Layout).

bt_lint(Params, Env) :-
    bt_dispatch(lint, _{v: 1, id: "t-l", verb: "lint", params: Params}, Env).

bt_export(Params, Env) :-
    bt_dispatch(export, _{v: 1, id: "t-e", verb: "export", params: Params}, Env).

% --- table-memo census (the C48 bounded-growth lock) --------------------------
% current_table/2 VARIANT-matches its goal argument, so an open pattern like
% pair_crossings(_,_,_) silently enumerates NOTHING (probe-verified on
% SWI 10.1.10; see core.pl reset_search_memos/0): enumerate unbound + filter.
bt_memo_table(Trie) :-
    current_table(M:G, Trie),
    M == crosswordsmith_core,
    functor(G, Name, Arity),
    memberchk(Name/Arity, [pair_crossings/3, answer_letters/2]).

bt_memo_census(Variants, Bytes) :-
    aggregate_all(count, bt_memo_table(_), Variants),
    aggregate_all(sum(S), ( bt_memo_table(T), trie_property(T, size(S)) ),
                  Bytes).

% Deterministic, globally-unique 5-letter vocabulary per set index (the "user
% edits their clue list between requests" story: every dispatch presents fresh
% (Letters, PLetters) pairs to the memo, the worst case for table growth).
% The first two letters encode SetIdx, so no word ever recurs across sets.
bt_vocab_word(SetIdx, WIdx, W) :-
    A is 0'A + SetIdx mod 26,
    B is 0'A + (SetIdx // 26) mod 26,
    C is 0'A + WIdx mod 26,
    D is 0'A + (SetIdx * 3 + WIdx * 5) mod 26,
    E is 0'A + (SetIdx * 7 + WIdx * 11) mod 26,
    string_codes(W, [A, B, C, D, E]).

bt_fresh_vocab_params(SetIdx, _{clues: Clues, size: 9, bestEffort: true}) :-
    numlist(1, 12, Is),
    maplist(bt_vocab_word(SetIdx), Is, Ws0),
    sort(Ws0, Ws),                       % dedupe: duplicate answers throw
    findall(_{answer: W}, member(W, Ws), Clues).

% load_clues/2 dispatches on the file extension, so the temp file must be .json.
bt_with_toy_input(Goal) :-
    tmp_file(bt_clues, Base),
    atom_concat(Base, '.json', File),
    bt_toy_params(P),
    setup_call_cleanup(
        setup_call_cleanup(open(File, write, S),
                           json_write_dict(S, _{clues: P.clues}),
                           close(S)),
        call(Goal, File),
        catch(delete_file(File), _, true)).


:- begin_tests(browser).

% --- envelope shape ----------------------------------------------------------

test(success_envelope_shape) :-
    bt_toy_params(P),
    bt_arrange(P, Env),
    get_dict(v, Env, 1),
    get_dict(id, Env, "t-1"),
    get_dict(verb, Env, "arrange"),
    bt_success_result(Env, R),
    get_dict(gridLength, R, 5),
    get_dict(words, R, Ws), length(Ws, 5),
    get_dict(diagnostics, R, _).

% An absent id echoes as JSON null (the literal, not the string "null" — the
% exact corruption JSON-both-ways exists to prevent).
test(absent_id_echoes_json_null) :-
    bt_toy_params(P),
    atom_json_dict(Payload, _{v: 1, verb: "arrange", params: P}, []),
    browser_dispatch(arrange, Payload, Json),
    once(sub_atom(Json, _, _, _, '"id":null')),
    \+ sub_atom(Json, _, _, _, '"id":"null"').

% id wrong-type shapes (request_id/2's atomic gate): atomic ids echo as
% themselves (number; boolean/null land as the atoms true/null after the
% parse, which serialize back to the JSON literals); structured ids
% (object/array) have no usable atomic echo -> JSON null, the absent-id rule.
test(id_shapes_echo_or_null) :-
    forall(member(IdJson-Expected,
                  [ '42'-'"id":42',
                    'true'-'"id":true',
                    '{"a":1}'-'"id":null',
                    '[1,2]'-'"id":null',
                    'null'-'"id":null' ]),
           ( format(atom(P), '{"v":1,"id":~w,"verb":"capabilities"}', [IdJson]),
             browser_dispatch(capabilities, P, Json),
             once(sub_atom(Json, _, _, _, Expected)) )).

% --- the malformed-params invariant fuzz (strategy §7's headline lock) -------
% Every case must yield a TYPED validation error — never a bare fail, never
% the generic masquerade. Each entry mutates the known-good toy params.

bt_bad_params(seed_neg1,        _{seed: -1}).
bt_bad_params(seed_neg,         _{seed: -5}).
bt_bad_params(seed_float,       _{seed: 3.5}).
bt_bad_params(seed_string,      _{seed: "42"}).
bt_bad_params(seed_1e21_float,  _{seed: 1.0e21}).
bt_bad_params(size_zero,        _{size: 0}).
bt_bad_params(size_neg,         _{size: -3}).
bt_bad_params(size_string,      _{size: "5"}).
bt_bad_params(size_float,       _{size: 2.5}).
bt_bad_params(mode_unknown,     _{mode: "diagonal"}).
bt_bad_params(mode_number,      _{mode: 5}).
bt_bad_params(besteffort_string, _{bestEffort: "yes"}).
bt_bad_params(besteffort_number, _{bestEffort: 1}).
bt_bad_params(seed_with_besteffort, _{seed: 7, bestEffort: true}).

test(fuzz_bad_params_typed_validation, [forall(bt_bad_params(_Name, Mut))]) :-
    bt_toy_params(P0),
    put_dict(Mut, P0, P),
    bt_arrange(P, Env),
    bt_error(Env, "validation").

% Clue-schema violations ride the same channel (core's existing throws).
test(fuzz_missing_clues) :-
    bt_arrange(_{size: 5}, Env), bt_error(Env, "validation").
test(fuzz_clues_not_list) :-
    bt_arrange(_{clues: "x", size: 5}, Env), bt_error(Env, "validation").
test(fuzz_answer_not_string) :-
    bt_arrange(_{clues: [_{answer: 42}], size: 5}, Env),
    bt_error(Env, "validation").
test(fuzz_duplicate_answers) :-
    bt_arrange(_{clues: [_{answer: "CAT"}, _{answer: "CAT"}], size: 5}, Env),
    bt_error(Env, "validation").
test(fuzz_missing_params) :-
    bt_dispatch(arrange, _{v: 1, id: "t-1", verb: "arrange"}, Env),
    bt_error(Env, "validation").
test(fuzz_params_not_object) :-
    bt_dispatch(arrange, _{v: 1, id: "t-1", verb: "arrange", params: "x"}, Env),
    bt_error(Env, "validation").

% The C56 pins: `size` is capped browser-side (max_browser_size/1 — the ONE
% deliberate divergence from the CLI's accepted domain: grid terms are
% size^2-arity compounds, and an uncapped size races the uncatchable WASM
% abort() on device). cap+1 -> typed validation naming the bound; an at-cap
% request still solves.
test(size_cap_plus_one_typed_validation) :-
    crosswordsmith_browser:max_browser_size(Cap),
    Over is Cap + 1,
    bt_toy_clues(Clues),
    bt_arrange(_{clues: Clues, size: Over}, Env),
    bt_error(Env, "validation"),
    get_dict(error, Env, E),
    get_dict(message, E, Msg),
    once(sub_string(Msg, _, _, _, "at most")).
test(size_at_cap_still_solves) :-
    crosswordsmith_browser:max_browser_size(Cap),
    bt_toy_clues(Clues),
    bt_arrange(_{clues: Clues, size: Cap}, Env),
    bt_success_result(Env, R),
    get_dict(gridLength, R, Cap).

% A huge seed that is still a JSON *integer* stays in the CLI's accepted
% domain (integer >= 0 — SWI reads bignums); only a float-degraded one is bad.
test(huge_integer_seed_accepted) :-
    bt_dispatch_payload(arrange,
        '{"v":1,"id":"t-1","verb":"arrange","params":{"clues":[{"answer":"CAT"},{"answer":"CAR"},{"answer":"ARC"},{"answer":"RAT"},{"answer":"TAR"}],"size":5,"seed":123456789012345678901234567890}}',
        Env),
    bt_success_result(Env, _).

% --- protocol errors ---------------------------------------------------------

% `fill` is deliberately NOT on the spine yet (strategy §6 phase 4) — it makes
% the perfect unknown-verb exemplar until it lands.
test(unknown_verb_typed) :-
    bt_dispatch(fill, _{v: 1, id: "t-1", verb: "fill", params: _{}}, Env),
    bt_error(Env, "unknown_verb").
test(unsupported_version_typed) :-
    bt_toy_params(P),
    bt_dispatch(arrange, _{v: 2, id: "t-1", verb: "arrange", params: P}, Env),
    bt_error(Env, "unsupported_version").
test(payload_not_json_typed) :-
    bt_dispatch_payload(arrange, 'this is not json', Env),
    bt_error(Env, "validation").
test(payload_not_object_typed) :-
    bt_dispatch_payload(arrange, '[1,2,3]', Env),
    bt_error(Env, "validation").

% The C55 pin: duplicate JSON keys (top-level or nested) are a client-input
% defect -> typed validation, with wire wording — never `internal` with the
% leaked dict_create/3 context the audit probe caught (p08).
test(duplicate_json_key_typed_validation) :-
    forall(member(Payload,
        [ '{"v":1,"v":2,"id":"t-d","verb":"arrange","params":{}}',
          '{"v":1,"id":"t-d","verb":"arrange","params":{"size":5,"size":7}}'
        ]),
        ( bt_dispatch_payload(arrange, Payload, Env),
          bt_error(Env, "validation"),
          get_dict(error, Env, E),
          get_dict(message, E, Msg),
          once(sub_string(Msg, _, _, _, "duplicate key")),
          \+ sub_string(Msg, _, _, _, "dict_create") )).

% --- protocol leniencies (the C59 pins) ---------------------------------------
% Four DELIBERATE lenient-reader behaviours (browser.pl parse_request/2's
% leniency note): the facade always sends a full well-formed envelope, so none
% can mask a facade bug — but the accepted shapes are contract, so lock them.

% (a) absent `v` is accepted and treated as v1.
test(lenient_absent_v_dispatches_as_v1) :-
    bt_toy_params(P),
    bt_dispatch(arrange, _{id: "t-v", verb: "arrange", params: P}, Env),
    bt_success_result(Env, _),
    get_dict(v, Env, 1).

% (b) trailing bytes after the one JSON value are ignored.
test(lenient_trailing_garbage_accepted) :-
    bt_dispatch_payload(arrange,
        '{"v":1,"id":"t-g","verb":"arrange","params":{"clues":[{"answer":"CAT"}],"size":5}} trailing garbage',
        Env),
    bt_success_result(Env, _).

% (c) `v` must be the JSON *integer* 1: 1.0 is unsupported_version
% (unreachable from a JS facade — JSON.stringify writes 1 — but hand-rolled
% clients can send it).
test(lenient_v_float_rejected) :-
    bt_dispatch_payload(arrange,
        '{"v":1.0,"id":"t-f","verb":"arrange","params":{"clues":[{"answer":"CAT"}],"size":5}}',
        Env),
    bt_error(Env, "unsupported_version").

% (d) the payload `verb` field is documentation only: the bound Verb argument
% dispatches and is echoed even when the field disagrees.
test(lenient_verb_field_not_cross_checked) :-
    bt_toy_params(P),
    bt_dispatch(arrange, _{v: 1, id: "t-m", verb: "lint", params: P}, Env),
    get_dict(verb, Env, "arrange"),
    bt_success_result(Env, R),
    get_dict(gridLength, R, 5).

% The C60 pin: an unbound Verb ARGUMENT must not silently dispatch (pre-fix it
% unified with the first verb clause and ran arrange). unknown_verb envelope
% with a stable '<unbound>' echo in the message; the envelope's verb field is
% JSON null (a variable's print name would be unstable bytes) and the caller's
% variable stays unbound.
test(unbound_verb_rejected_not_dispatched) :-
    bt_toy_params(P),
    atom_json_dict(Payload, _{v: 1, id: "t-u", verb: "arrange", params: P}, []),
    browser_dispatch(UnboundVerb, Payload, Json),
    var(UnboundVerb),
    atom_json_dict(Json, Env, [default_tag(json)]),
    bt_error(Env, "unknown_verb"),
    get_dict(verb, Env, null),
    get_dict(id, Env, "t-u"),
    get_dict(error, Env, E),
    get_dict(message, E, Msg),
    once(sub_string(Msg, _, _, _, "unbound")).

% Bound non-atom Verb arguments (string/number/compound — a facade defect,
% probe p25) fall to the unknown_verb catch-all: the flattened term is
% echoed, nothing dispatches. Pinned alongside C60's unbound case.
test(non_atom_verb_shapes_unknown_verb) :-
    forall(member(V, ["arrange", 42, f(x)]),
           ( browser_dispatch(V,
                 '{"v":1,"id":"t-nv","verb":"arrange","params":{}}', Json),
             atom_json_dict(Json, Env, [default_tag(json)]),
             bt_error(Env, "unknown_verb"),
             get_dict(id, Env, "t-nv") )).

% Non-text payloads (number/compound/unbound — Worker-impossible, the
% defensive path): the engine-level parse throws (open_string type error /
% instantiation_error) are honestly labelled `internal` — but still a tagged
% envelope, never a throw out of the spine.
test(non_text_payload_internal_envelope) :-
    forall(member(Payload, [42, f(x)]),
           ( browser_dispatch(arrange, Payload, Json),
             atom_json_dict(Json, Env, [default_tag(json)]),
             bt_error(Env, "internal") )),
    browser_dispatch(arrange, _UnboundPayload, Json2),
    atom_json_dict(Json2, Env2, [default_tag(json)]),
    bt_error(Env2, "internal").

% --- edge-case matrix (asserting the ENVELOPE, not "contains placed") --------

test(edge_empty_clues_success_empty_layout) :-
    bt_arrange(_{clues: [], size: 5}, Env),
    bt_success_result(Env, R),
    get_dict(words, R, []),
    get_dict(gridLength, R, 5).
test(edge_single_word) :-
    bt_arrange(_{clues: [_{answer: "CAT"}], size: 5}, Env),
    bt_success_result(Env, R),
    get_dict(words, R, [W]),
    get_dict(answer, W, "CAT").
test(edge_all_unplaceable_failure_with_words) :-
    bt_arrange(_{clues: [_{answer: "AAA"}, _{answer: "BBB"}], size: 5}, Env),
    bt_failure(Env, "unplaceable_words"),
    get_dict(detail, Env, D),
    get_dict(words, D, Ws),
    msort(Ws, ["AAA", "BBB"]).
test(edge_word_longer_than_grid_strict) :-
    bt_arrange(_{clues: [_{answer: "ABCDE"}], size: 3}, Env),
    bt_failure(Env, "unplaceable_words").
test(edge_nothing_placeable_best_effort) :-
    bt_arrange(_{clues: [_{answer: "ABCDE"}], size: 3, bestEffort: true}, Env),
    bt_failure(Env, "grid_too_small"),
    % grid_too_small carries no culprit list — `words` only rides
    % unplaceable_words (the failure.detail over-promise fix, strategy §4.1)
    get_dict(detail, Env, D),
    \+ get_dict(words, D, _).
test(edge_unicode_lowercase_answers) :-
    bt_arrange(_{clues: [_{answer: "café"}, _{answer: "férale"}], size: 6}, Env),
    bt_success_result(Env, R),
    get_dict(words, R, Ws),
    maplist([W]>>get_dict(answer, W, _), Ws),
    length(Ws, 2).
% The C62 pin: an empty-string answer is accepted by validation (core's
% doc_to_words domain, identical CLI-side — browser parity holds) and rides
% the wire as failure/unplaceable_words with words:[""]. Pinned AS-IS: kinder
% handling (a json_invalid_answer for empty/whitespace) is a core-lane
% decision, not a browser fix.
test(edge_empty_string_answer_current_shape) :-
    bt_arrange(_{clues: [_{answer: ""}], size: 5}, Env),
    bt_failure(Env, "unplaceable_words"),
    get_dict(detail, Env, D),
    get_dict(words, D, [""]).

test(edge_best_effort_drops_ride_diagnostics) :-
    bt_arrange(_{clues: [_{answer: "CAT"}, _{answer: "CAR"}, _{answer: "ZZZ"}],
                 size: 5, bestEffort: true}, Env),
    bt_success_result(Env, R),
    get_dict(diagnostics, R, Diag),
    get_dict(arrange, Diag, A),
    get_dict(dropped, A, ["ZZZ"]).

% --- same-instance determinism (the reset lock, strategy §3.2) ----------------

test(determinism_seedless_twice_identical) :-
    bt_toy_params(P),
    bt_arrange(P, Env1), bt_success_result(Env1, R1),
    bt_arrange(P, Env2), bt_success_result(Env2, R2),
    R1 == R2.

test(determinism_seed_42_twice_identical) :-
    bt_toy_params(P0), put_dict(_{seed: 42}, P0, P),
    bt_arrange(P, Env1), bt_success_result(Env1, R1),
    bt_arrange(P, Env2), bt_success_result(Env2, R2),
    R1 == R2.

% The leak lock: a seeded request must not perturb the NEXT seedless one.
test(determinism_seedless_after_seeded_unchanged) :-
    bt_toy_params(P0),
    bt_arrange(P0, EnvA), bt_success_result(EnvA, Baseline),
    put_dict(_{seed: 42}, P0, PSeed),
    bt_arrange(PSeed, _),
    bt_arrange(P0, EnvB), bt_success_result(EnvB, After),
    Baseline == After.

% Seed provenance rides diagnostics (json-output-spec §6.4), and seedless
% output must carry none.
test(seed_provenance_in_diagnostics) :-
    bt_toy_params(P0), put_dict(_{seed: 42}, P0, P),
    bt_arrange(P, Env), bt_success_result(Env, R),
    get_dict(seed, R.diagnostics.arrange, 42).
test(seedless_has_no_seed_provenance) :-
    bt_toy_params(P),
    bt_arrange(P, Env), bt_success_result(Env, R),
    \+ get_dict(seed, R.diagnostics.arrange, _).

% White-box: dispatch entry unconditionally clears all FOUR instance globals
% ({engine:true} does NOT — strategy §3.2's table). prng_state/1 is the
% fourth (C61): set_search_seed(7) in the setup asserts it alongside
% search_seed/1, and reset's set_search_seed(-1) must clear it transitively.
test(reset_clears_all_four_globals, [
        setup(( set_search_seed(7),
                set_check_target(2),
                set_verbose(true) )),
        cleanup(( set_search_seed(-1),
                  set_check_target(-1),
                  set_verbose(false) ))]) :-
    once(crosswordsmith_core:prng_state(_)),   % the 4th fact exists pre-reset
    crosswordsmith_browser:reset_request_state,
    \+ crosswordsmith_core:search_seed(_),
    \+ crosswordsmith_core:prng_state(_),
    \+ crosswordsmith_arrange:check_target_override(_),
    \+ crosswordsmith_core:verbose_mode.

% White-box: a seeded dispatch leaves its seed installed only until the next
% dispatch's reset (the fact itself gates the perturbed path).
test(dispatch_resets_previous_seed) :-
    bt_toy_params(P0), put_dict(_{seed: 42}, P0, PSeed),
    bt_arrange(PSeed, _),
    once(crosswordsmith_core:search_seed(42)),
    bt_arrange(P0, _),
    \+ crosswordsmith_core:search_seed(_).

% The C48 bounded-growth lock: fresh-vocab best-effort dispatches in THIS
% single engine must not accumulate memo tables. core's reset_search_memos/0
% seam (run at every search entry) bounds the store at ONE request's residue,
% flushed by the next request's search entry - pre-fix the greedy path never
% abolished and this exact scenario grew monotonically (98 KB -> 1.09 MB ->
% 5.59 MB over 1/10/50 dispatches, ~112 KB per 12-word request; audit C48).
% Post-fix, the census after dispatch N is single-request scale for every N.
% The factor-2 margin absorbs vocab-to-vocab residue variation (measured
% 77-94 variants across 50 distinct sets); pre-fix, five dispatches sit at
% ~5x and blow through it.
test(table_growth_bounded_across_dispatches) :-
    bt_fresh_vocab_params(1, P1),
    bt_dispatch(arrange, _{v: 1, id: "m-1", verb: "arrange", params: P1},
                Env1),
    bt_success_result(Env1, _),
    bt_memo_census(V1, B1),
    V1 > 0,
    forall(between(2, 5, I),
           ( bt_fresh_vocab_params(I, P),
             bt_dispatch(arrange,
                         _{v: 1, id: "m-n", verb: "arrange", params: P},
                         Env),
             bt_success_result(Env, _) )),
    bt_memo_census(V5, B5),
    assertion(V5 =< 2 * V1),
    assertion(B5 =< 2 * B1).

% The spine det lock (§10: exactly one solution, no choicepoint left — the
% main audit added deterministic/1 probes elsewhere; this is the browser
% spine's), across all three envelope arms: success, error, failure.
test(dispatch_det_lock_all_arms) :-
    bt_toy_params(P),
    atom_json_dict(Happy, _{v: 1, id: "d-1", verb: "arrange", params: P}, []),
    browser_dispatch(arrange, Happy, _),
    deterministic(D1), D1 == true,
    browser_dispatch(arrange, 'not json', _),
    deterministic(D2), D2 == true,
    atom_json_dict(Failing,
                   _{v: 1, id: "d-2", verb: "arrange",
                     params: _{clues: [_{answer: "AAA"}, _{answer: "BBB"}],
                               size: 5}}, []),
    browser_dispatch(arrange, Failing, _),
    deterministic(D3), D3 == true.

% E2E resource_exhausted + post-blowout worker health (the census probes,
% pinned): a genuine stack blowout inside an IN-CAP dispatch rides the wire
% as the typed resource_exhausted envelope — single-line message (C53's
% first-line rule) — and the very NEXT dispatch is byte-identical to its
% baseline: nothing about the blowout leaks forward. The blowout is forced
% by shrinking stack_limit around one dispatch (post-C56-cap a
% user-reachable blowout is too slow for a unit test; 0.5 MB fails the
% size-100 grid allocation in ~1 ms).
test(e2e_resource_exhausted_and_worker_health) :-
    bt_toy_params(P),
    atom_json_dict(Toy, _{v: 1, id: "rh", verb: "arrange", params: P}, []),
    browser_dispatch(arrange, Toy, Baseline),
    current_prolog_flag(stack_limit, Old),
    setup_call_cleanup(
        set_prolog_flag(stack_limit, 500_000),
        browser_dispatch(arrange,
            '{"v":1,"id":"rr","verb":"arrange","params":{"clues":[{"answer":"CAT"}],"size":100}}',
            Json),
        set_prolog_flag(stack_limit, Old)),
    atom_json_dict(Json, Env, [default_tag(json)]),
    bt_error(Env, "resource_exhausted"),
    get_dict(id, Env, "rr"),
    get_dict(error, Env, E),
    get_dict(message, E, Msg),
    once(sub_string(Msg, _, _, _, "Stack limit")),
    \+ sub_string(Msg, _, _, _, "\n"),
    browser_dispatch(arrange, Toy, After),
    Baseline == After.

% --- native value-parity vs the committed CLI goldens (DEC-8) ----------------
% The nested `result` is a live dict: correctness is value-equality of the
% parsed structures, never a byte diff (nesting re-indents / retabs). These
% prove the browser path IS the arrange verb — same fixtures, same goldens the
% CLI regression uses, fixed AND the max-crop framing.

test(value_parity_bundled_17_fixed) :-
    bt_fixture_params('fixtures/bundled_17_clues.pl', _{size: 17}, Params),
    bt_arrange(Params, Env),
    bt_success_result(Env, R),
    bt_golden('tests/golden/arrange_bundled_17_fixed.json', Golden),
    R == Golden.

test(value_parity_toc_demo_max_crop) :-
    bt_fixture_params('fixtures/toc_demo.pl', _{size: 25, mode: "max"}, Params),
    bt_arrange(Params, Env),
    bt_success_result(Env, R),
    bt_golden('tests/golden/arrange_toc_demo_max.json', Golden),
    R == Golden.

% --- resolver parity vs the CLI binary (strategy §7 / §10) --------------------
% The slice duplicates the CLI's option resolvers (throwing variants); this
% matrix locks the duplication: whatever the CLI rejects at the flag layer,
% the browser rejects as a typed validation error — and the CLI's accepted
% domain stays accepted. (JSON has true absence, so the CLI's -1 sentinel is
% deliberately NOT part of the browser domain — seed:-1 is bad_seed here.)

bt_parity_reject(seed_neg,   ['--seed', '-5'],           _{seed: -5}).
bt_parity_reject(seed_float, ['--seed', '3.5'],          _{seed: 3.5}).
bt_parity_reject(size_zero,  ['--size', '0'],            _{size: 0}).
bt_parity_reject(size_neg,   ['--size', '-3'],           _{size: -3}).
bt_parity_reject(seed_be,    ['--seed', '7', '--best-effort'],
                                                          _{seed: 7, bestEffort: true}).

test(cli_parity_rejects, [forall(bt_parity_reject(_Name, CliArgs, Mut))]) :-
    bt_with_toy_input([File]>>(
        ( memberchk('--size', CliArgs) -> SizeArgs = [] ; SizeArgs = ['--size', '5'] ),
        append([[arrange, '--input', File], SizeArgs, CliArgs], Args),
        bt_cli_exit(Args, Code),
        Code =\= 0
    )),
    bt_toy_params(P0), put_dict(Mut, P0, P),
    bt_arrange(P, Env),
    bt_error(Env, "validation").

bt_parity_accept(seed_zero, ['--seed', '0'], _{seed: 0}).
bt_parity_accept(seed_big,  ['--seed', '123456789012345678901234567890'],
                            _{seed: 123456789012345678901234567890}).

test(cli_parity_accepts, [forall(bt_parity_accept(_Name, CliArgs, Mut))]) :-
    bt_with_toy_input([File]>>(
        append([arrange, '--input', File, '--size', '5'], CliArgs, Args),
        bt_cli_exit(Args, Code),
        Code =:= 0
    )),
    bt_toy_params(P0), put_dict(Mut, P0, P),
    bt_arrange(P, Env),
    bt_success_result(Env, _).

% --- the lint + export verbs (strategy §6 phases 2-3) --------------------------
% Pure layout transforms on the spine: no engine, no seed, success-only happy
% paths (a FAIL verdict is still a successful lint — the CLI's verdict-as-exit-
% code is a shell convention, not part of the answer).

% Value parity: the browser report deep-equals the committed CLI golden
% (generated from the same arrange golden by `make regen-goldens`).
test(lint_value_parity_toc_golden) :-
    bt_layout(Layout),
    bt_lint(_{layout: Layout, profile: "toc"}, Env),
    get_dict(verb, Env, "lint"),
    bt_success_result(Env, R),
    bt_golden('tests/golden/lint_bundled_17_toc.json', Golden),
    R == Golden.

% ipuz is itself a JSON document, so the result IS the ipuz object.
test(export_ipuz_value_parity_golden) :-
    bt_layout(Layout),
    bt_export(_{layout: Layout, to: "ipuz"}, Env),
    bt_success_result(Env, R),
    bt_golden('tests/golden/export_bundled_17.ipuz', Golden),
    R == Golden.

% Exolve is plain text -> {format:"text", body}; body byte-equals the CLI's
% stdout (the golden file).
test(export_exolve_value_parity_golden) :-
    bt_layout(Layout),
    bt_export(_{layout: Layout, to: "exolve"}, Env),
    bt_success_result(Env, R),
    get_dict(format, R, "text"),
    get_dict(body, R, Body),
    read_file_to_string('tests/golden/export_bundled_17.exolve', GoldenText, []),
    Body == GoldenText.

% allowAsymmetry is honoured (echoed in the report + the report differs from
% the strict-symmetry golden, whose symmetry rule WARNs on this layout).
test(lint_allow_asymmetry_honoured) :-
    bt_layout(Layout),
    bt_lint(_{layout: Layout, profile: "toc", allowAsymmetry: true}, Env),
    bt_success_result(Env, R),
    get_dict(allowAsymmetry, R, true),
    bt_golden('tests/golden/lint_bundled_17_toc.json', Golden),
    R \== Golden.

% The composition the SDK exists for: arrange -> lint the result -> export it,
% all in one instance, each result feeding the next verb verbatim.
test(arrange_lint_export_composition) :-
    bt_toy_params(P),
    bt_arrange(P, AEnv),
    bt_success_result(AEnv, Layout),
    bt_lint(_{layout: Layout, profile: "american"}, LEnv),
    bt_success_result(LEnv, Report),
    get_dict(verdict, Report, V),
    memberchk(V, ["PASS", "WARN", "FAIL"]),
    get_dict(summary, Report, _),
    bt_export(_{layout: Layout, to: "exolve"}, EEnv),
    bt_success_result(EEnv, X),
    get_dict(body, X, Body),
    once(sub_string(Body, 0, _, _, "exolve-begin")).

% Malformed lint/export params -> typed validation envelopes, never bare-fail.
test(lint_export_param_fuzz) :-
    bt_layout(Layout),
    forall(member(V-P,
        [ lint-_{profile: "toc"},                          % layout absent
          lint-_{layout: 42, profile: "toc"},              % layout not an object
          lint-_{layout: Layout},                          % profile absent
          lint-_{layout: Layout, profile: "bogus"},        % unknown profile
          lint-_{layout: Layout, profile: 7},              % profile wrong type
          lint-_{layout: Layout, profile: "toc", allowAsymmetry: "yes"},
          export-_{to: "ipuz"},                            % layout absent
          export-_{layout: Layout},                        % to absent
          export-_{layout: Layout, to: "pdf"},             % unknown format
          export-_{layout: _{}, to: "ipuz"},               % fails the shape gate
          % the two C50 probe payloads: gate-passing garbage used to reach the
          % transforms (non-list grid row -> mislabelled internal "engine bug";
          % schema-less words entry -> degenerate SUCCESS ipuz), now rejected
          % by the deepened shared gate as typed validation errors
          export-_{layout: _{gridLength: 1, grid: [5], words: []}, to: "ipuz"},
          export-_{layout: _{gridLength: 1, grid: [5], words: []}, to: "exolve"},
          export-_{layout: _{gridLength: 3, grid: [], words: [_{bogus: 1}]}, to: "ipuz"},
          % obvious same-class shapes
          export-_{layout: _{gridLength: 2, grid: [[null, null], "xy"],
                             words: []}, to: "ipuz"},      % one row a string
          export-_{layout: _{gridLength: 3, grid: [], words: [7]}, to: "ipuz"},
          export-_{layout: _{gridLength: 3, grid: [],    % word: cells not a list
                             words: [_{answer: "CAT", direction: "across",
                                       cells: "garbage"}]}, to: "exolve"}
        ]),
        ( bt_dispatch(V, _{v: 1, id: "f-le", verb: V, params: P}, Env),
          bt_error(Env, "validation") )).

% The C50 regression locks, asserting the exact pre-fix mislabels can't return:
% the non-list-grid-row payload must NOT map to internal (the never-plain-fail
% backstop's "engine bug" wording), and the schema-less-word payload must NOT
% be blessed as success (the degenerate empty ipuz).
test(export_gate_nonlist_grid_row_not_internal) :-
    forall(member(To, ["ipuz", "exolve"]),
           ( bt_export(_{layout: _{gridLength: 1, grid: [5], words: []},
                         to: To}, Env),
             bt_error(Env, "validation"),
             get_dict(error, Env, E),
             get_dict(message, E, Msg),
             \+ sub_string(Msg, _, _, _, "engine bug") )).
test(export_gate_schemaless_word_not_success) :-
    bt_export(_{layout: _{gridLength: 3, grid: [], words: [_{bogus: 1}]},
                to: "ipuz"}, Env),
    \+ get_dict(status, Env, "success"),
    bt_error(Env, "validation").

% lint's own layout gate (lint_dict_layout/3's lint_* formals) maps to typed
% validation envelopes with the project's hook-rendered messages.
test(lint_layout_gate_fuzz) :-
    forall(member(L,
        [ _{},                                                     % no gridLength
          _{gridLength: 5},                                        % no words array
          _{gridLength: 5, words: [_{direction: "across"}]},       % word: no answer
          _{gridLength: 5, words: [_{answer: "CAT"}]},             % word: no direction
          _{gridLength: 5, words: [_{answer: "CAT", direction: "diagonal", cells: [[0, 0]]}]},
          _{gridLength: 5, words: [_{answer: "CAT", direction: "across"}]},  % no cells
          _{gridLength: 5, words: [_{answer: "CAT", direction: "across",
                                     cells: [[9, 9], [0, 1], [0, 2]]}]}      % cell off-grid
        ]),
        ( bt_lint(_{layout: L, profile: "toc"}, Env),
          bt_error(Env, "validation") )).

% Resolver parity with the CLI binary: the same rejections reject there.
test(cli_parity_lint_export_rejects) :-
    forall(member(Args,
        [ [lint, '--profile', bogus, 'tests/golden/arrange_bundled_17_fixed.json'],
          [lint, 'tests/golden/arrange_bundled_17_fixed.json'],
          [export, '--to', pdf, 'tests/golden/arrange_bundled_17_fixed.json'],
          [export, 'tests/golden/arrange_bundled_17_fixed.json']
        ]),
        ( bt_cli_exit(Args, Code), Code =\= 0 )).

% --- the verb registry (the C49 lock) ------------------------------------------
% verb/1 is THE verb list: white-box, every verb/1 fact must have a dispatch/3
% clause keyed on it and every verb-keyed dispatch clause must be registered
% (clause/2 leaves the head arg unbound for the invalid-request and unknown-verb
% catch-alls, so atom/1 selects exactly the per-verb clauses). Capabilities and
% the unknown_verb message render the same registry, so with this bijection
% locked none of the three can drift.
test(verb_registry_dispatch_bijection) :-
    findall(V, crosswordsmith_browser:verb(V), Registry),
    Registry \== [],
    is_set(Registry),
    findall(V,
            ( clause(crosswordsmith_browser:dispatch(V, _, _), _),
              atom(V) ),
            DispatchVerbs),
    msort(Registry, Sorted),
    msort(DispatchVerbs, Sorted).

% The wire copies source the registry: capabilities returns exactly the
% registered verbs (in registry order) and the unknown-verb diagnostic renders
% every one of them.
test(verb_registry_sources_capabilities_and_message) :-
    findall(V, crosswordsmith_browser:verb(V), Registry),
    maplist([A, S]>>atom_string(A, S), Registry, RegistryStrings),
    bt_dispatch(capabilities, _{v: 1, id: "c-r", verb: "capabilities"}, CEnv),
    bt_success_result(CEnv, R),
    get_dict(verbs, R, RegistryStrings),
    bt_dispatch(nope, _{v: 1, id: "u-r", verb: "nope", params: _{}}, UEnv),
    bt_error(UEnv, "unknown_verb"),
    get_dict(error, UEnv, E),
    get_dict(message, E, Msg),
    forall(member(VS, RegistryStrings), sub_string(Msg, _, _, _, VS)).

% --- the capabilities verb (engine-sourced, OQ-4) -----------------------------

test(capabilities_engine_sourced) :-
    bt_dispatch(capabilities, _{v: 1, id: "c-1", verb: "capabilities"}, Env),
    bt_success_result(Env, R),
    get_dict(verbs, R, Verbs),
    memberchk("arrange", Verbs),
    memberchk("lint", Verbs),
    memberchk("export", Verbs),
    get_dict(engine, R, E),
    get_dict(swipl, E, Swipl),
    string(Swipl).

% dispatch clause-1 ordering: an unparseable payload beats even the
% payload-less capabilities verb — the parse-error envelope answers (verb
% still echoed, id null), never a capabilities success built from garbage.
test(capabilities_with_unparseable_payload_is_parse_error) :-
    bt_dispatch_payload(capabilities, 'not json', Env),
    bt_error(Env, "validation"),
    get_dict(verb, Env, "capabilities"),
    get_dict(id, Env, null).

% --- classifier totality (white-box) ------------------------------------------
% outcome_envelope/4 must map EVERY term: the outcome atoms, mapped throws,
% unmapped throws -> internal, and the never-plain-fail backstop.

test(classifier_not_proven_is_budget_exceeded) :-
    crosswordsmith_browser:outcome_envelope(arrange, "t", not_proven, Env),
    Env.status == error,
    Env.error.type == budget_exceeded.
test(classifier_resource_error_mapped) :-
    crosswordsmith_browser:outcome_envelope(
        arrange, "t", thrown(error(resource_error(stack), _)), Env),
    Env.status == error,
    Env.error.type == resource_exhausted.
test(classifier_unmapped_throw_is_internal) :-
    crosswordsmith_browser:outcome_envelope(
        arrange, "t", thrown(error(weird_never_seen(x), _)), Env),
    Env.status == error,
    Env.error.type == internal.
test(classifier_bare_throw_is_internal) :-
    crosswordsmith_browser:outcome_envelope(arrange, "t", thrown(kaboom), Env),
    Env.status == error,
    Env.error.type == internal.
test(classifier_failed_dispatch_is_internal) :-
    crosswordsmith_browser:outcome_envelope(arrange, "t", failed(arrange), Env),
    Env.status == error,
    Env.error.type == internal.
test(classifier_unknown_outcome_is_internal) :-
    crosswordsmith_browser:outcome_envelope(arrange, "t", frobnicated(9), Env),
    Env.status == error,
    Env.error.type == internal.

% The C54 lock: even a SERIALIZATION failure yields a tagged envelope.
% White-box at the emit seam, forcing the impossible-today failures the
% armour exists for — a success outcome whose result dict smuggles (a) a
% non-JSON compound and (b) an unbound variable past the total classifier.
% Both must come back as the minimal internal envelope with id and verb
% still echoed — never a raw throw out of the spine's last step.
test(emit_serialization_failure_yields_internal_envelope) :-
    forall(member(Bad, [ _{x: foo(bar)}, _{x: _} ]),
           ( crosswordsmith_browser:emit_envelope(arrange, "t-54",
                                                  placed(Bad), Json),
             atom(Json),
             atom_json_dict(Json, Env, [default_tag(json)]),
             bt_error(Env, "internal"),
             get_dict(id, Env, "t-54"),
             get_dict(verb, Env, "arrange"),
             get_dict(error, Env, E),
             get_dict(message, E, Msg),
             once(sub_string(Msg, _, _, _, "serialization failed")) )).

% --- the arrange_outcome/5 seam (the exported wrapper, DEC-7) -----------------

test(arrange_outcome_placed_tag) :-
    bt_toy_clues(Clues),
    findall([A, _{}],
            ( member(D, Clues), get_dict(answer, D, S), atom_string(A, S) ),
            Words),
    arrange_outcome(Words, 5, strict, fixed, placed(Dict)),
    get_dict(gridLength, Dict, 5).
test(arrange_outcome_infeasible_tag) :-
    arrange_outcome([['AAA', _{}], ['BBB', _{}]], 5, strict, fixed,
                    infeasible(Bad)),
    msort(Bad, ['AAA', 'BBB']).
test(arrange_outcome_nothing_placeable_tag) :-
    arrange_outcome([['ABCDE', _{}]], 3, best_effort, fixed, nothing_placeable).
test(arrange_outcome_duplicate_throws, [throws(error(duplicate_answer('CAT'), _))]) :-
    arrange_outcome([['CAT', _{}], ['CAT', _{}]], 5, strict, fixed, _).

:- end_tests(browser).
