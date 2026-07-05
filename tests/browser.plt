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
:- use_module(library(http/json)).
:- use_module(library(process)).

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

% A huge seed that is still a JSON *integer* stays in the CLI's accepted
% domain (integer >= 0 — SWI reads bignums); only a float-degraded one is bad.
test(huge_integer_seed_accepted) :-
    bt_dispatch_payload(arrange,
        '{"v":1,"id":"t-1","verb":"arrange","params":{"clues":[{"answer":"CAT"},{"answer":"CAR"},{"answer":"ARC"},{"answer":"RAT"},{"answer":"TAR"}],"size":5,"seed":123456789012345678901234567890}}',
        Env),
    bt_success_result(Env, _).

% --- protocol errors ---------------------------------------------------------

test(unknown_verb_typed) :-
    bt_dispatch(lint, _{v: 1, id: "t-1", verb: "lint", params: _{}}, Env),
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

% White-box: dispatch entry unconditionally clears all three instance globals
% ({engine:true} does NOT — strategy §3.2's table).
test(reset_clears_all_three_globals, [
        setup(( set_search_seed(7),
                set_check_target(2),
                set_verbose(true) )),
        cleanup(( set_search_seed(-1),
                  set_check_target(-1),
                  set_verbose(false) ))]) :-
    crosswordsmith_browser:reset_request_state,
    \+ crosswordsmith_core:search_seed(_),
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

% --- the capabilities verb (engine-sourced, OQ-4) -----------------------------

test(capabilities_engine_sourced) :-
    bt_dispatch(capabilities, _{v: 1, id: "c-1", verb: "capabilities"}, Env),
    bt_success_result(Env, R),
    get_dict(verbs, R, Verbs),
    memberchk("arrange", Verbs),
    get_dict(engine, R, E),
    get_dict(swipl, E, Swipl),
    string(Swipl).

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
