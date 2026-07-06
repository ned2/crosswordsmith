% browser.pl - the browser/WASM dispatch spine (wasm-sdk-strategy §3/§4).
%
% One generic entry point, browser_dispatch/3: JSON-string request in, JSON
% envelope (a Prolog atom) out. The Worker binds Verb and the payload as
% forEach variables (never string-concatenated) and JSON is used in BOTH
% directions because SWI's native WASM data translation is lossy for our
% schemas (JS strings arrive as atoms; Prolog null/true/false return as the JS
% strings "null"/"true"/"false" — see wasm/client/solve_browser.pl's history).
%
% Wire contract (strategy §4.1) — every response is one of:
%   {v:1, id, verb, status:"success", result: {...}}   a live layout dict
%   {v:1, id, verb, status:"failure", detail: {reason, words?}}
%   {v:1, id, verb, status:"error",   error: {type, message}}
%
% The cross-cutting invariant (strategy §10): every dispatch path must unify a
% TAGGED outcome or throw a term the classifier MAPS — never plain-fail. A bare
% fail is what produced the spike's generic "no layout (search failed)"
% masquerade; browser_dispatch/3 therefore also backstops a failing dispatch/3
% with an `internal` error envelope rather than failing itself. In particular
% the option resolvers below validate-and-THROW (bad_size/bad_mode/bad_seed…)
% instead of leaning on a builtin failing (set_search_seed/1 FAILS on a bad
% seed — strategy §4.2's correction).
%
% Adding a verb = one verb/1 fact + one dispatch/3 clause (+ classifier rows if
% it introduces new outcome tags) — never new message wiring: the capabilities
% result and the unknown_verb diagnostic both read the verb/1 registry, so
% neither can drift from what actually dispatches. `arrange` is the first verb;
% `capabilities` rides along so the JS facade's capabilities() is
% engine-sourced, not a JS hardcode that lies after a qlf swap (OQ-4); `lint`
% and `export` (strategy §6 phases 2–3) prove the spine generalises: both are
% pure layout transforms — no engine, no seed, success-only happy paths.

:- module(crosswordsmith_browser,
          [ browser_dispatch/3
          ]).

% library(json), NOT the legacy library(http/json) alias: this module is
% qcompiled into the WASM image, where the alias does not resolve at all
% ("source_sink library(http/json) does not exist" load noise — wasm/README.md)
% and the predicates only worked by autoloading from core library(json).
% Explicit import list per the project imports below (qsave autoload(false)
% discipline).
:- use_module(library(json),
              [ atom_json_dict/3
              ]).
:- use_module(crosswordsmith(core),
              [ doc_to_words/2,
                set_verbose/1,
                set_search_seed/1
              ]).
:- use_module(crosswordsmith(arrange),
              [ arrange_outcome/5,
                set_check_target/1
              ]).
:- use_module(crosswordsmith(lint),
              [ lint_dict_layout/3,
                lint_run/5,
                lint_known_profile/1
              ]).
:- use_module(crosswordsmith(export),
              [ export_layout_dict/2,
                layout_to_ipuz/2,
                layout_to_exolve/2
              ]).

% browser_dispatch(+Verb:atom, +Payload, -JsonEnvelope:atom)
% Payload is the request envelope {v, id, verb, params, meta?} as a JSON
% string/atom (what JSON.stringify hands the Worker). JsonEnvelope is returned
% as an ATOM: an atom crosses the WASM boundary as a clean JS String, whereas a
% Prolog string arrives wrapped as a Prolog.String instance (strategy §4.3).
browser_dispatch(Verb, Payload, JsonEnvelope) :-
    reset_request_state,
    catch(parse_request(Payload, Req), ParseErr, Req = invalid(ParseErr)),
    request_id(Req, Id),
    catch(( dispatch(Verb, Req, Outcome0)
          ->  Outcome = Outcome0
          ;   Outcome = failed(Verb)          % the never-plain-fail backstop
          ),
          Err,
          Outcome = thrown(Err)),
    envelope_verb(Verb, VerbA),
    outcome_envelope(VerbA, Id, Outcome, Envelope),
    % width(0): a single-line wire envelope — nothing downstream reads the
    % bytes (JSON.parse both ways, DEC-8), so pretty-printing only costs
    % transport size and breaks line-oriented harnesses. as(atom) makes the
    % §4.3 atom-not-string contract explicit (it is also the default).
    atom_json_dict(JsonEnvelope, Envelope, [width(0), as(atom)]).

% --- per-request state reset (strategy §3.2) ---------------------------------
% Unconditionally clear the mutable module globals before every dispatch: a
% Worker reuses one Prolog instance across requests and the per-solve
% {engine:true} engine does NOT clear module dynamics, so a previous request's
% seed/checkTarget/verbose would otherwise leak forward. All three setters run
% defensively even though the slice only exposes `seed`; set_search_seed(-1)
% also clears core's prng_state/1 (the fourth dynamic fact) transitively.
%
% Tabled memos (core's pair_crossings/3 + answer_letters/2) are deliberately
% NOT abolished here: every arrange search entry runs core's
% reset_search_memos/0 seam (C48), so a single-engine embedding looping
% browser_dispatch/3 holds at most ONE request's table residue (~100 KB for a
% 12-word request), flushed by the next request's search entry - growth is
% bounded, no longer one memo variant per word pair ever seen. An abolish here
% would be redundant with that seam, and because this reset runs on ENTRY (not
% exit) it could not shrink the between-request footprint either. The bound is
% locked by browser.plt's bounded-growth test.
reset_request_state :-
    set_search_seed(-1),
    set_check_target(-1),
    set_verbose(false).

% --- request parsing ---------------------------------------------------------

% Parse the request JSON. atom_json_dict lands every value with the right
% Prolog type (strings stay strings, so core's `answer` guard works; numbers,
% null and booleans survive; non-object JSON parses to a non-dict, so the
% is_dict gate below still fires). Throws syntax_error(_) on malformed JSON,
% bad_request(_) on a non-object, unsupported_version(_) on a v we don't speak
% — all mapped by the classifier, so a garbage payload still gets a typed
% envelope.
parse_request(Payload, Req) :-
    atom_json_dict(Payload, Req0, [default_tag(json)]),
    (   is_dict(Req0)
    ->  Req = Req0
    ;   throw(error(bad_request(Req0), _))
    ),
    (   get_dict(v, Req, V), V \== 1
    ->  throw(error(unsupported_version(V), _))
    ;   true
    ).

% The response echoes the request's `id` so the facade can route each reply to
% its pending promise (concurrent same-verb calls — strategy §4.1). No usable
% id (unparseable request, absent/non-atomic id) => JSON null.
request_id(invalid(_), null) :- !.
request_id(Req, Id) :-
    (   get_dict(id, Req, Id0), atomic(Id0)
    ->  Id = Id0
    ;   Id = null
    ).

% The verb is echoed for self-describing logs. It arrives as an atom from the
% Worker's bound forEach variable; anything else (a defensive impossibility) is
% flattened so json_write_dict cannot choke on the envelope itself.
envelope_verb(Verb, Verb) :- atom(Verb), !.
envelope_verb(Verb, VerbA) :- term_to_atom(Verb, VerbA).

% --- the verb registry ---------------------------------------------------------
% THE one verb list (there used to be three hand-maintained copies — dispatch
% clauses, the capabilities literal, the unknown_verb message text — a
% capabilities-that-lies hazard, exactly what OQ-4 built the verb to prevent).
% engine_capabilities/1 collects these facts and the unknown_verb hook renders
% them at message time; dispatch/3 keys one clause per fact, and browser.plt's
% registry test locks the fact/clause bijection. Source order is wire order
% (the capabilities envelope is plt-locked).
verb(arrange).
verb(lint).
verb(export).
verb(capabilities).

% --- dispatch ----------------------------------------------------------------
% dispatch(+Verb, +Req, -Outcome): one clause per verb/1 fact, unknown-verb
% catch-all LAST. Every clause unifies a tagged Outcome or throws a mapped term.
dispatch(_Verb, invalid(Err), thrown(Err)) :- !.   % request never parsed
dispatch(arrange, Req, Outcome) :- !,
    arrange_browser(Req, Outcome).
dispatch(lint, Req, lint_report(Report)) :- !,
    lint_browser(Req, Report).
dispatch(export, Req, exported(Doc)) :- !,
    export_browser(Req, Doc).
dispatch(capabilities, _Req, capabilities(Caps)) :- !,
    engine_capabilities(Caps).
dispatch(Verb, _Req, thrown(error(unknown_verb(Verb), _))).

% --- the arrange verb (strategy §7) -------------------------------------------
% Exposed params (DEC-6): size (default 15, the CLI's bare-invocation default),
% mode ("fixed"|"max", default fixed), bestEffort (default false), seed
% (absent = the deterministic search). Fragment/candidates/enumerate are out of
% the slice; checkTarget is deferred (the reset still clears it defensively).
arrange_browser(Req, Outcome) :-
    request_params(Req, Params),
    doc_to_words(Params, Words),          % throws json_no_clues_array / …
    resolve_size(Params, GridLen),
    resolve_mode(Params, SizeMode),
    resolve_best_effort(Params, Drop),
    resolve_seed(Params, Seed),
    guard_seed_drop(Seed, Drop),
    apply_seed(Seed),
    arrange_outcome(Words, GridLen, Drop, SizeMode, Outcome).

request_params(Req, Params) :-
    (   get_dict(params, Req, P)
    ->  (   is_dict(P)
        ->  Params = P
        ;   throw(error(bad_params(P), _))
        )
    ;   Params = _{}      % absent params -> doc_to_words throws (typed)
    ).

% Local validate-and-throw resolvers (strategy §3.1: the CLI's resolvers report
% via cli_error->fail, which a browser path must not reuse; the parity matrix
% in tests/browser.plt locks the accepted-for-now duplication). Each mirrors
% the CLI's accepted domain: size ~ validate_size_flag (positive integer;
% absent -> resolve_size's 15), seed ~ validate_seed (integer >= 0), the
% seed+bestEffort guard ~ guard_seed_combos. JSON has true absence, so the
% CLI's -1 "not given" sentinels are NOT part of the browser domain: seed:-1 is
% a bad seed here, not "unset".

resolve_size(Params, GridLen) :-
    (   get_dict(size, Params, S)
    ->  (   integer(S), S > 0
        ->  GridLen = S
        ;   throw(error(bad_size(S), _))
        )
    ;   GridLen = 15
    ).

% mode:"max" = build on size N, crop to the tight enclosing square — the CLI's
% --max-size framing, routed through arrange_diag_layout_dict/5's SizeMode.
% (This un-blocks the spike's browser_max_mode_unsupported: the crop path is
% now golden-locked, wasm/test/value_golden.sh.)
resolve_mode(Params, SizeMode) :-
    (   get_dict(mode, Params, M)
    ->  (   M == "fixed" -> SizeMode = fixed
        ;   M == "max"   -> SizeMode = max
        ;   throw(error(bad_mode(M), _))
        )
    ;   SizeMode = fixed
    ).

resolve_best_effort(Params, Drop) :-
    (   get_dict(bestEffort, Params, B)
    ->  (   B == true  -> Drop = best_effort
        ;   B == false -> Drop = strict
        ;   throw(error(bad_best_effort(B), _))
        )
    ;   Drop = strict
    ).

resolve_seed(Params, Seed) :-
    (   get_dict(seed, Params, N)
    ->  (   integer(N), N >= 0
        ->  Seed = seed(N)
        ;   throw(error(bad_seed(N), _))
        )
    ;   Seed = none
    ).

% Randomisation perturbs the strict search only (matches the CLI's
% guard_seed_combos): best-effort never reaches the perturbed seam.
guard_seed_drop(seed(_), best_effort) :- !,
    throw(error(seed_with_best_effort, _)).
guard_seed_drop(_, _).

% Seed already validated integer >= 0, so set_search_seed/1 cannot fail here.
apply_seed(none).
apply_seed(seed(N)) :- set_search_seed(N).

% --- the lint verb (strategy §6 phase 2) ---------------------------------------
% Params: layout (a canonical arrange result), profile (required — same domain
% as the CLI's --profile), allowAsymmetry (default false). The report dict IS
% the success result: a FAIL verdict is still a successful lint (the CLI's
% verdict-as-exit-code is a shell convention, not part of the answer).
% lint_dict_layout/3 owns the layout validation and throws the lint_* formals
% (mapped below), so a malformed layout gets a typed validation envelope.
lint_browser(Req, Report) :-
    request_params(Req, Params),
    resolve_layout(Params, Layout),
    resolve_profile(Params, Profile),
    resolve_allow_asymmetry(Params, Allow),
    lint_dict_layout(Layout, GridLen, PlacedWords),
    lint_run(PlacedWords, GridLen, Profile, Allow, Report).

% --- the export verb (strategy §6 phase 3) --------------------------------------
% Params: layout + to ("ipuz" | "exolve"). Shape validation goes through the
% SAME export_layout_dict/2 the CLI's file load uses (no drifting twin). ipuz
% is itself a JSON document -> the result verbatim; Exolve is plain text ->
% wrapped as {format:"text", body} (§6) since an envelope result must be JSON.
export_browser(Req, Doc) :-
    request_params(Req, Params),
    resolve_layout(Params, Layout0),
    resolve_export_to(Params, To),
    export_layout_dict(Layout0, Layout),
    export_transform(To, Layout, Doc).

export_transform(ipuz, Layout, Ipuz) :-
    layout_to_ipuz(Layout, Ipuz).
export_transform(exolve, Layout, _{format: text, body: Body}) :-
    layout_to_exolve(Layout, Body).

% Shared by lint + export: the layout param must be present and a JSON object;
% finer shape checks belong to the verb's own gate (lint_dict_layout/3,
% export_layout_dict/2) so the browser stays a thin adapter.
resolve_layout(Params, Layout) :-
    (   get_dict(layout, Params, L)
    ->  (   is_dict(L)
        ->  Layout = L
        ;   throw(error(bad_layout(L), _))
        )
    ;   throw(error(missing_layout, _))
    ).

resolve_profile(Params, Profile) :-
    (   get_dict(profile, Params, P0)
    ->  (   string(P0), atom_string(P, P0), lint_known_profile(P)
        ->  Profile = P
        ;   throw(error(bad_profile(P0), _))
        )
    ;   throw(error(missing_profile, _))
    ).

resolve_allow_asymmetry(Params, Allow) :-
    (   get_dict(allowAsymmetry, Params, A)
    ->  (   A == true  -> Allow = true
        ;   A == false -> Allow = false
        ;   throw(error(bad_allow_asymmetry(A), _))
        )
    ;   Allow = false
    ).

resolve_export_to(Params, To) :-
    (   get_dict(to, Params, T0)
    ->  (   T0 == "ipuz"   -> To = ipuz
        ;   T0 == "exolve" -> To = exolve
        ;   throw(error(bad_to(T0), _))
        )
    ;   throw(error(missing_to, _))
    ).

% --- the capabilities verb -----------------------------------------------------
% Engine-sourced (OQ-4): the verb list comes from the loaded engine's verb/1
% registry, so a stale cached qlf answers for itself instead of a JS constant
% (or a Prolog literal) lying about it.
engine_capabilities(_{ verbs: Verbs,
                       engine: _{ swipl: Swipl } }) :-
    findall(V, verb(V), Verbs),
    current_prolog_flag(version_data, swi(Major, Minor, Patch, _)),
    format(atom(Swipl), "~w.~w.~w", [Major, Minor, Patch]).

% --- the outcome classifier (strategy §4.2) ------------------------------------
% outcome_envelope(+Verb, +Id, +Outcome, -Envelope) is TOTAL: every tagged
% outcome maps to a status, every thrown term maps to an error.type through
% thrown_type/2 with an `internal` fallback — never crash, never mislabel.

% success: the result is a LIVE dict (DEC-8) — correctness downstream is
% value-equality after JSON.parse, never a byte diff of the nested rendering.
outcome_envelope(Verb, Id, placed(Dict), Env) :- !,
    success_envelope(Verb, Id, Dict, Env).
outcome_envelope(Verb, Id, best_effort(Dict), Env) :- !,
    success_envelope(Verb, Id, Dict, Env).
outcome_envelope(Verb, Id, capabilities(Dict), Env) :- !,
    success_envelope(Verb, Id, Dict, Env).
outcome_envelope(Verb, Id, lint_report(Dict), Env) :- !,
    success_envelope(Verb, Id, Dict, Env).
outcome_envelope(Verb, Id, exported(Dict), Env) :- !,
    success_envelope(Verb, Id, Dict, Env).
% budget_exceeded is synthesised from the not_proven ATOM — the engine already
% separates budget from infeasible (arrange_best_layout/6); there is no throw.
outcome_envelope(Verb, Id, not_proven, Env) :- !,
    error_envelope(Verb, Id, budget_exceeded,
                   'not proven within budget (search did not complete; feasibility unresolved)',
                   Env).
% failure = a legitimate "no layout" answer, not an error. `words` is present
% only when there are genuinely-uncrossable answers (the typical infeasible
% case has none — the grid is just too tight for a full interlock).
outcome_envelope(Verb, Id, infeasible(Bad), Env) :- !,
    (   Bad = [_|_]
    ->  Detail = _{reason: unplaceable_words, words: Bad}
    ;   Detail = _{reason: no_interlock}
    ),
    failure_envelope(Verb, Id, Detail, Env).
outcome_envelope(Verb, Id, nothing_placeable, Env) :- !,
    failure_envelope(Verb, Id, _{reason: grid_too_small}, Env).
outcome_envelope(Verb, Id, thrown(Err), Env) :- !,
    thrown_envelope(Verb, Id, Err, Env).
outcome_envelope(Verb, Id, failed(_), Env) :- !,
    error_envelope(Verb, Id, internal,
                   'dispatch failed without tagging an outcome (engine bug)',
                   Env).
outcome_envelope(Verb, Id, Other, Env) :-        % totality net
    term_to_atom(Other, OtherA),
    format(atom(Msg), "unrecognised outcome term: ~w", [OtherA]),
    error_envelope(Verb, Id, internal, Msg, Env).

thrown_envelope(Verb, Id, error(Formal, _), Env) :-
    thrown_type(Formal, Type),
    !,
    formal_message(Formal, Msg),
    error_envelope(Verb, Id, Type, Msg, Env).
thrown_envelope(Verb, Id, Err, Env) :-           % unmatched formal / bare throw
    term_to_atom(Err, Msg),
    error_envelope(Verb, Id, internal, Msg, Env).

% The throw -> error.type table over the REAL terms this codebase throws.
% Anything unlisted falls through to `internal` above.
thrown_type(json_no_clues_array,    validation).
thrown_type(json_invalid_answer(_), validation).
thrown_type(json_invalid_meta(_),   validation).
thrown_type(duplicate_answer(_),    validation).
thrown_type(bad_request(_),         validation).
thrown_type(bad_params(_),          validation).
thrown_type(bad_size(_),            validation).
thrown_type(bad_mode(_),            validation).
thrown_type(bad_seed(_),            validation).
thrown_type(bad_best_effort(_),     validation).
thrown_type(seed_with_best_effort,  validation).
thrown_type(bad_layout(_),          validation).
thrown_type(missing_layout,         validation).
thrown_type(bad_profile(_),         validation).
thrown_type(missing_profile,        validation).
thrown_type(bad_allow_asymmetry(_), validation).
thrown_type(bad_to(_),              validation).
thrown_type(missing_to,             validation).
% lint's own layout gate (lint.pl, hooks alongside) + export's shape gate
thrown_type(lint_no_grid_length,        validation).
thrown_type(lint_no_words_array,        validation).
thrown_type(lint_invalid_word(_),       validation).
thrown_type(lint_invalid_direction(_),  validation).
thrown_type(lint_no_cells(_),           validation).
thrown_type(lint_invalid_cell(_, _),    validation).
thrown_type(export_invalid_layout,      validation).
thrown_type(syntax_error(_),        validation).    % request JSON didn't parse
thrown_type(resource_error(_),      resource_exhausted).
thrown_type(unknown_verb(_),        unknown_verb).
thrown_type(unsupported_version(_), unsupported_version).

% Render a formal through the same prolog:error_message//1 hooks the CLI uses
% (core.pl et al. + the browser-local clauses below), so the envelope's message
% matches the native rendering; unhooked formals degrade to term_to_atom.
% syntax_error is special-cased (ISO formals are rendered by SWI's own message
% machinery, not the project hook — and only the envelope wants this wording,
% so no global hook clause that would restyle the CLI's rendering).
formal_message(syntax_error(What), Msg) :- !,
    format(atom(Msg), "request is not valid JSON (~w)", [What]).
formal_message(Formal, Msg) :-
    (   catch(phrase(prolog:error_message(Formal), Lines), _, fail),
        catch(with_output_to(atom(Msg0),
                             print_message_lines(current_output, '', Lines)),
              _, fail)
    ->  % print_message_lines appends a newline the wire message must not carry
        atom_string(Msg0, S0),
        split_string(S0, "", " \n\t", [S1]),
        atom_string(Msg, S1)
    ;   term_to_atom(Formal, Msg)
    ).

:- multifile prolog:error_message//1.
prolog:error_message(bad_request(V)) -->
    [ 'request: expected a JSON object envelope (got ~q)'-[V] ].
prolog:error_message(bad_params(V)) -->
    [ 'params: expected a JSON object (got ~q)'-[V] ].
prolog:error_message(bad_size(V)) -->
    [ 'size must be a positive integer (got ~q)'-[V] ].
prolog:error_message(bad_mode(V)) -->
    [ 'mode must be "fixed" or "max" (got ~q)'-[V] ].
prolog:error_message(bad_seed(V)) -->
    [ 'seed must be a non-negative integer (got ~q)'-[V] ].
prolog:error_message(bad_best_effort(V)) -->
    [ 'bestEffort must be true or false (got ~q)'-[V] ].
prolog:error_message(seed_with_best_effort) -->
    [ 'seed perturbs the strict search only; not compatible with bestEffort' ].
prolog:error_message(bad_layout(V)) -->
    [ 'layout must be a JSON object (a canonical arrange result; got ~q)'-[V] ].
prolog:error_message(missing_layout) -->
    [ 'requires "layout": a canonical arrange result' ].
prolog:error_message(bad_profile(V)) -->
    [ 'unknown profile ~q; choose toc, blocked-uk, american, or barred-ximenean'-[V] ].
prolog:error_message(missing_profile) -->
    [ 'requires "profile": toc, blocked-uk, american, or barred-ximenean' ].
prolog:error_message(bad_allow_asymmetry(V)) -->
    [ 'allowAsymmetry must be true or false (got ~q)'-[V] ].
prolog:error_message(bad_to(V)) -->
    [ 'unknown format ~q; choose ipuz or exolve'-[V] ].
prolog:error_message(missing_to) -->
    [ 'requires "to": ipuz or exolve' ].
prolog:error_message(unknown_verb(V)) -->
    % The supported list is rendered from the verb/1 registry at message time
    % ({}/1), so this diagnostic cannot go stale when a verb lands.
    { findall(Verb, verb(Verb), Verbs),
      atomic_list_concat(Verbs, ', ', Supported) },
    [ 'unknown verb ~q; this engine supports: ~w'-[V, Supported] ].
prolog:error_message(unsupported_version(V)) -->
    [ 'unsupported envelope version ~q; this engine speaks v1'-[V] ].

% --- envelope builders --------------------------------------------------------
% Atom values serialize as JSON strings under json_write_dict (and the atoms
% null/true/false as the JSON literals), so `status`/`reason`/`type` atoms and
% an absent id (null) land exactly right on the wire.
success_envelope(Verb, Id, Result,
                 _{v: 1, id: Id, verb: Verb, status: success, result: Result}).
failure_envelope(Verb, Id, Detail,
                 _{v: 1, id: Id, verb: Verb, status: failure, detail: Detail}).
error_envelope(Verb, Id, Type, Msg,
               _{v: 1, id: Id, verb: Verb, status: error,
                 error: _{type: Type, message: Msg}}).
