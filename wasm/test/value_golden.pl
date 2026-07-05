% value_golden.pl - dual-mode driver for the wasm value-correctness golden
% (wasm-sdk-strategy §7; supersedes the spike's golden_parity.pl).
%
% DEC-8: the envelope nests `result` as a live dict, so correctness against
% the CLI reference is VALUE-equality of the parsed structures, never a byte
% diff (nesting re-indents, flips spaces->tabs, drops the trailing newline).
% The orchestrator (value_golden.sh) runs this file in two modes:
%
%   emit-request OUT CLUES EXTRA   (native swipl)
%       Build the v1 request envelope for the clue file CLUES (.json or .pl —
%       whatever load_clues/2 reads, so the SAME fixture feeds both the CLI
%       reference and the browser path) with the EXTRA params JSON (e.g.
%       '{"size":17,"seed":42}') merged in, and write it to OUT.
%
%   dispatch REQFILE...            (node build.wasm/src/swipl.js — the wasm VM)
%       Run each request file through browser_dispatch/3 on THIS ONE Prolog
%       instance, printing one envelope JSON per line. Running several
%       requests through one instance is the point: repeated/interleaved
%       seeded+seedless requests prove the per-request state reset under the
%       REAL wasm VM (same-instance determinism), not just under native.

:- set_prolog_flag(verbose, silent).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root, [relative_to(Dir), file_type(directory)]),
   directory_file_path(Root, 'load.pl', Load),
   ensure_loaded(Load).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, Argv),
    (   run(Argv)
    ->  halt(0)
    ;   format(user_error, "value_golden.pl: bad usage or failed: ~q~n", [Argv]),
        halt(2)
    ).

run(['emit-request', Out, CluesFile, ExtraJson]) :-
    load_clues(CluesFile, Words),
    findall(C, ( member(W, Words), word_clue(W, C) ), Clues),
    atom_json_dict(ExtraJson, Extra, [default_tag(json)]),
    put_dict(Extra, _{clues: Clues}, Params),
    Req = _{v: 1, id: "vg", verb: "arrange", params: Params},
    setup_call_cleanup(open(Out, write, S),
                       json_write_dict(S, Req),
                       close(S)).
run(['dispatch'|Files]) :-
    Files = [_|_],
    forall(member(F, Files), dispatch_file(F)).

% Fixture entries are [Answer, Meta] or bare [Answer]; answers go over the
% wire as JSON strings exactly as a JS caller would send them.
word_clue([A, M], _{answer: S, meta: M}) :- !, atom_string(A, S).
word_clue([A], _{answer: S}) :- atom_string(A, S).

% Dispatch on the request's own verb — exactly what the Worker does (it binds
% request.verb as the forEach Verb variable).
dispatch_file(F) :-
    read_file_to_string(F, Payload, []),
    atom_json_dict(Payload, Req, [default_tag(json)]),
    get_dict(verb, Req, VerbS),
    atom_string(Verb, VerbS),
    browser_dispatch(Verb, Payload, Json),
    format("~w~n", [Json]).
