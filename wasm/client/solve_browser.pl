% solve_browser.pl - file-free browser entry points for `arrange`.
%
% The Prolog side of the WASM browser client (see wasm/README.md and
% docs/plans/wasm-browser-deployment.md). Runs the solver with NO filesystem I/O.
% Still spike-grade: production promotes this to an exported `browser.pl` module
% (see the "NB" below and plan §9).
%
% Three entry points; the Worker calls solve_browser_str/2:
%
%   solve_browser_str/2   JSON string -> JSON string   (what the Worker calls)
%   solve_browser_json/2  InputDict   -> JSON string    (dict in, JSON out)
%   solve_browser/2       InputDict   -> ResultDict      (dict in, dict out; illustrative)
%
% JSON IN, JSON OUT is the right pattern for THIS app - learned the hard way in
% the browser. The tempting "query binding" round-trip (pass a JS object straight
% in, get a JS object back via SWI's WASM data translation, §13.2.3) is LOSSY in
% BOTH directions here: JS strings arrive as Prolog ATOMS (breaking core's
% string("answer") check) and Prolog null/true/false come back as the JS strings
% "null"/"true"/"false" (corrupting the layout schema's empty cells). Going
% through JSON both ways lands every value with the right type and stays
% byte-identical to the golden CLI output. No open/3, no MEMFS, no
% json_read_dict-from-file either way.
%
% NB (spike scope): this reaches into module internals as `Module:Pred(...)`
% (doc_to_words/2, arrange_best_layout/5, arrange_diag_layout_dict/5 are not
% exported). The PRODUCTION version should add a proper exported entry point to a
% small `browser.pl` module (or core) so the public API is explicit - see plan
% §9. Grid size is taken straight from the input dict; the CLI's --max-size
% cropping, fragment and seed machinery are not wired in - and `mode:"max"` is
% rejected outright by size_mode/2 (below) so an untested-but-reachable crop path
% can't ship. Only mode:"fixed" (the default) is supported.

% Load the whole implementation via load.pl, resolved relative to this file so
% it works from any cwd and under `node src/swipl.js` in the WASM build tree.
:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, '../../load.pl', LoadRel),
   absolute_file_name(LoadRel, LoadFile),
   ensure_loaded(LoadFile).

:- use_module(library(http/json)).

% --- the browser API --------------------------------------------------------

% solve_browser(+InputDict, -ResultDict)
% InputDict is the JS payload as a Prolog dict, e.g.
%   _{clues: [_{answer:"CAT"}, _{answer:"CAR"}, ...], size: 5, mode: fixed}
% ResultDict is the canonical arrange layout dict (identical shape to the JSON
% the CLI writes) - returned straight back to JS as an object.
solve_browser(InputDict, ResultDict) :-
    crosswordsmith_core:doc_to_words(InputDict, Words),
    grid_size(InputDict, GridLen),
    size_mode(InputDict, SizeMode),
    crosswordsmith_arrange:arrange_best_layout(Words, GridLen, Numbered, _Reward, _Outcome),
    crosswordsmith_arrange:arrange_diag_layout_dict(Numbered, Words, GridLen, SizeMode, ResultDict).

% solve_browser_json(+InputDict, -Json)
% Same input, but returns the exact JSON the CLI would write to --out, as a
% Prolog ATOM. JS does JSON.parse(Json).
%
% For THIS app this is the recommended OUTPUT path, not the dict above: the
% layout schema is full of JSON `null` (every empty cell) and `across`/`down`
% values, and per the WASM reverse-translation table Prolog's `null`/`true`/
% `false` atoms come back as the JS *strings* "null"/"true"/"false", not JS
% null/booleans. Serialising with json_write_dict and JSON.parse-ing on the JS
% side preserves null/bool fidelity and stays byte-identical to the golden CLI
% output. (Return an ATOM, not a string: an atom -> a clean JS String, whereas
% a Prolog string comes back wrapped as a Prolog.String instance.)
%
% ‼️ STATE-RESET INVARIANT (owed before any state-bearing key is added): this does
% NO per-request reset. It is safe ONLY because no `seed`/`checkTarget`/`verbose`
% key is read here, so the reused instance's search_seed/1 & check_target_override/1
% stay unset across solves ({engine:true} does NOT clear those module globals). The
% instant such a key is wired in, add an unconditional reset at the TOP of this
% predicate (set_search_seed(-1), set_check_target(-1), set_verbose(false)) and a
% same-worker determinism test - see plan §9.1 / §10.2.
solve_browser_json(InputDict, Json) :-
    crosswordsmith_core:doc_to_words(InputDict, Words),
    grid_size(InputDict, GridLen),
    drop_contract_of(InputDict, Drop),
    size_mode(InputDict, SizeMode),
    with_output_to(atom(Json),
                   crosswordsmith_arrange:arrange_solve(Words, GridLen, Drop, SizeMode)).

% solve_browser_str(+PayloadAtom, -Json)
% The Worker hands us the input as a single JSON string (JSON.stringify), NOT a
% JS object. We parse it here so every value lands with the RIGHT Prolog type.
% Rationale (found the hard way in the browser): the JS->Prolog query-binding
% path delivers JS strings as ATOMS, so a clue's "answer" fails core's
% string(RawAnswer) check ("every entry needs a string answer"). That is the
% exact input-side twin of the output-side quirk where Prolog null/true/false
% come back as JS strings. Going through JSON both ways is symmetric and
% translation-proof: json_read_dict yields strings, numbers, null and bools with
% correct types. (json_read_dict/atom_json_dict autoload from library(json) under
% wasm even though `use_module(library(http/json))` cannot resolve the compat
% shim there - see wasm/README.md "Browser gotchas".)
solve_browser_str(Payload, Json) :-
    setup_call_cleanup(
        open_string(Payload, In),
        json_read_dict(In, InputDict, [default_tag(json)]),
        close(In)),
    solve_browser_json(InputDict, Json).

% --- input helpers (spike-simple; production reuses the CLI's resolvers) -----

grid_size(InputDict, GridLen) :-
    ( get_dict(size, InputDict, N), integer(N), N > 0
    -> GridLen = N
    ;  throw(error(browser_missing_size, _))
    ).

size_mode(InputDict, Mode) :-
    ( get_dict(mode, InputDict, M)
    ->  ( M == max
        ->  % --max-size cropping IS reachable through arrange_solve/4, but it is
            % not golden-tested against the CLI yet (plan §9.1 Part B). Reject it
            % explicitly rather than silently shipping an unverified path.
            throw(error(browser_max_mode_unsupported, _))
        ;   Mode = fixed        % 'fixed' (or any unknown mode) -> fixed
        )
    ;   Mode = fixed
    ).

drop_contract_of(InputDict, Drop) :-
    ( get_dict(bestEffort, InputDict, true) -> Drop = best_effort ; Drop = strict ).

% --- self-test (native swipl): proves the file-free round-trip ---------------
% Run:  swipl -g browser_selftest -t halt wasm/client/solve_browser.pl

browser_selftest :-
    In = _{ clues: [ _{answer:"CAT"}, _{answer:"CAR"}, _{answer:"ARC"},
                     _{answer:"RAT"}, _{answer:"TAR"} ],
            size: 5, mode: fixed, bestEffort: true },
    ( solve_browser(In, Out)
    ->  get_dict(gridLength, Out, GL),
        get_dict(words, Out, WordObjs),
        length(WordObjs, NW),
        format("dict path OK: gridLength=~w, ~w placed words~n", [GL, NW]),
        ( solve_browser_json(In, Json)
        ->  string_length(Json, L),
            format("json path OK: ~w chars~n", [L])
        ;   format("json path FAILED~n", []), fail
        ),
        % str path: exactly what the Worker sends (a JSON string, answers as
        % strings) - proves the JSON-in adapter + json autoload.
        atom_json_dict(Payload, In, []),
        ( solve_browser_str(Payload, Json2)
        ->  string_length(Json2, L2),
            format("str path OK: ~w chars~n", [L2])
        ;   format("str path FAILED~n", []), fail
        )
    ;   format("dict path FAILED~n", []), fail
    ).
