% solve_browser.pl - the WASM client's Prolog load root (a thin adapter).
%
% This file is what wasm/build/build-wasm.sh step 4 qcompiles into
% crosswordsmith.qlf (include(user) folds the whole dependency closure in), so
% its NAME is load-bearing for the build; keep it. The actual browser API is
% crosswordsmith_browser:browser_dispatch/3 (prolog/crosswordsmith/browser.pl):
% JSON-string request in, discriminated JSON envelope out, per
% docs/plans/wasm-sdk-strategy.md §3/§4. load.pl imports that module's export
% into `user`, so the Worker's forEach goal
%
%     browser_dispatch(Verb, Payload, Json)
%
% (Verb + Payload ride as BOUND variables, never string-concatenated) resolves
% here without any Module:Pred reach — the spike's internals leak is gone.
%
% JSON IN, JSON OUT is mandatory for THIS app - learned the hard way in the
% browser. The tempting "query binding" round-trip (pass a JS object straight
% in, get a JS object back via SWI's WASM data translation, §13.2.3) is LOSSY
% in BOTH directions here: JS strings arrive as Prolog ATOMS (breaking core's
% string("answer") check) and Prolog null/true/false come back as the JS
% strings "null"/"true"/"false" (corrupting the layout schema's empty cells).
% (atom_json_dict lives in core library(json) — imported explicitly below. The
% legacy library(http/json) compat shim does not resolve under wasm at all -
% see wasm/README.md "Browser gotchas"; C6 swapped every import to the core
% library.)

% Load the whole implementation via load.pl, resolved relative to this file so
% it works from any cwd and under `node src/swipl.js` in the WASM build tree.
% directory_file_path/3 is autoload-only (library(filesex)); explicit so this
% root also loads under autoload(false) (P11/C5, matching load.pl).
:- use_module(library(filesex), [directory_file_path/3]).
% library(json), NOT the legacy library(http/json) alias (C6): the selftest
% below round-trips the dispatch envelope through atom_json_dict/3.
:- use_module(library(json), [atom_json_dict/3]).
:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, '../../load.pl', LoadRel),
   absolute_file_name(LoadRel, LoadFile),
   ensure_loaded(LoadFile).

% --- self-test (native swipl): proves the file-free dispatch round-trip ------
% Run:  swipl -g browser_selftest -t halt wasm/client/solve_browser.pl

browser_selftest :-
    atom_json_dict(Payload,
                   _{ v: 1, id: "selftest", verb: "arrange",
                      params: _{ clues: [ _{answer: "CAT"}, _{answer: "CAR"},
                                          _{answer: "ARC"}, _{answer: "RAT"},
                                          _{answer: "TAR"} ],
                                 size: 5 } },
                   []),
    browser_dispatch(arrange, Payload, Json),
    atom_json_dict(Json, Env, [default_tag(json)]),
    (   get_dict(status, Env, "success"),
        get_dict(result, Env, R),
        get_dict(gridLength, R, GL),
        get_dict(words, R, Ws), length(Ws, NW)
    ->  format("dispatch OK: status success, gridLength=~w, ~w placed words~n",
               [GL, NW])
    ;   format("dispatch FAILED: ~w~n", [Json]), fail
    ).
