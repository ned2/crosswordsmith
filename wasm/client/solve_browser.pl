% solve_browser.pl - the WASM client's Prolog load root (browser-specific).
%
% This file is what wasm/build/build-wasm.sh step 4 qcompiles into
% crosswordsmith.qlf (include(user) folds the whole dependency closure in), so
% its NAME is load-bearing for the build; keep it. The actual browser API is
% crosswordsmith_browser:browser_dispatch/3 (prolog/crosswordsmith/browser.pl):
% JSON-string request in, discriminated JSON envelope out, per
% docs/plans/wasm-sdk-strategy.md §3/§4. The imports below land that module's
% export in `user`, so the Worker's forEach goal
%
%     browser_dispatch(Verb, Payload, Json)
%
% (Verb + Payload ride as BOUND variables, never string-concatenated) resolves
% here without any Module:Pred reach — the spike's internals leak is gone.
%
% This root deliberately does NOT go through the general load.pl: that file
% owns the FULL product load order (CLI, tests, benchmarks) including
% `stockgrid` and `fill`, but the browser contract exposes only arrange, lint,
% export, and capabilities (fill is unbrowserified — strategy §6 phase 4). A
% browser-specific root keeps fill-only baggage (library(sha), library(fastrw),
% dictionary/stock-grid handling) OUT of the qlf, and — more importantly — it
% DEFINES the dependency closure the reduced preload/package profiles are
% derived from (docs/plans/wasm-payload-performance.md Phase 1: this list is
% the input to Phases 2-3). Widening the browser surface = adding the module
% here + its dispatch clause in browser.pl; load.pl stays the record for
% everything else.
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

% library(json), NOT the legacy library(http/json) alias (C6): the selftest
% below round-trips the dispatch envelope through atom_json_dict/3.
:- use_module(library(json), [atom_json_dict/3]).

% The `crosswordsmith` file-search alias, resolved relative to THIS file
% (cwd-independent, works under `node src/swipl.js` in the WASM build tree and
% under qcompile). The modules resolve each other through the alias, so it
% must exist before the first use_module below. atom_concat/3 on purpose,
% TWICE over: (a) unlike load.pl's directory_file_path/3 it needs no
% library(filesex), keeping file-utility baggage out of the qlf closure
% (payload plan §5); (b) unlike absolute_file_name/3 with file_type(directory)
% it never PROBES the path — qlf-embedded directives re-execute at load time,
% and in the browser that probe becomes a 404'd URL fetch + an unhandled-
% exception error line per worker boot (the alias is dead there anyway: every
% module is already compiled into the qlf, so it must merely assert quietly,
% exactly as load.pl's concat-only helper does).
% Guarded like load.pl so loading both roots in one process stays harmless.
:- prolog_load_context(directory, Dir),
   atom_concat(Dir, '/../../prolog/crosswordsmith', LibDir),
   (   user:file_search_path(crosswordsmith, LibDir)
   ->  true
   ;   assertz(user:file_search_path(crosswordsmith, LibDir))
   ).

% Exactly the browser surface (payload plan §5): the shared substrate, the
% metric predicates, the three browserified verbs' engines, and the dispatch
% spine. NO stockgrid, NO fill.
:- use_module(crosswordsmith(core)).
:- use_module(crosswordsmith(metrics)).
:- use_module(crosswordsmith(arrange)).
:- use_module(crosswordsmith(lint)).
:- use_module(crosswordsmith(export)).
:- use_module(crosswordsmith(browser)).

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
