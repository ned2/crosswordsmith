% solve_browser.pl - SPIKE: a file-free browser entry point for `arrange`.
%
% Throwaway exploration for the WASM/browser milestone (see
% docs/plans/wasm-browser-deployment.md). It demonstrates the two documented
% JS<->Prolog I/O patterns without touching the filesystem:
%
%   solve_browser/2       InputDict  -> ResultDict   (query-binding round-trip;
%                                                      the RECOMMENDED pattern)
%   solve_browser_json/2  InputDict  -> JSON string  (byte-identical to the CLI
%                                                      output; the safe fallback)
%
% In the browser these are called as, e.g.
%     Prolog.query("solve_browser(In, Out)",
%                  {In: {clues:[{answer:"CAT"}, ...], size:5}}).once().Out
% where the JS object becomes a Prolog dict on the way in and the Prolog dict
% comes back as a JS object (SWI's WASM data translation, §13.2.3). No `open/3`,
% no MEMFS, no json_read_dict-from-file.
%
% NB (spike scope): this reaches into module internals as `Module:Pred(...)`
% (doc_to_words/2, arrange_best_layout/5, arrange_diag_layout_dict/5 are not
% exported). The PRODUCTION version should add a proper exported
% `solve_browser/2` to a small `browser.pl` module (or core) so the public
% API is explicit - see the plan doc. Grid size and mode are taken straight
% from the input dict here; the CLI's --max-size cropping / fragment / seed
% machinery is intentionally out of scope for the spike.

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
solve_browser_json(InputDict, Json) :-
    crosswordsmith_core:doc_to_words(InputDict, Words),
    grid_size(InputDict, GridLen),
    drop_contract_of(InputDict, Drop),
    size_mode(InputDict, SizeMode),
    with_output_to(atom(Json),
                   crosswordsmith_arrange:arrange_solve(Words, GridLen, Drop, SizeMode)).

% --- input helpers (spike-simple; production reuses the CLI's resolvers) -----

grid_size(InputDict, GridLen) :-
    ( get_dict(size, InputDict, N), integer(N), N > 0
    -> GridLen = N
    ;  throw(error(browser_missing_size, _))
    ).

size_mode(InputDict, Mode) :-
    ( get_dict(mode, InputDict, M), memberchk(M, [fixed, max]) -> Mode = M ; Mode = fixed ).

drop_contract_of(InputDict, Drop) :-
    ( get_dict(bestEffort, InputDict, true) -> Drop = best_effort ; Drop = strict ).

% --- self-test (native swipl): proves the file-free round-trip ---------------
% Run:  swipl -g browser_selftest -t halt spikes/wasm-browser/solve_browser.pl

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
        )
    ;   format("dict path FAILED~n", []), fail
    ).
