% golden_parity.pl - browser side of the golden-parity check (native, no wasm).
%
% Reads a JSON payload file (argv[1]) exactly as the Worker would receive it and
% prints solve_browser_json/2's output to stdout. The driver golden_parity.sh
% diffs this against the bytes `crosswordsmith arrange` writes for the SAME
% request - proving the browser entry produces byte-identical output to the CLI
% (i.e. it really is the arrange verb underneath, and the spike's hand-rolled
% size/mode/drop resolution agrees with the CLI's flag resolution).
%
% Run (via the driver, not usually by hand):
%   swipl -q wasm/test/golden_parity.pl -- PAYLOAD.json

:- set_prolog_flag(verbose, silent).

:- prolog_load_context(directory, Dir),
   directory_file_path(Dir, '../client/solve_browser.pl', Rel),
   absolute_file_name(Rel, F),
   ensure_loaded(F).

:- initialization(main, main).

main :-
    current_prolog_flag(argv, [PayloadFile|_]),
    setup_call_cleanup(open(PayloadFile, read, S),
                       json_read_dict(S, In, [default_tag(json)]),
                       close(S)),
    solve_browser_json(In, Json),
    write(Json).
