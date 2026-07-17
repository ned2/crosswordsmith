#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   directory_file_path(Dir, 'ad2_delta.pl', AD2), use_module(AD2),
   directory_file_path(Dir, 'closeout_direct.pl', Direct), use_module(Direct).

:- initialization(main, main).

control(easy, 'fixtures/bundled_17_clues.pl', 17, 6, topleft_across, none).
control(seeded, 'fixtures/ladder_15x15_12w.pl', 15, 12, topright, 42).
control(dense, 'fixtures/ladder_21x21_80w.pl', 21, 80, topleft_across, none).

main :- catch(run, Error, (print_message(error, Error), halt(2))).

run :-
    forall(control(Name, Fixture, Grid, Count, Corner, Seed),
           verify(Name, Fixture, Grid, Count, Corner, Seed)),
    format('closeout-direct: 3 replay controls passed~n').

verify(Name, Fixture, Grid, Count, Corner, Seed) :-
    format(user_error, 'closeout-direct test=~w phase=start~n', [Name]),
    probe_arrange:load_fixture(Fixture, Count, Words),
    product_result(Seed, Words, Grid, Corner, Product),
    probe_arrange_closeout_direct:closeout_corner(
        Words, Grid, Corner, Seed, First, Stats, Trace),
    probe_arrange_closeout_direct:closeout_corner(
        Words, Grid, Corner, Seed, Second, Stats2, Trace2),
    probe_arrange_ad2:ad2_corner(
        Words, Grid, Corner, Seed, Differential, _Summary, _Counts, AD2Trace),
    require_identical(product_result, Product, First),
    require_identical(duplicate_result, First, Second),
    require_identical(differential_result, First, Differential),
    require_identical(duplicate_stats, Stats, Stats2),
    require_identical(duplicate_trace, Trace, Trace2),
    require_identical(differential_trace, Trace, AD2Trace),
    Stats.legal_decisions =:= Stats.places,
    format(user_error, 'closeout-direct test=~w phase=done~n', [Name]).

product_result(Seed, Words, Grid, Corner, Result) :-
    ( Seed == none
    -> probe_arrange_closeout_direct:product_corner(Words, Grid, Corner, Result)
    ; setup_call_cleanup(
          crosswordsmith_core:set_search_seed(Seed),
          ( crosswordsmith_core:reset_search_memos,
            crosswordsmith_core:init_gs(Grid, GS),
            crosswordsmith_core:start_loc(Corner, Grid, Start, Dir),
            once(crosswordsmith_core:assign_words_direct(
                     Words, [], Grid, Start, Dir, GS, _, Placed)) ),
          crosswordsmith_core:set_search_seed(-1)),
      crosswordsmith_core:assign_clue_numbers(Placed, Numbered),
      crosswordsmith_arrange:layout_reward(5, 1, Numbered, Reward),
      Result = result(Numbered, Reward)
    ).

require_identical(_, A, B) :- A == B, !.
require_identical(Label, A, B) :-
    variant_sha1(A, HashA),
    variant_sha1(B, HashB),
    throw(error(closeout_direct_mismatch(Label,HashA,HashB), _)).
