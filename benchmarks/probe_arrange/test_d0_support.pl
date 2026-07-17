#!/usr/bin/env swipl
:- set_prolog_flag(verbose, silent).
:- use_module(library(plunit)).

:- prolog_load_context(directory, Dir),
   absolute_file_name('../..', Root,
                      [relative_to(Dir), file_type(directory), access(read)]),
   directory_file_path(Root, 'load.pl', Load), consult(Load),
   directory_file_path(Dir, 'probe_arrange.pl', Probe), use_module(Probe),
   directory_file_path(Dir, 'd0_support.pl', D0), use_module(D0).

:- initialization(main, main).

main :-
    ( run_tests -> halt(0) ; halt(1) ).

:- begin_tests(probe_arrange_d0).

test(shadow_replays_product_and_decision_order) :-
    fixture(Words),
    probe_arrange:authority_corner(
        Words, 9, topleft_across, 500000000, null, Authority, _),
    probe_arrange_d0:d0_corner(
        Words, 9, topleft_across, null, false,
        Base, _BaseSummary, BaseDecisions, _),
    probe_arrange_d0:d0_corner(
        Words, 9, topleft_across, null, true,
        Observed, Summary, ObservedDecisions, _),
    equivalent(Authority, Base),
    equivalent(Authority, Observed),
    BaseDecisions == ObservedDecisions,
    Summary.classified_refreshes =:= Summary.verified_classifications,
    Summary.proposed_candidate_checks =< Summary.exact_candidate_checks.

test(duplicate_counter_passes_are_exact) :-
    fixture(Words),
    probe_arrange_d0:d0_corner(
        Words, 9, topright, null, true, _Warm, _WarmSummary, _WarmDecisions, _),
    probe_arrange_d0:d0_corner(
        Words, 9, topright, null, true, R1, S1, D1, T1),
    probe_arrange_d0:d0_corner(
        Words, 9, topright, null, true, R2, S2, D2, T2),
    R1 == R2,
    S1 =@= S2,
    D1 == D2,
    get_dict(inferences, T1, Inferences1),
    get_dict(inferences, T2, Inferences2),
    Inferences1 =:= Inferences2.

test(proof_multiplicity_is_observed) :-
    fixture(Words),
    probe_arrange_d0:d0_corner(
        Words, 9, topleft_across, null, true, _Result, Summary, _Decisions, _),
    Summary.proof_geometry_divergences > 0,
    Summary.exact_proofs > Summary.exact_unique_geometries.

test(two_corner_operation_replays_product) :-
    fixture(Words),
    probe_arrange:authority_operation(
        Words, 9, 500000000, null, Authority, _),
    probe_arrange_d0:d0_operation(
        Words, 9, null, false, Base, _BaseSummary, BaseDecisions, _),
    probe_arrange_d0:d0_operation(
        Words, 9, null, true, Observed, _Summary, ObservedDecisions, _),
    equivalent(Authority, Base),
    equivalent(Authority, Observed),
    BaseDecisions == ObservedDecisions.

fixture(Words) :-
    probe_arrange:load_fixture('fixtures/ladder_09x09_08w.pl', 8, Words).

equivalent(result(placed,ok,false,_,A,Reward),
           result(placed,ok,false,_,B,Reward)) :-
    probe_arrange:layout_signature(A, Signature),
    probe_arrange:layout_signature(B, Signature).

:- end_tests(probe_arrange_d0).
