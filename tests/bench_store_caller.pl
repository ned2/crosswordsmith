:- module(bench_store_test_caller, []).

:- use_module('../benchmarks/bench_store.pl').

build_recorded(Baseline, Run, Recorded) :-
    build_recorded_baseline(Baseline, Run, row_key, row_spec, Recorded).

row_key(Row, Row.id).

row_spec(Row, existing(Old), Spec) :-
    Spec = Old.put(value, Row.value).
row_spec(Row, new, _{value:Row.value, tier:Row.tier}).

replace_rejected(Path, Dict) :-
    replace_json_dict(Path, Dict, 80, reject_readback).

reject_readback(_) :-
    fail.
