
## 6.14 Transparent predicates: definition and context module

*The‘module-transparent’mechanism is still underlying the actual implementation. Direct usage by programmers is deprecated. Please use [meta_predicate/1](metapred.html#meta_predicate/1) to deal with meta-predicates.*

The qualification of module-sensitive arguments described in [section 6.5](metapred.html#sec:6.5) is realised using *transparent* predicates. It is now deprecated to use this mechanism directly. However, studying the underlying mechanism helps to understand SWI-Prolog's modules. In some respect, the transparent mechanism is more powerful than meta-predicate declarations.

Each predicate of the program is assigned a module, called its *definition module*. The definition module of a predicate is always the module in which the predicate was originally defined. Each active goal in the Prolog system has a *context module* assigned to it.

The context module is used to find predicates for a Prolog term. By default, the context module is the definition module of the predicate running the goal. For transparent predicates, however, this is the context module of the goal inherited from the parent goal. Below, we implement [maplist/3](apply.html#maplist/3) using the transparent mechanism. The code of [maplist/3](apply.html#maplist/3) and maplist\_/3 is the same as in [section 6.5](metapred.html#sec:6.5), but now we must declare both the main predicate and the helper as transparent to avoid changing the context module when calling the helper.

``` code
:- module(maplist, maplist/3).

:- module_transparent
        maplist/3,
        maplist_/3.

maplist(Goal, L1, L2) :-
        maplist_(L1, L2, G).

maplist_([], [], _).
maplist_([H0|T0], [H|T], Goal) :-
        call(Goal, H0, H),
        maplist_(T0, T, Goal).
```

Note that *any* call that translates terms into predicates is subject to the transparent mechanism, not just the terms passed to module-sensitive arguments. For example, the module below counts the number of unique atoms returned as bindings for a variable. It works as expected. If we use the directive `:- module_transparent ``count_atom_results/3.` instead, atom_result/2 is called wrongly in the module *calling* count_atom_results/3 . This can be solved using [strip_module/3](ctxmodule.html#strip_module/3) to create a qualified goal and a non-transparent helper predicate that is defined in the same module.

``` code
:- module(count_atom_results,
          [ count_atom_results/3
          ]).
:- meta_predicate count_atom_results(-,0,-).

count_atom_results(A, Goal, Count) :-
        setof(A, atom_result(A, Goal), As), !,
        length(As, Count).
count_atom_results(_, _, 0).

atom_result(Var, Goal) :-
        call(Goal),
        atom(Var).
```

The following predicates support the module-transparent interface:

:- **module_transparent**(`+Preds`)  
`Preds` is a comma-separated list of name/arity pairs (like [dynamic/1](dynamic.html#dynamic/1)). Each goal associated with a transparent-declared predicate will inherit the *context module* from its parent goal.

**context_module**(`-Module`)  
Unify `Module` with the context module of the current goal. [context_module/1](ctxmodule.html#context_module/1) itself is, of course, transparent.

**strip_module**(`+Term, -Module, -Plain`)  
Used in module-transparent predicates or meta-predicates to extract the referenced module and plain term. If `Term` is a module-qualified term, i.e. of the format `Module`:`Plain`, `Module` and `Plain` are unified to these values. Otherwise, `Plain` is unified to `Term` and `Module` to the context module.
