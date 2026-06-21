
## A.31 library(occurs): Finding and counting sub-terms

See also  
`library(terms)` provides similar predicates and is probably more wide-spread than this library.

This is a SWI-Prolog implementation of the corresponding Quintus library, based on the generalised [arg/3](manipterm.html#arg/3) predicate of SWI-Prolog.

\[semidet\]**contains_term**(`+Sub, +Term`)  
Succeeds if `Sub` is contained in `Term` (=, deterministically)

\[semidet\]**contains_var**(`+Sub, +Term`)  
Succeeds if `Sub` is contained in `Term` (==, deterministically)

\[semidet\]**free_of_term**(`+Sub, +Term`)  
Succeeds of `Sub` does not unify to any subterm of `Term`

\[semidet\]**free_of_var**(`+Sub, +Term`)  
Succeeds of `Sub` is not equal (`==`) to any subterm of `Term`

\[det\]**occurrences_of_term**(`@SubTerm, @Term, ?Count`)  
`Count` the number of SubTerms in `Term` that *unify* with `SubTerm`. As this predicate is implemented using backtracking, `SubTerm` and `Term` are not further instantiated. Possible constraints are enforced. For example, we can count the integers in `Term` using

``` code
?- freeze(S, integer(S)), occurrences_of_term(S, f(1,2,a), C).
C = 2,
freeze(S, integer(S)).
```

See also  
[occurrences_of_var/3](occurs.html#occurrences_of_var/3) for an equality ([==/2](compare.html#==/2)) based variant.

\[det\]**occurrences_of_var**(`@SubTerm, @Term, ?Count`)  
`Count` the number of SubTerms in `Term` that are *equal* to `SubTerm`. Equality is tested using [==/2](compare.html#==/2). Can be used to count the occurrences of a particular variable in `Term`.

See also  
[occurrences_of_term/3](occurs.html#occurrences_of_term/3) for a unification ([=/2](compare.html#=/2)) based variant.

**sub_term**(`-Sub, +Term`)  
Generates (on backtracking) all subterms of `Term`.

**sub_var**(`-Sub, +Term`)  
Generates (on backtracking) all subterms (`==`) of `Term`.

\[det\]**sub_term_shared_variables**(`+Sub, +Term, -Vars`)  
If `Sub` is a sub term of `Term`, `Vars` is bound to the list of variables in `Sub` that also appear outside `Sub` in `Term`. Note that if `Sub` appears twice in `Term`, its variables are all considered shared.

An example use-case is refactoring a large clause body by introducing intermediate predicates. This predicate can be used to find the arguments that must be passed to the new predicate.
