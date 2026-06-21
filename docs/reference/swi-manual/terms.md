
## A.60 library(terms): Term manipulation

Compatibility  
YAP, SICStus, Quintus. Not all versions of this library define exactly the same set of predicates, but defined predicates are compatible.

Compatibility library for term manipulation predicates. Most predicates in this library are provided as SWI-Prolog built-ins.

\[det\]**term_size**(`@Term, -Size`)  
True if `Size` is the size in *cells* occupied by `Term` on the global (term) stack. A *cell* is 8 bytes. The calculation does take *sharing* into account and handles *cycles* correctly. For example:

``` code
?- A = a(1,2,3), term_size(A,S).
S = 4.
?- A = a(1,2,3), term_size(a(A,A),S).
S = 7.
?- term_size(a(a(1,2,3), a(1,2,3)), S).
S = 11.
```

Note that small objects such as atoms and small integers have a size 0. Space is allocated for floats, large integers, strings and compound terms.

\[semidet\]**variant**(`@Term1, @Term2`)  
Same as SWI-Prolog `Term1 =@= Term2`.

**subsumes_chk**(`@Generic, @Specific`)  
True if `Generic` can be made equivalent to `Specific` without changing `Specific`.

deprecated  
Replace by [subsumes_term/2](compare.html#subsumes_term/2).

**subsumes**(`+Generic, @Specific`)  
True if `Generic` is unified to `Specific` without changing `Specific`.

deprecated  
It turns out that calls to this predicate almost always should have used [subsumes_term/2](compare.html#subsumes_term/2). Also the name is misleading. In case this is really needed, one is adviced to follow [subsumes_term/2](compare.html#subsumes_term/2) with an explicit unification.

\[det\]**term_subsumer**(`+Special1, +Special2, -General`)  
`General` is the most specific term that is a generalisation of `Special1` and `Special2`. The implementation can handle cyclic terms.

author  
Inspired by LOGIC.PRO by Stephen Muggleton

Compatibility  
SICStus

**term_factorized**(`+Term, -Skeleton, -Substiution`)  
Is true when `Skeleton` is `Term` where all subterms that appear multiple times are replaced by a variable and Substitution is a list of Var=Value that provides the subterm at the location Var. I.e., After unifying all substitutions in Substiutions, `Term` `==` `Skeleton`. `Term` may be cyclic. For example:

``` code
?- X = a(X), term_factorized(b(X,X), Y, S).
Y = b(_G255, _G255),
S = [_G255=a(_G255)].
```

**mapargs**(`:Goal, ?Term1, ?Term2`)  
`Term1` and `Term2` have the same functor (name/arity) and for each matching pair of arguments `call(Goal, A1, A2)` is true.

\[det\]**mapsubterms**(`:Goal, +Term1, -Term2`)  
\[det\]**mapsubterms_var**(`:Goal, +Term1, -Term2`)  
Recursively map sub terms of `Term1` into subterms of `Term2` for every pair for which `call(Goal, ST1, ST2)` succeeds. Procedurably, the mapping for each (sub) term pair `T1/T2` is defined as:

- If `T1` is a variable
  - [mapsubterms/3](terms.html#mapsubterms/3) unifies `T2` with `T1`.
  - [mapsubterms_var/3](terms.html#mapsubterms_var/3) treats variables as other terms.
- If `call(Goal, T1, T2)` succeeds we are done. Note that the mapping does not continue in `T2`. If this is desired, `Goal` must call [mapsubterms/3](terms.html#mapsubterms/3) explicitly as part of its conversion.
- If `T1` is a dict, map all values, i.e., the *tag* and *keys* are left untouched.
- If `T1` is a list, map all elements, i.e., the list structure is left untouched.
- If `T1` is a compound, use [same_functor/3](terms.html#same_functor/3) to instantiate `T2` and recurse over the term arguments left to right.
- Otherwise `T2` is unified with `T1`.

Both predicates are implemented using [foldsubterms/5](terms.html#foldsubterms/5).

\[semidet\]**foldsubterms**(`:Goal3, +Term1, +State0, -State`)  
\[semidet\]**foldsubterms**(`:Goal4, +Term1, ?Term2, +State0, -State`)  
The predicate [foldsubterms/5](terms.html#foldsubterms/5) calls `call(Goal4, SubTerm1, SubTerm2, StateIn, StateOut)` for each subterm, including variables, in `Term1`. If this call fails, `StateIn` and `StateOut` are the same. This predicate may be used to map subterms in a term while collecting state about the mapped subterms. The [foldsubterms/4](terms.html#foldsubterms/4) variant does not map the term.

\[semidet\]**same_functor**(`?Term1, ?Term2`)  
\[semidet\]**same_functor**(`?Term1, ?Term2, -Arity`)  
\[semidet\]**same_functor**(`?Term1, ?Term2, ?Name, ?Arity`)  
True when `Term1` and `Term2` are terms that have the same functor (`Name`/`Arity`). The arguments must be sufficiently instantiated, which means either `Term1` or `Term2` must be bound or both `Name` and `Arity` must be bound.

If `Arity` is 0, `Term1` and `Term2` are unified with `Name` for compatibility.

Compatibility  
SICStus
