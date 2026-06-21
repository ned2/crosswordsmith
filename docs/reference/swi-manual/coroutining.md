
## 8.2 Coroutining

Coroutining allows us to delay the execution of Prolog goals until their truth can be safely decided.

Among the most important coroutining predicates is [dif/2](coroutining.html#dif/2), which expresses *disequality* of terms in a sound way. The actual test is delayed until the terms are sufficiently different, or have become identical. For example:

``` code
?- dif(X, Y), X = a, Y = b.
X = a,
Y = b.

?- dif(X, Y), X = a, Y = a.
false.
```

There are also lower-level coroutining predicates that are intended as building blocks for higher-level constraints. For example, we can use [freeze/2](coroutining.html#freeze/2) to define a variable that can only be assigned an atom:

``` code
?- freeze(X, atom(X)), X = a.
X = a.
```

In this case, calling [atom/1](typetest.html#atom/1) earlier causes the whole query to fail:

``` code
?- atom(X), X = a.
false.
```

If available, domain-specific constraints should be used in such cases. For example, to state that a variable can only assume even integers, use the CLP(FD) constraint [\#=/2](clpfd.html##=/2):

``` code
?- X mod 2 #= 0.
X mod 2#=0.
```

Importantly, domain-specific constraints can apply stronger propagation by exploiting logical properties of their respective domains. For example:

``` code
?- X mod 2 #= 0, X in 1..3.
X = 2.
```

Remaining constraints, such as `X mod 2#=0` in the example above, are called *residual* goals. They are said to *flounder*, because their truth is not yet decided. Declaratively, the query is only true if all residual goals are satisfiable. Use [call_residue_vars/2](coroutining.html#call_residue_vars/2) to collect all variables that are involved in constraints.

**dif**(`@A, @B`)  
The [dif/2](coroutining.html#dif/2) predicate is a *constraint* that is true if and only if `A` and `B` are different terms. If `A` and `B` can never unify, [dif/2](coroutining.html#dif/2) succeeds deterministically. If `A` and `B` are identical, it fails immediately. Finally, if `A` and `B` can unify, goals are delayed that prevent `A` and `B` to become equal. It is this last property that makes [dif/2](coroutining.html#dif/2) a more general and more declarative alternative for [\\/2](compare.html#\=/2) and related predicates.

This predicate behaves as if defined by `dif(X, Y) :- when(?=(X,Y), X \== Y)`. See also [?=/2](compare.html#?=/2). The implementation can deal with cyclic terms.

The [dif/2](coroutining.html#dif/2) predicate is realised using attributed variables associated with the module `dif`. It is an autoloaded predicate that is defined in the library `library(dif)`.

**freeze**(`+Var, :Goal`)  
Delay the execution of `Goal` until `Var` is bound (i.e., is not a variable or attributed variable). If `Var` is bound on entry [freeze/2](coroutining.html#freeze/2) is equivalent to [call/1](metacall.html#call/1). The [freeze/2](coroutining.html#freeze/2) predicate is realised using an attributed variable associated with the module `freeze`. See also [frozen/2](coroutining.html#frozen/2).

\[det\]**frozen**(`@Term, -Goal`)  
Unify `Goal` with the goal or conjunction of goals delayed on some attributed variable in `Term`. If `Term` is free of attributed variables, `Goal` is unified to `true`. Note that [frozen/2](coroutining.html#frozen/2) reports all delayed goals, not only those delayed due to [freeze/2](coroutining.html#freeze/2). The goals are extracted using [copy_term/3](attvar.html#copy_term/3).^(194Versions prior to 8.3.7 only report goals delayed using [freeze/2](coroutining.html#freeze/2) on a plain variable. The new behaviour is compatible with SICStus.) See also [term_attvars/2](attvar.html#term_attvars/2) and [call_residue_vars/2](coroutining.html#call_residue_vars/2).

**when**(`@Condition, :Goal`)  
Execute `Goal` when `Condition` becomes true. `Condition` is one of `?=(X, Y)`, `nonvar(X)`, `ground(X)`, `,``(Cond1, Cond2)` or `;``(Cond1, Cond2)`. See also [freeze/2](coroutining.html#freeze/2) and [dif/2](coroutining.html#dif/2). The implementation can deal with cyclic terms in `X` and `Y`.

The [when/2](coroutining.html#when/2) predicate is realised using attributed variables associated with the module `when`. It is defined in the autoload library `library(when)`.

**call_residue_vars**(`:Goal, -Vars`)  
Find residual attributed variables left by `Goal`. This predicate is intended for reasoning about and debugging programs that use coroutining or constraints. To see why this predicate is necessary, consider a predicate that poses contradicting constraints on a variable, and where that variable does not appear in any argument of the predicate and hence does not yield any residual goals on the toplevel when the predicate is invoked. Such programs should fail, but sometimes succeed because the constraint solver is too weak to detect the contradiction. Ideally, delayed goals and constraints are all executed at the end of the computation. The meta predicate [call_residue_vars/2](coroutining.html#call_residue_vars/2) finds variables that are given attributes or whose attributes are modified by `Goal`, regardless of whether or not these variables are reachable from the arguments of `Goal`.^(195The implementation of [call_residue_vars/2](coroutining.html#call_residue_vars/2) is completely redone in version 7.3.2 (7.2.1) after discussion with Bart Demoen. The current implementation no longer performs full scans of the stacks. The overhead is proportional to the number of attributed variables on the stack, dead or alive.).
