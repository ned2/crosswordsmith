
## 4.40 Debugging and declaring determinism

A common issue with Prolog programs of a *procedural* nature is to guarantee deterministic behaviour and debug possible problems with determinism. SWI-Prolog provides several mechanisms to make writing, debugging and maintaining deterministic code easier. One of them is *Single Sided Unification* using **=\>/2** rules as described in [section 5.6](ssu.html#sec:5.6). This section deals with annotating your program.

If a program does not behave according to these annotations it raises an **error/2** exception where the formal term is `determinism_error(Pred, Declared, Observed, DeclType)`, where `Declared` is currently always `det`, `Observed` is one of `fail` or `nondet` and `DeclType` is one of `property` ([det/1](debug-determinism.html#det/1)), `guard` ([\$/0](debug-determinism.html#$/0)) or `goal` ([\$/1](debug-determinism.html#$/1)). Using [trap/1](prologdebug.html#trap/1) or gtrap/1 we can ask Prolog to start the debugger in such events using

``` code
?- gtrap(determinism_error(_,_,_,_)).
```

> **WARNING:** The primitives in this section are experimental. The naming and exact semantics may change. If you are interested in this, please follow and contribute to discussion on the Discourse forum.

\[directive\]**det**(`+PredicateIndicators`)  
Declare a number of predicates as `det` (*deterministic*). As a result, both failure and success with a choicepoint is considered an error. The behaviour if the declaration is violated is controlled with the Prolog flag [determinism_error](flags.html#flag:determinism_error). The default is to raise an exception (`error`). Consider the following program:

``` code
:- det(p/1).

p(1).
p(2).
```

Now, a call `?- p(1).` behaves normally. However:

``` code
?- p(X).
ERROR: Deterministic procedure p/1 succeeded with a choicepoint
ERROR: In:
ERROR:   [10] p(1)

?- p(a).
ERROR: Deterministic procedure p/1 failed
ERROR: In:
ERROR:   [10] p(a)
```

Violations throw an **error/2** exception `determinism_error(Pred, Declared, Observed, property)`.

The [trap/1](prologdebug.html#trap/1) (cli) or gtrap/1 (gui) predicate can be used to make the debugger stop near the error. For example:

``` code
?- gtrap(determinism_error(_,_,_,_)).
```

\[experimental\]**\$**  
The [\$/0](debug-determinism.html#$/0) constructs acts similar to the [!/0](control.html#!/0), but in addition declares that the remainder of the clause body shall succeed deterministically. It exploits the same underlying mechanism as the [det/1](debug-determinism.html#det/1) declaration. See also [\$/1](debug-determinism.html#$/1).

Violations throw an **error/2** exception `determinism_error(Pred, Declared, Observed, guard)`.

\[experimental\]**\$**(`:Goal`)  
Verify that `Goal` succeeds deterministically. This predicate has no effect if `Goal` succeeds without a choicepoint. Otherwise the result depends on the Prolog flag [determinism_error](flags.html#flag:determinism_error):

**silent**  
Act as [once/1](metacall.html#once/1).

**warning**  
Print a warning and act as [once/1](metacall.html#once/1).

**error**  
Raise a `determinism_error` exception.

Note that if [\$/1](debug-determinism.html#$/1) is used for the last call, last call optimization is not effective. This behaviour ensures consistent errors or warnings. Last call optimization with determinism checking can be realised using `..., $, Last.`, i.e. by executing [\$/0](debug-determinism.html#$/0) before the last call rather than wrapping the last call in [\$/1](debug-determinism.html#$/1).

Violations throw an **error/2** exception `determinism_error(Pred, Declared, Observed, goal)`.

A deterministic predicate may call normal predicates. No error is triggered as long as the deterministic predicate either ignores a possible failure, e.g., using [\\/1](control.html#\+/1) and prunes possible choice points created by called predicates. If the last predicate is a normal predicate the requirement to succeed deterministically is transferred to the new goal. As last-call optimization causes the information which predicate initially claimed to be deterministic to be lost, the error is associated with the called predicate. Debug mode (see [debug/0](debugger.html#debug/0) or the Prolog flag [debug](flags.html#flag:debug)) may be used to avoid last call optimization and find the call stack that causes the issue.
