
# 8 Constraint Logic Programming

This chapter describes the extensions primarily designed to support **constraint logic programming** (CLP), an important declarative programming paradigm with countless practical applications.

CLP(`X`) stands for constraint logic programming over the domain `X`. Plain Prolog can be regarded **as~CLP(`H`)**, where `H` stands for *Herbrand terms*. Over this domain, [=/2](compare.html#=/2) and [dif/2](coroutining.html#dif/2) are the most important constraints that express, respectively, equality and disequality of terms. Plain Prolog can thus be regarded as a special case of CLP.

There are dedicated constraint solvers for several important domains:

- CLP(FD) for **integers** ([section A.9](clpfd.html#sec:A.9))
- CLP(B) for **Boolean** variables ([section A.8](clpb.html#sec:A.8))
- CLP(Q) for **rational** numbers ([section A.10](clpqr.html#sec:A.10))
- CLP(R) for **floating point** numbers ([section A.10](clpqr.html#sec:A.10)).

In addition, CHR ([chapter 9](chr.html#sec:9)) provides a general purpose constraint handling language to reason over user-defined constraints.

Constraints blend in naturally into Prolog programs, and behave exactly like plain Prolog predicates in those cases that can also be expressed without constraints. However, there are two key differences between constraints and plain Prolog predicates:

- Constraints can *delay* checks until their truth can be safely decided. This feature can significantly increase the *generality* of your programs, and preserves their relational nature.
- Constraints can take into account everything you state about the entities you reason about, independent of the order in which you state it, both *before* and also *during* any search for concrete solutions. Using available information to prune parts of the search space is called constraint *propagation*, and it is performed automatically by the available constraint solvers for their respective domains. This feature can significantly increase the *performance* of your programs.

Due to these two key advantages over plain Prolog, CLP has become an extremely important declarative programming paradigm in practice.

Among its most important and typical instances is CLP(FD), constraint logic programming over *integers*. For example, using constraints, you can state in the most general way that a variable `X` is an integer greater than 0. If, later, `X` is bound to a concrete integer, the constraint solver automatically ensures this. If you in addition constrain `X` to integers less than 3, the constraint solver combines the existing knowledge to infer that `X` is either 1 or 2 (see below). To obtain concrete values for `X`, you can ask the solver to *label* `X` and produce 1 and 2 on backtracking. See [section A.9](clpfd.html#sec:A.9).

``` code
?- use_module(library(clpfd)).
...
true.

?- X #> 0, X #< 3.
X in 1..2.

?- X #> 0, X #< 3, indomain(X).
X = 1 ;
X = 2.
```

Contrast this with plain Prolog, which has no efficient means to deal with (integer) `X > 0` and `X < 3`. At best it could translate `X > 0` to `between(1, infinite, X)` and a similar primitive for `X < 3`. If the two are combined it has no choice but to generate and test over this infinite two-dimensional space.

Using constraints therefore makes your program more *declarative* in that it frees you from some procedural aspects and limitations of Prolog.

When working with constraints, keep in mind the following:

- As with plain Prolog, [!/0](control.html#!/0) also destroys the declarative semantics of constraints. A cut after a goal that is delayed may prematurely prune the search space, because the truth of delayed goals is not yet established. There are several ways to avoid cuts in constraint logic programs, retaining both generality and determinism of your programs. See for example [zcompare/3](clpfd.html#zcompare/3).
- Term-copying operations ([assertz/1](db.html#assertz/1), [retract/1](db.html#retract/1), [findall/3](allsolutions.html#findall/3), [copy_term/2](manipterm.html#copy_term/2), etc.) generally also copy constraints. The effect varies from ok, silent copying of huge constraint networks to violations of the internal consistency of constraint networks. As a rule of thumb, copying terms holding attributes must be deprecated. If you need to reason about a term that is involved in constraints, use [copy_term/3](attvar.html#copy_term/3) to obtain the constraints as Prolog goals, and use these goals for further processing.

All of the mentioned constraint solvers are implemented using the attributed variables interface described in [section 8.1](attvar.html#sec:8.1). These are lower-level predicates that are mainly intended for library authors, not for typical Prolog programmers.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[8.1 Attributed variables](attvar.html)

[8.1.1 Attribute manipulation predicates](attvar.html#sec:8.1.1)

[8.1.2 Attributed variable hooks](attvar.html#sec:8.1.2)

[8.1.3 Operations on terms with attributed variables](attvar.html#sec:8.1.3)

[8.1.4 Special purpose predicates for attributes](attvar.html#sec:8.1.4)

[8.2 Coroutining](coroutining.html)
