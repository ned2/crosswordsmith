
## 5.6 Single Sided Unification rules

For the execution of a normal Prolog clause, the goal term is unified with the head of the clause. This allows us to write facts such as below and use this relation in all four possible *modes*. This is the basis of SLD resolution that turns Prolog into a logic programming language.

``` code
parent('Bob', 'Susan').
```

In practice though, Prolog is both a logic programming language and a language for expressing computations in a near *procedural* style. The first is used to solve (notably) combinatorial problems while the latter is used for I/O, data transformation and the many non-logical operations that are involved in many applications.

Many Prolog programmers experience writing procedural style Prolog as fighting non-determinism and dealing with hard to debug silent failures because no clause matches some goal. Below are two typical queries on library predicates that have a procedural nature, i.e., are *single moded*.

``` code
?- sum_list(a, X).
false.

?- sum_list([1|T], X).
T = [],
X = 1 ;
ERROR: Arguments are not sufficiently instantiated
```

The definition of [sum_list/2](lists.html#sum_list/2) is it appears in library(lists) is below. This implementation can be considered elegant. Note that [sum_list/2](lists.html#sum_list/2) has only one meaningful mode: (+,-). A general (logical) implementation would allow for a partial list or a list holding one or more variables, With a proper list that holds a single variable we can still make a sound logical implementation. In all other cases the number of solutions is infinite and even *uncountable* for a partial list, making the predicate useless as a *generator* of solutions.

``` code
sum_list(Xs, Sum) :-
    sum_list(Xs, 0, Sum).

sum_list([], Sum, Sum).
sum_list([X|Xs], Sum0, Sum) :-
    Sum1 is Sum0 + X,
    sum_list(Xs, Sum1, Sum).
```

If we want to avoid the above dubious behaviour we have two options. First, we can verify that the first argument is a list before entering the recursion, changing the first clause as below. The disadvantage is that we process the list twice.

``` code
sum_list(Xs, Sum) :-
    must_be(list, Xs),
    sum_list(Xs, 0, Sum).
```

Alternatively, we can rewrite the second clause to verify the list on the fly. That leads to the code below. Most likely the overhead of this alternative compared to the above is even worse in many Prolog implementations. Most people would also consider this code rather inelegant.

``` code
sum_list(Var, _, _) :-
    var(Var),
    instantiation_error(Var).
sum_list([], Sum, Sum) :-
    !.
sum_list([X|Xs], Sum0, Sum) :-
    !,
    Sum1 is Sum0 + X,
    sum_list(Xs, Sum1, Sum).
sum_list(NoList, _, _) :-
    type_error(list, NoList).
```

Another example is a relation max/3 , expressing the maximum of two numbers. A classical textbook definition could be as below. This code has two drawbacks. First it leaves an open choice points in most Prolog implementations if `X` is the largest and second it compares the two numbers twice. Some Prolog systems detect this particular case, but in general it needs two know that one test is the strict negation of the other.

``` code
max(X,Y,X) :- X >= Y.
max(X,Y,Y) :- Y > X.
```

As a result people use a cut and might come up with the **wrong** solution below. Consider the query `?- max(5,2,2).` to see why this code is broken.

``` code
max(X,Y,X) :- X >= Y, !.
max(_,Y,Y).
```

A correct solution is below, *delaying* binding the output until after the cut.

``` code
max(X,Y,M) :- X >= Y, !, M = X.
max(_,Y,Y).
```

Some people may prefer using if-then-else as below. This is arguable the cleanest efficient solution in standard Prolog.

``` code
max(X,Y,M) :- ( X >= Y -> M = X ; M = Y ).
```

As we have seen from these examples, writing procedural code in Prolog requires us to follow the two basic principles below. Both principles have been properly described in *The Craft of Prolog* [O'Keefe, 1990](Bibliography.html#Keefe:90).

- Structure every clause as `Head :- Guard, !, Body`. Every clause has the cut as early as possible. `Guard` can be empty. The last clause often does not need a cut.
- Avoid that the head unification binds values in the goal term. We see this may lead to undesirable results such as sum_list(L,S) binding `L` to‘\[\]\` and `S` to‘0\` as well as loss of *steadfastness*, causing max(5,2,2) to succeed. The first requires additional [var/1](typetest.html#var/1) or [nonvar/1](typetest.html#nonvar/1) tests. The second requires delaying unification until after the cut.

[Picat](http://picat-lang.org/) provides the **=\>/2** alternative for the Prolog *neck* (**:-/2**) to force the above practices. A Picat rule has the following shape:

``` code
Head, Guard => Body.
```

This is semantically equivalent to the Prolog clause below. The [subsumes_term/2](compare.html#subsumes_term/2) guarantees the clause head is more *generic* than the goal term and thus unifying the two does not affect any of the arguments of the goal. This implies all output unification \_must\_ be done after the head unification.

``` code
p(V1,V2,...,Vn) :-
    Pattern = p(A1,A2,...,An),
    Args = p(V1,V2,...,Vn),
    subsumes_term(Pattern, Args),
    Pattern = Args,
    Guard,
    !,
    Body.
```

SWI-Prolog as of version 8.3.19 support **=\>/2** as an alternative to normal Prolog clauses. The construct comes with the following properties.

- A predicate either uses **:-/2** for all its clauses or **=\>/2**. Mixing is **not** allowed and raises a permission error for a clause that does not use the same *neck* as the first clause.
- Unlike Picat, it is an error if no clause matches.

Given **=\>/2** rules, we can rewrite [sum_list/2](lists.html#sum_list/2) as below. The first clause can be written using **:-/2** or **=\>/2**. As the head is the most general head and there is only one clause these are fully equivalent. The sum_list/3 helper needs a small modification: we need to delay the unification against `Sum` to the body. The last clause is equivalent.

``` code
sum_list(Xs, Sum) =>
    sum_list(Xs, 0, Sum).

sum_list([], Sum0, Sum) =>
    Sum = Sum0.
sum_list([X|Xs], Sum0, Sum) =>
    Sum1 is Sum0 + X,
    sum_list(Xs, Sum1, Sum).
```

Given this definition, `sum_list(L,S)` no longer matches a rule and neither does e.g., `sum_list(a,S)`. Both raise an error. Currently the error is defined as below.

``` code
existence_error(matching_rule, Head)
```

Should silent failure be desired if no rule matches, this is easily encoding by adding a rule at the end using the most general head and [fail/0](control.html#fail/0) as body:

``` code
sum_list(_,_,_) => fail.
```

### 5.6.1 Single Sided Unification Guards

Using the construction `Head, Guard => Body`, the `Guard` is executed *after* the single sided head unification. If the `Guard` succeeds the clause executes a cut ([!/0](control.html#!/0)) and proceeds normally. There are no restrictions on the guard code. A well behaved guard is a *test*. Notably:

- Though not enforced^(182We do not know about an efficient way to enforce unification against head arguments), guard code shall not instantiate variables in the `Head` because this breaks the promise of SSU and may make the node non-steadfast.
- It is bad style (but again, not enforced) to have any type of side effects (output, database change, etc.)
- Typically, guard calls are `semidet`. Non-deterministic calls are allowed. If the guard succeeds with choicepoints these are pruned before the body is entered.

As a special exception, explicit unification against a variable in the head is moved into the head. See [section 2.17.3](jitindex.html#sec:2.17.3). In the example below, the `X = f(I)` is moved into the head and (thus) is executed using single sided unification.

``` code
p(X), X = f(I), integer(I) => q(X).
```

> **Warning** Moving the guard unification into the head changes the semantics of the unification. This may be defended by the rules above that claim one should not unify against the head arguments in the guard. Future versions may use a dedicated operator to indicate that the unification may be moved into the head.

### 5.6.2 Consequences of `=>` single sided unification rules

The **=\>/2** construct is handled by the low-level compiler if no *guard* is present. If a guard is present it is currently compiled into the construct below. The Picat **?=\>/2** *neck* operator is like **=\>/2**, but does not *commit* to this rule. We are not yet sure whether or not SWI-Prolog will remain supporting **?=\>/2**.^(183**?=\>/2** is currently implemented but not defined as an operator.)

``` code
Head ?=> Guard, !, Body.
```

The main consequence is that [clause/2](examineprog.html#clause/2) cannot distinguish between a normal clause and a **=\>/2** clause. In the current implementation it operates on both without distinguishing the two. This implies e.g., *cross referencing* still works. Meta interpretation however does not work. In future versions [clause/2](examineprog.html#clause/2) may fail on these rules. As an alternative we provide [rule/2](ssu.html#rule/2),3.

**rule**(`:Head, -Rule`)  
**rule**(`:Head, -Rule, -Ref`)  
True when `Rule` is a rule/clause that implements `Head`. `Rule` is a complete rule term. For a normal clause this is a term `Head :- Body` and for a single sided unification rule it is a term `Head ``=>`` Body`.

### 5.6.3 Single sided unification for Definite Clause Grammars

Single sided unification is attractive for *generative DCG rules*, i.e., DCG rules that are used to *serialize* some term. In that context they avoid unwanted matching on variables and provide better error messages in case not all possible terms are described by the grammar. Single sided unification has no practical use for parsing because the arguments are typically *output* arguments.

If the head of an SSU DCG rules is a term `Head, Extra`, `Extra` is interpreted as a *push back list* if it is a list and as an SSU *guard* otherwise. The guard is *not* subject to DCG expansion, i.e., it is interpreted as if enclosed by `{}`.

### 5.6.4 SSU: Future considerations

The current implementation is a rather simple. Single sided unification is achieved doing normal head unification and backtrack if this unification bound variables in the goal term. Future versions are likely to backtrack as soon as we find a variable in the goal that needs to be unified.

It is likely that in due time significant parts of the libraries will be migrated to use SSU rules, turning many silent failures on type errors into errors.
