
## 4.6 Comparison and Unification of Terms

Although unification is mostly done implicitly while matching the head of a predicate, it is also provided by the predicate =/2.

\[ISO\]`?Term1` **=** `?Term2`  
Unify `Term1` with `Term2`. True if the unification succeeds. It acts as if defined by the following fact:

``` code
=(Term, Term).
```

For behaviour on cyclic terms see the Prolog flag [occurs_check](flags.html#flag:occurs_check). Calls to [=/2](compare.html#=/2) in a clause body are compiled and may be (re)moved depending on the Prolog flag [optimise_unify](flags.html#flag:optimise_unify). See also [section 2.17.3](jitindex.html#sec:2.17.3).

\[ISO\]`@Term1` **\\** `@Term2`  
Equivalent to `\+``Term1 = Term2`.

This predicate is logically sound if its arguments are sufficiently instantiated. In other cases, such as `?- X ``\=`` Y.`, the predicate fails although there are solutions. This is due to the incomplete nature of [\\/1](control.html#\+/1).

To make your programs work correctly also in situations where the arguments are not yet sufficiently instantiated, use [dif/2](coroutining.html#dif/2) instead.

### 4.6.1 Standard Order of Terms

Comparison and unification of arbitrary terms. Terms are ordered in the so-called “standard order” . This order is defined as follows:

1.  `Variables`` < ``Numbers`` < ``Strings`` < ``Atoms`` < ``Compound Terms`
2.  Variables are sorted by address.
3.  `Numbers` are compared by value. Mixed rational/float are compared using [cmpr/2](arith.html#f-cmpr/2).^(68Up to 9.1.4, comparison was done as float.) NaN is considered smaller than all numbers, including `-inf`. If the comparison is equal, the float is considered the smaller value. If the Prolog flag [iso](flags.html#flag:iso) is defined, all floating point numbers precede all rationals.
4.  `Strings` are compared alphabetically.
5.  `Atoms` are compared alphabetically.
6.  `Compound` terms are first checked on their arity, then on their functor name (alphabetically) and finally recursively on their arguments, leftmost argument first.

Although variables are ordered, there are some unexpected properties one should keep in mind when relying on variable ordering. This applies to the predicates below as to predicate such as [sort/2](builtinlist.html#sort/2) as well as libraries that reply on ordering such as library `library(assoc)` and library `library(ordsets)`. Obviously, an established relation `A` `@<` `B` no longer holds if `A` is unified with e.g., a number. Also unifying `A` with `B` invalidates the relation because they become equivalent (==/2) after unification.

As stated above, variables are sorted by address, which implies that they are sorted by‘age’, where‘older’variables are ordered before‘newer’variables. If two variables are unified their‘shared’age is the age of oldest variable. This implies we can examine a list of sorted variables with‘newer’(fresh) variables without invalidating the order. Attaching an *attribute*, see [section 8.1](attvar.html#sec:8.1), turns an‘old’variable into a‘new’one as illustrated below. Note that the first always succeeds as the first argument of a term is always the oldest. This only applies for the *first* attribute, i.e., further manipulation of the attribute list does *not* change the‘age’.

``` code
?- T = f(A,B), A @< B.
T = f(A, B).

?- T = f(A,B), put_attr(A, name, value), A @< B.
false.
```

The above implies you *can* use e.g., an assoc (from library `library(assoc)`, implemented as an AVL tree) to maintain information about a set of variables. You must be careful about what you do with the attributes though. In many cases it is more robust to use attributes to register information about variables.

Note that the standard order is not well defined on *rational trees*, also known as *cyclic terms*. This [issue was identified](https://swi-prolog.discourse.group/t/how-to-compare-3-without-surprises-on-non-ground-terms/6386/42=jan) by Mats Carlsson, quoted below. See also [issue#1162 on GitHub](https://github.com/SWI-Prolog/swipl-devel/issues/1162).

> Consider the terms `A` and `B` defined by the equations
>
> ``` code
> [1] A = s(B,0).
> [2] B = s(A,1).
> ```
>
> - Clearly, `A` and `B` are not identical, so either `A @< B` or `A @> B` must hold.
> - Assume `A @< B`. But then, `s(A,1) @> s(B,0)` i.e., `B @< A`. Contradiction.
> - Assume `A @> B`. But then, `s(A,1) @< s(B,0)` i.e., `B @< A`. Contradiction.

\[ISO\]`@Term1` **==** `@Term2`  
True if `Term1` is equivalent to `Term2`. A variable is only identical to a sharing variable.

\[ISO\]`@Term1` **\\=** `@Term2`  
Equivalent to `\+``Term1 == Term2`.

\[ISO\]`@Term1` **@\<** `@Term2`  
True if `Term1` is before `Term2` in the standard order of terms.

\[ISO\]`@Term1` **@=\<** `@Term2`  
True if both terms are equal ([==/2](compare.html#==/2)) or `Term1` is before `Term2` in the standard order of terms.

\[ISO\]`@Term1` **@\>** `@Term2`  
True if `Term1` is after `Term2` in the standard order of terms.

\[ISO\]`@Term1` **@\>=** `@Term2`  
True if both terms are equal ([==/2](compare.html#==/2)) or `Term1` is after `Term2` in the standard order of terms.

\[ISO\]**compare**(`?Order, @Term1, @Term2`)  
Determine or test the `Order` between two terms in the standard order of terms. `Order` is one of `<`, `>` or `=`, with the obvious meaning.

### 4.6.2 Special unification and comparison predicates

This section describes special purpose variations on Prolog unification. The predicate [unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2) provides sound unification and is part of the ISO standard. The predicate [subsumes_term/2](compare.html#subsumes_term/2) defines‘one-sided unification’and is part of the ISO proposal established in Edinburgh (2010). Finally, [unifiable/3](compare.html#unifiable/3) is a‘what-if’version of unification that is often used as a building block in constraint reasoners.

\[ISO\]**unify_with_occurs_check**(`+Term1, +Term2`)  
As [=/2](compare.html#=/2), but using *sound unification*. That is, a variable only unifies to a term if this term does not contain the variable itself. To illustrate this, consider the two queries below.

``` code
1 ?- A = f(A).
A = f(A).
2 ?- unify_with_occurs_check(A, f(A)).
false.
```

The first statement creates a *cyclic term*, also called a *rational tree*. The second executes logically sound unification and thus fails. Note that the behaviour of unification through [=/2](compare.html#=/2) as well as implicit unification in the head can be changed using the Prolog flag [occurs_check](flags.html#flag:occurs_check).

The SWI-Prolog implementation of [unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2) is cycle-safe and only guards against *creating* cycles, not against cycles that may already be present in one of the arguments. This is illustrated in the following two queries:

``` code
?- X = f(X), Y = X, unify_with_occurs_check(X, Y).
X = Y, Y = f(Y).
?- X = f(X), Y = f(Y), unify_with_occurs_check(X, Y).
X = Y, Y = f(Y).
```

Some other Prolog systems interpret [unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2) as if defined by the clause below, causing failure on the above two queries. Direct use of [acyclic_term/1](typetest.html#acyclic_term/1) is portable and more appropriate for such applications.

``` code
unify_with_occurs_check(X,X) :- acyclic_term(X).
```

`+Term1` **=@=** `+Term2`  
True if `Term1` is a *variant* of (or *structurally equivalent* to) `Term2`. Testing for a variant is weaker than equivalence ([==/2](compare.html#==/2)), but stronger than unification ([=/2](compare.html#=/2)). Two terms `A` and `B` are variants iff there exists a renaming of the variables in `A` that makes `A` equivalent (==) to `B` and vice versa.^(69Row 7 and 8 of this table may come as a surprise, but row 8 is satisfied by (left-to-right) `A→C`, `B→A` and (right-to-left) `C→A`, `A→B`. If the same variable appears in different locations in the left and right term, the variant relation can be broken by consistent binding of both terms. E.g., after binding the first argument in row 8 to a value, both terms are no longer variant.) Examples:

> |     |                     |       |
> |----:|:-------------------:|:-----:|
> |   1 |      `a =@= A`      | false |
> |   2 |      `A =@= B`      | true  |
> |   3 | `x(A,A) =@= x(B,C)` | false |
> |   4 | `x(A,A) =@= x(B,B)` | true  |
> |   5 | `x(A,A) =@= x(A,B)` | false |
> |   6 | `x(A,B) =@= x(C,D)` | true  |
> |   7 | `x(A,B) =@= x(B,A)` | true  |
> |   8 | `x(A,B) =@= x(C,A)` | true  |

A term is always a variant of a copy of itself. Term copying takes place in, e.g., [copy_term/2](manipterm.html#copy_term/2), [findall/3](allsolutions.html#findall/3) or proving a clause added with [asserta/1](db.html#asserta/1). In the pure Prolog world (i.e., without attributed variables), [=@=/2](compare.html#=@=/2) behaves as if defined below. With attributed variables, variant of the attributes is tested rather than trying to satisfy the constraints.

``` code
A =@= B :-
        copy_term(A, Ac),
        copy_term(B, Bc),
        numbervars(Ac, 0, N),
        numbervars(Bc, 0, N),
        Ac == Bc.
```

The SWI-Prolog implementation is cycle-safe and can deal with variables that are shared between the left and right argument. Its performance is comparable to [==/2](compare.html#==/2), both on success and (early) failure. ^(70The current implementation is contributed by Kuniaki Mukai.)

This predicate is known by the name [variant/2](terms.html#variant/2) in some other Prolog systems. Be aware of possible differences in semantics if the arguments contain attributed variables or share variables.^(71In many systems variant is implemented using two calls to [subsumes_term/2](compare.html#subsumes_term/2).)

`+Term1` **\\@=** `+Term2`  
Equivalent to `‘``\+``Term1 =@= Term2’`. See [=@=/2](compare.html#=@=/2) for details.

\[ISO\]**subsumes_term**(`@Generic, @Specific`)  
True if `Generic` can be made equivalent to `Specific` by only binding variables in `Generic`. The current implementation performs the unification and ensures that the variable set of `Specific` is not changed by the unification. On success, the bindings are undone.^(72This predicate is often named [subsumes_chk/2](terms.html#subsumes_chk/2) in older Prolog dialects. The current name was established in the ISO WG17 meeting in Edinburgh (2010). The `chk` postfix was considered to refer to determinism as in e.g., [memberchk/2](builtinlist.html#memberchk/2).) This predicate respects constraints.

See [section 5.6](ssu.html#sec:5.6) for defining clauses whose head is unified using *single sided unification*.

**term_subsumer**(`+Special1, +Special2, -General`)  
`General` is the most specific term that is a generalisation of `Special1` and `Special2`. The implementation can handle cyclic terms.

**unifiable**(`@X, @Y, -Unifier`)  
If `X` and `Y` can unify, unify `Unifier` with a list of `Var` = `Value`, representing the bindings required to make `X` and `Y` equivalent.^(73This predicate was introduced for the implementation of [dif/2](coroutining.html#dif/2) and [when/2](coroutining.html#when/2) after discussion with Tom Schrijvers and Bart Demoen. None of us is really happy with the name and therefore suggestions for a new name are welcome.) This predicate can handle cyclic terms. Attributed variables are handled as normal variables. Associated hooks are *not* executed.

**?=**(`@Term1, @Term2`)  
Succeeds if the syntactic equality of `Term1` and `Term2` can be decided safely, i.e. if the result of `Term1 == Term2` will not change due to further instantiation of either term. It behaves as if defined by `?=(X,Y) :- \+ unifiable(X,Y,[_|_]).`
