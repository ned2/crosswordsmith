
## 4.21 Analysing and Constructing Terms

\[ISO\]**functor**(`?Term, ?Name, ?Arity`)  
True when `Term` is a term with functor `Name`/`Arity`. If `Term` is a variable it is unified with a new term whose arguments are all different variables (such a term is called a skeleton). If `Term` is atomic, `Arity` will be unified with the integer 0, and `Name` will be unified with `Term`. Raises `instantiation_error()` if `Term` is unbound and `Name`/`Arity` is insufficiently instantiated.

SWI-Prolog also supports terms with arity 0, as in **`a()`** (see [section 5](extensions.html#sec:5)). Such terms must be processed using [functor/4](manipterm.html#functor/4) or [compound_name_arity/3](manipterm.html#compound_name_arity/3). The predicate [functor/3](manipterm.html#functor/3) and [=../2](manipterm.html#=../2) raise a `domain_error` when faced with these terms. Without this precaution a *round trip* of a term with arity 0 over [functor/3](manipterm.html#functor/3) would create an atom.

**functor**(`?Term, ?Name, ?Arity, ?Type`)  
As [functor/3](manipterm.html#functor/3), but designed to work with zero-arity terms (e.g., **`a()`**, see [section 5](extensions.html#sec:5)). `Type` is one of `atom`, `compound`, `callable` or `atomic`. `Type` *must* be instantiated if `Name` is an atom and `Arity` is 0 (zero). In other cases `Type` may be a variable. This predicate is true if `Term` (either initially or after having been created from `Name` and `Type`) and `Type` are related as below

- If `Term` is compound (including zero-arity compounds), `Type` must be `compound` or `callable`. If `Type` is unbound is is unified with `compound`.
- If `Term` is an atom, `Type` must be `atom` or `callable`. If `Type` is unbound is is unified with `atom`.
- Else `Type` is unified with `atomic`.

This predicate provides a safe *round trip* for zero-arity compounds and atoms. It can also be used as a variant of [functor/3](manipterm.html#functor/3) that only processes compound or callable terms. See also [compound/1](typetest.html#compound/1), [callable/1](typetest.html#callable/1) and [compound_name_arity/3](manipterm.html#compound_name_arity/3).

\[ISO\]**arg**(`?Arg, +Term, ?Value`)  
`Term` should be instantiated to a term, `Arg` to an integer between 1 and the arity of `Term`. `Value` is unified with the `Arg`-th argument of `Term`. `Arg` may also be unbound. In this case `Value` will be unified with the successive arguments of the term. On successful unification, `Arg` is unified with the argument number. Backtracking yields alternative solutions.^(115The instantiation pattern (-, +, ?) is an extension to‘standard’Prolog. Some systems provide genarg/3 that covers this pattern.) The predicate [arg/3](manipterm.html#arg/3) fails silently if `Arg`` = 0` or `Arg`` > `*`arity`* and raises the exception `domain_error(not_less_than_zero, ``Arg``)` if `Arg`` < 0`.

\[ISO\]`?Term` **=..** `?List`  
`List` is a list whose head is the functor of `Term` and the remaining arguments are the arguments of the term. Either side of the predicate may be a variable, but not both. This predicate is called‘Univ’.

``` code
?- foo(hello, X) =.. List.
List = [foo, hello, X]

?- Term =.. [baz, foo(1)].
Term = baz(foo(1))
```

SWI-Prolog also supports terms with arity 0, as in **`a()`** (see [section 5](extensions.html#sec:5)). Such terms must be processed using [compound_name_arguments/3](manipterm.html#compound_name_arguments/3). This predicate raises a domain error as shown below. See also [functor/3](manipterm.html#functor/3).

``` code
?- a() =.. L.
ERROR: Domain error: `compound_non_zero_arity' expected, found `a()'
```

**compound_name_arity**(`?Compound, ?Name, ?Arity`)  
Version of [functor/3](manipterm.html#functor/3) that only works for compound terms and can examine and create compound terms with zero arguments (e.g, **`name()`**). See also [compound_name_arguments/3](manipterm.html#compound_name_arguments/3). See also [functor/4](manipterm.html#functor/4).

**compound_name_arguments**(`?Compound, ?Name, ?Arguments`)  
Rationalized version of [=../2](manipterm.html#=../2) that can compose and decompose compound terms with zero arguments. See also [compound_name_arity/3](manipterm.html#compound_name_arity/3).

**numbervars**(`+Term, +Start, -End`)  
Unify the free variables in `Term` with a term `$VAR(N)`, where `N` is the number of the variable. Counting starts at `Start`. `End` is unified with the number that should be given to the next variable.^(bugOnly *tagged integers* are supported (see the Prolog flag [max_tagged_integer](flags.html#flag:max_tagged_integer)). This suffices to count all variables that can appear in the largest term that can be represented, but does not support arbitrary large integer values for `Start`. On overflow, a `representation_error(tagged_integer)` exception is raised.) The example below illustrates this. Note that the toplevel prints `'$VAR'(0)` as `A` due to the `numbervars(true)` option used to print answers.

``` code
?- Term = f(X,Y,X),
   numbervars(Term, 0, End, [singleton(true)]),
   write_canonical(Term), nl.
f('$VAR'(0),'$VAR'('_'),'$VAR'(0))
Term = f(A, _, A),
X = A,
Y = B,
End = 2.
```

See also the `numbervars` option to [write_term/3](termrw.html#write_term/3) and [numbervars/4](manipterm.html#numbervars/4).

**numbervars**(`+Term, +Start, -End, +Options`)  
As [numbervars/3](manipterm.html#numbervars/3), providing the following options:

**functor_name**(`+Atom`)  
Name of the functor to use instead of `$VAR`.

**attvar**(`+Action`)  
What to do if an attributed variable is encountered. Options are `skip`, which causes [numbervars/3](manipterm.html#numbervars/3) to ignore the attributed variable, `bind` which causes it to treat it as a normal variable and assign the next `'$VAR'`(N) term to it, or (default) `error` which raises a `type_error` exception.^(116This behaviour was decided after a long discussion between David Reitter, Richard O'Keefe, Bart Demoen and Tom Schrijvers.)

**singletons**(`+Bool`)  
If `true` (default `false`), [numbervars/4](manipterm.html#numbervars/4) does singleton detection. Singleton variables are unified with `'$VAR'('_')`, causing them to be printed as `_` by [write_term/2](termrw.html#write_term/2) using the numbervars option. This option is exploited by [portray_clause/2](listing.html#portray_clause/2) and [write_canonical/2](termrw.html#write_canonical/2).^(bugCurrently this option is ignored for cyclic terms.)

**var_number**(`@Term, -VarNumber`)  
True if `Term` is numbered by [numbervars/3](manipterm.html#numbervars/3) and `VarNumber` is the number given to this variable. This predicate avoids the need for unification with `'$VAR'(X)` and opens the path for replacing this valid Prolog term by an internal representation that has no textual equivalent.

\[ISO\]**term_variables**(`+Term, -List`)  
Unify `List` with a list of variables, each sharing with a unique variable of `Term`.^(117This predicate used to be called free_variables/2 . The name [term_variables/2](manipterm.html#term_variables/2) is more widely used. The old predicate is still available from the library `library(backcomp)`.) The variables in `List` are ordered in order of appearance traversing `Term` depth-first and left-to-right. See also [term_variables/3](manipterm.html#term_variables/3) and [nonground/2](manipterm.html#nonground/2). For example:

``` code
?- term_variables(a(X, b(Y, X), Z), L).
L = [X, Y, Z].
```

\[semidet\]**nonground**(`+Term, -Var`)  
True when `Var` is a variable in `Term`. Fails if `Term` is *ground* (see [ground/1](typetest.html#ground/1)). This predicate is intended for coroutining to trigger a wakeup if `Term` becomes ground, e.g., using [when/2](coroutining.html#when/2). The current implementation always returns the first variable in depth-first left-right search. Ideally it should return a random member of the set of variables (see [term_variables/2](manipterm.html#term_variables/2)) to realise logarithmic complexity for the ground trigger. Compatible with ECLiPSe and hProlog.

**term_variables**(`+Term, -List, ?Tail`)  
Difference list version of [term_variables/2](manipterm.html#term_variables/2). That is, `Tail` is the tail of the variable list `List`.

**term_singletons**(`+Term, -List`)  
Unify `List` with a list of variables, each sharing with a variable that appears only once in `Term`.^(bugIn the current implementation `Term` must be acyclic. If not, a `representation_error` is raised.) Note that, if a variable appears in a shared subterm, it is *not* considered singleton. Thus, `A` is *not* a singleton in the example below. See also the `singleton` option of [numbervars/4](manipterm.html#numbervars/4).

``` code

?- S = a(A), term_singletons(t(S,S), L).
L = [].
```

**is_most_general_term**(`@Term`)  
True if `Term` is a callable term where all arguments are non-sharing variables or `Term` is a list whose members are all non-sharing variables. This predicate is used to reason about call subsumption for tabling and is compatible with XSB. See also [subsumes_term/2](compare.html#subsumes_term/2). Examples:

> |     |                                |       |
> |----:|--------------------------------|:-----:|
> |   1 | `is_most_general_term(1)`      | false |
> |   2 | `is_most_general_term(p)`      | true  |
> |   3 | `is_most_general_term(p(_))`   | true  |
> |   4 | `is_most_general_term(p(_,a))` | false |
> |   5 | `is_most_general_term(p(X,X))` | false |
> |   6 | `is_most_general_term([])`     | true  |
> |   7 | `is_most_general_term([_|_])`  | false |
> |   8 | `is_most_general_term([_,_])`  | true  |
> |   9 | `is_most_general_term([X,X])`  | false |

\[ISO\]**copy_term**(`+In, -Out`)  
Create a version of `In` with renamed (fresh) variables and unify it to `Out`. Attributed variables (see [section 8.1](attvar.html#sec:8.1)) have their attributes copied. The implementation of [copy_term/2](manipterm.html#copy_term/2) can deal with infinite trees (cyclic terms). As pure Prolog cannot distinguish a ground term from another ground term with exactly the same structure, ground sub-terms are *shared* between `In` and `Out`. Sharing ground terms does affect [setarg/3](manipterm.html#setarg/3). SWI-Prolog provides [duplicate_term/2](manipterm.html#duplicate_term/2) to create a true copy of a term.

**copy_term**(`+VarsIn, +In, -VarsOut, -Out`)  
Similar to [copy_term/2](manipterm.html#copy_term/2), but only rename the variables in `VarsIn` that appear in `In`.^(118This predicate is based on a similar predicate in s(CASP) by Joaquin Arias.) Variables in `In` that do not appear in `VarsIn` are *shared* between `In` and `Out`. Sub terms that only contain such shared variables are shared as a whole between `In` and `Out`. `VarsIn` is often a list, but can be an arbitrary term. For example:

``` code
?- copy_term([X], q(X,Y), Vars, Term).
Vars = [_A],
Term = q(_A, Y).
```

Note that if `VarsIn` and `In` do not share any variables, `Out` is equivalent to `In` and `VarsOut` is a copy (as [copy_term/2](manipterm.html#copy_term/2)) of `VarsIn`. If `In` does not contain any variables not in `VarsIn` the result is the same as `copy_term(VarsIn-In, VarsOut-Out`).

**copy_term_nat**(`+VarsIn, +In, -VarsOut, -Out`)  
As [copy_term/4](manipterm.html#copy_term/4), using the attributed variable semantics of [copy_term_nat/2](attvar.html#copy_term_nat/2). This implies that attributed variables that appear in `VarsIn` appear as renamed plain variables in `VarsOut` and `Out`. Attributed variables in `In` that do *not* appear in `VarsIn` are shared between `In` and `Out`.

### 4.21.1 Non-logical operations on terms

Prolog is not able to *modify* instantiated parts of a term. Lacking that capability makes the language much safer, but unfortunately there are problems that suffer severely in terms of time and/or memory usage. Always try hard to avoid the use of these primitives, but they can be a good alternative to using dynamic predicates. See also [section 4.33](gvar.html#sec:4.33), discussing the use of global variables.

**setarg**(`+Arg, +Term, +Value`)  
Extra-logical predicate. Assigns the `Arg`-th argument of the compound term `Term` with the given `Value`. The assignment is undone if backtracking brings the state back into a position before the [setarg/3](manipterm.html#setarg/3) call. If the designated argument of `Term` is a variable, this variable is unified with `Value` using normal unification, i.e., [setarg/3](manipterm.html#setarg/3) behaves as [arg/3](manipterm.html#arg/3) in this case. Note that this may produce a cyclic term if `Value` contains this variable. See also [nb_setarg/3](manipterm.html#nb_setarg/3).

This predicate may be used for destructive assignment to terms, using them as an extra-logical storage bin. Always try hard to avoid the use of [setarg/3](manipterm.html#setarg/3) as it is not supported by many Prolog systems and one has to be very careful about unexpected copying as well as unexpected noncopying of terms. A good practice to improve somewhat on this situation is to make sure that terms whose arguments are subject to [setarg/3](manipterm.html#setarg/3) have one unused and unshared variable in addition to the used arguments. This variable avoids unwanted sharing in, e.g., [copy_term/2](manipterm.html#copy_term/2), and causes the term to be considered as non-ground. An alternative is to use [put_attr/3](attvar.html#put_attr/3) to attach information to attributed variables (see [section 8.1](attvar.html#sec:8.1)).

**nb_setarg**(`+Arg, +Term, +Value`)  
Assigns the `Arg`-th argument of the compound term `Term` with the given `Value` as [setarg/3](manipterm.html#setarg/3), but on backtracking the assignment is *not* reversed. If `Value` is not atomic, it is duplicated using [duplicate_term/2](manipterm.html#duplicate_term/2). This predicate uses the same technique as [nb_setval/2](gvar.html#nb_setval/2). We therefore refer to the description of [nb_setval/2](gvar.html#nb_setval/2) for details on non-backtrackable assignment of terms. This predicate is compatible with GNU-Prolog `setarg(A,T,V,false)`, removing the type restriction on `Value`. Below is an example for counting the number of solutions of a goal. Note that this implementation is thread-safe, reentrant and capable of handling exceptions. Realising these features with a traditional implementation based on assert/retract or [flag/3](db.html#flag/3) is much more complicated.

``` code
:- meta_predicate
        succeeds_n_times(0, -).

succeeds_n_times(Goal, Times) :-
        Counter = counter(0),
        (   Goal,
            arg(1, Counter, N0),
            N is N0 + 1,
            nb_setarg(1, Counter, N),
            fail
        ;   arg(1, Counter, Times)
        ).
```

See also [nb_linkarg/3](manipterm.html#nb_linkarg/3) and [foldall/4](aggregate.html#foldall/4).

**nb_linkarg**(`+Arg, +Term, +Value`)  
As [nb_setarg/3](manipterm.html#nb_setarg/3), but like [nb_linkval/2](gvar.html#nb_linkval/2) it does *not* duplicate `Value`. Use with extreme care and consult the documentation of [nb_linkval/2](gvar.html#nb_linkval/2) before use.

**duplicate_term**(`+In, -Out`)  
Version of [copy_term/2](manipterm.html#copy_term/2) that also copies ground terms and therefore ensures that destructive modification using [setarg/3](manipterm.html#setarg/3) does not affect the copy. See also [nb_setval/2](gvar.html#nb_setval/2), [nb_linkval/2](gvar.html#nb_linkval/2), [nb_setarg/3](manipterm.html#nb_setarg/3) and [nb_linkarg/3](manipterm.html#nb_linkarg/3).

\[semidet\]**same_term**(`@T1, @T2`)  
True if `T1` and `T2` are equivalent and will remain equivalent, even if [setarg/3](manipterm.html#setarg/3) is used on either of them. This means `T1` and `T2` are the same variable, equivalent atomic data or a compound term allocated at the same address.
