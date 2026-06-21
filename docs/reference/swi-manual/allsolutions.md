
## 4.30 Finding all Solutions to a Goal

\[ISO\]**findall**(`+Template, :Goal, -Bag`)  
Create a list of the instantiations `Template` gets successively on backtracking over `Goal` and unify the result with `Bag`. Succeeds with an empty list if `Goal` has no solutions.  

[findall/3](allsolutions.html#findall/3) is equivalent to [bagof/3](allsolutions.html#bagof/3) with all *free* variables appearing in `Goal` scoped to the `Goal` with an existential (caret) operator (**`^`**), except that [bagof/3](allsolutions.html#bagof/3) fails when `Goal` has no solutions.

**findall**(`+Template, :Goal, -Bag, +Tail`)  
As [findall/3](allsolutions.html#findall/3), but returns the result as the difference list `Bag`-`Tail`. The 3-argument version is defined as

``` code
findall(Templ, Goal, Bag) :-
        findall(Templ, Goal, Bag, [])
```

\[nondet\]**findnsols**(`+N, @Template, :Goal, -List`)  
\[nondet\]**findnsols**(`+N, @Template, :Goal, -List, ?Tail`)  
As [findall/3](allsolutions.html#findall/3) and [findall/4](allsolutions.html#findall/4), but generates at most `N` solutions. If `N` solutions are returned, this predicate succeeds with a choice point if `Goal` has a choice point. Backtracking returns the next chunk of (at most) `N` solutions. In addition to passing a plain integer for `N`, a term of the form `count(N)` is accepted. Using `count(N)`, the size of the next chunk can be controlled using [nb_setarg/3](manipterm.html#nb_setarg/3). The non-deterministic behaviour used to implement the *chunk* option in `library(pengines)`. Based on Ciao, but the Ciao version is deterministic. Portability can be achieved by wrapping the goal in [once/1](metacall.html#once/1). Below are three examples. The first illustrates standard chunking of answers. The second illustrates that the chunk size can be adjusted dynamically and the last illustrates that no choice point is left if `Goal` leaves no choice-point after the last solution.

``` code
?- findnsols(5, I, between(1, 12, I), L).
L = [1, 2, 3, 4, 5] ;
L = [6, 7, 8, 9, 10] ;
L = [11, 12].

?- State = count(2),
   findnsols(State, I, between(1, 12, I), L),
   nb_setarg(1, State, 5).
State = count(5), L = [1, 2] ;
State = count(5), L = [3, 4, 5, 6, 7] ;
State = count(5), L = [8, 9, 10, 11, 12].

?- findnsols(4, I, between(1, 4, I), L).
L = [1, 2, 3, 4].
```

\[ISO\]**bagof**(`+Template, :Goal, -Bag`)  
Unify `Bag` with the alternatives of `Template`. If `Goal` has free variables besides the one sharing with `Template`, [bagof/3](allsolutions.html#bagof/3) will backtrack over the alternatives of these free variables, unifying `Bag` with the corresponding alternatives of `Template`. The construct `+``Var``^``Goal` tells [bagof/3](allsolutions.html#bagof/3) not to bind `Var` in `Goal`. [bagof/3](allsolutions.html#bagof/3) fails if `Goal` has no solutions.

The example below illustrates [bagof/3](allsolutions.html#bagof/3) and the **`^`** operator. The variable bindings are printed together on one line to save paper.

``` code
2 ?- listing(foo).
foo(a, b, c).
foo(a, b, d).
foo(b, c, e).
foo(b, c, f).
foo(c, c, g).
true.

3 ?- bagof(C, foo(A, B, C), Cs).
A = a, B = b, C = G308, Cs = [c, d] ;
A = b, B = c, C = G308, Cs = [e, f] ;
A = c, B = c, C = G308, Cs = [g].

4 ?- bagof(C, A^foo(A, B, C), Cs).
A = G324, B = b, C = G326, Cs = [c, d] ;
A = G324, B = c, C = G326, Cs = [e, f, g].

5 ?-
```

\[ISO\]**setof**(`+Template, +Goal, -Set`)  
Equivalent to [bagof/3](allsolutions.html#bagof/3), but sorts the result using [sort/2](builtinlist.html#sort/2) to get a sorted list of alternatives without duplicates.
