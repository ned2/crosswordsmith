
## 4.31 Forall

\[semidet\]**forall**(`:Cond, :Action`)  
For all alternative bindings of `Cond`, `Action` can be proven. The example verifies that all arithmetic statements in the given list are correct. It does not say which is wrong if one proves wrong.

``` code
?- forall(member(Result = Formula, [2 = 1 + 1, 4 = 2 * 2]),
                 Result =:= Formula).
```

The predicate [forall/2](forall2.html#forall/2) is implemented as `\+ ( Cond, \+ Action)`, i.e., *There is no instantiation of `Cond` for which `Action` is false.*. The use of double negation implies that [forall/2](forall2.html#forall/2) *does not change any variable bindings*. It proves a relation. The [forall/2](forall2.html#forall/2) control structure can be used for its side-effects. E.g., the following asserts relations in a list into the dynamic database:

``` code
?- forall(member(Child-Parent, ChildPairs),
          assertz(child_of(Child, Parent))).
```

Using [forall/2](forall2.html#forall/2) as `forall(Generator, SideEffect)` is preferred over the classical *failure driven loop* as shown below because it makes it explicit which part of the construct is the generator and which part creates the side effects. Also, unexpected failure of the side effect causes the construct to fail. Failure makes it evident that there is an issue with the code, while a failure driven loop would succeed with an erroneous result.

``` code
        ...,
        (   Generator,
            SideEffect,
            fail
        ;   true
        )
```

If your intent is to create variable bindings, the [forall/2](forall2.html#forall/2) control structure is inadequate. Possibly you are looking for [maplist/2](apply.html#maplist/2), [findall/3](allsolutions.html#findall/3) or [foreach/2](aggregate.html#foreach/2).
