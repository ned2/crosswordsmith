
## A.13 library(dcg/high_order): High order grammar operations

This library provides facilities comparable [maplist/3](apply.html#maplist/3), [ignore/1](metacall.html#ignore/1) and [foreach/2](aggregate.html#foreach/2) for DCGs.

**STATUS**: This library is experimental. The interface and implementation may change based on feedback. Please send feedback to the mailinglist or the issue page of the `swipl-devel.git` repository.

\[nondet\]**sequence**(`:Element, ?List`)`//`  
Match or generate a sequence of `Element`. This predicate is deterministic if `List` is fully instantiated and `Element` is deterministic. When parsing, this predicate is *gready* and does not prune choice points. For example:

``` code
?- phrase(sequence(digit, Digits), `123a`, L).
Digits = "123",
L = [97] ;
Digits = [49, 50],
L = [51, 97] ;
...
```

\[nondet\]**sequence**(`:Element, :Sep, ?List`)`//`  
Match or generate a sequence of `Element` where each pair of elements is separated by `Sep`. When *parsing*, a matched `Sep` *commits*. The final element is *not* committed. More formally, it matches the following sequence:

``` code
(Element, (Sep,Element)*)?
```

See also [sequence//5](highorder.html#sequence//5).

\[semidet\]**sequence**(`:Start, :Element, :Sep, :End, ?List`)`//`  
Match or generate a sequence of `Element` enclosed by `Start` end `End`, where each pair of elements is separated by `Sep`. More formally, it matches the following sequence:

``` code
Start, (Element, (Sep,Element)*)?, End
```

The example below matches a Prolog list of integers:

``` code
?- phrase(sequence(("[",blanks),
                   number, (",",blanks),
                   (blanks,"]"), L),
                   `[1, 2, 3 ] a`, Tail).
L = [1, 2, 3],
Tail = [32, 97].
```

\[det\]**optional**(`:Match, :Default`)`//`  
Perform an optional match, executing `Default` if `Match` is not matched. This is comparable to [ignore/1](metacall.html#ignore/1). Both `Match` and `Default` are DCG body terms. `Default` is typically used to instantiate the output variables of `Match`, but may also be used to match a default representation. Using `[]` for `Default` succeeds without any additional actions if `Match` fails. For example:

``` code
?- phrase(optional(number(X), {X=0}), `23`, Tail).
X = 23,
Tail = [].
?- phrase(optional(number(X), {X=0}), `aap`, Tail).
X = 0,
Tail = `aap`.
```

\[det\]**foreach**(`:Generator, :Element`)`//`  
\[det\]**foreach**(`:Generator, :Element, :Sep`)`//`  
Generate a list from the solutions of `Generator`. This predicate collects all solutions of `Generator`, applies `Element` for each solution and `Sep` *between* each pair of solutions. For example:

``` code
?- phrase(foreach(between(1,5,X), number(X), ", "), L).
L = "1, 2, 3, 4, 5".
```
