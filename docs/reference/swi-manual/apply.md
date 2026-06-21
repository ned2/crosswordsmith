
## A.3 library(apply): Apply predicates on a list

See also  
\- `apply_macros.pl` provides compile-time expansion for part of this library.  
- [http://www.cs.otago.ac.nz/staffpriv/ok/pllib.htm](http://www.cs.otago.ac.nz/staffpriv/ok/pllib.htm)  
- Unit test code in `tests/library/test_apply.pl`

To be done  
Add include/4, include/5, exclude/4, exclude/5

This module defines meta-predicates that apply a predicate on all members of a list.

All predicates support partial application in the Goal argument. This means that these calls are identical:

``` code
?- maplist(=, [foo, foo], [X, Y]).
?- maplist(=(foo), [X, Y]).
```

\[det\]**include**(`:Goal, +List1, ?List2`)  
Filter elements for which `Goal` succeeds. True if `List2` contains those elements Xi of `List1` for which `call(Goal, Xi)` succeeds.

See also  
[exclude/3](apply.html#exclude/3), [partition/4](apply.html#partition/4), [convlist/3](apply.html#convlist/3).

Compatibility  
Older versions of SWI-Prolog had sublist/3 with the same arguments and semantics.

\[det\]**exclude**(`:Goal, +List1, ?List2`)  
Filter elements for which `Goal` fails. True if `List2` contains those elements Xi of `List1` for which `call(Goal, Xi)` fails.

See also  
[include/3](apply.html#include/3), [partition/4](apply.html#partition/4)

\[det\]**partition**(`:Pred, +List, ?Included, ?Excluded`)  
Filter elements of `List` according to `Pred`. True if `Included` contains all elements for which `call(Pred, X)` succeeds and `Excluded` contains the remaining elements.

See also  
[include/3](apply.html#include/3), [exclude/3](apply.html#exclude/3), [partition/5](apply.html#partition/5).

\[semidet\]**partition**(`:Pred, +List, ?Less, ?Equal, ?Greater`)  
Filter `List` according to `Pred` in three sets. For each element Xi of `List`, its destination is determined by `call(Pred, Xi, Place)`, where Place must be unified to one of `<`, `=` or `>`. `Pred` must be deterministic.

See also  
[partition/4](apply.html#partition/4)

**maplist**(`:Goal, ?List1`)  
**maplist**(`:Goal, ?List1, ?List2`)  
**maplist**(`:Goal, ?List1, ?List2, ?List3`)  
**maplist**(`:Goal, ?List1, ?List2, ?List3, ?List4`)  
True if `Goal` is successfully applied on all matching elements of the list. The maplist family of predicates is defined as:

``` code
maplist(G, [X_11, ..., X_1n],
           [X_21, ..., X_2n],
           ...,
           [X_m1, ..., X_mn]) :-
   call(G, X_11, ..., X_m1),
   call(G, X_12, ..., X_m2),
   ...
   call(G, X_1n, ..., X_mn).
```

This family of predicates is deterministic iff `Goal` is deterministic and `List1` is a proper list, i.e., a list that ends in `[]`.

\[det\]**convlist**(`:Goal, +ListIn, -ListOut`)  
Similar to [maplist/3](apply.html#maplist/3), but elements for which `call(Goal, ElemIn, _)` fails are omitted from `ListOut`. For example (using `library(yall)`):

``` code
?- convlist([X,Y]>>(integer(X), Y is X^2),
            [3, 5, foo, 2], L).
L = [9, 25, 4].
```

Compatibility  
Also appears in YAP `library(maplist)` and SICStus `library(lists)`.

**foldl**(`:Goal, +List, +V0, -V`)  
**foldl**(`:Goal, +List1, +List2, +V0, -V`)  
**foldl**(`:Goal, +List1, +List2, +List3, +V0, -V`)  
**foldl**(`:Goal, +List1, +List2, +List3, +List4, +V0, -V`)  
Fold an ensemble of *m* (0 `<=` *m* `<=` 4) lists of length *n* head-to-tail ("fold-left"), using columns of *m* list elements as arguments for `Goal`. The `foldl` family of predicates is defined as follows, with `V0` an initial value and `V` the final value of the folding operation:

``` code
foldl(G, [X_11, ..., X_1n],
         [X_21, ..., X_2n],
         ...,
         [X_m1, ..., X_mn], V0, V) :-
   call(G, X_11, ..., X_m1, V0, V1),
   call(G, X_12, ..., X_m2, V1, V2),
   ...
   call(G, X_1n, ..., X_mn, V<n-1>, V).
```

No implementation for a corresponding `foldr` is given. A `foldr` implementation would consist in first calling [reverse/2](lists.html#reverse/2) on each of the *m* input lists, then applying the appropriate `foldl`. This is actually more efficient than using a properly programmed-out recursive algorithm that cannot be tail-call optimized.

**scanl**(`:Goal, +List, +V0, -Values`)  
**scanl**(`:Goal, +List1, +List2, +V0, -Values`)  
**scanl**(`:Goal, +List1, +List2, +List3, +V0, -Values`)  
**scanl**(`:Goal, +List1, +List2, +List3, +List4, +V0, -Values`)  
Scan an ensemble of *m* (0 `<=` *m* `<=` 4) lists of length *n* head-to-tail ("scan-left"), using columns of *m* list elements as arguments for `Goal`. The `scanl` family of predicates is defined as follows, with `V0` an initial value and `V` the final value of the scanning operation:

``` code
scanl(G, [X_11, ..., X_1n],
         [X_21, ..., X_2n],
         ...,
         [X_m1, ..., X_mn], V0, [V0, V1, ..., Vn] ) :-
   call(G, X_11, ..., X_m1, V0, V1),
   call(G, X_12, ..., X_m2, V1, V2),
   ...
   call(G, X_1n, ..., X_mn, V<n-1>, Vn).
```

`scanl` behaves like a `foldl` that collects the sequence of values taken on by the `Vx` accumulator into a list.
