
## 7.3 Answer subsumption or mode directed tabling

Tabling as defined above has a serious limitation. Although the definition of connection/2 from section [section 7.2](tabling-non-termination.html#sec:7.2) can compute the transitive closure of connected cities, it cannot provide you with a route to travel. The reason is that there are infinitely many routes if there are cycles in the network and each new route found will be added to the answer table and cause the tabled execution's completion algorithm to search for more routes, eventually running out of memory.

The solution to this problem is called *mode directed tabling* or *answer subsumption*.^(187The term *answer subsumption* is used by XSB and *mode directed tabling* by YAP and B-Prolog. The idea is that some arguments are considered‘outputs’, where multiple values for the same‘input’are combined. Possibly *answer aggregation* would have been a better name.) In this execution model one or more arguments are *not* added to the table. Instead, we remember a single *aggregated* value for these arguments. The example below is derived from [section 7.2](tabling-non-termination.html#sec:7.2) and returns the connection as a list of cities. This argument is defined as a *moded* argument using the `lattice(PI)` mode.^(188This mode is compatible to XSB Prolog.) This causes the tabling engine each time that it finds an new path to call shortest/3 and keep the shortest route.

``` code
:- table
    connection(_,_,lattice(shortest/3)).

shortest(P1, P2, P):-
    length(P1, L1),
    length(P2, L2),
    (   L1 < L2
    ->  P = P1
    ;   P = P2
    ).

connection(X, Y, [X,Y]) :-
    connection(X, Y).
connection(X, Y, P) :-
    connection(X, Z, P0),
    connection(Z, Y),
    append(P0, [Y], P).
```

The mode declaration scheme is equivalent to XSB with partial compatibility support for YAP and B-Prolog. The `lattice(PI)` mode is the most general mode. The YAP `all` (B-Prolog `@`) mode is not yet supported. The list below describes the supported modes and indicates the portability.

**Var**  
**`+`**  
**index**  
A variable (XSB), the atom `index` (YAP) or a `+` (B-Prolog, YAP) declare that the argument is tabled normally.

**lattice**(`Pred`)  
`Pred` denotes a predicate with arity 3. It may be specified as an predicate indicator (`Name`/3), plain predicate name (`Name`) or a head term `Name(_,_,_)`. On each answer, `PI` is called with three arguments: the current aggregated answer and the new answer are inputs. The last argument must be unified with a term that represents the new aggregated answer.

**po**(`PI`)  
*Partial Ordering*. The new answer is added iff `call(PI, +Old, +Answer)` succeeds. For example, `po('<'/2)` accumulates the smallest result. In SWI-Prolog the arity (2) may be omitted, resulting in `po(<)`.

**`-`**  
**first**  
The atom `-` (B-Prolog, YAP) and `first` (YAP) declare to keep the first answer for this argument.

**last**  
The atom `last` (YAP) declares to keep the last answer.

**min**  
The atom `min` (YAP) declares to keep the smallest answer according to the standard order of terms (see [@\</2](compare.html#@%3C/2)). Note that in SWI-Prolog the standard order of terms orders numbers by value.

**max**  
The atom `max` (YAP) declares to keep the largest answer according to the standard order of terms (see [@\>/2](compare.html#@%3E/2)). Note that in SWI-Prolog the standard order of terms orders numbers by value.

**sum**  
The atom `sum` (YAP) declares to sum numeric answers.
