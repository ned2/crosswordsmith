
## 2.9 Reuse of top-level bindings

Bindings resulting from the successful execution of a top-level goal are asserted in a database *if they are not too large* (as defined by the Prolog flag [toplevel_var_size](flags.html#flag:toplevel_var_size)). These values may be reused in further top-level queries as \$Var. If the same variable name is used in a subsequent query the system associates the variable with the latest binding. Example:

``` code
1 ?- maplist(plus(1), `hello`, X).
X = [105,102,109,109,112].

2 ?- format('~s~n', [$X]).
ifmmp
true.

3 ?-
```

**Figure 1 :** Reusing top-level bindings

Note that variables may be set by executing [=/2](compare.html#=/2):

``` code
6 ?- X = statistics.
X = statistics.

7 ?- $X.
% Started at Fri Aug 24 16:42:53 2018
% 0.118 seconds cpu time for 456,902 inferences
% 7,574 atoms, 4,058 functors, 2,912 predicates, 56 modules, 109,791 VM-codes
%
%                     Limit   Allocated      In use
% Local  stack:           -       20 Kb    1,888  b
% Global stack:           -       60 Kb       36 Kb
% Trail  stack:           -       30 Kb    4,112  b
%        Total:    1,024 Mb      110 Kb       42 Kb
%
% 3 garbage collections gained 178,400 bytes in 0.000 seconds.
% 2 clause garbage collections gained 134 clauses in 0.000 seconds.
% Stack shifts: 2 local, 2 global, 2 trail in 0.000 seconds
% 2 threads, 0 finished threads used 0.000 seconds
true.
```
