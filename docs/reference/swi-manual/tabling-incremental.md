
## 7.7 Incremental tabling

Incremental tabling maintains the consistency of a set of tabled predicates that depend on a set of dynamic predicates. Both the tabled and dynamic predicates must have the property `incremental` set. See [dynamic/1](dynamic.html#dynamic/1) and [table/1](tabling-preds.html#table/1).

Incremental tabling causes the engine to connect the *answer tries* and incremental dynamic predicates in an *Incremental Dependency Graph* (IDG). Modifications ([asserta/1](db.html#asserta/1), [retract/1](db.html#retract/1), [retractall/1](db.html#retractall/1) and friends) of an incremental dynamic predicate mark all depending tables as invalid. Subsequent usage of these tables forces re-evaluation.

Re-evaluation of invalidated tables happens on demand, i.e., on access to an invalid table. First the dependency graph of invalid tables that lead to dynamic predicates is established. Next, tables are re-evaluated in *bottom-up* order. For each re-evaluated table the system determines whether the new table has changed. If the table has not changed, this event is propagated to the *affected* nodes of the IDG and no further re-evaluation may be needed. Consider the following program:

``` code
:- table (p/1, q/1) as incremental.
:- dynamic([d/1], [incremental(true)]).

p(X) :- q(X).
q(X) :- d(X), X < 10.

d(1).
```

Executing this program creates tables for `X=1` for p/1 and q/1 . After calling `assert(d(100))` the tables for p/1 and q/1 have an *invalid count* of `1`. Re-running `p(X)` first re-evaluates q/1 (bottom-up) which results to the same table as `X=100` does not lead to a new answer. Re-evaluation clears the invalid count for q/1 and, because the q/1 tables is not changed, decrements the invalid count of affected tables. This sets the *invalid count* for p/1 to zero, completing the re-evaluation.

Note that invalidating and re-evaluation is done at the level of tables. Notably asserting a clause invalidates all affected tables and may lead to re-evaluating of all these tables. Incremental tabling automates manual abolishing of invalid tables in a changing world and avoids unnecessary re-evaluation if indirectly affected tables prove unaffected because the answer set of dependent tables is unaffected by the change. This is the same policy as implemented in XSB [Swift, 2014](Bibliography.html#DBLP:journals/tplp/Swift14). Future versions may implement a more fine grained approach.
