
## 7.12 Tabling predicate reference

:- **table**(`:Specification`)  
Prepare the predicates specified by `Specification` for tabled execution. `Specification` is a *comma-list*, each member specifying tabled execution for a specific predicate. The individual specification is either a *predicate indicator* (name/arity or name//arity) or head specifying tabling with *answer subsumption*.

Although [table/1](tabling-preds.html#table/1) is normally used as a directive, SWI-Prolog allows calling it as a runtime predicate to prepare an existing predicate for tabled execution. The predicate [untable/1](tabling-preds.html#untable/1) can be used to remove the tabling instrumentation from a predicate.

The example below prepares the predicate edge/2 and the non-terminal statement//1 for tabled execution.

``` code
:- table edge/2, statement//1.
```

Below is an example declaring a predicate to use tabling with *answer subsumption*. Answer subsumption or *mode directed tabling* is discussed in [section 7.3](tabling-mode-directed.html#sec:7.3).

``` code
:- table connection(_,_,min).
```

Additional tabling options can be provided using a term **as/2**, which can be applied to a single specification or a comma list of specifications. The options themselves are a comma-list of one or more of the following atoms:

**variant**  
Default. Create a table for each call variant.

**subsumptive**  
Instead of creating a new table for each call variant, check whether there is a completed table for a more general goal and if this is the case extract the answers from this table. See [section 7.5](tabling-subsumptive.html#sec:7.5).

**shared**  
Declare that the table shall be shared between threads. See [section 7.9](tabling-shared.html#sec:7.9)

**private**  
Declare that the table shall be local to the calling thread. See [section 7.9](tabling-shared.html#sec:7.9)

**incremental**  
Declare that the table depends on other tables and *incremental* dynamic predicates. See [section 7.7](tabling-incremental.html#sec:7.7).

**dynamic**  
Declare that the predicate is dynamic. Often used together with `incremental`.

This syntax is closely related to the table declarations used in XSB Prolog. Where in XSB `as` is an operator with priority above the priority of the comma, it is an operator with priority below the comma in SWI-Prolog. Therefore, multiple predicates or options must be enclosed in parenthesis. For example:

``` code
:- table p/1 as subsumptive.
:- table (q/1, r/2) as subsumptive.
```

**tnot**(`:Goal`)  
The [tnot/1](tabling-preds.html#tnot/1) predicate implements *tabled negation*. This predicate realises *Well Founded Semantics*. See [section 7.6](WFS.html#sec:7.6) for details.

**not_exists**(`:Goal`)  
Handles tabled negation for non-ground (*floundering*) `Goal` as well as non tabled goals. If `Goal` is ground and tabled [not_exists/1](tabling-preds.html#not_exists/1) calls [tnot/1](tabling-preds.html#tnot/1). Otherwise it used `tabled_call(Goal)` to create a table and subsequently uses [tnot/1](tabling-preds.html#tnot/1) on the created table.

Logically, `not_exists(p(X))` is defined as tnot(`&exist``X`(p(`X`)))

Note that each `Goal` variant populates a table for [tabled_call/1](tabling-preds.html#tabled_call/1). Applications may need to abolish such tables to limit memory usage or guarantee consistency‘after the world changed’.

**tabled_call**(`:Goal`)  
Helper predicate for [not_exists/1](tabling-preds.html#not_exists/1). Defined as below. The helper is public because application may need to abolish its tables.

``` code
:- table tabled_call/1 as variant.
tabled_call(Goal) :- call(Goal).
```

**current_table**(`:Variant, -Trie`)  
True when `Trie` is the answer table for `Variant`.

**untable**(`:Specification`)  
Remove the tables and tabling instrumentation for the specified predicates. `Specification` is compatible with [table/1](tabling-preds.html#table/1), although tabling with *answer subsumption* may be removed using a name/arity specification. The [untable/1](tabling-preds.html#untable/1) predicate is first of all intended for examining the effect of various tabling scenarios on a particular program interactively from the toplevel.

Note that although using [untable/1](tabling-preds.html#untable/1) followed by [table/1](tabling-preds.html#table/1) may be used to flush all tables associated with the given predicate(s), flushing tables should be done using one of the table abolish predicates both for better performance and compatibility with other Prolog implementations: [abolish_all_tables/0](tabling-preds.html#abolish_all_tables/0), [abolish_private_tables/0](tabling-preds.html#abolish_private_tables/0), [abolish_shared_tables/0](tabling-preds.html#abolish_shared_tables/0), [abolish_module_tables/1](tabling-preds.html#abolish_module_tables/1) or [abolish_table_subgoals/1](tabling-preds.html#abolish_table_subgoals/1). For example, to remove all tables for p/3 , run the goal below. The predicate [functor/3](manipterm.html#functor/3) may be used to create a *head term* from a given name and arity.

``` code
?- abolish_table_subgoals(p(_,_,_)).
```

**abolish_all_tables**  
Remove all tables, both *private* and *shared* (see [section 7.9](tabling-shared.html#sec:7.9)). Since the introduction of *incremental tabling* (see [section 7.7](tabling-incremental.html#sec:7.7)) abolishing tables is rarely required to maintain consistency of the tables with a changed environment. Tables may be abolished regardless of the current state of the table. *Incomplete* tables are flagged for destruction when they are completed. See [section 7.9.1](tabling-shared.html#sec:7.9.1) for the semantics of destroying shared tables and the following predicates for destroying a subset of the tables: [abolish_private_tables/0](tabling-preds.html#abolish_private_tables/0), [abolish_shared_tables/0](tabling-preds.html#abolish_shared_tables/0), [abolish_table_subgoals/1](tabling-preds.html#abolish_table_subgoals/1) and [abolish_module_tables/1](tabling-preds.html#abolish_module_tables/1).

**abolish_private_tables**  
Abolish all tables that are private to this thread.

**abolish_shared_tables**  
Abolish all tables that are shared between threads. See also [section 7.9.1](tabling-shared.html#sec:7.9.1)

**abolish_table_subgoals**(`:Subgoal`)  
Abolish all tables that unify with `SubGoal`. Tables that have undefined answers that depend of the abolished table are abolished as well (recursively). For example, given the program below, `abolish_table_subgoals(und)` will also abolish the table for p/0 because its answer refers to und/0 .

``` code
p :- und.
und :- tnot(und).
```

**abolish_module_tables**(`+Module`)  
Remove all tables that belong to predicates in `Module`.

**abolish_nonincremental_tables**  
**abolish_nonincremental_tables**(`+Options`)  
Similar to [abolish_all_tables/0](tabling-preds.html#abolish_all_tables/0), but does not abolish *incremental* tables as their consistency is maintained by the system. Options:

**on_incomplete**(`Action`)  
`Action` is one of `skip` or `error`. If `Action` is `skip`, do not delete the table.^(bugXSB marks such tables for deletion after completion. That is not yet implemented.)
