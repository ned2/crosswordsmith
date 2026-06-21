
## A.21 library(increval): Incremental dynamic predicate modification

Compatibility  
XSB

This module emulates the XSB module `increval`. This module serves two goals: (1) provide alternatives for the dynamic clause manipulation predicates that propagate into the incremental tables and (2) query the dynamically maintained *Incremental Depency Graph* (IDG).

The change propagation for incremental dynamic predicates. SWI-Prolog relies in [prolog_listen/2](prolog-event.html#prolog_listen/2) to forward any change to dynamic predicates to the table IDG and incr_assert/1 and friends thus simply call the corresponding database update.

\[nondet\]**is_incremental_subgoal**(`?SubGoal`)  
This predicate non-deterministically unifies Subgoal with incrementally tabled subgoals that are currently table entries.

\[nondet\]**incr_directly_depends**(`:Goal1, :Goal2`)  
True if `Goal1` depends on `Goal2` in the IDG.

Compatibility  
: In XSB, at least one of Goal 1 or Goal 2 must be bound. This implementation may be used with both arguments unbound.

\[nondet\]**incr_trans_depends**(`:Goal1, Goal2`)  
True for each pair in the transitive closure of `incr_directly_depends(G1, G2)`.

\[det\]**incr_invalid_subgoals**(`-List`)  
`List` is a sorted list (set) of the incremental subgoals that are currently invalid.

\[semidet\]**incr_is_invalid**(`:Subgoal`)  
True when `Subgoal`’s table is marked as invalid.

\[det\]**incr_invalidate_calls**(`:Goal`)  
Invalidate all tables for subgoals of `Goal` as well as tables that are affected by these.

\[det\]**incr_invalidate_call**(`:Goal`)  
This is the XSB name, but the manual says [incr_invalidate_calls/1](increval.html#incr_invalidate_calls/1) and the comment with the code suggests this is misnamed.

deprecated  
Use [incr_invalidate_calls/1](increval.html#incr_invalidate_calls/1).

**incr_table_update**  
Updated all invalid tables

\[det\]**incr_propagate_calls**(`:Answer`)  
Activate the monotonic answer propagation similarly to when a new fact is asserted for a monotonic dynamic predicate. The `Answer` term must match a monotonic dynamic predicate.
