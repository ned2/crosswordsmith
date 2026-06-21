
# 7 Tabled execution (SLG resolution)

This chapter describes SWI-Prolog's support for *Tabled execution* for one or more Prolog predicates, also called *SLG resolution*. Tabling a predicate provides two properties:

1.  Re-evaluation of a tabled predicate is avoided by *memoizing* the answers. This can realise huge performance enhancements as illustrated in [section 7.1](tabling-memoize.html#sec:7.1). It also comes with two downsides: the memoized answers are not automatically updated or invalidated if the world (set of predicates on which the answers depend) changes and the answer tables must be stored (in memory).
2.  *Left recursion*, a goal calling a *variant* of itself recursively and thus *looping* under the normal Prolog SLD resolution is avoided by suspending the variant call and resuming it with answers from the table. This is illustrated in [section 7.2](tabling-non-termination.html#sec:7.2).

Tabling is particularly suited to simplify inference over a highly entangled set of predicates that express axioms and rules in a static (not changing) world. When using SLD resolution for such problems, it is hard to ensure termination and avoid frequent recomputation of intermediate results. A solution is to use Datalog style bottom-up evaluation, i.e., applying rules on the axioms and derived facts until a fixed point is reached. However, bottom-up evaluation typically derives many facts that are never used. Tabling provides a *goal oriented* resolution strategy for such problems and is enabled simply by adding a [table/1](tabling-preds.html#table/1) directive to the program.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[7.1 Example 1: using tabling for memoizing](tabling-memoize.html)

[7.2 Example 2: avoiding non-termination](tabling-non-termination.html)

[7.3 Answer subsumption or mode directed tabling](tabling-mode-directed.html)

[7.4 Tabling for impure programs](tnotpure.html)

[7.5 Variant and subsumptive tabling](tabling-subsumptive.html)

[7.6 Well Founded Semantics](WFS.html)

[7.6.1 Well founded semantics and the toplevel](WFS.html#sec:7.6.1)

[7.7 Incremental tabling](tabling-incremental.html)

[7.8 Monotonic tabling](tabling-monotonic.html)

[7.8.1 Eager and lazy monotonic tabling](tabling-monotonic.html#sec:7.8.1)

[7.8.2 Tracking new answers to monotonic tables](tabling-monotonic.html#sec:7.8.2)

[7.8.3 Monotonic tabling with external data](tabling-monotonic.html#sec:7.8.3)

[7.9 Shared tabling](tabling-shared.html)

[7.9.1 Abolishing shared tables](tabling-shared.html#sec:7.9.1)

[7.9.2 Status and future of shared tabling](tabling-shared.html#sec:7.9.2)

[7.10 Tabling and constraints](tabling-constraints.html)

[7.11 Tabling restraints: bounded rationality and tripwires](tabling-restraints.html)

[7.11.1 Restraint subgoal size](tabling-restraints.html#sec:7.11.1)

[7.11.2 Restraint answer size](tabling-restraints.html#sec:7.11.2)

[7.11.3 Restraint answer count](tabling-restraints.html#sec:7.11.3)

[7.12 Tabling predicate reference](tabling-preds.html)

[7.13 About the tabling implementation](tabling-about.html)

[7.13.1 Status of tabling](tabling-about.html#sec:7.13.1)
