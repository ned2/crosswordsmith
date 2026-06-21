
## A.59 library(tableutil): Table inspection and statistics utilities

This library provides tools intended to be used at the interactive toplevel for inspecting tables, table dependencies and table statistics.

\[nondet\]**table_statistics**(`?Stat, -Value`)  
\[nondet\]**table_statistics**(`?Variant, ?Stat, -Value`)  
Give summary statistics for the tables associated with all subgoals of `Variant`. The [table_statistics/2](tableutil.html#table_statistics/2) version considers all tables.

The values for `Stat` are:

**tables**  
Total number of answer tries

**answers**  
Total number of answers in the combined tries

**duplicate_ratio**  
Ratio of generated (and thus ignored) duplicate answers. `1` means no duplicates. `2` means for every answer there was (on everage) a duplicate generated.

**space_ratio**  
Number of nodes with a value divided by the total number of nodes in a trie. The maximum is `1`. A low number implies that a large amount of differently shaped data is included in the answer tries.

**complete_call**  
Number of times answers are generated from a completed table, i.e., times answers are *reused*.

**invalidated**  
Number of times an incremental table was invalidated.

**reevaluated**  
Number of times an invalidated table wa reevaluated. If lower than `invalidated` this implies that dependent nodes of the IDG were reevaluated to the same answer set.

**space**  
Summed memory usage of the answer tries in bytes.

**compiled_space**  
Summed size for the compiled representation of completed tables.

**table_statistics**  
Print a summary of statistics relevant to tabling.

See also  
[table_statistics/2](tableutil.html#table_statistics/2) for an explanation

**table_statistics**(`:Variant`)  
Print a summary for the statistics of all tables for subgoals of `Variant`. See [table_statistics/2](tableutil.html#table_statistics/2) for an explanation.

\[det\]**table_statistics_by_predicate**  
\[det\]**table_statistics_by_predicate**(`+Options`)  
Print statistics on memory usage and lookups per predicate. The version without options dumps all predicates without ordering. `Options`:

**order_by**(`+Key`)  
Order the predicates according to `Key`. Default is `tables`, the number of answer tables. See [table_statistics/2](tableutil.html#table_statistics/2) for a list of values for `Key`.

**top**(`N`)  
Only show the top `N` predicates.

**module**(`Module`)  
Limit the results to predicates of the given module.

**tstat**(`?Value, ?Top`)  
**tstat**(`?Variant, ?Value, ?Top`)  
Print the top-N (for positive `Top`) or bottom-N (for negative `Top`) for `Stat` for all tabled subgoals of `Variant` (or all tabled subgoals for [tstat/2](tableutil.html#tstat/2)). Stat is one of

**answers**  
**duplicate_ratio**  
**space_ratio**  
**gen**(`call`)  
**space**  
**compiled_space**  
See [table_statistics/2](tableutil.html#table_statistics/2).

**deadlock**  
Times this table was involved in a deadlock (cycle of threads waiting for each others table to complete)

**wait**  
Times a thread waited for this table.

**variables**  
The number of variables in the variant. The tabling logic adds a term `ret(...)` to the table for each answer, where each variable is an argument of the `ret(...)` term. The arguments are placed in depth-first lef-right order they appear in the variant. Optimal behaviour of the trie is achieved if the variance is as much as possible to the rightmost arguments. Poor allocation shows up as a low `space_ratio` statistics.

Below are some examples

\[det\]**tdump**  
\[det\]**tdump**(`:Goal`)  
\[det\]**tdump**(`:Goal, +Options`)  
Dump all tables and their status that *unify* with `Goal`. `Options`:

**scope**(`Scope`)  
Limit displayed tables to `local` or `global`.

**limit**(`Count`)  
Limit the number of answers displayed to `Count`

**reset**(`Boolean`)  
If `true`, also show reset (fresh) global tables. These are tables that have been abolished.

\[det\]**tidg**  
\[det\]**tidg**(`:Goal`)  
Dump the incremental dependency graph. idg/1 dumps the graph around a given node

\[det\]**summarize_idg**  
\[det\]**summarize_idg**(`+TopN`)  
Implements XSB's `statistics(summarize_idg)`
