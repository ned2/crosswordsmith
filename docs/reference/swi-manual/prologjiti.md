
## A.42 library(prolog_jiti): Just In Time Indexing (JITI) utilities

To be done  
Use [print_message/2](printmsg.html#print_message/2) and dynamically figure out the column width.

This module provides utilities to examine just-in-time indexes created by the system and can help diagnosing space and performance issues.

\[det\]**jiti_list**  
\[det\]**jiti_list**(`:Spec`)  
List the JITI (Just In Time Indexes) of selected predicates. The predicate [jiti_list/0](prologjiti.html#jiti_list/0) list all just-in-time indexed predicates. The predicate [jiti_list/1](prologjiti.html#jiti_list/1) takes one of the patterns below. All parts except for Name can be variables. The last pattern takes an arbitrary number of arguments.

- Module:Head
- Module:Name/Arity
- Module:Name

The columns use the following notation:

- The *Indexed* column describes the `argument(s)` indexed:
  - A plain integer refers to a 1-based argument number
  - `A+B` is a multi-argument index on the arguments `A` and `B`.
  - `P:L` is a deep-index `L` on sub-argument `P`. For example, `1/2:2+3` is an index of the 2nd and 3rd argument of the 2nd argument of a compound on the first argument of the predicate. This implies `x` and `y` in the head `p(f(_,g(_,x,y)))`
- The `Buckets` specifies the number of buckets of the hash table
- The `Speedup` specifies the selectivity of the index
- The `Flags` describes additional properties, currently:
  - `L` denotes that the index contains multiple compound terms with the same name/arity that may be used to create deep indexes. The deep indexes themselves are created as just-in-time indexes.
  - `V` denotes the index is *virtual*, i.e., it has not yet been materialized.

\[det\]**jiti_suggest_modes**  
\[det\]**jiti_suggest_modes**(`:Spec`)  
Propose modes for the predicates referenced by `Spec`. This utility may be executed *after* a clean load of your program and after running the program. It searches for static predicates that have been called and (thus) have been examined for candidate indexes. If candidate indexes have not been materialized this implies that the predicate was never called with a nonvar value for the corresponding argument. Adding a [mode/1](dynamic.html#mode/1) declaration may be used to inform the system thereof. The system will never examine arguments for indexing that have been declared as mode `-`.

**Note:** This predicate merely detects that some predicate is never called with instantiated specific arguments **during this run**. The user should verify whether the suggested `-` arguments are correct and typically complete the mode by changing `?` into `+` (or `-`) where applicable. Currently, in SWI-Prolog, [mode/1](dynamic.html#mode/1) declarations have no effect on the semantics of the code. In particular, a predicate that declares some argument as `-` may be called with this argument instantiated. This may change in the future.

|  |  |
|----|----|
| `Spec` | uses the same conventions as [jiti_list/1](prologjiti.html#jiti_list/1). |
