
# A The SWI-Prolog library

This chapter documents the SWI-Prolog library. As SWI-Prolog provides auto-loading, there is little difference between library predicates and built-in predicates. Part of the library is therefore documented in the rest of the manual. Library predicates differ from built-in predicates in the following ways:

- User definition of a built-in leads to a permission error, while using the name of a library predicate is allowed.
- If autoloading is disabled explicitly or because trapping unknown predicates is disabled (see [unknown/2](debugger.html#unknown/2) and [current_prolog_flag/2](flags.html#current_prolog_flag/2)), library predicates must be loaded explicitly.
- Using libraries reduces the footprint of applications that don't need them.

> *The documentation of the library has just started. Material from the standard packages should be moved here, some material from other parts of the manual should be moved too and various libraries are not documented at all.*

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[A.1 library(aggregate): Aggregation operators on backtrackable predicates](aggregate.html)

[A.2 library(ansi_term): Print decorated text to ANSI consoles](ansiterm.html)

[A.3 library(apply): Apply predicates on a list](apply.html)

[A.4 library(assoc): Association lists](assoc.html)

[A.4.1 Introduction](assoc.html#sec:A.4.1)

[A.4.2 Creating association lists](assoc.html#sec:A.4.2)

[A.4.3 Querying association lists](assoc.html#sec:A.4.3)

[A.4.4 Modifying association lists](assoc.html#sec:A.4.4)

[A.4.5 Conversion predicates](assoc.html#sec:A.4.5)

[A.4.6 Reasoning about association lists and their elements](assoc.html#sec:A.4.6)

[A.5 library(broadcast): Broadcast and receive event notifications](broadcast.html)

[A.6 library(charsio): I/O on Lists of Character Codes](charsio.html)

[A.7 library(check): Consistency checking](check.html)

[A.8 library(clpb): CLP(B): Constraint Logic Programming over Boolean Variables](clpb.html)

[A.8.1 Introduction](clpb.html#sec:A.8.1)

[A.8.2 Boolean expressions](clpb.html#sec:A.8.2)

[A.8.3 Interface predicates](clpb.html#sec:A.8.3)

[A.8.4 Examples](clpb.html#sec:A.8.4)

[A.8.5 Obtaining BDDs](clpb.html#sec:A.8.5)

[A.8.6 Enabling monotonic CLP(B)](clpb.html#sec:A.8.6)

[A.8.7 Example: Pigeons](clpb.html#sec:A.8.7)

[A.8.8 Example: Boolean circuit](clpb.html#sec:A.8.8)

[A.8.9 Acknowledgments](clpb.html#sec:A.8.9)

[A.8.10 CLP(B) predicate index](clpb.html#sec:A.8.10)

[A.9 library(clpfd): CLP(FD): Constraint Logic Programming over Finite Domains](clpfd.html)

[A.9.1 Introduction](clpfd.html#sec:A.9.1)

[A.9.2 Arithmetic constraints](clpfd.html#sec:A.9.2)

[A.9.3 Declarative integer arithmetic](clpfd.html#sec:A.9.3)

[A.9.4 Example: Factorial relation](clpfd.html#sec:A.9.4)

[A.9.5 Combinatorial constraints](clpfd.html#sec:A.9.5)

[A.9.6 Domains](clpfd.html#sec:A.9.6)

[A.9.7 Example: Sudoku](clpfd.html#sec:A.9.7)

[A.9.8 Residual goals](clpfd.html#sec:A.9.8)

[A.9.9 Core relations and search](clpfd.html#sec:A.9.9)

[A.9.10 Example: Eight queens puzzle](clpfd.html#sec:A.9.10)

[A.9.11 Optimisation](clpfd.html#sec:A.9.11)

[A.9.12 Reification](clpfd.html#sec:A.9.12)

[A.9.13 Enabling monotonic CLP(FD)](clpfd.html#sec:A.9.13)

[A.9.14 Custom constraints](clpfd.html#sec:A.9.14)

[A.9.15 Applications](clpfd.html#sec:A.9.15)

[A.9.16 Acknowledgments](clpfd.html#sec:A.9.16)

[A.9.17 CLP(FD) predicate index](clpfd.html#sec:A.9.17)

[A.9.17.1 Arithmetic constraints](clpfd.html#sec:A.9.17.1)

[A.9.17.2 Membership constraints](clpfd.html#sec:A.9.17.2)

[A.9.17.3 Enumeration predicates](clpfd.html#sec:A.9.17.3)

[A.9.17.4 Global constraints](clpfd.html#sec:A.9.17.4)

[A.9.17.5 Reification predicates](clpfd.html#sec:A.9.17.5)

[A.9.17.6 Reflection predicates](clpfd.html#sec:A.9.17.6)

[A.9.17.7 FD set predicates](clpfd.html#sec:A.9.17.7)

[A.9.17.8 FD miscellaneous predicates](clpfd.html#sec:A.9.17.8)

[A.9.18 Closing and opening words about CLP(FD)](clpfd.html#sec:A.9.18)

[A.10 library(clpqr): Constraint Logic Programming over Rationals and Reals](clpqr.html)

[A.10.1 Solver predicates](clpqr.html#sec:A.10.1)

[A.10.2 Syntax of the predicate arguments](clpqr.html#sec:A.10.2)

[A.10.3 Use of unification](clpqr.html#sec:A.10.3)

[A.10.4 Non-linear constraints](clpqr.html#sec:A.10.4)

[A.10.5 Status and known problems](clpqr.html#sec:A.10.5)

[A.11 library(csv): Process CSV (Comma-Separated Values) data](csv.html)

[A.12 library(dcg/basics): Various general DCG utilities](basics.html)

[A.13 library(dcg/high_order): High order grammar operations](highorder.html)

[A.14 library(debug): Print debug messages and test assertions](debug.html)

[A.15 library(dicts): Dict utilities](dicts.html)

[A.16 library(error): Error generating support](error.html)

[A.17 library(exceptions): Exception classification](exceptions.html)

[A.18 library(fastrw): Fast reading and writing of terms](fastrw.html)

[A.19 library(gensym): Generate unique symbols](gensym.html)

[A.20 library(heaps): heaps/priority queues](heaps.html)

[A.21 library(increval): Incremental dynamic predicate modification](increval.html)

[A.22 library(intercept): Intercept and signal interface](intercept.html)

[A.23 library(iostream): Utilities to deal with streams](iostream.html)

[A.24 library(listing): List programs and pretty print clauses](listing.html)

[A.25 library(lists): List Manipulation](lists.html)

[A.26 library(macros): Macro expansion](macros.html)

[A.26.1 Defining and using macros](macros.html#sec:A.26.1)

[A.26.2 Implementation details](macros.html#sec:A.26.2)

[A.26.3 Predicates](macros.html#sec:A.26.3)

[A.27 library(main): Provide entry point for scripts](main.html)

[A.28 library(nb_set): Non-backtrackable set](nb_set.html)

[A.29 library(writef): Old-style formatted write](writef.html)

[A.30 library(www_browser): Open a URL in the users browser](wwwbrowser.html)

[A.31 library(occurs): Finding and counting sub-terms](occurs.html)

[A.32 library(option): Option list processing](option.html)

[A.33 library(optparse): command line parsing](optparse.html)

[A.33.1 Notes and tips](optparse.html#sec:A.33.1)

[A.34 library(ordsets): Ordered set manipulation](ordsets.html)

[A.35 library(pairs): Operations on key-value lists](pairs.html)

[A.36 library(persistency): Provide persistent dynamic predicates](persistency.html)

[A.37 library(pio): Pure I/O](pio.html)

[A.37.1 library(pure_input): Pure Input from files and streams](pio.html#sec:A.37.1)

[A.38 library(portray_text): Portray text](portraytext.html)

[A.39 library(predicate_options): Declare option-processing of predicates](predicate_options.html)

[A.39.1 The strength and weakness of predicate options](predicate_options.html#sec:A.39.1)

[A.39.2 Options as arguments or environment?](predicate_options.html#sec:A.39.2)

[A.39.3 Improving on the current situation](predicate_options.html#sec:A.39.3)

[A.39.3.1 Options as types](predicate_options.html#sec:A.39.3.1)

[A.39.3.2 Reflective access to options](predicate_options.html#sec:A.39.3.2)

[A.40 library(prolog_coverage): Coverage analysis tool](prologcoverage.html)

[A.40.1 Coverage collection and threads](prologcoverage.html#sec:A.40.1)

[A.40.2 Combining coverage data from multiple runs](prologcoverage.html#sec:A.40.2)

[A.40.3 Predicate reference](prologcoverage.html#sec:A.40.3)

[A.41 library(prolog_debug): User level debugging tools](prologdebug.html)

[A.42 library(prolog_jiti): Just In Time Indexing (JITI) utilities](prologjiti.html)

[A.43 library(prolog_trace): Print access to predicates](prologtrace.html)

[A.44 library(prolog_versions): Demand specific (Prolog) versions](prologversions.html)

[A.45 library(prolog_xref): Prolog cross-referencer data collection](prologxref.html)

[A.46 library(quasi_quotations): Define Quasi Quotation syntax](quasiquotations.html)

[A.47 library(random): Random numbers](random.html)

[A.48 library(rbtrees): Red black trees](rbtrees.html)

[A.49 library(readutil): Read utilities](readutil.html)

[A.50 library(record): Access named fields in a term](record.html)

[A.51 library(registry): Manipulating the Windows registry](registry.html)

[A.52 library(rwlocks): Read/write locks](rwlocks.html)

[A.53 library(settings): Setting management](settings.html)

[A.54 library(statistics): Get information about resource usage](statistics.html)

[A.55 library(strings): String utilities](strings.html)

[A.56 library(simplex): Solve linear programming problems](simplex.html)

[A.56.1 Introduction](simplex.html#sec:A.56.1)

[A.56.2 Delayed column generation](simplex.html#sec:A.56.2)

[A.56.3 Solving LPs with special structure](simplex.html#sec:A.56.3)

[A.56.4 Examples](simplex.html#sec:A.56.4)

[A.56.4.1 Example 1](simplex.html#sec:A.56.4.1)

[A.56.4.2 Example 2](simplex.html#sec:A.56.4.2)

[A.56.4.3 Example 3](simplex.html#sec:A.56.4.3)

[A.57 library(solution_sequences): Modify solution sequences](solutionsequences.html)

[A.58 library(tables): XSB interface to tables](tables.html)

[A.59 library(tableutil): Table inspection and statistics utilities](tableutil.html)

[A.60 library(terms): Term manipulation](terms.html)

[A.61 library(thread): High level thread primitives](thread.html)

[A.62 library(thread_pool): Resource bounded thread management](threadpool.html)

[A.63 library(ugraphs): Graph manipulation library](ugraphs.html)

[A.64 library(url): Analysing and constructing URL](url.html)

[A.65 library(varnumbers): Utilities for numbered terms](varnumbers.html)

[A.66 library(yall): Lambda expressions](yall.html)
