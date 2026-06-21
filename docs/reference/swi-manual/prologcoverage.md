
## A.40 library(prolog_coverage): Coverage analysis tool

The purpose of this module is to find which part of the program has been used by a certain goal. Usage is defined in terms of clauses for which the *head unification* succeeded. For each clause we count how often it succeeded and how often it failed. In addition we track all *call sites*, creating goal-by-goal annotated clauses.

The result is represented as a list of clause-references. As the references to clauses of dynamic predicates cannot be guaranteed, these are omitted from the result.

Using [coverage/2](prologcoverage.html#coverage/2) with the option `annotate(true)`, implied by `ext(Ext)` or `dir(Dir)`, the analysis creates a line-by-line copy of the source files that is annotated with how many times this line was executed and with what logical results. These annotations rely on relating executable code to source locations which is shared by the source level debugger. Source level rewrites due to term or goal expansion may harm the results.

The typical usage is to load the program and run the query below to get a report by file with percentages and a directory `cov` holding annotated files that provide line-by-line annotations. See [show_coverage/1](prologcoverage.html#show_coverage/1) for details.

`?-` `coverage(Goal, [dir(cov)])`.

### A.40.1 Coverage collection and threads

The coverage collect data structure is shared by threads created from the thread that is collecting coverage data. Currently, this thread should be *joined* before we can operate on the coverage data.

### A.40.2 Combining coverage data from multiple runs

The coverage tools allow both combining data from running multiple queries as combining data from multiple Prolog processes.

For multiple queries in the same process, coverage data may be collected using [coverage/1](prologcoverage.html#coverage/1) which, unlike [coverage/2](prologcoverage.html#coverage/2), does not change the non-deterministic semantics of the `Goal` and adds to the already collected data. If no current collection is in progress, the currently collected data can be displayed using [show_coverage/1](prologcoverage.html#show_coverage/1).

Coverage data may be saved to a file using [cov_save_data/2](prologcoverage.html#cov_save_data/2). Saved data can be reloaded using [cov_load_data/2](prologcoverage.html#cov_load_data/2). Data from multiple Prolog runs can be combined in the same file using [cov_save_data/2](prologcoverage.html#cov_save_data/2) with the `append(true)` option. When possible, file locking is used to ensure that concurrect processes can safely use the same data file. The result can be shown by loading the code that was relevant to all runs, use [cov_load_data/2](prologcoverage.html#cov_load_data/2) and show the result using [show_coverage/1](prologcoverage.html#show_coverage/1).

Note that saving an loading the coverage data saves and restores references to the clauses as the Nth clause of a predicate defined in a specific file. This implies that the program must be loaded in exactly the same way, including optimization level, term/goal expansion and order of *multifile* predicates.

### A.40.3 Predicate reference

**coverage**(`:Goal`)  
As `call(Goal)`, collecting coverage information while `Goal` is running. If `Goal` succeeds with a choice point, coverage collection is suspended and resumed if we backtrack into `Goal`. Calls to [coverage/1](prologcoverage.html#coverage/1) may be nested.

\[semidet\]**coverage**(`:Goal, +Options`)  
Collect and optionally report coverage by `Goal`. `Goal` is executed as in [once/1](metacall.html#once/1). `Options` processed:

**show**(`+Boolean`)  
When `true` (default), call [show_coverage/1](prologcoverage.html#show_coverage/1) passing `Options` to show the collected coverage data and reset the data. When `false`, collect the data but do not reset it. If there is already existing data the new data is added.

\[det\]**show_coverage**(`+Options`)  
Show collected coverage data. By default it reports the percentage of called and failed clauses related to covered files. Using `dir(Dir)`, detailed line-by-line annotated files are created in the directory Dir. Other options control the level of detail.

**all**(`+Boolean`)  
When true, report on any file in which some predicate was called.

**modules**(`+Modules`)  
Only report on files that implement one of the given `Modules`.

**roots**(`+Directories`)  
Only report on files below one of the given roots. Each directory in `Directories` can be a specification for [absolute_file_name/3](files.html#absolute_file_name/3).

**annotate**(`+Bool`)  
Create an annotated file for the detailed results. This is implied if the `ext` or `dir` option are specified.

**ext**(`+Ext`)  
Extension to use for the annotated file. Default is‘.cov\`.

**dir**(`+Dir`)  
Dump the annotations in the given directory. If not given, the annotated files are created in the same directory as the source file. Each clause that is related to a physical line in the file is annotated with one of:

> |       |                                                  |
> |-------|--------------------------------------------------|
> | \###  | Clause was never executed.                       |
> | ++N   | Clause was entered N times and always succeeded  |
> | --N   | Clause was entered N times and never succeeded   |
> | +N-M  | Clause has succeeded N times and failed M times  |
> | +N\*M | Clause was entered N times and succeeded M times |

All *call sites* are annotated using the same conventions, except that `---` is used to annotate subgoals that were never called.

**line_numbers**(`Boolean`)  
If `true` (default), add line numbers to the annotated file.

**color**(`Boolean`)  
Controls using ANSI escape sequences to color the output in the annotated source. Default is `true`.

**width**(`+Columns`)  
Presumed width of the output window. A value of 40 is considered the minimum. Smaller values are handled as 40.

For example, run a goal and create annotated files in a directory `cov` using:

``` code
?- show_coverage([dir(cov)]).
```

bug  
Color annotations are created using ANSI escape sequences. On most systems these are displayed if the file is printed on the terminal. On most systems `less` may be used with the `-r` flag. Alternatively, programs such as `ansi2html` (Linux) may be used to convert the files to HTML. It would probably be better to integrate the output generation with `library(pldoc/doc_htmlsrc)`.

\[semidet,multifile\]**report_hook**(`+Succeeded, +Failed`)  
This hook is called after the data collection. It is passed a list of objects that have succeeded as well as a list of objects that have failed. The objects are one of

**`ClauseRef`**  
The specified clause

**call_site**(`ClauseRef, PC`)  
A call was make in `ClauseRef` at the given program counter.

\[det\]**cov_save_data**(`+File, +Options`)  
Save the coverage information to `File`. `Options`:

**append**(`true`)  
Append to `File` rather than truncating the data if the file exists.

The `File` is opened using `lock(exclusive)`, which implies that, provided the OS and file system implements file locking, multiple processes may save coverage data to the same file.

The saved data is highly specific to the setup in which it has been created. It can typically only be reloaded using [cov_load_data/2](prologcoverage.html#cov_load_data/2) in the same Prolog executable using the same options and with all relevant source file unmodified at the same location.

Reproducibility can be improved by using‘.qlf\` files or *saved states*.

\[det\]**cov_load_data**(`+File, +Options`)  
Reload coverage data from `File`. `Options`:

**load**(`true`)  
If specified and the file in which a clauses is expected to exist, load the file using [load_files/2](consulting.html#load_files/2) with the same options as used to initially load the file.

**silent**(`+Boolean`)  
When `true`, do not emit messages on not loaded source files.

Data is assumed to be reliable if the Nth-clause of a predicate is loaded from the same file at the same line number and has the same size. Unreliable data is ignored, silently if `silent(true)` is used.

\[det\]**cov_reset**  
Discard all collected coverage data. This predicate raises a permission error if coverage collection is in progress.

**cov_property**(`?Property`)  
True when coverage analysis satisfies `Property`. Currently defined properties are:

**active**(`?Nesting`)  
True when coverage data is being collected. `Nesting` expresses the nesting of [coverage/1](prologcoverage.html#coverage/1) calls and is normally 1 (one).
