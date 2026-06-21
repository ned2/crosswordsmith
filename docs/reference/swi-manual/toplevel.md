
## 4.37 User Top-level Manipulation

**break**  
Recursively start a new Prolog top level. This Prolog top level shares everything from the environment it was started in. Debugging is switched off on entering a break and restored on leaving one. The break environment is terminated by typing the system's end-of-file character (control-D). If that is somehow not functional, the term `end_of_file.` can be entered to return from the break environment. If the **-t** `toplevel` command line option is given, this goal is started instead of entering the default interactive top level ([prolog/0](toplevel.html#prolog/0)).

Notably the GUI based versions (**swipl-win** on Windows and MacOS) provide the menu **Run/New thread** that opens a new toplevel that runs concurrently with the initial toplevel. The concurrent toplevel can be used to examine the program, in particular global dynamic predicates. It can not access *global variables* or thread-local dynamic predicates (see [thread_local/1](threadcom.html#thread_local/1)) of the main thread.

**abort**  
Abort the Prolog execution and restart the top level. If the **-t** `toplevel` command line option is given, this goal is restarted instead of entering the default interactive top level.

Aborting is implemented by throwing the reserved exception `unwind(abort)`. This exception can be caught using [catch/3](exception.html#catch/3), but the recovery goal is wrapped with a predicate that prunes the choice points of the recovery goal (i.e., as [once/1](metacall.html#once/1)) and re-throws the exception. This is illustrated in the example below, where we press control-C and‘a’. See also [section 4.10.2](exception.html#sec:4.10.2).

``` code
?- catch((repeat,fail), E, true).
^CAction (h for help) ? abort
% Execution Aborted
```

\[ISO\]**halt**  
Terminate Prolog execution with default exit code using [halt/1](toplevel.html#halt/1). The default exit code is normally 0, but can be 1 if one of the Prolog flags [on_error](flags.html#flag:on_error) or [on_warning](flags.html#flag:on_warning) is set to `status` and there have been errors or warnings.

\[ISO\]**halt**(`+Status`)  
Terminate Prolog execution with `Status`. When possible, raise the exception `unwind(``halt(Status)``)`. Currently, this is used when [halt/1](toplevel.html#halt/1) is called in the `main` thread and there is no intermediate C function on the stack that called [PL_next_solution()](foreigninclude.html#PL_next_solution()) without the `PL_Q_PASS_EXCEPTION` flag. Future versions may also use signal based exit from threads.

After the exception bubbled up to the top or if the halt exception could not be raised, system termination starts. System termination (see also [PL_halt()](foreigninclude.html#PL_halt())) preforms the following steps:

1.  Set the Prolog flag [exit_status](flags.html#flag:exit_status) to `Status`.
2.  Call all hooks registered using [at_halt/1](consulting.html#at_halt/1). If `Status` equals 0 (zero), any of these hooks calls [cancel_halt/1](consulting.html#cancel_halt/1), termination is cancelled.
3.  Call all hooks registered using **PL_at_halt()**. In the future, if any of these hooks returns non-zero, termination will be cancelled. Currently, this only prints a warning.
4.  Perform the following system cleanup actions:
    - Raise `unwind(``halt(Status)``)` in all running threads.
    - Wait for a maximum of 1 second for all threads to respond to this exception. Registered [thread_at_exit/1](threadcreate.html#thread_at_exit/1) hooks are executed. Threads not responding within 1 second are cancelled forcefully.
    - Flush I/O and close all streams except for standard I/O.
    - Reset the terminal if its properties were changed.
    - Remove temporary files and incomplete compilation output.
    - Reclaim memory.
5.  Call exit(Status) to terminate the process

[halt/1](toplevel.html#halt/1) has been extended in SWI-Prolog to accept the arg `abort`. This performs as [halt/1](toplevel.html#halt/1) above except that:

- Termination cannot be cancelled with [cancel_halt/1](consulting.html#cancel_halt/1).
- **abort()** is called instead of exit(Status).

In addition to an integer status name we also allow passing a *signal name*. This is similar to `abort`, blocking halt cancellation and set the termination code to 128+signum. For example, using `halt(term)` the system exits with status 143 (= 128+15) on Linux.

**prolog**  
This goal starts the default interactive top level. Queries are read from the stream `user_input`. See also the Prolog flag **history**. The [prolog/0](toplevel.html#prolog/0) predicate is terminated (succeeds) by typing the end-of-file character (typically control-D).

The following hooks allow for expanding queries and handling the result of a query. These hooks are used by the top level variable expansion mechanism described in [section 2.9](topvars.html#sec:2.9).

**user:expand_query**(`+Query, -Expanded, +Bindings, -ExpandedBindings`)  
Hook in module `user`, normally not defined. `Query` and `Bindings` represents the query read from the user and the names of the free variables as obtained using [read_term/3](termrw.html#read_term/3). If this predicate succeeds, it should bind `Expanded` and `ExpandedBindings` to the query and bindings to be executed by the top level. This predicate is used by the top level ([prolog/0](toplevel.html#prolog/0)). See also expand_answer/2 and [term_expansion/2](consulting.html#term_expansion/2).

**prolog:expand_answer**(`+Goal, +Bindings, -ExpandedBindings`)  
Hook in module `prolog`, normally not defined. Expand the result of a successfully executed top-level query. `Bindings` is the query `<``Name``>=<``Value``>` binding list from the query. `ExpandedBindings` must be unified with the bindings the top level should print. `Goal` provides the instantiated query. This hook supersedes [user:expand_answer/2](toplevel.html#user:expand_answer/2).

\[deprecated\]**user:expand_answer**(`+Bindings, -ExpandedBindings`)  
Hook in module `user`, normally not defined. This hook provides backward compatibility and is superseded by [prolog:expand_answer/3](toplevel.html#prolog:expand_answer/3).
