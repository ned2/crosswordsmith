
## A.61 library(thread): High level thread primitives

author  
Jan Wielemaker

This module defines simple to use predicates for running goals concurrently. Where the core multi-threaded API is targeted at communicating long-living threads, the predicates here are defined to run goals concurrently without having to deal with thread creation and maintenance explicitely.

Note that these predicates run goals concurrently and therefore these goals need to be thread-safe. As the predicates in this module also abort branches of the computation that are no longer needed, predicates that have side-effect must act properly. In a nutshell, this has the following consequences:

- Nice clean Prolog code without side-effects (but with cut) works fine.
- Side-effects are bad news. If you really need assert to store intermediate results, use the [thread_local/1](threadcom.html#thread_local/1) declaration. This also guarantees cleanup of left-over clauses if the thread is cancelled. For other side-effects, make sure to use [call_cleanup/2](metacall.html#call_cleanup/2) to undo them should the thread be cancelled.
- Global variables are ok as they are thread-local and destroyed on thread cancellation. Note however that global variables in the calling thread are **not** available in the threads that are created. You have to pass the value as an argument and initialise the variable in the new thread.
- Thread-cancellation uses [thread_signal/2](threadcom.html#thread_signal/2). Using this code with long-blocking foreign predicates may result in long delays, even if another thread asks for cancellation.

\[semidet\]**concurrent**(`+N, :Goals, +Options`)  
Run `Goals` in parallel using `N` threads. This call blocks until all work has been done. The `Goals` must be independent. They should not communicate using shared variables or any form of global data. All `Goals` must be thread-safe.

Execution succeeds if all goals have succeeded. If one goal fails or throws an exception, other workers are abandoned as soon as possible and the entire computation fails or re-throws the exception. Note that if multiple goals fail or raise an error it is not defined which error or failure is reported.

On successful completion, variable bindings are returned. Note however that threads have independent stacks and therefore the goal is copied to the worker thread and the result is copied back to the caller of [concurrent/3](thread.html#concurrent/3).

Choosing the right number of threads is not always obvious. Here are some scenarios:

- If the goals are CPU intensive and normally all succeeding, typically the number of CPUs is the optimal number of threads. Less does not use all CPUs, more wastes time in context switches and also uses more memory.
- If the tasks are I/O bound the number of threads is typically higher than the number of CPUs.
- If one or more of the goals may fail or produce an error, using a higher number of threads may find this earlier.

|  |  |
|----|----|
| `N` | Number of worker-threads to create. Using 1, no threads are created. If `N` is larger than the number of `Goals` we create exactly as many threads as there are `Goals`. |
| `Goals` | List of callable terms. |
| `Options` | Passed to [thread_create/3](threadcreate.html#thread_create/3) for creating the workers. Only options changing the stack-sizes can be used. In particular, do not pass the detached or alias options. |

See also  
In many cases, [concurrent_maplist/2](thread.html#concurrent_maplist/2) and friends is easier to program and is tractable to program analysis.

\[semidet\]**concurrent_forall**(`:Generate, :Action`)  
\[semidet\]**concurrent_forall**(`:Generate, :Action, +Options`)  
True when `Action` is true for all solutions of `Generate`. This has the same semantics as [forall/2](forall2.html#forall/2), but the `Action` goals are executed in multiple threads. Notable a failing `Action` or a `Action` throwing an exception signals the calling thread which in turn aborts all workers and fails or re-throws the generated error. `Options`:

**threads**(`+Count`)  
Number of threads to use. The default is determined by the Prolog flag [cpu_count](flags.html#flag:cpu_count).

To be done  
Ideally we would grow the set of workers dynamically, similar to dynamic scheduling of HTTP worker threads. This would avoid creating threads that are never used if `Generate` is too slow or does not provide enough answers and would further raise the number of threads if `Action` is I/O bound rather than CPU bound.

**concurrent_and**(`:Generator, :Test`)  
**concurrent_and**(`:Generator, :Test, +Options`)  
Concurrent version of `(Generator,Test)`. This predicate creates a thread providing solutions for `Generator` that are handed to a pool of threads that run `Test` for the different instantiations provided by `Generator` concurrently. The predicate is logically equivalent to a simple conjunction except for two aspects: (1) terms are *copied* from `Generator` to the test `Test` threads while answers are copied back to the calling thread and (2) answers may be produced out of order.

If the evaluation of some `Test` raises an exception, [concurrent_and/2](thread.html#concurrent_and/2),3 is terminated with this exception. If the caller commits after a given answer or raises an exception while [concurrent_and/2](thread.html#concurrent_and/2),3 is active with pending choice points, all involved resources are reclaimed.

`Options`:

**threads**(`+Count`)  
Create a worker pool holding `Count` threads. The default is the Prolog flag [cpu_count](flags.html#flag:cpu_count).

This predicate was proposed by Jan Burse as `balance((Generator,Test))`.

\[semidet\]**concurrent_maplist**(`:Goal, +List`)  
\[semidet\]**concurrent_maplist**(`:Goal, +List1, +List2`)  
\[semidet\]**concurrent_maplist**(`:Goal, +List1, +List2, +List3`)  
Concurrent version of [maplist/2](apply.html#maplist/2). This predicate uses [concurrent/3](thread.html#concurrent/3), using multiple *worker* threads. The number of threads is the minimum of the list length and the number of cores available. The number of cores is determined using the prolog flag `cpu_count`. If this flag is absent or 1 or `List` has less than two elements, this predicate calls the corresponding maplist/N version using a wrapper based on [once/1](metacall.html#once/1). Note that all goals are executed as if wrapped in [once/1](metacall.html#once/1) and therefore these predicates are *semidet*.

Note that the the overhead of this predicate is considerable and therefore `Goal` must be fairly expensive before one reaches a speedup.

\[semidet\]**first_solution**(`-X, :Goals, +Options`)  
Try alternative solvers concurrently, returning the first answer. In a typical scenario, solving any of the goals in `Goals` is satisfactory for the application to continue. As soon as one of the tried alternatives is successful, all the others are killed and [first_solution/3](thread.html#first_solution/3) succeeds.

For example, if it is unclear whether it is better to search a graph breadth-first or depth-first we can use:

``` code
search_graph(Grap, Path) :-
         first_solution(Path, [ breadth_first(Graph, Path),
                                depth_first(Graph, Path)
                              ],
                        []).
```

`Options` include thread stack-sizes passed to thread_create, as well as the options `on_fail` and `on_error` that specify what to do if a solver fails or triggers an error. By default execution of all solvers is terminated and the result is returned. Sometimes one may wish to continue. One such scenario is if one of the solvers may run out of resources or one of the solvers is known to be incomplete.

**on_fail**(`Action`)  
If `stop` (default), terminate all threads and stop with the failure. If `continue`, keep waiting.

**on_error**(`Action`)  
As above, re-throwing the error if an error appears.

bug  
[first_solution/3](thread.html#first_solution/3) cannot deal with non-determinism. There is no obvious way to fit non-determinism into it. If multiple solutions are needed wrap the solvers in [findall/3](allsolutions.html#findall/3).

\[semidet\]**call_in_thread**(`+Thread, :Goal`)  
\[semidet\]**call_in_thread**(`+Thread, :Goal, +Options`)  
Run `Goal` as an interrupt in the context of `Thread`. This is based on [thread_signal/2](threadcom.html#thread_signal/2). If waiting times out, we inject a `stop(Reason)` exception into `Goal`. Interrupts can be nested, i.e., it is allowed to run a [call_in_thread/2](thread.html#call_in_thread/2) while the target thread is processing such an interrupt.

`Options` are passed to [thread_get_message/3](threadcom.html#thread_get_message/3) and notably allow for specifying a timeout. If a timeout is reached, this predicate will attempt to kill `Goal` in `Thread` and act according to the option `on_timeout`.

**on_timeout**(`:Goal`)  
If waiting terminates due to a `timeout(Time)`, or `deadline(Stamp)` option, call `Goal`. The default is `throw(time_limit_exceeded)`.

This predicate is primarily intended for debugging and inspection tasks.
