
## 10.2 Monitoring threads

Normal multithreaded applications should not need the predicates from this section because almost any usage of these predicates is unsafe. For example checking the existence of a thread before signalling it is of no use as it may vanish between the two calls. Catching exceptions using [catch/3](exception.html#catch/3) is the only safe way to deal with thread-existence errors.

These predicates are provided for diagnosis and monitoring tasks. See also [section 10.5](threadutil.html#sec:10.5), describing more high-level primitives.

**is_thread**(`@Term`)  
True if `Term` is a handle to an existing thread.

**thread_property**(`?Id, ?Property`)  
True if thread `Id` has `Property`. Either or both arguments may be unbound, enumerating all relations on backtracking. Calling [thread_property/2](thmonitor.html#thread_property/2) does not influence any thread. See also [thread_join/2](threadcreate.html#thread_join/2). For threads that have an alias name, this name is returned in `Id` instead of the opaque thread identifier. Defined properties are:

**alias**(`Alias`)  
`Alias` is the alias name of thread `Id`.

**debug**(`Boolean`)  
Whether the thread is a *debug* thread. No-debug threads may be created using the `debug(false)` option of [thread_create/3](threadcreate.html#thread_create/3).

**debug_mode**(`Boolean`)  
Obtain the value of the Prolog flag [debug](flags.html#flag:debug) for the thread.

**class**(`Atom`)  
True when the thread's class is set to `Atom`. See [thread_create/3](threadcreate.html#thread_create/3).

**detached**(`Boolean`)  
Current detached status of the thread.

**id**(`Integer`)  
Integer identifier for the thread. Can be used as argument to the thread predicates, but applications must be aware that these references are reused.

**status**(`Status`)  
Current status of the thread. `Status` is one of:

**running**  
The thread is running. This is the initial status of a thread. Please note that threads waiting for something are considered running too.

**suspended**  
Only if the thread is an engine (see [section 11](engines.html#sec:11)). Indicates that the engine is currently not associated with an OS thread.

**false**  
The `Goal` of the thread has been completed and failed.

**true**  
The `Goal` of the thread has been completed and succeeded.

**exited**(`Term`)  
The `Goal` of the thread has been terminated using [thread_exit/1](threadcreate.html#thread_exit/1) with `Term` as argument. If the underlying native thread has exited (using **pthread_exit()**) `Term` is unbound.

**exception**(`Term`)  
The `Goal` of the thread has been terminated due to an uncaught exception (see [throw/1](exception.html#throw/1) and [catch/3](exception.html#catch/3)).

**engine**(`Boolean`)  
If the thread is an engine (see [chapter 11](engines.html#sec:11)), `Boolean` is `true`. Otherwise the property is not present.

**thread**(`ThreadId`)  
If the thread is an engine that is currently attached to a thread, `ThreadId` is the thread that executes the engine.

**size**(`Bytes`)  
The amount of memory associated with this thread. This includes the thread structure, its stacks, its default message queue, its clauses in its thread local dynamic predicates (see [thread_local/1](threadcom.html#thread_local/1)) and memory used for representing thread-local answer tries (see [section 7](tabling.html#sec:7)).

**system_thread_id**(`Integer`)  
Thread identifier used by the operating system for the calling thread. Not available on all OSes. This is the same as the Prolog flag [system_thread_id](flags.html#flag:system_thread_id) for the calling thread. Access to the system thread identifier can, on some systems, be used to gain additional control over or information about Prolog threads.

See also [thread_statistics/3](thmonitor.html#thread_statistics/3) to obtain resource usage information and [message_queue_property/2](threadcom.html#message_queue_property/2) to get the number of queued messages for a thread.

**thread_statistics**(`+Id, +Key, -Value`)  
Obtains statistical information on thread `Id` as [statistics/2](builtin-statistics.html#statistics/2) does in single-threaded applications. This call supports all keys of [statistics/2](builtin-statistics.html#statistics/2), although only stack sizes, `cputime`, `inferences`, `epoch`, `errors` and `warnings` yield different values for each thread. For `errors` and `warnings` [statistics/2](builtin-statistics.html#statistics/2) gives the global process count and this predicate gives the counts for the calling thread.^(200There is no portable interface to obtain thread-specific CPU time and some operating systems provide no access to this information at all. On such systems the total process CPU is returned. Thread CPU time is supported on MS-Windows, Linux and MacOSX.)

**mutex_statistics**  
Print usage statistics on internal mutexes and mutexes associated with dynamic predicates. For each mutex two numbers are printed: the number of times the mutex was acquired and the number of *collisions*: the number of times the calling thread has to wait for the mutex. The output is written to `current_output` and can thus be redirected using [with_output_to/2](IO.html#with_output_to/2).
