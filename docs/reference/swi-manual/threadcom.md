
## 10.3 Thread communication

### 10.3.1 Message queues

Prolog threads can exchange data using dynamic predicates, database records, and other globally shared data. These provide no suitable means to wait for data or a condition as they can only be checked in an expensive polling loop. *Message queues* provide a means for threads to wait for data or conditions without using the CPU.

Each thread has a message queue attached to it that is identified by the thread. Additional queues are created using [message_queue_create/1](threadcom.html#message_queue_create/1). Explicitly created queues come in two flavours. When given an *alias*, they must be destroyed by the user. *Anonymous* message queues are identified by a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)) and subject to garbage collection.

**thread_send_message**(`+QueueOrThreadId, +Term`)  
Place `Term` in the given queue or default queue of the indicated thread (which can even be the message queue of itself, see [thread_self/1](threadcreate.html#thread_self/1)). Any term can be placed in a message queue, but note that the term is copied to the receiving thread and variable bindings are thus lost. This call returns immediately.

If more than one thread is waiting for messages on the given queue and at least one of these is waiting with a partially instantiated `Term`, the waiting threads are *all* sent a wake-up signal, starting a rush for the available messages in the queue. This behaviour can seriously harm performance with many threads waiting on the same queue as all-but-the-winner perform a useless scan of the queue. If there is only one waiting thread or all waiting threads wait with an unbound variable, an arbitrary thread is restarted to scan the queue.^(201See the documentation for the POSIX thread functions **pthread_cond_signal()** v.s. **pthread_cond_broadcast()** for background information.)

\[semidet\]**thread_send_message**(`+Queue, +Term, +Options`)  
As [thread_send_message/2](threadcom.html#thread_send_message/2), but providing additional `Options`. These are to deal with the case that the queue has a finite maximum size and is full: whereas [thread_send_message/2](threadcom.html#thread_send_message/2) will block until the queue has drained sufficiently to accept a new message, [thread_send_message/3](threadcom.html#thread_send_message/3) can accept a time-out or deadline analogously to [thread_get_message/3](threadcom.html#thread_get_message/3). The options are:

**deadline**(`+AbsTime`)  
The call fails (silently) if no space has become available before `AbsTime`. See [get_time/1](system.html#get_time/1) for the representation of absolute time. If `AbsTime` is earlier then the current time, [thread_send_message/3](threadcom.html#thread_send_message/3) fails immediately. Both resolution and maximum wait time is platform-dependent.^(202The implementation uses **MsgWaitForMultipleObjects()** on MS-Windows and **pthread_cond_timedwait()** on other systems.)

**timeout**(`+Time`)  
`Time` is a float or integer and specifies the maximum time to wait in seconds. This is a relative-time version of the `deadline` option. If both options are provided, the earlier time is effective.

If `Time` is 0 or 0.0, [thread_send_message/3](threadcom.html#thread_send_message/3) examines the queue and sends the message if space is available, but does not suspend if no space is available, failing immediately instead.

If `Time` `< 0`, [thread_send_message/3](threadcom.html#thread_send_message/3) fails immediately without sending the message.

**signals**(`+BoolOrTime`)  
Whether or not signals (see [thread_signal/2](threadcom.html#thread_signal/2)) are processed while waiting. As the underlying implementation does not handle signals on most platforms, the implementation by default (`true`) times out every 0.25 seconds and checks for signals. If `false`, signals are not checked. If a number is specified, we check for signals every `Time` seconds. Smaller times may be used to improved responsiveness to signals. Larger times may be used to reduce CPU usage.

**thread_get_message**(`?Term`)  
Examines the thread message queue and if necessary blocks execution until a term that unifies to `Term` arrives in the queue. After a term from the queue has been unified to `Term`, the term is deleted from the queue.

Please note that non-unifying messages remain in the queue. After the following has been executed, thread 1 has the term `b(gnu)` in its queue and continues execution using `A` = `gnat`.

``` code
   <thread 1>
   thread_get_message(a(A)),

   <thread 2>
   thread_send_message(Thread_1, b(gnu)),
   thread_send_message(Thread_1, a(gnat)),
```

`Term` may contain attributed variables (see [section 8](clp.html#sec:8)), in which case only terms for which the constraints successfully execute are returned. Handle constraints applies for all predicates that extract terms from message queues. For example, we can get the even numbers from a queue using this code:

``` code
get_matching_messages(Queue, Pattern, [H|T]) :-
    copy_term(Pattern, H),
    thread_get_message(Queue, H, [timeout(0)]),
    !,
    get_matching_messages(Queue, Pattern, T).
get_matching_messages(_, _, []).

even_numbers(Q, List) :-
    freeze(Even, Even mod 2 =:= 0),
    get_matching_messages(Q, Even, List).
```

See also [thread_peek_message/1](threadcom.html#thread_peek_message/1).

**thread_peek_message**(`?Term`)  
Examines the thread message queue and compares the queued terms with `Term` until one unifies or the end of the queue has been reached. In the first case the call succeeds, possibly instantiating `Term`. If no term from the queue unifies, this call fails. I.e., [thread_peek_message/1](threadcom.html#thread_peek_message/1) never waits and does not remove any term from the queue. See also [thread_get_message/3](threadcom.html#thread_get_message/3).

**message_queue_create**(`?Queue`)  
Equivalent to `message_queue_create(Queue,[])`. For compatibility, calling `message_queue_create(+Atom)` is equivalent to `message_queue_create(Queue, [alias(Atom)])`. New code should use [message_queue_create/2](threadcom.html#message_queue_create/2) to create a named queue.

**message_queue_create**(`-Queue, +Options`)  
Create a message queue from `Options`. Defined options are:

**alias**(`+Alias`)  
Create a message queue that is identified by the atom `Alias`. Message queues created this way must be explicitly destroyed by the user. If the alias option is omitted, an *Anonymous* queue is created that is identified by a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)) and subject to garbage collection.^(203Garbage collecting anonymous message queues is not part of the ISO proposal and most likely not a widely implemented feature.)

**max_size**(`+Size`)  
Maximum number of terms in the queue. If this number is reached, [thread_send_message/2](threadcom.html#thread_send_message/2) will suspend until the queue is drained. The option can be used if the source, sending messages to the queue, is faster than the drain, consuming the messages.

\[det\]**message_queue_destroy**(`+Queue`)  
Destroy a message queue created with [message_queue_create/1](threadcom.html#message_queue_create/1). A permission error is raised if `Queue` refers to (the default queue of) a thread. Other threads that are waiting for `Queue` using [thread_get_message/2](threadcom.html#thread_get_message/2) receive an existence error.

\[semidet\]**is_message_queue**(`@Term`)  
True if `Term` refers to an existing message queue. This predicate can not block and has no error conditions. Note that message queues may be destroyed asynchronously by another thread and *anonymous* message queues may be garbage collected asynchronously.

\[det\]**thread_get_message**(`+Queue, ?Term`)  
As [thread_get_message/1](threadcom.html#thread_get_message/1), operating on a given queue. It is allowed (but not advised) to get messages from the queue of other threads. This predicate raises an existence error exception if `Queue` doesn't exist or is destroyed using [message_queue_destroy/1](threadcom.html#message_queue_destroy/1) while this predicate is waiting.

\[semidet\]**thread_get_message**(`+Queue, ?Term, +Options`)  
As [thread_get_message/2](threadcom.html#thread_get_message/2), but providing additional `Options`:

**deadline**(`+AbsTime`)  
The call fails (silently) if no message has arrived before `AbsTime`. See [get_time/1](system.html#get_time/1) for the representation of absolute time. If `AbsTime` is earlier then the current time, [thread_get_message/3](threadcom.html#thread_get_message/3) fails immediately. Both resolution and maximum wait time is platform-dependent.^(204The implementation uses **MsgWaitForMultipleObjects()** on MS-Windows and **pthread_cond_timedwait()** on other systems.)

**timeout**(`+Time`)  
`Time` is a float or integer and specifies the maximum time to wait in seconds. This is a relative-time version of the `deadline` option. If both options are provided, the earlier time is effective.

If `Time` is 0 or 0.0, [thread_get_message/3](threadcom.html#thread_get_message/3) examines the queue but does not suspend if no matching term is available. Note that unlike [thread_peek_message/2](threadcom.html#thread_peek_message/2), a matching term is removed from the queue.

If `Time` `< 0`, [thread_get_message/3](threadcom.html#thread_get_message/3) fails immediately without removing any message from the queue.

**signals**(`+BoolOrTime`)  
Whether or not signals (see [thread_signal/2](threadcom.html#thread_signal/2)) are processed while waiting. As the underlying implementation does not handle signals on most platforms, the implementation by default (`true`) times out every 0.25 seconds and checks for signals. If `false`, signals are not checked. If a number is specified, we check for signals every `Time` seconds. Smaller times may be used to improved responsiveness to signals. Larger times may be used to reduce CPU usage.

\[semidet\]**thread_peek_message**(`+Queue, ?Term`)  
As [thread_peek_message/1](threadcom.html#thread_peek_message/1), operating on a given queue. It is allowed to peek into another thread's message queue, an operation that can be used to check whether a thread has swallowed a message sent to it.

**message_queue_property**(`?Queue, ?Property`)  
True if `Property` is a property of `Queue`. Defined properties are:

**alias**(`Alias`)  
Queue has the given alias name.

**max_size**(`Size`)  
Maximum number of terms that can be in the queue. See [message_queue_create/2](threadcom.html#message_queue_create/2). This property is not present if there is no limit (default).

**size**(`Size`)  
Queue currently contains `Size` terms. Note that due to concurrent access the returned value may be outdated before it is returned. It can be used for debugging purposes as well as work distribution purposes.

**waiting**(`-Count`)  
Number of threads waiting for this queue. This property is not present if no threads waits for this queue.

The `size(Size)` property is always present and may be used to enumerate the created message queues. Note that this predicate does *not enumerate* threads, but can be used to query the properties of the default queue of a thread.

**message_queue_set**(`+Queue, +Property`)  
Set a property on the queue. Supported properties are:

**max_size**(`+Size`)  
Change the number of terms that may appear in the message queue before the next [thread_send_message/\[2,3\]](threadcom.html#thread_send_message/2) blocks on it. If the value is higher then the current maximum and the queue has writers waiting, unblock the writers. The value can be lower than the current number of terms in the queue. In that case writers will block until the queue is drained below the new maximum.

Explicit message queues are designed with the *worker-pool* model in mind, where multiple threads wait on a single queue and pick up the first goal to execute. Below is a simple implementation where the workers execute arbitrary Prolog goals. Note that this example provides no means to tell when all work is done. This must be realised using additional synchronisation.

``` code
%%      create_workers(?Id, +N)
%
%       Create a pool with Id and number of workers.
%       After the pool is created, post_job/1 can be used to
%       send jobs to the pool.

create_workers(Id, N) :-
        message_queue_create(Id),
        forall(between(1, N, _),
               thread_create(do_work(Id), _, [])).

do_work(Id) :-
        repeat,
          thread_get_message(Id, Goal),
          (   catch(Goal, E, print_message(error, E))
          ->  true
          ;   print_message(error, goal_failed(Goal, worker(Id)))
          ),
        fail.

%%      post_job(+Id, +Goal)
%
%       Post a job to be executed by one of the pool's workers.

post_job(Id, Goal) :-
        thread_send_message(Id, Goal).
```

### 10.3.2 Waiting for events

While message queues realizes communicating *agents* sharing the same program and optionally dynamic data, the predicate [thread_wait/2](threadcom.html#thread_wait/2) facilitates agents that communicate based on a *shared blackboard*. An important difference is were message queues require the sender and receiver to know about the queue used to communicate and every message can unblock at most one thread, the blackboard model allows any number (including zero) of threads to *listen* to changes on the blackboard. Any module can act as a blackboard. The blackboard can be updated using the standard Prolog database update predicates ([assert/1](db.html#assert/1), [retract/1](db.html#retract/1) and friends).

Waiting is implemented using a POSIX *condition variable* and matching *mutex*. On a matching database change the condition variable is signalled using a *broadcast*, unblocking all threads waiting in [thread_wait/2](threadcom.html#thread_wait/2). Multiple database updates can be grouped and cause a single unblock event using [thread_update/2](threadcom.html#thread_update/2). This predicate also allows signalling the module condition variable without updating the database and controlling whether all or a single thread is activated.

The blackboard architecture is a good match for an intelligent agent system that has to react on a changing world. Input threads gather sensor data from the world and update a shared world view in a set of dynamic predicates in one or more modules. Agent threads listen to this data or a subset thereof and trigger actions. This is notably a good match with *tabling*, in particular incremental tabling (see [section 7.7](tabling-incremental.html#sec:7.7)) and *Well Founded Semantics* (see [section 7.6](WFS.html#sec:7.6)).^(205Future versions may provide additional triggers, for example to learn about invalidated tables. Please share your experience.)

**thread_wait**(`:Goal, :Options`)  
Block execution of the calling thread until `Goal` becomes true. The application must be prepared to handle spurious calls to `Goal`, i.e., more calls than asked for based on the `Options` list. A possible exception in `Goal` is propagated and thus terminates [thread_wait/2](threadcom.html#thread_wait/2).

The wait is associated with a module. This module is derived from the `Options` argument.

The `Options` list specifies when `Goal` is re-evaluated and optionally when the call terminates due to a timeout.

**deadline**(`+AbsTime`)  
**timeout**(`+Time`)  
Timeout and deadline handling. See [thread_get_message/3](threadcom.html#thread_get_message/3) for details. This predicate fails when it terminates due to one of these options.

**retry_every**(`+Time`)  
Retry goal every `Time` seconds regardless of whether an event happened. The default is 1 second. This ensures that signals (see [thread_signal/2](threadcom.html#thread_signal/2)) and time limits are respected with an optional delay.^(206Some operating systems process such signals immediately, while others only check for such events synchronously.)

**db**(`+Boolean`)  
Wake up on arbitrary changes to any dynamic predicate that is defined in the associated module. This is the default if `wait_preds(+Preds)` is not provided.

**wait_preds**(`+List`)  
Only call `Goal` if at least one of the predicates in `List` has been modified. Each element of `List` is a *predicate indicator* (*Name/Arity* or *Name//Arity* that is resolved to a predicate in the module this wait is associated with. If the element is `+``(PI)`^(207Note that `+``p/1` is read as /(+(p),1).), `Goal` is only triggered if a clause was added ([assert/1](db.html#assert/1)). If the element is `-``(PI)`, `Goal` is only triggered if a clause was retracted ([retract/1](db.html#retract/1) or [erase/1](db.html#erase/1)). Default is to wake up on both assert and retract.

**modified**(`-List`)  
The `List` variable normally also appears in `Goal` and is unified with a list of predicates from the `wait_preds` option that have been modified. `List` must be unbound at entry.

**module**(`+Module`)  
Specifies the module to act on explicitly.

The execution of `Goal` is synchronized between all threads calling this predicate on the same module, changes to dynamic predicates in this module and calls to [thread_update/2](threadcom.html#thread_update/2) on the same module.

This predicate raises a `permision_error` exception when called recursively or called from inside a transaction. See [section 4.14.1.4](db.html#sec:4.14.1.4) for details about interaction with transactions.

**thread_update**(`:Goal, :Options`)  
Update a module (typically using [assert/1](db.html#assert/1) and/or [retract/1](db.html#retract/1) and friends) and on completion signal threads waiting for this module using [thread_wait/2](threadcom.html#thread_wait/2) to reevaluate their `Goal`. `Goal` is synchronized between updating and waiting threads. `Options`:

**module**(`+Module`)  
Determines the module to operate on. Default is the context module associated with the `Options` argument.

**notify**(`+Atom`)  
Determines whether all waiting threads are activated (`broadcast`, default) or a single thread (`signal`).

*Compatibility* The [thread_wait/2](threadcom.html#thread_wait/2) predicate is modelled after the [Qu-Prolog](http://staff.itee.uq.edu.au/pjr/HomePages/QuPrologHome.html) predicate thread_wait_on_goal/2. It is largely compatible. Our current implementation does not support predicate time stamps.^(208See [predicate_property/2](examineprog.html#predicate_property/2), property `generation`.) We made this predicate act on a specific module rather than the entire database. The timeout specification follows that of the other thread waiting predicates and may be combined with the `retry_every` option. The default retry-time is also 1 second rather than *infinite*.

### 10.3.3 Signalling threads

The predicates in this section provide *signalling* between threads. A thread signal inserts any goal as an *interrupt* into the control flow of any target thread. The target thread processes the goal at the first safe opportunity. The mechanism was introduced with two goals in mind: (1) running a goal inside a thread for debugging purposes such as enabling the status or get access thread-specific data and (2) force a thread to abort its current goal by inserting an exception into its control flow.

Over time, more complicated use cases have been identified that may result in multiple signals that occur (nearly) simultaneous. As of version 8.5.1 the interface has been extended and the interaction with other built-in predicates has been specified in much more detail.

\[det\]**thread_signal**(`+ThreadId, :Goal`)  
Make thread `ThreadId` execute `Goal` at the first opportunity. The predicate [thread_signal/2](threadcom.html#thread_signal/2) itself places `Goal` into the signalled thread's signal queue and returns immediately.

`ThreadId` executes `Goal` as an *interrupt* at the first opportunity. Defined opportunities are:

- At the *call port* of any predicate except for predicates with the property `sig_atomic`. Currently this only applies to [sig_atomic/1](threadcom.html#sig_atomic/1).
- Before retrying a foreign predicate.
- Before backtracking to the next clause of a Prolog predicate.
- When a foreign predicate calls [PL_handle_signals()](foreigninclude.html#PL_handle_signals()). Foreign predicates that take long to complete should call [PL_handle_signals()](foreigninclude.html#PL_handle_signals()) regularly and return with `FALSE` after [PL_handle_signals()](foreigninclude.html#PL_handle_signals()) returned -1, indicating an exception was raised.
- Foreign predicates calling *blocking system calls* should attempt to make these system calls interruptible. To enable this on POSIX systems, SWI-Prolog sends a `SIGUSR2` to the signalled thread while the handler is an empty function. This causes most blocking system calls to return with `EINTR`. See also the commandline option **--sig-alert**. On Windows, [PL_handle_signals()](foreigninclude.html#PL_handle_signals()) is called when the user processes Windows messages.
- For some blocking (thread) APIs we use a timed version with a 0.25 sec timeout to achieve a *polling loop*.

If one or more signals are queued, the queue is processed. Processing the queue skips signals blocked due to [sig_block/1](threadcom.html#sig_block/1) and stops after the queue does not contain any more non-blocked signals or processing a signal results in an exception. After an exception, other signals remain in the queue and will be processed after unwinding to the matching [catch/3](exception.html#catch/3). Typically these queued signals will be processed during the `Recover` goal of the [catch/3](exception.html#catch/3). Note that [sig_atomic/1](threadcom.html#sig_atomic/1) may be used to protect the recovery goal.

The [thread_signal/2](threadcom.html#thread_signal/2) mechanism is primarily used by the system to insert debugging goals into the target thread ([tspy/1](threadutil.html#tspy/1), [tbacktrace/1](threadutil.html#tbacktrace/1), etc.) or to interrupt a thread using e.g., `thread_signal(Thread, abort)`. Predicates from library `library(thread)` use signals to stop workers for e.g. [concurrent_maplist/2](thread.html#concurrent_maplist/2) if some call fails. Applications may use it, typically for similar purposes such as asynchronously stopping tasks or inspecting the status of a task. Below we describe the behaviour of thread signalling in more detail. The following notes apply for `Goal` executing in `ThreadId`

- The execution is protected by [sig_atomic/1](threadcom.html#sig_atomic/1) and thus signal execution is *not nested*.
- If `Goal` *succeeds*, possible choice points are discarded. Changes to the Prolog stacks such as changes to backtrackable global variables remain.
- If `Goal` *fails*, no action is taken, i.e., failure is not considered a special condition.
- If `Goal` *raises an exception* the exception is propagated into the environment. This allows for forcefully stopping the target thread. The system uses this to implement [abort/0](toplevel.html#abort/0) and call_with_time_limit/2.
- Code into which signals may be injected must make sure to use [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) and friends to ensure proper cleanup in the case of an exception. This is good practice anyway to guard against unpredictable exceptions such as resource exhaustion.
- `Goal` may use stack inspection such as [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3) to determine what the thread is doing.

If `Goal` is an integer, it is taken as a signal number and this signal is raised in the target `ThreadId`. Signal numbers range between 1 and 32. The predicate [current_signal/3](signal.html#current_signal/3) can be used to map signal names to signal numbers and inspect the handler.

\[det\]**sig_pending**(`-List`)  
True when `List` contains all signals submitted using [thread_signal/2](threadcom.html#thread_signal/2) that are not yet processed. This includes signals blocked by [sig_block/1](threadcom.html#sig_block/1).

\[det\]**sig_remove**(`:Pattern, -List`)  
Remove all signals that unify with `Pattern` from the signal queue and make the removed signals available in `List`

\[det\]**sig_block**(`:Pattern`)  
Block thread signals queued using [thread_signal/2](threadcom.html#thread_signal/2) that match `Pattern`.

\[det\]**sig_unblock**(`:Pattern`)  
Remove any effect of [sig_block/1](threadcom.html#sig_block/1) for patterns that are more specific (see [subsumes_term/2](compare.html#subsumes_term/2)). If any patterns are removed, reschedule blocked signals. Note that [sig_unblock/1](threadcom.html#sig_unblock/1) normally causes all unblocked signals to be executed immediately.

\[semidet\]**sig_atomic**(`:Goal`)  
Execute `Goal` as [once/1](metacall.html#once/1) while blocking both thread signals (see [thread_signal/2](threadcom.html#thread_signal/2)) and OS signals (see [on_signal/3](signal.html#on_signal/3)). The system executes some goals while blocking signals. These are:

- The goal injected using [thread_signal/2](threadcom.html#thread_signal/2), i.e., signals do not interrupt a running signal handler.
- The `Setup` call of [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) and friends.
- The `Cleanup` call of [call_cleanup/2](metacall.html#call_cleanup/2) and friends.
- Compiling a file or loading a *quick load file*.

The call port of [sig_atomic/1](threadcom.html#sig_atomic/1) does not handle signals. This may notably be used to prevent interruption of the [catch/3](exception.html#catch/3) `Recover` goal. For example, we may ensure the recovery goal of a timeout is called using the code below. Without this precaution another signal may run before [writeln/1](termrw.html#writeln/1) and raise an exception to prevent its execution. Note that [catch/3](exception.html#catch/3) should generally *not* be used for cleanup of resources in case of an exception and thus it is typically fine if its `Recover` goal is interrupted. Use [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) or one of the other predicates from the [call_cleanup/2](metacall.html#call_cleanup/2) family for cleanup.

``` code
    ...,
    catch(call_with_time_limit(Time, Goal),
          time_limit_exceeded,
          sig_atomic(writeln('Time limit exceeded'))).
```

### 10.3.4 Threads and dynamic predicates

Besides queues ([section 10.3.1](threadcom.html#sec:10.3.1)) threads can share and exchange data using dynamic predicates. The multithreaded version knows about two types of dynamic predicates. By default, a predicate declared *dynamic* (see [dynamic/1](dynamic.html#dynamic/1)) is shared by all threads. Each thread may assert, retract and run the dynamic predicate. Synchronisation inside Prolog guarantees the consistency of the predicate. Updates are *logical*: visible clauses are not affected by assert/retract after a query started on the predicate. In many cases primitives from [section 10.4](threadsync.html#sec:10.4) should be used to ensure that application invariants on the predicate are maintained.

Besides shared predicates, dynamic predicates can be declared with the [thread_local/1](threadcom.html#thread_local/1) directive. Such predicates share their attributes, but the clause list is different in each thread.

**thread_local** `+Functor/+Arity, ...`  
This directive is related to the [dynamic/1](dynamic.html#dynamic/1) directive. It tells the system that the predicate may be modified using [assert/1](db.html#assert/1), [retract/1](db.html#retract/1), etc., during execution of the program. Unlike normal shared dynamic data, however, each thread has its own clause list for the predicate. As a thread starts, this clause list is empty. If there are still clauses when the thread terminates, these are automatically reclaimed by the system (see also [volatile/1](saved-states.html#volatile/1)). The thread_local property implies the properties *dynamic* and *volatile*.

Thread-local dynamic predicates are intended for maintaining thread-specific state or intermediate results of a computation.

It is not recommended to put clauses for a thread-local predicate into a file, as in the example below, because the clause is only visible from the thread that loaded the source file. All other threads start with an empty clause list.

``` code
:- thread_local
        foo/1.

foo(gnat).
```

**DISCLAIMER** Whether or not this declaration is appropriate in the sense of the proper mechanism to reach the goal is still debated. If you have strong feelings in favour or against, please share them in the SWI-Prolog mailing list.
