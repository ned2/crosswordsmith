
## 10.1 Creating and destroying Prolog threads

**thread_create**(`:Goal, -Id`)  
Shorthand for `thread_create(Goal, Id, [])`. See [thread_create/3](threadcreate.html#thread_create/3).

**thread_create**(`:Goal, -Id, +Options`)  
Create a new Prolog thread (and underlying operating system thread) and start it by executing `Goal`. If the thread is created successfully, the thread identifier of the created thread is unified to `Id`.

`Id` is the *alias* name if the option `alias(name)` is given. Otherwise it is a *blob* of type `thread`. The anonymous blobs are subject to atom garbage collection. If a thread handle is garbage collected and the thread is not *detached*, it is *joined* if it has already terminated (see [thread_join/2](threadcreate.html#thread_join/2)) and detached otherwise (see [thread_detach/1](threadcreate.html#thread_detach/1)).^(198Up to version 7.3.23, anonymous thread handles were integers. Using integers did not allow for safe checking of the thread's status as the thread may have died and the handle may have been reused and did not allow for garbage collection to take care of forgotten threads.) The thread identifier blobs are printed as `<thread>(``I`,`Ptr``)`, where `I` is the internal thread identifier and `Ptr` is the unique address of the identifier. The `I` is accepted as input argument for all thread APIs that accept a thread identifier for convenient interaction from the toplevel. See also [thread_property/2](thmonitor.html#thread_property/2).

`Options` is a list of options. The currently defined options are below. Stack size options can also take the value `inf` or `infinite`, which is mapped to the maximum stack size supported by the platform.

**affinity**(`+CpuSet`)  
Specify that the thread should only run on the specified CPUs (cores). `CpuSet` is a list of integers between 0 (zero) and the known number of CPUs (see [cpu_count](flags.html#flag:cpu_count)). If `CpuSet` is empty a `domain_error` is raised. Referring to CPUs equal to or higher than the known number of CPUs returns an `existence_error`.

This option is currently implemented for systems that provide **pthread_attr_setaffinity_np()**. The option is silently ignored on other systems.^(bugThere is currently no way to discover whether this option is supported.)

**alias**(`AliasName`)  
Associate an‘alias name’with the thread. This name may be used to refer to the thread and remains valid until the thread is joined (see [thread_join/2](threadcreate.html#thread_join/2)). If the OS supports it (e.g., Linux), the operating system thread is named as well.

**at_exit**(`:AtExit`)  
Register `AtExit` as using [thread_at_exit/1](threadcreate.html#thread_at_exit/1) before entering the thread goal. Unlike calling [thread_at_exit/1](threadcreate.html#thread_at_exit/1) as part of the normal `Goal`, this *ensures* the `AtExit` is called. Using [thread_at_exit/1](threadcreate.html#thread_at_exit/1), the thread may be signalled or run out of resources before [thread_at_exit/1](threadcreate.html#thread_at_exit/1) is reached. See [thread_at_exit/1](threadcreate.html#thread_at_exit/1) for details.

**debug**(`+Bool`)  
Enable/disable debugging the new thread. If `false` (default `true`), the new thread is created with the property `debug(false)` and debugging is disabled before the new thread is started. The thread debugging predicates such as [tspy/1](threadutil.html#tspy/1) and [tdebug/0](threadutil.html#tdebug/0) do not signal threads with the debug property set to `false`.^(199Currently, the flag is only used as a hint for the various debugging primitives, i.e., the system does not really enforce that the target thread stays in `nodebug` mode.)

**detached**(`Bool`)  
If `false` (default), the thread can be waited for using [thread_join/2](threadcreate.html#thread_join/2). [thread_join/2](threadcreate.html#thread_join/2) must be called on this thread to reclaim all resources associated with the thread. If `true`, the system will reclaim all associated resources automatically after the thread finishes. Please note that thread identifiers are freed for reuse after a detached thread finishes or a normal thread has been joined. See also [thread_join/2](threadcreate.html#thread_join/2) and [thread_detach/1](threadcreate.html#thread_detach/1).

If a detached thread dies due to failure or exception of the initial goal, the thread prints a message using [print_message/2](printmsg.html#print_message/2). If such termination is considered normal, the code must be wrapped using [ignore/1](metacall.html#ignore/1) and/or [catch/3](exception.html#catch/3) to ensure successful completion.

**class**(`+Atom`)  
Register the thread as belonging to the given class. This is used by the debugger to control *debug mode*.

**inherit_from**(`+ThreadId`)  
Inherit defaults from the given `ThreadId` instead of the calling thread. This option was added to ensure that the `__thread_pool_manager` (see [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4)), which is created lazily, has a predictable state. The following properties are inherited:

- The prompt (see [prompt/2](termrw.html#prompt/2))
- The *typein* module (see [module/1](mtoplevel.html#module/1))
- The standard streams (`user_input`, etc.).
- `current_input` and `current_output` are bound to `user_input` and `user_output`.
- The default encoding (see [encoding](flags.html#flag:encoding))
- The default locale (see [set_locale/1](locale.html#set_locale/1))
- All Prolog flags
- The stack limit (see Prolog flag [stack_limit](flags.html#flag:stack_limit)).

**queue_max_size**(`Size`)  
Enforces a maximum to the number of terms in the input queue. See [message_queue_create/2](threadcom.html#message_queue_create/2) with the `max_size(o)`ption for details.

**stack_limit**(`Bytes`)  
Set the size limit for the Prolog stacks. See the Prolog flag [stack_limit](flags.html#flag:stack_limit). The default is inherited from the calling thread or the thread specified using `inherit_from(ThreadId)`.

**c_stack**(`Bytes`)  
Set the limit to which the C stack of this thread may grow. The default, minimum and maximum values are system-dependent. The value is rounded up to the system page size and SWI-Prolog enforces a minimum of 64 K-bytes.

The `Goal` argument is *copied* to the new Prolog engine. This implies that further instantiation of this term in either thread does not have consequences for the other thread: Prolog threads do not share data from their stacks.

**thread_self**(`-Id`)  
Get the Prolog thread identifier of the running thread. If the thread has an alias, the alias name is returned.

**thread_join**(`+Id`)  
Calls [thread_join/2](threadcreate.html#thread_join/2) and succeeds if thread `Id` terminated with success. Otherwise the exception `error(``thread_error(Id, Status)``, _)` is raised, where `Status` is the status as returned by [thread_join/2](threadcreate.html#thread_join/2).

**thread_join**(`+Id, -Status`)  
Wait for the termination of the thread with the given `Id`. Then unify the result status of the thread with `Status`. After this call, `Id` becomes invalid and all resources associated with the thread are reclaimed. It is not allowed for two threads to join the same thread and the thread being joined cannot be *detached* (see the `detached(true)` option for [thread_create/3](threadcreate.html#thread_create/3) and [thread_detach/1](threadcreate.html#thread_detach/1)).

A thread that has been completed without [thread_join/2](threadcreate.html#thread_join/2) being called on it is partly reclaimed: the Prolog stacks are released and the C thread is destroyed. A small data structure representing the exit status of the thread is retained until [thread_join/2](threadcreate.html#thread_join/2) is called on the thread. Defined values for `Status` are:

**true**  
The goal has been proven successfully.

**false**  
The goal has failed.

**exception**(`Term`)  
The thread is terminated on an exception. See [print_message/2](printmsg.html#print_message/2) to turn system exceptions into readable messages.

**exited**(`Term`)  
The thread is terminated on [thread_exit/1](threadcreate.html#thread_exit/1) using the argument `Term`.

Note that the pthread primitive **pthread_join()** cannot be interrupted. Some systems provide **pthread_timedjoin_np()**. If this is provided [thread_join/2](threadcreate.html#thread_join/2) is implemented as a loop of timed joins with a 0.25 sec timeout while signals are being tested after each timeout. Such systems allow combining [thread_join/2](threadcreate.html#thread_join/2) with call_with_time_limit/2 as well as signalling threads blocked in [thread_join/2](threadcreate.html#thread_join/2) using [thread_signal/2](threadcom.html#thread_signal/2).

**set_thread**(`+Thread, +Property`)  
Set a property for `Thread`. The currently defined properties are below.

**alias**(`+Atom`)  
Set the alias name of the thread. An error is raised if `Thread` already has an alias or `Alias` is in use for a thread or message queue.

**debug**(`+Bool`)  
Set whether `Thread` can be debugged. See also [thread_create/3](threadcreate.html#thread_create/3) and [thread_property/2](thmonitor.html#thread_property/2).

**debug_mode**(`+Bool`)  
Set or clear the Prolog flag [debug](flags.html#flag:debug) in `Thread`

**class**(`+Atom`)  
Set the *class* of `Thread`. This may enable or disable the debug mode of the thread depending on whether or not debug mode is enabled for the specified thread class.

**thread_detach**(`+Id`)  
Switch thread into detached state (see `detached(Bool)` option at [thread_create/3](threadcreate.html#thread_create/3)) at runtime. `Id` is the identifier of the thread placed in detached state. This may be the result of [thread_self/1](threadcreate.html#thread_self/1).

One of the possible applications is to simplify debugging. Threads that are created as *detached* leave no traces if they crash. For non-detached threads the status can be inspected using [thread_property/2](thmonitor.html#thread_property/2). Threads nobody is waiting for may be created normally and detach themselves just before completion. This way they leave no traces on normal completion and their reason for failure can be inspected.

**thread_exit**(`+Term`)  
Terminates the thread immediately, leaving `exited(Term)` as result state for [thread_join/2](threadcreate.html#thread_join/2). If the thread has the attribute `detached(true)` it terminates, but its exit status cannot be retrieved using [thread_join/2](threadcreate.html#thread_join/2), making the value of `Term` irrelevant. The Prolog stacks and C thread are reclaimed.

The current implementation is based on the reserved `unwind(``thread_exit(Term)``)` exception. This implies that, unlike the previous implementation that was based on the C **pthread_exit()** function, the implementation is safe from the Prolog point of view. However, it is limited by the semantics of the *unwind exceptions*. See [section 4.10.1](exception.html#sec:4.10.1) for details.

This predicate raises a `permission_error` if it is known that the thread cannot handle this case.

**thread_initialization**(`:Goal`)  
Run `Goal` when thread is started. This predicate is similar to [initialization/1](consulting.html#initialization/1), but is intended for initialization operations of the runtime stacks, such as setting global variables as described in [section 4.33](gvar.html#sec:4.33). `Goal` is run on four occasions: at the call to this predicate, after loading a saved state, on starting a new thread and on creating a Prolog engine through the C interface. On loading a saved state, `Goal` is executed *after* running the [initialization/1](consulting.html#initialization/1) hooks.

**thread_at_exit**(`:Goal`)  
Run `Goal` just before releasing the thread resources. This is to be compared to [at_halt/1](consulting.html#at_halt/1), but only for the current thread. These hooks are run regardless of why the execution of the thread has been completed. When these hooks are run, the return code is already available through [thread_property/2](thmonitor.html#thread_property/2) using the result of [thread_self/1](threadcreate.html#thread_self/1) as thread identifier. Note that there are two scenarios for using exit hooks. Using [thread_at_exit/1](threadcreate.html#thread_at_exit/1) is typically used if the thread creates a side-effect that must be reverted if the thread dies. Another scenario is where the creator of the thread wants to be informed when the thread ends. That cannot be guaranteed by means of [thread_at_exit/1](threadcreate.html#thread_at_exit/1) because it is possible that the thread cannot be created or dies almost instantly due to a signal or resource error. The `at_exit(Goal)` option of [thread_create/3](threadcreate.html#thread_create/3) is designed to deal with this scenario.

The `Goal` is executed with signal processing disabled. This avoids that e.g., `thread_signal(Thread, abort)` kills the exit handler rather than the thread in the case the body of `Thread` has just finished when the signal arrives.

**thread_setconcurrency**(`-Old, +New`)  
Determine the concurrency of the process, which is defined as the maximum number of concurrently active threads.‘Active’here means they are using CPU time. This option is provided if the thread implementation provides **pthread_setconcurrency()**. On other systems this predicate unifies `Old` to 0 (zero) and succeeds silently.

**thread_affinity**(`+ThreadID, -Current, +New`)  
True when `Current` is unified with the current thread affinity and the thread affinity is successfully set to `New`. The *thread affinity* specifies the set of CPUs on which this thread is allowed to run. The affinity is represented as a list of non-negative integers. See also the option `affinity(+Affinity)` of [thread_create/3](threadcreate.html#thread_create/3).

This predicate is only present if this functionality can be supported and has been ported to the target operating system. Currently, only Linux support is provided.
