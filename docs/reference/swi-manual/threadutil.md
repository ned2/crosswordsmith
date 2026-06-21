
## 10.5 library(threadutil): Interactive thread utilities

This library provides utilities that are primarily intended for interactive usage in a threaded Prolog environment. It allows for inspecting threads, manage I/O of background threads (depending on the environment) and manipulating the debug status of threads.

**threads**  
List currently known threads with their status. For each thread it lists the *id*, *class*, *status*, *debug status*, *CPU time* and current *stack usage*. If a thread is listed with *Debug* set to‘X’, it cannot be debugged. If the *Debug* is show as‘V’, it is running in debug mode (see [debug/0](debugger.html#debug/0)) and responds to *spy points* and *break points*.

**join_threads**  
Join all terminated threads. For normal applications, dealing with terminated threads must be part of the application logic, either detaching the thread before termination or making sure it will be joined. The predicate [join_threads/0](threadutil.html#join_threads/0) is intended for interactive sessions to reclaim resources from threads that died unexpectedly during development.

\[det\]**with_stopped_threads**(`:Goal, Options`)  
Stop all threads except the caller while running `once(Goal)`. Note that this is in the thread user utilities as this is not something that should be used by normal applications. Notably, this may *deadlock* if the current thread requires input from some other thread to complete `Goal` or one of the stopped threads has a lock. `Options`:

**stop_nodebug_threads**(`+Boolean`)  
If `true` (default `false`), also stop threads created with the `debug(false)` option.

**except**(`+List`)  
Do not stop threads from this list.

bug  
Note that the threads are stopped when they process signals. As signal handling may be delayed, this implies they need not be stopped before `Goal` starts.

\[semidet\]**thread_has_console**  
True when the calling thread has an attached console.

See also  
[attach_console/0](threadutil.html#attach_console/0)

\[det\]**attach_console**  
\[det\]**attach_console**(`+Title`)  
Create a new console and make the standard Prolog streams point to it. If not provided, the title is built using the thread id. Does nothing if the current thread already has a console attached.

\[det\]**tspy**(`:Spec`)  
\[det\]**tspy**(`:Spec, +ThreadOrClass`)  
Trap the graphical debugger on reaching `Spec`. The predicate tspy/0 enabled debug mode in all threads using [tdebug/0](threadutil.html#tdebug/0) while [tspy/1](threadutil.html#tspy/1) enables debug mode using [tdebug/1](threadutil.html#tdebug/1).

\[det\]**tdebug**  
\[det\]**tdebug**(`+ThreadOrClass`)  
\[det\]**tnodebug**  
\[det\]**tnodebug**(`+ThreadOrClass`)  
Enable or disable a thread or group of threads for debugging using the graphical tracer. A group of threads is addressed based on the `class` property of a thread set by [thread_create/3](threadcreate.html#thread_create/3) or [set_thread/2](threadcreate.html#set_thread/2). This implies loading the graphical tracer and switching the thread to debug mode using [debug/0](debugger.html#debug/0). New threads created inherit their debug mode from the thread that created them.

Thread classes have been introduced in SWI-Prolog 10.0.2/10.1.5. This allows for more selective debugging as well as ensuring debugging works in newly created threads. For example, the HTTP server creates all its *worker threads* in the class `http`. Using query below, we reliable make sure spy points are trapped in HTTP handler threads, regardless of whether the worker existed or is lazily created and regardless of whether the user switched to *nodebug* mode while tracing a previous event (see debug_reset_from_class/0).

``` code
?- tdebug(http).
```

\[det\]**tbacktrace**(`+Thread`)  
\[det\]**tbacktrace**(`+Thread, +Options`)  
Print a backtrace for `Thread` to the stream `user_error` of the calling thread. This is achieved by inserting an interrupt into `Thread` using [call_in_thread/2](thread.html#call_in_thread/2). `Options`:

**depth**(`+MaxFrames`)  
Number of stack frames to show. Default is the current Prolog flag `backtrace_depth` or 20.

Other options are passed to get_prolog_backtrace/3.

bug  
[call_in_thread/2](thread.html#call_in_thread/2) may not process the event.

\[det\]**tprofile**(`+Thread`)  
Profile the operation of `Thread` until the user hits a key.

\[det\]**thread_alias**(`+Alias`)  
Set the alias for a thread.

deprecated  
Use [set_thread/2](threadcreate.html#set_thread/2) using `alias(Alias)`.
