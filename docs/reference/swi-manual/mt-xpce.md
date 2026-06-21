
## 10.7 Multithreading and the XPCE graphics system

GUI applications written in XPCE can benefit from Prolog threads if they need to do expensive computations that would otherwise block the UI. The XPCE message passing system is guarded with a single *mutex*, which synchronises both access from Prolog and activation through the GUI. In MS-Windows, GUI events are processed by the thread that created the window in which the event occurred, whereas in Unix/X11 they are processed by the thread that dispatches messages. In practice, the most feasible approach to graphical Prolog implementations is to control XPCE from a single thread and deploy other threads for (long) computations.

Traditionally, XPCE runs in the foreground (`main`) thread. We are working towards a situation where XPCE can run comfortably in a separate thread. A separate XPCE thread can be created using [pce_dispatch/1](mt-xpce.html#pce_dispatch/1). It is also possible to create this thread as the `library(()`pce) is loaded by setting the **xpce_threaded** to `true`.

Threads other than the thread in which XPCE runs are provided with two predicates to communicate with XPCE.

\[det\]**in_pce_thread**(`:Goal`)  
Assuming XPCE is running in the foreground thread, this call gives background threads the opportunity to make calls to the XPCE thread. A call to [in_pce_thread/1](mt-xpce.html#in_pce_thread/1) succeeds immediately, copying `Goal` to the XPCE thread. `Goal` is added to the XPCE event queue and executed synchronous to normal user events like typing and clicking.

\[semidet\]**in_pce_thread_sync**(`:Goal`)  
Same as [in_pce_thread/1](mt-xpce.html#in_pce_thread/1), but wait for `Goal` to be completed. Success depends on the success of executing `Goal`. Variable bindings inside `Goal` are visible to the caller, but it should be noted that the values are being *copied*. If `Goal` throws an exception, this exception is re-thrown by [in_pce_thread/1](mt-xpce.html#in_pce_thread/1). If the calling thread is the‘pce thread’, [in_pce_thread_sync/1](mt-xpce.html#in_pce_thread_sync/1) executes a direct meta-call. See also [in_pce_thread/1](mt-xpce.html#in_pce_thread/1).

Note that [in_pce_thread_sync/1](mt-xpce.html#in_pce_thread_sync/1) is expensive because it requires copying and thread communication. For example, `in_pce_thread_synctrue` runs at approximately 50,000 calls per second (AMD Phenom 9600B, Ubuntu 11.04).

**pce_dispatch**(`+Options`)  
Create a Prolog thread with the alias name `pce` for XPCE event handling. In the X11 version this call creates a thread that executes the X11 event-dispatch loop. In MS-Windows it creates a thread that executes a windows event-dispatch loop. The XPCE event-handling thread has the alias `pce`. `Options` specifies the thread attributes as [thread_create/3](threadcreate.html#thread_create/3).
