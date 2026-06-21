
## 10.6 Multithreaded mixed C and Prolog applications

All foreign code linked to the multithreading version of SWI-Prolog should be thread-safe (*reentrant*) or guarded in Prolog using [with_mutex/2](threadsync.html#with_mutex/2) from simultaneous access from multiple Prolog threads. If you want to write mixed multithreaded C and Prolog applications you should first familiarise yourself with writing multithreaded applications in C (C++).

If you are using SWI-Prolog as an embedded engine in a multithreaded application you can access the Prolog engine from multiple threads by creating an *engine* in each thread from which you call Prolog. Without creating an engine, a thread can only use functions that do *not* use the `term_t` type (for example [PL_new_atom()](foreigninclude.html#PL_new_atom())).

The system supports two models. [Section 10.6.1](foreignthread.html#sec:10.6.1) describes the original one-to-one mapping. In this schema a native thread attaches a Prolog thread if it needs to call Prolog and detaches it when finished, as opposed to the model from [section 10.6.2](foreignthread.html#sec:10.6.2), where threads temporarily use a Prolog engine.

### 10.6.1 A Prolog thread for each native thread (one-to-one)

In the one-to-one model, the thread that called [PL_initialise()](foreigninclude.html#PL_initialise()) has a Prolog engine attached. If another C thread in the system wishes to call Prolog it must first attach an engine using [PL_thread_attach_engine()](foreignthread.html#PL_thread_attach_engine()) and call [PL_thread_destroy_engine()](foreignthread.html#PL_thread_destroy_engine()) after all Prolog work is finished. This model is especially suitable with long running threads that need to do Prolog work regularly. See [section 10.6.2](foreignthread.html#sec:10.6.2) for the alternative many-to-many model.

`int` **PL_thread_self**()  
Returns the integer Prolog identifier of the engine or -1 if the calling thread has no Prolog engine. This function is also provided in the single-threaded version of SWI-Prolog, where it returns -2.

`int` **PL_unify_thread_id**(`term_t t, int i`)  
Unify `t` with the Prolog thread identifier for thread `i`. Thread identifiers are normally returned from [PL_thread_self()](foreignthread.html#PL_thread_self()). Returns -1 if the thread does not exist or the unification fails.

`int` **PL_thread_attach_engine**(`const PL_thread_attr_t *attr`)  
Creates a new Prolog engine in the calling thread. If the calling thread already has an engine the reference count of the engine is incremented. The `attr` argument can be `NULL` to create a thread with default attributes. Otherwise it is a pointer to a structure as defined below. The structure must be fully initialized, *including* hidden fields. For any field with value‘0’, the default is used. The `cancel` field may be filled with a pointer to a function that is called when [PL_cleanup()](foreigninclude.html#PL_cleanup()) terminates the running Prolog engines. The function is called with the thread id (see [PL_thread_self()](foreignthread.html#PL_thread_self())) as argument and must return `PL_THREAD_CANCEL_JOINED` if the thread was reclaimed successfully, `PL_THREAD_CANCEL_MUST_JOIN` if the thread as cancelled, but must still be joined or `PL_THREAD_CANCEL_FAILED` if the request cannot be honoured. If this function is not present or returns `PL_THREAD_CANCEL_FAILED` **pthread_cancel()** is used. The new thread inherits is properties from Prolog's `main` thread. The `flags` field defines the following flags:

**PL_THREAD_NO_DEBUG**  
If this flag is present, the thread starts in normal no-debug status. By default, the debug status is inherited from the main thread.

**PL_THREAD_NOT_DETACHED**  
By default the new thread is created in *detached* mode. With this flag it is created normally, allowing Prolog to *join* the thread.

**PL_THREAD_CUR_STREAMS**  
By default the `current_input` and `current_output` are set to `user_input` and `user_output` of the main thread. Using this flag, these streams are copied from the main thread. See also the `inherited_from` option of [thread_create/3](threadcreate.html#thread_create/3).

``` code
typedef struct
{ size_t    stack_limit;                /* Total stack limit (bytes) */
  size_t    table_space;                /* Total tabling space limit (bytes) */
  char *    alias;                      /* alias name */
  int       (*cancel)(int thread);      /* cancel function */
  intptr_t  flags;                      /* PL_THREAD_* flags */
  size_t    max_queue_size;             /* Max size of associated queue */
  char *    thread_class;               /* Class property of the thread */
} PL_thread_attr_t;
```

The structure may be destroyed after [PL_thread_attach_engine()](foreignthread.html#PL_thread_attach_engine()) has returned. On success it returns the Prolog identifier for the thread (as returned by [PL_thread_self()](foreignthread.html#PL_thread_self())). If an error occurs, -1 is returned. If this Prolog is not compiled for multithreading, -2 is returned.

`bool` **PL_thread_destroy_engine**()  
Destroy the Prolog engine in the calling thread. Only takes effect if [PL_thread_destroy_engine()](foreignthread.html#PL_thread_destroy_engine()) is called as many times as [PL_thread_attach_engine()](foreignthread.html#PL_thread_attach_engine()) in this thread. Returns `TRUE` on success and `FALSE` if the calling thread has no engine or this Prolog does not support threads.

Please note that construction and destruction of engines are relatively expensive operations. Only destroy an engine if performance is not critical and memory is a critical resource.

`bool` **PL_thread_at_exit**(`void (*function)(void *), void *closure, bool global`)  
Register a handle to be called as the Prolog engine is destroyed. The handler function is called with one `void *` argument holding `closure`. If `global` is `true`, the handler is installed *for all threads*. Globally installed handlers are executed after the thread-local handlers. If the handler is installed local for the current thread only (`global` == `false`) it is stored in the same FIFO queue as used by [thread_at_exit/1](threadcreate.html#thread_at_exit/1).

### 10.6.2 Using Prolog engines from C

Prolog engines live as entities that are independent from threads. They are always supported in the multi-threaded version and may be enabled in the single threaded version by providing `-DENGINES=ON` during the **cmake** configuration. Multiple threads may use a pool of engines for handling calls to Prolog. A single thread may use multiple engines to achieve *coroutining*. Engines are suitable in the following identified cases:

- *Many native threads with infrequent Prolog work*  
  Prolog threads are expensive in terms of memory and time to create and destroy them. For systems that use a large number of threads that only infrequently need to call Prolog, it is better to take an engine from a pool and return it there.
- *Prolog status must be handed to another thread*  
  This situation has been identified by Uwe Lesta when creating a .NET interface for SWI-Prolog. .NET distributes work for an active internet connection over a pool of threads. If a Prolog engine contains the state for a connection, it must be possible to detach the engine from a thread and re-attach it to another thread handling the same connection.
- *Achieving coroutines*  
  A single thread may use engines to implement *coroutining*. This is notably interesting when combined with *yielding* as described in [section 12.4.1.2](foreigninclude.html#sec:12.4.1.2).

`PL_engine_t` **PL_current_engine**()  
Returns the current engine of the calling thread or `NULL` if the thread has no Prolog engine.

`PL_engine_t` **PL_create_engine**(`PL_thread_attr_t *attributes`)  
Create a new Prolog engine. `attributes` is described with [PL_thread_attach_engine()](foreignthread.html#PL_thread_attach_engine()). Any thread can make this call after [PL_initialise()](foreigninclude.html#PL_initialise()) returns success. The returned engine is not attached to any thread and lives until [PL_destroy_engine()](foreignthread.html#PL_destroy_engine()) is used on the returned handle.

In the single-threaded version this call always returns `NULL`, indicating failure.

`bool` **PL_destroy_engine**(`PL_engine_t e`)  
Destroy the given engine. Destroying an engine is only allowed if the engine is not attached to any thread or attached to the calling thread. On success this function returns `TRUE`, on failure the return value is `FALSE`.

`int` **PL_set_engine**(`PL_engine_t engine, PL_engine_t *old`)  
Make the calling thread ready to use `engine`. If `old` is non-`NULL` the current engine associated with the calling thread is stored at the given location. If `engine` equals `PL_ENGINE_MAIN` the initial engine is attached to the calling thread. If `engine` is `PL_ENGINE_CURRENT` the engine is not changed. This can be used to query the current engine. This call returns `PL_ENGINE_SET` if the engine was switched successfully, `PL_ENGINE_INVAL` if `engine` is not a valid engine handle and `PL_ENGINE_INUSE` if the engine is currently in use by another thread.

Engines can be changed at any time. For example, it is allowed to select an engine to initiate a Prolog goal, detach it and at a later moment execute the goal from another thread. Note, however, that the `term_t`, `qid_t` and `fid_t` types are interpreted relative to the engine for which they are created. Behaviour when passing one of these types from one engine to another is undefined. The engine to which a query belongs can be requested using [PL_query_engine()](foreigninclude.html#PL_query_engine())

In versions that do not support engines this call only succeeds if `engine` refers to the main engine.

`void` **PL_WITH_ENGINE**(`PL_engine_t e`)  
This macro implements a C *for-loop* where the body is executed with the engine `e` as current engine. The body is executed exactly once. After completion of the body the current engine of the calling thread is restored to old state (either the old current engine or no engine). The user may use `break` to terminate the body early. The user *may not* use `return`. Using `return` does not reset the old engine.
