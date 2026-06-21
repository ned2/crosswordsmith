
## A.62 library(thread_pool): Resource bounded thread management

See also  
http_handler/3 and http_spawn/2.

The module `library(thread_pool)` manages threads in pools. A pool defines properties of its member threads and the maximum number of threads that can coexist in the pool. The call [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4) allocates a thread in the pool, just like [thread_create/3](threadcreate.html#thread_create/3). If the pool is fully allocated it can be asked to wait or raise an error.

The library has been designed to deal with server applications that receive a variety of requests, such as HTTP servers. Simply starting a thread for each request is a bit too simple minded for such servers:

- Creating many CPU intensive threads often leads to a slow-down rather than a speedup.
- Creating many memory intensive threads may exhaust resources
- Tasks that require little CPU and memory but take long waiting for external resources can run many threads.

Using this library, one can define a pool for each set of tasks with comparable characteristics and create threads in this pool. Unlike the worker-pool model, threads are not started immediately. Depending on the design, both approaches can be attractive.

The library is implemented by means of a manager thread with the fixed thread id `__thread_pool_manager`. All state is maintained in this manager thread, which receives and processes requests to create and destroy pools, create threads in a pool and handle messages from terminated threads. Thread pools are *not* saved in a saved state and must therefore be recreated using the [initialization/1](consulting.html#initialization/1) directive or otherwise during startup of the application.

\[det\]**thread_pool_create**(`+Pool, +Size, +Options`)  
Create a pool of threads. A pool of threads is a declaration for creating threads with shared properties (stack sizes) and a limited number of threads. Threads are created using [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4). If all threads in the pool are in use, the behaviour depends on the `wait` option of [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4) and the `backlog` option described below. `Options` are passed to [thread_create/3](threadcreate.html#thread_create/3), except for

**backlog**(`+MaxBackLog`)  
Maximum number of requests that can be suspended. Default is `infinite`. Otherwise it must be a non-negative integer. Using `backlog(0)` will never delay thread creation for this pool.

The pooling mechanism does *not* interact with the `detached` state of a thread. Threads can be created both `detached` and normal and must be joined using [thread_join/2](threadcreate.html#thread_join/2) if they are not detached.

\[det\]**thread_pool_destroy**(`+Name`)  
Destroy the thread pool named `Name`.

Errors  
`existence_error(thread_pool, Name)`.

\[nondet\]**current_thread_pool**(`?Name`)  
True if `Name` refers to a defined thread pool.

\[nondet\]**thread_pool_property**(`?Name, ?Property`)  
True if `Property` is a property of thread pool `Name`. Defined properties are:

**options**(`Options`)  
Thread creation options for this pool

**free**(`Size`)  
Number of free slots on this pool

**size**(`Size`)  
Total number of slots on this pool

**members**(`ListOfIDs`)  
`ListOfIDs` is the list or threads running in this pool

**running**(`Running`)  
Number of running threads in this pool

**backlog**(`Size`)  
Number of delayed thread creations on this pool

\[det\]**thread_create_in_pool**(`+Pool, :Goal, -Id, +Options`)  
Create a thread in `Pool`. `Options` overrule default thread creation options associated to the pool. In addition, the following option is defined:

**wait**(`+Boolean`)  
If `true` (default) and the pool is full, wait until a member of the pool completes. If `false`, throw a resource_error.

Errors  
\- `resource_error(threads_in_pool(Pool))` is raised if wait is `false` or the backlog limit has been reached.  
- `existence_error(thread_pool, Pool)` if `Pool` does not exist.

**worker_exitted**(`+PoolName, +WorkerId, :AtExit`)  
It is possible that’\_\_thread_pool_manager’no longer exists while closing down the process because the manager was killed before the worker.

To be done  
Find a way to discover that we are terminating Prolog.

\[semidet,multifile\]**create_pool**(`+PoolName`)  
Hook to create a thread pool lazily. The hook is called if [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4) discovers that the thread pool does not exist. If the hook succeeds, [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4) retries creating the thread. For example, we can use the following declaration to create threads in the pool `media`, which holds a maximum of 20 threads.

``` code
:- multifile thread_pool:create_pool/1.

thread_pool:create_pool(media) :-
    thread_pool_create(media, 20, []).
```
