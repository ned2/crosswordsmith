
## 4.43 Memory Management

### 4.43.1 Garbage collection

**garbage_collect**  
Invoke the global and trail stack garbage collector. Normally the garbage collector is invoked automatically if necessary. Explicit invocation might be useful to reduce the need for garbage collections in time-critical segments of the code. After the garbage collection [trim_stacks/0](memory.html#trim_stacks/0) is invoked to release the collected memory resources.

**garbage_collect_atoms**  
Reclaim unused atoms. Normally invoked after [agc_margin](flags.html#flag:agc_margin) (a Prolog flag) atoms have been created. On multithreaded versions the actual collection is delayed until there are no threads performing normal garbage collection. In this case [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0) returns immediately. Note that there is no guarantee it will *ever* happen, as there may always be threads performing garbage collection.

**garbage_collect_clauses**  
Reclaim retracted clauses. During normal operation, retracting a clause implies setting the *erased generation* to the current *generation* of the database and increment the generation. Keeping the clause around is both needed to realise the *logical update view* and deal with the fact that other threads may be executing the clause. Both static and dynamic code is processed this way.^(165Up to version 7.3.11, dynamic code was handled using *reference counts*.).

The clause garbage collector (CGC) scans the environment stacks of all threads for referenced dirty predicates and at which generation this reference accesses the predicate. It then removes the references for clauses that have been retracted before the oldest access generation from the clause list as well as the secondary clauses indexes of the predicate. If the clause list is not being scanned, the clause references and ultimately the clause itself is reclaimed.

The clause garbage collector is called under three conditions, (1) after *reloading* a source file, (2) if the memory occupied by retracted but not yet reclaimed clauses exceeds 12.5% of the program store, or (3) if skipping dead clauses in the clause lists becomes too costly. The cost of clause garbage collection is proportional with the total size of the local stack of all threads (the scanning phase) and the number of clauses in all‘dirty’predicates (the reclaiming phase).

**set_prolog_gc_thread**(`+Status`)  
Control whether or not atom and clause garbage collection are executed in a dedicated thread. The default is `true`. Values for `Status` are `true`, `false` and `stop`. The latter stops the `gc` thread but allows is to be recreated lazily. This is use by e.g., fork/1 to avoid forking a multi-threaded application. See also [gc_thread](flags.html#flag:gc_thread).

**trim_stacks**  
Release stack memory resources that are not in use at this moment, returning them to the operating system. It can be used to release memory resources in a backtracking loop, where the iterations require typically seconds of execution time and very different, potentially large, amounts of stack space. Such a loop can be written as follows:

``` code
loop :-
        generator,
            trim_stacks,
            potentially_expensive_operation,
        stop_condition, !.
```

The Prolog top-level loop is written this way, reclaiming memory resources after every user query. See also [trim_heap/0](memory.html#trim_heap/0) and [thread_idle/2](memory.html#thread_idle/2).

**set_prolog_stack**(`+Stack, +KeyValue`)  
Set a parameter for one of the Prolog runtime stacks. `Stack` is one of `local`, `global` or `trail`. The table below describes the `Key`(`Value`) pairs.

Current settings can be retrieved with [prolog_stack_property/2](memory.html#prolog_stack_property/2).

**min_free**(`+Cells`)  
Minimum amount of free space after trimming or shifting the stack. Setting this value higher can reduce the number of garbage collections and stack-shifts at the cost of higher memory usage. The amount is reported and specified in *cells*. A cell is 4 bytes in the 32-bit version and 8 bytes on the 64-bit version. See [address_bits](flags.html#flag:address_bits). See also [trim_stacks/0](memory.html#trim_stacks/0) and [debug/0](debugger.html#debug/0).

**low**(`+Cells`)  
**factor**(`+Number`)  
These two figures determine whether, if the stacks are low, a stack *shift* (expansion) or garbage collection is performed. This depends on these two parameters, the current stack usage and the amount of stack used after the last garbage collection. A garbage collection is started if `used > factor × lastused + low`.

**spare**(`+Cells`)  
All stacks trigger overflow before actually reaching the limit, so the resulting error can be handled gracefully. The spare stack is used for [print_message/2](printmsg.html#print_message/2) from the garbage collector and for handling exceptions. The default suffices, unless the user redefines related hooks. Do **not** specify large values for this because it reduces the amount of memory available for your real task.

Related hooks are [message_hook/3](printmsg.html#message_hook/3) (redefining GC messages), [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4) and prolog_exception_hook/5.

**prolog_stack_property**(`?Stack, ?KeyValue`)  
True if `KeyValue` is a current property of `Stack`. See [set_prolog_stack/2](memory.html#set_prolog_stack/2) for defined properties.

The total space limit for all stacks is controlled using the prolog flag [stack_limit](flags.html#flag:stack_limit).

### 4.43.2 Heap memory (malloc)

SWI-Prolog's memory management is based on the C runtime **malloc()** function and related functions. The characteristics of the **malloc()** implementation may affect performance and overall memory usage of the system. For most Prolog programs the performance impact of the allocator is small.^(166Multi-threaded applications may suffer from allocators that do not effectively avoid *false sharing* that affect CPU cache behaviour or operate using a single lock to provide thread safety. Such allocators should be rare in modern OSes.) The impact on total memory usage can be significant though, in particular for multi-threaded applications. This is due to two aspects of SWI-Prolog memory management:

- The Prolog stacks are allocated using **malloc()**. The stacks can be extremely large. SWI-Prolog assumes **malloc()** will use a mechanism that allows returning this memory to the OS. Most todays allocators satisfy this requirement.
- Atoms and clauses are allocated by the thread that requires them, but this memory is freed by the thread running the atom or clause garbage collector (see [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0) and [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0)). Normally these run in the thread `gc`, which means that all deallocation happens in this thread. Notably the [ptmalloc](http://www.malloc.de/en/) implementation used by the GNU C library (glibc) seems to handle this poorly.

Starting with version 8.1.27, SWI-Prolog by default links against [tcmalloc](https://github.com/google/tcmalloc) when available. Note that changing the allocator can only be done by linking the main executable (**swipl**) to an alternative library. When embedded (see [section 12.4.25](foreigninclude.html#sec:12.4.25)) the main program that embeds `libswipl` must be linked with tcmalloc. On ELF based systems (Linux), this effect can also be achieved using the environment variable `LD_PRELOAD`:

``` code
% LD_PRELOAD=/path/to/libtcmalloc.so swipl ...
```

SWI-Prolog attempts to detect the currently active allocator and sets the Prolog flag [malloc](flags.html#flag:malloc) if the detection succeeds. regardless of the malloc implementation, [trim_heap/0](memory.html#trim_heap/0) is provided.

\[det\]**trim_heap**  
his predicate attempts to return heap memory to the operating system. There is no portable way of doing so. If the system detects tcmalloc it calls **MallocExtension_ReleaseFreeMemory()**. If the system detects ptmalloc as provided by the GNU runtime library it calls **malloc_trim()**. In other cases this predicate simply succeeds. See also [trim_stacks/0](memory.html#trim_stacks/0)

#### 4.43.2.1 TCMalloc control predicates

If SWI-Prolog core detects that tcmalloc is the current allocator and provides the following additional predicates.

\[nondet\]**malloc_property**(`?Property`)  
True when `Property` is a property of the current allocator. The properties are defined by the allocator. The properties of tcmalloc are defined in `gperftools/malloc_extension.h`:^(167Documentation copied from the header.)

**’generic.current_allocated_bytes’**(`-Int`)  
Number of bytes currently allocated by application.

**’generic.heap_size’**(`-Int`)  
Number of bytes in the heap (= current_allocated_bytes + fragmentation + freed memory regions).

**’tcmalloc.max_total_thread_cache_bytes’**(`-Int`)  
Upper limit on total number of bytes stored across all thread caches.

**’tcmalloc.current_total_thread_cache_bytes’**(`-Int`)  
Number of bytes used across all thread caches.

**’tcmalloc.central_cache_free_bytes’**(`-Int`)  
Number of free bytes in the central cache that have been assigned to size classes. They always count towards virtual memory usage, and unless the underlying memory is swapped out by the OS, they also count towards physical memory usage.

**’tcmalloc.transfer_cache_free_bytes’**(`-Int`)  
Number of free bytes that are waiting to be transferred between the central cache and a thread cache. They always count towards virtual memory usage, and unless the underlying memory is swapped out by the OS, they also count towards physical

**’tcmalloc.thread_cache_free_bytes’**(`-Int`)  
Number of free bytes in thread caches. They always count towards virtual memory usage, and unless the underlying memory is swapped out by the OS, they also count towards physical memory usage.

**’tcmalloc.pageheap_free_bytes’**(`-Int`)  
Number of bytes in free, mapped pages in page heap. These bytes can be used to fulfill allocation requests. They always count towards virtual memory usage, and unless the underlying memory is swapped out by the OS, they also count towards physical memory usage. This property is not writable.

**’tcmalloc.pageheap_unmapped_bytes’**(`-Int`)  
Number of bytes in free, unmapped pages in page heap. These are bytes that have been released back to the OS, possibly by one of the MallocExtension "Release" calls. They can be used to fulfill allocation requests, but typically incur a page fault. They always count towards virtual memory usage, and depending on the OS, typically do not count towards physical memory usage.

\[det\]**set_malloc**(`+Property`)  
Set properties described in [malloc_property/1](memory.html#malloc_property/1). Currently the only writable property is `tcmalloc.max_total_thread_cache_bytes`. Setting an unknown property raises a `domain_error` and setting a read-only property raises a `permission_error` exception.

\[semidet\]**thread_idle**(`:Goal, +Duration`)  
Indicates to the system that the calling thread will idle for some time while calling `Goal` as [once/1](metacall.html#once/1). This call releases resources to the OS to minimise the footprint of the calling thread while it waits. Despite the name this predicate is always provided, also if the system is not configured with tcmalloc or is single threaded. `Duration` is one of

**short**  
Calls [trim_stacks/0](memory.html#trim_stacks/0) and, if tcmalloc is used, calls **MallocExtension_MarkThreadTemporarilyIdle()** which empties the thread's malloc cache but preserves the cache itself.

**long**  
Calls [garbage_collect/0](memory.html#garbage_collect/0) and [trim_stacks/0](memory.html#trim_stacks/0) and, if tcmalloc is used, calls **MallocExtension_MarkThreadIdle()** which releases all thread-specific allocation data structures.
