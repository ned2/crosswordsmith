
## 10.4 Thread synchronisation

All internal Prolog operations are thread-safe. This implies that two Prolog threads can operate on the same dynamic predicate without corrupting the consistency of the predicate. This section deals with user-level *mutexes* (called *monitors* in ADA or *critical sections* by Microsoft). A mutex is a **MUT**ual **EX**clusive device, which implies that at most one thread can *hold* a mutex.

Mutexes are used to realise related updates to the Prolog database. With‘related’, we refer to the situation where a‘transaction’implies two or more changes to the Prolog database. For example, we have a predicate address/2 , representing the address of a person and we want to change the address by retracting the old and asserting the new address. Between these two operations the database is invalid: this person has either no address or two addresses, depending on the assert/retract order.

The code below provides a solution to this problem based on [with_mutex/2](threadsync.html#with_mutex/2). It also illustrates the problem of mutexes. The predicate [with_mutex/2](threadsync.html#with_mutex/2) behaves as [once/1](metacall.html#once/1) with respect to the guarded goal. This means that our predicate address/2 is no longer a nice logical non-deterministic relation. This could be solved by explicit locking and unlocking a mutex using [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3), but at the risk of deadlocking the program if the choice point is left open by accident.

``` code
change_address(Id, Address) :-
        with_mutex(addressbook,
                   ( retractall(address_db(Id, _)),
                     asserta(address_db(Id, Address))
                   )).

address(Id, Address) :-
        with_mutex(addressbook,
                   address_db(Id, Address)).
```

Message queues (see [message_queue_create/2](threadcom.html#message_queue_create/2)) often provide simpler and more robust ways for threads to communicate. Still, mutexes can be a sensible solution and are therefore provided.

**mutex_create**(`?MutexId`)  
Create a mutex. If `MutexId` is an atom, a *named* mutex is created. If it is a variable, an anonymous mutex reference is returned. Anonymous mutexes are subject to (atom) garbage collection.

**mutex_create**(`-MutexId, +Options`)  
Create a mutex using options. Defined options are:

**alias**(`Alias`)  
Set the alias name. Using `mutex_create(X, [alias(name)])` is preferred over the equivalent `mutex_create(name)`.

**mutex_destroy**(`+MutexId`)  
Destroy a mutex. If the mutex is not locked, it is destroyed and further access yields an `existence_error` exception. As of version 7.1.19, this behaviour is reliable. If the mutex is locked, the mutex is scheduled for *delayed destruction*: it will be destroyed when it becomes unlocked.

**with_mutex**(`+MutexId, :Goal`)  
Execute `Goal` while holding `MutexId`. If `Goal` leaves choice points, these are destroyed (as in [once/1](metacall.html#once/1)). The mutex is unlocked regardless of whether `Goal` succeeds, fails or raises an exception. An exception thrown by `Goal` is re-thrown after the mutex has been successfully unlocked. See also [mutex_create/1](threadsync.html#mutex_create/1) and [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3).

Although described in the thread section, this predicate is also available in the single-threaded version, where it behaves simply as [once/1](metacall.html#once/1).

**mutex_lock**(`+MutexId`)  
Lock the mutex. Prolog mutexes are *recursive* mutexes: they can be locked multiple times by the same thread. Only after unlocking it as many times as it is locked does the mutex become available for locking by other threads. If another thread has locked the mutex the calling thread is suspended until the mutex is unlocked.

If `MutexId` is an atom, and there is no current mutex with that name, the mutex is created automatically using [mutex_create/1](threadsync.html#mutex_create/1). This implies named mutexes need not be declared explicitly.

Please note that locking and unlocking mutexes should be paired carefully. Especially make sure to unlock mutexes even if the protected code fails or raises an exception. For most common cases, use [with_mutex/2](threadsync.html#with_mutex/2), which provides a safer way for handling Prolog-level mutexes. The predicate [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) is another way to guarantee that the mutex is unlocked while retaining non-determinism.

**mutex_trylock**(`+MutexId`)  
As [mutex_lock/1](threadsync.html#mutex_lock/1), but if the mutex is held by another thread, this predicates fails immediately.

**mutex_unlock**(`+MutexId`)  
Unlock the mutex. This can only be called if the mutex is held by the calling thread. If this is not the case, a `permission_error` exception is raised.

\[deprecated\]**mutex_unlock_all**  
Unlock all mutexes held by the current thread. This predicate should not be needed if mutex unlocking is guaranteed with [with_mutex/2](threadsync.html#with_mutex/2) or [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3).^(209The also deprecated [thread_exit/1](threadcreate.html#thread_exit/1) bypasses the automatic cleanup.)

**mutex_property**(`?MutexId, ?Property`)  
True if `Property` is a property of `MutexId`. Defined properties are:

**alias**(`Alias`)  
Mutex has the defined alias name. See [mutex_create/2](threadsync.html#mutex_create/2) using the‘alias’option.

**status**(`Status`)  
Current status of the mutex. One of `unlocked` if the mutex is currently not locked, or `locked(Owner, Count)` if mutex is locked `Count` times by thread `Owner`. Note that unless `Owner` is the calling thread, the locked status can change at any time. There is no useful application of this property, except for diagnostic purposes.^(bugAs `Owner` and `Count` are fetched separately from the mutex, the values may be inconsistent.)
