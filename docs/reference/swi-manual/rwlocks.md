
## A.52 library(rwlocks): Read/write locks

This library implements *read/write* locks on top of [with_mutex/2](threadsync.html#with_mutex/2). *Read/write* locks are synchronization objects that allow for multiple readers or a single writer to be active.

**with_rwlock**(`+LockId, :Goal, +ModeSpec`)  
**with_rwlock**(`+LockId, :Goal, +ModeSpec, +Options`)  
Run `Goal`, synchronized with `LockId` in `ModeSpec`. `ModeSpec` is one of `read`, `write`, `read(Priority)` or `write(Priority)`. The default `read` priority is 100 and the default `write` priority is 200. These values prioritize writers over readers. `Goal` may start if

- If there is no goal waiting with higher priority **and**
  - It is a read goal and no write goal is running **or**
  - It is a write goal and no other goal is running.

If `Goal` may not start immediately the thread waits using [thread_wait/2](threadcom.html#thread_wait/2). The `Options` `timeout` and `deadline` are passed to [thread_wait/2](threadcom.html#thread_wait/2). If the time limit is exceeded an exception is raised.

*Read/write* locks are widely critized for their poor behaviour on several workloads. They perform well in scenarios where read operations take long, and write operations are relatively fast and occur only occasionally. *Transactions*, as implemented by [transaction/1](db.html#transaction/1),2 are often a better alternative.

This predicate uses a normal mutex and a flag with the same name. See [with_mutex/2](threadsync.html#with_mutex/2) and [flag/3](db.html#flag/3). Neither the mutex nor the flag should be used directly.

throws  
`time_limit_exceeded(rwlock)` if a timeout or deadline is specified and this is exceeded.

bug  
The current implementation is written in Prolog and comes with significant overhead. It is intended to synchronize slow operations.
