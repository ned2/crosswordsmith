
## B.8 Prolog events

Version 8.1.9 introduces a uniform mechanism to listen to events that happen in the Prolog engine. It replaces and generalises prolog_event_hook/1 , a hook that was introduced to support the graphical debugger. The current implementation deals with debug, thread and dynamic database events. We expect this mechanism to deal with more hooks in the future.

**prolog_listen**(`+Channel, :Closure`)  
**prolog_listen**(`+Channel, :Closure, +Options`)  
Call `Closure` if an event that matches `Channel` happens inside Prolog. Possible choice points are pruned as by [once/1](metacall.html#once/1). Possible failure is ignored, but exceptions are propagated into the environment. Multiple closures can be associated with the same channel. Execution of the list of closures may be terminated by an exception. Options:

**as**(`+Location`)  
`Location` is one of `first` (default) or `last` and determines whether the new handler is expected as first or last.

**name**(`+Atom`)  
Give the handler a name. A new registration using the same name replaces the existing handler rather than adding a new handler. Names are local to the `Channel`, i.e., different channels can use the same name.

Defined channels are described below. The `Channel` argument is the name of the term listed below. The arguments are added as additional arguments to the given `Closure`.

**abort**  
Called by [abort/0](toplevel.html#abort/0).

**erase**(`DbRef`)  
Called on an erased recorded database reference or clause. Note that a retracted clauses is not immediately removed. Clauses are reclaimed by [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0), which is normally executed automatically in the `gc` thread. This specific channel is used by clause_info/5 to reclaim source layout of reclaimed clauses. User applications should typically use the `PredicateIndicator` channel.

**break**(`Action, ClauseRef, PCOffset`)  
Traps events related to Prolog break points. See library `library(prolog_breakpoints)`

**frame_finished**(`FrameRef`)  
Indicates that a stack frame that has been examined using [prolog_current_frame/1](manipstack.html#prolog_current_frame/1), [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3) and friends has been deleted. Used by the source level debugger to avoid that the stack view references non-existing frames.

**thread_exit**(`Thread`)  
Globally registered channel that is called by any thread just before the thread is terminated.

**thread_start**(`Thread`)  
Globally registered channel that is called by any thread after the thread initialization and before running the thread's goal.

**this_thread_exit**  
Thread local version of the `thread_exit` channel that is also used by the `at_exit(Closure)` option of [thread_create/3](threadcreate.html#thread_create/3).

**PredicateIndicator**(`Action, Context`)  
Track changes to a predicate. This notably allows tracking modifications to dynamic predicates. The channel also allows tracking changes to *monotonic* tables ([section 7.8](tabling-monotonic.html#sec:7.8)). Both monotonic and incremental tabling use this to track changes to `incremental` and `monotonic` dynamic predicates. Below is an example illustrating events from changing a dynamic predicate.

``` code
:- dynamic p/1.
:- prolog_listen(p/1, updated(p/1)).

updated(Pred, Action, Context) :-
    format('Updated ~p: ~p ~p~n', [Pred, Action, Context]).
```

``` code
?- assert(p(a)).
Updated p/1: assertz <clause>(0x55db261709d0)
?- retractall(p(_)).
Updated p/1: retractall start(user:p(_12294))
Updated p/1: retract <clause>(0x55db261719c0)
Updated p/1: retractall end(user:p(_12294))
```

**asserta**  
**assertz**  
A new clauses has been added as first (last) for the given predicate. `Context` is a clause reference. The hook is called after the clause has been added. If the hook fails the clause is removed.

**retract**  
A clause was retracted from the given predicate using either [retract/1](db.html#retract/1), [erase/1](db.html#erase/1) or [retractall/1](db.html#retractall/1). `Context` is a clause reference. The hook is called before the clause is removed. If the hook fails, the clause is not removed.

**retractall**  
The beginning and end of [retractall/1](db.html#retractall/1) is indicated with the `Action` `retractall`. The context argument is `start(Head)` or `end(Head)`.

**rollback**(`Action`)  
Issued when rolling back (discarding) a transaction. `Action` is the local action being reverted and is one of `asserta`, `assertz` or `retract`. Context is the involved clause. See [transaction/1](db.html#transaction/1) and [snapshot/1](db.html#snapshot/1).

**new_answer**  
A new answer was added to a tabled predicate. The context is the answer term. Currently implemented for *monotonic* tabling only. Future versions may also implement this for normal tabling. See [section 7.8.2](tabling-monotonic.html#sec:7.8.2).

**prolog_unlisten**(`+Channel, :Closure`)  
Remove matching closures registered with [prolog_listen/3](prolog-event.html#prolog_listen/3).
