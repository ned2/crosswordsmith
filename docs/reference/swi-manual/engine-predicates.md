
## 11.3 Engine predicate reference

This section documents the built-in predicates that deal with engines. In addition to these, most predicates dealing with threads and message queue can be used to access engines.

\[det\]**engine_create**(`+Template, :Goal, ?Engine`)  
\[det\]**engine_create**(`+Template, :Goal, -Engine, +Options`)  
Create a new engine and unify `Engine` with a handle to it. `Template` and `Goal` form a pair similar to [findall/3](allsolutions.html#findall/3): the instantiation of `Template` becomes available through [engine_next/2](engine-predicates.html#engine_next/2) after `Goal` succeeds. `Options` is a list of the following options. See [thread_create/3](threadcreate.html#thread_create/3) for details.

**alias**(`+Name`)  
Give the engine a name. `Name` must be an atom. If this option is provided, `Engine` is unified with `Name`. The name space for engines is shared with threads and mutexes.

**stack**(`+Bytes`)  
Set the stack limit for the engine. The default is inherited from the calling thread.

The `Engine` argument of [engine_create/3](engine-predicates.html#engine_create/3) may be instantiated to an atom, creating an engine with the given alias.

\[det\]**engine_destroy**(`+Engine`)  
Destroy `Engine`.

\[semidet\]**engine_next**(`+Engine, -Term`)  
Ask the engine `Engine` to produce a next answer. On this first call on a specific engine, the `Goal` of the engine is started. If a previous call returned an answer through completion, this causes the engine to backtrack and finally, if the engine produces a previous result using [engine_yield/1](engine-predicates.html#engine_yield/1), execution proceeds after the [engine_yield/1](engine-predicates.html#engine_yield/1) call.

\[det\]**engine_next_reified**(`+Engine, -Term`)  
Similar to [engine_next/2](engine-predicates.html#engine_next/2), but instead of success, failure or raising an exception, `Term` is unified with one of terms below. This predicate is provided primarily for compatibility with Lean Prolog.

**the**(`Answer`)  
Goal succeeded with `Template` bound to `Answer` or Goal yielded with a term `Answer`.

**no**  
Goal failed.

**throw**(`Exception`)  
Goal raised `Exception`.

\[det\]**engine_post**(`+Engine, +Term`)  
Make `Term` available to [engine_fetch/1](engine-predicates.html#engine_fetch/1) inside the `Engine`. This call must be followed by a call to [engine_next/2](engine-predicates.html#engine_next/2) and the engine must call [engine_fetch/1](engine-predicates.html#engine_fetch/1).

\[det\]**engine_post**(`+Engine, +Term, -Reply`)  
Combines [engine_post/2](engine-predicates.html#engine_post/2) and [engine_next/2](engine-predicates.html#engine_next/2).

\[det\]**engine_yield**(`+Term`)  
Called from within the engine, causing [engine_next/2](engine-predicates.html#engine_next/2) in the caller to return with `Term`. A subsequent call to [engine_next/2](engine-predicates.html#engine_next/2) causes [engine_yield/1](engine-predicates.html#engine_yield/1) to‘return’. This predicate can only be called if the engine is not involved in a callback from C, i.e., when the engine calls a predicate defined in C that calls back Prolog it is not possible to use this predicate. Trying to do so results in a `permission_error` exception.

\[det\]**engine_fetch**(`-Term`)  
Called from within the engine to fetch the term made available through [engine_post/2](engine-predicates.html#engine_post/2) or [engine_post/3](engine-predicates.html#engine_post/3). If no term is available an existence_error exception is raised.

\[det\]**engine_self**(`-Engine`)  
Called from within the engine to get access to the handle to the engine itself.

\[semidet\]**is_engine**(`@Term`)  
True if `Term` is a reference to or the alias name of an existing engine.

\[nondet\]**current_engine**(`-Engine`)  
True when `Engine` is an existing engine.
