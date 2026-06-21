
## B.7 Hooks using the exception predicate

This section describes the predicate [exception/3](exception3.html#exception/3), which can be defined by the user in the module `user` as a multifile predicate. Unlike the name suggests, this is actually a *hook* predicate that has no relation to Prolog exceptions as defined by the ISO predicates [catch/3](exception.html#catch/3) and [throw/1](exception.html#throw/1).

The predicate [exception/3](exception3.html#exception/3) is called by the kernel on a couple of events, allowing the user to‘fix’errors just-in-time. The mechanism allows for *lazy* creation of objects such as predicates.

**exception**(`+Exception, +Context, -Action`)  
Dynamic predicate, normally not defined. Called by the Prolog system on run-time exceptions that can be repaired‘just-in-time’. The values for `Exception` are described below. See also [catch/3](exception.html#catch/3) and [throw/1](exception.html#throw/1).

If this hook predicate succeeds it must instantiate the `Action` argument to the atom `fail` to make the operation fail silently, `retry` to tell Prolog to retry the operation or `error` to make the system generate an exception. The action `retry` only makes sense if this hook modified the environment such that the operation can now succeed without error.

**undefined_predicate**  
`Context` is instantiated to a predicate indicator (\[module\]:\<`name`\>/\<`arity`\>). If the predicate fails, Prolog will generate an `existence_error` exception. The hook is intended to implement alternatives to the built-in autoloader, such as autoloading code from a database. Do *not* use this hook to suppress existence errors on predicates. See also [unknown](flags.html#flag:unknown) and [section 2.14](autoload.html#sec:2.14).

**undefined_global_variable**  
`Context` is instantiated to the name of the missing global variable. The hook must call [nb_setval/2](gvar.html#nb_setval/2) or [b_setval/2](gvar.html#b_setval/2) before returning with the action `retry`. See also [nb_current/2](gvar.html#nb_current/2).
