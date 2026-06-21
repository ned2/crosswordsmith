
## 4.33 Global variables

Global variables are associations between names (atoms) and terms. They differ in various ways from storing information using [assert/1](db.html#assert/1) or [recorda/3](db.html#recorda/3).

- The value lives on the Prolog (global) stack. This implies that lookup time is independent of the size of the term. This is particularly interesting for large data structures such as parsed XML documents or the CHR global constraint store.
- They support both global assignment using [nb_setval/2](gvar.html#nb_setval/2) and backtrackable assignment using [b_setval/2](gvar.html#b_setval/2).
- Only one value (which can be an arbitrary complex Prolog term) can be associated to a variable at a time.
- Their value cannot be shared among threads. Each thread has its own namespace and values for global variables.
- Currently global variables are scoped globally. We may consider module scoping in future versions.

Both [b_setval/2](gvar.html#b_setval/2) and [nb_setval/2](gvar.html#nb_setval/2) implicitly create a variable if the referenced name does not already refer to a variable.

Global variables may be initialised from directives to make them available during the program lifetime, but some considerations are necessary for saved states and threads. Saved states do not store global variables, which implies they have to be declared with [initialization/1](consulting.html#initialization/1) to recreate them after loading the saved state. Each thread has its own set of global variables, starting with an empty set. Using [thread_initialization/1](threadcreate.html#thread_initialization/1) to define a global variable it will be defined, restored after reloading a saved state and created in all threads that are created *after* the registration. Finally, global variables can be initialised using the exception hook [exception/3](exception3.html#exception/3). See also [nb_current/2](gvar.html#nb_current/2) and [nb_delete/1](gvar.html#nb_delete/1).

**b_setval**(`+Name, +Value`)  
Associate the term `Value` with the atom `Name` or replace the currently associated value with `Value`. On backtracking the assignment is reversed. If the variable `Name` did not exist before calling [b_setval/2](gvar.html#b_setval/2), backtracking causes the variable to be deleted.^(147Prior to version 8.3.28 backtracking over the variable creation caused the variable to get the value `[]`, i.e., the empty list. If this is desirable use `nb_setval(Var, [])` before [b_setval/2](gvar.html#b_setval/2).)

**b_getval**(`+Name, -Value`)  
Get the value associated with the global variable `Name` and unify it with `Value`. Note that this unification may further instantiate the value of the global variable. If this is undesirable the normal precautions (double negation or [copy_term/2](manipterm.html#copy_term/2)) must be taken. The [b_getval/2](gvar.html#b_getval/2) predicate generates errors if `Name` is not an atom or the requested variable does not exist.

&nbsp;

**nb_setval**(`+Name, +Value`)  
Associates a copy of `Value` created with [duplicate_term/2](manipterm.html#duplicate_term/2) with the atom `Name`. Note that this can be used to set an initial value other than `[]` prior to backtrackable assignment. Starting with version 9.3.18, if the new `Value` contains the old value, the old value is **not** copied. This implies that push/1 below has complexity O(1), regardless of the length of the list `Old`.

``` code
push(Var, Value) :-
    nb_getval(Var, Old),
    nb_setval(Var, [Value|Old]).
```

**nb_getval**(`+Name, -Value`)  
The [nb_getval/2](gvar.html#nb_getval/2) predicate is a synonym for [b_getval/2](gvar.html#b_getval/2), introduced for compatibility and symmetry. As most scenarios will use a particular global variable using either non-backtrackable or backtrackable assignment, using [nb_getval/2](gvar.html#nb_getval/2) can be used to document that the variable is non-backtrackable. Raises `existence_error(variable, Name)` if the variable does not exist. Alternatively, [nb_current/2](gvar.html#nb_current/2) can used to query a global variable. This version *fails* if the variable does not exist rather than raising an exception.

**nb_linkval**(`+Name, +Value`)  
Associates the term `Value` with the atom `Name` without copying it. This is a fast special-purpose variation of [nb_setval/2](gvar.html#nb_setval/2) intended for expert users only because the semantics on backtracking to a point before creating the link are poorly defined for compound terms. The principal term is always left untouched, but backtracking behaviour on arguments is undone if the original assignment was *trailed* and left alone otherwise, which implies that the history that created the term affects the behaviour on backtracking. Consider the following example:

``` code
demo_nb_linkval :-
        T = nice(N),
        (   N = world,
            nb_linkval(myvar, T),
            fail
        ;   nb_getval(myvar, V),
            writeln(V)
        ).
```

**nb_current**(`?Name, ?Value`)  
Enumerate all defined variables with their value. The order of enumeration is undefined. Note that [nb_current/2](gvar.html#nb_current/2) can be used as an alternative for [nb_getval/2](gvar.html#nb_getval/2) to request the value of a variable and fail silently if the variable does not exists. Note that if the variable is not defined, [exception/3](exception3.html#exception/3) is called attempting to define it. As of version 8.3.28, a failure of [exception/3](exception3.html#exception/3) to define the variable causes the variable to be defined with a reserved valued to avoid subsequent calls to [exception/3](exception3.html#exception/3).

**nb_delete**(`+Name`)  
Delete the named global variable. Succeeds also if the named variable does not exist. Deleting a global variable ensures the variable is associated to a reserved value to avoid subsequent calls to [exception/3](exception3.html#exception/3). Note that this implies that the resources associated with a global variable are never fully reclaimed.

### 4.33.1 Compatibility of SWI-Prolog Global Variables

Global variables have been introduced by various Prolog implementations recently. The implementation of them in SWI-Prolog is based on hProlog by Bart Demoen. In discussion with Bart it was decided that the semantics of hProlog [nb_setval/2](gvar.html#nb_setval/2), which is equivalent to [nb_linkval/2](gvar.html#nb_linkval/2), is not acceptable for normal Prolog users as the behaviour is influenced by how built-in predicates that construct terms ([read/1](termrw.html#read/1), =../2, etc.) are implemented.

GNU-Prolog provides a rich set of global variables, including arrays. Arrays can be implemented easily in SWI-Prolog using [functor/3](manipterm.html#functor/3) and [setarg/3](manipterm.html#setarg/3) due to the unrestricted arity of compound terms.
