
## A.17 library(exceptions): Exception classification

Prolog [catch/3](exception.html#catch/3) selects errors based on unification. This is problematic for two reasons. First, one typically wants the exception term to be more specific than the term passed to the 2nd (`Ball`) argument of [catch/3](exception.html#catch/3). Second, in many situations one wishes to select multiple errors that may be raised by some operations, but let the others pass. Unification is often not suitable for this. For example, [open/3](IO.html#open/3) can raise an *existence_error* or a *permission_error* (and a couple more), but *existence_error* are also raised on, for example, undefined procedures. This is very hard to specify, Below is an attempt that still assumes nothing throws `error(_,_)`.

``` code
    catch(open(...), error(Formal,ImplDefined),
          (   ( Formal = existence_error(source_sink,_)
              ; Formal = permission_error(open, source_sink, _)
              )
          ->  <handle>
          ;   throw(Formal, ImplDefined)
          )),
    ...
```

Besides being hard to specify, actual Prolog systems define a large number of additional error terms because there is no reasonable ISO exception defined. For example, SWI-Prolog [open/3](IO.html#open/3) may raise `resource_error(max_files)` if the maximum number of file handles of the OS is exceeded.

As a result, we see a lot of Prolog code in the wild that simply uses the construct below to simply fail. But, this may fail for lack of stack space, a programmer error that causes a type error, etc. This both makes it much harder to debug the code and provide meaningful feedback to the user of the application.

``` code
    catch(Goal, _, fail)
```

Many programing languages have their exceptions organised by a (class) hierarchy. Prolog has no hierarchy of terms. We introduce [exception/2](exceptions.html#exception/2) as exception(+Type, ?Term), which can both be used as a type test for an exception term and as a *constraint* for the `Ball` of [catch/3](exception.html#catch/3). Using a predicate we can express abstractions over concrete exception terms with more flexibility than a hierarchy. Using a *multifile* predicate, libraries can add their exceptions to defined types or introduce new types.

The predicate [catch/4](exceptions.html#catch/4) completes the interface.

**catch**(`:Goal, +ExceptionType, ?Ball, :Recover`)  
As [catch/3](exception.html#catch/3), only catching exceptions for which `exception(ErrorType,Ball)` is true. See error/2. For example, the code below properly informs the user some file could not be processed due do some issue with `File`, while propagating on all other reasons while process/1 could not be executed.

``` code
    catch(process(File), file_error, Ball,
          file_not_processed(File, Ball))

file_not_processed(File, Ball) :-
    message_to_string(Ball, Msg),
    format(user_error, 'Could not process ~p: ~s', [File, Msg]).
```

\[det\]**exception**(`:Type, --Ball`)  
\[semidet\]**exception**(`:Type, +Ball`)  
If `Ball` is unbound, adds a delayed goal that tests the error belongs to `Type` when `Ball` is instantiated (by [catch/3](exception.html#catch/3)). Else succeed is error is of the specified `Type`.

Note that the delayed goal is added using [freeze/2](coroutining.html#freeze/2) and therefore the stepwise instantiation of `Ball` does not work, e.g. `exception(file_error, error(Formal,_))` immediately fails.

Error types may be defined or extended (e.g., by libraries) by adding clauses to the multifile predicates [error_term/2](exceptions.html#error_term/2) and [exception_term/2](exceptions.html#exception_term/2). *Modules* may (re-)define local error types using the [exception_type/2](exceptions.html#exception_type/2) directive.

\[nondet,multifile\]**error_term**(`?Type, ?Term`)  
Describe the formal part of `error(Formal,ImplDefined)` exceptions.

\[nondet,multifile\]**exception_term**(`?Type, ?Term`)  
Describe exceptions that are not `error(Formal, _)` terms.

**exception_type**(`+Type, +Term`)  
Declare all exceptions subsumed by `Term` to be an exception of `Type`. This declaration is module specific.
