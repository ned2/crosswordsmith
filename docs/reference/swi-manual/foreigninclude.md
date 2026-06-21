
## 12.4 The Foreign Include File

### 12.4.1 Argument Passing and Control

If Prolog encounters a foreign predicate at run time it will call a function specified in the predicate definition of the foreign predicate. The arguments `1, ... , <``arity``>` pass the Prolog arguments to the goal as Prolog terms. Foreign functions should be declared of type `foreign_t`.

All the arguments to a foreign predicate must be of type `term_t`. The only operation that is allowed with an argument to a foreign predicate is unification; for anything that might over-write the term, you must use a copy created by [PL_copy_term_ref()](foreigntypes.html#PL_copy_term_ref()). For an example, see [PL_unify_list()](foreigninclude.html#PL_unify_list()).

Deterministic foreign functions return with either `TRUE` (success) or `FALSE` (failure).^(215`SWI-Prolog.h` defines the macros `PL_succeed` and `PL_fail` to return with success or failure. These macros should be considered deprecated.) The foreign function may raise an exception using [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) or one of the shorthands for commonly used exceptions such as [PL_type_error()](foreigninclude.html#PL_type_error()). Note that the C language does not provide exception handling and therefore the functions that raise an exception return (with the value `FALSE`). Functions that raise an exception *must* return `FALSE`.

#### 12.4.1.1 Non-deterministic Foreign Predicates

By default foreign predicates are deterministic. Using the `PL_FA_NONDETERMINISTIC` attribute (see [PL_register_foreign()](foreigninclude.html#PL_register_foreign())) it is possible to register a predicate as a non-deterministic predicate. Writing non-deterministic foreign predicates is slightly more complicated as the foreign function needs context information for generating the next solution. Note that the same foreign function should be prepared to be simultaneously active in more than one goal. Suppose the natural_number_below_n/2 is a non-deterministic foreign predicate, backtracking over all natural numbers lower than the first argument. Now consider the following predicate:

``` code
quotient_below_n(Q, N) :-
        natural_number_below_n(N, N1),
        natural_number_below_n(N, N2),
        Q =:= N1 / N2, !.
```

In this predicate the function natural_number_below_n/2 simultaneously generates solutions for both its invocations.

Non-deterministic foreign functions should be prepared to handle three different calls from Prolog:

- *Initial call (`PL_FIRST_CALL`)*  
  Prolog has just created a frame for the foreign function and asks it to produce the first answer.
- *Redo call (`PL_REDO`)*  
  The previous invocation of the foreign function associated with the current goal indicated it was possible to backtrack. The foreign function should produce the next solution.
- *Terminate call (`PL_PRUNED`)*  
  The choice point left by the foreign function has been destroyed by a cut or exception. The foreign function is given the opportunity to clean the environment. The context handle is the only meaningful argument -- the term arguments to the call are `(term_t)0`.

Both the context information and the type of call is provided by an argument of type `control_t` appended to the argument list for deterministic foreign functions. The macro [PL_foreign_control()](foreigninclude.html#PL_foreign_control()) extracts the type of call from the control argument. The foreign function can pass a context handle using the `PL_retry*()` macros and extract the handle from the extra argument using the `PL_foreign_context*()` macro.

`(return) foreign_t` **PL_retry**(`intptr_t value`)  
The foreign function succeeds while leaving a choice point. On backtracking over this goal the foreign function will be called again, but the control argument now indicates it is a‘Redo’call and the macro [PL_foreign_context()](foreigninclude.html#PL_foreign_context()) returns the handle passed via [PL_retry()](foreigninclude.html#PL_retry()). This handle is a signed value two bits smaller than a pointer, i.e., 30 or 62 bits (two bits are used for status indication). Defined as `return _`[`PL_retry(n)`](foreigninclude.html#PL_retry()). See also **PL_succeed()**.

`(return) foreign_t` **PL_retry_address**(`void *`)  
As [PL_retry()](foreigninclude.html#PL_retry()), but ensures an address as returned by **malloc()** is correctly recovered by [PL_foreign_context_address()](foreigninclude.html#PL_foreign_context_address()). Defined as `return _`[`PL_retry_address(n)`](foreigninclude.html#PL_retry_address()). See also **PL_succeed()**.

`int` **PL_foreign_control**(`control_t`)  
Extracts the type of call from the control argument. The return values are described above. Note that the function should be prepared to handle the `PL_PRUNED` case and should be aware that the other arguments are not valid in this case.

`intptr_t` **PL_foreign_context**(`control_t`)  
Extracts the context from the context argument. If the call type is `PL_FIRST_CALL` the context value is 0L. Otherwise it is the value returned by the last [PL_retry()](foreigninclude.html#PL_retry()) associated with this goal (both if the call type is `PL_REDO` or `PL_PRUNED`).

`void *` **PL_foreign_context_address**(`control_t`)  
Extracts an address as passed in by [PL_retry_address()](foreigninclude.html#PL_retry_address()).

`predicate_t` **PL_foreign_context_predicate**(`control_t`)  
Fetch the Prolog predicate that is executing this function. Note that if the predicate is imported, the returned predicate refers to the final definition rather than the imported predicate; i.e., the module reported by [PL_predicate_info()](foreigninclude.html#PL_predicate_info()) is the module in which the predicate is defined rather than the module where it was called. See also [PL_predicate_info()](foreigninclude.html#PL_predicate_info()).

Note: If a non-deterministic foreign function returns using **PL_succeed()** or **PL_fail()**, Prolog assumes the foreign function has cleaned its environment. **No** call with control argument `PL_PRUNED` will follow.

The code of [figure 5](foreigninclude.html#fig:nondetermf) shows a skeleton for a non-deterministic foreign predicate definition.

``` code
typedef struct                  /* define a context structure */
{ ...
} context;

foreign_t
my_function(term_t a0, term_t a1, control_t handle)
{ struct context * ctxt;

  switch( PL_foreign_control(handle) )
  { case PL_FIRST_CALL:
        if ( !(ctxt = malloc(sizeof *ctxt)) )
          return PL_resource_error("memory");
        <initialize ctxt>
        break;
    case PL_REDO:
        ctxt = PL_foreign_context_address(handle);
        break;
    case PL_PRUNED:
        ctxt = PL_foreign_context_address(handle);
        ...
        free(ctxt);
        return TRUE;
  }

  <find first/next solution from ctxt>
  ...
  // We fail */
  if ( <no_solution> )
  { free(ctx);
    return FALSE;
  }
  // We succeed without a choice point */
  if ( <last_solution> )
  { free(ctx);
    return TRUE;
  }
  // We succeed with a choice point */
  PL_retry_address(ctxt);
}
```

**Figure 5 :** Skeleton for non-deterministic foreign functions

#### 12.4.1.2 Yielding from foreign predicates

Starting with SWI-Prolog 8.5.5 we provide an experimental interface that allows using a SWI-Prolog engine for asynchronous processing. The idea is that an engine that calls a foreign predicate which would need to block may be suspended and later resumed. For example, consider an application that listens to a large number of network connections (sockets). SWI-Prolog offers three scenarios to deal with this:

1.  Using a thread per connection. This model fits Prolog well as it allows to keep state in e.g. a DCG using [phrase_from_stream/2](pio.html#phrase_from_stream/2). Maintaining an operating system thread per connection uses a significant amount of resources though.
2.  Using [wait_for_input/3](streamstat.html#wait_for_input/3) a single thread can wait for many connections. Each time input arrives we must associate this with a *state engine* and advance this engine using a chunk of input of unknown size. Programming a state engine in Prolog is typically a tedious job. Although we can use *delimited continuations* (see [section 4.9](delcont.html#sec:4.9)) in some scenarios this is not a universal solution.
3.  Using the primitives from this section we can create an *engine* (see **PL_engine_create()**) to handle a connection with the same benefits as using threads. When the engine calls a foreign predicate that would need to block it calls [PL_yield_address()](foreigninclude.html#PL_yield_address()) to suspend the engine. An overall scheduler watches for ready connections and calls [PL_next_solution()](foreigninclude.html#PL_next_solution()) to resume the suspended engine. This approach allows processing many connections on the same operating system thread.

As is, these features can only used through the foreign language interface. It was added after discussion with with Mattijs van Otterdijk aiming for using SWI-Prolog together with Rust's [asynchronous programming](https://rust-lang.github.io/async-book/01_getting_started/01_chapter.html) support. Note that this feature is related to the engine API as described in [section 11](engines.html#sec:11). It is different though. Where the Prolog engine API allows for communicating with a Prolog engine, the facilities of this section merely allow an engine to suspend, to be resumed later.

To prepare a query for asynchronous usage we first create an engine using [PL_create_engine()](foreignthread.html#PL_create_engine()). Next, we create a query in the engine using [PL_open_query()](foreigninclude.html#PL_open_query()) with the flags `PL_Q_ALLOW_YIELD` and `PL_Q_EXT_STATUS`. A foreign predicate that needs to be capable of suspending must be registered using [PL_register_foreign()](foreigninclude.html#PL_register_foreign()) and the flags `PL_FA_VARARGS` and `PL_FA_NONDETERMINISTIC`; i.e., only non-det predicates can yield. This is no restriction as non-det predicate can always return `TRUE` to indicate deterministic success. Finally, [PL_yield_address()](foreigninclude.html#PL_yield_address()) allows the predicate to yield control, preparing to resume similar to [PL_retry_address()](foreigninclude.html#PL_retry_address()) does for non-deterministic results. [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns `PL_S_YIELD` if a predicate calls [PL_yield_address()](foreigninclude.html#PL_yield_address()) and may be resumed by calling [PL_next_solution()](foreigninclude.html#PL_next_solution()) using the same query id (`qid`). We illustrate the above using some example fragments.

First, let us create a predicate that can read the available input from a Prolog stream and yield if it would block. Note that our predicate *must* the `PL_FA_VARARGS` interface, which implies the first argument is in `a0`, the second in `a0+1`, etc.^(216the other foreign interfaces do not support the yield API.)

``` code
/** read_or_block(+Stream, -String) is det.
*/

#define BUFSIZE 4096

static foreign_t
read_or_block(term_t a0, int arity, void *context)
{ IOSTREAM *s;

  switch(PL_foreign_control(context))
  { case PL_FIRST_CALL:
      if ( PL_get_stream(a0, &s, SIO_INPUT) )
      { Sset_timeout(s, 0);
        break;
      }
      return FALSE;
    case PL_RESUME:
      s = PL_foreign_context_address(context);
      break;
    case PL_PRUNED:
      return PL_release_stream(s);
    default:
      assert(0);
      return FALSE;
  }

  char buf[BUFSIZE];

  size_t n = Sfread(buf, sizeof buf[0], sizeof buf / sizeof buf[0], s);
  if ( n == 0 )                     // timeout or error
  { if ( (s->flags&SIO_TIMEOUT) )
      PL_yield_address(s);        // timeout: yield
    else
      return PL_release_stream(s);  // raise error
  } else
  { return ( PL_release_stream(s) &&
             PL_unify_chars(a0+1, PL_STRING|REP_ISO_LATIN_1, n, buf) );
  }
}
```

This function must be registered using [PL_register_foreign()](foreigninclude.html#PL_register_foreign()):

``` code
  PL_register_foreign("read_or_block", 2, read_or_block,
                      PL_FA_VARARGS|PL_FA_NONDETERMINISTIC);
```

Next, create an engine to run handle_connection/1 on a Prolog stream. Note that we omitted most of the error checking for readability. Also note that we must make our engine *current* using [PL_set_engine()](foreignthread.html#PL_set_engine()) before we can interact with it.

``` code
qid_t
start_connection(IOSTREAM *c)
{ predicate_t p = PL_predicate("handle_connection", 1, "user");

  PL_engine_t e = PL_create_engine(NULL);
  qid_t q = NULL;
  PL_WITH_ENGINE(e)
  { term_t av = PL_new_term_refs(1);
    PL_unify_stream(av+0, c);
    q = PL_open_query(e, NULL,
                      PL_Q_CATCH_EXCEPTION|
                      PL_Q_ALLOW_YIELD|
                      PL_Q_EXT_STATUS,
                      p, av);
  }
  return q;
}
```

Finally, our foreign code must manage this engine. Normally it will do so together with many other engines. First, we write a function that runs a query in the engine to which it belongs.^(217Possibly, future versions of [PL_next_solution()](foreigninclude.html#PL_next_solution()) may do that although the value is in general limited because interacting with the arguments of the query requires the query's engine to be current anyway.)

``` code
int
PL_engine_next_solution(qid_t qid)
{ int rc = FALSE;

  PL_WITH_ENGINE(PL_query_engine(qid))
  { rc = PL_next_solution(qid);
  }

  return rc;
}
```

Now we can simply handle a connection using the loop below which restarts the query as long as it yields. Realistic code manages multiple queries and will (in this case) use the POSIX **poll()** or **select()** interfaces to activate the next query that can continue without blocking.

``` code
  int rc;
  do
  { rc = PL_engine_next_solution(qid);
  } while( rc == PL_S_YIELD );
```

After the query completes it must be closed using [PL_close_query()](foreigninclude.html#PL_close_query()) or [PL_cut_query()](foreigninclude.html#PL_cut_query()). The engine may be destroyed using **PL_engine_destroy()** or reused for a new query.

`(return) foreign_t` **PL_yield_address**(`void *`)  
Cause [PL_next_solution()](foreigninclude.html#PL_next_solution()) of the active query to return with `PL_S_YIELD`. A subsequent call to [PL_next_solution()](foreigninclude.html#PL_next_solution()) on the same query calls the foreign predicate again with the control status set to `PL_RESUME`, after which [PL_foreign_context_address()](foreigninclude.html#PL_foreign_context_address()) retrieves the address passed to this function. The state of the Prolog engine is maintained, including `term_t` handles. If the passed address needs to be invalidated the predicate must do so when returning either `TRUE` or `FALSE`. If the engine terminates the predicate the predicate is called with status `PL_PRUNED`, in which case the predicate must cleanup.

`bool` **PL_can_yield**(`void`)  
Returns `TRUE` when called from inside a foreign predicate if the query that (indirectly) calls this foreign predicate can yield using [PL_yield_address()](foreigninclude.html#PL_yield_address()). Returns `FALSE` when either there is no current query or the query cannot yield.

**Discussion**

Asynchronous processing has become popular with modern programming languages, especially those aiming at network communication. Asynchronous processing uses fewer resources than threads while avoiding most of the complications associated with thread synchronization if only a single thread is used to manage the various states. The lack of good support for destructive state updates in Prolog makes it attractive to use threads for dealing with multiple inputs. The fact that Prolog discourages using shared global data such as dynamic predicates typically makes multithreaded code easy to manage.

It is not clear how much scalability we gain using Prolog engines instead of Prolog threads. The only difference between the two is the operating system task. Prolog engines are still rather memory intensive, mainly depending on the stack sizes. Global garbage collection (atoms and clauses) need to process all the stacks of all the engines and thus limit scalability.

One possible future direction is to allow all (possibly) blocking Prolog predicates to use the yield facility and provide a Prolog API to manage sets of engines that use this type of yielding. As is, these features are designed to allow SWI-Prolog for cooperating with languages that provide asynchronous functions.

#### 12.4.1.3 Implementing a yield based debugger

Starting with version 9.3.21, the SWI-Prolog foreign interface allows you to implement a debugger based on *yielding*, similar to [section 12.4.1.2](foreigninclude.html#sec:12.4.1.2). This implies that instead of using the built-in debugger or the [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4) hook that is used for e.g., driving the GUI debugger, [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns when trace interaction is required. The embedding system can interact with the user and resume [PL_next_solution()](foreigninclude.html#PL_next_solution()). This is notably required by the WASM described in [section 13](wasm-version.html#sec:13) version, where we need to return to the browser event loop to interact with the user. Yielding on behalf of the debugger is enabled using the `PL_Q_TRACE_WITH_YIELD` flag of [PL_open_query()](foreigninclude.html#PL_open_query()). The skeleton or using these facilities is below. [PL_get_trace_context()](foreigninclude.html#PL_get_trace_context()) retrieves a Prolog term that provides information similar to [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4).

``` code
  qid_t qid = PL_open_query(module,
                            PL_Q_EXT_STATUS|
                            PL_Q_ALLOW_YIELD|PL_Q_TRACE_WITH_YIELD,
                            predicate, argv);
  for(;;)
  { int rc = PL_next_solution(qid);
    switch(rc)
    { case PL_S_YIELD_DEBUG:
        PL_get_trace_context(msg);
        trace_reply(msg, action);
        PL_set_trace_action(action);
        break;
      case ...
    }
  }
```

The **trace_reply()** is a user defined function that typically calls Prolog to display the current goal and requests the user how to continue. **PL_set_trace_action()** processes the same term as output argument of [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4), telling Prolog how to continue. One option for **trace_reply()** is to call a Prolog predicate. Below is a partial implementation for such a predicate. Here, trace_reply/3 uses the character typed by the user to derive the trace continuation action.

``` code
read_trace_reply(Msg, Action) :-
    print_message(debug, Msg),
    format(user_error, '? ', []),
    get_single_char(Code),
    char_code(Char, Code),
    trace_reply(Char, Msg, Action).
```

`bool` **PL_get_trace_context**(`term_t -msg`)  
Unify `msg` with a Prolog term of the following shape: `frame(Frame, Choice, Port, PC)`. This provides the same information as passed to [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4). The additional `PC` argument provides the program counter in the executing clause. The `msg` term may be passed to [print_message/2](printmsg.html#print_message/2) to print the current port and goal, similarly to the commandline debugger.

`bool` **PL_set_trace_action, term_t +action**(`S`)  
et the continuation action. Supported values are described with the `Action` argument of [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4).

### 12.4.2 Atoms and functors

The following functions provide for communication using atoms and functors.

`atom_t` **PL_new_atom**(`const char *`)  
Return an atom handle for the given C-string. This function always succeeds. The returned handle is valid as long as the atom is referenced (see [section 12.4.2.1](foreigninclude.html#sec:12.4.2.1)). Currently aborts the process with a *fatal error* on failure. Future versions may raise a resource exception and return `(atom_t)0`.

The following atoms are provided as macros, giving access to the empty list symbol and the name of the list constructor. Prior to version 7, `ATOM_nil` is the same as [`PL_new_atom("[]")`](foreigninclude.html#PL_new_atom()) and `ATOM_dot` is the same as [`PL_new_atom(".")`](foreigninclude.html#PL_new_atom()). This is no longer the case in SWI-Prolog version 7.

`atom_t` **ATOM_nil**(`A`)  
tomic constant that represents the empty list. It is advised to use [PL_get_nil()](foreigninclude.html#PL_get_nil()), [PL_put_nil()](foreigninclude.html#PL_put_nil()) or [PL_unify_nil()](foreigninclude.html#PL_unify_nil()) where applicable.

`atom_t` **ATOM_dot**(`A`)  
tomic constant that represents the name of the list constructor. The list constructor itself is created using [`PL_new_functor(ATOM_dot,2)`](foreigninclude.html#PL_new_functor()). It is advised to use [PL_get_list()](foreigninclude.html#PL_get_list()), [PL_put_list()](foreigninclude.html#PL_put_list()) or [PL_unify_list()](foreigninclude.html#PL_unify_list()) where applicable.

`atom_t` **PL_new_atom_mbchars**(`int rep, size_t len, const char *s`)  
This function generalizes [PL_new_atom()](foreigninclude.html#PL_new_atom()) and [PL_new_atom_nchars()](foreigninclude.html#PL_new_atom_nchars()) while allowing for multiple encodings. The `rep` argument is one of `REP_ISO_LATIN_1`, `REP_UTF8` or `REP_MB`. If `len` is `(size_t)-1`, it is computed from `s` using **strlen()**. Raises an exception if `s` violates `rep` and returns `(atom_t)0`. For other error conditions, see [PL_new_atom()](foreigninclude.html#PL_new_atom()).

`bool` **PL_atom_mbchars**(`atom_t atom, size_t len, char *s, unsigned int flags`)  
This function generalizes fetching the text associated with an atom. The encoding depends on the flags `REP_UTF8`, `REP_MB` or `REP_ISO_LATIN_1`. Storage is defined by the `BUF_*` flags as described with [PL_get_chars()](foreigninclude.html#PL_get_chars()). The flag `CVT_EXCEPTION` defines whether or not the function fails silently or raises a Prolog exception. This function may fail because `atom` is not a text atom but a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)), conversion to the requested encoding is not possible or a resource error occurs.

`const char*` **PL_atom_chars**(`atom_t atom`)  
Deprecated. This function returns a pointer to the content represented by the atom or blob regardless of its type. New code that uses blobs should use the blob functions such as [PL_blob_data()](foreigninclude.html#PL_blob_data()) to get a pointer to the content, the size of the content, and the type of the content. Most applications that need to get text from a `term_t` handle should use [PL_atom_nchars()](foreigninclude.html#PL_atom_nchars()), [PL_atom_wchars()](foreigninclude.html#PL_atom_wchars()), or [PL_atom_mbchars()](foreigninclude.html#PL_atom_mbchars()). If it is *known* that `atom` is a classical Prolog text atom, one can use [PL_atom_nchars()](foreigninclude.html#PL_atom_nchars()) to obtain the C string and its length (for ISO-Latin-1 atoms) or [PL_atom_wchars()](foreigninclude.html#PL_atom_wchars()) to obtain a C wide string (`wchar_t`).

`size_t` **PL_atom_index**(`atom_t atom`)  
Extract the *index* of an atom. This is a relatively small integer. Atoms are numbered sequentially, starting at one (1). Note that the sequence may have holes due to atom garbage collection. The released index may later be reused for a new atom. The index may be used as a compact identifier for the atom. Extracting the index has no impact on the lifetime of the atom, i.e., the index is valid as long as the `atom_t` is valid.

`atom_t` **PL_atom_from_index**(`size_t index`)  
Recover an atom from its index as obtained by [PL_atom_index()](foreigninclude.html#PL_atom_index()).

`functor_t` **PL_new_functor**(`atom_t name, int arity`)  
Returns a *functor identifier*, a handle for the name/arity pair. The returned handle is valid for the entire Prolog session. Future versions may garbage collect functors as part of atom garbage collection. Currently aborts the process with a *fatal error* on failure. Future versions may raise a resource exception and return `(atom_t)0`.

`atom_t` **PL_functor_name**(`functor_t f`)  
Return an atom representing the name of the given functor.

`size_t` **PL_functor_arity**(`functor_t f`)  
Return the arity of the given functor.

#### 12.4.2.1 Atoms and atom garbage collection

With the introduction of atom garbage collection in version 3.3.0, atoms no longer live as long as the process. Instead, their lifetime is guaranteed only as long as they are referenced. In the single-threaded version, atom garbage collections are only invoked at the *call-port*. In the multithreaded version (see [chapter 10](threads.html#sec:10)), they appear asynchronously, except for the invoking thread.

For dealing with atom garbage collection, two additional functions are provided:

`void` **PL_register_atom**(`atom_t atom`)  
Increment the reference count of the atom by one. [PL_new_atom()](foreigninclude.html#PL_new_atom()) performs this automatically, returning an atom with a reference count of at least one.^(218Otherwise asynchronous atom garbage collection might destroy the atom before it is used.)

`void` **PL_unregister_atom**(`atom_t atom`)  
Decrement the reference count of the atom. If the reference count drops below zero, an assertion error is raised.

Please note that the following two calls are different with respect to atom garbage collection:

``` code
PL_unify_atom_chars(t, "text");
PL_unify_atom(t, PL_new_atom("text"));
```

The latter increments the reference count of the atom `text`, which effectively ensures the atom will never be collected. It is advised to use the \***\_chars()** or \***\_nchars()** functions whenever applicable.

### 12.4.3 Input and output

For input and output, `SWI-Stream.h` defines a set of functions that are similar to the C library functions, except prefixed by **S**, e.g. [Sfprintf()](foreign-streams.html#Sfprintf()). They differ from the C functions in following ways:

- Instead of returning the number of bytes written and a negative value for error, they return the number of characters written and a negative value for error.
- Instead of a `FILE`, they access the Prolog streams, using `IOSTREAM*`. In particular, `Scurrent_output` accesses the current output stream and works well with [with_output_to/2](IO.html#with_output_to/2). Similarly, there are `Scurrent_intput`, `Suser_output`, `Suser_error`, and `Suser_input`.
- If you wish to directly use the operating system's `stdin`, `stdout`, `stderr`, you can use `Sinput`, `Soutput`, `Serror`. These are not affected by predicates such as [with_output_to/2](IO.html#with_output_to/2).

In general, if a stream is acquired via [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()), an error is raised when [PL_release_stream()](foreign-streams.html#PL_release_stream()) is called, so in that situation, there's no need to check the return codes from the IO functions. Blob write callbacks are also called in the context of an acquired stream, so there is no need to check the return codes from its IO function calls. However, if you use one of the standard streams such as `Scurrent_output`, you should check the return code and return `FALSE` from the foreign predicate, at which point an error will be raised. Not all IO functions follow this, because they need to return other information, so you should check the details with each one (e.g., [Sputcode()](foreign-streams.html#Sputcode()) returns -1 on error).

For more details, including formatting extensions for printing terms, see [section 12.9](foreign-streams.html#sec:12.9).

### 12.4.4 Analysing Terms via the Foreign Interface

Each argument of a foreign function (except for the control argument) is of type `term_t`, an opaque handle to a Prolog term. Three groups of functions are available for the analysis of terms. The first just validates the type, like the Prolog predicates [var/1](typetest.html#var/1), [atom/1](typetest.html#atom/1), etc., and are called `PL_is_*()`. The second group attempts to translate the argument into a C primitive type. These predicates take a `term_t` and a pointer to the appropriate C type and return `TRUE` or `FALSE` depending on successful or unsuccessful translation. If the translation fails, the pointed-to data is never modified.

#### 12.4.4.1 Testing the type of a term

`int` **PL_term_type**(`term_t`)  
Obtain the type of a term, which should be a term returned by one of the other interface predicates or passed as an argument. The function returns the type of the Prolog term. The type identifiers are listed below. Note that the extraction functions `PL_get_*()` also validate the type and thus the two sections below are equivalent.

``` code
        if ( PL_is_atom(t) )
        { char *s;

          PL_get_atom_chars(t, &s);
          ...;
        }

or

        char *s;
        if ( PL_get_atom_chars(t, &s) )
        { ...;
        }
```

**Version 7** added `PL_NIL`, `PL_BLOB`, `PL_LIST_PAIR` and `PL_DICT`. Older versions classify `PL_NIL` and `PL_BLOB` as `PL_ATOM`, `PL_LIST_PAIR` as `PL_TERM` and do not have dicts.

|  |  |
|----|----|
| `PL_VARIABLE` | A variable or attributed variable |
| `PL_ATOM` | A Prolog atom |
| `PL_NIL` | The constant `[]` |
| `PL_BLOB` | A blob (see [section 12.4.10.2](foreigninclude.html#sec:12.4.10.2)) |
| `PL_STRING` | A string (see [section 5.2](string.html#sec:5.2)) |
| `PL_INTEGER` | A integer |
| `PL_RATIONAL` | A rational number |
| `PL_FLOAT` | A floating point number |
| `PL_TERM` | A compound term |
| `PL_LIST_PAIR` | A list cell (`[H|T]`) |
| `PL_DICT` | A dict (see [section 5.4](bidicts.html#sec:5.4))) |

The functions PL_is\_\<`type`\> are an alternative to [PL_term_type()](foreigninclude.html#PL_term_type()). The test [`PL_is_variable(term)`](foreigninclude.html#PL_is_variable()) is equivalent to [`PL_term_type(term)`](foreigninclude.html#PL_term_type())` == PL_VARIABLE`, but the first is considerably faster. On the other hand, using a switch over [PL_term_type()](foreigninclude.html#PL_term_type()) is faster and more readable then using an if-then-else using the functions below. All these functions return either `TRUE` or `FALSE`.

`bool` **PL_is_variable**(`term_t`)  
Returns non-zero if `term` is a variable.

`bool` **PL_is_ground**(`term_t`)  
Returns non-zero if `term` is a ground term. See also [ground/1](typetest.html#ground/1). This function is cycle-safe.

`bool` **PL_is_atom**(`term_t`)  
Returns non-zero if `term` is an atom.

`bool` **PL_is_string**(`term_t`)  
Returns non-zero if `term` is a string.

`bool` **PL_is_integer**(`term_t`)  
Returns non-zero if `term` is an integer.

`bool` **PL_is_rational**(`term_t`)  
Returns non-zero if `term` is a rational number (`P/Q`). Note that all integers are considered rational and this test thus succeeds for any term for which [PL_is_integer()](foreigninclude.html#PL_is_integer()) succeeds. See also [PL_get_mpq()](foreigninclude.html#PL_get_mpq()) and [PL_unify_mpq()](foreigninclude.html#PL_unify_mpq()).

`bool` **PL_is_float**(`term_t`)  
Returns non-zero if `term` is a float. Note that the corresponding [PL_get_float()](foreigninclude.html#PL_get_float()) converts rationals (and thus integers).

`bool` **PL_is_callable**(`term_t`)  
Returns non-zero if `term` is a callable term. See [callable/1](typetest.html#callable/1) for details.

`bool` **PL_is_compound**(`term_t`)  
Returns non-zero if `term` is a compound term.

`bool` **PL_is_functor**(`term_t, functor_t`)  
Returns non-zero if `term` is compound and its functor is `functor`. This test is equivalent to [PL_get_functor()](foreigninclude.html#PL_get_functor()), followed by testing the functor, but easier to write and faster.

`bool` **PL_is_list**(`term_t`)  
Returns non-zero if `term` is a compound term using the list constructor or the list terminator. See also [PL_is_pair()](foreigninclude.html#PL_is_pair()) and [PL_skip_list()](foreigninclude.html#PL_skip_list()).

`bool` **PL_is_pair**(`term_t`)  
Returns non-zero if `term` is a compound term using the list constructor. See also [PL_is_list()](foreigninclude.html#PL_is_list()) and [PL_skip_list()](foreigninclude.html#PL_skip_list()).

`bool` **PL_is_dict**(`term_t`)  
Returns non-zero if `term` is a dict. See also [PL_put_dict()](foreigninclude.html#PL_put_dict()) and [PL_get_dict_key()](foreigninclude.html#PL_get_dict_key()).

`bool` **PL_is_atomic**(`term_t`)  
Returns non-zero if `term` is atomic (not a variable or compound).

`bool` **PL_is_number**(`term_t`)  
Returns non-zero if `term` is an rational (including integers) or float.

`bool` **PL_is_acyclic**(`term_t`)  
Returns non-zero if `term` is acyclic (i.e. a finite tree).

#### 12.4.4.2 Reading data from a term

The functions `PL_get_*()` read information from a Prolog term. Most of them take two arguments. The first is the input term and the second is a pointer to the output value or a term reference. The return value is `TRUE` or `FALSE`, indicating the success of the "get" operation. Most functions have a related "\_ex" function that raises an error if the argument is the operation cannot be completed. If the Prolog term is not suitable, this is a type, domain or instantiation error. If the receiving C type cannot represent the value this is a representation error.

For integers an alternative interface exists, which helps deal with the various integer types in C and C++. They are convenient for use with `_Generic` selection or C++ overloading.

`bool` **PL_get_atom**(`term_t +t, atom_t *a`)  
If `t` is an atom, store the unique atom identifier over `a`. See also [PL_atom_chars()](foreigninclude.html#PL_atom_chars()) and [PL_new_atom()](foreigninclude.html#PL_new_atom()). If there is no need to access the data (characters) of an atom, it is advised to manipulate atoms using their handle. As the atom is referenced by `t`, it will live at least as long as `t` does. If longer lifetime is required, the atom should be locked using [PL_register_atom()](foreigninclude.html#PL_register_atom()).

`bool` **PL_get_atom_chars**(`term_t +t, char **s`)  
If `t` is an atom, store a pointer to a 0-terminated C-string in `s`. It is explicitly **not** allowed to modify the contents of this string. Some built-in atoms may have the string allocated in read-only memory, so‘temporary manipulation’can cause an error.

`int` **PL_get_string_chars**(`term_t +t, char **s, size_t *len`)  
If `t` is a string object, store a pointer to a 0-terminated C-string in `s` and the length of the string in `len`. Note that this pointer is invalidated by backtracking, garbage collection and stack-shifts, so generally the only safe operations are to pass it immediately to a C function that doesn't involve Prolog.

`bool` **PL_get_chars**(`term_t +t, char **s, unsigned flags`)  
Convert the argument term `t` to a 0-terminated C-string. `flags` is a bitwise disjunction from two groups of constants. The first specifies which term types should be converted and the second how the argument is stored. Below is a specification of these constants. `BUF_STACK` implies, if the data is not static (as from an atom), that the data is pushed on a stack. If `BUF_MALLOC` is used, the data must be freed using [PL_free()](foreignnotes.html#PL_free()) when no longer needed.

With the introduction of wide characters (see [section 2.18.1](widechars.html#sec:2.18.1)), not all atoms can be converted into a `char*`. This function fails if `t` is of the wrong type, but also if the text cannot be represented. See the `REP_*` flags below for details. See also [PL_get_wchars()](foreigninclude.html#PL_get_wchars()) and [PL_get_nchars()](foreigninclude.html#PL_get_nchars()).

The first set of flags (`CVT_ATOM` through `CVT_VARIABLE`, if set, are tested in order, using the first that matches. If none of these match, then a check is made for one of `CVT_WRITE`, `CVT_WRITE_CANONICAL`, `CVT_WRITEQ` being set. If none of the “CVT_WRITE\*” flags are set, then a `type_error` is raised.

**CVT_ATOM**  
Convert if term is an atom.

**CVT_STRING**  
Convert if term is a string.

**CVT_LIST**  
Convert if term is a list of characters (atoms of length 1) or character codes (integers representing Unicode code points).

**CVT_INTEGER**  
Convert if term is an integer.

**CVT_RATIONAL**  
Convert if term is a rational number (including integers). Non-integral numbers are written as \<`num`\>r\<`den`\>.

**CVT_XINTEGER**  
Convert if term is an integer to hexadecimal notation. May be combined with `CVT_RATIONAL` to represent rational numbers using hexadecimal notation. Hexadecimal notation is notably useful for transferring big integers to other programming environments if the target system can read hexadecimal notation because the result is both more compact and faster to write and read.

**CVT_FLOAT**  
Convert if term is a float. The characters returned are the same as [write/1](termrw.html#write/1) would write for the floating point number.

**CVT_NUMBER**  
Convert if term is an integer, rational number or float. Equivalent to `CVT_RATIONAL``|``CVT_FLOAT`. Note that `CVT_INTEGER` is implied by `CVT_RATIONAL`.

**CVT_ATOMIC**  
Convert if term is atomic. Equivalent to `CVT_NUMBER``|``CVT_ATOM``|``CVT_STRING`.

**CVT_ALL**  
Convert if term is any of the above. Integers and rational numbers are written as decimal (i.e., `CVT_XINTEGER` is *not* implied). Note that this does not include variables or terms (with the exception of a list of characters/codes). Equivalent to `CVT_ATOMIC``|``CVT_LIST`.

**CVT_VARIABLE**  
Convert variable to print-name (e.g., `_3290`).

**CVT_WRITE**  
Convert any term that is not converted by any of the other flags using [write/1](termrw.html#write/1). If no `BUF_*` is provided, `BUF_STACK` is implied.

**CVT_WRITEQ**  
As `CVT_WRITE`, but using [writeq/2](termrw.html#writeq/2).

**CVT_WRITE_CANONICAL**  
As `CVT_WRITE`, but using [write_canonical/2](termrw.html#write_canonical/2).

**CVT_EXCEPTION**  
If conversion fails due to a type error, raise a Prolog type error exception in addition to failure.

**BUF_DISCARDABLE**  
Data must be copied immediately. Use with caution - you must be certain that no PL\_\*() function is called before you use the memory because if garbage collection occurs, the value of `s` won't be properly updated - `BUF_STACK` is safe, but slower. Also, some flags can cause `BUF_DISCARDABLE` to be treated as `BUF_STACK`. See [section 12.4.14](foreigninclude.html#sec:12.4.14).

**BUF_STACK**  
Data is stored on a stack. The older `BUF_RING` is an alias for `BUF_STACK`. `BUF_STACK` is the default if `BUF_MALLOC` is not specified. `BUF_STACK` should not be used inside a loop, unless it is within a (`PL_STRINGS_MARK`-`PL_STRINGS_RELEASE`). See [section 12.4.14](foreigninclude.html#sec:12.4.14).

**BUF_MALLOC**  
Data is copied to a new buffer returned by **PL_malloc**(3). When no longer needed the user must call [PL_free()](foreignnotes.html#PL_free()) on the data. If `BUF_MALLOC` is specified, `BUF_STACK` is ignored. See [section 12.4.14](foreigninclude.html#sec:12.4.14).

**REP_ISO_LATIN_1**  
Text is in ISO Latin-1 encoding and the call fails if text cannot be represented. This flag has the value 0 and is thus the default.

**REP_UTF8**  
Convert the text to a UTF-8 string. This works for all text.

**REP_MB**  
Convert to default locale-defined 8-bit string. Success depends on the locale. Conversion is done using the **wcrtomb()** C library function.

`bool` **PL_get_list_chars**(`+term_t l, char **s, unsigned flags`)  
Same as [`PL_get_chars(``l``, ``s``, CVT_LIST|``flags``)`](foreigninclude.html#PL_get_chars()), provided `flags` contains none of the `CVT_*` flags.

`bool` **PL_get_integer**(`+term_t t, int *i`)  
If `t` is a Prolog integer, assign its value over `i`. On 32-bit machines, this is the same as [PL_get_long()](foreigninclude.html#PL_get_long()), but avoids a warning from the compiler. See also [PL_get_long()](foreigninclude.html#PL_get_long()) and [PL_get_integer_ex()](foreigninclude.html#PL_get_integer_ex()).

`bool` **PL_get_long**(`term_t +t, long *i`)  
If `t` is a Prolog integer that can be represented as a long, assign its value over `i`. If `t` is an integer that cannot be represented by a C long, this function returns `FALSE`. If `t` is a floating point number that can be represented as a long, this function succeeds as well. See also [PL_get_int64()](foreigninclude.html#PL_get_int64()) and [PL_get_long_ex()](foreigninclude.html#PL_get_long_ex()).

`bool` **PL_get_int64**(`term_t +t, int64_t *i`)  
If `t` is a Prolog integer or float that can be represented as a `int64_t`, assign its value over `i`. See also [PL_get_int64_ex()](foreigninclude.html#PL_get_int64_ex()).

`bool` **PL_get_uint64**(`term_t +t, uint64_t *i`)  
If `t` is a Prolog integer that can be represented as a `uint64_t`, assign its value over `i`. Note that this requires GMP support for representing `uint64_t` values with the high bit set. See also [PL_get_uint64_ex()](foreigninclude.html#PL_get_uint64_ex()).

`bool` **PL_get_intptr**(`term_t +t, intptr_t *i`)  
Get an integer that is at least as wide as a pointer. On most platforms this is the same as [PL_get_long()](foreigninclude.html#PL_get_long()), but on Win64 pointers are 8 bytes and longs only 4. Unlike [PL_get_pointer()](foreigninclude.html#PL_get_pointer()), the value is not modified.

`bool` **PL_get_bool**(`term_t +t, int *val`)  
If `t` has the value `true`, `false`, set `val` to the C constant `TRUE` or `FALSE` and return success, otherwise return failure. The values `on`, `1`, `off`, const0 and are also accepted.

`bool` **PL_get_pointer**(`term_t +t, void **ptr`)  
Together with [PL_put_pointer()](foreigninclude.html#PL_put_pointer()) and [PL_unify_pointer()](foreigninclude.html#PL_unify_pointer()), these functions allow representing a C pointer as a Prolog integer. The integer value is derived from the pointer, but not equivalent. The translation aims at producing smaller integers that fit more often in the *tagged* integer range. Representing C pointers as integers is *unsafe*. The *blob* API described in [section 12.4.10](foreigninclude.html#sec:12.4.10) provides a safe way for handling foreign resources that cooperates with Prolog garbage collection.

`bool` **PL_get_float**(`term_t +t, double *f`)  
If `t` is a float, integer or rational number, its value is assigned over `f`. Note that if `t` is an integer or rational conversion may fail because the number cannot be represented as a float.

`bool` **PL_get_functor**(`term_t +t, functor_t *f`)  
If `t` is compound or an atom, the Prolog representation of the name-arity pair will be assigned over `f`. See also [PL_get_name_arity()](foreigninclude.html#PL_get_name_arity()) and [PL_is_functor()](foreigninclude.html#PL_is_functor()).

`bool` **PL_get_name_arity**(`term_t +t, atom_t *name, size_t *arity`)  
If `t` is compound or an atom, the functor name will be assigned over `name` and the arity over `arity` (either or both may be NULL). See also [PL_get_compound_name_arity()](foreigninclude.html#PL_get_compound_name_arity()), [PL_get_functor()](foreigninclude.html#PL_get_functor()) and [PL_is_functor()](foreigninclude.html#PL_is_functor()).

`bool` **PL_get_compound_name_arity**(`term_t +t, atom_t *name, size_t *arity`)  
If `t` is compound term, the functor name will be assigned over `name` and the arity over `arity` (either or both may be `NULL`). This is the same as [PL_get_name_arity()](foreigninclude.html#PL_get_name_arity()), but this function fails if `t` is an atom.

`bool` **PL_get_module**(`term_t +t, module_t *module`)  
If `t` is an atom, the system will look up or create the corresponding module and assign an opaque pointer to it over *module*.

`bool` **PL_get_arg**(`size_t index, term_t +t, term_t -a`)  
If `t` is compound and index is between 1 and arity (inclusive), assign `a` with a term reference to the argument. Returns `FALSE` if `t` is not a compound or `index` is out of range (zero or higher than the arity of the compound). This function never raises a Prolog exception.

`int` **\_PL_get_arg**(`size_t index, term_t +t, term_t -a`)  
Same as [PL_get_arg()](foreigninclude.html#PL_get_arg()), but no checking is performed, neither whether `t` is actually a term nor whether `index` is a valid argument index. This is faster, but invalid usage leads to undefined behaviour.

`bool` **PL_get_dict_key**(`atom_t key, term_t +dict, term_t -value`)  
If `dict` is a dict, get the associated value in `value`. Fails silently if `key` does not appear in `dict` or if if `dict` is not a dict.

#### 12.4.4.3 Exchanging text using length and string

All internal text representation in SWI-Prolog is represented using `char *` plus length and allow for *0-bytes* in them. The foreign library supports this by implementing a \***\_nchars()** function for each applicable \***\_chars()** function. Below we briefly present the signatures of these functions. For full documentation consult the \***\_chars()** function.

`bool` **PL_get_atom_nchars**(`term_t t, size_t *len, char **s`)  
See [PL_get_atom_chars()](foreigninclude.html#PL_get_atom_chars()).

`bool` **PL_get_list_nchars**(`term_t t, size_t *len, char **s`)  
See [PL_get_list_chars()](foreigninclude.html#PL_get_list_chars()).

`bool` **PL_get_nchars**(`term_t t, size_t *len, char **s, unsigned int flags`)  
See [PL_get_chars()](foreigninclude.html#PL_get_chars()). The `len` pointer may be `NULL`.

`bool` **PL_put_atom_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_put_atom_chars()](foreigninclude.html#PL_put_atom_chars()).

`bool` **PL_put_string_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_put_string_chars()](foreigninclude.html#PL_put_string_chars()).

`bool` **PL_put_list_ncodes**(`term_t t, size_t len, const char *s`)  
See **PL_put_list_codes()**.

`bool` **PL_put_list_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_put_list_chars()](foreigninclude.html#PL_put_list_chars()).

`bool` **PL_unify_atom_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_unify_atom_chars()](foreigninclude.html#PL_unify_atom_chars()).

`bool` **PL_unify_string_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_unify_string_chars()](foreigninclude.html#PL_unify_string_chars()).

`bool` **PL_unify_list_ncodes**(`term_t t, size_t len, const char *s`)  
See **PL_unify_codes()**.

`bool` **PL_unify_list_nchars**(`term_t t, size_t len, const char *s`)  
See [PL_unify_list_chars()](foreigninclude.html#PL_unify_list_chars()).

In addition, the following functions are available for creating and inspecting atoms:

`atom_t` **PL_new_atom_nchars**(`size_t len, const char *s`)  
Create a new atom as [PL_new_atom()](foreigninclude.html#PL_new_atom()), but using the given length and characters. If `len` is `(size_t)-1`, it is computed from `s` using **strlen()**. See [PL_new_atom()](foreigninclude.html#PL_new_atom()) for error handling.

`const char *` **PL_atom_nchars**(`atom_t a, size_t *len`)  
Extract the text and length of an atom. If you do not need the length, pass NULL as the value of `len`. If [PL_atom_nchars()](foreigninclude.html#PL_atom_nchars()) is called for a blob, NULL is returned.

#### 12.4.4.4 Wide-character versions

Support for exchange of wide-character strings is still under consideration. The functions dealing with 8-bit character strings return failure when operating on a wide-character atom or Prolog string object. The functions below can extract and unify both 8-bit and wide atoms and string objects. Wide character strings are represented as C arrays of objects of the type `pl_wchar_t`, which is guaranteed to be the same as `wchar_t` on platforms supporting this type. For example, on MS-Windows, this represents a 16-bit UTF-16 string, while using the GNU C library (glibc) this represents 32-bit UCS4 characters.

`atom_t` **PL_new_atom_wchars**(`size_t len, const pl_wchar_t *s`)  
Create atom from wide-character string as [PL_new_atom_nchars()](foreigninclude.html#PL_new_atom_nchars()) does for ISO-Latin-1 strings. If `s` only contains ISO-Latin-1 characters a normal byte-array atom is created. If `len` is `(size_t)-1`, it is computed from `s` using **wcslen()**. See [PL_new_atom()](foreigninclude.html#PL_new_atom()) for error handling.

`const pl_wchar_t*` **PL_atom_wchars**(`atom_t atom, size_t *len`)  
Extract characters from a wide-character atom. Succeeds on any atom marked as‘text’. If the underlying atom is a wide-character atom, the returned pointer is a pointer into the atom structure. If the atom is represented as an ISO-Latin-1 string, the returned pointer comes from Prolog's‘buffer stack’(see [section 12.4.14](foreigninclude.html#sec:12.4.14)).

`bool` **PL_get_wchars**(`term_t t, size_t *len, pl_wchar_t **s, unsigned flags`)  
Wide-character version of [PL_get_chars()](foreigninclude.html#PL_get_chars()). The `flags` argument is the same as for [PL_get_chars()](foreigninclude.html#PL_get_chars()). Note that this operation may return a pointer into Prolog's‘buffer stack’(see [section 12.4.14](foreigninclude.html#sec:12.4.14)).

`bool` **PL_put_wchars**(`term_t -t, int type, size_t len, const pl_wchar_t *s`)  
*Put* text from a wide character array in `t`. Arguments are the same as [PL_unify_wchars()](foreigninclude.html#PL_unify_wchars()).^(219The current implementation uses [PL_put_variable()](foreigninclude.html#PL_put_variable()) followed by [PL_unify_wchars()](foreigninclude.html#PL_unify_wchars()).)

`bool` **PL_unify_wchars**(`term_t +t, int type, size_t len, const pl_wchar_t *s`)  
Unify `t` with a textual representation of the C wide-character array `s`. The `type` argument defines the Prolog representation and is one of `PL_ATOM`, `PL_STRING`, `PL_CODE_LIST` or `PL_CHAR_LIST`.

`bool` **PL_unify_wchars_diff**(`term_t +t, term_t -tail, int type, size_t len, const pl_wchar_t *s`)  
Difference list version of [PL_unify_wchars()](foreigninclude.html#PL_unify_wchars()), only supporting the types `PL_CODE_LIST` and `PL_CHAR_LIST`. It serves two purposes. It allows for returning very long lists from data read from a stream without the need for a resizing buffer in C. Also, the use of difference lists is often practical for further processing in Prolog. Examples can be found in `packages/clib/readutil.c` from the source distribution.

#### 12.4.4.5 Reading a list

The functions from this section are intended to read a Prolog list from C. Suppose we expect a list of atoms; the code below will print the atoms, each on a line. Please note the following:

- We need a `term_t` *term reference* for the elements (`head`). This reference is reused for each element.
- We walk over the list using [PL_get_list_ex()](foreigninclude.html#PL_get_list_ex()) which overwrites the list `term_t`. As it is not allowed to overwrite the `term_t` passed in as arguments to a predicate, we must *copy* the argument `term_t`.
- SWI-Prolog atoms are Unicode objects. The [PL_get_chars()](foreigninclude.html#PL_get_chars()) returns a `char*`. We want it to convert atoms, return the result as a *multibyte* string (`REP_UTF8` may also be used) and finally we want an exception on type, instantiation or representation errors (if the system's default encoding cannot represent some characters of the Unicode atom). This may create temporary copies of the atom text - [PL_STRINGS_MARK()](foreigninclude.html#PL_STRINGS_MARK()) `...` [PL_STRINGS_RELEASE()](foreigninclude.html#PL_STRINGS_RELEASE()) handles that.
- The \***\_ex()** API functions are functionally the same as the ones without the `_ex` suffix, but they raise type, domain, or instantiation errors when the input is invalid; whereas the plain version may only raise resource exceptions if the request cannot be fulfilled due to resource exhaustion.
- [PL_get_nil_ex()](foreigninclude.html#PL_get_nil_ex()) is designed to propagate an already raised exception.

``` code
foreign_t
pl_write_atoms(term_t l)
{ term_t head = PL_new_term_ref();   /* the elements */
  term_t tail = PL_copy_term_ref(l); /* copy (we modify tail) */
  int rc = TRUE;

  while( rc && PL_get_list_ex(tail, head, tail) )
  { PL_STRINGS_MARK();
      char *s;
      if (rc=PL_get_chars(head, &s, CVT_ATOM|REP_MB|CVT_EXCEPTION)) )
        rc = Sfprintf(Scurrent_output, "%s\n", s);
    PL_STRINGS_RELEASE();
  }

  return rc && PL_get_nil_ex(tail); /* test end for [] */
}
```

Note that as of version 7, lists have a new representation unless the option **--traditional** is used. see [section 5.1](ext-lists.html#sec:5.1).

`bool` **PL_get_list**(`term_t +l, term_t -h, term_t -t`)  
If `l` is a list and not the empty list, assign a term reference to the head to `h` and to the tail to `t`.

`bool` **PL_get_head**(`term_t +l, term_t -h`)  
If `l` is a list and not the empty list, assign a term reference to the head to `h`.

`bool` **PL_get_tail**(`term_t +l, term_t -t`)  
If `l` is a list and not the empty list, assign a term reference to the tail to `t`.

`bool` **PL_get_nil**(`term_t +l`)  
Succeeds if `l` represents the list termination constant.

`int` **PL_skip_list**(`term_t +list, term_t -tail, size_t *len`)  
This is a multi-purpose function to deal with lists. It allows for finding the length of a list, checking whether something is a list, etc. The reference `tail` is set to point to the end of the list, `len` is filled with the number of list-cells skipped, and the return value indicates the status of the list:

**PL_LIST**  
The list is a‘proper’list: one that ends in the list terminator constant and `tail` is filled with the terminator constant.

**PL_PARTIAL_LIST**  
The list is a‘partial’list: one that ends in a variable and `tail` is a reference to this variable.

**PL_CYCLIC_TERM**  
The list is cyclic (e.g. X = \[a\|X\]). `tail` points to an arbitrary cell of the list and `len` is at most twice the cycle length of the list.

**PL_NOT_A_LIST**  
The term `list` is not a list at all. `tail` is bound to the non-list term and `len` is set to the number of list-cells skipped.

It is allowed to pass 0 for `tail` and `NULL` for `len`.

#### 12.4.4.6 Processing option lists and dicts

`bool` **PL_scan_options**(`term_t options, int flags, const char* opttype, PL_option_t specs[], ...`)  
Process an *option list* as we find with, e.g., [write_term/2](termrw.html#write_term/2) and many other builtin predicates. This function takes an option list (or dict) and in the variadic argument list pointers that receive the option values. [PL_scan_options()](foreigninclude.html#PL_scan_options()) takes care of validating the list, ensuring the list is not cyclic, validating the option type and storing the converted values using the supplied pointers.

Below is an example. While `PL_option_t` is a struct, its members are initialised using the **PL_OPTION()** macro. The data structure is not constant because [PL_scan_options()](foreigninclude.html#PL_scan_options()) adds the option names as *atoms* to speed up option processing. The macro PL_OPTIONS_END terminates the option list.

``` code
static PL_option_t mypred_options[] =
{ PL_OPTION("quoted",   OPT_BOOL),
  PL_OPTION("length",   OPT_SIZE),
  PL_OPTION("callback", OPT_TERM),
  PL_OPTIONS_END
};

static foreign_t
mypred(term_t a1, term_t options)
{ int    quoted   = FALSE;
  size_t length   = 10;
  term_t callback = 0;

  if ( !PL_scan_options(options, 0, "mypred_options", mypred_options,
                        &quoted, &length, &callback) )
    return FALSE;

  <implement mypred>
}
```

The only defined value for `flags` is currently `OPT_ALL`, which causes this function to raise a domain error if an option is passed that is not in `specs`. Default in SWI-Prolog is to silently ignore unknown options, unless the Prolog flag [iso](flags.html#flag:iso) is `true`. The `opttype` argument defines the type (group) of the options, e.g., `"write_option"`. Option *types* are defined by the ISO standard. SWI-Prolog only uses this if `OPT_ALL` is specified, to raise a `domain_error` of the indicated type if some option is unused. The type name is normally the name of the predicate followed by `_option` or the name of a representative of a group of predicates to which the options apply.

Defined option types and their corresponding pointer type are described below.

**`OPT_BOOL` `int`**  
Convert the option value to a bool. This converts the values described by [PL_get_bool()](foreigninclude.html#PL_get_bool()). In addition, an option without a value (i.e., a plain atom that denotes the option name) can act as a boolean `TRUE`.

**`OPT_INT` `int`**  
**`OPT_INT64` `int64_t`**  
**`OPT_UINT64` `uint64_t`**  
**`OPT_SIZE` `size_t`**  
**`OPT_DOUBLE` `double`**  
Numeric values of various types. Raises an error if the Prolog value cannot be represented by the C type.

**`OPT_STRING` `char*`**  
Uses [PL_get_chars()](foreigninclude.html#PL_get_chars()) using the flags `CVT_ALL|REP_UTF8|BUF_STACK|CVT_EXCEPTION`. The buffered string must be guarded using [PL_STRINGS_MARK()](foreigninclude.html#PL_STRINGS_MARK()) and [PL_STRINGS_RELEASE()](foreigninclude.html#PL_STRINGS_RELEASE()).

**`OPT_ATOM` `atom_t`**  
Accepts an atom. Note that if the C function that implements the predicate wishes to keep hold of the atom after it returns it must use [PL_register_atom()](foreigninclude.html#PL_register_atom()).

**`OPT_TERM` `term_t`**  
Accepts an arbitrary Prolog term. The term handle is scoped by the foreign predicate invocation. Terms can be preserved using [PL_record()](foreigninclude.html#PL_record()).

The ISO standard demands that if an option is repeated the *last* occurrence holds. This implies that [PL_scan_options()](foreigninclude.html#PL_scan_options()) must scan the option list to the end.

#### 12.4.4.7 An example: defining write/1 in C

[Figure 6](foreigninclude.html#fig:pl-display) shows a simplified definition of [write/1](termrw.html#write/1) to illustrate the described functions. This simplified version does not deal with operators. It is called display/1, because it mimics closely the behaviour of this Edinburgh predicate.

``` code
foreign_t
pl_display(term_t t)
{ functor_t functor;
  int arity, len, n;
  char *s;

  switch( PL_term_type(t) )
  { case PL_VARIABLE:
    case PL_ATOM:
    case PL_INTEGER:
    case PL_FLOAT:
      PL_get_chars(t, &s, CVT_ALL);
      if (! Sfprintf(Scurrent_output, "%s", s) )
        PL_fail;
      break;
    case PL_STRING:
      if ( !PL_get_string_chars(t, &s, &len) &&
           !Sfprintf(Scurrent_output, "\"%s\"", s) )
        PL_fail;
      break;
    case PL_TERM:
    { term_t a = PL_new_term_ref();

      if ( !PL_get_name_arity(t, &name, &arity) &&
           !Sfprintf(Scurrent_output, "%s(", PL_atom_chars(name)) )
        PL_fail
      for(n=1; n<=arity; n++)
      { if ( ! PL_get_arg(n, t, a) )
          PL_fail;
        if ( n > 1 )
          if ( ! Sfprintf(Scurrent_output, ", ") )
            PL_fail;
        if ( !pl_display(a) )
          PL_fail;
      }
      if ( !Sfprintf(Scurrent_output, ")") )
        PL_fail;
      break;
    default:
      PL_fail;                          /* should not happen */
    }
  }

  PL_succeed;
}
```

**Figure 6 :** A Foreign definition of display/1

### 12.4.5 Constructing Terms

Terms can be constructed using functions from the `PL_put_*()` and `PL_cons_*()` families. This approach builds the term‘inside-out’, starting at the leaves and subsequently creating compound terms. Alternatively, terms may be created‘top-down’, first creating a compound holding only variables and subsequently unifying the arguments. This section discusses functions for the first approach. This approach is generally used for creating arguments for [PL_call()](foreigninclude.html#PL_call()) and [PL_open_query()](foreigninclude.html#PL_open_query()).

`bool` **PL_put_variable**(`term_t -t`)  
Put a fresh variable in the term, resetting the term reference to its initial state.^(220Older versions created a variable on the global stack.)

`bool` **PL_put_atom**(`term_t -t, atom_t a`)  
Put an atom in the term reference from a handle. See also [PL_new_atom()](foreigninclude.html#PL_new_atom()) and [PL_atom_chars()](foreigninclude.html#PL_atom_chars()).

`bool` **PL_put_bool**(`term_t -t, int val`)  
Put one of the atoms `true` or `false` in the term reference See also [PL_put_atom()](foreigninclude.html#PL_put_atom()), [PL_unify_bool()](foreigninclude.html#PL_unify_bool()) and [PL_get_bool()](foreigninclude.html#PL_get_bool()).

`bool` **PL_put_chars**(`term_t -t, int flags, size_t len, const char *chars`)  
New function to deal with setting a term from a `char*` with various encodings. The `flags` argument is a bitwise *or* specifying the Prolog target type and the encoding of `chars`. A Prolog type is one of `PL_ATOM`, `PL_STRING`, `PL_CODE_LIST` or `PL_CHAR_LIST`. A representation is one of `REP_ISO_LATIN_1`, `REP_UTF8` or `REP_MB`. See [PL_get_chars()](foreigninclude.html#PL_get_chars()) for a definition of the representation types. If `len` is `-1` `chars` must be zero-terminated and the length is computed from `chars` using **strlen()**.

`bool` **PL_put_atom_chars**(`term_t -t, const char *chars`)  
Put an atom in the term reference constructed from the zero-terminated string. The string itself will never be referenced by Prolog after this function.

`bool` **PL_put_string_chars**(`term_t -t, const char *chars`)  
Put a zero-terminated string in the term reference. The data will be copied. See also [PL_put_string_nchars()](foreigninclude.html#PL_put_string_nchars()).

`bool` **PL_put_string_nchars**(`term_t -t, size_t len, const char *chars`)  
Put a string, represented by a length/start pointer pair in the term reference. The data will be copied. This interface can deal with 0-bytes in the string. See also [section 12.4.24](foreigninclude.html#sec:12.4.24).

`bool` **PL_put_list_chars**(`term_t -t, const char *chars`)  
Put a list of ASCII values in the term reference.

`bool` **PL_put_integer**(`term_t -t, long i`)  
Put a Prolog integer in the term reference.

`bool` **PL_put_int64**(`term_t -t, int64_t i`)  
Put a Prolog integer in the term reference.

`bool` **PL_put_uint64**(`term_t -t, uint64_t i`)  
Put a Prolog integer in the term reference. Note that unbounded integer support is required for `uint64_t` values with the highest bit set to 1. Without unbounded integer support, too large values raise a `representation_error` exception.

`bool` **PL_put_pointer**(`term_t -t, void *ptr`)  
Put a Prolog integer in the term reference. Provided `ptr` is in the‘**malloc()**-area’, [PL_get_pointer()](foreigninclude.html#PL_get_pointer()) will get the pointer back.

`bool` **PL_put_float**(`term_t -t, double f`)  
Put a floating-point value in the term reference.

`bool` **PL_put_functor**(`term_t -t, functor_t functor`)  
Create a new compound term from `functor` and bind `t` to this term. All arguments of the term will be variables. To create a term with instantiated arguments, either instantiate the arguments using the `PL_unify_*()` functions or use [PL_cons_functor()](foreigninclude.html#PL_cons_functor()).

`bool` **PL_put_list**(`term_t -l`)  
As [PL_put_functor()](foreigninclude.html#PL_put_functor()), using the list-cell functor. Note that on classical Prolog systems or in SWI-Prolog using the option **--traditional**, this is `.``/2`, while on SWI-Prolog version 7 this is `[|]``/2`.

`bool` **PL_put_nil**(`term_t -l`)  
Put the list terminator constant in `l`. Always returns `TRUE`. Note that in classical Prolog systems or in SWI-Prolog using the option **--traditional**, this is the same as [`PL_put_atom_chars("[]")`](foreigninclude.html#PL_put_atom_chars()). See [section 5.1](ext-lists.html#sec:5.1).

`bool` **PL_put_term**(`term_t -t1, term_t +t2`)  
Make `t1` point to the same term as `t2`. Under the unusual condition that `t2` is a fresh term reference this function requires a global stack cell and may thus return `FALSE` and leave a resource exception in the environment.

`bool` **PL_cons_functor**(`term_t -h, functor_t f, ...`)  
Create a term whose arguments are filled from a variable argument list holding the same number of `term_t` objects as the arity of the functor. To create the term `animal(gnu, 50)`, use:

``` code
{ term_t a1 = PL_new_term_ref();
  term_t a2 = PL_new_term_ref();
  term_t t  = PL_new_term_ref();
  functor_t animal2;

  /* animal2 is a constant that may be bound to a global
     variable and re-used
  */
  animal2 = PL_new_functor(PL_new_atom("animal"), 2);

  PL_put_atom_chars(a1, "gnu");
  PL_put_integer(a2, 50);
  PL_cons_functor(t, animal2, a1, a2);
}
```

After this sequence, the term references `a1` and `a2` may be used for other purposes.

`bool` **PL_cons_functor_v**(`term_t -h, functor_t f, term_t a0`)  
Create a compound term like [PL_cons_functor()](foreigninclude.html#PL_cons_functor()), but `a0` is an array of term references as returned by [PL_new_term_refs()](foreigntypes.html#PL_new_term_refs()). The length of this array should match the number of arguments required by the functor.

`bool` **PL_cons_list**(`term_t -l, term_t +h, term_t +t`)  
Create a list (cons-) cell in `l` from the head `h` and tail `t`. As with [PL_cons_functor()](foreigninclude.html#PL_cons_functor()), the term references `h` and `t` may be used for other purposes after the call to [PL_cons_list()](foreigninclude.html#PL_cons_list()). The code below creates a list of atoms from a `char **`. The list is built tail-to-head. The `PL_unify_*()` functions can be used instead to build a list head-to-tail.

``` code
void
put_list(term_t l, int n, char **words)
{ term_t a = PL_new_term_ref();

  PL_put_nil(l);
  while( --n >= 0 )
  { PL_put_atom_chars(a, words[n]);
    PL_cons_list(l, a, l);
  }
}
```

`bool` **PL_put_dict**(`term_t -h, atom_t tag, size_t len, const atom_t *keys, term_t values`)  
Create a dict from a `tag` and vector of atom-value pairs and put the result in `h`. The dict's key is set by `tag`, which may be `0` to leave the tag unbound. The `keys` vector is a vector of atoms of at least `len` long. The `values` is a term vector allocated using [PL_new_term_refs()](foreigntypes.html#PL_new_term_refs()) of at least `len` long. This function returns `true` on success. On failure it returns `false` and leaves an exception in the environment. This is either a resource exception or a `duplicate_key(Key)` exception.^(221Versions up to 10.1.2 returned -1 if the `tag` or one of the `keys` was invalid and -2 if there are duplicate keys. Later versions terminate with a fatal ABI error if `tag` or one of the `keys` is invalid and raises a normal Prolog exception on duplicate keys.)

### 12.4.6 Unifying data

The functions of this section *unify* terms with other terms or translated C data structures. Except for [PL_unify()](foreigninclude.html#PL_unify()), these functions are specific to SWI-Prolog. They have been introduced because they shorten the code for returning data to Prolog and at the same time make this more efficient by avoiding the need to allocate temporary term references and reduce the number of calls to the Prolog API. Consider the case where we want a foreign function to return the host name of the machine Prolog is running on. Using the `PL_get_*()` and `PL_put_*()` functions, the code becomes:

``` code
foreign_t
pl_hostname(term_t name)
{ char buf[100];

  if ( gethostname(buf, sizeof buf) )
  { term_t tmp = PL_new_term_ref();

    PL_put_atom_chars(tmp, buf);
    return PL_unify(name, tmp);
  }

  PL_fail;
}
```

Using [PL_unify_atom_chars()](foreigninclude.html#PL_unify_atom_chars()), this becomes:

``` code
foreign_t
pl_hostname(term_t name)
{ char buf[100];

  if ( gethostname(buf, sizeof buf) )
    return PL_unify_atom_chars(name, buf);

  PL_fail;
}
```

Note that unification functions that perform multiple bindings may leave part of the bindings in case of failure. See [PL_unify()](foreigninclude.html#PL_unify()) for details.

`bool` **PL_unify**(`term_t ?t1, term_t ?t2`)  
Unify two Prolog terms and return `TRUE` on success. [PL_unify()](foreigninclude.html#PL_unify()) does not evaluate *attributed variables* (see [section 8.1](attvar.html#sec:8.1)), it merely schedules the goals associated with the attributes to be executed *after* the foreign predicate succeeds.^(222Goal associated with attributes may be non-deterministic, which we cannot handle from a callback. A callback could also result in deeply nested mutual recursion between C and Prolog and eventually trigger a C stack overflow.)

Care is needed if [PL_unify()](foreigninclude.html#PL_unify()) returns `FALSE` and the foreign function does not *immediately* return to Prolog with `FALSE`. Unification may perform multiple changes to either `t1` or `t2`. A failing unification may have created bindings before failure is detected. *Already created bindings are not undone*. For example, calling [PL_unify()](foreigninclude.html#PL_unify()) on `a(X, a)` and `a(c,b)` binds `X` to `c` and fails when trying to unify `a` to `b`. If control remains in C or if we want to return success to Prolog, we *must* undo such bindings. In addition, [PL_unify()](foreigninclude.html#PL_unify()) may have failed on an **exception**, typically a resource (stack) overflow. This can be tested using [PL_exception()](foreigninclude.html#PL_exception()), passing 0 (zero) for the query-id argument. Foreign functions that encounter an exception must return `FALSE` to Prolog as soon as possible or call [PL_clear_exception()](foreigninclude.html#PL_clear_exception()) if they wish to ignore the exception. Note that there can only be an exception if [PL_unify()](foreigninclude.html#PL_unify()) returned `FALSE`.

In some scenarios we need to undo *partial unifications*. Suppose we have a database that contains Prolog terms and we run a query over this database. We must succeed on the first successful unification. If a unification is not successful, we must stop if there is an exception or undo the partial unification and try again. Suppose our database contains `f(a,1)` and `f(b,2)` and our query is `f(A,2)`. This should succeed with `A` = `b`, but the first unification binds `A` to `a` before failing to unify 1 with 2.

``` code
static foreign_t
find_in_db(term_t target)
{ fid_t fid = PL_open_foreign_frame();
  term_t candidate = PL_new_term_ref();

  while(get_from_my_database(candidate))
  { if ( PL_unify(candidate, target) ) /* found */
    { PL_close_foreign_frame(fid);
      return TRUE;
    } else if ( PL_exception(0) )      /* error */
    { PL_close_foreign_frame(fid);
      return FALSE;
    }

    PL_rewind_foreign_frame(fid);      /* try next */
  }
  PL_close_foreign_frame(fid);         /* not found */
  return FALSE;
}
```

This code is only needed if the foreign predicate does not return immediately to Prolog when [PL_unify()](foreigninclude.html#PL_unify()) fails - there is an implicit frame around the entire predicate, and returning `FALSE` undoes all bindings when that frame is closed.

`bool` **PL_unify_atom**(`term_t ?t, atom_t a`)  
Unify `t` with the atom `a` and return non-zero on success.

`bool` **PL_unify_bool**(`term_t ?t, int a`)  
Unify `t` with either `false` or `true`, according to whether `a` is zero or non-zero. If `t` is instantiated, `off` and `on` are also accepted.

`bool` **PL_unify_chars**(`term_t ?t, int flags, size_t len, const char *chars`)  
New function to deal with unification of `char*` with various encodings to a Prolog representation. The `flags` argument is a bitwise *or* specifying the Prolog target type and the encoding of `chars`. A Prolog type is one of `PL_ATOM`, `PL_STRING`, `PL_CODE_LIST` or `PL_CHAR_LIST`. A representation is one of `REP_ISO_LATIN_1`, `REP_UTF8` or `REP_MB`. See [PL_get_chars()](foreigninclude.html#PL_get_chars()) for a definition of the representation types. If `len` is `-1` `chars` must be zero-terminated and the length is computed from `chars` using **strlen()**.

If `flags` includes `PL_DIFF_LIST` and type is one of `PL_CODE_LIST` or `PL_CHAR_LIST`, the text is converted to a *difference list*. The tail of the difference list is `t+1`.

`bool` **PL_unify_atom_chars**(`term_t ?t, const char *chars`)  
Unify `t` with an atom created from `chars` and return non-zero on success.

`bool` **PL_unify_list_chars**(`term_t ?t, const char *chars`)  
Unify `t` with a list of ASCII characters constructed from `chars`.

`bool` **PL_unify_string_chars**(`term_t ?t, const char *chars`)  
Unify `t` with a Prolog string object created from the zero-terminated string `chars`. The data will be copied. See also [PL_unify_string_nchars()](foreigninclude.html#PL_unify_string_nchars()).

`bool` **PL_unify_integer**(`term_t ?t, intptr_t n`)  
Unify `t` with a Prolog integer from `n`.

`bool` **PL_unify_int64**(`term_t ?t, int64_t n`)  
Unify `t` with a Prolog integer from `n`.

`bool` **PL_unify_uint64**(`term_t ?t, uint64_t n`)  
Unify `t` with a Prolog integer from `n`. Note that unbounded integer support is required if `n` does not fit in a *signed* `int64_t`. If unbounded integers are not supported a `representation_error` is raised.

`bool` **PL_unify_float**(`term_t ?t, double f`)  
Unify `t` with a Prolog float from `f`.

`bool` **PL_unify_pointer**(`term_t ?t, void *ptr`)  
Unify `t` with a Prolog integer describing the pointer. See also [PL_put_pointer()](foreigninclude.html#PL_put_pointer()) and [PL_get_pointer()](foreigninclude.html#PL_get_pointer()).

`bool` **PL_unify_functor**(`term_t ?t, functor_t f`)  
If `t` is a compound term with the given functor, just succeed. If it is unbound, create a term and bind the variable, else fail. Note that this function does not create a term if the argument is already instantiated. If `f` is a functor with arity 0, `t` is unified with an atom. See also [PL_unify_compound()](foreigninclude.html#PL_unify_compound()).

`bool` **PL_unify_compound**(`term_t ?t, functor_t f`)  
If `t` is a compound term with the given functor, just succeed. If it is unbound, create a term and bind the variable, else fail. Note that this function does not create a term if the argument is already instantiated. If `f` is a functor with arity 0, `t` is unified with compound without arguments. See also [PL_unify_functor()](foreigninclude.html#PL_unify_functor()).

`bool` **PL_unify_list**(`term_t ?l, term_t -h, term_t -t`)  
Unify `l` with a list-cell (`./2`). If successful, write a reference to the head of the list into `h` and a reference to the tail of the list into `t`. This reference to `h` may be used for subsequent calls to this function. Suppose we want to return a list of atoms from a `char **`. We could use the example described by [PL_cons_list()](foreigninclude.html#PL_cons_list()), followed by a call to [PL_unify()](foreigninclude.html#PL_unify()), or we can use the code below. If the predicate argument is unbound, the difference is minimal (the code based on [PL_cons_list()](foreigninclude.html#PL_cons_list()) is probably slightly faster). If the argument is bound, the code below may fail before reaching the end of the word list, but even if the unification succeeds, this code avoids a duplicate (garbage) list and a deep unification.

Note that [PL_unify_list()](foreigninclude.html#PL_unify_list()) is not used with `env` but with `tail`, which is a copy of `env`. [PL_copy_term_ref()](foreigntypes.html#PL_copy_term_ref()) creates a copy `term_t` holding the same Prolog term, i.e., *not* a copy of the Prolog term. The only thing that is allowed to be done with an argument to a foreign predicate (such as `env`) is unification; for anything that might over-write the term, you must use a copy created by [PL_copy_term_ref()](foreigntypes.html#PL_copy_term_ref()). The name [PL_unify_list()](foreigninclude.html#PL_unify_list()) is slightly misleading - it unifies the first argument (`l` but *overwrites* the second (`h`) and third (`t`) arguments.

``` code
foreign_t
pl_get_environ(term_t env)
{ term_t tail = PL_copy_term_ref(env);
  term_t item = PL_new_term_ref();
  extern char **environ;

  for(const char **e = environ; *e; e++)
  { if ( !PL_unify_list(tail, item, tail) ||
         !PL_unify_atom_chars(item, *e) )
      PL_fail;
  }

  return PL_unify_nil(tail);
}
```

In this example, `item` is initialized outside the loop. This allocates a single new reference to a term, which is used as a temporary inside the loop - there is no need to allocate a new reference each time around the loop because the `item` term reference can be reused and the call to [PL_unify_list()](foreigninclude.html#PL_unify_list()) copies a reference to the new list cell's head into the the term referenced by `item`.

`bool` **PL_unify_nil**(`term_t ?l`)  
Unify `l` with the atom `[]`.

`bool` **PL_unify_arg**(`int index, term_t ?t, term_t ?a`)  
Unifies the *index-th* argument (1-based) of `t` with `a`.

`bool` **PL_unify_term**(`term_t ?t, ...`)  
Unify `t` with a (normally) compound term. The remaining arguments are a sequence of a type identifier followed by the required arguments. This predicate is an extension to the Quintus and SICStus foreign interface from which the SWI-Prolog foreign interface has been derived, but has proved to be a powerful and comfortable way to create compound terms from C. Due to the vararg packing/unpacking and the required type-switching this interface is slightly slower than using the primitives. Please note that some bad C compilers have fairly low limits on the number of arguments that may be passed to a function.

Special attention is required when passing numbers. C‘promotes’any integral smaller than `int` to `int`. That is, the types `char`, `short` and `int` are all passed as `int`. In addition, on most 32-bit platforms `int` and `long` are the same. Up to version 4.0.5, only `PL_INTEGER` could be specified, which was taken from the stack as `long`. Such code fails when passing small integral types on machines where `int` is smaller than `long`. It is advised to use `PL_SHORT`, `PL_INT` or `PL_LONG` as appropriate. Similarly, C compilers promote `float` to `double` and therefore `PL_FLOAT` and `PL_DOUBLE` are synonyms.

The type identifiers are:

**`PL_VARIABLE` `none`**  
No op. Used in arguments of `PL_FUNCTOR`.

**`PL_BOOL` `int`**  
Unify the argument with `true` or `false`.

**`PL_ATOM` `atom_t`**  
Unify the argument with an atom, as in [PL_unify_atom()](foreigninclude.html#PL_unify_atom()).

**`PL_CHARS` `const char *`**  
Unify the argument with an atom constructed from the C `char *`, as in [PL_unify_atom_chars()](foreigninclude.html#PL_unify_atom_chars()).

**`PL_NCHARS` `size_t, const char *`**  
Unify the argument with an atom constructed from length and `char*` as in [PL_unify_atom_nchars()](foreigninclude.html#PL_unify_atom_nchars()).

**`PL_UTF8_CHARS` `const char *`**  
Create an atom from a UTF-8 string.

**`PL_UTF8_STRING` `const char *`**  
Create a packed string object from a UTF-8 string.

**`PL_MBCHARS` `const char *`**  
Create an atom from a multi-byte string in the current locale.

**`PL_MBCODES` `const char *`**  
Create a list of character codes from a multi-byte string in the current locale.

**`PL_MBSTRING` `const char *`**  
Create a packed string object from a multi-byte string in the current locale.

**`PL_NWCHARS` `size_t, const wchar_t *`**  
Create an atom from a length and a wide character pointer.

**`PL_NWCODES` `size_t, const wchar_t *`**  
Create a list of character codes from a length and a wide character pointer.

**`PL_NWSTRING` `size_t, const wchar_t *`**  
Create a packed string object from a length and a wide character pointer.

**`PL_SHORT` `short`**  
Unify the argument with an integer, as in [PL_unify_integer()](foreigninclude.html#PL_unify_integer()). As `short` is promoted to `int`, `PL_SHORT` is a synonym for `PL_INT`.

**`PL_INTEGER` `long`**  
Unify the argument with an integer, as in [PL_unify_integer()](foreigninclude.html#PL_unify_integer()).

**`PL_INT` `int`**  
Unify the argument with an integer, as in [PL_unify_integer()](foreigninclude.html#PL_unify_integer()).

**`PL_LONG` `long`**  
Unify the argument with an integer, as in [PL_unify_integer()](foreigninclude.html#PL_unify_integer()).

**`PL_INT64` `int64_t`**  
Unify the argument with a 64-bit integer, as in [PL_unify_int64()](foreigninclude.html#PL_unify_int64()).

**`PL_INTPTR` `intptr_t`**  
Unify the argument with an integer with the same width as a pointer. On most machines this is the same as `PL_LONG`. but on 64-bit MS-Windows pointers are 64 bits while longs are only 32 bits.

**`PL_DOUBLE` `double`**  
Unify the argument with a float, as in [PL_unify_float()](foreigninclude.html#PL_unify_float()). Note that, as the argument is passed using the C vararg conventions, a float must be casted to a double explicitly.

**`PL_FLOAT` `double`**  
Unify the argument with a float, as in [PL_unify_float()](foreigninclude.html#PL_unify_float()).

**`PL_POINTER` `void *`**  
Unify the argument with a pointer, as in [PL_unify_pointer()](foreigninclude.html#PL_unify_pointer()).

**`PL_STRING` `const char *`**  
Unify the argument with a string object, as in [PL_unify_string_chars()](foreigninclude.html#PL_unify_string_chars()).

**`PL_TERM` `term_t`**  
Unify a subterm. Note this may be the return value of a [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()) call to get access to a variable.

**`PL_FUNCTOR` `functor_t, ...`**  
Unify the argument with a compound term. This specification should be followed by exactly as many specifications as the number of arguments of the compound term.

**`PL_FUNCTOR_CHARS` `const char *name, int arity, ...`**  
Create a functor from the given name and arity and then behave as `PL_FUNCTOR`.

**`PL_LIST` `int length, ...`**  
Create a list of the indicated length. The remaining arguments contain the elements of the list.

For example, to unify an argument with the term `language(dutch)`, the following skeleton may be used:

``` code
static functor_t FUNCTOR_language1;

static void
init_constants()
{ FUNCTOR_language1 = PL_new_functor(PL_new_atom("language"),1);
}

foreign_t
pl_get_lang(term_t r)
{ return PL_unify_term(r,
                       PL_FUNCTOR, FUNCTOR_language1,
                           PL_CHARS, "dutch");
}

install_t
install()
{ PL_register_foreign("get_lang", 1, pl_get_lang, 0);
  init_constants();
}
```

`bool` **PL_chars_to_term**(`const char *chars, term_t -t`)  
Parse the string `chars` and put the resulting Prolog term into `t`. `chars` may or may not be closed using a Prolog full-stop (i.e., a dot followed by a blank). Returns `FALSE` if a syntax error was encountered and `TRUE` after successful completion. In addition to returning `FALSE`, the exception-term is returned in `t` on a syntax error. See also [term_to_atom/2](manipatom.html#term_to_atom/2).

The following example builds a goal term from a string and calls it.

``` code
int
call_chars(const char *goal)
{ fid_t fid = PL_open_foreign_frame();
  term_t g = PL_new_term_ref();
  BOOL rval;

  if ( PL_chars_to_term(goal, g) )
    rval = PL_call(goal, NULL);
  else
    rval = FALSE;

  PL_discard_foreign_frame(fid);
  return rval;
}
  ...
  call_chars("consult(load)");
  ...
```

[PL_chars_to_term()](foreigninclude.html#PL_chars_to_term()) is defined using [PL_put_term_from_chars()](foreigninclude.html#PL_put_term_from_chars()) which can deal with not null-terminated strings as well as strings using different encodings:

``` code
int
PL_chars_to_term(const char *s, term_t t)
{ return PL_put_term_from_chars(t, REP_ISO_LATIN_1, (size_t)-1, s);
}
```

`bool` **PL_wchars_to_term**(`const pl_wchar_t *chars, term_t -t`)  
Wide character version of [PL_chars_to_term()](foreigninclude.html#PL_chars_to_term()).

`char *` **PL_quote**(`int chr, const char *string`)  
Return a quoted version of `string`. If `chr` is `'\''`, the result is a quoted atom. If `chr` is `'"'`, the result is a string. The result string is stored in the same ring of buffers as described with the `BUF_STACK` argument of [PL_get_chars()](foreigninclude.html#PL_get_chars());

In the current implementation, the string is surrounded by `chr` and any occurrence of `chr` is doubled. In the future the behaviour will depend on the [character_escapes](flags.html#flag:character_escapes) Prolog flag.

`int` **PL_for_dict**(`term_t dict, int (*func)(term_t key, term_t value, void *closure), void *closure, int flags`)  
Iterates over `dict`, calling `func` for each item. In each call, `key` and `value` are the processed item's key-value pair and the `closure` argument is passed from the call to [PL_for_dict()](foreigninclude.html#PL_for_dict()). If `func` returns a non-`0` value, the iteration stops and [PL_for_dict()](foreigninclude.html#PL_for_dict()) returns that value; otherwise, all pairs are processed and [PL_for_dict()](foreigninclude.html#PL_for_dict()) returns `0`. If `flags` contains `PL_FOR_DICT_SORTED`, the key-value pairs are processed in the standard order of terms; otherwise the processing order is unspecified.

### 12.4.7 Convenient functions to generate Prolog exceptions

The typical implementation of a foreign predicate first uses the PL_get\_\*() functions to extract C data types from the Prolog terms. Failure of any of these functions is normally because the Prolog term is of the wrong type. The \***\_ex()** family of functions are wrappers around (mostly) the PL_get\_\*() functions, such that we can write code in the style below and get proper exceptions if an argument is uninstantiated or of the wrong type. [Section 12.4.8](foreigninclude.html#sec:12.4.8) documents an alternative API to fetch values for the C basic types.

``` code
/** set_size(+Name:atom, +Width:int, +Height:int) is det.

static foreign_t
set_size(term_t name, term_t width, term_t height)
{ char *n;
  int w, h;

  if ( !PL_get_chars(name, &n, CVT_ATOM|CVT_EXCEPTION) ||
       !PL_get_integer_ex(with, &w) ||
       !PL_get_integer_ex(height, &h) )
    return FALSE;

  ...

}
```

`bool` **PL_get_atom_ex**(`term_t t, atom_t *a`)  
As [PL_get_atom()](foreigninclude.html#PL_get_atom()), but raises a type or instantiation error if `t` is not an atom.

`bool` **PL_get_integer_ex**(`term_t t, int *i`)  
As [PL_get_integer()](foreigninclude.html#PL_get_integer()), but raises a type or instantiation error if `t` is not an integer, or a representation error if the Prolog integer does not fit in a C `int`.

`bool` **PL_get_long_ex**(`term_t t, long *i`)  
As [PL_get_long()](foreigninclude.html#PL_get_long()), but raises a type or instantiation error if `t` is not an atom, or a representation error if the Prolog integer does not fit in a C `long`.

`bool` **PL_get_int64_ex**(`term_t t, int64_t *i`)  
As [PL_get_int64()](foreigninclude.html#PL_get_int64()), but raises a type or instantiation error if `t` is not an integer, or a representation error if the Prolog integer does not fit in a C `int64_t`.

`bool` **PL_get_uint64_ex**(`term_t t, uint64_t *i`)  
As [PL_get_uint64()](foreigninclude.html#PL_get_uint64()), but raises a type, domain or instantiation error if `t` is not an integer or `t` is less than zero, or a representation error if the Prolog integer does not fit in a C `int64_t`.

`bool` **PL_get_intptr_ex**(`term_t t, intptr_t *i`)  
As [PL_get_intptr()](foreigninclude.html#PL_get_intptr()), but raises a type or instantiation error if `t` is not an atom, or a representation error if the Prolog integer does not fit in a C `intptr_t`.

`bool` **PL_get_size_ex**(`term_t t, size_t *i`)  
As [PL_get_intptr()](foreigninclude.html#PL_get_intptr()), but raises a type or instantiation error if `t` is not an integer, or a representation error if the Prolog integer does not fit in a C `size_t`.

`bool` **PL_get_bool_ex**(`term_t t, int *i`)  
As [PL_get_bool()](foreigninclude.html#PL_get_bool()), but raises a type or instantiation error if `t` is not a valid boolean value (`true`, `false`, `on`, constoff, `1` or `0`). Note that the pointer is to an `int` because C has no `bool` type.

`bool` **PL_get_float_ex**(`term_t t, double *f`)  
As [PL_get_float()](foreigninclude.html#PL_get_float()), but raises a type or instantiation error if `t` is not a float.

`bool` **PL_get_char_ex**(`term_t t, int *p, int eof`)  
Get a character code from `t`, where `t` is either an integer or an atom with length one. If `eof` is `TRUE` and `t` is -1, `p` is filled with -1. Raises an appropriate error if the conversion is not possible.

`bool` **PL_get_pointer_ex**(`term_t t, void **addrp`)  
As [PL_get_pointer()](foreigninclude.html#PL_get_pointer()), but raises a type or instantiation error if `t` is not a pointer.

`bool` **PL_get_list_ex**(`term_t l, term_t h, term_t t`)  
As [PL_get_list()](foreigninclude.html#PL_get_list()), but raises a type or instantiation error if `t` is not a list.

`bool` **PL_get_nil_ex**(`term_t l`)  
As [PL_get_nil()](foreigninclude.html#PL_get_nil()), but raises a type or instantiation error if `t` is not the empty list. Because [PL_get_nil_ex()](foreigninclude.html#PL_get_nil_ex()) is commonly used after a `while` loop over [PL_get_list_ex()](foreigninclude.html#PL_get_list_ex()), it fails immediately if there is an exception pending (from [PL_get_list_ex()](foreigninclude.html#PL_get_list_ex())).

`bool` **PL_unify_list_ex**(`term_t l, term_t h, term_t t`)  
As [PL_unify_list()](foreigninclude.html#PL_unify_list()), but raises a type error if `t` is not a variable, list-cell or the empty list.

`bool` **PL_unify_nil_ex**(`term_t l`)  
As [PL_unify_nil()](foreigninclude.html#PL_unify_nil()), but raises a type error if `t` is not a variable, list-cell or the empty list.

`bool` **PL_unify_bool_ex**(`term_t t, int val`)  
As [PL_unify_bool()](foreigninclude.html#PL_unify_bool()), but raises a type error if `t` is not a variable or a boolean.

The second family of functions in this section simplifies the generation of ISO compatible error terms. Any foreign function that calls this function must return to Prolog with the return code of the error function or the constant `FALSE`. If available, these error functions add the name of the calling predicate to the error context. See also [PL_raise_exception()](foreigninclude.html#PL_raise_exception()).

`bool` **PL_instantiation_error**(`term_t culprit`)  
Raise `instantiation_error`. `Culprit` is ignored, but should be bound to the term that is insufficiently instantiated. See [instantiation_error/1](error.html#instantiation_error/1).

`bool` **PL_uninstantiation_error**(`term_t culprit`)  
Raise `uninstantiation_error(culprit)`. This should be called if an argument that must be unbound at entry is bound to `culprit`. This error is typically raised for a pure output arguments such as a newly created stream handle (e.g., the third argument of [open/3](IO.html#open/3)).

`bool` **PL_representation_error**(`const char *resource`)  
Raise `representation_error(resource)`. See [representation_error/1](error.html#representation_error/1).

`bool` **PL_type_error**(`const char *expected, term_t culprit`)  
Raise `type_error(expected, culprit)`. See [type_error/2](error.html#type_error/2).

`bool` **PL_domain_error**(`const char *expected, term_t culprit`)  
Raise `domain_error(expected, culprit)`. See [domain_error/2](error.html#domain_error/2).

`bool` **PL_existence_error**(`const char *type, term_t culprit`)  
Raise `existence_error(type, culprit)`. See [existence_error/2](error.html#existence_error/2).

`bool` **PL_permission_error**(`const char *operation, const char *type, term_t culprit`)  
Raise `permission_error(operation, type, culprit)`. See [permission_error/3](error.html#permission_error/3).

`bool` **PL_resource_error**(`const char *resource`)  
Raise `resource_error(resource)`. See [resource_error/1](error.html#resource_error/1).

`bool` **PL_syntax_error**(`const char *message, IOSTREAM *in`)  
Raise `syntax_error(message)`. If `in` is not `NULL`, add information about the current position of the input stream.

### 12.4.8 Foreign language wrapper support functions

In addition to the functions described in [section 12.4.4.2](foreigninclude.html#sec:12.4.4.2), there is a family of functions that is used for automatic generation of wrapper functions, for example using the Prolog library `library(qpforeign)` that provides a Quintus/SICStus compatible foreign language interface.

The PL_cvt_i\_\*() family of functions is suitable for use with a `_Generic` selector or C++ overloading.^(223`_Generic` needs to take into account that there's no `bool` type in C but there is in C++. An overloaded **integer()** method is provided in the C++ interface.)

Note that the documentation on this API is incomplete. Also note that many of these functions are equivalent to the PL_get\_\***\_ex()** functions described in [section 12.4.7](foreigninclude.html#sec:12.4.7).

`bool` **PL_cvt_i_bool**(`term_t p, int *c`)  
Equivalent to [PL_get_bool_ex()](foreigninclude.html#PL_get_bool_ex()). Note that the pointer is to an `int` because C has no `bool` type. The return value is either `0` or `1`.

`bool` **PL_cvt_i_char**(`term_t p, char *c`)  
`bool` **PL_cvt_i_schar**(`term_t p, signed char *c`)  
`bool` **PL_cvt_i_uchar**(`term_t p, unsigned char *c`)  
`bool` **PL_cvt_i_short**(`term_t p, short *s`)  
`bool` **PL_cvt_i_ushort**(`term_t p, unsigned short *s`)  
`bool` **PL_cvt_i_int**(`term_t p, int *c`)  
`bool` **PL_cvt_i_uint**(`term_t p, unsigned int *c`)  
`bool` **PL_cvt_i_long**(`term_t p, long *c`)  
`bool` **PL_cvt_i_ulong**(`term_t p, unsigned long *c`)  
`bool` **PL_cvt_i_llong**(`term_t p, long long *c`)  
`bool` **PL_cvt_i_ullong**(`term_t p, unsigned long long *c`)  
`bool` **PL_cvt_i_int32**(`term_t p, int32_t *c`)  
`bool` **PL_cvt_i_uint32**(`term_t p, uint32_t *c`)  
`bool` **PL_cvt_i_int64**(`term_t p, int64_t *c`)  
`bool` **PL_cvt_i_uint64**(`term_t p, uint64_t *c`)  
`bool` **PL_cvt_i_size_t**(`term_t p, size_t *c`)  
Convert a Prolog integer into a C integer of the specified size. Generate an exception and return `FALSE` if the conversion is impossible because the Prolog term is not an integer or the C type cannot represent the value of the Prolog integer.

### 12.4.9 Serializing and deserializing Prolog terms

`bool` **PL_put_term_from_chars**(`term_t t, int flags, size_t len, const char *s`)  
Parse the text from the C-string `s` holding `len` bytes and put the resulting term in `t`. `len` can be `(size_t)-1`, assuming a 0-terminated string. The `flags` argument controls the encoding and is currently one of `REP_UTF8` (string is UTF8 encoded), `REP_MB` (string is encoded in the current locale) or 0 (string is encoded in ISO latin 1). The string may, but is not required, to be closed by a full stop (.).

If parsing produces an exception; the behaviour depends on the `CVT_EXCEPTION` flag. If present, the exception is propagated into the environment. Otherwise, the exception is placed in `t` and the return value is `FALSE`.^(224The `CVT_EXCEPTION` was added in version 8.3.12).

### 12.4.10 BLOBS: Using atoms to store arbitrary binary data

SWI-Prolog atoms as well as strings can represent arbitrary binary data of arbitrary length. This facility is attractive for storing foreign data such as images in an atom. An atom is a unique handle to this data and the atom garbage collector is able to destroy atoms that are no longer referenced by the Prolog engine. This property of atoms makes them attractive as a handle to foreign resources, such as Java atoms, Microsoft's COM objects, etc., providing safe combined garbage collection.

To exploit these features safely and in an organised manner, the SWI-Prolog foreign interface allows creating‘atoms’with additional type information. The type is represented by a structure holding C function pointers that tell Prolog how to handle releasing the atom, writing it, sorting it, etc. Two atoms created with different types can represent the same sequence of bytes. Atoms are first ordered on the rank number of the type and then on the result of the [compare()](foreigninclude.html#compare()) function. Rank numbers are assigned when the type is registered. This implies that the results of inequality comparisons between blobs of different types is undefined and can change if the program is run twice (the ordering within a blob type will not change, of course).

While the blob is alive, neither its handle nor the location of the contents (see [PL_blob_data()](foreigninclude.html#PL_blob_data())) change. If the blob's type has the `PL_BLOB_UNIQUE` feature, the content of the blob must remain unmodified. If the blob's type does not have the `PL_BLOB_UNIQUE` feature multiple instances of this blob type may contain the same data. The blob *handle* (`atom_t`) is reclaimed *only* by the atom garbage collector. The blob's *content* (data) is normally reclaimed when the garbage collector reclaims the blob. If the blob's type defines the [release()](foreigninclude.html#release()) function, this function is called. This hook may deal with side effects and is responsible of releasing the data if the blob's type has the `PL_BLOB_NOCOPY` flag. The content of a `PL_BLOB_NOCOPY` blob may be released before the blob itself can be garbage collected using [PL_free_blob()](foreigninclude.html#PL_free_blob()). This immediately triggers the [release()](foreigninclude.html#release()) function. After [PL_free_blob()](foreigninclude.html#PL_free_blob()) has reclaimed the content, this function will not be called when the `atom_t` handle is reclaimed. An `atom_t` handle may be reused for a new atom or blob after it has been garbage collected.

If foreign code stores the `atom_t` handle in some permanent location it must make sure the handle is *registered* to prevent it from being garbage collected. If the handle is obtained from a `term_t` object it is **not** registered because it is protected by the `term_t` object. This applies to e.g., [PL_get_atom()](foreigninclude.html#PL_get_atom()). Functions that create a handle from data, such as [PL_new_atom()](foreigninclude.html#PL_new_atom()), return a registered handle to prevent the asynchronous atom garbage collector from reclaiming it immediately. Note that many of the API functions create an atom or blob handle and use this to fill a `term_t` object, e.g., [PL_unify_blob()](foreigninclude.html#PL_unify_blob()), [PL_unify_chars()](foreigninclude.html#PL_unify_chars()), etc. In this scenario the handle is protected by the `term_t` object. Registering and unregistering `atom_t` handles is done by [PL_register_atom()](foreigninclude.html#PL_register_atom()) and [PL_unregister_atom()](foreigninclude.html#PL_unregister_atom()).

Note that during program shutdown using [PL_cleanup()](foreigninclude.html#PL_cleanup()), all atoms and blobs are reclaimed as described above. **These objects are reclaimed regardless of their registration count. The order in which the atoms or blobs are reclaimed under [PL_cleanup()](foreigninclude.html#PL_cleanup()) is undefined.** However, when these objects are reclaimed using [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0), registration counts are taken into account.

#### 12.4.10.1 Defining a BLOB type

The type `PL_blob_t` represents a structure with the layout displayed below. The structure contains additional fields at the ... for internal bookkeeping as well as future extensions.

``` code
typedef struct PL_blob_t
{ uintptr_t     magic;          /* PL_BLOB_MAGIC */
  uintptr_t     flags;          /* Bitwise or of PL_BLOB_* */
  const char *  name;           /* name of the type */
  int           (*release)(atom_t a);
  int           (*compare)(atom_t a, atom_t b);
  int           (*write)(IOSTREAM *s, atom_t a, int flags);
  void          (*acquire)(atom_t a);
  int           (*save)(atom_t a, IOSTREAM *s);
  atom_t        (*load)(IOSTREAM *s);
  ...
} PL_blob_t;
```

For each type, exactly one such structure should be allocated and must not be moved because the address of the structure determines the blob's "type". Its first field must be initialised to `PL_BLOB_MAGIC`. If a blob type is registered from a loadable object (shared object or DLL) the blob type must be deregistered using [PL_unregister_blob_type()](foreigninclude.html#PL_unregister_blob_type()) before the object may be released.

The `flags` is a bitwise *or* of the following constants:

**PL_BLOB_TEXT**  
If specified, the blob is assumed to contain text and is considered a normal Prolog atom. The (currently) two predefined blob types that represent atoms have this flag set. User-defined blobs may not specify this, even if they contain only text. Applications should *not* use the blob API to create normal text atoms or get access to the text represented by normal text atoms. Most applications should use [PL_get_nchars()](foreigninclude.html#PL_get_nchars()) and [PL_unify_chars()](foreigninclude.html#PL_unify_chars()) to get text from Prolog terms or create Prolog terms that represent text.

**PL_BLOB_UNIQUE**  
If specified the system ensures that the blob-handle is a unique reference for a blob with the given type, length and content. If this flag is not specified, each lookup creates a new blob. Uniqueness is determined by comparing the bytes in the blobs unless `PL_BLOB_NOCOPY` is also specified, in which case the pointers are compared. Note that the lookup does *not* use the blob's compare function when testing for equality, but only tests the bytes; this means that terms from the recorded database or C++-style strings will typically not compare as equal when doing blob lookup.

**PL_BLOB_NOCOPY**  
By default the content of the blob is copied. Using this flag the blob references the external data directly. The user must ensure the provided pointer is valid as long as the atom lives. If `PL_BLOB_UNIQUE` is also specified, uniqueness is determined by comparing the pointer rather than the data pointed at. Using `PL_BLOB_UNIQUE``|``PL_BLOB_NOCOPY` can be used to make a blob reference an arbitrary pointer where the pointer data may be reclaimed in the [release()](foreigninclude.html#release()) handler.

**PL_BLOB_WCHAR**  
If `PL_BLOB_TEXT` is also set, then the text is made up of `pl_wchar_t` items and the blob's lenght is the number of bytes (that is, the number of characters times `sizeof(pl_wchar_t)`). As `PL_BLOB_TEXT`, this flag should not be set in user-defined blobs.

The `name` field represents the type name as available to Prolog. See also [current_blob/2](examineprog.html#current_blob/2). The other fields are function pointers that must be initialised to proper functions or `NULL` to get the default behaviour of built-in atoms. Below are the defined member functions:

`void` **acquire**(`atom_t a`)  
Called if a new blob of this type is created through [PL_put_blob()](foreigninclude.html#PL_put_blob()), [PL_unify_blob()](foreigninclude.html#PL_unify_blob()), or [PL_new_blob()](foreigninclude.html#PL_new_blob()). Note this this call is done as part of creating the blob. The call to [PL_unify_blob()](foreigninclude.html#PL_unify_blob()) may fail if the unification fails or cannot be completed due to a resource error. [PL_put_blob()](foreigninclude.html#PL_put_blob()) has no such error conditions. This callback is typically used to store the `atom_t` handle into the content of the blob. Given a pointer to the content, we can now use [PL_unify_atom()](foreigninclude.html#PL_unify_atom()) to bind a Prolog term with a reference to the pointed to object. If the content of the blob can be modified (`PL_BLOB_UNIQUE` is not present) this is the only way to get access to the `atom_t` handle that belongs to this blob. If `PL_BLOB_UNIQUE` is provided and respected, [PL_unify_blob()](foreigninclude.html#PL_unify_blob()) given the same pointer and length will produce the same `atom_t` handle.

`int` **release**(`atom_t a`)  
The blob (atom) `a` is about to be released. The [release()](foreigninclude.html#release()) function is called when the atom is reclaimed by the atom garbage collector, when an explicit call to [PL_free_blob()](foreigninclude.html#PL_free_blob()) is made or during shutdown of Prolog. This function can retrieve the data of the blob using [PL_blob_data()](foreigninclude.html#PL_blob_data()). If the [release()](foreigninclude.html#release()) function returns `FALSE`, the atom garbage collector will *not* reclaim the atom. For critical resources such as file handles or significant memory resources, it may be desirable to have an explicit call to dispose (most of) the resources. For example, [close/1](IO.html#close/1) reclaims the file handle and most of the resources associated with a stream, leaving only a tiny bit of content to the garbage collector. See also [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3).

The [release()](foreigninclude.html#release()) callback is called in the context of the thread executing the atom garbage collect, the thread executing [PL_free_blob()](foreigninclude.html#PL_free_blob()) or the thread initiating the shutdown. Normally the thread `gc` runs all atom and clause garbage collections. The [release()](foreigninclude.html#release()) function may not call any of the PL\_\*() functions except for [PL_blob_data()](foreigninclude.html#PL_blob_data()) or [PL_unregister_atom()](foreigninclude.html#PL_unregister_atom()) to unregister other atoms that are part data associated to the blob. Calling any of the other PL\_\* functions may result in deadlocks or crashes. The [release()](foreigninclude.html#release()) function should not call any potentially slow or blocking functions as this may cause serious slowdowns in the rest of the system.

Blobs that require cleanup that is slow, blocking or requires calling Prolog must pass the data to be cleaned to another thread. Be aware that if the blob uses `PL_BLOB_NOCOPY` the user is responsible for discarding the data, otherwise the atom garbage collector will free the data.

As SWI-Prolog atom garbage collector is *conservative*, there is no guarantee that the [release()](foreigninclude.html#release()) function will ever be called. If it is important to clean up some resource, there should be an explicit predicate for doing that, and calling that predicate should be guaranteed by using [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) or some a process finalization hook such as [at_halt/1](consulting.html#at_halt/1).

Normally, Prolog does not clean memory during shutdown. It does so on an explicit call to [PL_cleanup()](foreigninclude.html#PL_cleanup()).^(225Or if the system is compiled with the **cmake** *build type* `Debug`.) In such a situation, there is no guarantee of the order in which atoms are released; if a blob contains an atom (or another blob), those atoms (or blobs) may have already been released. See also [PL_blob_data()](foreigninclude.html#PL_blob_data()).

`int` **compare**(`atom_t a, atom_t b`)  
Compare the blobs `a` and `b`, both of which are of the type associated to this blob type. Return values are as **memcmp()**: `< 0` if `a` is less than `b`, `= 0` if both are equal, and `> 0` otherwise. The default implementation is a bitwise comparison of the blobs’contents. This default implementation suffices if `PL_BLOB_UNIQUE` is set and the blob follows the requirement that its contents do not change, although it might give an unexpected ordering, and the ordering may change if the blob is saved and restored using save_program/2.

If the [compare()](foreigninclude.html#compare()) function is defined, the [sort/2](builtinlist.html#sort/2) predicate uses that to determine if two blobs are equal and only keeps one of them. This can cause unexpected results with blobs that are actually different; if you cannot guarantee that the blobs all have unique contents, then you should incorporate the blob address (the system guarantees that blobs are not shifted in memory after they are allocated). This function should not call any PL\_\*() functions other than [PL_blob_data()](foreigninclude.html#PL_blob_data()).

The following minimal compare function gives a stable total ordering:

``` code
static int
compare_my_blob(atom_t a, atom_t b)
{ const struct my_blob_data *blob_a = PL_blob_data(a, NULL, NULL);
  const struct my_blob_data *blob_b = PL_blob_data(b, NULL, NULL);
  return (blob_a < blob_b) ? -1 : (blob_a > blob_b) ? 1 : 0;
}
```

`int` **write**(`IOSTREAM *s, atom_t a, int flags`)  
Write the content of the blob `a` to the stream `s` respecting the `flags`. The return value is `TRUE` or `FALSE` and does *not* follow the Unix convention of the number of bytes (where zero is possible) and negative for errors. Any I/O operations to `s` are in the context of a [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()); upon return, the [PL_release_stream()](foreign-streams.html#PL_release_stream()) handles any errors, so it is safe to not check return codes from [Sfprintf()](foreign-streams.html#Sfprintf()), etc.

In general, the output from the [write()](foreigninclude.html#write()) callback should be minimal. If you wish to output more debug information, it is suggested that you either add a debug option to your "open" predicate to output more information, or provide a "properties" predicate. A typical implementation is:

``` code
static int write_my_blob(IOSTREAM *s, atom_t symbol, int flags)
{ (void)flags; /* unused */
  Sfprintf(s, "<my_blob>(%p)", PL_blob_data(symbol, NULL, NULL));
  return TRUE;
}
```

The `flags` are a bitwise *or* of zero or more of the `PL_WRT_*` flags that were passed in to the calling [PL_write_term()](foreign-streams.html#PL_write_term()) that called [write()](foreigninclude.html#write()), and are defined in `SWI-Prolog.h`. The `flags` do not have the `PL_WRT_NEWLINE` bit set, so it is safe to call [PL_write_term()](foreign-streams.html#PL_write_term()) and there is no need for writing a trailing newline. This prototype is available if the `SWI-Stream.h` is included *before* `SWI-Prolog.h`. This function can retrieve the data of the blob using [PL_blob_data()](foreigninclude.html#PL_blob_data()).

Most blobs reference some external data identified by a pointer and the [write()](foreigninclude.html#write()) function writes `<``type``>(``address``)`. If this function is not provided, [write/1](termrw.html#write/1) emits the content of the blob for blobs of type `PL_BLOB_TEXT` or a string of the format `<#`*hex data*`>` for binary blobs.

&nbsp;

`int` **save**(`atom_t a, IOSTREAM *s`)  
Write the blob to stream `s`, in an opaque form that is known only to the blob. If a “save” function is not provided (that is, the field is `NULL`), the default implementation saves and restores the blob as if it is an array of bytes which may contain null (`’`  
`0’`) bytes.

`SWI-Stream.h` defines a number of PL_qlf_put\_\*() functions that write data in a machine-independent form that can be read by the corresponding PL_qlf_get\_\*() functions.

If the “save” function encounters an error, it should call [PL_warning()](foreigninclude.html#PL_warning()), raise an exception (see [PL_raise_exception()](foreigninclude.html#PL_raise_exception())), and return `FALSE`.^(226Details are subject to change.) Note that failure to save/restore a blob makes it impossible to compile a file that contains such a blob using [qcompile/2](consulting.html#qcompile/2) as well as creating a *saved state* from a program that contains such a blob impossible. Here, *contains* means that the blob appears in a clause or directive.

`atom_t` **load**(`IOSTREAM *s`)  
Read the blob from its saved form as written by the “save” function of the same blob type. If this cannot be done (e.g., a stream read failure or a corrupted external form), the “load” function should call [PL_warning()](foreigninclude.html#PL_warning()), then [PL_fatal_error()](foreigninclude.html#PL_fatal_error()), and return constFALSE.^(227Details are subject to change; see the “save” function.) If a “load” function is not provided (that is, the field is `NULL`, the default implementation assumes that the blob was written by the default “save” - that is, as an array of bytes

`SWI-Stream.h` defines a number of PL_qlf_get\_\*() functions that read data in a machine-independent form, as written by the by the corresponding PL_qlf_put\_\*() functions.

The atom that the “load” function returns can be created using [PL_new_blob()](foreigninclude.html#PL_new_blob()).

`bool` **PL_unregister_blob_type**(`PL_blob_t *type`)  
Unlink the blob type from the registered type and transform the type of possible living blobs to `unregistered`, avoiding further reference to the type structure, functions referred by it, as well as the data. This function returns `TRUE` if no blobs of this type existed and `FALSE` otherwise. [PL_unregister_blob_type()](foreigninclude.html#PL_unregister_blob_type()) is intended for the **uninstall()** hook of foreign modules, avoiding further references to the module.

`int` **PL_register_blob_type**(`PL_blob_t *type`)  
This function does not need to be called explicitly. It is called if needed when a blob is created by [PL_unify_blob()](foreigninclude.html#PL_unify_blob()), [PL_put_blob()](foreigninclude.html#PL_put_blob()), or [PL_new_blob()](foreigninclude.html#PL_new_blob()).

#### 12.4.10.2 Accessing blobs

The blob access functions are similar to the atom accessing functions. Blobs being atoms, the atom functions operate on blobs and vice versa. For clarity and possible future compatibility issues, however, it is not advised to rely on this.

`bool` **PL_is_blob**(`term_t t, PL_blob_t **type`)  
Succeeds if `t` refers to a blob, in which case `type` is filled with the type of the blob.

`bool` **PL_unify_blob**(`term_t t, void *blob, size_t len, PL_blob_t *type`)  
Unify `t` to a blob constructed from the given data and associated with the given type. This performs the following steps:

1.  If the `type` has `PL_BLOB_UNIQUE` set, search the blob database for a blob of the same type with the same content. If found, unify `t` with the existing handle.
2.  If not found or `PL_BLOB_UNIQUE` is not set, create a new blob handle. If `PL_BLOB_NOCOPY` is set, associate it to the given memory; else, copy the memory to a new area owned by the blob. Call the [acquire()](foreigninclude.html#acquire()) function of the `type`.
3.  Unify `t` with the existing or new handle. This succeeds if `t` is already bound to the existing blob handle. If `t` is a variable, it succeeds *if* sufficient resources are available to perform the unification; if `t` is bound to something else, this fails.

It is possible that a blob referencing critical resources is created after which the unification fails. Typically these resources are eventually reclaimed because the new blob is not referenced and reclaimed by the atom garbage collector. As described with the [release()](foreigninclude.html#release()) function, it can be desirable to reclaim the critical resources after the failing [PL_unify_blob()](foreigninclude.html#PL_unify_blob()) call.

`bool` **PL_put_blob**(`term_t t, void *blob, size_t len, PL_blob_t *type`)  
Store the described blob in `t`. The return value indicates whether a new blob was allocated (`FALSE`) or the blob is a reference to an existing blob (`TRUE`). Reporting new/existing can be used to deal with external objects having their own reference counts. If the return is `TRUE` this reference count must be incremented, and it must be decremented on blob destruction callback. See also [PL_put_atom_nchars()](foreigninclude.html#PL_put_atom_nchars()).

`atom_t` **PL_new_blob**(`void *blob, size_t len, PL_blob_t *type`)  
Create a blob from its internal opaque form. This function is intended for the “load” function of a blob.

`bool` **PL_get_blob**(`term_t t, void **blob, size_t *len, PL_blob_t **type`)  
If `t` holds a blob or atom, get the data and type and return `TRUE`. Otherwise return `FALSE`. Each result pointer may be `NULL`, in which case the requested information is ignored.

`void *` **PL_blob_data**(`atom_t a, size_t *len, PL_blob_t **type`)  
Get the data and type associated to a blob. This function is mainly used from the callback functions described in [section 12.4.10.1](foreigninclude.html#sec:12.4.10.1). Note that if the [release()](foreigninclude.html#release()) hook is called from [PL_cleanup()](foreigninclude.html#PL_cleanup()), blobs are released regardless of whether or not they are referenced and the order in which blobs are released is undefined (the order depends on the ordering in the atom hash table). [PL_blob_data()](foreigninclude.html#PL_blob_data()) may be called safely on a blob that has already been released. If this happens during [PL_cleanup()](foreigninclude.html#PL_cleanup()) the return value is guaranteed to be `NULL`. During normal execution it may return the content of a newly allocated blob that reuses the released handle.

`bool` **PL_free_blob**(`atom_t blob`)  
New in 9.1.12. This function may be used on blobs with the `PL_BLOB_NOCOPY` flag set and the blob type implements the [release()](foreigninclude.html#release()) callback. It causes the [release()](foreigninclude.html#release()) callback to be called, after which the data and size are set to 0 if the [release()](foreigninclude.html#release()) returns `TRUE`. After this sequence, the [release()](foreigninclude.html#release()) for this blob is never called again. The related `atom_t` handle remains valid until it is no longer referenced and reclaimed by the atom garbage collector. If the blob data is accessed using e.g., [PL_get_blob()](foreigninclude.html#PL_get_blob()) it returns `NULL` for the data and 0 for the size.^(228This means that any predicates or callbacks that use the blob must check the result of [PL_blob_data()](foreigninclude.html#PL_blob_data()).) If the [release()](foreigninclude.html#release()) function is not called, or if it returns `FALSE`, `FALSE` is returned.

[PL_free_blob()](foreigninclude.html#PL_free_blob()) may be called multiple times on the same `atom_t`, provided the handle is still valid. Subsequent calls after a successful call have no effect and return `FALSE`.

#### 12.4.10.3 Considerations for non-C code

The blob API assumes that Prolog will take care of memory management, using the [release(c)](foreigninclude.html#release())allback to handle any cleanup.

Other programming languages have their own memory management, which might not fit nicely with the Prolog memory management. For more details on blobs written with C++, see [C++ interface to SWI-Prolog (Version 2)](https://www.swi-prolog.org/pldoc/man?section=cpp2).

### 12.4.11 Exchanging GMP numbers

If SWI-Prolog is linked with the GNU Multiple Precision Arithmetic Library (GMP, used by default), the foreign interface provides functions for exchanging numeric values to GMP types. To access these functions the header `<gmp.h>` must be included *before* `<SWI-Prolog.h>`. Foreign code using GMP linked to SWI-Prolog asks for some considerations.

- SWI-Prolog normally rebinds the GMP allocation functions using **mp_set_memory_functions()**. This means SWI-Prolog must be initialised before the foreign code touches any GMP function. You can call `PL_action(PL_GMP_SET_ALLOC_FUNCTIONS, TRUE)` to force Prolog's GMP initialization without doing the rest of the Prolog initialization. If you do not want Prolog rebinding the GMP allocation, call `PL_action(PL_GMP_SET_ALLOC_FUNCTIONS, FALSE)` *before* initializing Prolog.
- On Windows, each DLL has its own memory pool. To make exchange of GMP numbers between Prolog and foreign code possible you must either let Prolog rebind the allocation functions (default) or you must recompile SWI-Prolog to link to a DLL version of the GMP library.

Here is an example exploiting the function **mpz_nextprime()**:

``` code
#include <gmp.h>
#include <SWI-Prolog.h>

static foreign_t
next_prime(term_t n, term_t prime)
{ mpz_t mpz;
  int rc;

  mpz_init(mpz);
  if ( PL_get_mpz(n, mpz) )
  { mpz_nextprime(mpz, mpz);

    rc = PL_unify_mpz(prime, mpz);
  } else
    rc = FALSE;

  mpz_clear(mpz);
  return rc;
}

install_t
install()
{ PL_register_foreign("next_prime", 2, next_prime, 0);
}
```

`bool` **PL_get_mpz**(`term_t t, mpz_t mpz`)  
If `t` represents an integer, `mpz` is filled with the value and the function returns `TRUE`. Otherwise `mpz` is untouched and the function returns `FALSE`. Note that `mpz` must have been initialised before calling this function and must be cleared using **mpz_clear()** to reclaim any storage associated with it.

`bool` **PL_get_mpq**(`term_t t, mpq_t mpq`)  
If `t` is an integer or rational number (term `rdiv/2`), `mpq` is filled with the *normalised* rational number and the function returns `TRUE`. Otherwise `mpq` is untouched and the function returns `FALSE`. Note that `mpq` must have been initialised before calling this function and must be cleared using **mpq_clear()** to reclaim any storage associated with it.

`bool` **PL_unify_mpz**(`term_t t, mpz_t mpz`)  
Unify `t` with the integer value represented by `mpz` and return `TRUE` on success. The `mpz` argument is not changed.

`bool` **PL_unify_mpq**(`term_t t, mpq_t mpq`)  
Unify `t` with a rational number represented by `mpq` and return `TRUE` on success. Note that `t` is unified with an integer if the denominator is 1. The `mpq` argument is not changed.

### 12.4.12 Calling Prolog from C

The Prolog engine can be called from C. There are two interfaces for this. For the first, a term is created that could be used as an argument to [call/1](metacall.html#call/1), and then [PL_call()](foreigninclude.html#PL_call()) is used to call Prolog. This system is simple, but does not allow to inspect the different answers to a non-deterministic goal and is relatively slow as the runtime system needs to find the predicate. The other interface is based on [PL_open_query()](foreigninclude.html#PL_open_query()), [PL_next_solution()](foreigninclude.html#PL_next_solution()), and [PL_cut_query()](foreigninclude.html#PL_cut_query()) or [PL_close_query()](foreigninclude.html#PL_close_query()). This mechanism is more powerful, but also more complicated to use.

#### 12.4.12.1 Predicate references

This section discusses the functions used to communicate about predicates. Though a Prolog predicate may be defined or not, redefined, etc., a Prolog predicate has a handle that is neither destroyed nor moved. This handle is known by the type `predicate_t`.

`predicate_t` **PL_pred**(`functor_t f, module_t m`)  
Return a handle to a predicate for the specified name/arity in the given module. If the module argument `m` is `NULL`, the current context module is used. If the target predicate does not exist a handle to a new *undefined* predicate is returned. The predicate may fail, returning `(predicate_t)0` after setting a resource exception, if the target module has a limit on the `program_space`, see [set_module/1](manipmodule.html#set_module/1). Currently aborts the process with a *fatal error* when out of memory. Future versions may raise a resource exception and return `(predicate_t)0`.

`predicate_t` **PL_predicate**(`const char *name, int arity, const char* module`)  
Same as [PL_pred()](foreigninclude.html#PL_pred()), but provides a more convenient interface to the C programmer. If the module argument `module` is `NULL`, the current context module is used. The `predicate_t` handle may be stored as global data and reused for future queries^(229[PL_predicate()](foreigninclude.html#PL_predicate()) involves 5 hash lookups (two to get the atoms, one to get the module, one to get the functor and the final one to get the predicate associated with the functor in the module)) as illustrated below.

``` code
static predicate_t p = 0;

  ...
  if ( !p )
    p = PL_predicate("is_a", 2, "database");
```

Note that [PL_cleanup()](foreigninclude.html#PL_cleanup()) invalidates the predicate handle. Foreign libraries that use the above mechanism must implement the module **uninstall()** function to clear the `predicate_t` global variable.

`bool` **PL_predicate_info**(`predicate_t p, atom_t *n, size_t *a, module_t *m`)  
Return information on the predicate `p`. The name is stored over `n`, the arity over `a`, while `m` receives the definition module. Note that the latter need not be the same as specified with [PL_predicate()](foreigninclude.html#PL_predicate()). If the predicate is imported into the module given to [PL_predicate()](foreigninclude.html#PL_predicate()), this function will return the module where the predicate is defined. Any of the arguments `n`, `a` and `m` can be `NULL`. Currently always returns `TRUE`.

#### 12.4.12.2 Initiating a query from C

This section discusses the functions for creating and manipulating queries from C. Note that a foreign context can have at most one active query. This implies that it is allowed to make strictly nested calls between C and Prolog (Prolog calls C, calls Prolog, calls C, etc.), but it is **not** allowed to open multiple queries and start generating solutions for each of them by calling [PL_next_solution()](foreigninclude.html#PL_next_solution()). Be sure to call [PL_cut_query()](foreigninclude.html#PL_cut_query()) or [PL_close_query()](foreigninclude.html#PL_close_query()) on any query you opened before opening the next or returning control back to Prolog. Failure to do so results in "undefined behavior" (typically, a crash).

`qid_t` **PL_open_query**(`module_t ctx, int flags, predicate_t p, term_t +t0`)  
Opens a query and returns an identifier for it. `ctx` is the *context module* of the goal. When `NULL`, the context module of the calling context will be used, or `user` if there is no calling context (as may happen in embedded systems). Note that the context module only matters for *meta-predicates*. See [meta_predicate/1](metapred.html#meta_predicate/1), [context_module/1](ctxmodule.html#context_module/1) and [module_transparent/1](ctxmodule.html#module_transparent/1). The term reference `t0` is the first of a vector of term references as returned by [PL_new_term_refs(n)](foreigntypes.html#PL_new_term_refs()). Raise a resource exception and returns `(qid_t)0` on failure.

Every use of [PL_open_query()](foreigninclude.html#PL_open_query()) must have a corresponding call to [PL_cut_query()](foreigninclude.html#PL_cut_query()) or [PL_close_query()](foreigninclude.html#PL_close_query()) before the foreign predicate returns either `TRUE` or `FALSE`.

The `flags` arguments provides some additional options concerning debugging and exception handling. It is a bitwise *or* of the following values below. Note that exception propagation is defined by the flags `PL_Q_NORMAL`, `PL_Q_CATCH_EXCEPTION` and `PL_Q_PASS_EXCEPTION`. Exactly one of these flags must be specified (if none of them is specified, the behavior is as if `PL_Q_NODEBUG` is specified)..

**`PL_Q_NORMAL`**  
Normal operation. It is named "normal" because it makes a call to Prolog behave as it did before exceptions were implemented, i.e., an error (now uncaught exception) triggers the debugger. See also the Prolog flag [debug_on_error](flags.html#flag:debug_on_error). This mode is still useful when calling Prolog from C if the C code is not willing to handle exceptions.

**`PL_Q_NODEBUG`**  
Switch off the debugger while executing the goal. This option is used by many calls to hook-predicates to avoid tracing the hooks. An example is [print/1](termrw.html#print/1) calling [portray/1](termrw.html#portray/1) from foreign code. This is the default if flags is `0`.

**`PL_Q_CATCH_EXCEPTION`**  
If an exception is raised while executing the goal, make it available by calling [`PL_exception(qid)`](foreigninclude.html#PL_exception()), where `qid` is the `qid_t` returned by [PL_open_query()](foreigninclude.html#PL_open_query()). The exception is implicitly cleared from the environment when the query is closed and the exception term returned from [`PL_exception(qid)`](foreigninclude.html#PL_exception()) becomes invalid. Use `PL_Q_PASS_EXCEPTION` if you wish to propagate the exception.

**`PL_Q_PASS_EXCEPTION`**  
As `PL_Q_CATCH_EXCEPTION`, making the exception on the inner environment available using [`PL_exception(0)`](foreigninclude.html#PL_exception()) in the parent environment. If [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns `FALSE`, you must call [PL_cut_query()](foreigninclude.html#PL_cut_query()) or [PL_close_query()](foreigninclude.html#PL_close_query()). After that you may verify whether failure was due to logical failure of the called predicate or an exception by calling [`PL_exception(0)`](foreigninclude.html#PL_exception()). If the predicate failed due to an exception you should return with `FALSE` from the foreign predicate or call [PL_clear_exception()](foreigninclude.html#PL_clear_exception()) to clear it. If you wish to process the exception in C, it is advised to use `PL_Q_CATCH_EXCEPTION` instead, but only if you have no need to raise an exception or re-raise the caught exception.

Note that `PL_Q_PASS_EXCEPTION` is used by the debugger to decide whether the exception is *caught*. If there is no matching [catch/3](exception.html#catch/3) call in the current query and the query was started using `PL_Q_PASS_EXCEPTION` the debugger searches the parent queries until it either finds a matching [catch/3](exception.html#catch/3), a query with `PL_Q_CATCH_EXCEPTION` (in which case it considers the exception handled by C) or the top of the query stack (in which case it considers the exception *uncaught*). Uncaught exceptions use the `library(library(prolog_stack))` to add a backtrace to the exception and start the debugger as soon as possible if the Prolog flag [debug_on_error](flags.html#flag:debug_on_error) is `true`.

**`PL_Q_ALLOW_YIELD`**  
Support the `I_YIELD` instruction for engine-based coroutining. See \$engine_yield/2 in `boot/init.pl` for details.

**`PL_Q_TRACE_WITH_YIELD`**  
Allows for implementing a *yield* based debugger. See [section 12.4.1.3](foreigninclude.html#sec:12.4.1.3)

**`PL_Q_EXT_STATUS`**  
Make [PL_next_solution()](foreigninclude.html#PL_next_solution()) return extended status. Instead of only `TRUE` or `FALSE` extended status as illustrated in the following table:

|  |  |  |
|----|----|----|
| **Extended** | **Normal** |  |
| PL_S_NOT_INNER | FALSE | [PL_next_solution()](foreigninclude.html#PL_next_solution()) may only be called on the innermost query |
| PL_S_EXCEPTION | FALSE | Exception available through [PL_exception()](foreigninclude.html#PL_exception()) |
| PL_S_FALSE | FALSE | Query failed |
| PL_S_TRUE | TRUE | Query succeeded with choicepoint |
| PL_S_LAST | TRUE | Query succeeded without choicepoint |
| PL_S_YIELD | n/a | Query was yielded. See [section 12.4.1.2](foreigninclude.html#sec:12.4.1.2) |
| PL_S_YIELD_DEBUG | n/a | Yielded on behalf of the debugger. See [section 12.4.1.3](foreigninclude.html#sec:12.4.1.3) |

[PL_open_query()](foreigninclude.html#PL_open_query()) can return the query identifier `0` if there is not enough space on the environment stack (and makes the exception available through [`PL_exception(0)`](foreigninclude.html#PL_exception())). This function succeeds, even if the referenced predicate is not defined. In this case, running the query using [PL_next_solution()](foreigninclude.html#PL_next_solution()) may return an existence_error. See [PL_exception()](foreigninclude.html#PL_exception()).

The example below opens a query to the predicate is_a/2 to find the ancestor of‘me’. The reference to the predicate is valid for the duration of the process or until [PL_cleanup()](foreigninclude.html#PL_cleanup()) is called (see [PL_predicate()](foreigninclude.html#PL_predicate()) for details) and may be cached by the client.

``` code
char *
ancestor(const char *me)
{ term_t a0 = PL_new_term_refs(2);
  static predicate_t p;

  if ( !p )
    p = PL_predicate("is_a", 2, "database");

  PL_put_atom_chars(a0, me);
  PL_open_query(NULL, PL_Q_PASS_EXCEPTION, p, a0);
  ...
}
```

`int` **PL_next_solution**(`qid_t qid`)  
Generate the first (next) solution for the given query. The return value is `TRUE` if a solution was found, or `FALSE` to indicate the query could not be proven. This function may be called repeatedly until it fails to generate all solutions to the query. The return value `PL_S_NOT_INNER` is returned if `qid` is not the innermost query.

If the [PL_open_query()](foreigninclude.html#PL_open_query()) had the flag `PL_Q_EXT_STATUS`, there are additional return values (see [section 12.4.1.2](foreigninclude.html#sec:12.4.1.2)).

`int` **PL_cut_query**(`qid_t qid`)  
Discards the query, but does not delete any of the data created by the query. It just invalidates `qid`, allowing for a new call to [PL_open_query()](foreigninclude.html#PL_open_query()) in this context. [PL_cut_query()](foreigninclude.html#PL_cut_query()) may invoke cleanup handlers (see [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3)) and therefore may experience exceptions. If an exception occurs the return value is `FALSE` and the exception is accessible through [`PL_exception(0)`](foreigninclude.html#PL_exception()).

An example of a handler that can trigger an exception in [PL_cut_query()](foreigninclude.html#PL_cut_query()) is:

``` code
test_setup_call_cleanup(X) :-
    setup_call_cleanup(
        true,
        between(1, 5, X),
        throw(error)).
```

where [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns `TRUE` on the first result and the `throw(error)` will only run when [PL_cut_query()](foreigninclude.html#PL_cut_query()) or [PL_close_query()](foreigninclude.html#PL_close_query()) is run. On the other hand, if the goal in [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) has completed (failure, exception, deterministic success), the cleanup handler will have done its work before control gets back to Prolog and therefore [PL_next_solution()](foreigninclude.html#PL_next_solution()) will have generated the exception. The return value `PL_S_NOT_INNER` is returned if `qid` is not the innermost query.

`int` **PL_close_query**(`qid_t qid`)  
As [PL_cut_query()](foreigninclude.html#PL_cut_query()), but all data and bindings created by the query are destroyed as if the query is called as `\+ \+ Goal`. This reduces the need for garbage collection, but also rewinds side effects such as setting global variables using [b_setval/2](gvar.html#b_setval/2). The return value `PL_S_NOT_INNER` is returned if `qid` is not the innermost query.

`qid_t` **PL_current_query**(`void`)  
Returns the query id of the current query or `0` if the current thread is not executing any queries.

`PL_engine_t` **PL_query_engine**(`qid_t qid`)  
Return the engine to which `qid` belongs. Note that interacting with a query or the Prolog terms associated with a query requires the engine to be *current*. See [PL_set_engine()](foreignthread.html#PL_set_engine()).

`term_t` **PL_query_arguments**(`qid_t qid`)  
Return a `term_t` handle to the first argument of the main goal of the query `qid`. This allows for enumerating a query and acting on the binding of one of the arguments without additional context. Note that the returned `term_t` is *not* the same handle that was used in [PL_open_query()](foreigninclude.html#PL_open_query()) to pass the arguments. The content of the returned vector, however, is the same.

`void*` **PL_set_query_data**(`qid_t qid, unsigned int offset, void* data`)  
`void*` **PL_query_data**(`qid_t qid, unsigned int offset`)  
Associate user data with a query. `offset` must be smaller than `PL_MAX_QUERY_DATA` (currently 2). [PL_set_query_data()](foreigninclude.html#PL_set_query_data()) returns the old value.

`bool` **PL_call_predicate**(`module_t m, int flags, predicate_t pred, term_t +t0`)  
Shorthand for [PL_open_query()](foreigninclude.html#PL_open_query()), [PL_next_solution()](foreigninclude.html#PL_next_solution()), [PL_cut_query()](foreigninclude.html#PL_cut_query()), generating a single solution. The arguments are the same as for [PL_open_query()](foreigninclude.html#PL_open_query()), the return value is the same as [PL_next_solution()](foreigninclude.html#PL_next_solution()).

`bool` **PL_call**(`term_t t, module_t m`)  
Call term `t` just like the Prolog predicate [once/1](metacall.html#once/1). `t` is called in the module `m`, or in the context module if `m` == NULL. Returns `TRUE` if the call succeeds, `FALSE` otherwise. If the goal raises an exception the return value is `FALSE` and the exception term is available using [PL_exception(0)](foreigninclude.html#PL_exception()).^(230Up to version 9.1.11 the debugger was started and the exception was not propagated.) [Figure 7](foreigninclude.html#fig:calling) shows an example to obtain the number of defined atoms.

### 12.4.13 Discarding Data

The Prolog data created and term references needed to set up the call and/or analyse the result can in most cases be discarded right after the call. [PL_close_query()](foreigninclude.html#PL_close_query()) allows for destroying the data, while leaving the term references. The calls below may be used to destroy term references and data. See [figure 7](foreigninclude.html#fig:calling) for an example.

`fid_t` **PL_open_foreign_frame**()  
Create a foreign frame, holding a mark that allows the system to undo bindings and destroy data created after it, as well as providing the environment for creating term references. Each call to a foreign predicate is wrapped in a [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()) and [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()) pair. This ensures that `term_t` handles created during the execution of a foreign predicate are scoped to this execution. Note that if the foreign predicate is *non-deterministic*, `term_t` handles are scoped to *each* activation of the foreign function.

The user may create explicit foreign frames to undo (backtrack) changes to Prolog terms. See [PL_unify()](foreigninclude.html#PL_unify()) for an example. An explicit foreign frame must also be used for creating a callback from C to Prolog (see [PL_open_query()](foreigninclude.html#PL_open_query())) to ensure the existence of such a frame and to scope the `term_t` handles needed to setup the call to Prolog.

On success, the stack has room for at least 10 `term_t` handles. This implies that foreign predicates as well as code inside an explicitly created foreign frame may use [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()) to create up to 10 `term_t` handles without checking the return status.

Returns `(fid_t)0` on failure. Failure is either lack of space on the stacks, in which case a resource exception is scheduled or atom-gc being in progress in the current thread, in which case no exception is scheduled. The latter is an exceptional case that prevents doing a callback on Prolog from *blob release handlers*.^(231Such a callback would *deadlock* if the callback creates new atoms or requires stack shifts or garbage collection.)

`void` **PL_close_foreign_frame**(`fid_t id`)  
Discard all term references created after the frame was opened. Prolog data referenced by the discarded term references is not affected.

`void` **PL_discard_foreign_frame**(`fid_t id`)  
Same as [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()), but also undo all bindings made since the open and destroy all Prolog data.

`void` **PL_rewind_foreign_frame**(`fid_t id`)  
Undo all bindings and discard all term references created since the frame was created, but do not pop the frame. That is, the same frame can be rewound multiple times, and must eventually be closed or discarded.

It is obligatory to call either of the two closing functions to discard a foreign frame. Foreign frames may be nested.

``` code
int
count_atoms()
{ fid_t fid = PL_open_foreign_frame();
  term_t goal  = PL_new_term_ref();
  term_t a1    = PL_new_term_ref();
  term_t a2    = PL_new_term_ref();
  functor_t s2 = PL_new_functor(PL_new_atom("statistics"), 2);
  int atoms;

  PL_put_atom_chars(a1, "atoms");
  PL_cons_functor(goal, s2, a1, a2);
  PL_call(goal, NULL);         /* call it in current module */

  PL_get_integer(a2, &atoms);
  PL_discard_foreign_frame(fid);

  return atoms;
}
```

**Figure 7 :** Calling Prolog

### 12.4.14 String buffering

Many of the functions of the foreign language interface involve strings. Some of these strings point into static memory like those associated with atoms. These strings are valid as long as the atom is protected against atom garbage collection, which generally implies the atom must be locked using [PL_register_atom()](foreigninclude.html#PL_register_atom()) or be part of an accessible term. Other strings are more volatile. Several functions provide a BUF\_\* flag that can be set to either `BUF_STACK` (default) or `BUF_MALLOC`.^(232`BUF_DISCARDABLE` is also defined - it can often avoid allocating a new string but in some situations it behaves the same as `BUF_STACK`.) Strings returned by a function accepting `BUF_MALLOC` **must** be freed using [PL_free()](foreignnotes.html#PL_free()). Strings returned using `BUF_STACK` are pushed on a stack that is cleared when a foreign predicate returns control back to Prolog More fine grained control may be needed if functions that return strings are called outside the context of a foreign predicate or a foreign predicate creates many strings during its execution. Temporary strings are scoped using these macros:

`void` **PL_STRINGS_MARK**()  
`void` **PL_STRINGS_RELEASE**()  
These macros must be paired and create a C *block* ({...}). Any string created using `BUF_STACK` after [PL_STRINGS_MARK()](foreigninclude.html#PL_STRINGS_MARK()) is released by the corresponding [PL_STRINGS_RELEASE()](foreigninclude.html#PL_STRINGS_RELEASE()). These macros should be used like below. Note that strings returned by any of the Prolog functions between this pair may be invalidated.

``` code
  ...
  PL_STRINGS_MARK();
    <operations involving strings>
  PL_STRINGS_RELEASE();
  ...
```

The Prolog flag [string_stack_tripwire](flags.html#flag:string_stack_tripwire) may be used to set a *tripwire* to help finding places where scoping strings may help reducing resources.

### 12.4.15 Foreign Code and Modules

Modules are identified via a unique handle. The following functions are available to query and manipulate modules.

`module_t` **PL_context**()  
Return the module identifier of the context module of the currently active foreign predicate. If there is no currently active predicate it returns a handle to the `user` module.

`bool` **PL_strip_module**(`term_t +raw, module_t *m, term_t -plain`)  
Utility function. If `raw` is a term, possibly holding the module construct \<`module`\>`:`\<`rest`\>, this function will make `plain` a reference to \<`rest`\> and fill `module *` with \<`module`\>. For further nested module constructs the innermost module is returned via `module *`. If `raw` is not a module construct, `raw` will simply be put in `plain`. The value pointed to by `m` must be initialized before calling [PL_strip_module()](foreigninclude.html#PL_strip_module()), either to the default module or to `NULL`. A `NULL` value is replaced by the current context module if `raw` carries no module. The following example shows how to obtain the plain term and module if the default module is the user module:

``` code
{ module m = PL_new_module(PL_new_atom("user"));
  term_t plain = PL_new_term_ref();

  PL_strip_module(term, &m, plain);
  ...
}
```

Returns `TRUE` on success and `FALSE` on error, leaving an exception. Currently the only exception condition is `raw` to be a cyclic term.

`atom_t` **PL_module_name**(`module_t module`)  
Return the name of `module` as an atom.

`module_t` **PL_new_module**(`atom_t name`)  
Find an existing module or create a new module with the name `name`. Currently aborts the process with a *fatal error* on failure. Future versions may raise a resource exception and return `(module_t)0`.

### 12.4.16 Prolog exceptions in foreign code

This section discusses [PL_exception()](foreigninclude.html#PL_exception()) and [PL_raise_exception()](foreigninclude.html#PL_raise_exception()), the interface functions to detect and generate Prolog exceptions from C code. [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) from the C interface registers the exception term and returns `FALSE`. If a foreign predicate returns `FALSE`, while an exception term is registered, a Prolog exception will be raised by the virtual machine. This implies for a foreign function that implements a predicate and wishes to raise an exception, the function must call [PL_raise_exception()](foreigninclude.html#PL_raise_exception()), perform any necessary cleanup, and return the return code of [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) or explicitly `FALSE`. Calling [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) outside the context of a function implementing a foreign predicate results in undefined behaviour.

Note that many of the C API functions may call [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) and return `FALSE`. The user must test for this, cleanup, and make the foreign function return `FALSE`.

[PL_exception()](foreigninclude.html#PL_exception()) may be used to inspect the currently registered exception. It is normally called after a call to [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns `FALSE`, and returns a term reference to an exception term if an exception is pending, and `(term_t)0` otherwise. It may also be called after, e.g., [PL_unify()](foreigninclude.html#PL_unify()) to distinguish a normal failing unification from a unification that raised an resource error exception. [PL_exception()](foreigninclude.html#PL_exception()) must only be called after a function such as [PL_next_solution()](foreigninclude.html#PL_next_solution()) or [PL_unify()](foreigninclude.html#PL_unify()) returns failure; if called elsewhere, the return value is undefined.

If a C function implementing a predicate that calls Prolog should use [PL_open_query()](foreigninclude.html#PL_open_query()) with the flag `PL_Q_PASS_EXCEPTION` and make the function return FALSE if [PL_next_solution()](foreigninclude.html#PL_next_solution()) returns `FALSE` and [PL_exception()](foreigninclude.html#PL_exception()) indicates an exception is pending.

Both for C functions implementing a predicate and when Prolog is called while the main control of the process is in C, user code should always check for exceptions. As explained above, C functions implementing a predicate should normally cleanup and return with `FALSE`. If the C function wishes to continue it may call [PL_clear_exception()](foreigninclude.html#PL_clear_exception()). Note that this may cause any exception to be ignored, including *time outs* and *abort*. Typically the user should check the exception details before ignoring an exception (using [`PL_exception(0)`](foreigninclude.html#PL_exception()) or [`PL_exception(qid)`](foreigninclude.html#PL_exception()) as appropriate). If the C code does not implement a predicate it normally prints the exception and calls [PL_clear_exception()](foreigninclude.html#PL_clear_exception()) to discard it. Exceptions may be printed by calling [print_message/2](printmsg.html#print_message/2) through the C interface.

`bool` **PL_raise_exception**(`term_t exception`)  
Generate an exception (as [throw/1](exception.html#throw/1)) and return `FALSE`. If there is already a pending exception, the most urgent exception is kept; and if both are of the same urgency, the new exception is kept. Urgency of exceptions is described in secrefurgentexceptions.

This function is rarely used directly. Instead, errors are typically raised using the functions in [section 12.4.7](foreigninclude.html#sec:12.4.7) or the C API functions that end in `_ex` such as [PL_get_atom_ex()](foreigninclude.html#PL_get_atom_ex()). Below we give an example returning an exception from a foreign predicate the verbose way. Note that the exception is raised in a sequence of actions connected using `&&`. This ensures that a proper exception is raised should any of the calls used to build or raise the exception themselves raise an exception. In this simple case [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()) is guaranteed to succeed because the system guarantees at least 10 available term references before entering the foreign predicate. [PL_unify_term()](foreigninclude.html#PL_unify_term()) however may raise a resource exception for the global stack.

``` code
foreign_t
pl_hello(term_t to)
{ char *s;

  if ( PL_get_atom_chars(to, &s) )
  { return Sfprintf(Scurrent_output, "Hello \"%s\"\n", s);
  } else
  { term_t except;

    return  ( (except=PL_new_term_ref()) &&
              PL_unify_term(except,
                            PL_FUNCTOR_CHARS, "type_error", 2,
                              PL_CHARS, "atom",
                              PL_TERM, to) &&
              PL_raise_exception(except) );
  }
}
```

For reference, the preferred implementation of the above is below. The `CVT_EXCEPTION` tells the system to generate an exception if the conversion fails. The other `CVT_` flags define the admissible types and `REP_MB` requests the string to be provided in the current *locale* representation. This implies that Unicode text is printed correctly if the current environment can represent it. If not, a `representation_error` is raised.

``` code
foreign_t
pl_hello(term_t to)
{ char *s;

  if ( PL_get_chars(to, &s, CVT_ATOM|CVT_STRING|CVT_EXCEPTION|REP_MB) )
  { return Sfprintf(Scurrent_output, "Hello \"%s\"\n", s);
  }

  return FALSE;
}
```

`bool` **PL_throw**(`term_t exception`)  
Similar to [PL_raise_exception()](foreigninclude.html#PL_raise_exception()), but returns using the C **longjmp()** function to the innermost [PL_next_solution()](foreigninclude.html#PL_next_solution()). This function is deprecated as it does not provide the opportunity to cleanup.

`term_t` **PL_exception**(`qid_t qid`)  
Return the *pending exception*. Exceptions may be raised by most of the API calls described in this chapter, a common possibility being `resource_error` exceptions. Some return `type_error` or `domain_error` exceptions. A call to Prolog using [PL_next_solution()](foreigninclude.html#PL_next_solution()) may return any exception, including those thrown by explicit calls to [throw/1](exception.html#throw/1). If no exception is pending this function returns `(term_t)0`.

Normally `qid` should be `0`. An explicit `qid` must be used after a call to [PL_next_solution()](foreigninclude.html#PL_next_solution()) that returns `FALSE` when the query was created using the `PL_Q_PASS_EXCEPTION` flag (see [PL_open_query()](foreigninclude.html#PL_open_query())).

Note that an API may only raise an exception when it fails; if the API call succeeds, the result of [`PL_exception(0)`](foreigninclude.html#PL_exception()) will be 0.^(233Provided no exception was pending before calling the API function. As clients must deal with exceptions immediately after an API call raises one, this can not happen in a well behaved client.) The implementation of a foreign predicate should normally cleanup and return `FALSE` after an exception is raised (and typically also after an API call failed for logical reasons; see [PL_unify()](foreigninclude.html#PL_unify()) for an elaboration on this topic). If the call to Prolog is not the implementation of a foreign predicate, e.g., when the overall process control is in some other language, exceptions may be printed by calling [print_message/2](printmsg.html#print_message/2) and should be discarded by calling [PL_clear_exception()](foreigninclude.html#PL_clear_exception()).

`void` **PL_clear_exception**(`void`)  
Tells Prolog that the encountered exception must be ignored. This function must be called if control remains in C after a previous API call fails with an exception.^(234This feature is non-portable. Other Prolog systems (e.g., YAP) have no facilities to ignore raised exceptions, and the design of YAP's exception handling does not support such a facility.) If there is no pending exception, [PL_clear_exception()](foreigninclude.html#PL_clear_exception()) does nothing.

### 12.4.17 Catching Signals (Software Interrupts)

SWI-Prolog offers both a C and Prolog interface to deal with software interrupts (signals). The Prolog mapping is defined in [section 4.12](signal.html#sec:4.12). This subsection deals with handling signals from C.

If a signal is not used by Prolog and the handler does not call Prolog in any way, the native signal interface routines may be used.

Any handler that wishes to call one of the Prolog interface functions should call [PL_sigaction()](foreigninclude.html#PL_sigaction()) to install the handler. [PL_signal()](foreigninclude.html#PL_signal()) provides a deprecated interface that is notably not capable of properly restoring the old signal status if the signal was previously handled by Prolog.

`int` **PL_sigaction**(`int sig, pl_sigaction_t *act, pl_sigaction_t *oldact`)  
Install or query the status for signal `sig`. The signal is an integer between 1 and 64, where the where the signals up to 32 are mapped to OS signals and signals above that are handled by Prolog's synchronous signal handling. The `pl_sigaction_t` is a struct with the following definition:

``` code
typedef struct pl_sigaction
{ void        (*sa_cfunction)(int);     /* traditional C function */
  predicate_t sa_predicate;             /* call a predicate */
  int         sa_flags;                 /* additional flags */
} pl_sigaction_t;
```

The `sa_flags` is a bitwise or of `PLSIG_THROW`, `PLSIG_SYNC` and `PLSIG_NOFRAME`. Signal handling is enabled if `PLSIG_THROW` is provided, `sa_cfunction` or `sa_predicate` is provided. `sa_predicate` is a predicate handle for a predicate with arity 1. If no action is provided the signal handling for this signal is restored to the default before [PL_initialise()](foreigninclude.html#PL_initialise()) was called.

Finally, 0 (zero) may be passed for `sig`. In that case the system allocates a free signal in the *Prolog range* (32 ... 64). Such signal handler are activated using **PL_thread_raise()**.

`void (*)()` **PL_signal**(`sig, func`)  
This function is equivalent to the BSD-Unix **signal()** function, regardless of the platform used. The signal handler is blocked while the signal routine is active, and automatically reactivated after the handler returns.

After a signal handler is registered using this function, the native signal interface redirects the signal to a generic signal handler inside SWI-Prolog. This generic handler validates the environment, creates a suitable environment for calling the interface functions described in this chapter and finally calls the registered user-handler.

By default, signals are handled asynchronously (i.e., at the time they arrive). It is inherently dangerous to call extensive code fragments, and especially exception related code from asynchronous handlers. The interface allows for *synchronous* handling of signals. In this case the native OS handler just schedules the signal using [PL_raise()](foreigninclude.html#PL_raise()), which is checked by [PL_handle_signals()](foreigninclude.html#PL_handle_signals()) at the call- and redo-port. This behaviour is realised by *or*-ing `sig` with the constant `PL_SIGSYNC`.^(235A better default would be to use synchronous handling, but this interface preserves backward compatibility.)

Signal handling routines may raise exceptions using [PL_raise_exception()](foreigninclude.html#PL_raise_exception()). The use of [PL_throw()](foreigninclude.html#PL_throw()) is not safe. If a synchronous handler raises an exception, the exception is delayed to the next call to [PL_handle_signals()](foreigninclude.html#PL_handle_signals());

`bool` **PL_raise**(`int sig`)  
Register `sig` for *synchronous* handling by Prolog. Synchronous signals are handled at the call-port or if foreign code calls [PL_handle_signals()](foreigninclude.html#PL_handle_signals()). See also [thread_signal/2](threadcom.html#thread_signal/2).

`int` **PL_handle_signals**(`void`)  
Handle any signals pending from [PL_raise()](foreigninclude.html#PL_raise()). [PL_handle_signals()](foreigninclude.html#PL_handle_signals()) is called at each pass through the call- and redo-port at a safe point. Exceptions raised by the handler using [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) are properly passed to the environment.

The user may call this function inside long-running foreign functions to handle scheduled interrupts. This routine returns the number of signals handled. If a handler raises an exception, the return value is -1 and the calling routine should return with `FALSE` as soon as possible.

`int` **PL_get_signum_ex**(`term_t t, int *sig`)  
Extract a signal specification from a Prolog term and store as an integer signal number in `sig`. The specification is an integer, a lowercase signal name without `SIG` or the full signal name. These refer to the same: `9`, `kill` and `SIGKILL`. Leaves a typed, domain or instantiation error if the conversion fails.

### 12.4.18 Miscellaneous

#### 12.4.18.1 Term Comparison

`int` **PL_compare**(`term_t t1, term_t t2`)  
Compares two terms using the standard order of terms and returns -1, 0 or 1. See also [compare/3](compare.html#compare/3).

`bool` **PL_same_compound**(`term_t t1, term_t t2`)  
Yields `TRUE` if `t1` and `t2` refer to physically the same compound term and `FALSE` otherwise.

#### 12.4.18.2 Recorded database

In some applications it is useful to store and retrieve Prolog terms from C code. For example, the XPCE graphical environment does this for storing arbitrary Prolog data as slot-data of XPCE objects.

Please note that the returned handles have no meaning at the Prolog level and the recorded terms are not visible from Prolog. The functions [PL_recorded()](foreigninclude.html#PL_recorded()) and [PL_erase()](foreigninclude.html#PL_erase()) are the only functions that can operate on the stored term.

Two groups of functions are provided. The first group ([PL_record()](foreigninclude.html#PL_record()) and friends) store Prolog terms on the Prolog heap for retrieval during the same session. These functions are also used by [recorda/3](db.html#recorda/3) and friends. The recorded database may be used to communicate Prolog terms between threads.

`record_t` **PL_record**(`term_t +t`)  
Record the term `t` into the Prolog database as [recorda/3](db.html#recorda/3) and return an opaque handle to the term. The returned handle remains valid until [PL_erase()](foreigninclude.html#PL_erase()) is called on it. [PL_recorded()](foreigninclude.html#PL_recorded()) is used to copy recorded terms back to the Prolog stack. Currently aborts the process with a *fatal error* on failure. Future versions may raise a resource exception and return `(record_t)0`.

`record_t` **PL_duplicate_record**(`record_t record`)  
Return a duplicate of `record`. As records are read-only objects this function merely increments the records reference count. Returns `(record_t)0` if the `record` is an *external record* (see [PL_record_external()](foreigninclude.html#PL_record_external())).

`bool` **PL_recorded**(`record_t record, term_t -t`)  
Copy a recorded term back to the Prolog stack. The same record may be used to copy multiple instances at any time to the Prolog stack. Returns `TRUE` on success, and `FALSE` if there is not enough space on the stack to accommodate the term. See also [PL_record()](foreigninclude.html#PL_record()) and [PL_erase()](foreigninclude.html#PL_erase()).

`void` **PL_erase**(`record_t record`)  
Remove the recorded term from the Prolog database, reclaiming all associated memory resources.

The second group (headed by [PL_record_external()](foreigninclude.html#PL_record_external())) provides the same functionality, but the returned data has properties that enable storing the data on an external device. It has been designed for fast and compact storage of Prolog terms in an external database. Here are the main features:

- *Independent of session*  
  Records can be communicated to another Prolog session and made visible using [PL_recorded_external()](foreigninclude.html#PL_recorded_external()).
- *Binary*  
  The representation is binary for maximum performance. The returned data may contain zero bytes.
- *Byte-order independent*  
  The representation can be transferred between machines with different byte order.
- *No alignment restrictions*  
  There are no memory alignment restrictions and copies of the record can thus be moved freely. For example, it is possible to use this representation to exchange terms using shared memory between different Prolog processes.
- *Compact*  
  It is assumed that a smaller memory footprint will eventually outperform slightly faster representations.
- *Stable*  
  The format is designed for future enhancements without breaking compatibility with older records.

`char *` **PL_record_external**(`term_t +t, size_t *len`)  
Similar to [PL_record()](foreigninclude.html#PL_record()), but the term is serialized such that it can be reloaded in another Prolog session. This implies that atoms and functors are stored by their content rather than their handle. As a result, [PL_record_external()](foreigninclude.html#PL_record_external()) fails (returning `NULL` if the term contains *blobs* that cannot be serialized, such as streams.

These functions are used to implement library `library(fastrw)` as well as for storing Prolog terms in external databases such as BerkeleyDB (library `library(bdb)`) or RocksDB. The representation is optimized for plain atoms and numbers.

Records that are used only in the same Prolog process should use [PL_record()](foreigninclude.html#PL_record()) as this can represent any term, is more compact and faster.

The returned string may be copied. Note that the string may contain null bytes and is not null terminated. The length in bytes is returned in `len`. After copying, the returned string may be discarded using [PL_erase_external()](foreigninclude.html#PL_erase_external()).

[PL_recorded_external()](foreigninclude.html#PL_recorded_external()) is used to copy the term represented in the data back to the Prolog stack. [PL_recorded_external()](foreigninclude.html#PL_recorded_external()) can be used on the returned string as well as on a copy.

`bool` **PL_recorded_external**(`const char *record, term_t -t`)  
Copy a recorded term back to the Prolog stack. The same record may be used to copy multiple instances at any time to the Prolog stack. See also [PL_record_external()](foreigninclude.html#PL_record_external()) and [PL_erase_external()](foreigninclude.html#PL_erase_external()).

`bool` **PL_erase_external**(`char *record`)  
Remove the recorded term from the Prolog database, reclaiming all associated memory resources.

#### 12.4.18.3 Database

`bool` **PL_assert**(`term_t t, module_t m, int flags`)  
Provides direct access to [asserta/1](db.html#asserta/1) and [assertz/1](db.html#assertz/1) by asserting `t` into the database in the module `m`. Defined flags are:

**PL_ASSERTZ**  
Add the new clause as last. Calls [assertz/1](db.html#assertz/1). This macros is defined as 0 and thus the default.

**PL_ASSERTA**  
Add the new clause as first. Calls [asserta/1](db.html#asserta/1).

**PL_CREATE_THREAD_LOCAL**  
If the predicate is not defined, create it as thread-local. See [thread_local/1](threadcom.html#thread_local/1).

**PL_CREATE_INCREMENTAL**  
If the predicate is not defined, create it as *incremental* see [table/1](tabling-preds.html#table/1) and [section 7.7](tabling-incremental.html#sec:7.7).

On success this function returns `TRUE`. On failure `FALSE` is returned and an exception is left in the environment that describes the reason of failure. See [PL_exception()](foreigninclude.html#PL_exception()).

This predicate bypasses creating a Prolog callback environment and is faster than setting up a call to [assertz/1](db.html#assertz/1). It may be used together with [PL_chars_to_term()](foreigninclude.html#PL_chars_to_term()), but the typical use case will create a number of clauses for the same predicate. The fastest way to achieve this is by creating a term that represents the invariable structure of the desired clauses using variables for the variable sub terms. Now we can loop over the data, binding the variables, asserting the term and undoing the bindings. Below is an example loading words from a file that contains a word per line.

``` code
#include <SWI-Prolog.h>
#include <stdio.h>
#include <string.h>

#define MAXWLEN 256

static foreign_t
load_words(term_t name)
{ char *fn;

  if ( PL_get_file_name(name, &fn, PL_FILE_READ) )
  { FILE *fd = fopen(fn, "r");
    char word[MAXWLEN];
    module_t m = PL_new_module(PL_new_atom("words"));
    term_t cl = PL_new_term_ref();
    term_t w  = PL_new_term_ref();
    fid_t fid;

    if ( !PL_unify_term(cl, PL_FUNCTOR_CHARS, "word", 1, PL_TERM, w) )
      return FALSE;

    if ( (fid = PL_open_foreign_frame()) )
    { while(fgets(word, sizeof word, fd))
      { size_t len;

        if ( (len=strlen(word)) )
        { word[len-1] = '\0';
          if ( !PL_unify_chars(w, PL_ATOM|REP_MB, (size_t)-1, word) ||
               !PL_assert(cl, m, 0) )
            return FALSE;
          PL_rewind_foreign_frame(fid);
        }
      }

      PL_close_foreign_frame(fid);
    }

    fclose(fd);
    return TRUE;
  }

  return FALSE;
}

install_t
install(void)
{ PL_register_foreign("load_words", 1, load_words, 0);
}
```

#### 12.4.18.4 Getting file names

The function [PL_get_file_name()](foreigninclude.html#PL_get_file_name()) provides access to Prolog filenames and its file-search mechanism described with [absolute_file_name/3](files.html#absolute_file_name/3). Its existence is motivated to realise a uniform interface to deal with file properties, search, naming conventions, etc., from foreign code.

`bool` **PL_get_file_name**(`term_t spec, char **name, int flags`)  
Translate a Prolog term into a file name. The name is stored in the buffer stack described with the [PL_get_chars()](foreigninclude.html#PL_get_chars()) option `BUF_STACK`, which is popped upon return from the foreign predicate to Prolog. Conversion from the internal UNICODE encoding is done using standard C library functions. `flags` is a bit-mask controlling the conversion process. On failure, `PL_FILE_NOERRORS` controls whether an exception is raised. Options are:

**`PL_FILE_ABSOLUTE`**  
Return an absolute path to the requested file.

**`PL_FILE_OSPATH`**  
Return the name using the hosting OS conventions. On MS-Windows, `\` is used to separate directories rather than the canonical `/`.

**`PL_FILE_SEARCH`**  
Invoke [absolute_file_name/3](files.html#absolute_file_name/3). This implies rules from [file_search_path/2](consulting.html#file_search_path/2) are used.

**`PL_FILE_EXIST`**  
Demand the path to refer to an existing entity.

**`PL_FILE_READ`**  
Demand read-access on the result.

**`PL_FILE_WRITE`**  
Demand write-access on the result.

**`PL_FILE_EXECUTE`**  
Demand execute-access on the result.

**`PL_FILE_NOERRORS`**  
Do not raise any exceptions.

`bool` **PL_get_file_nameW**(`term_t spec, wchar_t **name, int flags`)  
Same as [PL_get_file_name()](foreigninclude.html#PL_get_file_name()), but returns the filename as a wide-character string. This is intended for Windows to access the Unicode version of the Win32 API. Note that the flag `PL_FILE_OSPATH` must be provided to fetch a filename in OS native (e.g., `C:\x\y`) notation.

#### 12.4.18.5 Dealing with Prolog flags from C

Foreign code can set or create Prolog flags using [PL_set_prolog_flag()](foreigninclude.html#PL_set_prolog_flag()). See [set_prolog_flag/2](flags.html#set_prolog_flag/2) and [create_prolog_flag/3](flags.html#create_prolog_flag/3). To retrieve the value of a flag you can use [PL_current_prolog_flag()](foreigninclude.html#PL_current_prolog_flag()).

`bool` **PL_set_prolog_flag**(`const char *name, int type, ...`)  
Set/create a Prolog flag from C. `name` is the name of the affected flag. `type` is one of the values below, which also dictates the type of the final argument. The function returns `TRUE` on success and `FALSE` on failure. This function can be called *before* [PL_initialise()](foreigninclude.html#PL_initialise()), making the flag available to the Prolog startup code.

**`PL_BOOL`**  
Create a boolean (`true` or `false`) flag. The argument must be an `int`.

**`PL_ATOM`**  
Create a flag with an atom as value. The argument must be of type `const char *`.

**`PL_INTEGER`**  
Create a flag with an integer as value. The argument must be of type `intptr_t *`.

&nbsp;

`bool` **PL_current_prolog_flag**(`atom_t name, int type, void *value`)  
Retrieve the value of a Prolog flag from C. `name` is the name of the flag as an `atom_t` (see [current_prolog_flag/2](flags.html#current_prolog_flag/2)). `type` specifies the kind of value to be retrieved, it is one of the values below. `value` is a pointer to a location where to store the value. The user is responsible for making sure this memory location is of the appropriate size/type (see the returned types below to determine the size/type). The function returns `TRUE` on success and `FALSE` on failure.

**`PL_ATOM`**  
Retrieve a flag whose value is an `atom`. The returned value is an atom handle of type `atom_t`.

**`PL_BOOL`**  
Retrieve a flag whose value is a `bool`. The returned value is a C boolean of type `bool`. Note that boolean flags are internally represented as one of the atoms `true` or `false`.

**`PL_INTEGER`**  
Retrieve a flag whose value is an `integer`. The returned value is an integer of type `int64_t`.

**`PL_FLOAT`**  
Retrieve a flag whose value is a `float`. The returned value is a floating point number of type `double`.

**`PL_TERM`**  
Retrieve a flag whose value is a `term`. The returned value is a term handle of type `term_t`.

#### 12.4.18.6 Foreign code and Well Founded Semantics

`bool` **PL_get_delay_list**(`term_t -dl`)  
Fetch the current *delay list*. If this list is not empty, the current answer is *undefined*. In the logical sense, this function always succeeds and sets `dl` to the delay list. It returns `FALSE` if the delay list is empty (and answer is well defined) and `TRUE` if the delay list is not empty. If `dl` is 0 no list is instantiated, while the return value is the same. This allows for testing that an answer is undefined as below.

``` code
  if ( PL_get_delay_list(0) )
    <undefined>
  else
    <normal answer>
```

For now, we consider the content of the list elements opaque. See `boot/tabling.pl` for examples.

### 12.4.19 Errors and warnings

[PL_warning()](foreigninclude.html#PL_warning()) prints a standard Prolog warning message to the standard error (`user_error`) stream. Please note that new code should consider using [PL_raise_exception()](foreigninclude.html#PL_raise_exception()) to raise a Prolog exception. See also [section 4.10](exception.html#sec:4.10).

`bool` **PL_warning**(`format, a1, ...`)  
Print an error message starting with‘`[WARNING: `’, followed by the output from `format`, followed by a‘`]`’and a newline. Then start the tracer. `format` and the arguments are the same as for **printf**(2). Always returns `FALSE`.

`int` **PL_fatal_error**(`format, a1, ...`)  
As [PL_warning()](foreigninclude.html#PL_warning()), but using `[FATAL ERROR: at <``time``> ...]` and terminates the process after cleanup using **abort()**. If the process is a Windows GUI application it uses a message box. This function should be used if an unrepairable error is detected. For example, Prolog uses it to signal it cannot find the compiled Prolog startup or memory allocation fails in a place from where we cannot gracefully generate an exception.^(236Currently most memory allocation except for most of the big allocations such as for the Prolog stacks.)

`int` **PL_system_error**(`format, a1, ...`)  
As [PL_fatal_error()](foreigninclude.html#PL_fatal_error()), but using `[ERROR: system error:]` and provides additional technical details such as the thread that trapped the error and backtrace of the C and Prolog stacks. This function should be used to when an unexpected and unrepairable error is detected. For example, Prolog uses this after it finds an inconsistency in the data during garbage collection.

`int` **PL_api_error**(`format, a1, ...`)  
As [PL_system_error()](foreigninclude.html#PL_system_error()), but using `[ERROR: API error:]` and provides additional technical details such as the thread that trapped the error and backtrace of the C and Prolog stacks. This function is used by the C API and may be used by other language bindings to report invalid use of the API. This function causes the process to be terminated.

`bool` **PL_print_message**(`atom_t severity, ...`)  
Calls [print_message/2](printmsg.html#print_message/2) using a term constructed from the remaining arguments that are passed to [PL_unify_term()](foreigninclude.html#PL_unify_term()). This is similar to setting up a call to [print_message/2](printmsg.html#print_message/2) using [PL_call_predicate()](foreigninclude.html#PL_call_predicate()), except that it saves and restores possibly pending exceptions and delayed goals. The `severity` argument must be valid for [print_message/2](printmsg.html#print_message/2).

### 12.4.20 Environment Control from Foreign Code

`bool` **PL_action**(`int, ...`)  
Perform some action on the Prolog system. `int` describes the action. Remaining arguments depend on the requested action. The actions are listed below:

**PL_ACTION_TRACE**  
Start Prolog tracer ([trace/0](debugger.html#trace/0)). Requires no arguments.

**PL_ACTION_DEBUG**  
Switch on Prolog debug mode ([debug/0](debugger.html#debug/0)). Requires no arguments.

**PL_ACTION_BACKTRACE**  
Print backtrace on current output stream. The argument (an `int`) is the number of frames printed.

**PL_ACTION_HALT**  
Halt Prolog execution. This action should be called rather than Unix **exit()** to give Prolog the opportunity to clean up. This call does not return. The argument (an `int`) is the exit code. See [halt/1](toplevel.html#halt/1).

**PL_ACTION_ABORT**  
Generate a Prolog abort ([abort/0](toplevel.html#abort/0)). This call does not return. Requires no arguments.

**PL_ACTION_BREAK**  
Create a standard Prolog break environment ([break/0](toplevel.html#break/0)). Returns after the user types the end-of-file character. Requires no arguments.

**PL_ACTION_GUIAPP**  
Windows: Used to indicate to the kernel that the application is a GUI application if the argument is not 0, and a console application if the argument is 0. If a fatal error occurs, the system uses a windows messagebox to report this on a GUI application, and otherwise simply prints the error and exits.

**PL_ACTION_TRADITIONAL**  
Same effect as using **--traditional**. Must be called *before* [PL_initialise()](foreigninclude.html#PL_initialise()).

**PL_ACTION_WRITE**  
Write the argument, a `char *` to the current output stream.

**PL_ACTION_FLUSH**  
Flush the current output stream. Requires no arguments.

**PL_ACTION_ATTACH_CONSOLE**  
Attach a console to a thread if it does not have one. See [attach_console/0](threadutil.html#attach_console/0).

**PL_GMP_SET_ALLOC_FUNCTIONS**  
Takes an integer argument. If `TRUE`, the GMP allocations are immediately bound to the Prolog functions. If `FALSE`, SWI-Prolog will never rebind the GMP allocation functions. See **mp_set_memory_functions()** in the GMP documentation. The action returns `FALSE` if there is no GMP support or GMP is already initialised.

`unsigned int` **PL_version_info**(`int key`)  
Query version information. This function may be called before [PL_initialise()](foreigninclude.html#PL_initialise()). If the key is unknown the function returns 0. See [section 2.21](abi-versions.html#sec:2.21) for a more in-depth discussion on binary compatibility. Versions up to SWI-Prolog 8.5.2 defined this function as **PL_version()**. It was renamed to avoid a conflict with Perl affecting [Yaswi](https://github.com/salva/p5-Language-Prolog-Yaswi). **PL_version()** is provided as a macro for compatibility. Defined keys are:

**PL_VERSION_SYSTEM**  
SWI-Prolog version as `10,000 × major + 100 × minor + patch`.

**PL_VERSION_FLI**  
Incremented if the foreign interface defined in this chapter changes in a way that breaks backward compatibility.

**PL_VERSION_REC**  
Incremented if the binary representation of terms as used by [PL_record_external()](foreigninclude.html#PL_record_external()) and [fast_write/2](IO.html#fast_write/2) changes.

**PL_VERSION_QLF**  
Incremented if the QLF file format changes.

**PL_VERSION_QLF_LOAD**  
Represents the oldest loadable QLF file format version.

**PL_VERSION_VM**  
A hash that represents the VM instructions and their arguments.

**PL_VERSION_BUILT_IN**  
A hash that represents the names, arities and properties of all built-in predicates defined in C. If this function is called before [PL_initialise()](foreigninclude.html#PL_initialise()) it returns 0.

### 12.4.21 Querying Prolog

`long` **PL_query**(`int`)  
Obtain status information on the Prolog system. The actual argument type depends on the information required. `int` describes what information is wanted.^(237Returning pointers and integers as a long is bad style. The signature of this function should be changed.) The options are given in [table 9](foreigninclude.html#tab:query).

> |  |  |
> |----|----|
> | `PL_QUERY_ARGC` | Return an integer holding the number of arguments given to Prolog from Unix. |
> | `PL_QUERY_ARGV` | Return a `char **` holding the argument vector given to Prolog from Unix. |
> | `PL_QUERY_SYMBOLFILE` | Return a `char *` holding the current symbol file of the running process. |
> | `PL_MAX_INTEGER` | Return a long, representing the maximal integer value represented by a Prolog integer. |
> | `PL_MIN_INTEGER` | Return a long, representing the minimal integer value. |
> | `PL_QUERY_VERSION` | Return a long, representing the version as `10,000 × M + 100 × m + p`, where `M` is the major, `m` the minor version number and `p` the patch level. For example, `20717` means `2.7.17`. |
> | `PL_QUERY_ENCODING` | Return the default stream encoding of Prolog (of type `IOENC`). |
> | `PL_QUERY_USER_CPU` | Get amount of user CPU time of the process in milliseconds. |

**Table 9 :** [PL_query()](foreigninclude.html#PL_query()) options

### 12.4.22 Registering Foreign Predicates

`bool` **PL_register_foreign_in_module**(`char *mod, char *name, int arity, foreign_t (*f)(), int flags, ...`)  
Register the C function `f` to implement a Prolog predicate. After this call returns successfully a predicate with name `name` (a `char *`) and arity `arity` (a C `int`) is created in module `mod`. If `mod` is `NULL`, the predicate is created in the module of the calling context, or if no context is present in the module `user`.

When called in Prolog, Prolog will call `function`. `flags` form a bitwise *or*’ed list of options for the installation. These are:

|  |  |
|----|----|
| `PL_FA_META` | Provide meta-predicate info (see below) |
| `PL_FA_TRANSPARENT` | Predicate is module transparent (deprecated) |
| `PL_FA_NONDETERMINISTIC` | Predicate is non-deterministic. See also [PL_retry()](foreigninclude.html#PL_retry()). |
| `PL_FA_NOTRACE` | Predicate cannot be seen in the tracer |
| `PL_FA_VARARGS` | Use alternative calling convention. |

If `PL_FA_META` is provided, [PL_register_foreign_in_module()](foreigninclude.html#PL_register_foreign_in_module()) takes one extra argument. This argument is of type `const char*`. This string must be exactly as long as the number of arguments of the predicate and filled with characters from the set `0-9:^-+?`. See [meta_predicate/1](metapred.html#meta_predicate/1) for details. `PL_FA_TRANSPARENT` is implied if at least one meta-argument is provided (`0-9:^`). Note that meta-arguments are *not always* passed as \<`module`\>:\<`term`\>. Always use [PL_strip_module()](foreigninclude.html#PL_strip_module()) to extract the module and plain term from a meta-argument.^(238It is encouraged to pass an additional `NULL` pointer for non-meta-predicates.)

Predicates may be registered either before or after [PL_initialise()](foreigninclude.html#PL_initialise()). When registered before initialisation the registration is recorded and executed after installing the system predicates and before loading the saved state.

Default calling (i.e., without `PL_FA_VARARGS`) `function` is passed the same number of `term_t` arguments as the arity of the predicate and, if the predicate is non-deterministic, an extra argument of type `control_t` (see [section 12.4.1.1](foreigninclude.html#sec:12.4.1.1)). If `PL_FA_VARARGS` is provided, `function` is called with three arguments. The first argument is a `term_t` handle to the first argument. Further arguments can be reached by adding the offset (see also [PL_new_term_refs()](foreigntypes.html#PL_new_term_refs())). The second argument is the arity, which defines the number of valid term references in the argument vector. The last argument is used for non-deterministic calls. It is currently undocumented and should be defined of type `void*`. Here is an example:

``` code
static foreign_t
atom_checksum(term_t a0, int arity, void* context)
{ char *s;

  if ( PL_get_atom_chars(a0, &s) )
  { int sum;

    for(sum=0; *s; s++)
      sum += *s&0xff;

    return PL_unify_integer(a0+1, sum&0xff);
  }

  return FALSE;
}

install_t
install(void)
{ PL_register_foreign("atom_checksum", 2,
                      atom_checksum, PL_FA_VARARGS);
}
```

Note that the `function` is documented as `foreign_t (*function)()`. The actual prototype uses `void*` as C11 and later do not allow for function types with unspecified arguments.

`bool` **PL_register_foreign**(`const char *name, int arity, foreign_t (*function)(), int flags, ...`)  
Same as [PL_register_foreign_in_module()](foreigninclude.html#PL_register_foreign_in_module()), passing `NULL` for the `module`.

`void` **PL_register_extensions_in_module**(`const char *module, PL_extension *e`)  
Register a series of predicates from an array of definitions of the type `PL_extension` in the given `module`. If `module` is `NULL`, the predicate is created in the module of the calling context, or if no context is present in the module `user`. The `PL_extension` type is defined as

``` code
typedef struct PL_extension
{ char          *predicate_name; /* Name of the predicate */
  short         arity;           /* Arity of the predicate */
  pl_function_t function;        /* Implementing functions */
  short         flags;           /* Or of PL_FA_... */
} PL_extension;
```

For details, see [PL_register_foreign_in_module()](foreigninclude.html#PL_register_foreign_in_module()). Here is an example of its usage:

``` code
static PL_extension predicates[] = {
{ "foo",        1,      pl_foo, 0 },
{ "bar",        2,      pl_bar, PL_FA_NONDETERMINISTIC },
{ NULL,         0,      NULL,   0 }
};

main(int argc, char **argv)
{ PL_register_extensions_in_module("user", predicates);

  if ( !PL_initialise(argc, argv) )
    PL_halt(1);

  ...
}
```

`void` **PL_register_extensions**(` PL_extension *e`)  
Same as [PL_register_extensions_in_module()](foreigninclude.html#PL_register_extensions_in_module()) using `NULL` for the `module` argument.

### 12.4.23 Foreign Code Hooks

For various specific applications some hooks are provided.

`PL_dispatch_hook_t` **PL_dispatch_hook**(`PL_dispatch_hook_t`)  
If this hook is not NULL, this function is called when reading from the terminal. It is supposed to dispatch events when SWI-Prolog is connected to a window environment. It can return two values: `PL_DISPATCH_INPUT` indicates Prolog input is available on file descriptor 0 or `PL_DISPATCH_TIMEOUT` to indicate a timeout. The old hook is returned. The type `PL_dispatch_hook_t` is defined as:

``` code
typedef int  (*PL_dispatch_hook_t)(void);
```

`void` **PL_abort_hook**(`PL_abort_hook_t`)  
Install a hook when [abort/0](toplevel.html#abort/0) is executed. SWI-Prolog [abort/0](toplevel.html#abort/0) is implemented using C **setjmp()**/**longjmp()** construct. The hooks are executed in the reverse order of their registration after the **longjmp()** took place and before the Prolog top level is reinvoked. The type `PL_abort_hook_t` is defined as:

``` code
typedef void (*PL_abort_hook_t)(void);
```

`bool` **PL_abort_unhook**(`PL_abort_hook_t`)  
Remove a hook installed with [PL_abort_hook()](foreigninclude.html#PL_abort_hook()). Returns `FALSE` if no such hook is found, `TRUE` otherwise.

`void` **PL_on_halt**(`int (*f)(int, void *), void *closure`)  
Register the function `f` to be called if SWI-Prolog is halted. The function is called with two arguments: the exit code of the process (0 if this cannot be determined) and the `closure` argument passed to the [PL_on_halt()](foreigninclude.html#PL_on_halt()) call. Handlers *must* return 0. Other return values are reserved for future use. See also [at_halt/1](consulting.html#at_halt/1).^(bugAlthough both [PL_on_halt()](foreigninclude.html#PL_on_halt()) and [at_halt/1](consulting.html#at_halt/1) are called in FIFO order, *all* [at_halt/1](consulting.html#at_halt/1) handlers are called before *all* [PL_on_halt()](foreigninclude.html#PL_on_halt()) handlers.) These handlers are called *before* system cleanup and can therefore access all normal Prolog resources. See also [PL_exit_hook()](foreigninclude.html#PL_exit_hook()).

`void` **PL_exit_hook**(`int (*f)(int, void *), void *closure`)  
Similar to [PL_on_halt()](foreigninclude.html#PL_on_halt()), but the hooks are executed by [PL_halt()](foreigninclude.html#PL_halt()) instead of [PL_cleanup()](foreigninclude.html#PL_cleanup()) just before calling **exit()**.

`PL_agc_hook_t` **PL_agc_hook**(`PL_agc_hook_t new`)  
Register a hook with the atom-garbage collector (see [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0)) that is called on any atom that is reclaimed. The old hook is returned. If no hook is currently defined, `NULL` is returned. The argument of the called hook is the atom that is to be garbage collected. The return value is an `int`. If the return value is zero, the atom is **not** reclaimed. The hook may invoke any Prolog predicate.

The example below defines a foreign library for printing the garbage collected atoms for debugging purposes.

``` code
#include <SWI-Stream.h>
#include <SWI-Prolog.h>

static int
atom_hook(atom_t a)
{ Sdprintf("AGC: deleting %s\n", PL_atom_chars(a));

  return TRUE;
}

static PL_agc_hook_t old;

install_t
install()
{ old = PL_agc_hook(atom_hook);
}

install_t
uninstall()
{ PL_agc_hook(old);
}
```

### 12.4.24 Storing foreign data

When combining foreign code with Prolog, it can be necessary to make data represented in the foreign language available to Prolog. For example, to pass it to another foreign function. At the end of this section, there is a partial implementation of using foreign functions to manage bit-vectors. Another example is the SGML/XML library that manages a‘parser’object, an object that represents the current state of the parser and that can be directed to perform actions such as parsing a document or make queries about the document content.

This section provides some hints for handling foreign data in Prolog. There are four options for storing such data:

- *Natural Prolog data*  
  Uses the representation one would choose if no foreign interface was required. For example, a bitvector representing a list of small integers can be represented as a Prolog list of integers.
- *Opaque packed data on the stacks*  
  It is possible to represent the raw binary representation of the foreign object as a Prolog string (see [section 5.2](string.html#sec:5.2)). Strings may be created from foreign data using [PL_put_string_nchars()](foreigninclude.html#PL_put_string_nchars()) and retrieved using [PL_get_string_chars()](foreigninclude.html#PL_get_string_chars()). It is good practice to wrap the string in a compound term with arity 1, so Prolog can identify the type. The hook [portray/1](termrw.html#portray/1) rules may be used to streamline printing such terms during development.
- *Opaque packed data in a blob*  
  Similar to the above solution, binary data can be stored in an atom. The blob interface ([section 12.4.10](foreigninclude.html#sec:12.4.10)) provides additional facilities to assign a type and hook-functions that act on creation and destruction of the underlying atom.
- *Natural foreign data, passed as a pointer*  
  An alternative is to pass a pointer to the foreign data. Again, the pointer is often wrapped in a compound term.

The choice may be guided using the following distinctions

- *Is the data opaque to Prolog*  
  With‘opaque’data, we refer to data handled in foreign functions, passed around in Prolog, but where Prolog never examines the contents of the data itself. If the data is opaque to Prolog, the selection will be driven solely by simplicity of the interface and performance.
- *What is the lifetime of the data*  
  With‘lifetime’we refer to how it is decided that the object is (or can be) destroyed. We can distinguish three cases:
  1.  The object must be destroyed on backtracking and normal Prolog garbage collection (i.e., it acts as a normal Prolog term). In this case, representing the object as a Prolog string (second option above) is the only feasible solution.
  2.  The data must survive Prolog backtracking. This leaves two options. One is to represent the object using a pointer and use explicit creation and destruction, making the programmer responsible. The alternative is to use the blob-interface, leaving destruction to the (atom) garbage collector.
  3.  The data lives as during the lifetime of a foreign function that implements a predicate. If the predicate is deterministic, foreign automatic variables are suitable. If the predicate is non-deterministic, the data may be allocated using **malloc()** and a pointer may be passed. See [section 12.4.1.1](foreigninclude.html#sec:12.4.1.1).

#### 12.4.24.1 Examples for storing foreign data

In this section, we outline some examples, covering typical cases. In the first example, we will deal with extending Prolog's data representation with integer sets, represented as bit-vectors. Then, we discuss the outline of the DDE interface.

**Integer sets** with not-too-far-apart upper- and lower-bounds can be represented using bit-vectors. Common set operations, such as union, intersection, etc., are reduced to simple *and*’ing and *or*’ing the bit-vectors. This can be done using Prolog's unbounded integers.

For really demanding applications, foreign representation will perform better, especially time-wise. Bit-vectors are naturally expressed using string objects. If the string is wrapped in `bitvector/1`, the lower-bound of the vector is 0 and the upper-bound is not defined; an implementation for getting and putting the sets as well as the union predicate for it is below.

``` code
#include <SWI-Prolog.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

static functor_t FUNCTOR_bitvector1;

static int
get_bitvector(term_t in, int *len, unsigned char **data)
{ if ( PL_is_functor(in, FUNCTOR_bitvector1) )
  { term_t a = PL_new_term_ref();

    PL_get_arg(1, in, a);
    return PL_get_string(a, (char **)data, len);
  }

  PL_fail;
}

static int
unify_bitvector(term_t out, int len, const unsigned char *data)
{ if ( PL_unify_functor(out, FUNCTOR_bitvector1) )
  { term_t a = PL_new_term_ref();

    PL_get_arg(1, out, a);

    return PL_unify_string_nchars(a, len, (const char *)data);
  }

  PL_fail;
}

static foreign_t
pl_bitvector_union(term_t t1, term_t t2, term_t u)
{ unsigned char *s1, *s2;
  int l1, l2;

  if ( get_bitvector(t1, &l1, &s1) &&
       get_bitvector(t2, &l2, &s2) )
  { int l = max(l1, l2);
    unsigned char *s3 = alloca(l);

    if ( s3 )
    { int n;
      int ml = min(l1, l2);

      for(n=0; n<ml; n++)
        s3[n] = s1[n] | s2[n];
      for( ; n < l1; n++)
        s3[n] = s1[n];
      for( ; n < l2; n++)
        s3[n] = s2[n];

      return unify_bitvector(u, l, s3);
    }

    return PL_warning("Not enough memory");
  }

  PL_fail;
}

install_t
install()
{ PL_register_foreign("bitvector_union", 3, pl_bitvector_union, 0);

  FUNCTOR_bitvector1 = PL_new_functor(PL_new_atom("bitvector"), 1);
}
```

**The DDE interface** (see [section 4.44](DDE.html#sec:4.44)) represents another common usage of the foreign interface: providing communication to new operating system features. The DDE interface requires knowledge about active DDE server and client channels. These channels contains various foreign data types. Such an interface is normally achieved using an open/close protocol that creates and destroys a *handle*. The handle is a reference to a foreign data structure containing the relevant information.

There are a couple of possibilities for representing the handle. The choice depends on responsibilities and debugging facilities. The simplest approach is to use [PL_unify_pointer()](foreigninclude.html#PL_unify_pointer()) and [PL_get_pointer()](foreigninclude.html#PL_get_pointer()). This approach is fast and easy, but has the drawbacks of (untyped) pointers: there is no reliable way to detect the validity of the pointer, nor to verify that it is pointing to a structure of the desired type. The pointer may be wrapped into a compound term with arity 1 (i.e., `dde_channel(<``Pointer``>)`), making the type-problem less serious.

Alternatively (used in the DDE interface), the interface code can maintain a (preferably variable length) array of pointers and return the index in this array. This provides better protection. Especially for debugging purposes, wrapping the handle in a compound is a good suggestion.

### 12.4.25 Embedding SWI-Prolog in other applications

With embedded Prolog we refer to the situation where the‘main’program is not the Prolog application. Prolog is sometimes embedded in C, C++, Java or other languages to provide logic based services in a larger application. Embedding loads the Prolog engine as a library to the external language. Prolog itself only provides for embedding in the C language (compatible with C++). Embedding in Java is achieved using JPL using a C-glue between the Java and Prolog C interfaces.

The most simple embedded program is below. The interface function [PL_initialise()](foreigninclude.html#PL_initialise()) **must** be called before any of the other SWI-Prolog foreign language functions described in this chapter, except for **PL_initialise_hook()**, [PL_new_atom()](foreigninclude.html#PL_new_atom()), [PL_new_functor()](foreigninclude.html#PL_new_functor()) and [PL_register_foreign()](foreigninclude.html#PL_register_foreign()). [PL_initialise()](foreigninclude.html#PL_initialise()) interprets all the command line arguments, except for the **-t** `toplevel` flag that is interpreted by [PL_toplevel()](foreigninclude.html#PL_toplevel()).

``` code
int
main(int argc, char **argv)
{ if ( !PL_initialise(argc, argv) )
    PL_halt(1);

  PL_halt(PL_toplevel() ? 0 : 1);
}
```

`bool` **PL_initialise**(`int argc, char **argv`)  
Initialises the SWI-Prolog heap and stacks, restores the Prolog state, loads the system and personal initialisation files, runs the [initialization/1](consulting.html#initialization/1) hooks and finally runs the initialization goals registered using **-g** `goal`.

Special consideration is required for `argv[0]`. On **Unix**, this argument passes the part of the command line that is used to locate the executable. Prolog uses this to find the file holding the running executable. The **Windows** version uses this to find a *module* of the running executable. If the specified module cannot be found, it tries the module `libswipl.dll`, containing the Prolog runtime kernel. In all these cases, the resulting file is used for two purposes:

- See whether a Prolog saved state is appended to the file. If this is the case, this state will be loaded instead of the default `boot.prc` file from the SWI-Prolog home directory. See also [qsave_program/\[1,2\]](saved-states.html#qsave_program/1) and [section 12.5](plld.html#sec:12.5).
- Find the Prolog home directory. This process is described in detail in [section 12.6](findhome.html#sec:12.6).

[PL_initialise()](foreigninclude.html#PL_initialise()) returns 1 if all initialisation succeeded and 0 otherwise.^(bugVarious fatal errors may cause [PL_initialise()](foreigninclude.html#PL_initialise()) to call [PL_halt(1)](foreigninclude.html#PL_halt()), preventing it from returning at all.)

In most cases, `argc` and `argv` will be passed from the main program. It is allowed to create your own argument vector, provided `argv[0]` is constructed according to the rules above. For example:

``` code
int
main(int argc, char **argv)
{ char *av[10];
  int ac = 0;

  av[ac++] = argv[0];
  av[ac++] = "-x";
  av[ac++] = "mystate";
  av[ac]   = NULL;

  if ( !PL_initialise(ac, av) )
    PL_halt(1);
  ...
}
```

Please note that the passed argument vector may be referred from Prolog at any time and should therefore be valid as long as the Prolog engine is used.

A good setup in Windows is to add SWI-Prolog's `bin` directory to your `PATH` and either pass a module holding a saved state, or `"libswipl.dll"` as `argv[0]`. If the Prolog state is attached to a DLL (see the **-dll** option of **swipl-ld**), pass the name of this DLL.

`bool` **PL_winitialise**(`int argc, wchar_t **argv`)  
Wide character version of [PL_initialise()](foreigninclude.html#PL_initialise()). Can be used in Windows combined with the **wmain()** entry point.

`bool` **PL_is_initialised**(`int *argc, char ***argv`)  
Test whether the Prolog engine is already initialised. Returns `FALSE` if Prolog is not initialised and `TRUE` otherwise. If the engine is initialised and `argc` is not `NULL`, the argument count used with [PL_initialise()](foreigninclude.html#PL_initialise()) is stored in `argc`. Same for the argument vector `argv`.

`bool` **PL_set_resource_db_mem**(`const unsigned char *data, size_t size`)  
This function must be called at most once and *before* calling [PL_initialise()](foreigninclude.html#PL_initialise()). The memory area designated by `data` and `size` must contain the resource data and be in the format as produced by [qsave_program/2](saved-states.html#qsave_program/2). The memory area is accessed by [PL_initialise()](foreigninclude.html#PL_initialise()) as well as calls to [open_resource/3](program-resources.html#open_resource/3).^(239This implies that the data must remain accessible during the lifetime of the process if [open_resource/3](program-resources.html#open_resource/3) is used. Future versions may provide a function to detach the resource database and cause [open_resource/3](program-resources.html#open_resource/3) to raise an exception.)

For example, we can include the bootstrap data into an embedded executable using the steps below. The advantage of this approach is that it is fully supported by any OS and you obtain a single file executable.

1.  Create a saved state using [qsave_program/2](saved-states.html#qsave_program/2) or

    ``` code
    % swipl -o state -c file.pl ...
    ```

2.  Create a C source file from the state using e.g., the Unix utility **xxd**(1):

    ``` code
    % xxd -i state > state.h
    ```

3.  Embed Prolog as in the example below. Instead of calling the toplevel you probably want to call your application code.

    ``` code
    #include <SWI-Prolog.h>
    #include "state.h"

    int
    main(int argc, char **argv)
    { if ( !PL_set_resource_db_mem(state, state_len) ||
           !PL_initialise(argc, argv) )
        PL_halt(1);

      return PL_toplevel();
    }
    ```

Alternative to **xxd**, it is possible to use inline assembler, e.g. the **gcc** `incbin` instruction. Code for **gcc** was provided by Roberto Bagnara on the SWI-Prolog mailinglist. Given the state in a file `state`, create the following assembler program:

``` code
        .globl _state
        .globl _state_end
_state:
        .incbin "state"
_state_end:
```

Now include this as follows:

``` code
#include <SWI-Prolog.h>

#if __linux
#define STATE _state
#define STATE_END _state_end
#else
#define STATE state
#define STATE_END state_end
#endif

extern unsigned char STATE[];
extern unsigned char STATE_END[];

int
main(int argc, char **argv)
{ if ( !PL_set_resource_db_mem(STATE, STATE_END - STATE) ||
       !PL_initialise(argc, argv) )
    PL_halt(1);
  return PL_toplevel();
}
```

As Jose Morales pointed at [https://github.com/graphitemaster/incbin](https://github.com/graphitemaster/incbin), which contains a portability layer on top of the above idea.

`bool` **PL_toplevel**()  
Runs the goal of the **-t** `toplevel` switch (default [prolog/0](toplevel.html#prolog/0)) and returns 1 if successful, 0 otherwise.

`int` **PL_cleanup**(`int status_and_flags`)  
This function may be called instead of [PL_halt()](foreigninclude.html#PL_halt()) to cleanup Prolog without exiting the process. It performs the reverse of [PL_initialise()](foreigninclude.html#PL_initialise()). It runs the [PL_on_halt()](foreigninclude.html#PL_on_halt()) and [at_halt/1](consulting.html#at_halt/1) handlers, closes all streams (except for the *standard I/O* streams, which are flushed only), restores all signal handlers and reclaims all memory unless asked not to. `status_and_flags` accepts the following flags:

**`PL_CLEANUP_NO_RECLAIM_MEMORY`**  
Do not reclaim memory. This is the default when called from [PL_halt()](foreigninclude.html#PL_halt()) for the release versions because the OS will do so anyway.

**`PL_CLEANUP_NO_CANCEL`**  
Do not allow Prolog and foreign *halt* hooks to cancel the cleanup.

The return value of [PL_cleanup()](foreigninclude.html#PL_cleanup()) is one of the following:

**`PL_CLEANUP_CANCELED`**  
A Prolog or foreign *halt* hook cancelled the cleanup. Note that some of the halt hooks may have been executed.

**`PL_CLEANUP_SUCCESS`**  
Cleanup completed successfully. Unless `PL_CLEANUP_NO_RECLAIM_MEMORY` was specified this implies most of the memory was reclaimed and Prolog may be reinitialized in the same process using [PL_initialise()](foreigninclude.html#PL_initialise()).

**`PL_CLEANUP_FAILED`**  
Cleanup failed. This happens if the user requested to reclaim all memory but this failed because the system was not able to *join* all Prolog threads and could therefore not reclaim the memory.

**`PL_CLEANUP_RECURSIVE`**  
[PL_cleanup()](foreigninclude.html#PL_cleanup()) was called from a hook called by the cleanup process.

[PL_cleanup()](foreigninclude.html#PL_cleanup()) allows deleting and restarting the Prolog system in the same process. In versions older than 8.5.9 this did not work. As of version 8.5.9, it works for the basic Prolog engine. Many of the plugins that contain foreign code do not implement a suitable *uninstall* handler and will leak memory and possibly other resources. Note that shutting Prolog down and renitializing it is slow. For almost all scenarios there are faster alternatives such as reloading modified code, using *temporary modules*, using *threads*, etc.

`void` **PL_cleanup_fork**()  
Stop intervaltimer that may be running on behalf of [profile/1](profile.html#profile/1). The call is intended to be used in combination with **fork()**:

``` code
    if ( (pid=fork()) == 0 )
    { PL_cleanup_fork();
      <some exec variation>
    }
```

The call behaves the same on Windows, though there is probably no meaningful application.

`bool` **PL_halt**(`int status`)  
Terminate the Prolog process as [halt/1](toplevel.html#halt/1). The `status` is used to set the return code from calling **`exit()`** - this is a value between 0 and 65535. Additional behaviour can be specified by OR-ing the status with these flags:

**PL_HALT_WITH_EXCEPTION**  
When specified, try to raise `unwind(``halt(status)``)`. If this is not possible, ignore this flag. Return `false`.

**PL_CLEANUP_NO_RECLAIM_MEMORY**  
Never reclaim the memory, leaving this task to the OS. This is the default unless the system is compiled for debugging.

**PL_CLEANUP_NO_CANCEL**  
Do not allow hooks to cancel halting the system.

Unless `PL_HALT_WITH_EXCEPTION` was specified and effective, Clean up the Prolog environment using [PL_cleanup()](foreigninclude.html#PL_cleanup()) and if successful call **exit()** with the status argument. Returns `true` if exit was cancelled by [PL_cleanup()](foreigninclude.html#PL_cleanup()).^(240Versions up to 9.3.12 returned `false` on cancel. If necessary, the two behaviours can be distinguished based on the existence of the `PL_HALT_WITH_EXCEPTION` macro.)

#### 12.4.25.1 Threading, Signals and embedded Prolog

This section applies to Unix-based environments that have signals or multithreading. The Windows version is compiled for multithreading, and Windows lacks proper signals.

We can distinguish two classes of embedded executables. There are small C/C++ programs that act as an interfacing layer around Prolog. Most of these programs can be replaced using the normal Prolog executable extended with a dynamically loaded foreign extension and in most cases this is the preferred route. In other cases, Prolog is embedded in a complex application that---like Prolog---wants to control the process environment. A good example is Java. Embedding Prolog is generally the only way to get these environments together in one process image. Java VMs, however, are by nature multithreaded and appear to do signal handling (software interrupts).

On Unix systems, SWI-Prolog installs handlers for the following signals:

**SIGUSR2**  
has an empty signal handler. This signal is sent to a thread after sending a thread-signal (see [thread_signal/2](threadcom.html#thread_signal/2)). It causes blocking system calls to return with `EINTR`, which gives them the opportunity to react to thread-signals.

In some cases the embedded system and SWI-Prolog may both use `SIGUSR2` without conflict. If the embedded system redefines `SIGUSR2` with a handler that runs quickly and no harm is done in the embedded system due to spurious wakeup when initiated from Prolog, there is no problem. If SWI-Prolog is initialised *after* the embedded system it will call the handler set by the embedded system and the same conditions as above apply. SWI-Prolog's handler is a simple function only chaining a possibly previously registered handler. SWI-Prolog can handle spurious `SIGUSR2` signals.

**SIGINT**  
is used by the top level to activate the tracer (typically bound to control-C). The first control-C posts a request for starting the tracer in a safe, synchronous fashion. If control-C is hit again before the safe route is executed, it prompts the user whether or not a forced interrupt is desired.

**SIGTERM, SIGABRT and SIGQUIT**  
are caught to cleanup before killing the process again using the same signal.

**SIGSEGV, SIGILL, SIGBUS, SIGFPE and SIGSYS**  
are caught by to print a backtrace before killing the process again using the same signal.

**SIGHUP**  
is caught and causes the process to exit with status 2 after cleanup.

The **--no-signals** option can be used to inhibit all signal processing except for `SIGUSR2`. The handling of `SIGUSR2` is vital for dealing with blocking system call in threads. The used signal may be changed using the **--sigalert=NUM** option or disabled using `--sigalert=0`.
