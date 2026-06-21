
## 4.12 Handling signals

As of version 3.1.0, SWI-Prolog is able to handle software interrupts (signals) in Prolog as well as in foreign (C) code (see [section 12.4.17](foreigninclude.html#sec:12.4.17)).

Signals are used to handle internal errors (execution of a non-existing CPU instruction, arithmetic domain errors, illegal memory access, resource overflow, etc.), as well as for dealing with asynchronous interprocess communication.

Signals are defined by the POSIX standard and part of all Unix machines. The MS-Windows Win32 provides a subset of the signal handling routines, lacking the vital functionality to raise a signal in another thread for achieving asynchronous interprocess (or interthread) communication (Unix **kill()** function).

**on_signal**(`+Signal, -Old, :New`)  
Determines how `Signal` is processed. `Old` is unified with the old behaviour, while the behaviour is switched to `New`. As with similar environment control predicates, the current value is retrieved using `on_signal(Signal, Current, Current)`.

The action description is an atom denoting the name of the predicate that will be called if `Signal` arrives. [on_signal/3](signal.html#on_signal/3) is a meta-predicate, which implies that \<`Module`\>:\<`Name`\> refers to \<`Name`\>/1 in module \<`Module`\>. The handler is called with a single argument: the name of the signal as an atom. The Prolog names for signals are explained below.

Four names have special meaning. `throw` implies Prolog will map the signal onto a Prolog exception as described in [section 4.10](exception.html#sec:4.10), `ignore` causes Prolog to ignore the signal entirely, `debug` specifies the debug interrupt prompt that is initially bound to `SIGINT` (Control-C) and `default` resets the handler to the settings active before SWI-Prolog manipulated the handler.

Signals bound to a foreign function through [PL_signal()](foreigninclude.html#PL_signal()) are reported using the term `’$foreign_function’(Address)`.

After receiving a signal mapped to `throw`, the exception raised has the following structure:

> `error(signal(<``SigName``>, <``SigNum``>), <``Context``>) `

The signal names are defined by the POSIX standard as symbols of the form `SIG`\<`SIGNAME`\>. The Prolog name for a signal is the lowercase version of \<`SIGNAME`\>. The predicate [current_signal/3](signal.html#current_signal/3) may be used to map between names and signals.

Initially, the following signals are handled unless the command line option **--no-signals** is specified:

**int**  
Prompts the user, allowing to inspect the current state of the process and start the tracer.

**usr2**  
Bound to an empty signal handler used to make blocking system calls return. This allows [thread_signal/2](threadcom.html#thread_signal/2) to interrupt threads blocked in a system call. See also [prolog_alert_signal/2](signal.html#prolog_alert_signal/2).

**pipe**  
Ignored.

**hup, term, abrt, quit**  
Causes normal Prolog cleanup (e.g., [at_halt/1](consulting.html#at_halt/1)) before terminating the process with the same signal.

**segv, ill, bus, sys**  
Dumps the C and Prolog stacks and runs cleanup before terminating the process with the same signal.

**fpe, alrm, xcpu, xfsz, vtalrm**  
Throw a Prolog exception (see above).

**current_signal**(`?Name, ?Id, ?Handler`)  
Enumerate the currently defined signal handling. `Name` is the signal name, `Id` is the numerical identifier and `Handler` is the currently defined handler (see [on_signal/3](signal.html#on_signal/3)).

**prolog_alert_signal**(`?Old, +New`)  
Query or set the signal used to unblock blocking system calls on Unix(-like) systems and process pending Prolog signals. The default is `SIGUSR2`. See also **--sigalert**. `New` can be a signal name or number. See [on_signal/3](signal.html#on_signal/3) for how the Prolog signal name is defined. The `Old` argument is unified to the signal name if known and the number otherwise. Notably the value 0 (zero) indicates that the system does not use an alarm signal. On POSIX systems, this implies that system calls are not interrupted by [thread_signal/2](threadcom.html#thread_signal/2).

This predicate is only defined on systems where the alert signal mechanism is used.

### 4.12.1 Notes on signal handling

Before deciding to deal with signals in your application, please consider the following:

- *Portability*  
  On MS-Windows, the signal interface is severely limited. Different Unix brands support different sets of signals, and the relation between signal name and number may vary. Currently, the system only supports signals numbered 1 to 32^(83TBD: the system should support the Unix realtime signals). Installing a signal outside the limited set of supported signals in MS-Windows crashes the application.

- *Safety*  
  Immediately delivered signals (see below) are unsafe. This implies that foreign functions called from a handler cannot safely use the SWI-Prolog API and cannot use C **longjmp()**. Handlers defined as `throw` are unsafe. Handlers defined to call a predicate are safe. Note that the predicate can call [throw/1](exception.html#throw/1), but the delivery is delayed until Prolog is in a safe state.

  The C-interface described in [section 12.4.17](foreigninclude.html#sec:12.4.17) provides the option `PL_SIGSYNC` to select either safe synchronous or unsafe asynchronous delivery.

- *Time of delivery*  
  Using `throw` or a foreign handler, signals are delivered immediately (as defined by the OS). When using a Prolog predicate, delivery is delayed to a safe moment. Blocking system calls or foreign loops may cause long delays. Foreign code can improve on that by calling [PL_handle_signals()](foreigninclude.html#PL_handle_signals()).

  Signals are blocked when the garbage collector is active.
