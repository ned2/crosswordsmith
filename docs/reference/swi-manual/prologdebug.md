
## A.41 library(prolog_debug): User level debugging tools

This library provides tools to control the Prolog debuggers. Traditionally this code was built-in. Because these tools are only required in (interactive) debugging sessions they have been moved into the library.

\[multifile\]prolog:**debug_control_hook**(`+Action`)  
Allow user-hooks in the Prolog debugger interaction. See the calls below for the provided hooks. We use a single predicate with action argument to avoid an uncontrolled poliferation of hooks.

\[det\]**spy**(`:Spec`)  
\[det\]**nospy**(`:Spec`)  
\[det\]**nospyall**  
Set/clear spy-points. A successfully set or cleared spy-point is reported using [print_message/2](printmsg.html#print_message/2), level `informational`, with one of the following terms, where `Spec` is of the form M:Head.

- `spy(Spec)`
- `nospy(Spec)`

See also  
[spy/1](debugger.html#spy/1) and [nospy/1](debugger.html#nospy/1) call the hook [prolog:debug_control_hook/1](prologdebug.html#prolog:debug_control_hook/1) to allow for alternative specifications of the thing to debug.

\[det\]**debugging**  
Report current status of the debugger.

\[multifile\]**debugging_hook**(`+DebugMode`)  
Multifile hook that is called as `forall(debugging_hook(DebugMode), true)` and that may be used to extend the information printed from other debugging libraries.

\[det\]**trap**(`+Formal`)  
\[det\]**notrap**(`+Formal`)  
Install a trap on `error(Formal, Context)` exceptions that unify. The tracer is started when a matching exception is raised. This predicate enables *debug mode* using [debug/0](debugger.html#debug/0) to get more context about the exception. Even with debug mode disabled exceptions are still trapped and thus one may call [nodebug/0](debugger.html#nodebug/0) to run in normal mode after installing a trap. Exceptions are trapped in any thread. Debug mode is only enabled in the calling thread. To enable debug mode in all threads use [tdebug/0](threadutil.html#tdebug/0).

Calling [debugging/0](debugger.html#debugging/0) lists the enabled traps. The predicate [notrap/1](prologdebug.html#notrap/1) removes matching (unifying) traps.

In many cases debugging an exception that is caught is as simple as below (assuming run/0 starts your program).

``` code
?- trap(_).
?- run.
```

The multifile hook [trap_alias/2](prologdebug.html#trap_alias/2) allow for defining short hands for commonly used traps. Currently this defines

**det**  
Trap determinism exceptions raised as a result of the [det/1](debug-determinism.html#det/1) directive.

**`=>`**  
Trap rule existence error exceptions.

See also  
\- gtrap/1 to trap using the graphical debugger.  
- *Edit exceptions* menu in PceEmacs and the graphical debugger that provide a graphical frontend to trap exceptions.

\[multifile\]**trap_alias**(`+Alias, -Error`)  
Define short hands for commonly used exceptions.

\[failure\]**exception_hook**(`+ExIn, -ExOut, +Frame, +Catcher, +DebugMode`)  
Trap exceptions and consider whether or not to start the tracer.
