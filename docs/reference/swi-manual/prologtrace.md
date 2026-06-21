
## A.43 library(prolog_trace): Print access to predicates

See also  
`library(debug)` for adding conditional print statements to a program.

This library prints accesses to specified predicates by wrapping the predicate.

\[det\]**trace**(`:Pred`)  
\[det\]**trace**(`:Pred, +PortSpec`)  
Print passes through *ports* of specified predicates. `Pred` is a, possible partial, specification of a predicate as it is also used be [spy/1](debugger.html#spy/1) and similar predicates. Where a full predicate specification is of the shape `Module:Name/Arity` (or‘`//`Arity for non-terminals), both the module and arity may be omitted in which case `Pred` refers to all matching predicates. `PortSpec` is either a single port (`call`, `exit`, `fail` or `redo`), preceded with `+` or `-` or a list of these. The predicate modifies the current trace specification and then installs a suitable wrapper for the predicate using wrap_predicate/4. For example:

``` code
?- trace(append).
%     lists:append/2: [all]
%     lists:append/3: [all]
%     append/1: [all]
true.

?- append([a,b], [c], L).
 T [10] Call: lists:append([a, b], [c], _18032)
 T [19] Call: lists:append([b], [c], _19410)
 T [28] Call: lists:append([], [c], _20400)
 T [28 +0.1ms] Exit: lists:append([], [c], [c])
 T [19 +0.2ms] Exit: lists:append([b], [c], [b, c])
 T [10 +0.5ms] Exit: lists:append([a, b], [c], [a, b, c])
L = [a, b, c].

?- trace(append, -all).
%     lists:append/2: Not tracing
%     lists:append/3: Not tracing
%     append/1: Not tracing
```

The text between `[]` indicates the call depth (first number) and for all ports except the `call` port the *wall* time since the start (call port) in milliseconds. Note that the instrumentation and print time is included in the time. In the example above the actual time is about 0.00001ms on todays hardware.

In addition, **conditions** may be specified. In this case the the specification takes the shape `trace(:Head, Port(Condition))`. For example:

``` code
?- trace(current_prolog_flag(Flag, Value), call(var(Flag))).
?- list_tracing.
% Trace points (see trace/1,2) on:
%     system:current_prolog_flag(A,_): [call(var(A))]
```

This specification will only print the goal if the registered condition succeeds. Note that we can use the condition for its side effect and then fail to avoid printing the event. Clearing the trace event on all relevant ports removes the condition. There is currently no way to modify the condition without clearing the trace point first.

**tracing**(`:Spec, -Ports`)  
True if `Spec` is traced using `Ports`. `Spec` is a fully qualified head term.

**list_tracing**  
List predicates we are currently tracing

\[det\]**notraceall**  
Remove all trace points
