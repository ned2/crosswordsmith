
## B.6 Adding context to errors: prolog_exception_hook

The hook [prolog:prolog_exception_hook/5](excepthook.html#prolog:prolog_exception_hook/5) has been introduced to provide dedicated exception handling facilities for application frameworks, for example non-interactive server applications that wish to provide extensive context for exceptions for offline debugging.

**prolog:prolog_exception_hook**()  
+ExceptionIn, -ExceptionOut, +Frame, +CatcherFrame, +DebugMode This hook predicate, if defined in the module `prolog`, is between raising an exception and handling it. It is intended to allow a program adding additional context to an exception to simplify diagnosing the problem. `ExceptionIn` is the exception term as raised by [throw/1](exception.html#throw/1) or one of the built-in predicates. The output argument `ExceptionOut` describes the exception that is actually raised. `Frame` is the innermost frame. See [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3) and the library `library(prolog_stack)` for getting information from this. `CatcherFrame` is a reference to the frame calling the matching [catch/3](exception.html#catch/3), `none` if the exception is not caught or `’C’` if the exception is caught in C calling Prolog using the flag `PL_Q_CATCH_EXCEPTION`. `DebugMode` contains the setting of the Prolog flag [debug](flags.html#flag:debug) from the calling context.

The hook is run in‘nodebug’mode. If it succeeds, `ExceptionOut` is considered the current exception. If it fails, `ExceptionIn` is used for further processing. The hook is *never* called recursively. The hook is *not* allowed to modify `ExceptionOut` in such a way that it no longer unifies with the catching frame.

Typically, [prolog:prolog_exception_hook/5](excepthook.html#prolog:prolog_exception_hook/5) is used to fill the second argument of `error(Formal, Context)` exceptions. `Formal` is defined by the ISO standard, while SWI-Prolog defines `Context` as a term `context(Location, Message)`. `Location` is bound to a term \<`name`\>/\<`arity`\> by the kernel. This hook can be used to add more information on the calling context, such as a full stack trace.

Applications that use exceptions as part of normal processing must do a quick test of the environment before starting expensive gathering information on the state of the program.

The hook can call [trace/0](debugger.html#trace/0) to enter trace mode immediately. For example, imagine an application performing an unwanted division by zero while all other errors are expected and handled. We can force the debugger using the hook definition below. Run the program in debug mode (see [debug/0](debugger.html#debug/0)) to preserve as much as possible of the error context.

``` code
user:prolog_exception_hook(
         error(evaluation_error(zero_divisor), _),
         _, _, _) :-
        trace, fail.
```

This hook is used by `library(prolog_stack)` to print stack traces on uncaught exceptions, [trap/1](prologdebug.html#trap/1) to debug after exceptions and the GUI exception editor that is part of the GUI debugger.
