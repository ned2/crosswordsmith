
## B.5 Breakpoint and watchpoint handling

SWI-Prolog support *breakpoints*. Breakpoints can be manipulated with the library `library(prolog_breakpoints)`. Setting a breakpoint replaces a virtual machine instruction with the `D_BREAK` instruction. If the virtual machine executes a `D_BREAK`, it performs a callback to decide on the action to perform. This section describes this callback, called [prolog:break_hook/7](breakpoint.html#prolog:break_hook/7).

\[hook,semidet\]**prolog:break_hook**(`+Clause, +PC, +FR, +BFR, +Expression, +Debug, -Action`)  
*Experimental* This hook is called if the virtual machine executes a `D_BREAK`, set using set_breakpoint/4. `Clause` and `PC` identify the breakpoint. `FR` and `BFR` provide the environment frame and current choicepoint. `Debug` is `true` if the system was in debug mode when the breakpoint was reached, otherwise `Debug` is `false`. `Expression` identifies the action that is interrupted, and is one of the following:

**call**(`Goal`)  
The instruction will call `Goal`. This is generated for nearly all instructions. Note that `Goal` is semantically equivalent to the compiled body term, but might differ syntactically. This is notably the case when arithmetic expressions are compiled in optimized mode (see [optimise](flags.html#flag:optimise)). In particular, the arguments of arithmetic expressions have already been evaluated. Thus, `A` is 3\*`B`, where `B` equals 3 results in a term `call(A is 9)` if the clause was compiled with optimization enabled.

**`!`**  
The instruction will call the cut. Because the semantics of metacalling the cut differs from executing the cut in its original context we do not wrap the cut in `call/1`.

**`:-`**  
The breakpoint is on the *neck* instruction, i.e., after performing the head unifications.

**exit**  
The breakpoint is on the *exit* instruction, i.e., at the end of the clause. Note that the exit instruction may not be reached due to last-call optimisation.

**unify_exit**  
The breakpoint is on the completion of an in-lined unification while the system is not in debug mode. If the system is in debug mode, inlined unification is returned as call(Var=Term).^(255This hack will disappear if we find a good solution for applying D_BREAK to inlined unification. Only option might be to place the break on both the unification start and end instructions.)

If [prolog:break_hook/7](breakpoint.html#prolog:break_hook/7) succeeds, it must unify `Action` with a value that describes how execution must continue. Possible values for `Action` are:

**continue**  
Just continue as if no breakpoint was present.

**debug**  
Continue in *debug mode*. See [debug/0](debugger.html#debug/0).

**trace**  
Continue in *trace mode*. See [trace/0](debugger.html#trace/0).

**call**(`Goal`)  
Execute `Goal` instead of the goal that would be executed. `Goal` is executed as [call/1](metacall.html#call/1), preserving (non-)determinism and exceptions.

If this hook throws an exception, the exception is propagated normally. If this hook is not defined or fails, the default action is executed. This implies that, if the thread is in debug mode, the tracer will be enabled (`trace`) and otherwise the breakpoint is ignored (`continue`).

This hook allows for injecting various debugging scenarios into the executable without recompiling. The hook can access variables of the calling context using the frame inspection predicates. Here are some examples.

- Create *conditional* breakpoints by imposing conditions before deciding the return `trace`.
- Watch variables at a specific point in the execution. Note that binding of these variables can be monitored using *attributed variables*, see [section 8.1](attvar.html#sec:8.1).
- Dynamically add *assertions* on variables using [assertion/1](debug.html#assertion/1).
- Wrap the `Goal` into a meta-call that traces progress of the `Goal`.
