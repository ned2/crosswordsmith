
## B.3 Intercepting the Tracer

**prolog_trace_interception**(`+Port, +Frame, +Choice, -Action`)  
Dynamic predicate, normally not defined. This predicate is called from the SWI-Prolog debugger just before it would show a port. If this predicate succeeds, the debugger assumes that the trace action has been taken care of and continues execution as described by `Action`. Otherwise the normal Prolog debugger actions are performed.

`Port` denotes the reason to activate the tracer (\`port’in the 4/5-port, but with some additions):

**call**  
Normal entry through the call port of the 4-port debugger.

**redo**(`PC`)  
Normal entry through the redo port of the 4-port debugger. The `redo` port signals resuming a predicate to generate alternative solutions. If `PC` is 0 (zero), clause indexing has found another clause that will be tried next. Otherwise, `PC` is the program counter in the current clause where execution continues. This implies we are dealing with an in-clause choice point left by, e.g., [;/2](control.html#;/2). Note that non-determinism in foreign predicates are also handled using an in-clause choice point.

**unify**  
The unify port represents the *neck* instruction, signalling the end of the head-matching process. This port is normally invisible. See [leash/1](debugger.html#leash/1) and [visible/1](debugger.html#visible/1).

**exit**  
The exit port signals the goal is proved. It is possible for the goal to have alternatives. See [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3) to examine the goal stack.

**fail**  
The fail port signals final failure of the goal.

**exception**(`Except`)  
An exception is raised and still pending. This port is activated on each parent frame of the frame generating the exception until the exception is caught or the user restarts normal computation using `retry`. `Except` is the pending exception term.

**cut_call**(`PC`)  
A cut is encountered at `PC`. This port is used by the graphical debugger to visualise the effect of the cut.

**cut_exit**(`PC`)  
A cut has been executed. See `cut_call(PC)` for more information.

`Frame` is a reference to the current local stack frame, which can be examined using [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3). `Choice` is a reference to the last choice point and can be examined using [prolog_choice_attribute/3](manipstack.html#prolog_choice_attribute/3). `Action` must be unified with a term that specifies how execution must continue. The following actions are defined:

**abort**  
Abort execution. See [abort/0](toplevel.html#abort/0).

**halt**  
Terminate execution. See [halt/0](toplevel.html#halt/0).

**continue**  
Continue (i.e., *creep* in the command line debugger).

**fail**  
Make the current goal fail.

**ignore**  
Step over the current goal without executing it.

**nodebug**  
Continue execution in normal nodebug mode. See [nodebug/0](debugger.html#nodebug/0).

**leap**  
Continue execution in normal debug mode. See [debug/0](debugger.html#debug/0).

**retry**  
Retry the current frame.

**retry**(`Frame`)  
Retry the given frame. This must be a parent of the current frame.

**skip**  
Skip over the current goal (i.e., *skip* in the command line debugger).

**skip**(`Frame`)  
Skip to the end the execution of `Frame`. This is used to implement *finish* on an arbitrary frame in the GUI debugger.

**up**  
Skip to the parent goal (i.e., *up* in the command line debugger). This is the same as `skip(Frame)` using the parent frame of the current frame.

Together with the predicates described in [section 4.39](debugger.html#sec:4.39) and the other predicates of this chapter, this predicate enables the Prolog user to define a complete new debugger in Prolog. Besides this, it enables the Prolog programmer to monitor the execution of a program. The example below records all goals trapped by the tracer in the database.

``` code
prolog_trace_interception(Port, Frame, _PC, continue) :-
        prolog_frame_attribute(Frame, goal, Goal),
        prolog_frame_attribute(Frame, level, Level),
        recordz(trace, trace(Port, Level, Goal)).
```

To trace the execution of‘go’this way the following query should be given:

``` code
?- trace, go, notrace.
```

As of version 9.1.12, unification against variables in the passed data as well as changes to backtrackable global variables persist. The hook should not unify variables in its arguments. One solution to this is to backtrace over the body of the interceptor. Note that the `Action` needs to be preserved.

``` code
user:prolog_trace_interception(Port, Frame, Choice, Action) :-
    State = state(0),
    (   my_trace_interception(Port, Frame, Choice, Action),
        nb_setarg(1, State, Action),
        fail
    ;   arg(1, State, Action)
    ).
```

**prolog_skip_level**(`-Old, +New`)  
Unify `Old` with the old value of‘skip level’and then set this level according to `New`. `New` is an integer, the atom `very_deep` (meaning don't skip) or the atom `skip_in_redo` (see prolog_skip_frame/1). The‘skip level’is a setting of each Prolog thread that disables the debugger on all recursion levels deeper than the level of the variable. See also prolog_skip_frame/1.
