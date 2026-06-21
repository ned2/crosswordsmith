
## B.1 Examining the Environment Stack

\[det\]**prolog_current_frame**(`-Frame`)  
Unify `Frame` with an integer providing a reference to the parent of the current local stack frame. A pointer to the current local frame cannot be provided as the predicate succeeds deterministically and therefore its frame is destroyed immediately after succeeding.

\[semidet\]**prolog_current_choice**(`-Choice`)  
Unify `Choice` with an integer provided a reference to the last choice point. Fails if the current environment has no choice points. See also [prolog_choice_attribute/3](manipstack.html#prolog_choice_attribute/3).

**prolog_frame_attribute**(`+Frame, +Key, :Value`)  
Obtain information about the local stack frame `Frame`. `Frame` is a frame reference as obtained through [prolog_current_frame/1](manipstack.html#prolog_current_frame/1), [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4) or this predicate. The key values are described below.

**alternative**  
`Value` is unified with an integer reference to the local stack frame in which execution is resumed if the goal associated with `Frame` fails. Fails if the frame has no alternative frame.

**has_alternatives**  
`Value` is unified with `true` if `Frame` still is a candidate for backtracking; `false` otherwise.

**goal**  
`Value` is unified with the goal associated with `Frame`. If the definition module of the active predicate is not the calling context, the goal is represented as `<``module``>:<``goal``>`. Do not instantiate variables in this goal unless you **know** what you are doing! Note that the returned term may contain references to the frame and should be discarded before the frame terminates.^(254The returned term is actually an illegal Prolog term that may hold references from the global to the local stack to preserve the variable names.)

**parent_goal**  
**parent_goal**(`-Parent`)  
If `Value` is instantiated to a callable term, find a frame executing the predicate described by `Value` and unify the arguments of `Value` to the goal arguments associated with the frame. This is intended to check the current execution context. The user must ensure the checked parent goal is not removed from the stack due to last-call optimisation and be aware of the slow operation on deeply nested calls.

The variant `parent_goal(-Parent)` unifies the frame reference of the parent of the found frame with `Parent`. That allows for finding frames higher up in the stack running the same goal.

**predicate_indicator**  
Similar to `goal`, but only returning the \[\<`module`\>:\]\<`name`\>/\<`arity`\> term describing the term, not the actual arguments. It avoids creating an illegal term as `goal` and is used by the library `library(prolog_stack)`.

**clause**  
`Value` is unified with a reference to the currently running clause. Fails if the current goal is associated with a foreign (C) defined predicate. See also [nth_clause/3](examineprog.html#nth_clause/3) and [clause_property/2](examineprog.html#clause_property/2).

**level**  
`Value` is unified with the recursion level of `Frame`. The top level frame is at level‘0’.

**parent**  
`Value` is unified with an integer reference to the parent local stack frame of `Frame`. Fails if `Frame` is the top frame.

**context_module**  
`Value` is unified with the name of the context module of the environment.

**top**  
`Value` is unified with `true` if `Frame` is the top Prolog goal from a recursive call back from the foreign language; `false` otherwise.

**hidden**  
`Value` is unified with `true` if the frame is hidden from the user, either because a parent has the hide-childs attribute (all system predicates), or the system has no trace-me attribute.

**skipped**  
`Value` is `true` if this frame was skipped in the debugger.

**pc**  
`Value` is unified with the program pointer saved on behalf of the parent goal if the parent goal is not owned by a foreign predicate or belongs to a compound meta-call (e.g., call((a,b))).

**argument**(`N`)  
`Value` is unified with the `N`-th slot of the frame. Argument 1 is the first argument of the goal. Arguments above the arity refer to local variables. Fails silently if `N` is out of range.

**prolog_choice_attribute**(`+ChoicePoint, +Key, -Value`)  
Extract attributes of a choice point. `ChoicePoint` is a reference to a choice point as passed to [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4) on the 3rd argument or obtained using [prolog_current_choice/1](manipstack.html#prolog_current_choice/1). `Key` specifies the requested information:

**parent**  
Requests a reference to the first older choice point.

**frame**  
Requests a reference to the frame to which the choice point refers.

**type**  
Requests the type. Defined values are `clause` (the goal has alternative clauses), `foreign` (non-deterministic foreign predicate), `jump` (clause internal choice point), `top` (first dummy choice point), `catch` ([catch/3](exception.html#catch/3) to allow for undo), `debug` (help the debugger), or `none` (has been deleted).

**pc**  
Requests the program counter to which the choice point refers. Only applicable for in-clause choice points.

**clause**  
Request the clause that will be tried if this choice point is activated. Only applicable for choice points of type `clause`.

This predicate is used for the graphical debugger to show the choice point stack.

**deterministic**(`-Boolean`)  
Unifies its argument with `true` if no choice point exists that is more recent than the entry of the clause in which it appears. There are few realistic situations for using this predicate. It is used by the [prolog/0](toplevel.html#prolog/0) top level to check whether Prolog should prompt the user for alternatives. Similar results can be achieved in a more portable fashion using [call_cleanup/2](metacall.html#call_cleanup/2).
