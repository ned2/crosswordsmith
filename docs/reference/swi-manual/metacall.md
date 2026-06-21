
## 4.8 Meta-Call Predicates

Meta-call predicates are used to call terms constructed at run time. The basic meta-call mechanism offered by SWI-Prolog is to use variables as a subclause (which should of course be bound to a valid goal at runtime). A meta-call is slower than a normal call as it involves actually searching the database at runtime for the predicate, while for normal calls this search is done at compile time.

\[ISO\]**call**(`:Goal`)  
Call `Goal`. This predicate is normally used for goals that are not known at compile time. For example, the Prolog toplevel essentially performs `read(Goal), call(Goal)`. Also a *meta* predicates such as [ignore/1](metacall.html#ignore/1) are defined using call:

``` code
ignore(Goal) :- call(Goal), !.
ignore(_).
```

Note that a plain variable as a body term acts as [call/1](metacall.html#call/1) and the above is equivalent to the code below. SWI-Prolog produces the same code for these two programs and [listing/1](listing.html#listing/1) prints the program above.

``` code
ignore(Goal) :- Goal, !.
ignore(_).
```

Note that [call/1](metacall.html#call/1) restricts the scope of the cut ([!/0](control.html#!/0)). A cut inside `Goal` only affects choice points created by `Goal`.

\[ISO\]**call**(`:Goal, +ExtraArg1, ...`)  
Append `ExtraArg1, ExtraArg2, ...` to the argument list of `Goal` and call the result. For example, `call(plus(1), 2, X)` will call `plus(1, 2, X)`, binding `X` to 3.

The call/\[2..\] construct is handled by the compiler. The predicates call/\[2-8\] are defined as real (meta-)predicates and are available to inspection through [current_predicate/1](examineprog.html#current_predicate/1), [predicate_property/2](examineprog.html#predicate_property/2), etc.^(76Arities 2..8 are demanded by ISO/IEC 13211-1:1995/Cor.2:2012.) Higher arities are handled by the compiler and runtime system, but the predicates are not accessible for inspection.^(77Future versions of the reflective predicate may fake the presence of call/9.. . Full logical behaviour, generating all these pseudo predicates, is probably undesirable and will become impossible if *max_arity* is removed.)

**autoload_call**(`:Goal`)  
As [call/1](metacall.html#call/1), but first trigger the autoloader if `Goal` is not defined. The autoloader is used, even if the flag [autoload](flags.html#flag:autoload) is `false`. In addition, using this predicate causes the dependency checker, e.g., [list_undefined/0](check.html#list_undefined/0), to ignore the argument. This predicate is particularly used to call into the development tools.

\[deprecated\]**apply**(`:Goal, +List`)  
Append the members of `List` to the arguments of `Goal` and call the resulting term. For example: `apply(plus(1), [2, X])` calls `plus(1, 2, X)`. New code should use call/\[2..\] if the length of `List` is fixed.

\[deprecated\]**not**(`:Goal`)  
True if `Goal` cannot be proven. Retained for compatibility only. New code should use [\\/1](control.html#\+/1).

\[ISO\]**once**(`:Goal`)  
Make a possibly *nondet* goal *semidet*, i.e., succeed at most once. Defined as:

``` code
once(Goal) :-
    call(Goal), !.
```

[once/1](metacall.html#once/1) can in many cases be replaced with [-\>/2](control.html#-%3E/2). The only difference is how the cut behaves (see [!/0](control.html#!/0)). The following two clauses below are identical. Be careful about the interaction with [;/2](control.html#;/2). The `library(apply_macros)` library defines an inline expansion of [once/1](metacall.html#once/1), mapping it to `(Goal->true;fail)`. Using the full if-then-else constructs prevents its semantics from being changed when embedded in a [;/2](control.html#;/2) disjunction.

``` code
1) a :- once((b, c)), d.
2) a :- b, c -> d.
```

**ignore**(`:Goal`)  
Calls `Goal` as [once/1](metacall.html#once/1), but succeeds, regardless of whether `Goal` succeeded or not. Defined as:

``` code
ignore(Goal) :-
        Goal, !.
ignore(_).
```

**call_with_depth_limit**(`:Goal, +Limit, -Result`)  
If `Goal` can be proven without recursion deeper than `Limit` levels, [call_with_depth_limit/3](metacall.html#call_with_depth_limit/3) succeeds, binding `Result` to the deepest recursion level used during the proof. Otherwise, `Result` is unified with `depth_limit_exceeded` if the limit was exceeded during the proof, or the entire predicate fails if `Goal` fails without exceeding `Limit`.

The depth limit is guarded by the internal machinery. This may differ from the depth computed based on a theoretical model. For example, [true/0](control.html#true/0) is translated into an inline virtual machine instruction. Also, [repeat/0](control.html#repeat/0) is not implemented as below, but as a non-deterministic foreign predicate.

``` code
repeat.
repeat :-
        repeat.
```

As a result, [call_with_depth_limit/3](metacall.html#call_with_depth_limit/3) may still loop infinitely on programs that should theoretically finish in finite time. This problem can be cured by using Prolog equivalents to such built-in predicates.

This predicate may be used for theorem provers to realise techniques like *iterative deepening*. See also [call_with_inference_limit/3](metacall.html#call_with_inference_limit/3). It was implemented after discussion with Steve Moyle [smoyle@ermine.ox.ac.uk](mailto:smoyle@ermine.ox.ac.uk).

**call_with_inference_limit**(`:Goal, +Limit, -Result`)  
Equivalent to `call(Goal)`, but limits the number of inferences *for each solution of `Goal`*.^(78This predicate was realised after discussion with Ulrich Neumerkel and Markus Triska.). Execution may terminate as follows:

- If `Goal` does *not* terminate before the inference limit is exceeded, `Goal` is aborted by injecting the exception `inference_limit_exceeded` into its execution. After termination of `Goal`, `Result` is unified with the atom `inference_limit_exceeded`. *Otherwise*,
- If `Goal` fails, [call_with_inference_limit/3](metacall.html#call_with_inference_limit/3) fails.
- If `Goal` succeeds *without a choice point*, `Result` is unified with `!`.
- If `Goal` succeeds *with a choice point*, `Result` is unified with `true`.
- If `Goal` throws an exception, [call_with_inference_limit/3](metacall.html#call_with_inference_limit/3) re-throws the exception.

An inference is defined as a call or redo on a predicate. Please note that some primitive built-in predicates are compiled to virtual machine instructions for which inferences are not counted. The execution of predicates defined in other languages (e.g., C, C++) count as a single inference. This includes potentially expensive built-in predicates such as [sort/2](builtinlist.html#sort/2).

Calls to this predicate may be nested. An inner call that sets the limit below the current is honoured. An inner call that would terminate after the current limit does not change the effective limit. See also [call_with_depth_limit/3](metacall.html#call_with_depth_limit/3) and call_with_time_limit/2.

**setup_call_cleanup**(`:Setup, :Goal, :Cleanup`)  
Calls `(once(Setup), Goal)`. If `Setup` succeeds, `Cleanup` will be called exactly once after `Goal` is finished: either on failure, deterministic success, commit, or an exception. The execution of `Setup` is protected from asynchronous interrupts like call_with_time_limit/2 (package clib) or [thread_signal/2](threadcom.html#thread_signal/2). In most uses, `Setup` will perform temporary side-effects required by `Goal` that are finally undone by `Cleanup`.

Success or failure of `Cleanup` is ignored, and choice points it created are destroyed (as [once/1](metacall.html#once/1)). If `Cleanup` throws an exception, this is executed as normal while it was not triggered as the result of an exception the exception is propagated as normal. If `Cleanup` was triggered by an exception the rules are described in [section 4.10.2](exception.html#sec:4.10.2)

Typically, this predicate is used to cleanup permanent data storage required to execute `Goal`, close file descriptors, etc. The example below provides a non-deterministic search for a term in a file, closing the stream as needed.

``` code
term_in_file(Term, File) :-
        setup_call_cleanup(open(File, read, In),
                           term_in_stream(Term, In),
                           close(In) ).

term_in_stream(Term, In) :-
        repeat,
        read(In, T),
        (   T == end_of_file
        ->  !, fail
        ;   T = Term
        ).
```

Note that it is impossible to implement this predicate in Prolog. The closest approximation would be to read all terms into a list, close the file and call [member/2](lists.html#member/2). Without [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) there is no way to gain control if the choice point left by [repeat/0](control.html#repeat/0) is removed by a cut or an exception.

[setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) can also be used to test determinism of a goal, providing a portable alternative to [deterministic/1](manipstack.html#deterministic/1):

``` code
?- setup_call_cleanup(true,(X=1;X=2), Det=yes).

X = 1 ;

X = 2,
Det = yes ;
```

This predicate is under consideration for inclusion into the ISO standard. For compatibility with other Prolog implementations see [call_cleanup/2](metacall.html#call_cleanup/2).

**setup_call_catcher_cleanup**(`:Setup, :Goal, +Catcher, :Cleanup`)  
Similar to `setup_call_cleanup(Setup, Goal, Cleanup)` with additional information on the reason for calling `Cleanup`. Prior to calling `Cleanup`, `Catcher` unifies with the termination code (see below). If this unification fails, `Cleanup` is *not* called.

**exit**  
`Goal` succeeded without leaving any choice points.

**fail**  
`Goal` failed.

**`!`**  
`Goal` succeeded with choice points and these are now discarded by the execution of a cut (or other pruning of the search tree such as if-then-else).

**exception**(`Exception`)  
`Goal` raised the given `Exception`.

**external_exception**(`Exception`)  
`Goal` succeeded with choice points and these are now discarded due to an exception. For example:

``` code
?- setup_call_catcher_cleanup(true, (X=1;X=2),
                              Catcher, writeln(Catcher)),
   throw(ball).
external_exception(ball)
ERROR: Unhandled exception: Unknown message: ball
```

**call_cleanup**(`:Goal, :Cleanup`)  
Same as `setup_call_cleanup(true, Goal, Cleanup)`. This is provided for compatibility with a number of other Prolog implementations only. Do not use [call_cleanup/2](metacall.html#call_cleanup/2) if you perform side-effects prior to calling that will be undone by `Cleanup`. Instead, use [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) with an appropriate first argument to perform those side-effects.

**undo**(`:Goal`)  
Add `Goal` to the *trail*. `Goal` is executed as [ignore/1](metacall.html#ignore/1) on the first opportunity after backtracking to a point before the call to `Goal`. This predicate is intended to make otherwise persistent changes to the database or created by foreign procedures backtrackable if it is possible to define a goal that reverts the effect of the initial call. A typical use case is to define a *backtrackable assert*.

``` code
b_assertz(Term) :-
    assertz(Term, Ref),
    undo(erase(Ref)).
```

Without [undo/1](metacall.html#undo/1) we can achieve something similar by leaving a choicepoint using the almost portable^(79[assertz/2](db.html#assertz/2) is not part of the ISO standard but supported by multiple systems.) alternative below.

``` code
b_assertz(Term) :-
    assertz(Term, Ref),
    (   true
    ;   erase(Ref),
        fail
    ).
```

The [undo/1](metacall.html#undo/1) based solution avoids leaving a choice point open and, more importantly, keeps undoing the assert also if the choice point from the second alternative is pruned.

Currently the following remarks apply

- `Goal` is *copied* when it is registered.
- “First opportunity” means after backtracking or at the first call port reached.
- Multiple undo goals may be scheduled that are executed as a batch. If multiple goals raise an exception, the most urgent is preserved after all goals have been executed.
- It is not allowed for `Goal` to call [undo/1](metacall.html#undo/1). An attempt to do so results in a `permission_error` exception.
- Note that an exception that is caught higher in the call stack backtracks and therefore ensures `Goal` is called.

See also [snapshot/1](db.html#snapshot/1) and [transaction/1](db.html#transaction/1).
