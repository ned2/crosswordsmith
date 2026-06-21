
## 7.6 Well Founded Semantics

According to [Wikipedia](https://en.wikipedia.org/wiki/Well-founded_semantics), "*Well Founded Semantics* is one definition of how we can make conclusions from a set of logical rules". Well Founded Semantics (WFS) defines a *three valued logic* representing *true*, *false* and something that is neither true or false. This latter value is often referred to as *bottom*, *undefined* or *unknown*. SWI-Prolog uses [undefined/0](WFS.html#undefined/0).

Well Founded Semantics allows for reasoning about programs with contradictions or multiple answer sets. It allows for obtaining true/false results for literals that do not depend on the sub program that has no unambiguous solution, propagating the notion of *undefined* to literals that cannot be resolved otherwise and obtaining the *residual* program that expresses why an answer is not unambiguous.

The notion of *Well Founded Semantics* is only relevant if the program uses *negation* as implemented by [tnot/1](tabling-preds.html#tnot/1). The argument of [tnot/1](tabling-preds.html#tnot/1), as the name implies, must be a goal associated with a tabled predicate (see [table/1](tabling-preds.html#table/1)). In a nutshell, resolving a goal that implies [tnot/1](tabling-preds.html#tnot/1) is implemented as follows:

Consider the following partial *body term*:

``` code
        ...,
        tnot(p),
        q.
```

1.  If `p` has an unconditional answer in its table, fail.
2.  Else, *delay* the negation. If an unconditional answer arrives at some time, resume with failure.
3.  If at the end of the traditional tabled evaluation we can still not decide on `p`, execute the *continuation* (`q` above) while maintaining the *delay list* set to `tnot(p)`. If executing the continuation results in an answer for some tabled predicate, record this answer as a *conditional* answer, in this case with the condition `tnot(p)`.
4.  If a conditional answer is added to a table, it is propagated to its *followers*, say `f`, adding the pair {`f`,answer} to the delay list. If this leads to an answer, the answer is conditional on this pair.
5.  After the continuations of all unresolved [tnot/1](tabling-preds.html#tnot/1) calls have been executed the various tables may have conditional answers in addition to normal answers.
6.  If there are negative literals that have neither conditional answers nor unconditional answers, the condition `tnot(g)` is true. This conclusion is propagated by simplifying the conditions for all answers that depend on `tnot(g)`. This may result in a definite *false* condition, in which case the answer is removed or a definite *true* condition in which case the answer is made unconditional. Both events can make other conditional answers definitely true or false, etc.
7.  At the end of the simplifying process some answers may still be conditional. A final *answer completion* step analyses the graph of depending nodes, eliminating *positive loops*, e.g., “`p` :- `q`. `q` :- `p`” . The answers in such a loop are removed, possibly leading to more simplification. This process is executed until *fixed point* is reached, i.e., no further positive loops exist and no further simplification is possible.

The above process may complete without any remaining conditional answers, in which case we are back in the normal Prolog world. It is also possible that some answers remain conditional. The most obvious case is represented by [undefined/0](WFS.html#undefined/0). The toplevel responds with **undefined** instead of **true** if an answer is conditional.

**undefined**  
Unknown represents neither `true` nor `false` in the well formed model. It is implemented as

``` code
:- table undefined/0.

undefined :- tnot(undefined).
```

Solving a set of predicates under well formed semantics results in a *residual program*. This program contains clauses for all tabled predicates with condition answers where each clause head represents and answer and each clause body its condition. The condition is a disjunction of conjunctions where each literal is either a tabled goal or [tnot/1](tabling-preds.html#tnot/1) of a tabled goal. The remaining model has at least a cycle through a negative literal ([tnot/1](tabling-preds.html#tnot/1)) and has no single solution in the *stable model semantics*, i.e., it either expresses a contradiction (as [undefined/0](WFS.html#undefined/0), i.e., there is no stable model) or a multiple stable models as in the program below, where both {p} and {q} are stable models.

``` code
:- table p/0, q/0.

p :- tnot(q).
q :- tnot(p).
```

Note that it is possible that some literals have the same truth value in all stable models but are still *undefined* under the stable model semantics.

The residual program is an explanation of why an answer is undefined. SWI-Prolog offers the following predicates to access the residual program.

**call_residual_program**(`:Goal, -Program`)  
True when `Goal` is an answer according to the Well Founded Semantics. If `Program` is the empty list, `Goal` is unconditionally true. Otherwise this is a program as described by [delays_residual_program/2](WFS.html#delays_residual_program/2).

**call_delays**(`:Goal, -Condition`)  
True when `Goal` is an answer that is true when Condition can be satisfied. If `Condition` is `true`, `Answer` is unconditional. Otherwise it is a conjunction of goals, each of which is associated with a tabled predicate.

**delays_residual_program**(`:Condition, -Program`)  
Program is a list of clauses that represents the connected program associated with `Condition`. Each clause head represents a conditional answer from a table and each corresponding clause body is the condition that must hold for this answer to be true. The body is a disjunction of conjunctions. Each leaf in this condition is either a term `tnot(Goal)` or a plain `Goal`. Each `Goal` is associated with a tabled predicate. The program always contains at least one cycle that involves [tnot/1](tabling-preds.html#tnot/1).

### 7.6.1 Well founded semantics and the toplevel

The toplevel supports two modes for reporting that it is undefined whether the current answer is true. The mode is selected by the Prolog flag [toplevel_list_wfs_residual_program](flags.html#flag:toplevel_list_wfs_residual_program). If `true`, the toplevel uses [call_delays/2](WFS.html#call_delays/2) and [delays_residual_program/2](WFS.html#delays_residual_program/2) to find the conditional answers used and the *residual* program associated with these answers. It then prints the residual program, followed by the answer and the conditional answers. For [undefined/0](WFS.html#undefined/0), this results in the following output:

``` code
?- undefined.
% WFS residual program
    undefined :-
        tnot(undefined).
undefined.
```

If the [toplevel_list_wfs_residual_program](flags.html#flag:toplevel_list_wfs_residual_program) is false, any undefined answer is a conjunction with [undefined/0](WFS.html#undefined/0). See the program and output below.

``` code
:- table p/0, q/0.

p :- tnot(q).
q :- tnot(p).
```

``` code
?- p.
% WFS residual program
    p :-
        tnot(q).
    q :-
        tnot(p).
p.

?- set_prolog_flag(toplevel_list_wfs_residual_program, false).
true.

?- p.
undefined.
```
