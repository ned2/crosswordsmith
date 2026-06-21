
## 4.7 Control Predicates

The predicates of this section implement control structures. Normally the constructs in this section, except for [repeat/0](control.html#repeat/0), are translated by the compiler. Please note that complex goals passed as arguments to meta-predicates such as [findall/3](allsolutions.html#findall/3) below cause the goal to be compiled to a temporary location before execution. It is faster to define a sub-predicate (i.e., one_character_atoms/1 in the example below) and make a call to this simple predicate. See also the Prolog flag [compile_meta_arguments](flags.html#flag:compile_meta_arguments).

``` code
one_character_atoms(As) :-
        findall(A, (current_atom(A), atom_length(A, 1)), As).
```

\[ISO\]**fail**  
Always fail. The predicate [fail/0](control.html#fail/0) is translated into a single virtual machine instruction.

\[ISO\]**false**  
Same as fail, but the name has a more declarative connotation.

\[ISO\]**true**  
Always succeed. The predicate [true/0](control.html#true/0) is translated into a single virtual machine instruction.

\[ISO\]**repeat**  
Always succeed, provide an infinite number of choice points.

\[ISO\]**!**  
Cut. Discard all choice points created since entering the predicate in which the cut appears. In other words, *commit* to the clause in which the cut appears *and* discard choice points that have been created by goals to the left of the cut in the current clause. Meta calling is opaque to the cut. This implies that cuts that appear in a term that is subject to meta-calling ([call/1](metacall.html#call/1)) only affect choice points created by the meta-called term. The following control structures are transparent to the cut: [;/2](control.html#;/2), [-\>/2](control.html#-%3E/2) and [\*-\>/2](control.html#*-%3E/2). Cuts appearing in the *condition* part of [-\>/2](control.html#-%3E/2) and [\*-\>/2](control.html#*-%3E/2) are opaque to the cut. The table below explains the scope of the cut with examples. *Prunes* here means “prunes `X` choice point created by `X`” .

|                                 |                       |
|---------------------------------|-----------------------|
| `t0 :- (a, !, b).`              | % prunes a/0 and t0/0 |
| `t1 :- (a, !, fail ; b).`       | % prunes a/0 and t1/0 |
| `t2 :- (a -> b, ! ; c).`        | % prunes b/0 and t2/0 |
| `t3 :- (a, !, b -> c ; d).`     | % prunes a/0          |
| `t4 :- call((a, !, fail ; b)).` | % prunes a/0          |
| `t5 :- ``\+``(a, !, fail).`     | % prunes a/0          |

\[ISO\]`:Goal1` **,** `:Goal2`  
Conjunction (*and*). True if both `Goal1` and `Goal2` are true.

\[ISO\]`:Goal1` **;** `:Goal2`  
Disjunction (*or*). True if either `Goal1` or `Goal2` succeeds. Note that the semantics change if `Goal1` contains [-\>/2](control.html#-%3E/2) or [\*-\>/2](control.html#*-%3E/2). [;/2](control.html#;/2) is transparent to cuts. See [!/0](control.html#!/0) for details. For example:

``` code
?- (between(1,2,X) ; X = a).
X = 1 ;
X = 2 ;
X = a.
```

It is strongly advised to *always* use parenthesis around disjunctions. Conjunctions inside a disjunction should not use parenthesis. Traditionally the `;` is placed at the start of the line rather than at the end because a `;` at the end of a line is easily overlooked. Below is an example of the preferred style used in SWI-Prolog.^(74Some users prefer a newline after the `;`.)

``` code
p :-
    a,
    (   b,
        c
    ;   d
    ).
```

Although [;/2](control.html#;/2) is a *control structure* that is normally handled by the compiler, SWI-Prolog implements [;/2](control.html#;/2) as a true predicate to support [call/2](metacall.html#call/2) and friends as well as to allow for querying predicate properties, for example to support code analysis.

`:Goal1` **\|** `:Goal2`  
Equivalent to [;/2](control.html#;/2). Retained for compatibility only. New code should use [;/2](control.html#;/2).

\[ISO\]`:Condition` **-\>** `:Action`  
If-then and If-Then-Else. The [-\>/2](control.html#-%3E/2) construct commits to the choices made at its left-hand side, destroying choice points created inside the clause (by [;/2](control.html#;/2)), or by goals called by this clause. Unlike [!/0](control.html#!/0), the choice point of the predicate as a whole (due to multiple clauses) is **not** destroyed. Disregarding the interaction with [!/0](control.html#!/0), the combination [;/2](control.html#;/2) and [-\>/2](control.html#-%3E/2) acts as if defined as:

``` code
If -> Then; _Else :- If, !, Then.
If -> _Then; Else :- !, Else.
If -> Then :- If, !, Then.
```

Please note that (If `->` Then) acts as (If `->` Then ; **fail**), making the construct *fail* if the condition fails. This unusual semantics is part of the ISO and all de-facto Prolog standards.

Please note that (if`->`then;else) is read as ((if`->`then);else) and that the *combined* semantics of this syntactic construct as defined above is *different* from the simple nesting of the two individual constructs, i.e., the semantics of [-\>/2](control.html#-%3E/2) *changes* when embedded in [;/2](control.html#;/2). See also [once/1](metacall.html#once/1).

As with [;/2](control.html#;/2), this construct is always nested in parenthesis. Here is an example of the preferred layout for SWI-Prolog.

``` code
p :-
    a,
    (   b,
        c
    ->  d,
        e
    ;   f
    ->  g
    ;   h
    ).
```

`:Condition` **\*-\>** `:Action ; :Else`  
This construct implements the so-called‘soft-cut’. The control is defined as follows: If `Condition` succeeds at least once, the semantics is the same as (`call(Condition)`, `Action`).^(75Note that the `Condition` is wrapped in [call/1](metacall.html#call/1), limiting the scope of the cut ([!/0](control.html#!/0)) If `Condition` does not succeed, the semantics is that of (`\+` `Condition`, `Else`). In other words, if `Condition` succeeds at least once, simply behave as the conjunction of `call(Condition)` and `Action`, otherwise execute `Else`. The construct is known under the name if/3 in some other Prolog implementations.

The construct `A` `*->` `B`, i.e., without an `Else` branch, the semantics is the same as (`call(A)`, `B`).

This construct is rarely used. An example use case is the implementation of optional in sparql. The optional construct should preserve all solutions if the argument succeeds as least once but still succeed otherwise. This is implemented as below.

``` code
optional(Goal) :-
    (   Goal
    *-> true
    ;   true
    ).
```

Now calling e.g., `optional(member(X, [a,b]))` has the solutions `X=a` and `X=b`, while `optional(member(X,[]))` succeeds without binding `X`.

\[ISO\]**\\** `:Goal`  
True if‘Goal’cannot be proven (mnemonic: `+` refers to *provable* and the backslash (`\`) is normally used to indicate negation in Prolog). In contrast to the ISO standard, but compatible with several other Prolog systems, SWI-Prolog implements [\\/1](control.html#\+/1) as a *control structure*. This implies that its argument is compiled as part of the enclosing clause and possible variables in goal positions are translated to [call/1](metacall.html#call/1). As a result, if such a variable is at runtime bound to a ([!/0](control.html#!/0)), the cut is scoped to the [call/1](metacall.html#call/1) call rather than the enclosing [\\/1](control.html#\+/1).

Many Prolog implementations (including SWI-Prolog) provide [not/1](metacall.html#not/1). The [not/1](metacall.html#not/1) alternative is deprecated due to its strong link to logical negation.
