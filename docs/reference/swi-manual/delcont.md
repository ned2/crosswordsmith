
## 4.9 Delimited continuations

The predicates [reset/3](delcont.html#reset/3) and [shift/1](delcont.html#shift/1) implement *delimited continuations* for Prolog. Delimited continuations for Prolog are described in [Schrijvers *et al.*, 2013](Bibliography.html#DBLP:journals/tplp/SchrijversDDW13) ([preprint PDF](https://www.swi-prolog.org/download/publications/iclp2013.pdf)). The mechanism allows for proper *coroutines*, two or more routines whose execution is interleaved, while they exchange data. Note that coroutines in this sense differ from coroutines realised using attributed variables as described in [chapter 8](clp.html#sec:8).

Note that [shift/1](delcont.html#shift/1) captures the *forward continuation*. It notably does not capture choicepoints. Choicepoints created before the continuation is captured remain open, while choicepoints created when the continuation is executed live their normal life. Unfortunately the consequences for *committing* a choicepoint is complicated. In general a cut ([!/0](control.html#!/0)) in the continuation does not have the expected result. Negation ([\\/1](control.html#\+/1)) and if-then(-else) ([-\>/2](control.html#-%3E/2)) behave as expected, *provided the continuation is called immediately*. This works because for [\\/1](control.html#\+/1) and [-\>/2](control.html#-%3E/2) the continuation contains a reference to the choicepoint that must be cancelled and this reference is restored when possible. If, as with tabling, the continuation is saved and called later, the commit has no effect. We illustrate the three scenarios using with the programs below.

``` code
t1 :-
    reset(gbad, ball, Cont),
    (   Cont == 0
    ->  true
    ;   writeln(resuming),
        call(Cont)
    ).

gbad :-
    n, !, fail.
gbad.

n :-
    shift(ball),
    writeln(n).
```

Here, the [!/0](control.html#!/0) has **no effect**:

``` code
?- t1.
resuming
n
true.
```

The second example uses [\\/1](control.html#\+/1), which is essentially `(G->fail;true)`.

``` code
t2 :-
    reset(gok, ball, Cont),
    (   Cont == 0
    ->  true
    ;   writeln(resuming),
        call(Cont)
    ).

gok :-
    \+ n.
```

In this scenario the normal semantics of [\\/1](control.html#\+/1) is preserved:

``` code
?- t1.
resuming
n
false.
```

In the last example we illustrate what happens if we assert the continuation to be executed later. We write the negation using if-then-else to make it easier to explain the behaviour.

``` code
:- dynamic cont/1.

t3 :-
    retractall(cont(_)),
    reset(gassert, ball, Cont),
    (   Cont == 0
    ->  true
    ;   asserta(cont(Cont))
    ).

c3 :-
    cont(Cont),
    writeln(resuming),
    call(Cont).

gassert :-
    (   n
    ->  fail
    ;   true
    ).
```

Now, t3/0 succeeds *twice*. This is because n/0 shifts, so the commit to the [fail/0](control.html#fail/0) branch is not executed and the [true/0](control.html#true/0) branch is evaluated normally. Calling the continuation later using c3/0 fails because the choicepoint that realised the if-then-else does not exist in the continuation and thus the effective continuation is the remainder of n/0 and [fail/0](control.html#fail/0) in gassert/0 .

``` code
?- t3.
true ;
true.

?- c3.
resuming
n
false.
```

The suspension mechanism provided by delimited continuations is used to implement *tabling* [Desouter *et al.*, 2015](Bibliography.html#DBLP:journals/tplp/DesouterDS15), ([available here](https://www.cambridge.org/core/journals/theory-and-practice-of-logic-programming/article/div-classtitletabling-as-a-library-with-delimited-controldiv/227B7C0227FD715CF159B6AF894DE96E)). See [section 7](tabling.html#sec:7).

**reset**(`:Goal, ?Ball, -Continuation`)  
Call `Goal`. If `Goal` calls [shift/1](delcont.html#shift/1) and the argument of [shift/1](delcont.html#shift/1) can be unified with `Ball`,^(80The argument order described in [Schrijvers *et al.*, 2013](Bibliography.html#DBLP:journals/tplp/SchrijversDDW13) is `reset(Goal,Continuation,Ball)`. We swapped the argument order for compatibility with [catch/3](exception.html#catch/3)) [shift/1](delcont.html#shift/1) causes [reset/3](delcont.html#reset/3) to return, unifying `Continuation` with a goal that represents the *continuation* after [shift/1](delcont.html#shift/1). In other words, meta-calling `Continuation` completes the execution where shift left it. If `Goal` does not call [shift/1](delcont.html#shift/1), `Continuation` are unified with the integer `0` (zero).^(81Note that older versions also unify `Ball` with `0`. Testing whether or not shift happened on `Ball` however is *always* ambiguous.)

**shift**(`+Ball`)  
Abandon the execution of the current goal, returning control to just *after* the matching [reset/3](delcont.html#reset/3) call. This is similar to [throw/1](exception.html#throw/1) except that (1) nothing is‘undone’and (2) the 3th argument of [reset/3](delcont.html#reset/3) is unified with the *continuation*, which allows the code calling [reset/3](delcont.html#reset/3) to *resume* the current goal.

\[experimental\]**shift_for_copy**(`+Ball`)  
Similar to [shift/1](delcont.html#shift/1). This version is intended for situations where it is assumed the continuation is copied and saved to be executed one or multiple times in a different context. This notably prevents restoring choice points saved for [\\/1](control.html#\+/1), *If`->`Then`;`Else*, etc.
