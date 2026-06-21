
## 5.3 Syntax changes since SWI-Prolog 7

### 5.3.1 Operators and quoted atoms

As of SWI-Prolog version 7, quoted atoms lose their operator property. This means that expressions such as `A = 'dynamic'/1` are valid syntax, regardless of the operator definitions. From questions on the mailinglist this is what people expect.^(174We believe that most users expect an operator declaration to define a new token, which would explain why the operator name is often quoted in the declaration, but not while the operator is used. We are afraid that allowing for this easily creates ambiguous syntax. Also, many development environments are based on tokenization. Having dynamic tokenization due to operator declarations would make it hard to support Prolog in such editors.) To accommodate for real quoted operators, a quoted atom that *needs* quotes can still act as an operator.^(175Suggested by Joachim Schimpf.) A good use-case for this is a unit library^(176[https://groups.google.com/d/msg/comp.lang.prolog/ozqdzI-gi_g/2G16GYLIS0IJ](https://groups.google.com/d/msg/comp.lang.prolog/ozqdzI-gi_g/2G16GYLIS0IJ)), which allows for expressions such as below.

``` code
?- Y isu 600kcal - 1h*200'W'.
Y = 1790400.0'J'.
```

### 5.3.2 Compound terms with zero arguments

As of SWI-Prolog version 7, the system supports compound terms that have no arguments. This implies that e.g., **`name()`** is valid syntax. This extension aims at functions on dicts (see [section 5.4](bidicts.html#sec:5.4)) as well as the implementation of domain specific languages (DSLs). To minimise the consequences, the classic predicates [functor/3](manipterm.html#functor/3) and [=../2](manipterm.html#=../2) have not been modified. The predicates [compound_name_arity/3](manipterm.html#compound_name_arity/3) and [compound_name_arguments/3](manipterm.html#compound_name_arguments/3) have been added. These predicates operate only on compound terms and behave consistently for compounds with zero arguments. Code that *generalises* a term using the sequence below should generally be changed to use [compound_name_arity/3](manipterm.html#compound_name_arity/3).

``` code
    ...,
    functor(Specific, Name, Arity),
    functor(General, Name, Arity),
    ...,
```

Replacement of [=../2](manipterm.html#=../2) by [compound_name_arguments/3](manipterm.html#compound_name_arguments/3) is typically needed to deal with code that follow the skeleton below.

``` code
    ...,
    Term0 =.. [Name|Args0],
    maplist(convert, Args0, Args),
    Term =.. [Name|Args],
    ...,
```

For predicates, goals and arithmetic functions (evaluable terms), \<`name`\> and \<`name`\>() are *equivalent*. Below are some examples that illustrate this behaviour.

``` code
go() :- format('Hello world~n').

?- go().
Hello world

?- go.
Hello world

?- Pi is pi().
Pi = 3.141592653589793.

?- Pi is pi.
Pi = 3.141592653589793.
```

Note that the *canonical* representation of predicate heads and functions without arguments is an atom. Thus, `clause(`**`go()`**`, Body)` returns the clauses for go/0 , but `clause(-Head, -Body, +Ref)` unifies `Head` with an atom if the clause specified by `Ref` is part of a predicate with zero arguments.

### 5.3.3 Block operators

Introducing curly bracket and array subscripting.^(177Introducing block operators was proposed by Jose Morales. It was discussed in the Prolog standardization mailing list, but there were too many conflicts with existing extensions (ECLiPSe and B-Prolog) and doubt about their need to reach an agreement. Increasing need to get to some solution resulted in what is documented in this section. These extensions are also implemented in recent versions of YAP.) The symbols `[]` and `{}` may be declared as an operator, which has the following effect:

**\[ \]**  
This operator is typically declared as a low-priority `yf` postfix operator, which allows for `array[index]` notation. This syntax produces a term `[]([index],array)`.

**{ }**  
This operator is typically declared as a low-priority `xf` postfix operator, which allows for `head(arg) { body }` notation. This syntax produces a term `{}({body},head(arg))`.

Below is an example that illustrates the representation of a typical‘curly bracket language’in Prolog.

``` code
?- op(100, xf, {}).
?- op(100, yf, []).
?- op(1100, yf, ;).

?- displayq(func(arg)
            { a[10] = 5;
              update();
            }).
{}({;(=([]([10],a),5),;(update()))},func(arg))
```
