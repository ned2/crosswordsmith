
## A.66 library(yall): Lambda expressions

author  
Paulo Moura and Jan Wielemaker

To be done  
Extend optimization support

Prolog realizes *high-order* programming with meta-calling. The core predicate of this is [call/1](metacall.html#call/1), which simply calls its argument. This can be used to define higher-order predicates such as [ignore/1](metacall.html#ignore/1) or [forall/2](forall2.html#forall/2). The call/N construct calls a *closure* with N-1 *additional arguments*. This is used to define higher-order predicates such as the [maplist/2](apply.html#maplist/2)-5 family or [foldl/4](apply.html#foldl/4)-7.

The *closure* concept used here is somewhat different from the closure concept from functional programming. The latter is a function that is always evaluated in the context that existed at function creation time. Here, a closure is a term of arity *0 `=<` L `=<` K*. The term's functor is the name of a predicate of arity *K* and the term's *L* arguments (where *L* could be 0) correspond to *L* leftmost arguments of said predicate, bound to parameter values. For example, a closure involving [atom_concat/3](manipatom.html#atom_concat/3) might be the term `atom_concat(prefix)`. In order of increasing *L*, one would have increasingly more complete closures that could be passed to call/3, all giving the same result:

``` code
call(atom_concat,prefix,suffix,R).
call(atom_concat(prefix),suffix,R).
call(atom_concat(prefix,suffix),R).
call(atom_concat(prefix,suffix,R)).
```

The problem with higher order predicates based on call/N is that the additional arguments are always added to the end of the closure's argument list. This often requires defining trivial helper predicates to get the argument order right. For example, if you want to add a common postfix to a list of atoms you need to apply `atom_concat(In,Postfix,Out)`, but `maplist(atom_concat(Postfix),ListIn,ListOut)` calls `atom_concat(Postfix,In,Out)`. This is where `library(yall)` comes in, where the module name, *yall*, stands for *Yet Another Lambda Library*.

The library allows us to write a lambda expression that *wraps around* the (possibly complex) goal to call:

``` code
?- maplist([In,Out]>>atom_concat(In,'_p',Out), [a,b], ListOut).
ListOut = [a_p, b_p].
```

A bracy list `{...}` specifies which variables are *shared* between the wrapped goal and the surrounding context. This allows us to write the code below. Without the `{Postfix}` a fresh variable would be passed to [atom_concat/3](manipatom.html#atom_concat/3).

``` code
add_postfix(Postfix, ListIn, ListOut) :-
    maplist({Postfix}/[In,Out]>>atom_concat(In,Postfix,Out),
            ListIn, ListOut).
```

This introduces the second application area of lambda expressions: the ability to confine variables to the called goal's context. This features shines when combined with [bagof/3](allsolutions.html#bagof/3) or [setof/3](allsolutions.html#setof/3) where one normally has to list those variables whose bindings one is *not* interested in using the `Var^Goal` construct (marking `Var` as existentially quantified and confining it to the called goal's context). Lambda expressions allow you to do the converse: specify the variables which one *is* interested in. These variables are common to the context of the called goal and the surrounding context.

Lambda expressions use the syntax below

``` code
{...}/[...]>>Goal.
```

The `{...}` optional part is used for lambda-free variables (the ones shared between contexts). The order of variables doesn't matter, hence the `{...}` set notation.

The `[...]` optional part lists lambda parameters. Here, order of variables matters, hence the list notation.

As `/` and `>>` are standard infix operators, no new operators are added by this library. An advantage of this syntax is that we can simply unify a lambda expression with `{Free}/[Parameters]>>Lambda` to access each of its components. Spaces in the lambda expression are not a problem although the goal may need to be written between’()'s. Goals that are qualified by a module prefix also need to be wrapped inside parentheses.

Combined with `library(apply_macros)`, `library(yall)` allows writing one-liners for many list operations that have the same performance as hand-written code.

This module implements [Logtalk's lambda expressions syntax](https://logtalk.org/manuals/refman/grammar.html\#lambda-expressions).

The development of this module was sponsored by Kyndi, Inc.

`+Parameters` **\>\>** `+Lambda`  
**\>\>**(`+Parameters, +Lambda, ?A1`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2, ?A3`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2, ?A3, ?A4`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2, ?A3, ?A4, ?A5`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2, ?A3, ?A4, ?A5, ?A6`)  
**\>\>**(`+Parameters, +Lambda, ?A1, ?A2, ?A3, ?A4, ?A5, ?A6, ?A7`)  
Calls a copy of `Lambda`. This is similar to `call(Lambda,A1,...)`, but arguments are reordered according to the list `Parameters`:

- The first `length(Parameters)` arguments from `A1`, ... are unified with (a copy of) `Parameters`, which *may* share them with variables in `Lambda`.
- Possible excess arguments are passed by position.

|  |  |
|----|----|
| `Parameters` | is either a plain list of parameters or a term `{Free}/List`. `Free` represents variables that are shared between the context and the `Lambda` term. This is needed for compiling `Lambda` expressions. |

`+Free` **/** `:Lambda`  
**/**(`+Free, :Lambda, ?A1`)  
**/**(`+Free, :Lambda, ?A1, ?A2`)  
**/**(`+Free, :Lambda, ?A1, ?A2, ?A3`)  
**/**(`+Free, :Lambda, ?A1, ?A2, ?A3, ?A4`)  
**/**(`+Free, :Lambda, ?A1, ?A2, ?A3, ?A4, ?A5`)  
**/**(`+Free, :Lambda, ?A1, ?A2, ?A3, ?A4, ?A5, ?A6`)  
**/**(`+Free, :Lambda, ?A1, ?A2, ?A3, ?A4, ?A5, ?A6, ?A7`)  
Shorthand for `Free/[]>>Lambda`. This is the same as applying call/N on `Lambda`, except that only variables appearing in `Free` are bound by the call. For example

``` code
p(1,a).
p(2,b).

?- {X}/p(X,Y).
X = 1;
X = 2.
```

This can in particularly be combined with [bagof/3](allsolutions.html#bagof/3) and [setof/3](allsolutions.html#setof/3) to *select* particular variables to be concerned rather than using existential quantification (^/2) to *exclude* variables. For example, the two calls below are equivalent.

``` code
setof(X, Y^p(X,Y), Xs)
setof(X, {X}/p(X,_), Xs)
```

\[semidet\]**is_lambda**(`@Term`)  
True if `Term` is a valid Lambda expression.

\[det\]**lambda_calls**(`+LambdaExpression, -Goal`)  
\[det\]**lambda_calls**(`+LambdaExpression, +ExtraArgs, -Goal`)  
`Goal` is the goal called if call/N is applied to `LambdaExpression`, where `ExtraArgs` are the additional arguments to call/N. `ExtraArgs` can be an integer or a list of concrete arguments. This predicate is used for cross-referencing and code highlighting.
