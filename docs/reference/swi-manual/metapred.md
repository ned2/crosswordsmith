
## 6.5 Defining a meta-predicate

A meta-predicate is a predicate that calls other predicates dynamically, modifies a predicate, or reasons about properties of a predicate. Such predicates use either a compound term or a *predicate indicator* to describe the predicate they address, e.g., `assert(name(jan))` or `abolish(``name/1)`. With modules, this simple schema no longer works as each module defines its own mapping from name+arity to predicate. This is resolved by wrapping the original description in a term \<`module`\>:\<`term`\>, e.g., `assert(person:name(jan))` or `abolish(``person:name/1)`.

Of course, when calling [assert/1](db.html#assert/1) from inside a module, we expect to assert to a predicate local to this module. In other words, we do not wish to provide this `:``/2` wrapper by hand. The [meta_predicate/1](metapred.html#meta_predicate/1) directive tells the compiler that certain arguments are terms that will be used to look up a predicate and thus need to be wrapped (qualified) with \<`module`\>:\<`term`\>, unless they are already wrapped.

In the example below, we use this to define [maplist/3](apply.html#maplist/3) inside a module. The argument‘2’in the meta_predicate declaration means that the argument is module-sensitive and refers to a predicate with an arity that is two more than the term that is passed in. The compiler only distinguishes the values 0..9 and `:`, which denote module-sensitive arguments, from `+`, `-` and `?`, which denote *modes*. The values 0..9 are used by the *cross-referencer* and syntax highlighting. Note that the helper predicate maplist\_/3 does not need to be declared as a meta-predicate because the [maplist/3](apply.html#maplist/3) wrapper already ensures that `Goal` is qualified as \<`module`\>:`Goal`. See the description of [meta_predicate/1](metapred.html#meta_predicate/1) for details.

``` code
:- module(maplist, [maplist/3]).
:- meta_predicate maplist(2, ?, ?).

%%      maplist(:Goal, +List1, ?List2)
%
%       True if Goal can successfully be applied to all
%       successive pairs of elements from List1 and List2.

maplist(Goal, L1, L2) :-
        maplist_(L1, L2, Goal).

maplist_([], [], _).
maplist_([H0|T0], [H|T], Goal) :-
        call(Goal, H0, H),
        maplist_(T0, T, Goal).
```

**meta_predicate** `+Head, ...`  
Define the predicates referenced by the comma-separated list `Head` as *meta-predicates*. Each argument of each head is a *meta argument specifier*. Defined specifiers are given below. Only 0..9, `:`, `^` and `//` are interpreted; the mode declarations `+`, `-`, `*` and `?` are ignored.

**0..9**  
The argument is a term that is used to reference a predicate with `N` more arguments than the given argument term. For example: `call(0)` or `maplist(1, +)`.

**`:`**  
The argument is module-sensitive, but does not directly refer to a predicate. For example: `consult(:)`.

**`^`**  
This extension is used to denote the possibly `^`-annotated goal of [setof/3](allsolutions.html#setof/3), [bagof/3](allsolutions.html#bagof/3), [aggregate/3](aggregate.html#aggregate/3) and [aggregate/4](aggregate.html#aggregate/4). It is processed similar to‘0’, but leaving the `^`/2 intact.

**`//`**  
The argument is a DCG body. See [phrase/3](DCG.html#phrase/3).

**`-`**  
**`?`**  
**`*`**  
**`+`**  
All these have the same semantics, declaring the argument to be not module sensitive. The `*` notation is an alias for `?` for compatibility with e.g., Logtalk. The specific mode has merely documentation value. See [section 4.1.1](preddesc.html#sec:4.1.1) for details.

Each argument that is module-sensitive (i.e., marked 0..9, `:` or `^`) is qualified with the context module of the caller if it is not already qualified. The implementation ensures that the argument is passed as \<`module`\>:\<`term`\>, where \<`module`\> is an atom denoting the name of a module and \<`term`\> itself is not a `:``/2` term where the first argument is an atom. Below is a simple declaration and a number of queries.

``` code
:- meta_predicate
        meta(0, +).

meta(Module:Term, _Arg) :-
        format('Module=~w, Term = ~q~n', [Module, Term]).
```

``` code
?- meta(test, x).
Module=user, Term = test
?- meta(m1:test, x).
Module=m1, Term = test
?- m2:meta(test, x).
Module=m2, Term = test
?- m1:meta(m2:test, x).
Module=m2, Term = test
?- meta(m1:m2:test, x).
Module=m2, Term = test
?- meta(m1:42:test, x).
Module=42, Term = test
```

The [meta_predicate/1](metapred.html#meta_predicate/1) declaration is the portable mechanism for defining meta-predicates and replaces the old SWI-Prolog specific mechanism provided by the deprecated predicates [module_transparent/1](ctxmodule.html#module_transparent/1), [context_module/1](ctxmodule.html#context_module/1) and [strip_module/3](ctxmodule.html#strip_module/3). See also [section 6.16](modulecompat.html#sec:6.16).
