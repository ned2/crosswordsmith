
## 4.5 Verify Type of a Term

Type tests are semi-deterministic predicates that succeed if the argument satisfies the requested type. Type-test predicates have no error condition and do not instantiate their argument. See also library `library(error)`.

\[ISO\]**var**(`@Term`)  
True if `Term` currently is a free variable.

\[ISO\]**nonvar**(`@Term`)  
True if `Term` currently is not a free variable.

\[ISO\]**integer**(`@Term`)  
True if `Term` is bound to an integer.

\[ISO\]**float**(`@Term`)  
True if `Term` is bound to a floating point number.

**rational**(`@Term`)  
True if `Term` is bound to a rational number. Rational numbers include integers.

**rational**(`@Term, -Numerator, -Denominator`)  
True if `Term` is a rational number with given `Numerator` and `Denominator`. The `Numerator` and `Denominator` are in canonical form, which means `Denominator` is a positive integer and there are no common divisors between `Numerator` and `Denominator`.

\[ISO\]**number**(`@Term`)  
True if `Term` is bound to a rational number (including integers) or a floating point number.

\[ISO\]**atom**(`@Term`)  
True if `Term` is bound to an atom.

**blob**(`@Term, ?Type`)  
True if `Term` is a *blob* of type `Type`. See [section 12.4.10](foreigninclude.html#sec:12.4.10).

**string**(`@Term`)  
True if `Term` is bound to a string. Note that string here refers to the built-in atomic type string as described in [section 5.2](string.html#sec:5.2). Starting with version 7, the syntax for a string object is text between double quotes, such as `"hello"`.^(65In traditional Prolog systems, double quoted text is often mapped to a list of *character codes*.) See also the Prolog flag [double_quotes](flags.html#flag:double_quotes).

\[ISO\]**atomic**(`@Term`)  
True if `Term` is bound (i.e., not a variable) and is not compound. Thus, atomic acts as if defined by:

``` code
atomic(Term) :-
        nonvar(Term),
        \+ compound(Term).
```

SWI-Prolog defines the following atomic datatypes: atom ([atom/1](typetest.html#atom/1)), string ([string/1](typetest.html#string/1)), integer ([integer/1](typetest.html#integer/1)), floating point number ([float/1](typetest.html#float/1)), rational ([rational/1](typetest.html#rational/1)) and blob ([blob/2](typetest.html#blob/2)). In addition, the symbol `[]` (empty list) is atomic, but not an atom. See [section 5.1](ext-lists.html#sec:5.1).

\[ISO\]**compound**(`@Term`)  
True if `Term` is bound to a compound term. See also [functor/3](manipterm.html#functor/3) =../2, [compound_name_arity/3](manipterm.html#compound_name_arity/3) and [compound_name_arguments/3](manipterm.html#compound_name_arguments/3).

\[ISO\]**callable**(`@Term`)  
True if `Term` is bound to an atom or a compound term. This was intended as a type-test for arguments to [call/1](metacall.html#call/1), [call/2](metacall.html#call/2) etc. Note that callable only tests the *surface term*. Terms such as (22,true) are considered callable, but cause [call/1](metacall.html#call/1) to raise a type error. Module-qualification of meta-argument (see [meta_predicate/1](metapred.html#meta_predicate/1)) using `:``/2` causes callable to succeed on any meta-argument.^(66We think that [callable/1](typetest.html#callable/1) should be deprecated and there should be two new predicates, one performing a test for callable that is minimally module aware and possibly consistent with type-checking in [call/1](metacall.html#call/1) and a second predicate that tests for atom or compound.) Consider the program and query below:

``` code
:- meta_predicate p(0).

p(G) :- callable(G), call(G).

?- p(22).
ERROR: Type error: `callable' expected, found `22'
ERROR: In:
ERROR:    [6] p(user:22)
```

\[ISO\]**ground**(`@Term`)  
True if `Term` holds no free variables. See also [nonground/2](manipterm.html#nonground/2) and [term_variables/2](manipterm.html#term_variables/2).

**cyclic_term**(`@Term`)  
True if `Term` contains cycles, i.e. is an infinite term. See also [acyclic_term/1](typetest.html#acyclic_term/1) and [section 2.16](cyclic.html#sec:2.16).^(67The predicates [cyclic_term/1](typetest.html#cyclic_term/1) and [acyclic_term/1](typetest.html#acyclic_term/1) are compatible with SICStus Prolog. Some Prolog systems supporting cyclic terms use is_cyclic/1 .)

\[ISO\]**acyclic_term**(`@Term`)  
True if `Term` does not contain cycles, i.e. can be processed recursively in finite time. See also [cyclic_term/1](typetest.html#cyclic_term/1) and [section 2.16](cyclic.html#sec:2.16).
