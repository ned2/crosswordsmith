
## 9.3 CHR in SWI-Prolog Programs

### 9.3.1 Embedding CHR in Prolog Programs

The CHR constraints defined in a `.pl` file are associated with a module. The default module is `user`. One should never load different `.pl` files with the same CHR module name.

### 9.3.2 CHR Constraint declaration

:- **chr_constraint**(`+Specifier`)  
Every constraint used in CHR rules has to be declared with a [chr_constraint/1](practical.html#chr_constraint/1) declaration by the *constraint specifier*. For convenience multiple constraints may be declared at once with the same [chr_constraint/1](practical.html#chr_constraint/1) declaration followed by a comma-separated list of constraint specifiers.

A constraint specifier is, in its compact form, `F``/``A` where `F` and `A` are respectively the functor name and arity of the constraint, e.g.:

``` code
:- chr_constraint foo/1.
:- chr_constraint bar/2, baz/3.
```

In its extended form, a constraint specifier is `c``(``A_1``, ... ,``A_n``)` where `c` is the constraint's functor, `n` its arity and the `A_i` are argument specifiers. An argument specifier is a mode, optionally followed by a type. Example:

``` code
:- chr_constraint get_value(+,?).
:- chr_constraint domain(?int, +list(int)),
                  alldifferent(?list(int)).
```

**Modes**

A mode is one of:

**`-`**  
The corresponding argument of every occurrence of the constraint is always unbound.

**`+`**  
The corresponding argument of every occurrence of the constraint is always ground.

**`?`**  
The corresponding argument of every occurrence of the constraint can have any instantiation, which may change over time. This is the default value.

**Types**

A type can be a user-defined type or one of the built-in types. A type comprises a (possibly infinite) set of values. The type declaration for a constraint argument means that for every instance of that constraint the corresponding argument is only ever bound to values in that set. It does not state that the argument necessarily has to be bound to a value.

The built-in types are:

**int**  
The corresponding argument of every occurrence of the constraint is an integer number.

**dense_int**  
The corresponding argument of every occurrence of the constraint is an integer that can be used as an array index. Note that if this argument takes values in `[0,n]`, the array takes `O(n)` space.

**float**  
... a floating point number.

**number**  
... a number.

**natural**  
... a positive integer.

**any**  
The corresponding argument of every occurrence of the constraint can have any type. This is the default value.

&nbsp;

:- **chr_type**(`+TypeDeclaration`)  
User-defined types are algebraic data types, similar to those in Haskell or the discriminated unions in Mercury. An algebraic data type is defined using [chr_type/1](practical.html#chr_type/1):

``` code
:- chr_type type ---> body.
```

If the type term is a functor of arity zero (i.e. one having zero arguments), it names a monomorphic type. Otherwise, it names a polymorphic type; the arguments of the functor must be distinct type variables. The body term is defined as a sequence of constructor definitions separated by semi-colons.

Each constructor definition must be a functor whose arguments (if any) are types. Discriminated union definitions must be transparent: all type variables occurring in the body must also occur in the type.

Here are some examples of algebraic data type definitions:

``` code
:- chr_type color ---> red ; blue ; yellow ; green.

:- chr_type tree --->  empty ; leaf(int) ; branch(tree, tree).

:- chr_type list(T) ---> [] ; [T | list(T)].

:- chr_type pair(T1, T2) ---> (T1 - T2).
```

Each algebraic data type definition introduces a distinct type. Two algebraic data types that have the same bodies are considered to be distinct types (name equivalence).

Constructors may be overloaded among different types: there may be any number of constructors with a given name and arity, so long as they all have different types.

Aliases can be defined using ==. For example, if your program uses lists of lists of integers, you can define an alias as follows:

``` code
:- chr_type lli == list(list(int)).
```

**Type Checking**

Currently two complementary forms of type checking are performed:

1.  Static type checking is always performed by the compiler. It is limited to CHR rule heads and CHR constraint calls in rule bodies.

    Two kinds of type error are detected. The first is where a variable has to belong to two types. For example, in the program:

    ``` code
    :-chr_type foo ---> foo.
    :-chr_type bar ---> bar.

    :-chr_constraint abc(?foo).
    :-chr_constraint def(?bar).

    foobar @ abc(X) <=> def(X).
    ```

    the variable `X` has to be of both type `foo` and `bar`. This is reported as a type clash error:

    ``` code
    CHR compiler ERROR:
        `--> Type clash for variable _ in rule foobar:
                    expected type foo in body goal def(_, _)
                    expected type bar in head def(_, _)
    ```

    The second kind of error is where a functor is used that does not belong to the declared type. For example in:

    ``` code
    :- chr_type foo ---> foo.
    :- chr_type bar ---> bar.

    :- chr_constraint abc(?foo).

    foo @ abc(bar) <=> true.
    ```

    `bar` appears in the head of the rule where something of type `foo` is expected. This is reported as:

    ``` code
    CHR compiler ERROR:
        `--> Invalid functor in head abc(bar) of rule foo:
                    found `bar',
                    expected type `foo'!
    ```

    No runtime overhead is incurred in static type checking.

2.  Dynamic type checking checks at runtime, during program execution, whether the arguments of CHR constraints respect their declared types. The [when/2](coroutining.html#when/2) coroutining library is used to delay dynamic type checks until variables are instantiated.

    The kind of error detected by dynamic type checking is where a functor is used that does not belong to the declared type. For example, for the program:

    ``` code
    :-chr_type foo ---> foo.

    :-chr_constraint abc(?foo).
    ```

    we get the following error in an erroneous query:

    ``` code
    ?- abc(bar).
    ERROR: Type error: `foo' expected, found `bar'
           (CHR Runtime Type Error)
    ```

    Dynamic type checking is weaker than static type checking in the sense that it only checks the particular program execution at hand rather than all possible executions. It is stronger in the sense that it tracks types throughout the whole program.

    Note that it is enabled only in debug mode, as it incurs some (minor) runtime overhead.

### 9.3.3 CHR Compilation

The SWI-Prolog CHR compiler exploits [term_expansion/2](consulting.html#term_expansion/2) rules to translate the constraint handling rules to plain Prolog. These rules are loaded from the library `library(chr)`. They are activated if the compiled file has the `.chr` extension or after finding a declaration in the following format:

``` code
:- chr_constraint ...
```

It is advised to define CHR rules in a module file, where the module declaration is immediately followed by including the library(chr) library as exemplified below:

``` code
:- module(zebra, [ zebra/0 ]).
:- use_module(library(chr)).

:- chr_constraint ...
```

Using this style, CHR rules can be defined in ordinary Prolog .pl files and the operator definitions required by CHR do not leak into modules where they might cause conflicts.
