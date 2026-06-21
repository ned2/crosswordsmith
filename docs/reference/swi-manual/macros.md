
## A.26 library(macros): Macro expansion

This library defines a macro expansion mechanism that operates on arbitrary terms. Unlike [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2), a term is explicitly designed for expansion using the term `#(Macro)`. Macros are first of all intended to deal with compile time constants. They can also be used to construct terms at compile time.

### A.26.1 Defining and using macros

Macros are defined for the current module using one of the three constructs below.

``` code
#define(Macro, Replacement).
#define(Macro, Replacement) :- Code.
#import(ModuleFile).
```

`Macro` is a *callable term*, not being `define(_,_)`, or `import(_)`. `Replacement` is an arbitrary Prolog term. `Code` is a Prolog *body term* that *must* succeed and can be used to dynamically generate (parts of) `Replacement`.

The `#import(ModuleFile)` definition makes all macros from the given module available for expansion in the module it appears. Normally this shall be appear after local macro definitions.

A macro is called using the term `#(Macro)`. `#` is defined as a low-priority (10) prefix operator to allow for `#Macro`. Macros can appear at the following places:

- An entire sentence (clause)
- Any argument of a compound. This implies also the head and body of a clause.
- Anywhere in a list, including as the tail of a list
- As a value for a dict key or as a dict key name.

Macros can **not** appear as name of a compound or tag of a dict. A term `#Macro` appearing in one of the allowed places **must** have a matching macro defined, i.e., `#Macro` is **always** expanded. An error is emitted if the expansion fails. Macro expansion is applied recursively and thus, macros may be passed to macro arguments and macro expansion may use other macros.

Macros are matched to terms using *Single Sided Unification* (SSU), implemented using `Head => Body` rules. This implies that the matching never instantiates variables in the term that is being expanded.

Below are some examples. The first line defines the macro and the indented line after show example usage of the macro.

``` code
#define(max_width, 100).
    W < #max_width

#define(calc(Expr), Value) :- Value is Expr.
    fact(#calc(#max_width*2)).

#define(pt(X,Y), point{x:X, y:Y}).
    reply_json(json{type:polygon,
                    points:[#pt(0,0), #pt(0,5), #pt(5,0)]}).
```

Macro expansion expands terms `#(Callable)`. If the argument to the \#-term is not a `callable`, the \#-term is not modified. This notably allows for `#(Var)` as used by `library(clpfd)` to indicate that a variable is constraint to be an (`clp(fd)`) integer.

### A.26.2 Implementation details

A macro `#define(Macro, Expanded) :- Body.` is, after some basic sanity checks, translated into a rule

``` code
'$macro'(Macro, Var), Body => Var = Expanded.
```

The `#import(File)` is translated into `:- use_module(File, [])` and a *link clause* that links the macro expansion from the module defined in `File` to the current module.

Macro expansion is realised by creating a clause for [term_expansion/2](consulting.html#term_expansion/2) in the current module. This clause results from expanding the first `#define` or `#import` definition. Thus, if macros are defined before any other local definition for [term_expansion/2](consulting.html#term_expansion/2) it is executed as the first step. The macro expansion fails if no macros were encounted in the term, allowing other term_expansion rules local to the module to take effect. In other words, a term holding macros is not subject to any other term expansion local to the module. It is subject to term expansion defined in module `user` and `system` that is performed after the local expansion is completed.

### A.26.3 Predicates

\[semidet\]**include_macros**(`+M, +Macro, -Expanded`)  
Include macros from another module. This predicate is a helper for `#import(File)`. It calls’\$macro’/2 in `M`, but fails silently in case `Macro` is not defined in `M` as it may be defined in another imported macro file or further down in the current file.

\[semidet\]**expand_macros**(`+Module, +TermIn, -TermOut, +PosIn, -PosOut`)  
Perform macro expansion on `TermIn` with layout `PosIn` to produce `TermOut` with layout `PosOut`. The transformation is performed if the current load context module is `Module` (see [prolog_load_context/2](consulting.html#prolog_load_context/2)).

This predicate is not intended for direct usage.

\[det\]**macro_position**(`-Position`)  
True when `Position` is the position of the macro. `Position` is a term `File:Line:LinePos`. If `File` is unknown it is unified with `-`. If Line and/or LinePos are unknown they are unified with 0. This predicate can be used in the body of a macro definition to provide the source location. The example below defines `#pp(Var)` to print a variable together with the variable name and source location.

``` code
#define(pp(Var), print_message(debug, dump_var(Pos, Name, Var))) :-
    (   var_property(Var, name(Name))
    ->  true
    ;   Name = 'Var'
    ),
    macro_position(Pos).

:- multifile prolog:message//1.
prolog:message(dump_var(Pos,Name,Var)) -->
    [ url(Pos), ': ',
      ansi([fg(magenta),bold], '~w', [Name]), ' = ',
      ansi(code, '~p', [Var])
    ].
```
