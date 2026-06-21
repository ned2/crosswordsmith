
## 4.15 Declaring predicate properties

This section describes directives which manipulate attributes of predicate definitions. The functors [dynamic/1](dynamic.html#dynamic/1), [multifile/1](dynamic.html#multifile/1), [discontiguous/1](dynamic.html#discontiguous/1) and [public/1](dynamic.html#public/1) are operators of priority 1150 (see [op/3](operators.html#op/3)), which implies that the list of predicates they involve can just be a comma-separated list:

``` code
:- dynamic
        foo/0,
        baz/2.
```

In SWI-Prolog all these directives are just predicates. This implies they can also be called by a program. Do not rely on this feature if you want to maintain portability to other Prolog implementations.

Notably with the introduction of tabling (see [section 7](tabling.html#sec:7)) it is common that a set of predicates require multiple options to be set. SWI-Prolog offers two mechanisms to cope with this. The predicate [dynamic/2](dynamic.html#dynamic/2) can be used to make a list of predicates dynamic and set additional options. In addition and for compatibility with XSB,^(92Note that `as` is in XSB a high-priority operator and in SWI a low-priority and therefore both the sets of predicate indicators as multiple options require parenthesis.) all the predicates below accept a term `as((:PredicateIndicator, ... ), (+Options))`, where `Options` is a *comma-list* of one of more of the following options:

**incremental**  
Include a dynamic predicate into the incremental tabling dependency graph. See [section 7.7](tabling-incremental.html#sec:7.7).

**opaque**  
Opposite of `incremental`. For XSB compatibility.^(93In XSB, `opaque` is distinct from the default in the sense that dynamic switching between `opaque` and `incremental` is allowed.)

**abstract**(`Level`)  
Used together with `incremental` to reduce the dependency graph. See [section 7.7](tabling-incremental.html#sec:7.7).

**volatile**  
Do not save this predicate. See [volatile/1](saved-states.html#volatile/1).

**multifile**  
Predicate may have clauses in multiple clauses. See [multifile/1](dynamic.html#multifile/1).

**discontiguous**  
Predicate clauses may not be contiguous in the file. See [discontiguous/1](dynamic.html#discontiguous/1).

**shared**  
Dynamic predicate is shared between all threads. This is currently the default.

**local**  
**private**  
Dynamic predicate has distinct set of clauses in each thread. See [thread_local/1](threadcom.html#thread_local/1).

Below are some examples, where the last two are semantically identical.

``` code
:- dynamic person/2 as incremental.
:- dynamic (person/2,organization/2) as (incremental, abstract(0)).
:- dynamic([ person/2,
             organization/2
           ],
           [ incremental(true),
             abstract(0)
           ]).
```

\[ISO\]**dynamic** `:PredicateIndicator, ...`  
Informs the interpreter that the definition of the predicate(s) may change during execution (using [assert/1](db.html#assert/1) and/or [retract/1](db.html#retract/1)). In the multithreaded version, the clauses of dynamic predicates are shared between the threads. The directive [thread_local/1](threadcom.html#thread_local/1) provides an alternative where each thread has its own clause list for the predicate. Dynamic predicates can be turned into static ones using [compile_predicates/1](dynamic.html#compile_predicates/1).

**dynamic**(`:ListOfPredicateIndicators, +Options`)  
As [dynamic/1](dynamic.html#dynamic/1), but allows for setting additional properties. This predicate allows for setting multiple properties on multiple predicates in a single call. SWI-Prolog also offers the XSB compatible `:- dynamic (``p/1) as (incremental,abstract(0)).` syntax. See the introduction of [section 4.15](dynamic.html#sec:4.15). Defined `Options` are:

**incremental**(`+Boolean`)  
Make the dynamic predicate signal depending *tables*. See [section 7.7](tabling-incremental.html#sec:7.7).

**abstract**(`0`)  
This option must be used together with `incremental`. The only supported value is `0`. With this option a call to the incremental dynamic predicate is recorded as the most generic term for the predicate rather than the specific variant.

**thread**(`+Local`)  
`Local` is one of `shared` (default) or `local`. See also [thread_local/1](threadcom.html#thread_local/1).

**multifile**(`+Boolean`)  
**discontiguous**(`+Boolean`)  
**volatile**(`+Boolean`)  
Set the corresponding property. See [multifile/1](dynamic.html#multifile/1), [discontiguous/1](dynamic.html#discontiguous/1) and [volatile/1](saved-states.html#volatile/1).

**mode** `+Head, ...`  
Define the *modes* for the predicates referenced by the comma-separated list `Head`. Each argument of each head is a *mode specifier*. Defined specifiers are `+`, `-` and `?`. Mode declarations have a long history in Prolog. This predicate uses the definitions derived from e.g., Quintus Prolog and used for documentation in the ISO standard. The specifiers declare whether or not an argument is *unbound* at *call time*. The `+` declares that the argument is [nonvar/1](typetest.html#nonvar/1), the `-` declares the argument is [var/1](typetest.html#var/1) and `?` makes no claim.

Several Prolog systems define `mode` as an operator using the declaration below. Currently we do not define the operator.

``` code
:- op(1150, fx, mode).
```

SWI-Prolog uses the mode information for its just-in-time clause indexing as described in [section 2.17](jitindex.html#sec:2.17). JIT only uses the `-` specifier to avoid examining that argument as a candidate for indexing.

**compile_predicates**(`:ListOfPredicateIndicators`)  
Compile a list of specified dynamic predicates (see [dynamic/1](dynamic.html#dynamic/1) and [assert/1](db.html#assert/1)) into normal static predicates. This call tells the Prolog environment the definition will not change anymore and further calls to [assert/1](db.html#assert/1) or [retract/1](db.html#retract/1) on the named predicates raise a permission error. This predicate is designed to deal with parts of the program that are generated at runtime but do not change during the remainder of the program execution.^(94The specification of this predicate is from Richard O'Keefe. The implementation is allowed to optimise the predicate. This is not yet implemented. In multithreaded Prolog, however, static code runs faster as it does not require synchronisation. This is particularly true on SMP hardware.)

\[ISO\]**multifile** `:PredicateIndicator, ...`  
Informs the system that the specified predicate(s) may be defined over more than one file. This stops [consult/1](consulting.html#consult/1) from redefining a predicate when a new definition is found.

\[ISO\]**discontiguous** `:PredicateIndicator, ...`  
Informs the system that the clauses of the specified predicate(s) might not be together in the source file. See also [style_check/1](debugger.html#style_check/1).

**public** `:PredicateIndicator, ...`  
Instructs the cross-referencer that the predicate can be called. It has no semantics.^(95This declaration is compatible with SICStus. In YAP, [public/1](dynamic.html#public/1) instructs the compiler to keep the source. As the source is always available in SWI-Prolog, our current interpretation also enhances the compatibility with YAP.) The public declaration can be queried using [predicate_property/2](examineprog.html#predicate_property/2). The [public/1](dynamic.html#public/1) directive does *not* export the predicate (see [module/1](mtoplevel.html#module/1) and [export/1](altmoduleapi.html#export/1)). The public directive is used for (1) direct calls into the module from, e.g., foreign code, (2) direct calls into the module from other modules, or (3) flag a predicate as being called if the call is generated by meta-calling constructs that are not analysed by the cross-referencer.

**non_terminal** `:PredicateIndicator, ...`  
Sets the `non_terminal` property on the predicate. This indicates that the predicate implements a *grammar rule*. See [predicate_property/2](examineprog.html#predicate_property/2). The `non_terminal` property is set for predicates exported as `Name`//`Arity` as well as predicates that have at least one clause written using the `-->``/2` notation.
