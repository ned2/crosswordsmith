
# D Glossary of Terms

**anonymous \[variable\]**  
The variable `_` is called the [anonymous](glossary.html#gloss:anonymou) variable. Multiple occurrences of `_` in a single [term](glossary.html#gloss:term) are not [shared](glossary.html#gloss:shared).

**arguments**  
Arguments are [terms](glossary.html#gloss:term) that appear in a [compound](glossary.html#gloss:compound) [term](glossary.html#gloss:term). `A1` and `a2` are the first and second argument of the term `myterm(A1, a2)`.

**arity**  
Argument count (= number of arguments) of a [compound](glossary.html#gloss:compound) [term](glossary.html#gloss:term).

**assert**  
Add a [clause](glossary.html#gloss:clause) to a [predicate](glossary.html#gloss:predicate). Clauses can be added at either end of the clause-list of a [predicate](glossary.html#gloss:predicate). See [asserta/1](db.html#asserta/1) and [assertz/1](db.html#assertz/1).

**atom**  
Textual constant. Used as name for [compound](glossary.html#gloss:compound) terms, to represent constants or text.

**backtracking**  
Search process used by Prolog. If a predicate offers multiple [clauses](glossary.html#gloss:clause) to solve a [goal](glossary.html#gloss:goal), they are tried one-by-one until one [succeeds](glossary.html#gloss:succeed). If a subsequent part of the proof is not satisfied with the resulting [variable](glossary.html#gloss:variable) [binding](glossary.html#gloss:binding), it may ask for an alternative [solution](glossary.html#gloss:solution) (= [binding](glossary.html#gloss:binding) of the [variables](glossary.html#gloss:variable)), causing Prolog to reject the previously chosen [clause](glossary.html#gloss:clause) and try the next one.

**binding \[of a variable\]**  
Current value of the [variable](glossary.html#gloss:variable). See also [backtracking](glossary.html#gloss:backtracking) and [query](glossary.html#gloss:query).

**built-in \[predicate\]**  
Predicate that is part of the Prolog system. Built-in predicates cannot be redefined by the user, unless this is overruled using [redefine_system_predicate/1](db.html#redefine_system_predicate/1).

**body**  
Part of a [clause](glossary.html#gloss:clause) behind the [neck](glossary.html#gloss:neck) operator (`:-`).

**choice point**  
A [choice point](glossary.html#gloss:choice-point) represents a choice in the search for a [solution](glossary.html#gloss:solution). Choice points are created if multiple clauses match a [query](glossary.html#gloss:query) or using disjunction ([;/2](control.html#;/2)). On [backtracking](glossary.html#gloss:backtracking), the execution state of the most recent [choice point](glossary.html#gloss:choice-point) is restored and search continues with the next alternative (i.e., next clause or second branch of [;/2](control.html#;/2)).

**clause**  
‘Sentence’of a Prolog program. A [clause](glossary.html#gloss:clause) consists of a [head](glossary.html#gloss:head) and [body](glossary.html#gloss:body) separated by the [neck](glossary.html#gloss:neck) operator (`:-`) or it is a [fact](glossary.html#gloss:fact). For example:

``` code
parent(X) :-
        father(X, _).
```

Expressed as “X is a parent if X is a father of someone” . See also [variable](glossary.html#gloss:variable) and [predicate](glossary.html#gloss:predicate).

**compile**  
Process where a Prolog [program](glossary.html#gloss:program) is translated to a sequence of instructions. See also [interpreted](glossary.html#gloss:interpreted). SWI-Prolog always compiles your program before executing it.

**compound \[term\]**  
Also called [structure](glossary.html#gloss:structure). It consists of a name followed by `N` [arguments](glossary.html#gloss:argument), each of which are [terms](glossary.html#gloss:term). `N` is called the [arity](glossary.html#gloss:arity) of the term.

**context module**  
If a [term](glossary.html#gloss:term) is referring to a [predicate](glossary.html#gloss:predicate) in a [module](glossary.html#gloss:module), the [context module](glossary.html#gloss:context-module) is used to find the target module. The context module of a [goal](glossary.html#gloss:goal) is the module in which the [predicate](glossary.html#gloss:predicate) is defined, unless this [predicate](glossary.html#gloss:predicate) is [module transparent](glossary.html#gloss:module-transparent), in which case the [context module](glossary.html#gloss:context-module) is inherited from the parent [goal](glossary.html#gloss:goal). See also [module_transparent/1](ctxmodule.html#module_transparent/1) and [meta-predicate](glossary.html#gloss:meta-predicate).

**dcg**  
Abbreviation for **Definite Clause Grammar**.

**det \[determinism\]**  
Short for [deterministic](glossary.html#gloss:deterministic).

**determinism**  
How many solutions a [goal](glossary.html#gloss:goal) can provide. Values are‘nondet’(zero to infinite),‘multi’(one to infinite),‘det’(exactly one) and‘semidet’(zero or one).

**deterministic**  
A [predicate](glossary.html#gloss:predicate) is [deterministic](glossary.html#gloss:deterministic) if it succeeds exactly one time without leaving a [choice point](glossary.html#gloss:choice-point).

**dynamic \[predicate\]**  
A [dynamic](glossary.html#gloss:dynamic) predicate is a predicate to which [clauses](glossary.html#gloss:clause) may be [assert](glossary.html#gloss:assert)ed and from which [clauses](glossary.html#gloss:clause) may be [retract](glossary.html#gloss:retract)ed while the program is running. See also [update view](glossary.html#gloss:update-view).

**exported \[predicate\]**  
A [predicate](glossary.html#gloss:predicate) is said to be [exported](glossary.html#gloss:exported) from a [module](glossary.html#gloss:module) if it appears in the [public list](glossary.html#gloss:public-list). This implies that the predicate can be [imported](glossary.html#gloss:imported) into another module to make it visible there. See also [use_module/\[1,2\]](import.html#use_module/1).

**fact**  
[Clause](glossary.html#gloss:clause) without a [body](glossary.html#gloss:body). This is called a fact because, interpreted as logic, there is no condition to be satisfied. The example below states `john` is a person.

``` code
person(john).
```

**fail**  
A [goal](glossary.html#gloss:goal) is said to have failed if it could not be [proven](glossary.html#gloss:prove).

**float**  
Computer's crippled representation of a real number. Represented as‘IEEE double’.

**foreign**  
Computer code expressed in languages other than Prolog. SWI-Prolog can only cooperate directly with the C and C++ computer languages.

**functor**  
Combination of name and [arity](glossary.html#gloss:arity) of a [compound](glossary.html#gloss:compound) term. The term `foo(a, b, c)` is said to be a term belonging to the functor foo/3 . foo/0 is used to refer to the [atom](glossary.html#gloss:atom) `foo`.

**goal**  
Question stated to the Prolog engine. A [goal](glossary.html#gloss:goal) is either an [atom](glossary.html#gloss:atom) or a [compound](glossary.html#gloss:compound) term. A [goal](glossary.html#gloss:goal) either succeeds, in which case the [variables](glossary.html#gloss:variable) in the [compound](glossary.html#gloss:compound) terms have a [binding](glossary.html#gloss:binding), or it [fails](glossary.html#gloss:fail) if Prolog fails to prove it.

**hashing**  
[Indexing](glossary.html#gloss:indexing) technique used for quick lookup.

**head**  
Part of a [clause](glossary.html#gloss:clause) before the [neck](glossary.html#gloss:neck) operator (`:-`). This is an [atom](glossary.html#gloss:atom) or [compound](glossary.html#gloss:compound) term.

**imported \[predicate\]**  
A [predicate](glossary.html#gloss:predicate) is said to be [imported](glossary.html#gloss:imported) into a [module](glossary.html#gloss:module) if it is defined in another [module](glossary.html#gloss:module) and made available in this [module](glossary.html#gloss:module). See also [chapter 6](modules.html#sec:6).

**indexing**  
Indexing is a technique used to quickly select candidate [clauses](glossary.html#gloss:clause) of a [predicate](glossary.html#gloss:predicate) for a specific [goal](glossary.html#gloss:goal). In most Prolog systems, indexing is done (only) on the first [argument](glossary.html#gloss:argument) of the [head](glossary.html#gloss:head). If this argument is instantiated to an [atom](glossary.html#gloss:atom), [integer](glossary.html#gloss:integer), [float](glossary.html#gloss:float) or [compound](glossary.html#gloss:compound) term with [functor](glossary.html#gloss:functor), [hashing](glossary.html#gloss:hashing) is used to quickly select all [clauses](glossary.html#gloss:clause) where the first argument may [unify](glossary.html#gloss:unify) with the first argument of the [goal](glossary.html#gloss:goal). SWI-Prolog supports just-in-time and multi-argument indexing. See [section 2.17](jitindex.html#sec:2.17).

**integer**  
Whole number. On all implementations of SWI-Prolog integers are at least 64-bit signed values. When linked to the GNU GMP library, integer arithmetic is unbounded. See also [current_prolog_flag/2](flags.html#current_prolog_flag/2), flags [bounded](flags.html#flag:bounded), [max_integer](flags.html#flag:max_integer) and [min_integer](flags.html#flag:min_integer).

**interpreted**  
As opposed to [compiled](glossary.html#gloss:compile), interpreted means the Prolog system attempts to prove a [goal](glossary.html#gloss:goal) by directly reading the [clauses](glossary.html#gloss:clause) rather than executing instructions from an (abstract) instruction set that is not or only indirectly related to Prolog.

**instantiation \[of an argument\]**  
To what extend a term is bound to a value. Typical levels are‘unbound’(a [variable](glossary.html#gloss:variable)),‘ground’(term without variables) or‘partially bound’(term with embedded variables).

**meta-predicate**  
A [predicate](glossary.html#gloss:predicate) that reasons about other [predicates](glossary.html#gloss:predicate), either by calling them, (re)defining them or querying [properties](glossary.html#gloss:property).

**mode \[declaration\]**  
Declaration of an argument [instantiation](glossary.html#gloss:instantiation) pattern for a [predicate](glossary.html#gloss:predicate), often accompanied with a [determinism](glossary.html#gloss:determinism).

**module**  
Collection of predicates. Each module defines a name-space for predicates. [built-in](glossary.html#gloss:built-in) predicates are accessible from all modules. Predicates can be published ([exported](glossary.html#gloss:exported)) and [imported](glossary.html#gloss:imported) to make their definition available to other modules.

**module transparent \[predicate\]**  
A [predicate](glossary.html#gloss:predicate) that does not change the [context module](glossary.html#gloss:context-module). Sometimes also called a [meta-predicate](glossary.html#gloss:meta-predicate).

**multi \[determinism\]**  
A [predicate](glossary.html#gloss:predicate) is said to have [determinism](glossary.html#gloss:determinism) multi if it generates at *least* one answer.

**multifile \[predicate\]**  
Predicate for which the definition is distributed over multiple source files. See [multifile/1](dynamic.html#multifile/1).

**neck**  
Operator (`:-`) separating [head](glossary.html#gloss:head) from [body](glossary.html#gloss:body) in a [clause](glossary.html#gloss:clause).

**nondet**  
Short for [non deterministic](glossary.html#gloss:non-deterministic).

**non deterministic**  
A [non deterministic](glossary.html#gloss:non-deterministic) predicate is a predicate that may fail or succeed any number of times.

**operator**  
Symbol ([atom](glossary.html#gloss:atom)) that may be placed before its [operand](glossary.html#gloss:operand) (prefix), after its [operand](glossary.html#gloss:operand) (postfix) or between its two [operands](glossary.html#gloss:operand) (infix).

In Prolog, the expression `a+b` is exactly the same as the canonical term `+(a,b)`.

**operand**  
[Argument](glossary.html#gloss:argument) of an [operator](glossary.html#gloss:operator).

**precedence**  
The [priority](glossary.html#gloss:priority) of an [operator](glossary.html#gloss:operator). Operator precedence is used to interpret `a+b*c` as `+(a, *(b,c))`.

**predicate**  
Collection of [clauses](glossary.html#gloss:clause) with the same [functor](glossary.html#gloss:functor) (name/[arity](glossary.html#gloss:arity)). If a [goal](glossary.html#gloss:goal) is proved, the system looks for a [predicate](glossary.html#gloss:predicate) with the same functor, then uses [indexing](glossary.html#gloss:indexing) to select candidate [clauses](glossary.html#gloss:clause) and then tries these [clauses](glossary.html#gloss:clause) one-by-one. See also [backtracking](glossary.html#gloss:backtracking).

**predicate indicator**  
Term of the form Name/Arity (traditional) or Name//Arity (ISO DCG proposal), where Name is an atom and Arity a non-negative integer. It acts as an *indicator* (or reference) to a predicate or [DCG](glossary.html#gloss:dcg) rule.

**priority**  
In the context of [operators](glossary.html#gloss:operator) a synonym for [precedence](glossary.html#gloss:precedence).

**program**  
Collection of [predicates](glossary.html#gloss:predicate).

**property**  
Attribute of an object. SWI-Prolog defines various *\*\_property* predicates to query the status of predicates, clauses. etc.

**prove**  
Process where Prolog attempts to prove a [query](glossary.html#gloss:query) using the available [predicates](glossary.html#gloss:predicate).

**public list**  
List of [predicates](glossary.html#gloss:predicate) exported from a [module](glossary.html#gloss:module).

**query**  
See [goal](glossary.html#gloss:goal).

**retract**  
Remove a [clause](glossary.html#gloss:clause) from a [predicate](glossary.html#gloss:predicate). See also [dynamic](glossary.html#gloss:dynamic), [update view](glossary.html#gloss:update-view) and [assert](glossary.html#gloss:assert).

**semidet**  
Shorthand for

**semi deterministic**  
.

**semi deterministic**  
A [predicate](glossary.html#gloss:predicate) that is [semi deterministic](glossary.html#gloss:semi-deterministic) either fails or succeeds exactly once without a [choice point](glossary.html#gloss:choice-point). See also [deterministic](glossary.html#gloss:deterministic).

**shared**  
Two [variables](glossary.html#gloss:variable) are called [shared](glossary.html#gloss:shared) after they are [unified](glossary.html#gloss:unify). This implies if either of them is [bound](glossary.html#gloss:binding), the other is bound to the same value:

``` code
?- A = B, A = a.
A = B, B = a.
```

**singleton \[variable\]**  
[Variable](glossary.html#gloss:variable) appearing only one time in a [clause](glossary.html#gloss:clause). SWI-Prolog normally warns for this to avoid you making spelling mistakes. If a variable appears on purpose only once in a clause, write it as `_` (see [anonymous](glossary.html#gloss:anonymou)). Rules for naming a variable and avoiding a warning are given in [section 2.15.1.10](syntax.html#sec:2.15.1.10).

**solution**  
[Bindings](glossary.html#gloss:binding) resulting from a successfully [prove](glossary.html#gloss:prove)n [goal](glossary.html#gloss:goal).

**structure**  
Synonym for [compound](glossary.html#gloss:compound) term.

**string**  
Used for the following representations of text: a packed array (see [section 5.2](string.html#sec:5.2), SWI-Prolog specific), a list of character codes or a list of one-character [atoms](glossary.html#gloss:atom).

**succeed**  
A [goal](glossary.html#gloss:goal) is said to have [succeeded](glossary.html#gloss:succeed) if it has been [proven](glossary.html#gloss:prove).

**term**  
Value in Prolog. A [term](glossary.html#gloss:term) is either a [variable](glossary.html#gloss:variable), [atom](glossary.html#gloss:atom), [integer](glossary.html#gloss:integer), [float](glossary.html#gloss:float) or [compound](glossary.html#gloss:compound) term. In addition, SWI-Prolog also defines the type [string](glossary.html#gloss:string).

**transparent**  
See [module transparent](glossary.html#gloss:module-transparent).

**unify**  
Prolog process to make two terms equal by assigning variables in one term to values at the corresponding location of the other term. For example:

``` code
?- foo(a, B) = foo(A, b).
A = a,
B = b.
```

Unlike assignment (which does not exist in Prolog), unification is not directed.

**update view**  
How Prolog behaves when a [dynamic](glossary.html#gloss:dynamic) [predicate](glossary.html#gloss:predicate) is changed while it is running. There are two models. In most older Prolog systems the change becomes immediately visible to the [goal](glossary.html#gloss:goal), in modern systems including SWI-Prolog, the running [goal](glossary.html#gloss:goal) is not affected. Only new [goals](glossary.html#gloss:goal) ‘see’the new definition.

**variable**  
A Prolog variable is a value that‘is not yet bound’. After [binding](glossary.html#gloss:binding) a variable, it cannot be modified. [Backtracking](glossary.html#gloss:backtracking) to a point in the execution before the variable was bound will turn it back into a variable:

``` code
?- A = b, A = c.
false.

?- (A = b; true; A = c).
A = b ;
true ;
A = c .
```

See also [unify](glossary.html#gloss:unify).
