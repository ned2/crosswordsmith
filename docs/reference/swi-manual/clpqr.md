
## A.10 library(clpqr): Constraint Logic Programming over Rationals and Reals

> Author: *Christian Holzbaur*, ported to SWI-Prolog by *Leslie De Koninck*, K.U. Leuven

This CLP(Q,R) system is a port of the CLP(Q,R) system of Sicstus Prolog by Christian Holzbaur: Holzbaur C.: OFAI clp(q,r) Manual, Edition 1.3.3, Austrian Research Institute for Artificial Intelligence, Vienna, TR-95-09, 1995.^(252http://www.ai.univie.ac.at/cgi-bin/tr-online?number+95-09) This manual is roughly based on the manual of the above mentioned CLP(Q,R) implementation.

The CLP(Q,R) system consists of two components: the CLP(Q) library for handling constraints over the rational numbers and the CLP(R) library for handling constraints over the real numbers (using floating point numbers as representation). Both libraries offer the same predicates (with exception of [bb_inf/4](clpqr.html#bb_inf/4) in CLP(Q) and [bb_inf/5](clpqr.html#bb_inf/5) in CLP(R)). It is allowed to use both libraries in one program, but using both CLP(Q) and CLP(R) constraints on the same variable will result in an exception.

Please note that the `library(clpqr)` library is *not* an *autoload* library and therefore this library must be loaded explicitly before using it:

``` code
:- use_module(library(clpq)).
```

or

``` code
:- use_module(library(clpr)).
```

### A.10.1 Solver predicates

The following predicates are provided to work with constraints:

**{}**(`+Constraints`)  
Adds the constraints given by `Constraints` to the constraint store.

**entailed**(`+Constraint`)  
Succeeds if `Constraint` is necessarily true within the current constraint store. This means that adding the negation of the constraint to the store results in failure.

**inf**(`+Expression, -Inf`)  
Computes the infimum of `Expression` within the current state of the constraint store and returns that infimum in `Inf`. This predicate does not change the constraint store.

**sup**(`+Expression, -Sup`)  
Computes the supremum of `Expression` within the current state of the constraint store and returns that supremum in `Sup`. This predicate does not change the constraint store.

**minimize**(`+Expression`)  
Minimizes `Expression` within the current constraint store. This is the same as computing the infimum and equating the expression to that infimum.

**maximize**(`+Expression`)  
Maximizes `Expression` within the current constraint store. This is the same as computing the supremum and equating the expression to that supremum.

**bb_inf**(`+Ints, +Expression, -Inf, -Vertex, +Eps`)  
This predicate is offered in CLP(R) only. It computes the infimum of `Expression` within the current constraint store, with the additional constraint that in that infimum, all variables in `Ints` have integral values. `Vertex` will contain the values of `Ints` in the infimum. `Eps` denotes how much a value may differ from an integer to be considered an integer. E.g. when `Eps` = 0.001, then X = 4.999 will be considered as an integer (5 in this case). `Eps` should be between 0 and 0.5.

**bb_inf**(`+Ints, +Expression, -Inf, -Vertex`)  
This predicate is offered in CLP(Q) only. It behaves the same as [bb_inf/5](clpqr.html#bb_inf/5) but does not use an error margin.

**bb_inf**(`+Ints, +Expression, -Inf`)  
The same as [bb_inf/5](clpqr.html#bb_inf/5) or [bb_inf/4](clpqr.html#bb_inf/4) but without returning the values of the integers. In CLP(R), an error margin of 0.001 is used.

**dump**(`+Target, +Newvars, -CodedAnswer`)  
Returns the constraints on `Target` in the list `CodedAnswer` where all variables of `Target` have been replaced by `NewVars`. This operation does not change the constraint store. E.g. in

``` code
dump([X,Y,Z],[x,y,z],Cons)
```

`Cons` will contain the constraints on X, Y and Z, where these variables have been replaced by atoms x, y and z.

### A.10.2 Syntax of the predicate arguments

The arguments of the predicates defined in the subsection above are defined in [table 10](clpqr.html#tab:clpqrbnf). Failing to meet the syntax rules will result in an exception.

|  |  |  |  |
|----|---:|----|----|
| \<`Constraints`\> | ::= | \<`Constraint`\> | single constraint |
|  | \| | \<`Constraint`\> , \<`Constraints`\> | conjunction |
|  | \| | \<`Constraint`\> ; \<`Constraints`\> | disjunction |
| \<`Constraint`\> | ::= | \<`Expression`\> `<` \<`Expression`\> | less than |
|  | \| | \<`Expression`\> `>` \<`Expression`\> | greater than |
|  | \| | \<`Expression`\> `=<` \<`Expression`\> | less or equal |
|  | \| | `<=`(\<`Expression`\>, \<`Expression`\>) | less or equal |
|  | \| | \<`Expression`\> `>=` \<`Expression`\> | greater or equal |
|  | \| | \<`Expression`\> `=\=` \<`Expression`\> | not equal |
|  | \| | \<`Expression`\> =:= \<`Expression`\> | equal |
|  | \| | \<`Expression`\> = \<`Expression`\> | equal |
| \<`Expression`\> | ::= | \<`Variable`\> | Prolog variable |
|  | \| | \<`Number`\> | Prolog number |
|  | \| | +\<`Expression`\> | unary plus |
|  | \| | -\<`Expression`\> | unary minus |
|  | \| | \<`Expression`\> + \<`Expression`\> | addition |
|  | \| | \<`Expression`\> - \<`Expression`\> | substraction |
|  | \| | \<`Expression`\> \* \<`Expression`\> | multiplication |
|  | \| | \<`Expression`\> / \<`Expression`\> | division |
|  | \| | abs(\<`Expression`\>) | absolute value |
|  | \| | sin(\<`Expression`\>) | sine |
|  | \| | cos(\<`Expression`\>) | cosine |
|  | \| | tan(\<`Expression`\>) | tangent |
|  | \| | exp(\<`Expression`\>) | exponent |
|  | \| | pow(\<`Expression`\>) | exponent |
|  | \| | \<`Expression`\> `^` \<`Expression`\> | exponent |
|  | \| | min(\<`Expression`\>, \<`Expression`\>) | minimum |
|  | \| | max(\<`Expression`\>, \<`Expression`\>) | maximum |

**Table 10 :** CLP(Q,R) constraint BNF

### A.10.3 Use of unification

Instead of using the [{}/1](clpqr.html#%7B%7D/1) predicate, you can also use the standard unification mechanism to store constraints. The following code samples are equivalent:

- *Unification with a variable*  

  ``` code
  {X =:= Y}
  {X = Y}
  X = Y
  ```

- *Unification with a number*  

  ``` code
  {X =:= 5.0}
  {X = 5.0}
  X = 5.0
  ```

### A.10.4 Non-linear constraints

The CLP(Q,R) system deals only passively with non-linear constraints. They remain in a passive state until certain conditions are satisfied. These conditions, which are called the isolation axioms, are given in [table 11](clpqr.html#tab:clpqraxioms).

|                 |                           |                            |
|-----------------|---------------------------|----------------------------|
| `A = B * C`     | B or C is ground          | A = 5 \* C or A = B \* 4   |
|                 | A and (B or C) are ground | 20 = 5 \* C or 20 = B \* 4 |
| `A = B / C`     | C is ground               | A = B / 3                  |
|                 | A and B are ground        | 4 = 12 / C                 |
| `X = min(Y,Z)`  | Y and Z are ground        | X = min(4,3)               |
| `X = max(Y,Z)`  | Y and Z are ground        | X = max(4,3)               |
| `X = abs(Y)`    | Y is ground               | X = abs(-7)                |
| `X = pow(Y,Z)`  | X and Y are ground        | 8 = 2 `^` Z                |
| `X = exp(Y,Z)`  | X and Z are ground        | 8 = Y `^` 3                |
| `X = Y` `^` `Z` | Y and Z are ground        | X = 2 `^` 3                |
| `X = sin(Y)`    | X is ground               | 1 = sin(Y)                 |
| `X = cos(Y)`    | Y is ground               | X = sin(1.5707)            |
| `X = tan(Y)`    |                           |                            |

**Table 11 :** CLP(Q,R) isolating axioms

### A.10.5 Status and known problems

The clpq and clpr libraries are‘orphaned’, i.e., they currently have no maintainer.

- *Top-level output*  
  The top-level output may contain variables not present in the original query:

  ``` code
  ?- {X+Y>=1}.
  {Y=1-X+_G2160, _G2160>=0}.

  ?-
  ```

  Nonetheless, for linear constraints this kind of answer means unconditional satisfiability.

- *Dumping constraints*  
  The first argument of [dump/3](clpqr.html#dump/3) has to be a list of free variables at call-time:

  ``` code
  ?- {X=1},dump([X],[Y],L).
  ERROR: Unhandled exception: Unknown message:
         instantiation_error(dump([1],[_G11],_G6),1)
  ?-
  ```
