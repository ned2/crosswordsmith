
## 2.17 Just-in-time clause indexing

SWI-Prolog provides‘just-in-time’indexing over multiple arguments.‘Just-in-time’means that clause indexes are not built by the compiler (or [asserta/1](db.html#asserta/1) for dynamic predicates), but on the first call to such a predicate where an index might help (i.e., a call where at least one argument is instantiated). This section describes the rules used by the indexing logic. Note that this logic is not‘set in stone’. The indexing capabilities of the system will change. Although this inevitably leads to some regressing on some particular use cases, we strive to avoid significant slowdowns.

The list below describes the clause selection process for various predicates and calls. The alternatives are considered in the order they are presented.

- *Special purpose code*  
  Currently two special cases are recognised by the compiler: static code with exactly one clause and static code with two clauses, one where some argument is the empty list (`[]`) and one where the same argument is a non-empty list (`[_|_]`). Note that if this argument is an *output argument*, semantics are preserved, but efficiency suffers slightly. This slow-down can be avoided using [mode/1](dynamic.html#mode/1) to declare the argument as `-`.^(36Up to version 9.3.18, only the first argument was considered.)

- *Linear scan on primary index argument*  
  The principal clause list maintains a *key*, normally for the first argument. An indexing key is either a constant or a functor (name/arity reference). Calls with an instantiated primary index argument and less than 10 clauses perform a linear scan for a possible matching clause using this index key. If the result is deterministic it is used. Otherwise the system looks for better indexes.^(37Up to 7.7.2 this result was used also when non-deterministic.). The primary index argument is the first argument for which at least one of the clauses has an indexable (nonvar) value for that argument.^(38Up to version 9.3.18, the primary index was fixed to the first argument.)

- *Hash lookup*  
  If none of the above applies, the system considers the available hash tables for which the corresponding argument is instantiated. If a table is found with acceptable characteristics, it is used. Otherwise it assesses the clauses for all instantiated arguments and selects the best candidate for creating a new hash table. If there is no single argument that provides an acceptable hash quality it will search for a combination of arguments.^(39The last step was added in SWI-Prolog 7.5.8.) Searching for index candidates is only performed on the first 254 arguments.

  If a single-argument index contains multiple compound terms with the same name and arity and at least one non-variable argument, a *list index* is created. A subsequent query where this argument is bound to a compound causes jiti indexing to be applied *recursively* on the arguments of the term. This is called *deep indexing*.^(40Deep indexing was added in version 7.7.4.) See also [section 2.17.1](jitindex.html#sec:2.17.1)

  Clauses that have a variable at an otherwise indexable argument must be linked into all hash buckets. Currently, predicates that have more than 10% such clauses for a specific argument are not considered for indexing on that argument.

  Disregarding variables, the suitability of an argument for hashing is expressed as the number of unique indexable values divided by the standard deviation of the number of duplicate values for each value plus one.^(41Earlier versions simply used the number of unique values, but poor distribution of values makes a table less suitable. This was analysed by Fabien Noth and Günter Kniesel.)

  The indexes of dynamic predicates are deleted if the number of clauses is doubled since its creation or reduced below 1/4th. The JIT approach will recreate a suitable index on the next call. Indexes of running predicates cannot be deleted. They are added to a‘removed index list’associated to the predicate. Outdated indexes of predicates are reclaimed by [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0). The clause garbage collector is scheduled automatically, based on time and space based heuristics. See [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0) for details.

The library `library(prolog_jiti)` provides [jiti_list/0](prologjiti.html#jiti_list/0),1 to list the characteristics of all or some of the created hash tables.

**Dynamic predicates** are indexed using the same rules as static predicates, except that the *special purpose* schemes are never applied. In addition, the JITI index is discarded if the number of clauses has doubled since the predicate was last assessed or shrinks below one fourth. A subsequent call reassesses the statistics of the dynamic predicate and, when applicable, creates a new index.

Jit indexing is controlled by a set of Prolog flags whose names starts with `ci_`, e.g., **ci_min_speedup**. See [current_prolog_flag/2](flags.html#current_prolog_flag/2).

### 2.17.1 Deep indexing

As introduced in [section 2.17](jitindex.html#sec:2.17), *deep indexing* creates hash tables distinguish clauses that share a compound with the same name and arity. Deep indexes allow for efficient lookup of arbitrary terms. Without it is advised to *flatten* the term, i.e., turn `F(X)` into two arguments for the fact, one argument denoting the functor `F` and the second the argument X. This works fine as long as the arity of the each of the terms is the same. Alternatively we can use [term_hash/2](db.html#term_hash/2) or [term_hash/4](db.html#term_hash/4) to add a column holding the hash of the term. That approach can deal with arbitrary arities, but requires us to know that the term is ground ([term_hash/2](db.html#term_hash/2)) or up to which depth we get sufficient selectivity ([term_hash/4](db.html#term_hash/4)).

Deep indexing does not require this knowledge and leads to efficient lookup regardless of the instantiation of the query and term. The current version does come with some limitations:

- The decision which index to use is taken independently at each level. Future versions may be smarter on this.
- Deep indexing only applies to a *single argument* indexes (on any argument).
- Currently, the depth of indexing is limited to 7 levels.

Note that, when compiling DCGs (see [section 4.13](DCG.html#sec:4.13)) and the first body term is a *literal*, it is included into the clause head. See for example the grammar and its plain Prolog representation below.

``` code
det(det(a), sg)  --> "a".
det(det(an), pl) --> "an".
det(det(the), _) --> "the".
```

``` code
?- listing(det).
det(det(a), sg, [97|A], A).
det(det(an), pl, [97, 110|A], A).
det(det(the), _, [116, 104, 101|A], A).
```

Deep argument indexing will create indexes for the 3rd list argument, providing speedup and making clause selection deterministic if all rules start with a literal and all literals are unique in the first 6 elements. Note that deep index creation stops as soon as a deterministic choice can be made or there are no two clauses that have the same name/arity combination.

### 2.17.2 Future directions

- The‘special cases’can be extended. This is notably attractive for static predicates with a relatively small number of clauses where a hash lookup is too costly.
- Create an efficient decision diagram for selecting between low numbers of static clauses.
- Implement a better judgements for selecting between deep and plain indexes.

### 2.17.3 Indexing for body code

The current SWI-Prolog versions only consider the head for generating clause indexing. This would make it impossible to examine a head argument and pass the argument in the body without copying the argument. Consider the two clauses below. Both have equal semantics under Prolog. The first version would loose clause indexing while the second creates a copy of the **f/1** argument. Neither is desirable.

``` code
p(X) :- X = f(I), integer(I), q(X).
p(f(I)) :- integer(I), q(f(X)).
```

As of SWI-Prolog 8.3.21, unifications against head arguments that happen before anything else in the body are compiled special. Effectively, the term unified too is moved into the head (providing indexing) and places where this term is used simply use the corresponding argument. The explicit unification is removed. Decompilation ([clause/2](examineprog.html#clause/2)) reverses this process, but may not produce exactly the same term. The re-inserted unifications are ordered according to the argument position and the variable is always on the left hand of the [=/2](compare.html#=/2). Thus,

``` code
p(X,Y) :- f(_) = Y, X = g(_), q(X,Y).
```

Is decompiled into the following equivalent clause.

``` code
p(X,Y) :- X = g(_), Y = f(_), q(X,Y).
```

Additional notes:

- This transformation is only performed on *static* code.
- The unifications must *immediately* follow the head in a *conjunction*.
- As sole exception, calls to [true/0](control.html#true/0) are skipped. This allows [goal_expansion/2](consulting.html#goal_expansion/2) to convert goals to `true` while preserving this optimization.
- If the head argument is not used the body unification is still moved into the head. The decompiler does not inverse the process in that case. Thus, `p(X) :- X = a.` is fully equivalent to `p(a).`
- Currently this optimisation is enabled regardless of the Prolog flag [optimise](flags.html#flag:optimise). As this optimization harms source-level debugging, this may not be desirable. On the other hand we do not want determinism to depend on optimization while this optimization affects determinism.

### 2.17.4 Indexing and portability

The base-line functionality of Prolog implementations provides indexing on constants and functor (name/arity) on the first argument. This must be your assumption if wide portability of your program is important. This can typically be achieved by exploiting [term_hash/2](db.html#term_hash/2) or [term_hash/4](db.html#term_hash/4) and/or maintaining multiple copies of a predicate with reordered arguments and wrappers that update all implementations (assert/retract) and selects the appropriate implementation (query).

YAP provides full JIT indexing, including indexing arguments of compound terms. YAP's indexing has been the inspiration for enhancing SWI-Prolog's indexing capabilities.
