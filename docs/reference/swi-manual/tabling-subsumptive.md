
## 7.5 Variant and subsumptive tabling

By default, SWI-Prolog (and other Prolog systems with tabling) create a table per call *variant*. A call (term) is a variant of another call (term) if there is a renaming of variables that makes the two terms equal. See [=@=/2](compare.html#=@=/2) for details. Consider the following program:

``` code
:- table p/1.

p(X) :- p(Y), Y < 10 000, X is Y+1.
p(1).
```

Calling `p(X)` creates a table for this variant with 10,000 answers. Calling `p(42)` creates a new table where the recursive call (`p(Y)`) uses the previously created table to enumerate all values `1 ... 41` before deriving `p(42)` is true. *Early completion* (see [section 7.4](tnotpure.html#sec:7.4)) in this case prevents enumerating all answers for `p(Y)` (`1 ... 10,000`). As a result, the query below runs in quadratic time and creates a 10,000 additional tables.

``` code
?- time(forall(between(1, 10 000, X), p(X))).
% 150,370,553 inferences, 13.256 CPU in 13.256 seconds
```

*Subsumptive* tabling answers a query using answers from a more general table. In this case, this means it uses basically [trie_gen/2](db.html#trie_gen/2) to get the answer `p(42)` from the table `p(_)`. This leads to the program and results shown below.

``` code
:- table p/1 as subsumptive.

p(X) :- p(Y), Y < 10 000, X is Y+1.
p(1).
```

``` code
?- time(p(_)).
% 140,066 inferences, 0.015 CPU in 0.015 seconds
?- time(t).
% 170,005 inferences, 0.016 CPU in 0.016 seconds
```

*Subsumptive* tabling can be activated in two ways. Per table assign the `... as subsumptive` option and globally by setting the [table_subsumptive](flags.html#flag:table_subsumptive) flag to `true`.

One may wonder why subsumptive tabling is not the default. There are also some drawbacks:

- Subsumptive tabling only provides correct support if instances (more specific) queries indeed provides answers that are consistent with the more general query. This is true for *pure programs*, but not guaranteed for arbitrary Prolog programs.
- Finding more generic tables is more complicated and expensive than finding the call variant table and extracting the subset of answers that match the more specific query can be expensive.
- Using subsumptive tables can create more dependencies in the call graph which can slow down the table completion process. Larger dependent components also negatively impact the issues discussed in [section 7.4](tnotpure.html#sec:7.4).
