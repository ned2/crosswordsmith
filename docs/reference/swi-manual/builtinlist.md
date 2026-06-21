
## 4.29 Built-in list operations

Most list operations are defined in the library `library(lists)` described in [section A.25](lists.html#sec:A.25). Some that are implemented with more low-level primitives are built-in and described here.

**is_list**(`+Term`)  
True if `Term` is bound to the empty list (`[]`) or a compound term with name‘`[|]`’^(136The traditional list functor name is the dot (`’.’`). This is still the case of the command line option **--traditional** is given. See also [section 5.1](ext-lists.html#sec:5.1).) and arity 2 and the second argument is a list.^(137In versions before 5.0.1, [is_list/1](builtinlist.html#is_list/1) just checked for `[]` or `[_|_]` and proper_list/1 had the role of the current [is_list/1](builtinlist.html#is_list/1). The current definition conforms to the de facto standard. Assuming proper coding standards, there should only be very few cases where a quick-and-dirty [is_list/1](builtinlist.html#is_list/1) is a good choice. Richard O'Keefe pointed at this issue.) This predicate acts as if defined by the definition below on *acyclic* terms. The implementation safely *fails* if `Term` represents a cyclic list.

``` code
is_list(X) :-
        var(X), !,
        fail.
is_list([]).
is_list([_|T]) :-
        is_list(T).
```

\[semidet\]**memberchk**(`?Elem, +List`)  
True when `Elem` is an element of `List`. This‘chk’variant of [member/2](lists.html#member/2) is semi deterministic and typically used to test membership of a list. Raises a `type` error if scanning `List` encounters a non-list. Note that [memberchk/2](builtinlist.html#memberchk/2) does *not* perform a full list typecheck. For example, `memberchk(a, [a|b])` succeeds without error. If `List` is cyclic and `Elem` is not a member of `List`, [memberchk/2](builtinlist.html#memberchk/2) eventually raises a `type` error.^(138*Eventually* here means it will scan as many elements as the longest list that may exist given the current stack usage before raising the exception.)

\[ISO\]**length**(`?List, ?Length`)  
True if `Length` represents the number of elements in `List`. This predicate is a true relation and can be used to find the length of a list or produce a list (holding variables) of length `Length`. The predicate is non-deterministic, producing lists of increasing length if `List` is a *partial list* and `Length` is a variable.

``` code
?- length(List,4).
List = [_27940,_27946,_27952,_27958].

?- length(List,Length).
List = [], Length = 0 ;
List = [_24698], Length = 1 ;
List = [_24698,_25826], Length = 2
...
```

It raises errors if `Length` is bound to a non-integer or a negative integer or if `List` is neither a list nor a partial list. This error condition includes cyclic lists:^(139ISO demands failure here. We think an error is more appropriate.)

``` code
?- A=[1,2,3|A], length(A,L).
ERROR: Type error: `list' expected ...
```

Covering an edge case, the predicate fails if the tail of `List` is equivalent to `Length`:^(140This is logically correct. An exception would be more appropriate, but to our best knowledge, current practice in Prolog does not describe a suitable candidate exception term.)

``` code
?- List=[1,2,3|Length],length(List,Length).
false.

?- length(Length,Length).
false.
```

\[ISO\]**sort**(`+List, -Sorted`)  
True if `Sorted` can be unified with a list holding the elements of `List`, sorted to the standard order of terms (see [section 4.6](compare.html#sec:4.6)). Duplicates are removed. The implementation is in C, using *natural merge sort*.^(141Contributed by Richard O'Keefe.) The [sort/2](builtinlist.html#sort/2) predicate can sort a cyclic list, returning a non-cyclic version with the same elements.

Note that `List` may contain non-ground terms. If `Sorted` is unbound at call-time, for each consecutive pair of elements in `Sorted`, the relation `E1 @< E2` will hold. However, unifying a variable in `Sorted` may cause this relation to become invalid, *even* unifying a variable in `Sorted` with another (older) variable. See also [section 4.6.1](compare.html#sec:4.6.1).

**sort**(`+Key, +Order, +List, -Sorted`)  
True when `Sorted` can be unified with a list holding the element of `List`. `Key` determines which part of each element in `List` is used for comparing two term and `Order` describes the relation between each set of consecutive elements in `Sorted`.^(142The definition of this predicate was established after discussion with Joachim Schimpf from the ECLiPSe team. ECLiPSe currently only accepts `<`, `=<`, `>` and `>=` for the `Order` argument but this is likely to change. SWI-Prolog extends this predicate to deal with dicts.)

If `Key` is the integer zero (0), the entire term is used to compare two elements. Using `Key`=0 can be used to sort arbitrary Prolog terms. Other values for `Key` can only be used with compound terms or dicts (see [section 5.4](bidicts.html#sec:5.4)). An integer key extracts the `Key`-th argument from a compound term. An integer or atom key extracts the value from a dict that is associated with the given key. A type_error is raised if the list element is of the wrong type and an existence_error is raised if the compound has not enough argument or the dict does not contain the requested key.

Deeper nested elements of structures can be selected by using a list of keys for the `Key` argument.

The `Order` argument is described in the table below:^(143For compatibility with ECLiPSe, the values `<`, `=<`, `>` and `>=` are allowed as synonyms.)

|       |            |                    |
|-------|------------|--------------------|
| Order | Ordering   | Duplicate handling |
| `@<`  | ascending  | remove             |
| `@=<` | ascending  | keep               |
| `@>`  | descending | remove             |
| `@>=` | descending | keep               |

The sort is *stable*, which implies that, if duplicates are kept, the order of duplicates is not changed. If duplicates are removed, only the first element of a sequence of duplicates appears in `Sorted`.

This predicate supersedes most of the other sorting primitives, for example:

``` code
sort(List, Sorted)     :- sort(0,  @<, List,  Sorted).
msort(List, Sorted)    :- sort(0, @=<, List,  Sorted).
keysort(Pairs, Sorted) :- sort(1, @=<, Pairs, Sorted).
```

The following example sorts a list of rows, for example resulting from [csv_read_file/2](csv.html#csv_read_file/2)) ascending on the 3rd column and descending on the 4th column (for sets of rows where the 3rd column is equal):

``` code
    sort(4, @>=, Rows0, Rows1),
    sort(3, @=<, Rows1, Sorted).
```

See also [sort/2](builtinlist.html#sort/2) (ISO), [msort/2](builtinlist.html#msort/2), [keysort/2](builtinlist.html#keysort/2), [predsort/3](builtinlist.html#predsort/3) and [order_by/2](solutionsequences.html#order_by/2).

**msort**(`+List, -Sorted`)  
Equivalent to [sort/2](builtinlist.html#sort/2), but does not remove duplicates. Raises a `type_error` if `List` is a cyclic list or not a list.

\[ISO\]**keysort**(`+List, -Sorted`)  
Sort a list of *pairs*. `List` must be a list of `Key``-``Value` pairs, terms whose principal functor is (-)/2. `List` is sorted on `Key` according to the standard order of terms (see [section 4.6.1](compare.html#sec:4.6.1)). Duplicates are *not* removed. Sorting is *stable* with regard to the order of the `Values`, i.e., the order of multiple elements that have the same `Key` is not changed.

The [keysort/2](builtinlist.html#keysort/2) predicate is often used together with library `library(pairs)`. It can be used to sort lists on different or multiple criteria. For example, the following predicates sorts a list of atoms according to their length, maintaining the initial order for atoms that have the same length.

``` code
:- use_module(library(pairs)).

sort_atoms_by_length(Atoms, ByLength) :-
        map_list_to_pairs(atom_length, Atoms, Pairs),
        keysort(Pairs, Sorted),
        pairs_values(Sorted, ByLength).
```

**predsort**(`+Pred, +List, -Sorted`)  
Sorts similar to [sort/2](builtinlist.html#sort/2), but determines the order of two terms by calling `Pred`(-`Delta`, +`E1`, +`E2`) . This call must unify `Delta` with one of `<`, `>` or `=`. Duplicates are removed (i.e. equivalence classes of elements as defined by `Pred` are collapsed to a single element in `Sorted`) If the built-in predicate [compare/3](compare.html#compare/3) is used, the result is the same as [sort/2](builtinlist.html#sort/2). See also [keysort/2](builtinlist.html#keysort/2).
