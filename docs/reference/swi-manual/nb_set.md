
## A.28 library(nb_set): Non-backtrackable set

The library `library(nb_set)` defines *non-backtrackable sets*, implemented as binary trees. The sets are represented as compound terms and manipulated using [nb_setarg/3](manipterm.html#nb_setarg/3). Non-backtrackable manipulation of data structures is not supported by a large number of Prolog implementations, but it has several advantages over using the database. It produces less garbage, is thread-safe, reentrant and deals with exceptions without leaking data.

Similar to the `library(assoc)` library, keys can be any Prolog term, but it is not allowed to instantiate or modify a term.

One of the ways to use this library is to generate unique values on backtracking *without* generating *all* solutions first, for example to act as a filter between a generator producing many duplicates and an expensive test routine, as outlined below:

``` code
generate_and_test(Solution) :-
        empty_nb_set(Set),
        generate(Solution),
        add_nb_set(Solution, Set, true),
        test(Solution).
```

**empty_nb_set**(`?Set`)  
True if `Set` is a non-backtrackable empty set.

**add_nb_set**(`+Key, !Set`)  
Add `Key` to `Set`. If `Key` is already a member of `Set`, [add_nb_set/3](nb_set.html#add_nb_set/3) succeeds without modifying `Set`.

**add_nb_set**(`+Key, !Set, ?New`)  
If `Key` is not in `Set` and `New` is unified to `true`, `Key` is added to `Set`. If `Key` is in `Set`, `New` is unified to `false`. It can be used for many purposes:

|                           |                            |
|---------------------------|----------------------------|
| `add_nb_set(+, +, false)` | Test membership            |
| `add_nb_set(+, +, true)`  | Succeed only if new member |
| `add_nb_set(+, +, Var)`   | Succeed, binding `Var`     |

**gen_nb_set**(`+Set, -Key`)  
Generate all members of `Set` on backtracking in the standard order of terms. To test membership, use [add_nb_set/3](nb_set.html#add_nb_set/3).

**size_nb_set**(`+Set, -Size`)  
Unify `Size` with the number of elements in `Set`.

**nb_set_to_list**(`+Set, -List`)  
Unify `List` with a list of all elements in `Set` in the standard order of terms (i.e., an *ordered list*).
