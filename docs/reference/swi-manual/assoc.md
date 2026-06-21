
## A.4 library(assoc): Association lists

Authors: *Richard A. O'Keefe, L.Damas, V.S.Costa and [Markus Triska](https://www.metalevel.at)*

### A.4.1 Introduction

An *association list* as implemented by this library is a collection of unique *keys* that are associated to *values*. Keys must be ground, values need not be.

An association list can be used to *fetch* elements via their keys and to *enumerate* its elements in ascending order of their keys.

This library uses **AVL trees** to implement association lists. This means that

- inserting a key
- changing an association
- fetching a single element

are all *O(`log(N)`)* *worst-case* (and expected) time operations, where *N* denotes the number of elements in the association list.

The logarithmic overhead is often acceptable in practice. Notable advantages of association lists over several other methods are:

- `library(assoc)` is written entirely in Prolog, making it portable to other systems
- the interface predicates fit the declarative nature of Prolog, avoiding destructive updates to terms
- AVL trees scale very predictably and can be used to represent sparse arrays efficiently.

### A.4.2 Creating association lists

An association list is *created* with one of the following predicates:

\[semidet\]**empty_assoc**(`?Assoc`)  
Is true if `Assoc` is the empty association list.

\[det\]**list_to_assoc**(`+Pairs, -Assoc`)  
Create an association from a list `Pairs` of Key-Value pairs. List must not contain duplicate keys.

Errors  
`domain_error(unique_key_pairs, List)` if List contains duplicate keys

\[det\]**ord_list_to_assoc**(`+Pairs, -Assoc`)  
`Assoc` is created from an ordered list `Pairs` of Key-Value pairs. The pairs must occur in strictly ascending order of their keys.

Errors  
`domain_error(key_ordered_pairs, List)` if pairs are not ordered.

### A.4.3 Querying association lists

An association list can be *queried* with:

\[semidet\]**get_assoc**(`+Key, +Assoc, -Value`)  
True if `Key`-`Value` is an association in `Assoc`.

\[semidet\]**get_assoc**(`+Key, +Assoc0, ?Val0, ?Assoc, ?Val`)  
True if `Key`-`Val0` is in `Assoc0` and `Key`-`Val` is in `Assoc`.

\[semidet\]**max_assoc**(`+Assoc, -Key, -Value`)  
True if `Key`-`Value` is in `Assoc` and `Key` is the largest key.

\[semidet\]**min_assoc**(`+Assoc, -Key, -Value`)  
True if `Key`-`Value` is in assoc and `Key` is the smallest key.

\[nondet\]**gen_assoc**(`?Key, +Assoc, ?Value`)  
True if `Key`-`Value` is an association in `Assoc`. Enumerates keys in ascending order on backtracking.

See also  
[get_assoc/3](assoc.html#get_assoc/3).

### A.4.4 Modifying association lists

Elements of an association list can be changed and inserted with:

\[det\]**put_assoc**(`+Key, +Assoc0, +Value, -Assoc`)  
`Assoc` is `Assoc0`, except that `Key` is associated with `Value`. This can be used to insert and change associations.

\[semidet\]**del_assoc**(`+Key, +Assoc0, ?Value, -Assoc`)  
True if `Key`-`Value` is in `Assoc0`. `Assoc` is `Assoc0` with `Key`-`Value` removed.

\[semidet\]**del_min_assoc**(`+Assoc0, ?Key, ?Val, -Assoc`)  
True if `Key`-Value is in `Assoc0` and `Key` is the smallest key. `Assoc` is `Assoc0` with `Key`-Value removed. Warning: This will succeed with *no* bindings for `Key` or `Val` if `Assoc0` is empty.

\[semidet\]**del_max_assoc**(`+Assoc0, ?Key, ?Val, -Assoc`)  
True if `Key`-Value is in `Assoc0` and `Key` is the greatest key. `Assoc` is `Assoc0` with `Key`-Value removed. Warning: This will succeed with *no* bindings for `Key` or `Val` if `Assoc0` is empty.

### A.4.5 Conversion predicates

Conversion of (parts of) an association list to *lists* is possible with:

\[det\]**assoc_to_list**(`+Assoc, -Pairs`)  
Translate `Assoc` to a list `Pairs` of Key-Value pairs. The keys in `Pairs` are sorted in ascending order.

\[det\]**assoc_to_keys**(`+Assoc, -Keys`)  
True if `Keys` is the list of keys in `Assoc`. The keys are sorted in ascending order.

\[det\]**assoc_to_values**(`+Assoc, -Values`)  
True if `Values` is the list of values in `Assoc`. `Values` are ordered in ascending order of the key to which they were associated. `Values` may contain duplicates.

### A.4.6 Reasoning about association lists and their elements

Further inspection predicates of an association list and its elements are:

\[semidet\]**is_assoc**(`+Assoc`)  
True if `Assoc` is an association list. This predicate checks that the structure is valid, elements are in order, and tree is balanced to the extent guaranteed by AVL trees. I.e., branches of each subtree differ in depth by at most 1. Does *not* validate that keys are sufficiently instantiated to ensure the tree remains valid if a key is further instantiated.

\[semidet\]**map_assoc**(`:Pred, +Assoc`)  
True if `Pred`(Value) is true for all values in `Assoc`.

\[semidet\]**map_assoc**(`:Pred, +Assoc0, ?Assoc`)  
Map corresponding values. True if `Assoc` is `Assoc0` with `Pred` applied to all corresponding pairs of of values.
