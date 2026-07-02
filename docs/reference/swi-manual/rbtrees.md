
## A.48 library(rbtrees): Red black trees

author  
Vitor Santos Costa, Jan Wielemaker, Samer Abdallah, Peter Ludemann.

See also  
\- `library(pairs)`, `library(assoc)`  
- "Introduction to Algorithms", Second Edition Cormen, Leiserson, Rivest, and Stein, MIT Press

Red-Black trees are balanced search binary trees. They are named because nodes can be classified as either red or black. The code we include is based on "Introduction to Algorithms", second edition, by Cormen, Leiserson, Rivest and Stein. The library includes routines to insert, lookup and delete elements in the tree.

A Red black tree is represented as a term `t(Nil, Tree)`, where Nil is the Nil-node, a node shared for each nil-node in the tree. Any node has the form `colour(Left, Key, Value, Right)`, where *colour* is one of `red` or `black`.

**Warning: instantiation of keys**

Red-Black trees depend on the Prolog *standard order of terms* to organize the keys as a (balanced) binary tree. This implies that any term may be used as a key. The tree may produce wrong results, such as not being able to find a key, if the ordering of keys changes after the key has been inserted into the tree. The user is responsible to ensure that variables used as keys or appearing in a term used as key that may affect ordering are not unified, with the exception of unification against new fresh variables. For this reason, *ground* terms are safe keys. When using non-ground terms, either make sure the variables appear in places that do not affect the standard order relative to other keys in the tree or make sure to not unify against these variables as long as the tree is being used.

\[det\]**rb_new**(`-Tree`)  
Create a new Red-Black tree `Tree`.

deprecated  
Use [rb_empty/1](rbtrees.html#rb_empty/1).

\[semidet\]**rb_empty**(`?Tree`)  
Succeeds if `Tree` is an empty Red-Black tree.

\[semidet\]**rb_lookup**(`+Key, -Value, +Tree`)  
True when `Value` is associated with `Key` in the Red-Black tree `Tree`. The given `Key` may include variables, in which case the RB tree is searched for a key with equivalent variables (using (`==`)/2). Time complexity is O(log N) in the number of elements in the tree.

See also  
[rb_in/3](rbtrees.html#rb_in/3) for backtracking over keys.

\[semidet\]**rb_min**(`+Tree, -Key, -Value`)  
`Key` is the minimum key in `Tree`, and is associated with Val.

\[semidet\]**rb_max**(`+Tree, -Key, -Value`)  
`Key` is the maximal key in `Tree`, and is associated with Val.

\[semidet\]**rb_next**(`+Tree, +Key, -Next, -Value`)  
`Next` is the next element after `Key` in `Tree`, and is associated with Val. Fails if `Key` isn't in `Tree` or if `Key` is the maximum key.

\[semidet\]**rb_previous**(`+Tree, +Key, -Previous, -Value`)  
`Previous` is the previous element after `Key` in `Tree`, and is associated with Val. Fails if `Key` isn't in `Tree` or if `Key` is the minimum key.

\[semidet\]**rb_update**(`+Tree, +Key, ?NewVal, -NewTree`)  
`Tree` `NewTree` is tree `Tree`, but with value for `Key` associated with `NewVal`. Fails if `Key` is not in `Tree` (using (`==`)/2). This predicate may fail or give unexpected results if `Key` is not sufficiently instantiated.

See also  
[rb_in/3](rbtrees.html#rb_in/3) for backtracking over keys.

\[semidet\]**rb_update**(`+Tree, +Key, -OldVal, ?NewVal, -NewTree`)  
Same as `rb_update(Tree, Key, NewVal, NewTree)` but also unifies `OldVal` with the value associated with `Key` in `Tree`.

\[semidet\]**rb_apply**(`+Tree, +Key, :G, -NewTree`)  
If the value associated with key `Key` is Val0 in `Tree`, and if `call(G,Val0,ValF)` holds, then `NewTree` differs from `Tree` only in that `Key` is associated with value ValF in tree `NewTree`. Fails if it cannot find `Key` in `Tree`, or if `call(G,Val0,ValF)` is not satisfiable.

\[nondet\]**rb_in**(`?Key, ?Value, +Tree`)  
True when `Key`-`Value` is a key-value pair in red-black tree `Tree`. Same as below, but does not materialize the pairs.

``` code
rb_visit(Tree, Pairs), member(Key-Value, Pairs)
```

Leaves a choicepoint even if `Key` is instantiated; to avoid a choicepoint, use [rb_lookup/3](rbtrees.html#rb_lookup/3).

\[det\]**rb_insert**(`+Tree, +Key, ?Value, -NewTree`)  
Add an element with key `Key` and `Value` to the tree `Tree` creating a new red-black tree `NewTree`. If `Key` is a key in `Tree`, the associated value is replaced by `Value`. See also [rb_insert_new/4](rbtrees.html#rb_insert_new/4). Does *not* validate that `Key` is sufficiently instantiated to ensure the tree remains valid if a key is further instantiated.

\[semidet\]**rb_insert_new**(`+Tree, +Key, ?Value, -NewTree`)  
Add a new element with key `Key` and `Value` to the tree `Tree` creating a new red-black tree `NewTree`. Fails if `Key` is a key in `Tree`. Does *not* validate that `Key` is sufficiently instantiated to ensure the tree remains valid if a key is further instantiated.

**rb_delete**(`+Tree, +Key, -NewTree`)  
Delete element with key `Key` from the tree `Tree`, returning the value Val associated with the key and a new tree `NewTree`. Fails if `Key` is not in `Tree` (using (`==`)/2).

See also  
[rb_in/3](rbtrees.html#rb_in/3) for backtracking over keys.

**rb_delete**(`+Tree, +Key, -Val, -NewTree`)  
Same as `rb_delete(Tree, Key, NewTree)`, but also unifies `Val` with the value associated with `Key` in `Tree`.

**rb_del_min**(`+Tree, -Key, -Val, -NewTree`)  
Delete the least element from the tree `Tree`, returning the key `Key`, the value `Val` associated with the key and a new tree `NewTree`. Fails if `Tree` is empty.

**rb_del_max**(`+Tree, -Key, -Val, -NewTree`)  
Delete the largest element from the tree `Tree`, returning the key `Key`, the value `Val` associated with the key and a new tree `NewTree`. Fails if `Tree` is empty.

\[det\]**rb_visit**(`+Tree, -Pairs`)  
`Pairs` is an infix visit of tree `Tree`, where each element of `Pairs` is of the form Key-Value.

\[det\]**rb_visit_range**(`+Tree, +Min, +Max, -Pairs`)  
Retrieves a range of pairs with keys between `Min` and `Max` (inclusive) from a `Tree` using standard term comparison.

\[semidet\]**rb_map**(`+Tree, :G, -NewTree`)  
For all nodes Key in the tree `Tree`, if the value associated with key Key is Val0 in tree `Tree`, and if `call(G,Val0,ValF)` holds, then the value associated with Key in `NewTree` is ValF. Fails if `call(G,Val0,ValF)` is not satisfiable for all Val0. If `G` is non-deterministic, [rb_map/3](rbtrees.html#rb_map/3) will backtrack over all possible values from `call(G,Val0,ValF)`. You should not depend on the order of tree traversal (currently: key order).

\[semidet\]**rb_map**(`+T, :Goal`)  
True if `call(Goal, Value)` is true for all nodes in `T`.

**rb_fold**(`:Goal, +Tree, +State0, -State`)  
Fold the given predicate over all the key-value pairs in `Tree`, starting with initial state `State0` and returning the final state `State`. Pred is called as

``` code
call(Pred, Key-Value, State1, State2)
```

Determinism depends on `Goal`.

\[det\]**rb_clone**(`+TreeIn, -TreeOut, -Pairs`)  
‘Clone’the red-back tree `TreeIn` into a new tree `TreeOut` with the same keys as the original but with all values set to unbound values. `Pairs` is a list containing all new nodes as pairs K-V.

**rb_partial_map**(`+Tree, +Keys, :G, -NewTree`)  
For all nodes Key in `Keys`, if the value associated with key Key is Val0 in tree `Tree`, and if `call(G,Val0,ValF)` holds, then the value associated with Key in `NewTree` is ValF, otherwise it is the value associated with the key in `Tree`. Fails if Key isn't in `Tree` or if `call(G,Val0,ValF)` is not satisfiable for all Val0 in `Keys`. Assumes keys are sorted and not repeated (fails if this is not true).

\[det\]**rb_keys**(`+Tree, -Keys`)  
`Keys` is unified with an ordered list of all keys in the Red-Black tree `Tree`.

\[det\]**list_to_rbtree**(`+List, -Tree`)  
`Tree` is the red-black tree corresponding to the mapping in `List`, which should be a list of Key-Value pairs. `List` should not contain more than one entry for each distinct key, but this is not validated by [list_to_rbtree/2](rbtrees.html#list_to_rbtree/2).

\[det\]**ord_list_to_rbtree**(`+List, -Tree`)  
`Tree` is the red-black tree corresponding to the mapping in list `List`, which should be a list of Key-Value pairs. `List` should not contain more than one entry for each distinct key, but this is not validated by [ord_list_to_rbtree/2](rbtrees.html#ord_list_to_rbtree/2). `List` is assumed to be sorted according to the standard order of terms.

\[det\]**rb_size**(`+Tree, -Size`)  
`Size` is the number of elements in `Tree`.

\[semidet\]**is_rbtree**(`@Term`)  
True if `Term` is a valid Red-Black tree. Processes the entire tree, checking the coloring of the nodes, the balance and the ordering of keys. Does *not* validate that keys are sufficiently instantiated to ensure the tree remains valid if a key is further instantiated.
