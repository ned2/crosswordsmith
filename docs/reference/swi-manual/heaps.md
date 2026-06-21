
## A.20 library(heaps): heaps/priority queues

author  
Lars Buitinck

Heaps are data structures that return the entries inserted into them in an ordered fashion, based on a priority. This makes them the data structure of choice for implementing priority queues, a central element of algorithms such as best-first or A\* search and Kruskal's minimum-spanning-tree algorithm.

This module implements min-heaps, meaning that key-value items are retrieved in ascending order of key. In other words, the key determines the priority. It was designed to be compatible with the SICStus Prolog library module of the same name. [merge_heaps/3](heaps.html#merge_heaps/3) and [singleton_heap/3](heaps.html#singleton_heap/3) are SWI-specific extension. The portray_heap/1 predicate is not implemented.

Although the values can be arbitrary Prolog terms, the keys determine the priority, so keys must be ordered by [@=\</2](compare.html#@=%3C/2). This means that variables can be used as keys, but binding them in between heap operations may change the ordering. It also means that rational terms (cyclic trees), for which standard order is not well-defined, cannot be used as keys.

The current version implements pairing heaps. These support insertion and merging both in constant time, deletion of the minimum in logarithmic amortized time (though delete-min, i.e., get_from_heap/3, takes linear time in the worst case).

\[semidet\]**add_to_heap**(`+Heap0, +Key, ?Value, -Heap`)  
Adds `Value` with priority `Key` to `Heap0`, constructing a new heap in `Heap`.

\[semidet\]**delete_from_heap**(`+Heap0, -Key, +Value, -Heap`)  
Deletes `Value` from `Heap0`, leaving its priority in `Key` and the resulting data structure in `Heap`. Fails if `Value` is not found in `Heap0`.

bug  
This predicate is extremely inefficient and exists only for SICStus compatibility.

\[semidet\]**empty_heap**(`?Heap`)  
True if `Heap` is an empty heap. Complexity: constant.

\[semidet\]**singleton_heap**(`?Heap, ?Key, ?Value`)  
True if `Heap` is a heap with the single element `Key`-`Value`.

Complexity: constant.

\[semidet\]**get_from_heap**(`?Heap0, ?Key, ?Value, -Heap`)  
Retrieves the minimum-priority pair `Key`-`Value` from `Heap0`. `Heap` is `Heap0` with that pair removed. Complexity: logarithmic (amortized), linear in the worst case.

\[det\]**heap_size**(`+Heap, -Size:int`)  
Determines the number of elements in `Heap`. Complexity: constant.

\[det\]**heap_to_list**(`+Heap, -List:list`)  
Constructs a list `List` of Key-Value terms, ordered by (ascending) priority. Complexity: O(N log N).

\[semidet\]**is_heap**(`+X`)  
Returns true if `X` is a heap. Validates the consistency of the entire heap. Complexity: linear.

\[det\]**list_to_heap**(`+List:list, -Heap`)  
If `List` is a list of Key-Value terms, constructs a heap out of `List`. Complexity: linear.

\[semidet\]**min_of_heap**(`+Heap, ?Key, ?Value`)  
Unifies `Value` with the minimum-priority element of `Heap` and `Key` with its priority value. Complexity: constant.

\[semidet\]**min_of_heap**(`+Heap, ?Key1, ?Value1, ?Key2, ?Value2`)  
Gets the two minimum-priority elements from `Heap`. Complexity: logarithmic (amortized).

bug  
This predicate is extremely inefficient and exists for compatibility with earlier implementations of this library and SICStus compatibility. It performs a linear amount of work in the worst case that a following get_from_heap has to re-do.

\[det\]**merge_heaps**(`+Heap0, +Heap1, -Heap`)  
Merge the two heaps `Heap0` and `Heap1` in `Heap`. Complexity: constant.
