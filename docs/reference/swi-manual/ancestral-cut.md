
## B.2 Ancestral cuts

**prolog_cut_to**(`+Choice`)  
Prunes all choice points created since `Choice`. Can be used together with [prolog_current_choice/1](manipstack.html#prolog_current_choice/1) to implement *ancestral* cuts. This predicate is in the hackers corner because it should not be used in normal Prolog code. It may be used to create new high level control structures, particularly for compatibility purposes.

Note that in the current implementation, the pruned choice points and environment frames are *not* reclaimed. As a consequence, where predicates that are deterministic due to clause indexing, normal cuts or `(if->then;else)` and and tail recursive run in bounded local stack space, predicates using [prolog_cut_to/1](ancestral-cut.html#prolog_cut_to/1) will run out of stack.
