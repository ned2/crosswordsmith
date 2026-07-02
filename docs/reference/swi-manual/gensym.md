
## A.19 library(gensym): Generate unique symbols

Gensym (*Generate Symbols*) is an old library for generating unique symbols (atoms). Such symbols are generated from a base atom which gets a sequence number appended. Of course there is no guarantee that `catch22` is not an already defined atom and therefore one must be aware these atoms are only unique in an isolated context.

The SWI-Prolog gensym library is thread-safe. The sequence numbers are global over all threads and therefore generated atoms are unique over all threads.

**gensym**(`+Base, -Unique`)  
Generate \<`Base`\>1, \<`Base`\>2, etc atoms on each subsequent call. Note that there is nothing that prevents other parts of the application to‘invent’the same identifier. The predicate [gensym/2](gensym.html#gensym/2) is thread-safe in the sense that two threads generating identifiers from the same `Base` will never generate the same identifier.

See also  
uuid/1, [term_hash/2](db.html#term_hash/2), [variant_sha1/2](db.html#variant_sha1/2) may be used to generate various unique or content-based identifiers safely.

**reset_gensym**  
Reset gensym for all registered keys. This predicate is available for compatibility only. New code is strongly advised to avoid the use of reset_gensym or at least to reset only the keys used by your program to avoid unexpected side effects on other components.

**reset_gensym**(`+Base`)  
Restart generation of identifiers from `Base` at \<`Base`\>1. Used to make sure a program produces the same results on subsequent runs. Use with care.
