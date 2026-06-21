
## 4.16 Examining the program

**current_atom**(`-Atom`)  
Successively unifies `Atom` with all atoms known to the system. Note that [current_atom/1](examineprog.html#current_atom/1) always succeeds if `Atom` is instantiated to an atom.

**current_blob**(`?Blob, ?Type`)  
Examine the type or enumerate blobs of the given `Type`. Typed blobs are supported through the foreign language interface for storing arbitrary BLOBs (Binary Large Object) or handles to external entities. See [section 12.4.10](foreigninclude.html#sec:12.4.10) for details.

**current_functor**(`?Name, ?Arity`)  
True when `Name`/`Arity` is a known functor. This means that at some point in time a term with name `Name` and `Arity` arguments was created. Functor objects are currently not subject to garbage collection. Due to timing, t/2 below with instantiated `Name` and `Arity` can theoretically fail, i.e., a functor may be visible in instantiated mode while it is not yet visible in unbound mode. Considering that the only practical value of [current_functor/2](examineprog.html#current_functor/2) we are aware of is to analyse resource usage we accept this impure behaviour.

``` code
t(Name, Arity) :-
    (   current_functor(Name, Arity)
    ->  current_functor(N, A), N == Name, A == Arity
    ;   true
    ).
```

**current_flag**(`-FlagKey`)  
Successively unifies `FlagKey` with all keys used for flags (see [flag/3](db.html#flag/3)).

**current_key**(`-Key`)  
Successively unifies `Key` with all keys used for records (see [recorda/3](db.html#recorda/3), etc.).

\[ISO\]**current_predicate**(`:PredicateIndicator`)  
True if `PredicateIndicator` is a currently defined predicate. A predicate is considered defined if it exists in the specified module, is imported into the module or is defined in one of the modules from which the predicate will be imported if it is called (see [section 6.10](importmodule.html#sec:6.10)). Note that [current_predicate/1](examineprog.html#current_predicate/1) does *not* succeed for predicates that can be *autoloaded* unless they are imported using [autoload/2](module-autoload.html#autoload/2). See also [current_predicate/2](examineprog.html#current_predicate/2) and [predicate_property/2](examineprog.html#predicate_property/2).

If `PredicateIndicator` is not fully specified, the predicate only generates values that are defined in or already imported into the target module. Generating all callable predicates therefore requires enumerating modules using [current_module/1](manipmodule.html#current_module/1). Generating predicates callable in a given module requires enumerating the import modules using [import_module/2](importmodule.html#import_module/2) and the autoloadable predicates using the [predicate_property/2](examineprog.html#predicate_property/2) `autoload`.

**current_predicate**(`?Name, :Head`)  
Classical pre-ISO implementation of [current_predicate/1](examineprog.html#current_predicate/1), where the predicate is represented by the head term. The advantage is that this can be used for checking the existence of a predicate before calling it without the need for [functor/3](manipterm.html#functor/3):

``` code
call_if_exists(G) :-
        current_predicate(_, G),
        call(G).
```

Because of this intended usage, [current_predicate/2](examineprog.html#current_predicate/2) also succeeds if the predicate can be autoloaded. Unfortunately, checking the autoloader makes this predicate relatively slow, in particular because a failed lookup of the autoloader will cause the autoloader to verify that its index is up-to-date.

**predicate_property**(`:Head, ?Property`)  
True when `Head` refers to a predicate that has property `Property`. With sufficiently instantiated `Head`, [predicate_property/2](examineprog.html#predicate_property/2) tries to resolve the predicate the same way as calling it would do: if the predicate is not defined it scans the default modules (see [default_module/2](importmodule.html#default_module/2)) and finally tries the autoloader. Unlike calling, failure to find the target predicate causes [predicate_property/2](examineprog.html#predicate_property/2) to fail silently. If `Head` is not sufficiently bound, only currently locally defined and already imported predicates are enumerated. See [current_predicate/1](examineprog.html#current_predicate/1) for enumerating all predicates. A common issue concerns *generating* all built-in predicates. This can be achieved using the code below:

``` code
generate_built_in(Name/Arity) :-
    predicate_property(system:Head, built_in),
    functor(Head, Name, Arity),
    \+ sub_atom(Name, 0, _, _, $).   % discard reserved names
```

The predicate [predicate_property/2](examineprog.html#predicate_property/2) is covered by part-II of the ISO standard (modules). Although we are not aware of any Prolog system that implements part-II of the ISO standard, [predicate_property/2](examineprog.html#predicate_property/2) is available in most systems. There is little consensus on the implemented properties though. SWI-Prolog's *auto loading* feature further complicate this predicate.

`Property` is one of:

**autoload**(`File`)  
True if the predicate can be autoloaded from the file `File`. Like `undefined`, this property is *not* generated.

**built_in**  
True if the predicate is locked as a built-in predicate. This implies it cannot be redefined in its definition module and it can normally not be seen in the tracer.

**defined**  
True if the predicate is defined. This property is aware of sources being *reloaded*, in which case it claims the predicate defined only if it is defined in another source or it has seen a definition in the current source. See [compile_aux_clauses/1](consulting.html#compile_aux_clauses/1).

**det**  
The predicate is defined to be deterministic using [det/1](debug-determinism.html#det/1).

**discontiguous**  
True after [discontiguous/1](dynamic.html#discontiguous/1) was used to flag that the clauses of the predicates may not be contiguous.

**dynamic**  
True if [assert/1](db.html#assert/1) and [retract/1](db.html#retract/1) may be used to modify the predicate. This property is set using [dynamic/1](dynamic.html#dynamic/1).

**exported**  
True if the predicate is in the public list of the context module.

**imported_from**(`Module`)  
Is true if the predicate is imported into the context module from module `Module`.

**file**(`FileName`)  
Unify `FileName` with the name of the source file in which the predicate is defined. See also [source_file/2](consulting.html#source_file/2) and the property `line_count`. Note that this reports the file of the first clause of a predicate. A more robust interface can be achieved using [nth_clause/3](examineprog.html#nth_clause/3) and [clause_property/2](examineprog.html#clause_property/2).

**foreign**  
True if the predicate is defined in the C language.

**implementation_module**(`-Module`)  
True when `Module` is the module in which `Head` is or will be defined. Resolving this property goes through the same search mechanism as when an undefined predicate is encountered, but does not perform any loading. It searches (1) the module inheritance hierarchy (see [default_module/2](importmodule.html#default_module/2)) and (2) the autoload index if the [unknown](flags.html#flag:unknown) flag is not set to `fail` in the target module.

**indexed**(`Indexes`)  
`Indexes` is a list of additional (hash) indexes on the predicate. Each element of the list is a dict holding the following keys:

**arguments**: `List`  
1-based argument numbers of the predicate or compound used for the index.

**position**: `List`  
1-based argument numbers to find the compound for a *deep* index. Empty list (`[]`) for predicate arguments. A list `[1,2]` implies we index on the 2nd argument of the 1st argument of the predicate, i.e., on the `I` in the head `p(f(_,I))`

**buckets**: `Count`  
Number of buckets of the hash

**list**: `Boolean`  
The hash is a *deep* index, i.e., it points at compounds and there may be sub-indexes on these compounds.

**realised**: `Boolean`  
If `false`, the index is associated with the predicate but not filled. It will be filled if the predicate is called such that this index is at least **ci_min_speedup** times better than the best already realised index.

**size**: `Bytes`  
Memory usage in bytes of the index.

**speedup**: `Float`  
Estimated speedup. The basic value is the number of unique values in the argument of the clauses to which this index applies. The value is reduced when there are clauses with a variable in this position and when the standard deviation for the sets of clauses with the same value is high, i.e., we prefer indexes where the number of candidate clauses is similar, regardless of the value used in the call.

**Note:** This predicate property should be used for analysis and statistics only. The exact representation of `Indexes` may change between versions. The utilities [jiti_list/0](prologjiti.html#jiti_list/0) [jiti_list/1](prologjiti.html#jiti_list/1) list the *jit* indexes of matching predicates in a user friendly way.

**interpreted**  
True if the predicate is defined in Prolog. We return true on this because, although the code is actually compiled, it is completely transparent, just like interpreted code.

**iso**  
True if the predicate is covered by the ISO standard (ISO/IEC 13211-1).

**line_count**(`LineNumber`)  
Unify `LineNumber` with the line number of the first clause of the predicate. Fails if the predicate is not associated with a file. See also [source_file/2](consulting.html#source_file/2). See also the `file` property above, notably the reference to [clause_property/2](examineprog.html#clause_property/2).

**multifile**  
True if there may be multiple (or no) files providing clauses for the predicate. This property is set using [multifile/1](dynamic.html#multifile/1).

**meta_predicate**(`Head`)  
If the predicate is declared as a meta-predicate using [meta_predicate/1](metapred.html#meta_predicate/1), unify `Head` with the head-pattern. The head-pattern is a compound term with the same name and arity as the predicate where each argument of the term is a meta-predicate specifier. See [meta_predicate/1](metapred.html#meta_predicate/1) for details.

**mode**(`Head`)  
If the mode for the predicate is defined using [mode/1](dynamic.html#mode/1). Head is a term as in the property `meta_predicate(Head)`, but the specifiers are limited to `+`, `-` and `?`.

**monotonic**  
True if the predicate is tabled or dynamic using monotonic propagation. See [section 7.8](tabling-monotonic.html#sec:7.8).

**nodebug**  
Details of the predicate are not shown by the debugger. This is the default for built-in predicates. User predicates can be compiled this way using the Prolog flag [generate_debug_info](flags.html#flag:generate_debug_info).

**non_terminal**  
True if the predicate implements a *grammar rule*. See [non_terminal/1](dynamic.html#non_terminal/1).

**notrace**  
Do not show ports of this predicate in the debugger.

**number_of_clauses**(`ClauseCount`)  
Unify `ClauseCount` to the number of clauses associated with the predicate. Fails for foreign predicates. This property respects the *logical update view* and counts visible clauses at the moment the predicate was started.

**number_of_rules**(`RuleCount`)  
Similar to `number_of_clauses(ClauseCount)`, but only counts *rules*. A *rule* is defined as a clauses that has a body that is not just `true` (i.e., a *fact*).

**last_modified_generation**(`Generation`)  
Database generation at which the predicate was modified for the last time. Intended to quickly assesses the validity of caches.

**opaque**  
This property applies to dynamic and tabled predicates. For dynamic predicates it (temporary) stops propagating updates to dependent incrementally or monotonic tabled predicates. For tabled predicates it is not an error for an opaque predicate to depend on incremental or monotonic dynamic or tabled predicates.

**public**  
Predicate is declared public using [public/1](dynamic.html#public/1). Note that without further definition, public predicates are considered undefined and this property is *not* reported.

**quasi_quotation_syntax**  
The predicate (with arity 4) is declared to provide quasi quotation syntax with [quasi_quotation_syntax/1](quasiquotations.html#quasi_quotation_syntax/1).

**size**(`Bytes`)  
Memory used for this predicate. This includes the memory of the predicate header, the combined memory of all clauses including erased but not yet garbage collected clauses (see [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0) and [clause_property/2](examineprog.html#clause_property/2)) and the memory used by clause indexes (see the `indexed(Indexes)` property. *Excluded* are *lingering* data structures. These are garbage data structures that have been detached from the predicate but cannot yet be reclaimed because they may be in use by some thread.

**ssu**  
The predicate has been defined using *single sided unification* rules. See [section 5.6](ssu.html#sec:5.6).

**static**  
The definition can *not* be modified using [assertz/1](db.html#assertz/1) and friends. This property is the opposite from `dynamic`, i.e., for each defined predicate, either `static` or `dynamic` is true but never both.

**tabled**  
True of the predicate is *tabled*. The `tabled(?Flag)` property can be used to obtain details about how the predicate is tabled.

**tabled**(`?Flag`)  
True of the predicate is *tabled* and `Flag` applies. Any tabled predicate has one of the mutually exclusive flags `variant` or `subsumptive`. In addition, tabled predicates may have one or more of the following flags

**shared**  
The table is shared between threads. See [section 7.9](tabling-shared.html#sec:7.9).

**incremental**  
The table is subject to *incremental tabling*. See [section 7.7](tabling-incremental.html#sec:7.7)

Use the `tabled` property to enumerate all tabled predicates. See [table/1](tabling-preds.html#table/1) for details.

**thread_local**  
If true (only possible on the multithreaded version) each thread has its own clauses for the predicate. This property is set using [thread_local/1](threadcom.html#thread_local/1).

**transparent**  
True if the predicate is declared transparent using the [module_transparent/1](ctxmodule.html#module_transparent/1) or [meta_predicate/1](metapred.html#meta_predicate/1) declaration. In the latter case the property `meta_predicate(Head)` is also provided. See [chapter 6](modules.html#sec:6) for details.

**undefined**  
True if a procedure definition block for the predicate exists, but there are no clauses for it and it is not declared dynamic or multifile. This is true if the predicate occurs in the body of a loaded predicate, an attempt to call it has been made via one of the meta-call predicates, the predicate has been declared as e.g., a meta-predicate or the predicate had a definition in the past. Originally used to find missing predicate definitions. The current implementation of [list_undefined/0](check.html#list_undefined/0) used cross-referencing. Deprecated.

**visible**  
True when predicate can be called without raising a predicate existence error. This means that the predicate is (1) defined, (2) can be inherited from one of the default modules (see [default_module/2](importmodule.html#default_module/2)) or (3) can be autoloaded. The behaviour is logically consistent iff the property `visible` is provided explicitly. If the property is left unbound, only defined predicates are enumerated.

**volatile**  
If true, the clauses are not saved into a saved state by [qsave_program/\[1,2\]](saved-states.html#qsave_program/1). This property is set using [volatile/1](saved-states.html#volatile/1).

**dwim_predicate**(`+Term, -Dwim`)  
‘Do What I Mean’(\`dwim’) support predicate. `Term` is a term, whose name and arity are used as a predicate specification. `Dwim` is instantiated with the most general term built from `Name` and the arity of a defined predicate that matches the predicate specified by `Term` in the‘Do What I Mean’sense. See [dwim_match/2](miscpreds.html#dwim_match/2) for‘Do What I Mean’string matching. Internal system predicates are not generated, unless the access level is `system` (see [access_level](flags.html#flag:access_level)). Backtracking provides all alternative matches.

\[ISO\]**clause**(`:Head, ?Body`)  
True if `Head` can be unified with a clause head and `Body` with the corresponding clause body. Gives alternative clauses on backtracking. For facts, `Body` is unified with the atom `true`. Note that SWI-Prolog allows [clause/2](examineprog.html#clause/2) to work on both dynamic and static code.^(96Using [clause/2](examineprog.html#clause/2) is disallowed if either the flag [iso](flags.html#flag:iso) or [protect_static_code](flags.html#flag:protect_static_code) is `true`.) Note that [clause/2](examineprog.html#clause/2) *decompiles* the actual clause and may return a clause that is different from the source or asserted clause, i.e., [clause/2](examineprog.html#clause/2) only promises *semantic equivalence*.

**clause**(`:Head, ?Body, ?Reference`)  
Equivalent to [clause/2](examineprog.html#clause/2), but unifies `Reference` with a unique reference to the clause (see also [assert/2](db.html#assert/2), [erase/1](db.html#erase/1)). If `Reference` is instantiated to a reference the clause's head and body will be unified with `Head` and `Body`. The `Reference` is a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)), which implies it is subject to *atom garbage collection*. The `Reference` provides safe access to the clause while it exists and generates a reliable existence_error exception after the clause has been erased.

**nth_clause**(`?Pred, ?Index, ?Reference`)  
Provides access to the clauses of a predicate using their index number. Counting starts at 1. If `Reference` is specified it unifies `Pred` with the most general term with the same name/arity as the predicate and `Index` with the index number of the clause. Otherwise the name and arity of `Pred` are used to determine the predicate. If `Index` is provided, `Reference` will be unified with the clause reference. If `Index` is unbound, backtracking will yield both the indexes and the references of all clauses of the predicate. The following example finds the 2nd clause of [append/3](lists.html#append/3):

``` code
?- use_module(library(lists)).
...
?- nth_clause(append(_,_,_), 2, Ref), clause(Head, Body, Ref).
Ref = <clause>(0x994290),
Head = lists:append([_G23|_G24], _G21, [_G23|_G27]),
Body = append(_G24, _G21, _G27).
```

**clause_property**(`+ClauseRef, -Property`)  
Queries properties of a clause. `ClauseRef` is a reference to a clause as produced by [clause/3](examineprog.html#clause/3), [nth_clause/3](examineprog.html#nth_clause/3) or [prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3). Unlike most other predicates that access clause references, [clause_property/2](examineprog.html#clause_property/2) may be used to get information about erased clauses that have not yet been reclaimed. `Property` is one of the following:

**file**(`FileName`)  
Unify `FileName` with the name of the file from which the clause is loaded. Fails if the clause was not created by loading a file (e.g., clauses added using [assertz/1](db.html#assertz/1)). See also `source`.

**line_count**(`LineNumber`)  
Unify `LineNumber` with the line number of the clause. Fails if the clause is not associated to a file.

**size**(`SizeInBytes`)  
True when `SizeInBytes` is the size that the clause uses in memory in bytes. The size required by a predicate also includes the predicate data record, a linked list of clauses, clause selection instructions and optionally one or more clause indexes.

**source**(`FileName`)  
Unify `FileName` with the name of the source file that created the clause. This is the same as the `file` property, unless the file is loaded from a file that is textually included into source using [include/1](consulting.html#include/1). In this scenario, `file` is the included file, while the `source` property refers to the *main* file.

**fact**  
True if the clause has no body.

**erased**  
True if the clause has been erased, but not yet reclaimed because it is referenced.

**predicate**(`PredicateIndicator`)  
`PredicateIndicator` denotes the predicate to which this clause belongs. This is needed to obtain information on erased clauses because the usual way to obtain this information using [clause/3](examineprog.html#clause/3) fails for erased clauses.

**module**(`Module`)  
`Module` is the context module used to execute the body of the clause. For normal clauses, this is the same as the module in which the predicate is defined. However, if a clause is compiled with a module qualified *head*, the clause belongs to the predicate with the qualified head, while the body is executed in the context of the module in which the clause was defined.
