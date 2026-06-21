
## 6.15 Module properties

The following predicates can be used to query the module system for reflexive programming:

\[nondet\]**current_module**(`?Module`)  
True if `Module` is a currently defined module. This predicate enumerates all modules, whether loaded from a file or created dynamically. Note that modules cannot be destroyed in the current version of SWI-Prolog.

**module_property**(`?Module, ?Property`)  
True if `Property` is a property of `Module`. Defined properties are:

**class**(`-Class`)  
True when `Class` is the class of the module. Defined classes are

**user**  
Default for user-defined modules.

**system**  
Module `system` and modules from `<``home``>/boot`.

**library**  
Other modules from the system directories.

**temporary**  
Module is temporary.

**test**  
Modules that create tests.

**development**  
Modules that only support the development environment.

**file**(`?File`)  
True if `Module` was loaded from `File`.

**line_count**(`-Line`)  
True if `Module` was loaded from the N-th line of file.

**exports**(`-ListOfPredicateIndicators`)  
True if `Module` exports the given predicates. Predicate indicators are in canonical form (i.e., always using name/arity and never the DCG form name//arity). Future versions may also use the DCG form. See also [predicate_property/2](examineprog.html#predicate_property/2). Succeeds with an empty list if the module exports no predicates.

**exported_operators**(`-ListOfOperators`)  
True if `Module` exports the given operators. Each exported operator is represented as a term `op(Pri,Assoc,Name)`. Succeeds with an empty list if the module exports no operators.

**size**(`-Bytes`)  
Total size in bytes used to represent the module. This includes the module itself, its (hash) tables and the summed size of all predicates defined in this module. See also the `size(Bytes)` property in [predicate_property/2](examineprog.html#predicate_property/2).

**program_size**(`-Bytes`)  
Memory (in bytes) used for storing the predicates of this module. This figure includes the predicate header and clauses.

**program_space**(`-Bytes`)  
If present, this number limits the `program_size`. See [set_module/1](manipmodule.html#set_module/1).

**last_modified_generation**(`-Generation`)  
Integer expression the last database generation where a clause was added or removed from a predicate that is implemented in this module. See also [predicate_property/2](examineprog.html#predicate_property/2).

**set_module**(`:Property`)  
Modify properties of the module. Currently, the following properties may be modified:

**base**(`+Base`)  
Set the default import module of the current module to `Module`. Typically, `Module` is one of `user` or `system`. See [section 6.10](importmodule.html#sec:6.10).

**class**(`+Class`)  
Set the class of the module. See [module_property/2](manipmodule.html#module_property/2).

**program_space**(`+Bytes`)  
Maximum amount of memory used to store the predicates defined inside the module. Raises a permission error if the current usage is above the requested limit. Setting the limit to 0 (zero) removes the limit. An attempt to assert clauses that causes the limit to be exceeded causes a `resource_error(program_space)` exception. See [assertz/1](db.html#assertz/1) and [module_property/2](manipmodule.html#module_property/2).
