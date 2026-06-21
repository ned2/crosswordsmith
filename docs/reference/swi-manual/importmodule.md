
## 6.10 Dynamic importing using import modules

Until now we discussed the public module interface that is, at least to some extent, portable between Prolog implementations with a module system that is derived from Quintus Prolog. The remainder of this chapter describes the underlying mechanisms that can be used to emulate other module systems or implement other code-reuse mechanisms.

In addition to built-in predicates, imported predicates and locally defined predicates, SWI-Prolog modules can also call predicates from its *import modules*. Each module has a (possibly empty) list of import modules. In the default setup, each new module has a single import module, which is `user` for all normal user modules and `system` for all system library modules. Module `user` imports from `system` where all built-in predicates reside. These special modules are described in more detail in [section 6.11](resmodules.html#sec:6.11).

In general, the import relations between modules form an acyclic directed graph. The import relation affects the following mechanisms:

**Predicate visibility**  
When looking for a specific predicate definition the system starts in the target module. If the predicate is undefined there it walks the module import relations depth-first left-to-right searching for a module that defines the predicate. The first encountered definition is used. Note that using the default setup this means it searches the `user` and `system` modules (in that order).

**Operators**  
Operators are also searched through the import relations. System operators are defined in the module `system`. The user may define operators in `user` to make them globally visible for compatibility with e.g., SICStus Prolog that has no local operators. Normally operators are defined in a module and, when applicable, exported using the [module/2](defmodule.html#module/2) module header.

**The [unknown](flags.html#flag:unknown) flag**  
This flag controls the response to encountering an undefined predicate in the target module.

**Term and goal expansion**  
The hooks [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2) (see [section 4.3.1](consulting.html#sec:4.3.1)) are *chained* over the import modules that define these hooks. This implies we collect all modules that provide definitions for these hook predicates by traversing the import module relation depth-first and left-to-right. Next, we perform the transformations in a *pipeline*, starting at the target module.

The list of import modules for a specific module can be manipulated and queried using the following predicates, as well as using [set_module/1](manipmodule.html#set_module/1).

\[nondet\]**import_module**(`+Module, -Import`)  
True if `Module` inherits directly from `Import`. All normal modules only import from `user`, which imports from `system`. The predicates [add_import_module/3](importmodule.html#add_import_module/3) and [delete_import_module/2](importmodule.html#delete_import_module/2) can be used to manipulate the import list. See also [default_module/2](importmodule.html#default_module/2).

\[multi\]**default_module**(`+Module, -Default`)  
True if predicates and operators in `Default` are visible in `Module`. Modules are returned in the same search order used for predicates and operators. That is, `Default` is first unified with `Module`, followed by the depth-first transitive closure of [import_module/2](importmodule.html#import_module/2).

**add_import_module**(`+Module, +Import, +StartOrEnd`)  
If `Import` is not already an import module for `Module`, add it to this list at the `start` or `end` depending on `StartOrEnd`. See also [import_module/2](importmodule.html#import_module/2) and [delete_import_module/2](importmodule.html#delete_import_module/2).

**delete_import_module**(`+Module, +Import`)  
Delete `Import` from the list of import modules for `Module`. Fails silently if `Import` is not in the list.
