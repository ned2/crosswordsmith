
## 6.4 Controlled autoloading for modules

SWI-Prolog by default support *autoloading* from its standard library. Autoloading implies that when a predicate is found missing during execution the library is searched and the predicate is imported lazily using [use_module/2](import.html#use_module/2). See [section 2.14](autoload.html#sec:2.14) for details.

The advantage of autoloading is that it requires less typing while it reduces the startup time and reduces the memory footprint of an application. It also allows moving old predicates or emulation thereof the module `library(backcomp)` without affecting existing code. This procedure keeps the libraries and system clean. We make sure that there are not two modules that provide the same predicate as autoload predicate.

Nevertheless, a disadvantage of this autoloader is that the dependencies of a module on the libraries are not explicit and tooling such as PceEmacs or [gxref/0](xref.html#gxref/0) are required to find these dependencies. Some users want explicit control over which library predicates are accessed from where, preferably by using [use_module/2](import.html#use_module/2) which explicitly states which predicates are imported from which library.^(185Note that built-in predicates still add predicates for general use to all name spaces.)

Large applications typically contain source files that are not immediately needed and often are not needed at all in many runs of the program. This can be solved by creating an application-specific autoload library, but with multiple parties providing autoloadable predicates the maintenance becomes fragile. For these two reasons we added [autoload/1](module-autoload.html#autoload/1) and [autoload/2](module-autoload.html#autoload/2) that behave similar to [use_module/\[1,2\]](import.html#use_module/1), but do not perform the actual loading. The generic autoloader now proceeds as follows if a missing predicate is encountered:

1.  Check [autoload/2](module-autoload.html#autoload/2) declarations. If one specifies the predicate, import it using [use_module/2](import.html#use_module/2).
2.  Check [autoload/1](module-autoload.html#autoload/1) declarations. If the specified file is loaded, check its export list. Otherwise read the module declaration of the target file to find the exports. If the target predicate is found, import it using [use_module/2](import.html#use_module/2).
3.  Perform autoloading from the library if the [autoload](flags.html#flag:autoload) is `true`.

**autoload**(`:File`)  
**autoload**(`:File, +Imports`)  
Declare that possibly missing predicates in the module in which this declaration occurs are to be resolved by using [use_module/2](import.html#use_module/2) on `File` to (possibly) load the file and make the target predicate available. The [autoload/2](module-autoload.html#autoload/2) variant is tried before [autoload/1](module-autoload.html#autoload/1). It is not allowed for two [autoload/2](module-autoload.html#autoload/2) declarations to provide the same predicate and it is not allowed to define a predicate provided in this way locally. See also [require/1](consulting.html#require/1), which allows specifying predicates for autoloading from their default location.

Predicates made available using [autoload/2](module-autoload.html#autoload/2) behave as defined predicates, which implies that any operation on them will perform autoloading if necessary. Notably [predicate_property/2](examineprog.html#predicate_property/2), [current_predicate/1](examineprog.html#current_predicate/1) and [clause/2](examineprog.html#clause/2) are supported.

Currently, neither the existence of `File`, nor whether it actually exports the given predicates ([autoload/2](module-autoload.html#autoload/2)) is verified when the file is loaded. Instead, the declarations are verified when searching for a missing predicate.

If the Prolog flag [autoload](flags.html#flag:autoload) is set to `false`, these declarations are interpreted as [use_module/\[1,2\]](import.html#use_module/1).
