
## 6.12 An alternative import/export interface

The [use_module/1](import.html#use_module/1) predicate from [section 6.3](import.html#sec:6.3) defines import and export relations based on the filename from which a module is loaded. If modules are created differently, such as by asserting predicates into a new module as described in [section 6.13](dynamic-modules.html#sec:6.13), this interface cannot be used. The interface below provides for import/export from modules that are not created using a module file.

**export**(`+PredicateIndicator, ...`)  
Add predicates to the public list of the context module. This implies the predicate will be imported into another module if this module is imported with [use_module/\[1,2\]](import.html#use_module/1). Note that predicates are normally exported using the directive [module/2](defmodule.html#module/2). [export/1](altmoduleapi.html#export/1) is meant to handle export from dynamically created modules.

**import**(`+PredicateIndicator, ...`)  
Import predicates `PredicateIndicator` into the current context module. `PredicateIndicator` must specify the source module using the \<`module`\>:\<`pi`\> construct. Note that predicates are normally imported using one of the directives [use_module/\[1,2\]](import.html#use_module/1). The [import/1](altmoduleapi.html#import/1) alternative is meant for handling imports into dynamically created modules. See also [export/1](altmoduleapi.html#export/1) and export_list/2.
