
## 6.3 Importing Predicates into a Module

Predicates can be added to a module by *importing* them from another module. Importing adds predicates to the namespace of a module. An imported predicate can be called exactly the same as a locally defined predicate, although its implementation remains part of the module in which it has been defined.

Importing the predicates from another module is achieved using the directives [use_module/1](import.html#use_module/1) or [use_module/2](import.html#use_module/2). Note that both directives take `filename(s)` as arguments. That is, modules are imported based on their filename rather than their module name.

**use_module**(`+Files`)  
Load the file(s) specified with `Files` just like [ensure_loaded/1](consulting.html#ensure_loaded/1). The files must all be module files. All exported predicates from the loaded files are imported into the module from which this predicate is called. This predicate is equivalent to [ensure_loaded/1](consulting.html#ensure_loaded/1), except that it raises an error if `Files` are not module files.

The imported predicates act as *weak symbols* in the module into which they are imported. This implies that a local definition of a predicate overrides (clobbers) the imported definition. If the flag [warn_override_implicit_import](flags.html#flag:warn_override_implicit_import) is `true` (default), a warning is printed. Below is an example of a module that uses library(lists), but redefines [flatten/2](lists.html#flatten/2), giving it a totally different meaning:

``` code
:- module(shapes, []).
:- use_module(library(lists)).

flatten(cube, square).
flatten(ball, circle).
```

Loading the above file prints the following message:

``` code
Warning: /home/janw/Bugs/Import/t.pl:5:
        Local definition of shapes:flatten/2
        overrides weak import from lists
```

This warning can be avoided by (1) using [use_module/2](import.html#use_module/2) to only import the predicates from the `lists` library that are actually used in the‘shapes’module, (2) using the `except([`[`flatten/2`](lists.html#flatten/2)`])` option of [use_module/2](import.html#use_module/2), (3) use `:- abolish(`[`flatten/2`](lists.html#flatten/2)`).` before the local definition or (4) setting [warn_override_implicit_import](flags.html#flag:warn_override_implicit_import) to `false`. Globally disabling this warning is only recommended if overriding imported predicates is common as a result of design choices or the program is ported from a system that silently overrides imported predicates.

Note that it is always an error to import two modules with [use_module/1](import.html#use_module/1) that export the same predicate. Such conflicts must be resolved with [use_module/2](import.html#use_module/2) as described above.

**use_module**(`+File, +ImportList`)  
Load `File`, which must be a module file, and import the predicates as specified by `ImportList`. `ImportList` is a list of predicate indicators specifying the predicates that will be imported from the loaded module. `ImportList` also allows for renaming or import-everything-except. See also the `import` option of [load_files/2](consulting.html#load_files/2). The first example below loads [member/2](lists.html#member/2) from the `lists` library and [append/2](lists.html#append/2) under the name `list_concat`, which is how this predicate is named in YAP. The second example loads all exports from library `option` except for [meta_options/3](option.html#meta_options/3). These renaming facilities are generally used to deal with portability issues with as few changes as possible to the actual code. See also [section C](dialect.html#sec:C) and [section 6.8](reexport.html#sec:6.8).

``` code
:- use_module(library(lists), [ member/2,
                                append/2 as list_concat
                              ]).
:- use_module(library(option), except([meta_options/3])).
```

In most cases a module is imported because some of its predicates are being used. However, sometimes a module is imported for other reasons, e.g., for its declarations. In such cases it is best practice to use [use_module/2](import.html#use_module/2) with empty ImportList. This distinguishes an imported module that is used, although not for its predicates, from a module that is needlessly imported.

The [module/2](defmodule.html#module/2), [use_module/1](import.html#use_module/1) and [use_module/2](import.html#use_module/2) directives are sufficient to partition a simple Prolog program into modules. The SWI-Prolog graphical cross-referencing tool [gxref/0](xref.html#gxref/0) can be used to analyse the dependencies between non-module files and propose module declarations for each file.
