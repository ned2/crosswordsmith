
## 6.2 Defining a Module

Modules are normally created by loading a *module file*. A module file is a file holding a [module/2](defmodule.html#module/2) directive as its first term. The [module/2](defmodule.html#module/2) directive declares the name and the public (i.e., externally visible) predicates of the module. The rest of the file is loaded into the module. Below is an example of a module file, defining [reverse/2](lists.html#reverse/2) and hiding the helper predicate rev/3. A module can use all built-in predicates and, by default, cannot redefine system predicates.

``` code
:- module(reverse, [reverse/2]).

reverse(List1, List2) :-
        rev(List1, [], List2).

rev([], List, List).
rev([Head|List1], List2, List3) :-
        rev(List1, [Head|List2], List3).
```

The module is named `reverse`. Typically, the name of a module is the same as the name of the file by which it is defined without the filename extension, but this naming is not enforced. Modules are organised in a single and flat namespace and therefore module names must be chosen with some care to avoid conflicts. As we will see, typical applications of the module system rarely use the name of a module explicitly in the source text.

:- **module**(`+Module, +PublicList`)  
This directive can only be used as the first term of a source file. It declares the file to be a *module file*, defining a module named `Module`. Note that a module name is an atom. The module exports the predicates of `PublicList`. `PublicList` is a list of predicate indicators (name/arity or name//arity pairs) or operator declarations using the format `op(Precedence, Type, Name)`. Operators defined in the export list are available inside the module as well as to modules importing this module. See also [section 4.25](operators.html#sec:4.25).

Compatible to Ciao Prolog, if `Module` is unbound, it is unified with the basename without extension of the file being loaded.

:- **module**(`+Module, +PublicList, +Dialect`)  
Same as [module/2](defmodule.html#module/2). The additional `Dialect` argument provides a list of *language options*. Each atom in the list `Dialect` is mapped to a [use_module/1](import.html#use_module/1) goal as given below. See also [section C](dialect.html#sec:C). The third argument is supported for compatibility with the [Prolog Commons project](http://prolog-commons.org/).

``` code
:- use_module(library(dialect/LangOption)).
```
