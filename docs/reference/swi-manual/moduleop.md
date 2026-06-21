
## 6.9 Operators and modules

Operators ([section 4.25](operators.html#sec:4.25)) are local to modules, where the initial table behaves as if it is copied from the module `user` (see [section 6.11](resmodules.html#sec:6.11)). A specific operator can be disabled inside a module using `:- op(0, Type, Name)`. Inheritance from the public table can be restored using `:- op(-1, Type, Name)`.

In addition to using the [op/3](operators.html#op/3) directive, operators can be declared in the [module/2](defmodule.html#module/2) directive as shown below. Such operator declarations are visible inside the module, and importing such a module makes the operators visible in the target module. Exporting operators is typically used by modules that implement sub-languages such as chr (see [chapter 9](chr.html#sec:9)). The example below is copied from the library `library(clpfd)`.

``` code
:- module(clpfd,
          [ op(760, yfx, #<==>),
            op(750, xfy, #==>),
            op(750, yfx, #<==),
            op(740, yfx, #\/),
            ...
            (#<==>)/2,
            (#==>)/2,
            (#<==)/2,
            (#\/)/2,
            ...
          ]).
```
