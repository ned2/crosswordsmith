
## 4.25 Operators

Operators are defined to improve the readability of source code. For example, without operators, to write `2*3+4*5` one would have to write `+(*(2,3),*(4,5))`. In Prolog, a number of operators have been predefined. All operators, except for the comma (,) can be redefined by the user.

Some care has to be taken before defining new operators. Defining too many operators might make your source‘natural’looking, but at the same time using many operators can make it hard to understand the limits of your syntax.

In SWI-Prolog, operators are local to the module in which they are defined. Operators can be exported from modules using a term `op(Precedence, Type, Name)` in the export list as specified by [module/2](defmodule.html#module/2). Many modern Prolog systems have module specific operators. Unfortunately, there is no established interface for exporting and importing operators. SWI-Prolog's convention has been adopted by YAP.

The module table of the module `user` acts as default table for all modules and can be modified explicitly from inside a module to achieve compatibility with other Prolog that do not have module-local operators:

``` code
:- module(prove,
          [ prove/1
          ]).

:- op(900, xfx, user:(=>)).
```

Although operators are module-specific and the predicates that define them ([op/3](operators.html#op/3)) or rely on them such as [current_op/3](operators.html#current_op/3), [read/1](termrw.html#read/1) and [write/1](termrw.html#write/1) are module sensitive, they are not proper meta-predicates. If they were proper meta predicates [read/1](termrw.html#read/1) and [write/1](termrw.html#write/1) would use the module from which they are called, breaking compatibility with other Prolog systems. The following rules apply:

1.  If the module is explicitly specified by qualifying the third argument ([op/3](operators.html#op/3), [current_op/3](operators.html#current_op/3)) or specifying a `module(Module)` option ([read_term/3](termrw.html#read_term/3), [write_term/3](termrw.html#write_term/3)), this module is used.
2.  While compiling, the module into which the compiled code is loaded applies.
3.  Otherwise, the *typein module* applies. This is normally `user` and may be changed using [module/1](mtoplevel.html#module/1).

In SWI-Prolog, a *quoted atom* never acts as an operator. Note that the portable way to stop an atom acting as an operator is to enclose it in parentheses like this: (myop). See also [section 5.3.1](ext-syntax.html#sec:5.3.1).

\[ISO\]**op**(`+Precedence, +Type, :Name`)  
Declare `Name` to be an operator of type `Type` with precedence `Precedence`. `Name` can also be a list of names, in which case all elements of the list are declared to be identical operators. `Precedence` is an integer between 0 and 1200. Precedence 0 removes the declaration. `Type` is one of: `xf`, `yf`, `xfx`, `xfy`, `yfx`, `fy` or `fx`. The‘`f`’indicates the position of the functor, while `x` and `y` indicate the position of the arguments.‘`y`’should be interpreted as “on this position a term with precedence lower or equal to the precedence of the functor should occur” . For‘`x`’the precedence of the argument must be strictly lower. The precedence of a term is 0, unless its principal functor is an operator, in which case the precedence is the precedence of this operator. A term enclosed in parentheses `( ... )` has precedence 0.

The predefined operators are shown in [table 5](operators.html#tab:operators). Operators can be redefined, unless prohibited by one of the limitations below. Applications must be careful with (re-)defining operators because changing operators may cause (other) files to be interpreted **differently**. Often this will lead to a syntax error. In other cases, text is read silently into a different term which may lead to subtle and difficult to track errors.

- It is not allowed to redefine the comma (`','`).
- The bar (`|`) can only be (re-)defined as infix operator with priority not less than 1001.

In SWI-Prolog, operators are *local* to a module (see also [section 6.9](moduleop.html#sec:6.9)). Keeping operators in modules and using controlled import/export of operators as described with the [module/2](defmodule.html#module/2) directive keep the issues manageable. The module `system` provides the operators from [table 5](operators.html#tab:operators) and these operators cannot be modified. Files that are loaded from the SWI-Prolog directories resolve operators and predicates from this `system` module rather than `user`, which makes the semantics of the library and development system modules independent of operator changes to the `user` module. See [section 4.25](operators.html#sec:4.25) for details about the relation between operators and modules.

|  |  |  |
|---:|----|----|
| 1200 | xfx | **`-->`**, **`:-`**, **`=>`**, **==\>** |
| 1200 | fx | **`:-`**, **`?-`** |
| 1150 | fx | **dynamic**, **discontiguous**, **initialization**, **meta_predicate**, **module_transparent**, **multifile**, **public**, **thread_local**, **thread_initialization**, **volatile** |
| 1105 | xfy | **`|`** |
| 1100 | xfy | **`;`** |
| 1050 | xfy | **`->`**, **`*->`** |
| 1000 | xfy | **`,`** |
| 990 | xfx | **:=** |
| 900 | fy | **`\+`** |
| 700 | xfx | **`<`**, **`=`**, **`=..`**, **`=@=`**, **`\=@=`**, **`=:=`**, **`=<`**, **`==`**, **`=\=`**, **`>`**, **`>=`**, **`@<`**, **`@=<`**, **`@>`**, **`@>=`**, **`\=`**, **`\==`**, **as**, **is**, **`>:<`**, **`:<`** |
| 600 | xfy | **`:`** |
| 500 | yfx | **`+`**, **`-`**, **`/\`**, **`\/`**, **xor** |
| 500 | fx | **`?`** |
| 400 | yfx | **`*`**, **`/`**, **`//`**, **div**, **rdiv**, **`<<`**, **`>>`**, **mod**, **rem** |
| 200 | xfx | **`**`** |
| 200 | xfy | **`^`** |
| 200 | fy | **`+`**, **`-`**, **`\`** |
| 100 | yfx | **`.`** |
| 1 | fx | **`$`** |

**Table 5 :** System operators

\[ISO\]**current_op**(`?Precedence, ?Type, ?:Name`)  
True if `Name` is currently defined as an operator of type `Type` with precedence `Precedence`. See also [op/3](operators.html#op/3). Note that an *unqualified* `Name` does **not** resolve to the calling context but, when compiling, to the compiler's target module and otherwise to the *typein module*. See [section 4.25](operators.html#sec:4.25) for details.
