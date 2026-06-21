
## 6.8 Composing modules from other modules

The predicates in this section are intended to create new modules from the content of other modules. Below is an example to define a *composite* module. The example exports all public predicates of `module_1`, `module_2` and `module_3`, pred/1 from `module_4`, all predicates from `module_5` except do_not_use/1 and all predicates from `module_6` while renaming pred/1 into mypred/1.

``` code
:- module(my_composite, []).
:- reexport([ module_1,
              module_2,
              module_3
            ]).
:- reexport(module_4, [ pred/1 ]).
:- reexport(module_5, except([do_not_use/1])).
:- reexport(module_6, except([pred/1 as mypred])).
```

**reexport**(`+Files`)  
Load and import predicates as [use_module/1](import.html#use_module/1) and re-export all imported predicates. The reexport declarations must immediately follow the module declaration.

**reexport**(`+File, +Import`)  
Import from `File` as [use_module/2](import.html#use_module/2) and re-export the imported predicates. All formats accepted by [use_module/2](import.html#use_module/2) for `Import` are accepted. The reexport declarations must immediately follow the module declaration.
