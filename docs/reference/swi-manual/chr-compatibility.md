
## 9.6 CHR compatibility

### 9.6.1 The Old SICStus CHR implementation

There are small differences between the current K.U.Leuven CHR system in SWI-Prolog, older versions of the same system, and SICStus’CHR system.

The current system maps old syntactic elements onto new ones and ignores a number of no longer required elements. However, for each a *deprecated* warning is issued. You are strongly urged to replace or remove deprecated features.

Besides differences in available options and pragmas, the following differences should be noted:

- *The constraints/1 declaration*  
  This declaration is deprecated. It has been replaced with the [chr_constraint/1](practical.html#chr_constraint/1) declaration.
- *The option/2 declaration*  
  This declaration is deprecated. It has been replaced with the [chr_option/2](chr-syntaxandsemantics.html#chr_option/2) declaration.
- *The handler/1 declaration*  
  In SICStus every CHR module requires a handler/1 declaration declaring a unique handler name. This declaration is valid syntax in SWI-Prolog, but will have no effect. A warning will be given during compilation.
- *The rules/1 declaration*  
  In SICStus, for every CHR module it is possible to only enable a subset of the available rules through the rules/1 declaration. The declaration is valid syntax in SWI-Prolog, but has no effect. A warning is given during compilation.
- *Guard bindings*  
  The `check_guard_bindings` option only turns invalid calls to unification into failure. In SICStus this option does more: it intercepts instantiation errors from Prolog built-ins such as [is/2](arith.html#is/2) and turns them into failure. In SWI-Prolog, we do not go this far, as we like to separate concerns more. The CHR compiler is aware of the CHR code, the Prolog system, and the programmer should be aware of the appropriate meaning of the Prolog goals used in guards and bodies of CHR rules.

### 9.6.2 The Old ECLiPSe CHR implementation

The old ECLiPSe CHR implementation features a label_with/1 construct for labeling variables in CHR constraints. This feature has long since been abandoned. However, a simple transformation is all that is required to port the functionality.

``` code
label_with Constraint1 if Condition1.
...
label_with ConstraintN if ConditionN.
Constraint1 :- Body1.
...
ConstraintN :- BodyN.
```

is transformed into

``` code
:- chr_constraint my_labeling/0.

my_labeling \ Constraint1 <=> Condition1 | Body1.
...
my_labeling \ ConstraintN <=> ConditionN | BodyN.
my_labeling <=> true.
```

Be sure to put this code after all other rules in your program! With my_labeling/0 (or another predicate name of your choosing) the labeling is initiated, rather than ECLiPSe's chr_labeling/0 .
