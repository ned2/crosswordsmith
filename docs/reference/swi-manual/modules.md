
# 6 Modules

A Prolog module is a collection of predicates which defines a public interface by means of a set of provided predicates and operators. Prolog modules are defined by an ISO standard. Unfortunately, the standard is considered a failure and, as far as we are aware, not implemented by any concrete Prolog implementation. The SWI-Prolog module system syntax is derived from the Quintus Prolog module system. The Quintus module system has been the starting point for the module systems of a number of mainstream Prolog systems, such as SICStus, Ciao and YAP. The underlying primitives of the SWI-Prolog module system differ from the mentioned systems. These primitives allow for multiple modules in a file, hierarchical modules, emulation of other modules interfaces, etc.

This chapter motivates and describes the SWI-Prolog module system. Novices can start using the module system after reading [section 6.2](defmodule.html#sec:6.2) and [section 6.3](import.html#sec:6.3). The primitives defined in these sections suffice for basic usage until one needs to export predicates that call or manage other predicates dynamically (e.g., use [call/1](metacall.html#call/1), [assert/1](db.html#assert/1), etc.). Such predicates are called *meta predicates* and are discussed in [section 6.5](metapred.html#sec:6.5). [Section 6.6](overrule.html#sec:6.6) to [section 6.9](moduleop.html#sec:6.9) describe more advanced issues. Starting with [section 6.10](importmodule.html#sec:6.10), we discuss more low-level aspects of the SWI-Prolog module system that are used to implement the visible module system, and can be used to build other code reuse mechanisms.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[6.1 Why Use Modules?](whymodules.html)

[6.2 Defining a Module](defmodule.html)

[6.3 Importing Predicates into a Module](import.html)

[6.4 Controlled autoloading for modules](module-autoload.html)

[6.5 Defining a meta-predicate](metapred.html)

[6.6 Overruling Module Boundaries](overrule.html)

[6.6.1 Explicit manipulation of the calling context](overrule.html#sec:6.6.1)

[6.7 Interacting with modules from the top level](mtoplevel.html)

[6.8 Composing modules from other modules](reexport.html)

[6.9 Operators and modules](moduleop.html)

[6.10 Dynamic importing using import modules](importmodule.html)

[6.11 Reserved Modules and using the‘user’module](resmodules.html)

[6.12 An alternative import/export interface](altmoduleapi.html)

[6.13 Dynamic Modules](dynamic-modules.html)

[6.14 Transparent predicates: definition and context module](ctxmodule.html)

[6.15 Module properties](manipmodule.html)

[6.16 Compatibility of the Module System](modulecompat.html)
