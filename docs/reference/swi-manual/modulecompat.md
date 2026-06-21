
## 6.16 Compatibility of the Module System

The SWI-Prolog module system is largely derived from the Quintus Prolog module system, which is also adopted by SICStus, Ciao and YAP. Originally, the mechanism for defining meta-predicates in SWI-Prolog was based on the [module_transparent/1](ctxmodule.html#module_transparent/1) directive and [strip_module/3](ctxmodule.html#strip_module/3). Since 5.7.4 it supports the de-facto standard [meta_predicate/1](metapred.html#meta_predicate/1) directive for implementing meta-predicates, providing much better compatibility.

The support for the [meta_predicate/1](metapred.html#meta_predicate/1) mechanism, however, is considerably different. On most systems, the *caller* of a meta-predicate is compiled differently to provide the required \<`module`\>:\<`term`\> qualification. This implies that the meta-declaration must be available to the compiler when compiling code that calls a meta-predicate. In practice, this implies that other systems pose the following restrictions on meta-predicates:

- Modules that provide meta-predicates for a module to be compiled must be loaded explicitly by that module.
- The meta-predicate directives of exported predicates must follow the [module/2](defmodule.html#module/2) directive immediately.
- After changing a meta-declaration, all modules that *call* the modified predicates need to be recompiled.

In SWI-Prolog, meta-predicates are also *module-transparent*, and qualifying the module-sensitive arguments is done inside the meta-predicate. As a result, the caller need not be aware that it is calling a meta-predicate and none of the above restrictions hold for SWI-Prolog. However, code that aims at portability must obey the above rules.

Other differences are listed below.

- If a module does not define a predicate, it is searched for in the *import modules*. By default, the import module of any user-defined module is the `user` module. In turn, the `user` module imports from the module `system` that provides all built-in predicates. The auto-import hierarchy can be changed using [add_import_module/3](importmodule.html#add_import_module/3) and [delete_import_module/2](importmodule.html#delete_import_module/2).

  This mechanism can be used to realise a simple object-oriented system or a hierarchical module system.

- Operator declarations are local to a module and may be exported. In Quintus and SICStus all operators are global. YAP and Ciao also use local operators. SWI-Prolog provides global operator declarations from within a module by explicitly qualifying the operator name with the `user` module. I.e., operators are inherited from the *import modules* (see above).

  ``` code
  :- op(precedence, type, user:(operatorname)).
  ```
