
## 6.11 Reserved Modules and using the‘user’module

As mentioned above, SWI-Prolog contains two special modules. The first one is the module `system`. This module contains all built-in predicates. Module `system` has no import module, i.e., is a *root* of the module graph. The second special module is the module `user`. This module forms the initial working space of the user. Initially it is empty.^(186Unfortunately some *hooks* are traditionally defined in the user module). The import module of module `user` is `system`, making all built-in predicates available.

All normal application modules import from the module `user`. This implies they can use all predicates imported into `user` without explicitly importing them. If an application loads all modules from the `user` module using [use_module/1](import.html#use_module/1), one achieves a scoping system similar to the C-language, where every module can access all exported predicates without any special precautions.

All *library* modules (see [module_property/2](manipmodule.html#module_property/2)) import directly from `system`. Library modules are modules loaded from the SWI-Prolog installation. As they import from `system`, the functionality of a library is not affected by operator or predicate definitions in the `user` module.
