
## 6.13 Dynamic Modules

So far, we discussed modules that were created by loading a module file. These modules have been introduced to facilitate the development of large applications. The modules are fully defined at load-time of the application and normally will not change during execution. Having the notion of a set of predicates as a self-contained world can be attractive for other purposes as well. For example, assume an application that can reason about multiple worlds. It is attractive to store the data of a particular world in a module, so we extract information from a world simply by invoking goals in this world.

Dynamic modules can easily be created. Any built-in predicate that tries to locate a predicate in a specific module will create this module as a side-effect if it did not yet exist. For example:

``` code
?- assert(world_a:consistent),
   set_prolog_flag(world_a:unknown, fail).
```

These calls create a module called‘world_a’and make the call‘world_a:consistent’succeed. Undefined predicates will not raise an exception for this module (see [unknown](flags.html#flag:unknown)).

Import and export from a dynamically created world can be achieved using [import/1](altmoduleapi.html#import/1) and [export/1](altmoduleapi.html#export/1) or by specifying the import module as described in [section 6.10](importmodule.html#sec:6.10).

``` code
?- world_b:export(solve/2).          % exports solve/2 from world_b
?- world_c:import(world_b:solve/2).  % and import it to world_c
```
