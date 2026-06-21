
## 6.7 Interacting with modules from the top level

Debugging often requires interaction with predicates that reside in modules: running them, setting spy points on them, etc. This can be achieved using the \<`module`\>:\<`term`\> construct explicitly as described above. In SWI-Prolog, you may also wish to omit the module qualification. Setting a spy point ([spy/1](debugger.html#spy/1)) on a plain predicate sets a spy point on any predicate with that name in any module. Editing ([edit/1](edit.html#edit/1)) or calling an unqualified predicate invokes the DWIM (Do What I Mean) mechanism, which generally suggests the correct qualified query.

Mainly for compatibility, we provide [module/1](mtoplevel.html#module/1) to switch the module with which the interactive top level interacts:

**module**(`+Module`)  
The call `module(``Module``)` may be used to switch the default working module for the interactive top level (see [prolog/0](toplevel.html#prolog/0)). This may be used when debugging a module. The example below lists the clauses of file_of_label/2 in the module `tex`.

``` code
1 ?- module(tex).
true.
tex: 2 ?- listing(file_of_label/2).
...
```
