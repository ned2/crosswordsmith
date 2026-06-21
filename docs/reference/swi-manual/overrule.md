
## 6.6 Overruling Module Boundaries

The module system described so far is sufficient to distribute programs over multiple modules. There are, however, cases in which we would like to be able to overrule this schema and explicitly call a predicate in some module or assert explicitly into some module. Calling in a particular module is useful for debugging from the user's top level or to access multiple implementations of the same interface that reside in multiple modules. Accessing the same interface from multiple modules cannot be achieved using importing because importing a predicate with the same name and arity from two modules results in a name conflict. Asserting in a different module can be used to create models dynamically in a new module. See [section 6.13](dynamic-modules.html#sec:6.13).

Direct addressing of modules is achieved using a `:``/2` explicitly in a program and relies on the module qualification mechanism described in [section 6.5](metapred.html#sec:6.5). Here are a few examples:

``` code
?- assert(world:done).   % asserts done/0 into module world
?- world:asserta(done).  % the same
?- world:done.           % calls done/0 in module world
```

Note that the second example is the same due to the Prolog flag [colon_sets_calling_context](flags.html#flag:colon_sets_calling_context). The system predicate [asserta/1](db.html#asserta/1) is called in the module `world`, which is possible because system predicates are *visible* in all modules. At the same time, the *calling context* is set to `world`. Because meta arguments are qualified with the calling context, the resulting call is the same as the first example.

### 6.6.1 Explicit manipulation of the calling context

Quintus’derived module systems have no means to separate the lookup module (for finding predicates) from the calling context (for qualifying meta arguments). Some other Prolog implementations (e.g., ECLiPSe and IF/Prolog) distinguish these operations, using `@/2` for setting the calling context of a goal. This is provided by SWI-Prolog, currently mainly to support compatibility layers.

**@**(`:Goal, +Module`)  
Execute `Goal`, setting the calling context to `Module`. Setting the calling context affects meta-predicates, for which meta arguments are qualified with `Module` and transparent predicates (see [module_transparent/1](ctxmodule.html#module_transparent/1)). It has no implications for other predicates.

For example, the code `asserta(done)@world` is the same as `asserta(world:done)`. Unlike in `world:asserta(done)`, [asserta/1](db.html#asserta/1) is resolved in the current module rather than the module `world`. This makes no difference for system predicates, but usually does make a difference for user predicates.

Not that SWI-Prolog does not define `@` as an operator. Some systems define this construct using `op(900, xfx, @)`.
