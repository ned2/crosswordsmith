
## 9.7 CHR Programming Tips and Tricks

In this section we cover several guidelines on how to use CHR to write constraint solvers and how to do so efficiently.

- *Check guard bindings yourself*  
  It is considered bad practice to write guards that bind variables of the head and to rely on the system to detect this at runtime. It is inefficient and obscures the working of the program.
- *Set semantics*  
  The CHR system allows the presence of identical constraints, i.e. multiple constraints with the same functor, arity and arguments. For most constraint solvers, this is not desirable: it affects efficiency and possibly termination. Hence appropriate simpagation rules should be added of the form: \[ constraint `\`constraint \<=\> true \]
- *Multi-headed rules*  
  Multi-headed rules are executed more efficiently when the constraints share one or more variables.
- *Mode and type declarations*  
  Provide mode and type declarations to get more efficient program execution. Make sure to disable debug (**--no-debug**) and enable optimization (**-O**).
- *Compile once, run many times*  
  Does consulting your CHR program take a long time in SWI-Prolog? Probably it takes the CHR compiler a long time to compile the CHR rules into Prolog code. When you disable optimizations the CHR compiler will be a lot quicker, but you may lose performance. Alternatively, you can just use SWI-Prolog's [qcompile/1](consulting.html#qcompile/1) to generate a `.qlf` file once from your `.pl` file. This `.qlf` contains the generated code of the CHR compiler (be it in a binary format). When you consult the `.qlf` file, the CHR compiler is not invoked and consultation is much faster.
- *Finding Constraints*  
  The [find_chr_constraint/1](chr-debugging.html#find_chr_constraint/1) predicate is fairly expensive. Avoid it, if possible. If you must use it, try to use it with an instantiated top-level constraint symbol.
