
## 3.8 Accessing the IDE from your program

Over the years a collection of IDE components have been developed, each with its own interface. In addition, some of these components require each other, and loading IDE components must be on demand to avoid the IDE being part of a saved state (see [qsave_program/2](saved-states.html#qsave_program/2)). For this reason, access to the IDE is concentrated on a single interface called [prolog_ide/1](idepreds.html#prolog_ide/1):

**prolog_ide**(`+Action`)  
This predicate ensures the IDE-enabling XPCE component is loaded, creates the XPCE class *prolog_ide* and sends `Action` to its one and only instance `@prolog_ide`. `Action` is one of the following:

**open_navigator**(`+Directory`)  
Open the Prolog Navigator (see [section 3.6](navigator.html#sec:3.6)) in the given `Directory`.

**open_debug_status**  
Open a window to edit spy and trace points.

**open_query_window**  
Open a little window to run Prolog queries from a GUI component.

**thread_monitor**  
Open a graphical window indicating existing threads and their status.

**debug_monitor**  
Open a graphical front-end for the `library(debug)` library that provides an overview of the topics and catches messages.

**xref**  
Open a graphical front-end for the cross-referencer that provides an overview of predicates and their callers.
