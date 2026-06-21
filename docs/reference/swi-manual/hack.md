
# B Hackers corner

This appendix describes a number of predicates which enable the Prolog user to inspect the Prolog environment and manipulate (or even redefine) the debugger. They can be used as entry points for experiments with debugging tools for Prolog. The predicates described here should be handled with some care as it is easy to corrupt the consistency of the Prolog system by misusing them.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[B.1 Examining the Environment Stack](manipstack.html)

[B.2 Ancestral cuts](ancestral-cut.html)

[B.3 Intercepting the Tracer](tracehook.html)

[B.4 Simulating a debugger interrupt](interrupt.html)

[B.5 Breakpoint and watchpoint handling](breakpoint.html)

[B.6 Adding context to errors: prolog_exception_hook](excepthook.html)

[B.7 Hooks using the exception predicate](exception3.html)

[B.8 Prolog events](prolog-event.html)

[B.9 Hooks for integrating libraries](intlibs.html)

[B.10 Hooks for loading files](loadfilehook.html)
