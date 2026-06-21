
# 4 Built-in Predicates

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[4.1 Notation of Predicate Descriptions](preddesc.html)

[4.1.1 The argument mode indicator](preddesc.html#sec:4.1.1)

[4.1.2 Predicate indicators](preddesc.html#sec:4.1.2)

[4.1.2.1 Non-terminal indicators](preddesc.html#sec:4.1.2.1)

[4.1.3 Predicate behaviour and determinism](preddesc.html#sec:4.1.3)

[4.2 Character representation](chars.html)

[4.3 Loading Prolog source files](consulting.html)

[4.3.1 Conditional compilation and program transformation](consulting.html#sec:4.3.1)

[4.3.1.1 Program transformation with source layout info](consulting.html#sec:4.3.1.1)

[4.3.1.2 Conditional compilation](consulting.html#sec:4.3.1.2)

[4.3.2 Reloading files, active code and threads](consulting.html#sec:4.3.2)

[4.3.2.1 Errors and warnings during compilation](consulting.html#sec:4.3.2.1)

[4.3.2.2 Compilation of mutually dependent code](consulting.html#sec:4.3.2.2)

[4.3.2.3 Compilation with multiple threads](consulting.html#sec:4.3.2.3)

[4.3.3 Quick load files](consulting.html#sec:4.3.3)

[4.4 Editor Interface](edit.html)

[4.4.1 Customizing the editor interface](edit.html#sec:4.4.1)

[4.5 Verify Type of a Term](typetest.html)

[4.6 Comparison and Unification of Terms](compare.html)

[4.6.1 Standard Order of Terms](compare.html#sec:4.6.1)

[4.6.2 Special unification and comparison predicates](compare.html#sec:4.6.2)

[4.7 Control Predicates](control.html)

[4.8 Meta-Call Predicates](metacall.html)

[4.9 Delimited continuations](delcont.html)

[4.10 Exception handling](exception.html)

[4.10.1 Unwind exceptions](exception.html#sec:4.10.1)

[4.10.2 Urgency of exceptions](exception.html#sec:4.10.2)

[4.10.3 Debugging and exceptions](exception.html#sec:4.10.3)

[4.10.4 The exception term](exception.html#sec:4.10.4)

[4.10.4.1 General form of the ISO standard exception term](exception.html#sec:4.10.4.1)

[4.10.4.2 Throwing exceptions from applications and libraries](exception.html#sec:4.10.4.2)

[4.11 Printing messages](printmsg.html)

[4.11.1 Printing from libraries](printmsg.html#sec:4.11.1)

[4.12 Handling signals](signal.html)

[4.12.1 Notes on signal handling](signal.html#sec:4.12.1)

[4.13 DCG Grammar rules](DCG.html)

[4.14 Database](db.html)

[4.14.1 Managing (dynamic) predicates](db.html#sec:4.14.1)

[4.14.1.1 Update view](db.html#sec:4.14.1.1)

[4.14.1.2 Indexing databases](db.html#sec:4.14.1.2)

[4.14.1.3 Transactions](db.html#sec:4.14.1.3)

[4.14.1.4 Impact of transactions](db.html#sec:4.14.1.4)

[4.14.2 The recorded database](db.html#sec:4.14.2)

[4.14.3 Flags](db.html#sec:4.14.3)

[4.14.4 Tries](db.html#sec:4.14.4)

[4.15 Declaring predicate properties](dynamic.html)

[4.16 Examining the program](examineprog.html)

[4.17 Input and output](IO.html)

[4.17.1 Predefined stream aliases](IO.html#sec:4.17.1)

[4.17.2 ISO Input and Output Streams](IO.html#sec:4.17.2)

[4.17.3 Edinburgh-style I/O](IO.html#sec:4.17.3)

[4.17.4 Switching between Edinburgh and ISO I/O](IO.html#sec:4.17.4)

[4.17.5 Adding IRI schemas](IO.html#sec:4.17.5)

[4.17.6 Write onto atoms, code-lists, etc.](IO.html#sec:4.17.6)

[4.17.7 Fast binary term I/O](IO.html#sec:4.17.7)

[4.18 Status of streams](streamstat.html)

[4.19 Primitive character I/O](chario.html)

[4.20 Term reading and writing](termrw.html)

[4.21 Analysing and Constructing Terms](manipterm.html)

[4.21.1 Non-logical operations on terms](manipterm.html#sec:4.21.1)

[4.22 Analysing and Constructing Atoms](manipatom.html)

[4.23 Localization (locale) support](locale.html)

[4.24 Character properties](chartype.html)

[4.24.1 Case conversion](chartype.html#sec:4.24.1)

[4.24.2 White space normalization](chartype.html#sec:4.24.2)

[4.24.3 Language-specific comparison](chartype.html#sec:4.24.3)

[4.25 Operators](operators.html)

[4.26 Character Conversion](charconv.html)

[4.27 Arithmetic](arith.html)

[4.27.1 Special purpose integer arithmetic](arith.html#sec:4.27.1)

[4.27.2 General purpose arithmetic](arith.html#sec:4.27.2)

[4.27.2.1 Arithmetic types](arith.html#sec:4.27.2.1)

[4.27.2.2 Rational number examples](arith.html#sec:4.27.2.2)

[4.27.2.3 Rational numbers or floats](arith.html#sec:4.27.2.3)

[4.27.2.4 IEEE 754 floating point arithmetic](arith.html#sec:4.27.2.4)

[4.27.2.5 Floating point arithmetic precision](arith.html#sec:4.27.2.5)

[4.27.2.6 Arithmetic Functions](arith.html#sec:4.27.2.6)

[4.28 Misc arithmetic support predicates](miscarith.html)

[4.29 Built-in list operations](builtinlist.html)

[4.30 Finding all Solutions to a Goal](allsolutions.html)

[4.31 Forall](forall2.html)

[4.32 Formatted Write](format.html)

[4.32.1 Programming Format](format.html#sec:4.32.1)

[4.33 Global variables](gvar.html)

[4.33.1 Compatibility of SWI-Prolog Global Variables](gvar.html#sec:4.33.1)

[4.34 Terminal Control](tty.html)

[4.35 Operating System Interaction](system.html)

[4.35.1 Windows-specific Operating System Interaction](system.html#sec:4.35.1)

[4.35.2 Apple specific Operating System Interaction](system.html#sec:4.35.2)

[4.35.3 Dealing with time and date](system.html#sec:4.35.3)

[4.35.3.1 Time and date data structures](system.html#sec:4.35.3.1)

[4.35.3.2 Time and date predicates](system.html#sec:4.35.3.2)

[4.35.4 Controlling the **swipl-win** (Epilog) console window](system.html#sec:4.35.4)

[4.36 File System Interaction](files.html)

[4.37 User Top-level Manipulation](toplevel.html)

[4.38 Creating a Protocol of the User Interaction](protocol.html)

[4.39 Debugging and Tracing Programs](debugger.html)

[4.40 Debugging and declaring determinism](debug-determinism.html)

[4.41 Obtaining Runtime Statistics](builtin-statistics.html)

[4.42 Execution profiling](profile.html)

[4.42.1 library(prolog_profile): Execution profiler](profile.html#sec:4.42.1)

[4.42.2 Visualizing profiling data](profile.html#sec:4.42.2)

[4.42.3 Information gathering](profile.html#sec:4.42.3)

[4.42.3.1 Profiling in the Windows Implementation](profile.html#sec:4.42.3.1)

[4.43 Memory Management](memory.html)

[4.43.1 Garbage collection](memory.html#sec:4.43.1)

[4.43.2 Heap memory (malloc)](memory.html#sec:4.43.2)

[4.43.2.1 TCMalloc control predicates](memory.html#sec:4.43.2.1)

[4.44 Windows DDE interface](DDE.html)

[4.44.1 DDE client interface](DDE.html#sec:4.44.1)

[4.44.2 DDE server mode](DDE.html#sec:4.44.2)

[4.45 Miscellaneous](miscpreds.html)
