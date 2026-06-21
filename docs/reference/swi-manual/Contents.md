
# Table of Contents

[1 Introduction](intro.html)

[1.1 Positioning SWI-Prolog](swiprolog.html)

[1.2 Status and releases](status.html)

[1.3 Should I be using SWI-Prolog?](swiorother.html)

[1.4 Support the SWI-Prolog project](sponsor.html)

[1.5 Implementation history](implhistory.html)

[1.6 Acknowledgements](acknowledge.html)

[2 Overview](overview.html)

[2.1 Getting started quickly](quickstart.html)

[2.1.1 Starting SWI-Prolog](quickstart.html#sec:2.1.1)

[2.1.1.1 Starting SWI-Prolog on Unix](quickstart.html#sec:2.1.1.1)

[2.1.1.2 Starting SWI-Prolog on Windows](quickstart.html#sec:2.1.1.2)

[2.1.2 Adding rules from the console](quickstart.html#sec:2.1.2)

[2.1.3 Executing a query](quickstart.html#sec:2.1.3)

[2.1.4 Examining and modifying your program](quickstart.html#sec:2.1.4)

[2.1.5 Stopping Prolog](quickstart.html#sec:2.1.5)

[2.2 The user's initialisation file](initfile.html)

[2.3 Initialisation files and goals](initgoal.html)

[2.4 Command line options](cmdline.html)

[2.4.1 Informational command line options](cmdline.html#sec:2.4.1)

[2.4.2 Command line options for running Prolog](cmdline.html#sec:2.4.2)

[2.4.3 Controlling the stack sizes](cmdline.html#sec:2.4.3)

[2.4.4 Running goals from the command line](cmdline.html#sec:2.4.4)

[2.4.5 Compilation options](cmdline.html#sec:2.4.5)

[2.4.6 Maintenance options](cmdline.html#sec:2.4.6)

[2.5 UI Themes](theme.html)

[2.5.1 Status of theme support](theme.html#sec:2.5.1)

[2.6 GNU Emacs Interface](gemacs.html)

[2.7 Online Help](online-help.html)

[2.7.1 library(help): Text based manual](online-help.html#sec:2.7.1)

[2.7.2 library(explain): Describe Prolog Terms](online-help.html#sec:2.7.2)

[2.8 Command line history](history.html)

[2.9 Reuse of top-level bindings](topvars.html)

[2.10 Overview of the Debugger](debugoverview.html)

[2.10.1 The Byrd Box Model And Ports](debugoverview.html#sec:2.10.1)

[2.10.2 Trace Mode Example](debugoverview.html#sec:2.10.2)

[2.10.3 Trace Mode Options: leash/1 and visible/1](debugoverview.html#sec:2.10.3)

[2.10.4 Trace Mode Commands When Paused](debugoverview.html#sec:2.10.4)

[2.10.4.1 Control Flow Commands](debugoverview.html#sec:2.10.4.1)

[2.10.4.2 Informational Commands](debugoverview.html#sec:2.10.4.2)

[2.10.4.3 Formatting Commands](debugoverview.html#sec:2.10.4.3)

[2.10.5 Trace Mode vs. Trace Point](debugoverview.html#sec:2.10.5)

[2.10.6 Spy Points and Debug Mode](debugoverview.html#sec:2.10.6)

[2.10.7 Breakpoints](debugoverview.html#sec:2.10.7)

[2.10.8 Command Line Debugger Summary](debugoverview.html#sec:2.10.8)

[2.10.8.1 Trace Mode](debugoverview.html#sec:2.10.8.1)

[2.10.8.2 Trace Points](debugoverview.html#sec:2.10.8.2)

[2.11 Loading and running projects](compilation.html)

[2.11.1 Running an application](compilation.html#sec:2.11.1)

[2.11.1.1 Using PrologScript](compilation.html#sec:2.11.1.1)

[2.11.1.2 Creating a shell script](compilation.html#sec:2.11.1.2)

[2.11.1.3 Creating a saved state](compilation.html#sec:2.11.1.3)

[2.11.1.4 Compilation using the -c command line option](compilation.html#sec:2.11.1.4)

[2.11.1.5 SWI-Prolog app scripts](compilation.html#sec:2.11.1.5)

[2.12 Environment Control (Prolog flags)](flags.html)

[2.13 An overview of hook predicates](hooks.html)

[2.14 Automatic loading of libraries](autoload.html)

[2.15 The SWI-Prolog syntax](syntax.html)

[2.15.1 ISO Syntax Support](syntax.html#sec:2.15.1)

[2.15.1.1 Processor Character Set](syntax.html#sec:2.15.1.1)

[2.15.1.2 Nested comments](syntax.html#sec:2.15.1.2)

[2.15.1.3 Character Escape Syntax](syntax.html#sec:2.15.1.3)

[2.15.1.4 Syntax for non-decimal numbers](syntax.html#sec:2.15.1.4)

[2.15.1.5 Using digit groups in large integers](syntax.html#sec:2.15.1.5)

[2.15.1.6 Rational number syntax](syntax.html#sec:2.15.1.6)

[2.15.1.7 NaN and Infinity floats and their syntax](syntax.html#sec:2.15.1.7)

[2.15.1.8 Force only underscore to introduce a variable](syntax.html#sec:2.15.1.8)

[2.15.1.9 Unicode Prolog source](syntax.html#sec:2.15.1.9)

[2.15.1.10 Singleton variable checking](syntax.html#sec:2.15.1.10)

[2.16 Rational trees (cyclic terms)](cyclic.html)

[2.17 Just-in-time clause indexing](jitindex.html)

[2.17.1 Deep indexing](jitindex.html#sec:2.17.1)

[2.17.2 Future directions](jitindex.html#sec:2.17.2)

[2.17.3 Indexing for body code](jitindex.html#sec:2.17.3)

[2.17.4 Indexing and portability](jitindex.html#sec:2.17.4)

[2.18 Wide character support](widechars.html)

[2.18.1 Wide character encodings on streams](widechars.html#sec:2.18.1)

[2.18.1.1 BOM: Byte Order Mark](widechars.html#sec:2.18.1.1)

[2.19 System limits](limits.html)

[2.19.1 Limits on memory areas](limits.html#sec:2.19.1)

[2.19.1.1 The heap](limits.html#sec:2.19.1.1)

[2.19.2 Other Limits](limits.html#sec:2.19.2)

[2.19.3 Reserved Names](limits.html#sec:2.19.3)

[2.20 SWI-Prolog and 32-bit machines](32bits.html)

[2.21 Binary compatibility](abi-versions.html)

[3 Initialising and Managing a Prolog Project](IDE.html)

[3.1 The project source files](projectfiles.html)

[3.1.1 File Names and Locations](projectfiles.html#sec:3.1.1)

[3.1.1.1 File Name Extensions](projectfiles.html#sec:3.1.1.1)

[3.1.1.2 Project Directories](projectfiles.html#sec:3.1.1.2)

[3.1.1.3 Sub-projects using search paths](projectfiles.html#sec:3.1.1.3)

[3.1.2 Project Special Files](projectfiles.html#sec:3.1.2)

[3.1.3 International source files](projectfiles.html#sec:3.1.3)

[3.2 Using modules](usingmodules.html)

[3.3 The test-edit-reload cycle](editreload.html)

[3.3.1 Locating things to edit](editreload.html#sec:3.3.1)

[3.3.2 Editing and incremental compilation](editreload.html#sec:3.3.2)

[3.4 Using the PceEmacs built-in editor](pceemacs.html)

[3.4.1 Activating PceEmacs](pceemacs.html#sec:3.4.1)

[3.4.2 Bluffing through PceEmacs](pceemacs.html#sec:3.4.2)

[3.4.2.1 Edit modes](pceemacs.html#sec:3.4.2.1)

[3.4.2.2 Frequently used editor commands](pceemacs.html#sec:3.4.2.2)

[3.4.3 Prolog Mode](pceemacs.html#sec:3.4.3)

[3.4.3.1 Finding your way around](pceemacs.html#sec:3.4.3.1)

[3.5 The Graphical Debugger](guitracer.html)

[3.5.1 Invoking the window-based debugger](guitracer.html#sec:3.5.1)

[3.6 The Prolog Navigator](navigator.html)

[3.7 Cross-referencer](xref.html)

[3.8 Accessing the IDE from your program](idepreds.html)

[3.9 Summary of the IDE](idesummary.html)

[4 Built-in Predicates](builtin.html)

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

[5 SWI-Prolog extensions](extensions.html)

[5.1 Lists are special](ext-lists.html)

[5.1.1 Motivating‘`[|]`’and `[]` for lists](ext-lists.html#sec:5.1.1)

[5.2 The string type and its double quoted syntax](string.html)

[5.2.1 Representing text: strings, atoms and code lists](string.html#sec:5.2.1)

[5.2.2 Predicates that operate on strings](string.html#sec:5.2.2)

[5.2.3 Why has the representation of double quoted text changed?](string.html#sec:5.2.3)

[5.2.4 Adapting code for double quoted strings](string.html#sec:5.2.4)

[5.2.5 Predicates to support adapting code for double quoted strings](string.html#sec:5.2.5)

[5.3 Syntax changes since SWI-Prolog 7](ext-syntax.html)

[5.3.1 Operators and quoted atoms](ext-syntax.html#sec:5.3.1)

[5.3.2 Compound terms with zero arguments](ext-syntax.html#sec:5.3.2)

[5.3.3 Block operators](ext-syntax.html#sec:5.3.3)

[5.4 Dicts: structures with named arguments](bidicts.html)

[5.4.1 Dict types and tags](bidicts.html#sec:5.4.1)

[5.4.2 Functions on dicts](bidicts.html#sec:5.4.2)

[5.4.2.1 User defined functions on dicts](bidicts.html#sec:5.4.2.1)

[5.4.2.2 Predefined functions on dicts](bidicts.html#sec:5.4.2.2)

[5.4.3 Predicates for managing dicts](bidicts.html#sec:5.4.3)

[5.4.3.1 Destructive assignment in dicts](bidicts.html#sec:5.4.3.1)

[5.4.4 When to use dicts?](bidicts.html#sec:5.4.4)

[5.4.5 A motivation for dicts as primary citizens](bidicts.html#sec:5.4.5)

[5.4.6 Implementation notes about dicts](bidicts.html#sec:5.4.6)

[5.5 Integration of strings and dicts in the libraries](ext-integration.html)

[5.5.1 Dicts and option processing](ext-integration.html#sec:5.5.1)

[5.5.2 Dicts in core data structures](ext-integration.html#sec:5.5.2)

[5.5.3 Dicts, strings and XML](ext-integration.html#sec:5.5.3)

[5.5.4 Dicts, strings and JSON](ext-integration.html#sec:5.5.4)

[5.5.5 Dicts, strings and HTTP](ext-integration.html#sec:5.5.5)

[5.6 Single Sided Unification rules](ssu.html)

[5.6.1 Single Sided Unification Guards](ssu.html#sec:5.6.1)

[5.6.2 Consequences of `=>` single sided unification rules](ssu.html#sec:5.6.2)

[5.6.3 Single sided unification for Definite Clause Grammars](ssu.html#sec:5.6.3)

[5.6.4 SSU: Future considerations](ssu.html#sec:5.6.4)

[5.7 Remaining issues](ext-issues.html)

[6 Modules](modules.html)

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

[7 Tabled execution (SLG resolution)](tabling.html)

[7.1 Example 1: using tabling for memoizing](tabling-memoize.html)

[7.2 Example 2: avoiding non-termination](tabling-non-termination.html)

[7.3 Answer subsumption or mode directed tabling](tabling-mode-directed.html)

[7.4 Tabling for impure programs](tnotpure.html)

[7.5 Variant and subsumptive tabling](tabling-subsumptive.html)

[7.6 Well Founded Semantics](WFS.html)

[7.6.1 Well founded semantics and the toplevel](WFS.html#sec:7.6.1)

[7.7 Incremental tabling](tabling-incremental.html)

[7.8 Monotonic tabling](tabling-monotonic.html)

[7.8.1 Eager and lazy monotonic tabling](tabling-monotonic.html#sec:7.8.1)

[7.8.2 Tracking new answers to monotonic tables](tabling-monotonic.html#sec:7.8.2)

[7.8.3 Monotonic tabling with external data](tabling-monotonic.html#sec:7.8.3)

[7.9 Shared tabling](tabling-shared.html)

[7.9.1 Abolishing shared tables](tabling-shared.html#sec:7.9.1)

[7.9.2 Status and future of shared tabling](tabling-shared.html#sec:7.9.2)

[7.10 Tabling and constraints](tabling-constraints.html)

[7.11 Tabling restraints: bounded rationality and tripwires](tabling-restraints.html)

[7.11.1 Restraint subgoal size](tabling-restraints.html#sec:7.11.1)

[7.11.2 Restraint answer size](tabling-restraints.html#sec:7.11.2)

[7.11.3 Restraint answer count](tabling-restraints.html#sec:7.11.3)

[7.12 Tabling predicate reference](tabling-preds.html)

[7.13 About the tabling implementation](tabling-about.html)

[7.13.1 Status of tabling](tabling-about.html#sec:7.13.1)

[8 Constraint Logic Programming](clp.html)

[8.1 Attributed variables](attvar.html)

[8.1.1 Attribute manipulation predicates](attvar.html#sec:8.1.1)

[8.1.2 Attributed variable hooks](attvar.html#sec:8.1.2)

[8.1.3 Operations on terms with attributed variables](attvar.html#sec:8.1.3)

[8.1.4 Special purpose predicates for attributes](attvar.html#sec:8.1.4)

[8.2 Coroutining](coroutining.html)

[9 CHR: Constraint Handling Rules](chr.html)

[9.1 Introduction to CHR](chr-intro.html)

[9.2 CHR Syntax and Semantics](chr-syntaxandsemantics.html)

[9.2.1 Syntax of CHR rules](chr-syntaxandsemantics.html#sec:9.2.1)

[9.2.2 Semantics of CHR](chr-syntaxandsemantics.html#sec:9.2.2)

[9.3 CHR in SWI-Prolog Programs](practical.html)

[9.3.1 Embedding CHR in Prolog Programs](practical.html#sec:9.3.1)

[9.3.2 CHR Constraint declaration](practical.html#sec:9.3.2)

[9.3.3 CHR Compilation](practical.html#sec:9.3.3)

[9.4 Debugging CHR programs](chr-debugging.html)

[9.4.1 CHR debug ports](chr-debugging.html#sec:9.4.1)

[9.4.2 Tracing CHR programs](chr-debugging.html#sec:9.4.2)

[9.4.3 CHR Debugging Predicates](chr-debugging.html#sec:9.4.3)

[9.5 CHR Examples](chr-examples.html)

[9.6 CHR compatibility](chr-compatibility.html)

[9.6.1 The Old SICStus CHR implementation](chr-compatibility.html#sec:9.6.1)

[9.6.2 The Old ECLiPSe CHR implementation](chr-compatibility.html#sec:9.6.2)

[9.7 CHR Programming Tips and Tricks](chr-guidelines.html)

[9.8 CHR Compiler Errors and Warnings](chr-warnings-and-errors.html)

[9.8.1 CHR Compiler Errors](chr-warnings-and-errors.html#sec:9.8.1)

[10 Multithreaded applications](threads.html)

[10.1 Creating and destroying Prolog threads](threadcreate.html)

[10.2 Monitoring threads](thmonitor.html)

[10.3 Thread communication](threadcom.html)

[10.3.1 Message queues](threadcom.html#sec:10.3.1)

[10.3.2 Waiting for events](threadcom.html#sec:10.3.2)

[10.3.3 Signalling threads](threadcom.html#sec:10.3.3)

[10.3.4 Threads and dynamic predicates](threadcom.html#sec:10.3.4)

[10.4 Thread synchronisation](threadsync.html)

[10.5 library(threadutil): Interactive thread utilities](threadutil.html)

[10.6 Multithreaded mixed C and Prolog applications](foreignthread.html)

[10.6.1 A Prolog thread for each native thread (one-to-one)](foreignthread.html#sec:10.6.1)

[10.6.2 Using Prolog engines from C](foreignthread.html#sec:10.6.2)

[10.7 Multithreading and the XPCE graphics system](mt-xpce.html)

[11 Coroutining using Prolog engines](engines.html)

[11.1 Examples using engines](engine-examples.html)

[11.1.1 Aggregation using engines](engine-examples.html#sec:11.1.1)

[11.1.2 State accumulation using engines](engine-examples.html#sec:11.1.2)

[11.1.3 Scalable many-agent applications](engine-examples.html#sec:11.1.3)

[11.2 Engine resource usage](engine-resources.html)

[11.3 Engine predicate reference](engine-predicates.html)

[12 Foreign Language Interface](foreign.html)

[12.1 Overview of the Interface](foreignoverview.html)

[12.2 Linking Foreign Modules](foreignlink.html)

[12.2.1 What linking is provided?](foreignlink.html#sec:12.2.1)

[12.2.2 What kind of loading should I be using?](foreignlink.html#sec:12.2.2)

[12.2.3 library(shlib): Utility library for loading foreign objects (DLLs, shared objects)](foreignlink.html#sec:12.2.3)

[12.2.4 Low-level operations on shared libraries](foreignlink.html#sec:12.2.4)

[12.2.5 Static Linking](foreignlink.html#sec:12.2.5)

[12.3 Interface Data Types](foreigntypes.html)

[12.3.1 Type `term_t`: a reference to a Prolog term](foreigntypes.html#sec:12.3.1)

[12.3.1.1 Interaction with the garbage collector and stack-shifter](foreigntypes.html#sec:12.3.1.1)

[12.3.2 Other foreign interface types](foreigntypes.html#sec:12.3.2)

[12.3.2.1 PL_ARITY_AS_SIZE](foreigntypes.html#sec:12.3.2.1)

[12.3.2.2 Notes on C API bool return values](foreigntypes.html#sec:12.3.2.2)

[12.4 The Foreign Include File](foreigninclude.html)

[12.4.1 Argument Passing and Control](foreigninclude.html#sec:12.4.1)

[12.4.1.1 Non-deterministic Foreign Predicates](foreigninclude.html#sec:12.4.1.1)

[12.4.1.2 Yielding from foreign predicates](foreigninclude.html#sec:12.4.1.2)

[12.4.1.3 Implementing a yield based debugger](foreigninclude.html#sec:12.4.1.3)

[12.4.2 Atoms and functors](foreigninclude.html#sec:12.4.2)

[12.4.2.1 Atoms and atom garbage collection](foreigninclude.html#sec:12.4.2.1)

[12.4.3 Input and output](foreigninclude.html#sec:12.4.3)

[12.4.4 Analysing Terms via the Foreign Interface](foreigninclude.html#sec:12.4.4)

[12.4.4.1 Testing the type of a term](foreigninclude.html#sec:12.4.4.1)

[12.4.4.2 Reading data from a term](foreigninclude.html#sec:12.4.4.2)

[12.4.4.3 Exchanging text using length and string](foreigninclude.html#sec:12.4.4.3)

[12.4.4.4 Wide-character versions](foreigninclude.html#sec:12.4.4.4)

[12.4.4.5 Reading a list](foreigninclude.html#sec:12.4.4.5)

[12.4.4.6 Processing option lists and dicts](foreigninclude.html#sec:12.4.4.6)

[12.4.4.7 An example: defining write/1 in C](foreigninclude.html#sec:12.4.4.7)

[12.4.5 Constructing Terms](foreigninclude.html#sec:12.4.5)

[12.4.6 Unifying data](foreigninclude.html#sec:12.4.6)

[12.4.7 Convenient functions to generate Prolog exceptions](foreigninclude.html#sec:12.4.7)

[12.4.8 Foreign language wrapper support functions](foreigninclude.html#sec:12.4.8)

[12.4.9 Serializing and deserializing Prolog terms](foreigninclude.html#sec:12.4.9)

[12.4.10 BLOBS: Using atoms to store arbitrary binary data](foreigninclude.html#sec:12.4.10)

[12.4.10.1 Defining a BLOB type](foreigninclude.html#sec:12.4.10.1)

[12.4.10.2 Accessing blobs](foreigninclude.html#sec:12.4.10.2)

[12.4.10.3 Considerations for non-C code](foreigninclude.html#sec:12.4.10.3)

[12.4.11 Exchanging GMP numbers](foreigninclude.html#sec:12.4.11)

[12.4.12 Calling Prolog from C](foreigninclude.html#sec:12.4.12)

[12.4.12.1 Predicate references](foreigninclude.html#sec:12.4.12.1)

[12.4.12.2 Initiating a query from C](foreigninclude.html#sec:12.4.12.2)

[12.4.13 Discarding Data](foreigninclude.html#sec:12.4.13)

[12.4.14 String buffering](foreigninclude.html#sec:12.4.14)

[12.4.15 Foreign Code and Modules](foreigninclude.html#sec:12.4.15)

[12.4.16 Prolog exceptions in foreign code](foreigninclude.html#sec:12.4.16)

[12.4.17 Catching Signals (Software Interrupts)](foreigninclude.html#sec:12.4.17)

[12.4.18 Miscellaneous](foreigninclude.html#sec:12.4.18)

[12.4.18.1 Term Comparison](foreigninclude.html#sec:12.4.18.1)

[12.4.18.2 Recorded database](foreigninclude.html#sec:12.4.18.2)

[12.4.18.3 Database](foreigninclude.html#sec:12.4.18.3)

[12.4.18.4 Getting file names](foreigninclude.html#sec:12.4.18.4)

[12.4.18.5 Dealing with Prolog flags from C](foreigninclude.html#sec:12.4.18.5)

[12.4.18.6 Foreign code and Well Founded Semantics](foreigninclude.html#sec:12.4.18.6)

[12.4.19 Errors and warnings](foreigninclude.html#sec:12.4.19)

[12.4.20 Environment Control from Foreign Code](foreigninclude.html#sec:12.4.20)

[12.4.21 Querying Prolog](foreigninclude.html#sec:12.4.21)

[12.4.22 Registering Foreign Predicates](foreigninclude.html#sec:12.4.22)

[12.4.23 Foreign Code Hooks](foreigninclude.html#sec:12.4.23)

[12.4.24 Storing foreign data](foreigninclude.html#sec:12.4.24)

[12.4.24.1 Examples for storing foreign data](foreigninclude.html#sec:12.4.24.1)

[12.4.25 Embedding SWI-Prolog in other applications](foreigninclude.html#sec:12.4.25)

[12.4.25.1 Threading, Signals and embedded Prolog](foreigninclude.html#sec:12.4.25.1)

[12.5 Linking embedded applications using swipl-ld](plld.html)

[12.5.1 A simple example](plld.html#sec:12.5.1)

[12.6 The Prolog‘home’directory](findhome.html)

[12.7 Example of Using the Foreign Interface](foreignxmp.html)

[12.8 Notes on Using Foreign Code](foreignnotes.html)

[12.8.1 Foreign debugging functions](foreignnotes.html#sec:12.8.1)

[12.8.2 Memory Allocation](foreignnotes.html#sec:12.8.2)

[12.8.3 Compatibility between Prolog versions](foreignnotes.html#sec:12.8.3)

[12.8.4 Foreign hash tables](foreignnotes.html#sec:12.8.4)

[12.8.5 Debugging and profiling foreign code (valgrind, asan)](foreignnotes.html#sec:12.8.5)

[12.8.6 Name Conflicts in C modules](foreignnotes.html#sec:12.8.6)

[12.8.7 Compatibility of the Foreign Interface](foreignnotes.html#sec:12.8.7)

[12.9 Foreign access to Prolog IO streams](foreign-streams.html)

[12.9.1 Get IO stream handles](foreign-streams.html#sec:12.9.1)

[12.9.2 Creating an IO stream](foreign-streams.html#sec:12.9.2)

[12.9.3 Interacting with foreign streams](foreign-streams.html#sec:12.9.3)

[12.9.3.1 Writing Prolog terms to foreign streams](foreign-streams.html#sec:12.9.3.1)

[12.9.4 Foreign stream error handling](foreign-streams.html#sec:12.9.4)

[12.9.5 Foreign stream encoding](foreign-streams.html#sec:12.9.5)

[12.9.6 Foreign stream line endings](foreign-streams.html#sec:12.9.6)

[12.9.7 Foreign stream position information](foreign-streams.html#sec:12.9.7)

[12.9.8 Support functions for blob save/load](foreign-streams.html#sec:12.9.8)

[13 Using SWI-Prolog in your browser (WASM)](wasm-version.html)

[13.1 Loading and initializing Prolog](wasm-loading.html)

[13.1.1 Loading Prolog files](wasm-loading.html#sec:13.1.1)

[13.2 Calling Prolog from JavaScript](wasm-calling.html)

[13.2.1 The JavaScript class Query](wasm-calling.html#sec:13.2.1)

[13.2.2 Using engines](wasm-calling.html#sec:13.2.2)

[13.2.3 Translating data between JavaScript and Prolog](wasm-calling.html#sec:13.2.3)

[13.2.3.1 Translating JavaScript data to Prolog](wasm-calling.html#sec:13.2.3.1)

[13.2.3.2 Translating Prolog data to JavaScript](wasm-calling.html#sec:13.2.3.2)

[13.3 Accessing JavaScript from Prolog](wasm-js-call.html)

[13.3.1 Asynchronous access to JavaScript from Prolog](wasm-js-call.html#sec:13.3.1)

[13.3.2 JavaScript Promise that can be aborted](wasm-js-call.html#sec:13.3.2)

[14 Deploying applications](runtime.html)

[14.1 Deployment options](deployment-options.html)

[14.2 Understanding saved states](saved-states.html)

[14.2.1 Creating a saved state](saved-states.html#sec:14.2.1)

[14.2.2 Limitations of qsave_program](saved-states.html#sec:14.2.2)

[14.2.3 Runtimes and Foreign Code](saved-states.html#sec:14.2.3)

[14.3 State initialization](state-initialization.html)

[14.4 Using program resources](program-resources.html)

[14.4.1 Resources as files](program-resources.html#sec:14.4.1)

[14.4.2 Access resources using open_resource](program-resources.html#sec:14.4.2)

[14.4.3 Declaring resources](program-resources.html#sec:14.4.3)

[14.4.4 Managing resource files](program-resources.html#sec:14.4.4)

[14.5 Debugging and updating deployed systems](debug-deployed-systems.html)

[14.6 Protecting your code](protect-code.html)

[14.6.1 Obfuscating code in saved states](protect-code.html#sec:14.6.1)

[14.7 Finding Application files](findappfile.html)

[15 Packs: community add-ons](packs.html)

[15.1 Installing packs](pack-install.html)

[15.2 Built-in predicates for attaching packs](pack-attach.html)

[15.3 library(prolog_pack): A package manager for Prolog](prologpack.html)

[15.4 Structure of a pack](pack-structure.html)

[15.5 Developing a pack](pack-devel.html)

[15.5.1 The pack meta data](pack-devel.html#sec:15.5.1)

[15.5.1.1 Pack requirements on Prolog](pack-devel.html#sec:15.5.1.1)

[15.5.2 Packs with foreign code](pack-devel.html#sec:15.5.2)

[15.5.2.1 Compiling a foreign extension using a simple Makefile](pack-devel.html#sec:15.5.2.1)

[15.5.2.2 Publishing a pack](pack-devel.html#sec:15.5.2.2)

[15.5.2.3 Compiling a foreign extension using CMake](pack-devel.html#sec:15.5.2.3)

[15.5.3 Updating a package](pack-devel.html#sec:15.5.3)

[A The SWI-Prolog library](libpl.html)

[A.1 library(aggregate): Aggregation operators on backtrackable predicates](aggregate.html)

[A.2 library(ansi_term): Print decorated text to ANSI consoles](ansiterm.html)

[A.3 library(apply): Apply predicates on a list](apply.html)

[A.4 library(assoc): Association lists](assoc.html)

[A.4.1 Introduction](assoc.html#sec:A.4.1)

[A.4.2 Creating association lists](assoc.html#sec:A.4.2)

[A.4.3 Querying association lists](assoc.html#sec:A.4.3)

[A.4.4 Modifying association lists](assoc.html#sec:A.4.4)

[A.4.5 Conversion predicates](assoc.html#sec:A.4.5)

[A.4.6 Reasoning about association lists and their elements](assoc.html#sec:A.4.6)

[A.5 library(broadcast): Broadcast and receive event notifications](broadcast.html)

[A.6 library(charsio): I/O on Lists of Character Codes](charsio.html)

[A.7 library(check): Consistency checking](check.html)

[A.8 library(clpb): CLP(B): Constraint Logic Programming over Boolean Variables](clpb.html)

[A.8.1 Introduction](clpb.html#sec:A.8.1)

[A.8.2 Boolean expressions](clpb.html#sec:A.8.2)

[A.8.3 Interface predicates](clpb.html#sec:A.8.3)

[A.8.4 Examples](clpb.html#sec:A.8.4)

[A.8.5 Obtaining BDDs](clpb.html#sec:A.8.5)

[A.8.6 Enabling monotonic CLP(B)](clpb.html#sec:A.8.6)

[A.8.7 Example: Pigeons](clpb.html#sec:A.8.7)

[A.8.8 Example: Boolean circuit](clpb.html#sec:A.8.8)

[A.8.9 Acknowledgments](clpb.html#sec:A.8.9)

[A.8.10 CLP(B) predicate index](clpb.html#sec:A.8.10)

[A.9 library(clpfd): CLP(FD): Constraint Logic Programming over Finite Domains](clpfd.html)

[A.9.1 Introduction](clpfd.html#sec:A.9.1)

[A.9.2 Arithmetic constraints](clpfd.html#sec:A.9.2)

[A.9.3 Declarative integer arithmetic](clpfd.html#sec:A.9.3)

[A.9.4 Example: Factorial relation](clpfd.html#sec:A.9.4)

[A.9.5 Combinatorial constraints](clpfd.html#sec:A.9.5)

[A.9.6 Domains](clpfd.html#sec:A.9.6)

[A.9.7 Example: Sudoku](clpfd.html#sec:A.9.7)

[A.9.8 Residual goals](clpfd.html#sec:A.9.8)

[A.9.9 Core relations and search](clpfd.html#sec:A.9.9)

[A.9.10 Example: Eight queens puzzle](clpfd.html#sec:A.9.10)

[A.9.11 Optimisation](clpfd.html#sec:A.9.11)

[A.9.12 Reification](clpfd.html#sec:A.9.12)

[A.9.13 Enabling monotonic CLP(FD)](clpfd.html#sec:A.9.13)

[A.9.14 Custom constraints](clpfd.html#sec:A.9.14)

[A.9.15 Applications](clpfd.html#sec:A.9.15)

[A.9.16 Acknowledgments](clpfd.html#sec:A.9.16)

[A.9.17 CLP(FD) predicate index](clpfd.html#sec:A.9.17)

[A.9.17.1 Arithmetic constraints](clpfd.html#sec:A.9.17.1)

[A.9.17.2 Membership constraints](clpfd.html#sec:A.9.17.2)

[A.9.17.3 Enumeration predicates](clpfd.html#sec:A.9.17.3)

[A.9.17.4 Global constraints](clpfd.html#sec:A.9.17.4)

[A.9.17.5 Reification predicates](clpfd.html#sec:A.9.17.5)

[A.9.17.6 Reflection predicates](clpfd.html#sec:A.9.17.6)

[A.9.17.7 FD set predicates](clpfd.html#sec:A.9.17.7)

[A.9.17.8 FD miscellaneous predicates](clpfd.html#sec:A.9.17.8)

[A.9.18 Closing and opening words about CLP(FD)](clpfd.html#sec:A.9.18)

[A.10 library(clpqr): Constraint Logic Programming over Rationals and Reals](clpqr.html)

[A.10.1 Solver predicates](clpqr.html#sec:A.10.1)

[A.10.2 Syntax of the predicate arguments](clpqr.html#sec:A.10.2)

[A.10.3 Use of unification](clpqr.html#sec:A.10.3)

[A.10.4 Non-linear constraints](clpqr.html#sec:A.10.4)

[A.10.5 Status and known problems](clpqr.html#sec:A.10.5)

[A.11 library(csv): Process CSV (Comma-Separated Values) data](csv.html)

[A.12 library(dcg/basics): Various general DCG utilities](basics.html)

[A.13 library(dcg/high_order): High order grammar operations](highorder.html)

[A.14 library(debug): Print debug messages and test assertions](debug.html)

[A.15 library(dicts): Dict utilities](dicts.html)

[A.16 library(error): Error generating support](error.html)

[A.17 library(exceptions): Exception classification](exceptions.html)

[A.18 library(fastrw): Fast reading and writing of terms](fastrw.html)

[A.19 library(gensym): Generate unique symbols](gensym.html)

[A.20 library(heaps): heaps/priority queues](heaps.html)

[A.21 library(increval): Incremental dynamic predicate modification](increval.html)

[A.22 library(intercept): Intercept and signal interface](intercept.html)

[A.23 library(iostream): Utilities to deal with streams](iostream.html)

[A.24 library(listing): List programs and pretty print clauses](listing.html)

[A.25 library(lists): List Manipulation](lists.html)

[A.26 library(macros): Macro expansion](macros.html)

[A.26.1 Defining and using macros](macros.html#sec:A.26.1)

[A.26.2 Implementation details](macros.html#sec:A.26.2)

[A.26.3 Predicates](macros.html#sec:A.26.3)

[A.27 library(main): Provide entry point for scripts](main.html)

[A.28 library(nb_set): Non-backtrackable set](nb_set.html)

[A.29 library(writef): Old-style formatted write](writef.html)

[A.30 library(www_browser): Open a URL in the users browser](wwwbrowser.html)

[A.31 library(occurs): Finding and counting sub-terms](occurs.html)

[A.32 library(option): Option list processing](option.html)

[A.33 library(optparse): command line parsing](optparse.html)

[A.33.1 Notes and tips](optparse.html#sec:A.33.1)

[A.34 library(ordsets): Ordered set manipulation](ordsets.html)

[A.35 library(pairs): Operations on key-value lists](pairs.html)

[A.36 library(persistency): Provide persistent dynamic predicates](persistency.html)

[A.37 library(pio): Pure I/O](pio.html)

[A.37.1 library(pure_input): Pure Input from files and streams](pio.html#sec:A.37.1)

[A.38 library(portray_text): Portray text](portraytext.html)

[A.39 library(predicate_options): Declare option-processing of predicates](predicate_options.html)

[A.39.1 The strength and weakness of predicate options](predicate_options.html#sec:A.39.1)

[A.39.2 Options as arguments or environment?](predicate_options.html#sec:A.39.2)

[A.39.3 Improving on the current situation](predicate_options.html#sec:A.39.3)

[A.39.3.1 Options as types](predicate_options.html#sec:A.39.3.1)

[A.39.3.2 Reflective access to options](predicate_options.html#sec:A.39.3.2)

[A.40 library(prolog_coverage): Coverage analysis tool](prologcoverage.html)

[A.40.1 Coverage collection and threads](prologcoverage.html#sec:A.40.1)

[A.40.2 Combining coverage data from multiple runs](prologcoverage.html#sec:A.40.2)

[A.40.3 Predicate reference](prologcoverage.html#sec:A.40.3)

[A.41 library(prolog_debug): User level debugging tools](prologdebug.html)

[A.42 library(prolog_jiti): Just In Time Indexing (JITI) utilities](prologjiti.html)

[A.43 library(prolog_trace): Print access to predicates](prologtrace.html)

[A.44 library(prolog_versions): Demand specific (Prolog) versions](prologversions.html)

[A.45 library(prolog_xref): Prolog cross-referencer data collection](prologxref.html)

[A.46 library(quasi_quotations): Define Quasi Quotation syntax](quasiquotations.html)

[A.47 library(random): Random numbers](random.html)

[A.48 library(rbtrees): Red black trees](rbtrees.html)

[A.49 library(readutil): Read utilities](readutil.html)

[A.50 library(record): Access named fields in a term](record.html)

[A.51 library(registry): Manipulating the Windows registry](registry.html)

[A.52 library(rwlocks): Read/write locks](rwlocks.html)

[A.53 library(settings): Setting management](settings.html)

[A.54 library(statistics): Get information about resource usage](statistics.html)

[A.55 library(strings): String utilities](strings.html)

[A.56 library(simplex): Solve linear programming problems](simplex.html)

[A.56.1 Introduction](simplex.html#sec:A.56.1)

[A.56.2 Delayed column generation](simplex.html#sec:A.56.2)

[A.56.3 Solving LPs with special structure](simplex.html#sec:A.56.3)

[A.56.4 Examples](simplex.html#sec:A.56.4)

[A.56.4.1 Example 1](simplex.html#sec:A.56.4.1)

[A.56.4.2 Example 2](simplex.html#sec:A.56.4.2)

[A.56.4.3 Example 3](simplex.html#sec:A.56.4.3)

[A.57 library(solution_sequences): Modify solution sequences](solutionsequences.html)

[A.58 library(tables): XSB interface to tables](tables.html)

[A.59 library(tableutil): Table inspection and statistics utilities](tableutil.html)

[A.60 library(terms): Term manipulation](terms.html)

[A.61 library(thread): High level thread primitives](thread.html)

[A.62 library(thread_pool): Resource bounded thread management](threadpool.html)

[A.63 library(ugraphs): Graph manipulation library](ugraphs.html)

[A.64 library(url): Analysing and constructing URL](url.html)

[A.65 library(varnumbers): Utilities for numbered terms](varnumbers.html)

[A.66 library(yall): Lambda expressions](yall.html)

[B Hackers corner](hack.html)

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

[C Compatibility with other Prolog dialects](dialect.html)

[C.1 Some considerations for writing portable code](portabilitystrategies.html)

[C.2 Notes on specific dialects](dialect-notes.html)

[C.2.1 Notes on specific dialects](dialect-notes.html#sec:C.2.1)

[C.2.1.1 Loading XSB source files](dialect-notes.html#sec:C.2.1.1)

[C.2.2 The XSB import directive](dialect-notes.html#sec:C.2.2)

[D Glossary of Terms](glossary.html)

[E SWI-Prolog License Conditions and Tools](license.html)

[E.1 Contributing to the SWI-Prolog project](contrib.html)

[E.2 Software support to keep track of license conditions](softlicense.html)

[E.3 License conditions inherited from used code](otherlicenses.html)

[E.3.1 Cryptographic routines](otherlicenses.html#sec:E.3.1)

[F Summary](summary.html)

[F.1 Predicates](predsummary.html)

[F.2 Library predicates](library.html)

[F.2.1 library(aggregate)](library.html#sec:F.2.1)

[F.2.2 library(ansi_term)](library.html#sec:F.2.2)

[F.2.3 library(apply)](library.html#sec:F.2.3)

[F.2.4 library(assoc)](library.html#sec:F.2.4)

[F.2.5 library(broadcast)](library.html#sec:F.2.5)

[F.2.6 library(charsio)](library.html#sec:F.2.6)

[F.2.7 library(check)](library.html#sec:F.2.7)

[F.2.8 library(clpb)](library.html#sec:F.2.8)

[F.2.9 library(clpfd)](library.html#sec:F.2.9)

[F.2.10 library(clpqr)](library.html#sec:F.2.10)

[F.2.11 library(csv)](library.html#sec:F.2.11)

[F.2.12 library(dcgbasics)](library.html#sec:F.2.12)

[F.2.13 library(dcghighorder)](library.html#sec:F.2.13)

[F.2.14 library(debug)](library.html#sec:F.2.14)

[F.2.15 library(dicts)](library.html#sec:F.2.15)

[F.2.16 library(error)](library.html#sec:F.2.16)

[F.2.17 library(exceptions)](library.html#sec:F.2.17)

[F.2.18 library(fastrw)](library.html#sec:F.2.18)

[F.2.19 library(explain)](library.html#sec:F.2.19)

[F.2.20 library(help)](library.html#sec:F.2.20)

[F.2.21 library(gensym)](library.html#sec:F.2.21)

[F.2.22 library(heaps)](library.html#sec:F.2.22)

[F.2.23 library(increval)](library.html#sec:F.2.23)

[F.2.24 library(intercept)](library.html#sec:F.2.24)

[F.2.25 library(iostream)](library.html#sec:F.2.25)

[F.2.26 library(listing)](library.html#sec:F.2.26)

[F.2.27 library(lists)](library.html#sec:F.2.27)

[F.2.28 library(macros)](library.html#sec:F.2.28)

[F.2.29 library(main)](library.html#sec:F.2.29)

[F.2.30 library(occurs)](library.html#sec:F.2.30)

[F.2.31 library(option)](library.html#sec:F.2.31)

[F.2.32 library(optparse)](library.html#sec:F.2.32)

[F.2.33 library(ordsets)](library.html#sec:F.2.33)

[F.2.34 library(persistency)](library.html#sec:F.2.34)

[F.2.35 library(portraytext)](library.html#sec:F.2.35)

[F.2.36 library(predicate_options)](library.html#sec:F.2.36)

[F.2.37 library(prologcoverage)](library.html#sec:F.2.37)

[F.2.38 library(prologdebug)](library.html#sec:F.2.38)

[F.2.39 library(prologjiti)](library.html#sec:F.2.39)

[F.2.40 library(prologpack)](library.html#sec:F.2.40)

[F.2.41 library(prologversions)](library.html#sec:F.2.41)

[F.2.42 library(prologtrace)](library.html#sec:F.2.42)

[F.2.43 library(prologxref)](library.html#sec:F.2.43)

[F.2.44 library(pairs)](library.html#sec:F.2.44)

[F.2.45 library(pio)](library.html#sec:F.2.45)

[F.2.45.1 library(pure_input)](library.html#sec:F.2.45.1)

[F.2.46 library(random)](library.html#sec:F.2.46)

[F.2.47 library(rbtrees)](library.html#sec:F.2.47)

[F.2.48 library(readutil)](library.html#sec:F.2.48)

[F.2.49 library(record)](library.html#sec:F.2.49)

[F.2.50 library(registry)](library.html#sec:F.2.50)

[F.2.51 library(rwlocks)](library.html#sec:F.2.51)

[F.2.52 library(settings)](library.html#sec:F.2.52)

[F.2.53 library(simplex)](library.html#sec:F.2.53)

[F.2.54 library(statistics)](library.html#sec:F.2.54)

[F.2.55 library(tableutils)](library.html#sec:F.2.55)

[F.2.56 library(terms)](library.html#sec:F.2.56)

[F.2.57 library(ugraphs)](library.html#sec:F.2.57)

[F.2.58 library(url)](library.html#sec:F.2.58)

[F.2.59 library(writef)](library.html#sec:F.2.59)

[F.2.60 library(www_browser)](library.html#sec:F.2.60)

[F.2.61 library(solution_sequences)](library.html#sec:F.2.61)

[F.2.62 library(thread)](library.html#sec:F.2.62)

[F.2.63 library(thread_pool)](library.html#sec:F.2.63)

[F.2.64 library(varnumbers)](library.html#sec:F.2.64)

[F.2.65 library(yall)](library.html#sec:F.2.65)

[F.3 Arithmetic Functions](funcsummary.html)

[F.4 Operators](opsummary.html)

[G Bibliography](Bibliography.html)
