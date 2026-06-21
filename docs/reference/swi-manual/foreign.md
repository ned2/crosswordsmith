
# 12 Foreign Language Interface

SWI-Prolog offers a powerful interface to C [Kernighan & Ritchie, 1978](Bibliography.html#Kernighan:78). The main design objectives of the foreign language interface are flexibility and performance. A foreign predicate is a C function that has the same number of arguments as the predicate represented. C functions are provided to analyse the passed terms, convert them to basic C types as well as to instantiate arguments using unification. Non-deterministic foreign predicates are supported, providing the foreign function with a handle to control backtracking.

C can call Prolog predicates, providing both a query interface and an interface to extract multiple solutions from a non-deterministic Prolog predicate. There is no limit to the nesting of Prolog calling C, calling Prolog, etc. It is also possible to write the‘main’in C and use Prolog as an embedded logical engine.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

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
