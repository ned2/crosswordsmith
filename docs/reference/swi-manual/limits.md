
## 2.19 System limits

### 2.19.1 Limits on memory areas

The SWI-Prolog engine uses three *stacks* the *local stack* (also called *environment stack*) stores the environment frames used to call predicates as well as choice points. The *global stack* (also called *heap*) contains terms, floats, strings and large integers. Finally, the *trail stack* records variable bindings and assignments to support *backtracking*. Except for available memory, there is no *hard limit* for the sizes of the stacks.^(43As of version 9.3.6. Older versions have a hard limit on 32-bit hardware of 128Mb for each stack.)

The combined stack size (per thread) has a *soft limit* implemented by the writeable flag [stack_limit](flags.html#flag:stack_limit) or the command line option **--stack-limit**. Currently the default limit is 1Gb. Considering portability, applications that need to modify the default limits are advised to do so using the Prolog flag [stack_limit](flags.html#flag:stack_limit).

[TABLE]

**Table 2 :** Memory areas

#### 2.19.1.1 The heap

With the heap, we refer to the memory area used by **malloc()** and friends. SWI-Prolog uses the area to store atoms, functors, predicates and their clauses, records and other dynamic data. No limits are imposed on the addresses returned by **malloc()** and friends.

### 2.19.2 Other Limits

**Clauses**  
The only limit on clauses is their arity (the number of arguments to the head), which is limited to 1024. Raising this limit is easy and relatively cheap; removing it is harder.

**Atoms and Strings**  
SWI-Prolog has no limits on the length of atoms or strings. The number of atoms is unlimited. Atoms are subject to garbage collection. See [section 12.4.2.1](foreigninclude.html#sec:12.4.2.1). Both atoms and strings can represent all Unicode code points, including 0 (`\u0000`). Currently, SWI-Prolog uses a separate representation for ISO Latin 1 text (code points `0 ... 255`) and text that includes higher code points. The latter is represented using the C `wchar_t` type. On most systems this implies UCS-4, i.e., 32-bit unsigned integers. On Windows `wchar_t` uses UTF-16, which implies that it cannot represent the code points reserved for *surrogate pairs* as single code points. Future versions may switch to using UTF-8 throughout.

**Nesting of terms**  
Most built-in predicates that process Prolog terms create an explicitly managed stack and perform optimization for processing the last argument of a term. This implies they can process deeply nested terms at constant and low usage of the C stack, and the system raises a resource error if no more stack can be allocated. Currently only [read/1](termrw.html#read/1) and [write/1](termrw.html#write/1) (and all variations thereof) still use the C stack and may cause the system to crash in an uncontrolled way (i.e., not mapped to a Prolog exception that can be caught).

**Integers**  
SWI-Prolog has two integer representations. *Tagged integers* are currently limited to 57 bits.^(44Before version 9.3.6, tagged integers on 32-bit systems had 25 bits and there was a third representation for 64 bit integers.) Unbounded integers are by default provided by the GNU GMP library. Alternatively, they may be provided by the bundled LibBf library. The system can be built without support for unbounded integers.

**Floating point numbers**  
Floating point numbers are represented as C-native double precision floats, 64-bit IEEE on most machines.

### 2.19.3 Reserved Names

The boot compiler (see **-b** option) does not support the module system. As large parts of the system are written in Prolog itself we need some way to avoid name clashes with the user's predicates, database keys, etc. Like Edinburgh C-Prolog [Pereira, 1986](Bibliography.html#CPROLOG:manual) all predicates, database keys, etc., that should be hidden from the user start with a dollar (`$`) sign.
