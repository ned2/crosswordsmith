
## 12.8 Notes on Using Foreign Code

### 12.8.1 Foreign debugging functions

The functions in this section are primarily intended for debugging foreign extensions or embedded Prolog. Violating the constraints of the foreign interface often leads to crashes in a subsequent garbage collection. If this happens, the system needs to be compiled for debugging using `cmake -DCMAKE_BUILD_TYPE=Debug`, after which all functions and predicates listed below are available to use from the debugger (e.g. **gdb**) or can be placed at critical location in your code or the system code.

`void` **PL_backtrace**(`int depth, int flags`)  
Dump a Prolog backtrace to the `user_error` stream. `Depth` is the number of frames to dump. `Flags` is a bitwise or of the following constants:

**PL_BT_SAFE**  
(0x1) Do not try to print *goals*. Instead, just print the predicate name and arity. This reduces the likelihood to crash if [PL_backtrace()](foreignnotes.html#PL_backtrace()) is called in a damaged environment.

**PL_BT_USER**  
(0x2) Only show‘user’frames. Default is to also show frames of hidden built-in predicates.

`char *` **PL_backtrace_string**(`int depth, int flags`)  
As [PL_backtrace()](foreignnotes.html#PL_backtrace()), but returns the stack as a string. The string uses UTF-8 encoding. The returned string must be freed using [PL_free()](foreignnotes.html#PL_free()). This function is was added to get stack traces from running servers where I/O is redirected or discarded. For example, using **gdb**, a stack trace is printed in the gdb console regardless of Prolog I/O redirection using the following command:

``` code
(gdb) printf "%s", PL_backtrace_string(25,0)
```

The source distribution provides the script `scripts/swipl-bt` that exploits **gdb** and [PL_backtrace_string()](foreignnotes.html#PL_backtrace_string()) to print stack traces in various formats for a SWI-Prolog process, given its process id.

`bool` **PL_check_data**(`term_t data`)  
Check the consistency of the term `data`. Returns `TRUE` this is actually implemented in the current version and `FALSE` otherwise. The actual implementation only exists if the system is compiled with the cflag `-DO_DEBUG` or `-DO_MAINTENANCE`. This is *not* the default.

`bool` **PL_check_stacks**()  
Check the consistency of the runtime stacks of the calling thread. Returns `TRUE` this is actually implemented in the current version and `FALSE` otherwise. The actual implementation only exists if the system is compiled with the cflag `-DO_DEBUG` or `-DO_MAINTENANCE`. This is *not* the default.

The Prolog kernel sources use the macro **DEBUG(Topic, Code)**. These macros are disabled in the production version and must be enabled by recompiling the system as described above. Specific topics can be enabled and disabled using the predicates [prolog_debug/1](foreignnotes.html#prolog_debug/1) and [prolog_nodebug/1](foreignnotes.html#prolog_nodebug/1). In addition, they can be activated from the commandline using commandline option `-d topics`, where `topics` is a comma-separated list of debug topics to enable. For example, the code below adds many consistency checks and prints messages if the Prolog signal handler dispatches signals.

``` code
$ swipl -d chk_secure,msg_signal
```

**prolog_debug**(`+Topic`)  
**prolog_nodebug**(`+Topic`)  
Enable/disable a debug topic. `Topic` is an atom that identifies the desired topic. The available topics are defined in `src/pl-debug.h`. Please search the sources to find out what is actually printed and when. We highlight one topic here:

**chk_secure**(`A`)  
dd many expensive consistency checks to the system. This should typically be used when the system crashes, notably in the garbage collector. Garbage collection crashes are in most cases caused by invalid data on the Prolog stacks. This debug topic may help locating how the invalid data was created.

These predicates require the system to be compiled for debugging using `cmake -DCMAKE_BUILD_TYPE=Debug`.

`int` **PL_prolog_debug**(`const char *topic`)  
`int` **PL_prolog_nodebug**(`const char *topic`)  
(De)activate debug topics. The `topics` argument is a comma-separated string of topics to enable or disable. Matching is case-insensitive. See also [prolog_debug/1](foreignnotes.html#prolog_debug/1) and [prolog_nodebug/1](foreignnotes.html#prolog_nodebug/1).

These functions require the system to be compiled for debugging using `cmake -DCMAKE_BUILD_TYPE=Debug`.

### 12.8.2 Memory Allocation

SWI-Prolog's heap memory allocation is based on the **malloc**(3) library routines. SWI-Prolog provides the functions below as a wrapper around **malloc()**. Allocation errors in these functions trap SWI-Prolog's fatal-error handler, in which case [PL_malloc()](foreignnotes.html#PL_malloc()) or [PL_realloc()](foreignnotes.html#PL_realloc()) do not return.

Portable applications must use [PL_free()](foreignnotes.html#PL_free()) to release strings returned by [PL_get_chars()](foreigninclude.html#PL_get_chars()) using the `BUF_MALLOC` argument. Portable applications may use both [PL_malloc()](foreignnotes.html#PL_malloc()) and friends or **malloc()** and friends but should not mix these two sets of functions on the same memory.

`void *` **PL_malloc**(`size_t bytes`)  
Allocate `bytes` of memory. On failure SWI-Prolog's fatal-error handler is called and [PL_malloc()](foreignnotes.html#PL_malloc()) does not return. Memory allocated using these functions must use [PL_realloc()](foreignnotes.html#PL_realloc()) and [PL_free()](foreignnotes.html#PL_free()) rather than **realloc()** and **free()**.

`void *` **PL_realloc**(`void *mem, size_t size`)  
Change the size of the allocated chunk, possibly moving it. The `mem` argument must be obtained from a previous [PL_malloc()](foreignnotes.html#PL_malloc()) or [PL_realloc()](foreignnotes.html#PL_realloc()) call.

`void` **PL_free**(`void *mem`)  
Release memory. The `mem` argument must be obtained from a previous [PL_malloc()](foreignnotes.html#PL_malloc()) or [PL_realloc()](foreignnotes.html#PL_realloc()) call.

### 12.8.3 Compatibility between Prolog versions

Great care is taken to ensure binary compatibility of foreign extensions between different Prolog versions. Only the much less frequently used stream interface has been responsible for binary incompatibilities.

Source code that relies on new features of the foreign interface can use the macro `PLVERSION` to find the version of `SWI-Prolog.h` and [PL_query()](foreigninclude.html#PL_query()) using the option `PL_QUERY_VERSION` to find the version of the attached Prolog system. Both follow the same numbering schema explained with [PL_query()](foreigninclude.html#PL_query()).

### 12.8.4 Foreign hash tables

As of SWI-Prolog 8.3.2 the foreign API provides access to the internal thread-safe and lock-free hash tables that associate pointers or objects that fit in a pointer such as atoms (`atom_t`). An argument against providing these functions is that they have little to do with Prolog. The argument in favor is that it is hard to implement efficient lock-free tables without low-level access to the underlying Prolog threads and exporting this interface has a low cost.

The functions below **can only be called if the calling thread is associated with a Prolog thread**. Failure to do so causes the call to be ignored or return the failure code where applicable.

`hash_table_t` **PL_new_hash_table**(`int size, void (*free_symbol)(void *n, void *v)`)  
Create a new table for `size` key-value pairs. The table is resized when needed. If you know the table will hold 10,000 key-value pairs, providing a suitable initial size avoids resizing. The `free_symbol` function is called whenever a key-value pair is removed from the table. This can be `NULL`.

`int` **PL_free_hash_table**(`hash_table_t table`)  
Destroy the hash table. First calls [PL_clear_hash_table()](foreignnotes.html#PL_clear_hash_table()).

`void*` **PL_lookup_hash_table**(`hash_table_t table, void *key`)  
Return the value matching `key` or `NULL` if `key` does not appear in the table.

`void*` **PL_add_hash_table**(`hash_table_t table, void *key, void *value, int flags`)  
Add `key`-`value` to the table. The behaviour if `key` is already in the table depends on `flags`. If `0`, this function returns the existing value without updating the table. If `PL_HT_UPDATE` the old `value` is *replaced* and the function returns the old value. If `PL_HT_NEW`, a message and backtrace are printed and the function returns `NULL` if `key` is already in the table.

`void*` **PL_del_hash_table**(`hash_table_t table, void *key`)  
Delete `key` from the table, returning the old associated value or `NULL`

`int` **PL_clear_hash_table**(`hash_table_t table`)  
Delete all key-value pairs from the table. Call `free_symbol` for each deleted pair.

`hash_table_enum_t` **PL_new_hash_table_enum**(`hash_table_t table`)  
Return a table *enumerator* (cursor) that can be used to enumerate all key-value pairs using [PL_advance_hash_table_enum()](foreignnotes.html#PL_advance_hash_table_enum()). The enumerator must be discarded using [PL_free_hash_table_enum()](foreignnotes.html#PL_free_hash_table_enum()). It is safe for another thread to add symbols while the table is being enumerated, but undefined whether or not these new symbols are visible. If another thread deletes a key that is not yet enumerated it will not be enumerated.

`void` **PL_free_hash_table_enum**(`hash_table_enum_t e`)  
Discard an enumerator object created using [PL_new_hash_table_enum()](foreignnotes.html#PL_new_hash_table_enum()). Failure to do so causes the table to use more and more memory on subsequent modifications.

`int` **PL_advance_hash_table_enum**(`hash_table_enum_t e, void **key, void **value`)  
Get the next key-value pair from a cursor.

### 12.8.5 Debugging and profiling foreign code (valgrind, asan)

This section is only relevant for Unix users on platforms supported by [valgrind](http://valgrind.org/). Valgrind is an excellent binary instrumentation platform. Unlike many other instrumentation platforms, valgrind can deal with code loaded through **dlopen()**.

The **callgrind** tool can be used to profile foreign code loaded under SWI-Prolog. Compile the foreign library adding **-g** option to **gcc** or **swipl-ld**. By setting the environment variable `VALGRIND` to `yes`, SWI-Prolog will *not* release loaded shared objects using **dlclose()**. This trick is required to get source information on the loaded library. Without, valgrind claims that the shared object has no debugging information.^(242Tested using valgrind version 3.2.3 on x64.) Here is the complete sequence using **bash** as login shell:

``` code
% VALGRIND=yes valgrind --tool=callgrind pl <args>
<prolog interaction>
% kcachegrind callgrind.out.<pid>
```

Instead of **valgrind**, you can use [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer). Here is a short example for building with *asan* enabled and then running the resulting binary. When you exit **swipl**, a message is printed and any leaks are printed. During execution, other messages may be printed out, such as freeing an address twice or using freed or unallocated memory.

``` code
% cd build.sanitize
% cmake -G Ninja -DCMAKE_BUILD_TYPE=Sanitize ..
% ninja
% ASAN_OPTIONS=detect_leaks=1 build.sanitize/src/swipl
<prolog interaction>
% halt
Running LSAN memory leak check (reclaim_memory=1)
No leaks detected
```

### 12.8.6 Name Conflicts in C modules

In the current version of the system all public C functions of SWI-Prolog are in the symbol table. This can lead to name clashes with foreign code. Someday I should write a program to strip all these symbols from the symbol table (why does Unix not have that?). For now I can only suggest you give your function another name. You can do this using the C preprocessor. If---for example---your foreign package uses a function **warning()**, which happens to exist in SWI-Prolog as well, the following macro should fix the problem:

``` code
#define warning warning_
```

Note that shared libraries do not have this problem as the shared library loader will only look for symbols in the main executable for symbols that are not defined in the library itself.

### 12.8.7 Compatibility of the Foreign Interface

The term reference mechanism was first used by Quintus Prolog version 3. SICStus Prolog version 3 is strongly based on the Quintus interface. The described SWI-Prolog interface is similar to using the Quintus or SICStus interfaces, defining all foreign-predicate arguments of type `+term`. SWI-Prolog explicitly uses type `functor_t`, while Quintus and SICStus use \<`name`\> and \<`arity`\>. As the names of the functions differ from Prolog to Prolog, a simple macro layer dealing with the names can also deal with this detail. For example:

``` code
#define QP_put_functor(t, n, a) \
        PL_put_functor(t, PL_new_functor(n, a))
```

The `PL_unify_*()` functions are lacking from the Quintus and SICStus interface. They can easily be emulated, or the put/unify approach should be used to write compatible code.

The [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame())/[PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()) combination is lacking from both other Prologs. SICStus has [PL_new_term_refs(0)](foreigntypes.html#PL_new_term_refs()), followed by [PL_reset_term_refs()](foreigntypes.html#PL_reset_term_refs()), that allows for discarding term references.

The Prolog interface for the graphical user interface package XPCE shares about 90% of the code using a simple macro layer to deal with different naming and calling conventions of the interfaces.
