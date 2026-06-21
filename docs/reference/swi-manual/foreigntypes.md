
## 12.3 Interface Data Types

### 12.3.1 Type `term_t`: a reference to a Prolog term

The principal data type is `term_t`. Type `term_t` is what Quintus calls `QP_term_ref`. This name indicates better what the type represents: it is a *handle* for a term rather than the term itself. Terms can only be represented and manipulated using this type, as this is the only safe way to ensure the Prolog kernel is aware of all terms referenced by foreign code and thus allows the kernel to perform garbage collection and/or stack-shifts while foreign code is active, for example during a callback from C.

A term reference is a C `uintptr_t`, representing the offset of a variable on the Prolog environment stack. A foreign function is passed term references for the predicate arguments, one for each argument. If references for intermediate results are needed, such references may be created using [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()) or [PL_new_term_refs()](foreigntypes.html#PL_new_term_refs()). These references normally live till the foreign function returns control back to Prolog. Their scope can be explicitly limited using [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()) and [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame())/[PL_discard_foreign_frame()](foreigninclude.html#PL_discard_foreign_frame()).

A `term_t` always refers to a valid Prolog term (variable, atom, integer, float or compound term). A term lives either until backtracking takes us back to a point before the term was created, the garbage collector has collected the term, or the term was created after a [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()) and [PL_discard_foreign_frame()](foreigninclude.html#PL_discard_foreign_frame()) has been called.

The foreign interface functions can either *read*, *unify* or *write* to term references. In this document we use the following notation for arguments of type `term_t`:

> |  |  |
> |----|----|
> | `term_t +t` | Accessed in read-mode. The‘+’indicates the argument is‘input’. |
> | `term_t -t` | Accessed in write-mode. |
> | `term_t ?t` | Accessed in unify-mode. |

**WARNING** Term references that are accessed in‘write’(-) mode will refer to an invalid term if the term is allocated on the global stack and backtracking takes us back to a point before the term was written.^(214This could have been avoided by *trailing* term references when data is written to them. This seriously hurts performance in some scenarios though. If this is desired, use [PL_put_variable()](foreigninclude.html#PL_put_variable()) followed by one of the PL_unify\_\*() functions.) Compound terms, dicts, large integers, rational numbers, floats and strings are all allocated on the global stack. Below is a typical scenario where this may happen. The first solution writes a term extracted from the solution into `a`. After the system backtracks due to [PL_next_solution()](foreigninclude.html#PL_next_solution()), `a` becomes a reference to a term that no longer exists.

``` code
term_t a = PL_new_term_ref();
...
query = PL_open_query(...);
while(PL_next_solution(query))
{ PL_get_arg(i, ..., a);
}
PL_close_query(query);
```

There are two solutions to this problem. One is to scope the term reference using [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()) and [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()) and makes sure it goes out of scope before backtracking happens. The other is to clear the term reference using [PL_put_variable()](foreigninclude.html#PL_put_variable()) before backtracking.

Term references are obtained in any of the following ways:

- *Passed as argument*  
  The C functions implementing foreign predicates are passed their arguments as term references. These references may be read or unified. Writing to these variables causes undefined behaviour.
- *Created by [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref())*  
  A term created by [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()) is normally used to build temporary terms or to be written by one of the interface functions. For example, [PL_get_arg()](foreigninclude.html#PL_get_arg()) writes a reference to the term argument in its last argument.
- *Created by [PL_new_term_refs(size_t n)](foreigntypes.html#PL_new_term_refs())*  
  This function returns a set of term references with the same characteristics as [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()). See [PL_open_query()](foreigninclude.html#PL_open_query()).
- *Created by [PL_copy_term_ref(term_t t)](foreigntypes.html#PL_copy_term_ref())*  
  Creates a new term reference to the same term as the argument. The term may be written to. See [figure 6](foreigninclude.html#fig:pl-display).

Term references can safely be copied to other C variables of type `term_t`, but all copies will always refer to the same term.

`term_t` **PL_new_term_ref**()  
Return a fresh reference to a term. The reference is allocated on the *local* stack. Allocating a term reference may trigger a stack-shift on machines that cannot use sparse memory management for allocation of the Prolog stacks. The returned reference describes a variable. Raise a resource exception and returns `(term_t)0` on failure.

`term_t` **PL_new_term_refs**(`size_t n`)  
Return `n` new term references. The first term reference is returned. The others are `t``+1`, `t``+2`, etc. Raise a resource exception and returns `(term_t)0` on failure. There are two reasons for using this function. [PL_open_query()](foreigninclude.html#PL_open_query()) and [PL_cons_functor()](foreigninclude.html#PL_cons_functor()) expect the arguments as a set of consecutive term references, and *very* time-critical code requiring a number of term references can be written as:

``` code
pl_mypredicate(term_t a0, term_t a1)
{ term_t t0 = PL_new_term_refs(2);
  term_t t1 = t0+1;

  ...
}
```

`term_t` **PL_copy_term_ref**(`term_t from`)  
Create a new term reference and make it point initially to the same term as `from`. This function is commonly used to copy a predicate argument to a term reference that may be written. Raise a resource exception and returns `(term_t)0` on failure. An example of its use is given below, in the sample code **`pl_write_atoms()`**.

`void` **PL_free_term_ref**(`term_t t`)  
Release a specific term reference. Normally all term references in a *scope* are discarded together or all term references created after a specific one are reclaimed using [PL_reset_term_refs()](foreigntypes.html#PL_reset_term_refs()). This function shrinks the current foreign frame if `t` is the last one in the frame. Else it marks `t` for reuse by [PL_new_term_ref()](foreigntypes.html#PL_new_term_ref()).

`void` **PL_reset_term_refs**(`term_t after`)  
Destroy all term references that have been created after `after`, including `after` itself. Any reference to the invalidated term references after this call results in undefined behaviour.

Note that returning from the foreign context to Prolog will reclaim all references used in the foreign context. This call is only necessary if references are created inside a loop that never exits back to Prolog. See also [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()), [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()) and [PL_discard_foreign_frame()](foreigninclude.html#PL_discard_foreign_frame()).

#### 12.3.1.1 Interaction with the garbage collector and stack-shifter

Prolog implements two mechanisms for avoiding stack overflow: garbage collection and stack expansion. On machines that allow for it, Prolog will use virtual memory management to detect stack overflow and expand the runtime stacks. On other machines Prolog will reallocate the stacks and update all pointers to them. To do so, Prolog needs to know which data is referenced by C code. As all Prolog data known by C is referenced through term references (`term_t`), Prolog has all the information necessary to perform its memory management without special precautions from the C programmer.

### 12.3.2 Other foreign interface types

**atom_t**  
The type `atom_t` actually represents a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)). Blobs are the super type of Prolog atoms, where atoms are blobs that represent textual content. Textual content is also represented by Prolog string (see [section 5.2](string.html#sec:5.2)), which makes the general notion of *string* in Prolog ambiguous. The core idea behind blobs/atoms is to represent arbitrary content using a *unique* handle, such that comparing the handles is enough to prove equivalence of the contents; i.e., given two different atom handles we know they represent different texts. This uniqueness feature allows the core engine to reason about atom equality and inequality without considering their content. Blobs without the `PL_BLOB_UNIQUE` feature are also tested for uniqueness without considering their content. Each time an atom or a `PL_BLOB_UNIQUE` blob is created, it must be looked up in the atom table; if a blob without `PL_BLOB_UNIQUE` is created, no lookup is done. *Strings* ([section 5.2](string.html#sec:5.2)) and blobs without the `PL_BLOB_UNIQUE` feature do *not* have this uniqueness property - to test for equality, the contents of the strings or blobs must be compared. For both atoms and strings, comparisons for ordering (e.g., used by [sort/2](builtinlist.html#sort/2) or @\</2) must use the contents; in the case of blobs, [compare()](foreigninclude.html#compare()) can be specified in the `PL_blob_t` structure to override the default bitwise comparison.

Because atoms are often used to represent (parts of) arbitrary input, intermediate results, and output of data processed by Prolog, it is necessary that atoms be subject to *garbage collection* (see [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0)). The garbage collection makes atoms ideal handles for arbitrary data structures, which are generalized as *blobs*. Blobs provide *safe* access to many internal Prolog data structures such as streams, clause references, etc.

**functor_t**  
A functor is the internal representation of a name/arity pair. They are used to find the name and arity of a compound term as well as to construct new compound terms. Like atoms they live for the whole Prolog session and are unique.

**predicate_t**  
Handle to a Prolog predicate. Predicate handles live forever (although they can lose their definition).

**qid_t**  
Query identifier. Used by [PL_open_query()](foreigninclude.html#PL_open_query()), [PL_next_solution()](foreigninclude.html#PL_next_solution()), [PL_cut_query()](foreigninclude.html#PL_cut_query()), and [PL_close_query()](foreigninclude.html#PL_close_query()) to handle calling Prolog from C.

**fid_t**  
Frame identifier. Used by [PL_open_foreign_frame()](foreigninclude.html#PL_open_foreign_frame()) and [PL_close_foreign_frame()](foreigninclude.html#PL_close_foreign_frame()).

**module_t**  
A module is a unique handle to a Prolog module. Modules are used only to call predicates in a specific module.

**foreign_t**  
Return type for a C function implementing a Prolog predicate.

**control_t**  
Passed as additional argument to non-deterministic foreign functions. See PL_retry\*() and PL_foreign_context\*().

**install_t**  
Type for the **install()** and **uninstall()** functions of shared or dynamic link libraries. See [section 12.2.3](foreignlink.html#sec:12.2.3).

**int64_t**  
Actually part of the C99 standard rather than Prolog. As of version 5.5.6, Prolog integers are 64-bit on all hardware. The C99 type `int64_t` is defined in the `stdint.h` standard header and provides platform-independent 64-bit integers. Portable code accessing Prolog should use this type to exchange integer values. Please note that [PL_get_long()](foreigninclude.html#PL_get_long()) can return `FALSE` on Prolog integers that cannot be represented as a C long. Robust code should not assume any of the integer fetching functions to succeed, *even* if the Prolog term is known to be an integer.

#### 12.3.2.1 PL_ARITY_AS_SIZE

As of SWI-Prolog 7.3.12, the arity of terms has changed from `int` to `size_t`. To deal with this transition, all affecting functions have two versions, where the old name exchanges the arity as `int` and a new function with name \***\_sz()** exchanges the arity as `size_t`. Up to 8.1.28, the default was to use the old `int` functions. As of 8.1.29/8.2.x, the default is to use `size_t` and the old behaviour can be restored by defining `PL_ARITY_AS_SIZE` to `0` (zero). This makes old code compatible, but the following warning is printed when compiling:

``` code
#warning "Term arity has changed from int to size_t."
#warning "Please update your code or use #define PL_ARITY_AS_SIZE 0."
```

To make the code compile silently again, change the types you use to represent arity from `int` to `size_t`. Please be aware that `size_t` is *unsigned*. At some point representing arity as `int` will be dropped completely.

#### 12.3.2.2 Notes on C API bool return values

Most of the SWI-Prolog C-API consists of C functions that return a Boolean result. Up to version 9.3.10, these functions are defined to return `int`. Later versions define these functions to return the `bool`. This type is provided by the standard header `stdbool.h` and will be supported as a native type starting with the C23 standard, which introduces the keywords `false`, `true` and `bool`. `SWI-Prolog.h` defines the constants `FALSE` and `TRUE`. These constants are consistent with `false`, and `true` and may be used interchangeably. Future versions will deprecate `FALSE` and `TRUE`. As of version 9.3.11 `SWI-Prolog.h` includes `stdbool.h` and thus provides the standard names.

The Boolean result `true` indicates success, while `false` may indicate an error or *logical failure*. Which of the two happened can be examined by calling [PL_exception(0)](foreigninclude.html#PL_exception()), which returns a `term_t` of value 0 if there was a logical failure. Otherwise the returned term reference is a handle to the Prolog exception. Typically there is no need to test whether or not there has been an exception. Instead, the implementation of a foreign predicate can often simply return `false` in case some API returned `false`. Prolog will map this to logical failure or raise the pending exception. The C API defines several groups of `bool` functions that behave consistently. Note that errors which as the Prolog term handle (`term_t`) not being a valid is not reported through the API. If this is detected [PL_api_error()](foreigninclude.html#PL_api_error()) is called, which aborts the process with a diagnostic message. If not detected, such errors lead to *undefined behaviour* (read: arbitrary crashes or wrong behaviour now or later).

**PL_is\_\*()**  
These are *type checking* functions. They have no side effects and no error conditions. Returning `false` implies the argument is not of the tested type.

**PL_get\_\*()**  
This group extracts C value from a Prolog term. If the term is not of the expected type or the C value cannot represent the value the function returns `false`. No exception is raised.

**PL_get\_\***\_ex()****  
This group is similar to PL_get\_\*(), but raises a Prolog exception. The exception is either an `instantiation_error` in case the term is unbound but should not be, a `type_error` in case the term is of the wrong type or a `representation_error` in case the C type cannot represent the Prolog value (e.g., a C `int` while the Prolog integer is out of reach for this type).

**PL_put\_\*()**  
This group converts C data to a Prolog term. Such a function returning `false` always raises a `resource_error`, indicating that Prolog does not have sufficient resources to store the result.

**PL_unify\_\*()**  
This group unifies a Prolog term to a converted C value. Here, the failure can be *logical* if the unification failed because the term was already bound to some other value or the failure may be the result of a resource error as with the PL_put\_\*() group.
