
## 12.2 Linking Foreign Modules

Foreign modules may be linked to Prolog in two ways. Using *static linking*, the extensions, a (short) file defining **main()** which attaches the extension calls to Prolog, and the SWI-Prolog kernel distributed as a C library, are linked together to form a new executable. Using *dynamic linking*, the extensions are linked to a shared library (`.so` file on most Unix systems) or dynamic link library (`.DLL` file on Microsoft platforms) and loaded into the running Prolog process.^(212The system also contains code to load `.o` files directly for some operating systems, notably Unix systems using the BSD `a.out` executable format. As the number of Unix platforms supporting this quickly gets smaller and this interface is difficult to port and slow, it is no longer described in this manual. The best alternative would be to use the dld package on machines that do not have shared libraries.)

### 12.2.1 What linking is provided?

The *static linking* schema can be used on all versions of SWI-Prolog. Whether or not dynamic linking is supported can be deduced from the Prolog flag [open_shared_object](flags.html#flag:open_shared_object) (see [current_prolog_flag/2](flags.html#current_prolog_flag/2)). If this Prolog flag yields `true`, [open_shared_object/2](foreignlink.html#open_shared_object/2) and related predicates are defined. See [section 12.2.3](foreignlink.html#sec:12.2.3) for a suitable high-level interface to these predicates.

### 12.2.2 What kind of loading should I be using?

All described approaches have their advantages and disadvantages. Static linking is portable and allows for debugging on all platforms. It is relatively cumbersome and the libraries you need to pass to the linker may vary from system to system, though the utility program **swipl-ld** described in [section 12.5](plld.html#sec:12.5) often hides these problems from the user.

Loading shared objects (DLL files on Windows) provides sharing and protection and is generally the best choice. If a saved state is created using [qsave_program/\[1,2\]](saved-states.html#qsave_program/1), an [initialization/1](consulting.html#initialization/1) directive may be used to load the appropriate library at startup.

Note that the definition of the foreign predicates is the same, regardless of the linking type used.

### 12.2.3 library(shlib): Utility library for loading foreign objects (DLLs, shared objects)

This section discusses the functionality of the (autoload) `library(shlib)`, providing an interface to manage shared libraries. We describe the procedure for using a foreign resource (DLL in Windows and shared object in Unix) called `mylib`.

First, one must assemble the resource and make it compatible to SWI-Prolog. The details for this vary between platforms. The `swipl-ld(1)` utility can be used to deal with this in a portable manner. The typical commandline is:

``` code
swipl-ld -shared -o mylib file.{c,o,cc,C} ...
```

Make sure that one of the files provides a global function `install_mylib()` that initialises the module using calls to PL_register_foreign(). Below is a simple example file `mylib.c`, which prints a "hello" message. Note that we use SWI-Prolog's Sprintf() rather than C standard `printf()` to print the outout through Prolog's `current_output` stream, making the example work in a windowed environment. The standard C `printf()` works in a console environment, but this bypasses Prolog's output redirection. Also note the use of the standard C `bool` type, which is supported in 9.2.x and more actively promoted in the 9.3.x development series.

``` code
#include <SWI-Prolog.h>
#include <SWI-Stream.h>
#include <stdbool.h>

static foreign_t
pl_say_hello(term_t to)
{ char *s;

  if ( PL_get_chars(to, &s, CVT_ALL|REP_UTF8) )
  { Sprintf("hello %Us", s);

    return true;
  }

  return false;
}

install_t
install_mylib(void)
{ PL_register_foreign("say_hello", 1, pl_say_hello, 0);
}
```

Now write a file `mylib.pl`:

``` code
:- module(mylib, [ say_hello/1 ]).
:- use_foreign_library(foreign(mylib)).
```

The file `mylib.pl` can be loaded as a normal Prolog file and provides the predicate defined in C. The generated `mylib.so` (or `.dll`, etc.) must be placed in a directory searched for using the Prolog search path `foreign` (see [absolute_file_name/3](files.html#absolute_file_name/3)). To load this from the current directory, we can use the `-p alias=dir` option:

``` code
swipl -p foreign=. mylib.pl
?- say_hello(world).
hello world
true.
```

\[det\]**use_foreign_library**(`+FileSpec`)  
\[det\]**use_foreign_library**(`+FileSpec, +Options:list`)  
Load and install a foreign library as [load_foreign_library/1](foreignlink.html#load_foreign_library/1),2 and register the installation using [initialization/2](consulting.html#initialization/2) with the option `now`. This is similar to using:

``` code
:- initialization(load_foreign_library(foreign(mylib))).
```

but using the [initialization/1](consulting.html#initialization/1) wrapper causes the library to be loaded *after* loading of the file in which it appears is completed, while [use_foreign_library/1](foreignlink.html#use_foreign_library/1) loads the library *immediately*. I.e. the difference is only relevant if the remainder of the file uses functionality of the C-library.

As of SWI-Prolog 8.1.22, [use_foreign_library/1](foreignlink.html#use_foreign_library/1),2 is in provided as a built-in predicate that, if necessary, loads `library(shlib)`. This implies that these directives can be used without explicitly loading `library(shlib)` or relying on demand loading.

\[semidet,multifile\]qsave:**compat_arch**(`Arch1, Arch2`)  
User definable hook to establish if `Arch1` is compatible with `Arch2` when running a shared object. It is used in saved states produced by [qsave_program/2](saved-states.html#qsave_program/2) to determine which shared object to load at runtime.

See also  
`foreign` option in [qsave_program/2](saved-states.html#qsave_program/2) for more information.

\[det\]**load_foreign_library**(`:FileSpec`)  
\[det\]**load_foreign_library**(`:FileSpec, +Options:list`)  
Load a *shared object* or *DLL*. After loading the Entry function is called without arguments. The default entry function is composed from =install\_=, followed by the file base-name. E.g., the load-call below calls the function `install_mylib()`. If the platform prefixes extern functions with =\_=, this prefix is added before calling. `Options` provided are below. Other options are passed to [open_shared_object/3](foreignlink.html#open_shared_object/3).

**install**(`+Function`)  
Installation function to use. Default is `default(install)`, which derives the function from `FileSpec`.

``` code
    ...
    load_foreign_library(foreign(mylib)),
    ...
```

|  |  |
|----|----|
| `FileSpec` | is a specification for [absolute_file_name/3](files.html#absolute_file_name/3). If searching the file fails, the plain name is passed to the OS to try the default method of the OS for locating foreign objects. The default definition of [file_search_path/2](consulting.html#file_search_path/2) searches `<`prolog home`>`/lib/`<`arch`>` on Unix and `<`prolog home`>`/bin on Windows. |

See also  
[use_foreign_library/1](foreignlink.html#use_foreign_library/1),2 are intended for use in directives.

\[det\]**unload_foreign_library**(`+FileSpec`)  
\[det\]**unload_foreign_library**(`+FileSpec, +Exit:atom`)  
Unload a *shared object* or *DLL*. After calling the `Exit` function, the shared object is removed from the process. The default exit function is composed from =uninstall\_=, followed by the file base-name.

**current_foreign_library**(`?File, ?Public`)  
Query currently loaded shared libraries.

**reload_foreign_libraries**  
Reload all foreign libraries loaded (after restore of a state created using [qsave_program/2](saved-states.html#qsave_program/2).

### 12.2.4 Low-level operations on shared libraries

The interface defined in this section allows the user to load shared libraries (`.so` files on most Unix systems, `.dll` files on Windows). This interface is portable to Windows as well as to Unix machines providing **dlopen**(2) (Solaris, Linux, FreeBSD, Irix and many more) or **shl_open**(2) (HP/UX). It is advised to use the predicates from [section 12.2.3](foreignlink.html#sec:12.2.3) in your application.

**open_shared_object**(`+File, -Handle`)  
`File` is the name of a *shared object* file (DLL in MS-Windows). This file is attached to the current process, and `Handle` is unified with a handle to the library. Equivalent to `open_shared_object(File, Handle, [])`. See also [open_shared_object/3](foreignlink.html#open_shared_object/3), [load_foreign_library/1](foreignlink.html#load_foreign_library/1) and [use_foreign_library/1](foreignlink.html#use_foreign_library/1).

On errors, an exception `shared_object(Action, Message)` is raised. `Message` is the return value from **dlerror()**.

**open_shared_object**(`+File, -Handle, +Options`)  
As [open_shared_object/2](foreignlink.html#open_shared_object/2), but allows for additional flags to be passed. Defined options are below. These options map to `RTLD_NOW`, `RTLD_LAZY`, `RTLD_GLOBAL`, `RTLD_NODELETE`, `RTLD_NOLOAD` and `RTLD_DEEPBIND` on systems where this predicate is implemented using **dlopen()** and these flags are supported. If the flag is not supported on the target OS, the corresponding option is silently ignored.

**resolve**(`Atom`)  
When to resolve symbols. Values are `lazy` (default) or `now`.

**visibility**(`Atom`)  
Visibility of the new symbols. Values are are `local` (default) or `global`, making the new symbols available to all subsequently loaded shared objects.

**now**(`Bool`)  
`now(true)` is the same as `resolve(now)`. Provided for backward compatibility.

**global**(`Bool`)  
`global(true)` is the same as `visibility(global)`. Provided for backward compatibility.

**delete**(`Bool`)  
If `false`, include `RTLD_NODELETE`.

**load**(`Bool`)  
if `false`, include `RTLD_NOLOAD`. This returns a handle to the object if it is already loaded and `NULL` otherwise. It causes this predicate to *fail silently* if the object is not loaded.

**deepbind**(`Bool`)  
if `true`, include `RTLD_DEEPBIND`.

Note that these flags may not be supported by your operating system. Check the documentation of **dlopen()** or equivalent on your operating system. Unsupported flags are silently ignored.

**close_shared_object**(`+Handle`)  
Detach the shared object identified by `Handle`.

**call_shared_object_function**(`+Handle, +Function`)  
Call the named function in the loaded shared library. The function is called without arguments and the return value is ignored. Normally this function installs foreign language predicates using calls to [PL_register_foreign()](foreigninclude.html#PL_register_foreign()).

### 12.2.5 Static Linking

Older versions of SWI-Prolog were shipped by default with a *static library*. In recent versions we no longer ship a static library because practically every OS properly supports dynamic linking without serious drawbacks and dynamic linking has several advantages. It is on many platforms required to be able to load SWI-Prolog foreign libraries (see [use_foreign_library/1](foreignlink.html#use_foreign_library/1)). Only on ELF based systems such as Linux we can load foreign libraries *if* the main executable is linked to export its global symbols (**gcc** `-rdynamic` option). Another advantage of dynamic libraries is that the user does not have to worry about libraries that this particular build of SWI-Prolog requires such as `libgmp` as well as OS specific libraries.

If one really wants a static library, use the **CMake** flag `-DSWIPL_STATIC_LIB=ON` while configuring a build from source. This causes building and installing `libswipl_static.a`. Note the `_static` postfix to avoid a name conflict on Windows between the *import library* and the *static library*.^(213As is, the Windows build is cross-compiled using MinGW which produces `libswipl_static.a`. This file can, as far as we know, *not* be used by MSVC.).
