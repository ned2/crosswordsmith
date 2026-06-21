
## 2.12 Environment Control (Prolog flags)

The predicates [current_prolog_flag/2](flags.html#current_prolog_flag/2) and [set_prolog_flag/2](flags.html#set_prolog_flag/2) allow the user to examine and modify the execution environment. It provides access to whether optional features are available on this version, operating system, foreign code environment, command line arguments, version, as well as runtime flags to control the runtime behaviour of certain predicates to achieve compatibility with other Prolog environments.

\[ISO\]**current_prolog_flag**(`?Key, -Value`)  
The predicate [current_prolog_flag/2](flags.html#current_prolog_flag/2) defines an interface to installation features: options compiled in, version, home, etc. With both arguments unbound, it will generate all defined Prolog flags. With `Key` instantiated, it unifies `Value` with the value of the Prolog flag or fails if the `Key` is not a Prolog flag.

Flags marked *changeable* can be modified by the user using [set_prolog_flag/2](flags.html#set_prolog_flag/2). Flag values are typed. Flags marked as `bool` can have the values `true` or `false`. The predicate [create_prolog_flag/3](flags.html#create_prolog_flag/3) may be used to create flags that describe or control behaviour of libraries and applications. The library `library(settings)` provides an alternative interface for managing notably application parameters.

Some Prolog flags are not defined in all versions, which is normally indicated in the documentation below as *“if present and true’*. A boolean Prolog flag is true iff the Prolog flag is present **and** the `Value` is the atom `true`. Tests for such flags should be written as below:

``` code
        (   current_prolog_flag(windows, true)
        ->  <Do MS-Windows things>
        ;   <Do normal things>
        )
```

Some Prolog flags are scoped to a source file. This implies that if they are set using a directive inside a file, the flag value encountered when loading of the file started is restored when loading of the file is completed. Currently, the following flags are scoped to the source file: [generate_debug_info](flags.html#flag:generate_debug_info) and [optimise](flags.html#flag:optimise).

A new thread (see [section 10](threads.html#sec:10)) *copies* all flags from the thread that created the new thread (its *parent*).^(15This is implemented using the copy-on-write technique.) As a consequence, modifying a flag inside a thread does not affect other threads.

**abi_version**(`dict`)  
The flag value is a dict with keys that describe the version of the various Application Binary Interface (ABI) components. See [section 2.21](abi-versions.html#sec:2.21) for details.

**access_level**(`atom`, changeable)  
This flag defines a normal‘user’view (`user`, default) or a‘system’view. In system view all system code is fully accessible as if it was normal user code. In user view, certain operations are not permitted and some details are kept invisible. We leave the exact consequences undefined, but, for example, system code can be traced using system access and system predicates can be redefined.

**address_bits**(`integer`)  
Address size of the hosting machine. Typically 32 or 64. Except for the maximum stack limit, this has few implications to the user. See also the Prolog flag [arch](flags.html#flag:arch).

**agc_close_streams**(`boolean`, changeable)  
When `true` (default `false`^(16Future versions are likely to change the default to `true`.)), that atom garbage collector streams that are garbage collected while being open. In addition, a warning is printed. Below is an example of such a warning.

``` code
WARNING: AGC: closed <stream>(0x560e29014400)
```

Note that closing I/O streams should not be left to the (atom) garbage collector because it may take long before the atom garbage collector runs and because that atom garbage collector is *conservative*, which implies that it is not guaranteed that all garbage atoms are reclaimed. Code that uses I/O streams should use [setup_call_cleanup/3](metacall.html#setup_call_cleanup/3) using the skeleton below, where process/1 is a predicate that reads from or writes to `Stream`.

``` code
    setup_call_cleanup(
        open(..., Stream),
        process(Stream),
        close(Stream)),
    ...
```

Note that the setting for this flag in the `main` thread applies.

**agc_margin**(`integer`, changeable)  
If this amount of atoms possible garbage atoms exist perform atom garbage collection at the first opportunity. Initial value is 10,000. May be changed. A value of 0 (zero) disables atom garbage collection. See also [PL_register_atom()](foreigninclude.html#PL_register_atom()).^(17Given that SWI-Prolog has no limit on the length of atoms, 10,000 atoms may still occupy a lot of memory. Applications using extremely large atoms may wish to call [garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0) explicitly or lower the margin.)

**allow_dot_in_atom**(`bool`, changeable)  
If `true` (default `false`), dots may be embedded into atoms that are not quoted and start with a letter. The embedded dot *must* be followed by an identifier continuation character (i.e., letter, digit or underscore). The dot is allowed in identifiers in many languages, which can make this a useful flag for defining DSLs. Note that this conflicts with cascading functional notation. For example, `Post.meta.author` is read as `.(Post,’meta.author’` if this flag is set to `true`.

**allow_variable_name_as_functor**(`bool`, changeable)  
If true (default is false), `Functor(arg)` is read as if it were written `’Functor’(arg)`. Some applications use the Prolog [read/1](termrw.html#read/1) predicate for reading an application-defined script language. In these cases, it is often difficult to explain to non-Prolog users of the application that constants and functions can only start with a lowercase letter. Variables can be turned into atoms starting with an uppercase atom by calling [read_term/2](termrw.html#read_term/2) using the option `variable_names` and binding the variables to their name. Using this feature, F(x) can be turned into valid syntax for such script languages. Suggested by Robert van Engelen. SWI-Prolog specific.

**android**(`bool`)  
If present and true, it indicates we are running on the Android OS. The flag is not present in other operating systems.

**android_api**(`integer`)  
If running on Android, it indicates the compile-time API Level defined by the C macro `__ANDROID_API__`. It is not defined if running on other operating systems. The API level may or may not match the API level of the running device, since it is the API level at compile time.

**answer_write_options**(`term`, changeable)  
This flag is used by the interactive toplevel to print the value if *bindings* (answers). The flag value is passed to [write_term/2](termrw.html#write_term/2) when printing an answer queries. Default is `[quoted(true), portray(true), max_depth(10), attributes(portray)]`.

**apple**(`bool`)  
If present and `true`, the operating system is MacOSX. Defined if the C compiler used to compile this version of SWI-Prolog defines `__APPLE__`. Note that the [unix](flags.html#flag:unix) is also defined for MacOSX.

**apple_universal_binary**(`bool`)  
If present and `true`, SWI-Prolog has been build as a *universal binary*. Universal binaries contain native executable code for multiple architectures. Currently the supported architectures are `x86_64` and `arm64`. The archirecture prefix for components is `fat-darwin` while the [arch](flags.html#flag:arch) depends on the actual CPU type.

**arch**(`atom`)  
Identifier for the hardware and operating system SWI-Prolog is running on. Used to select foreign files for the right architecture. See also [section 12.2.3](foreignlink.html#sec:12.2.3) and [file_search_path/2](consulting.html#file_search_path/2). For Apple, see also [apple_universal_binary](flags.html#flag:apple_universal_binary).

**argv**(`list`, changeable)  
List is a list of atoms representing the application command line arguments. Application command line arguments are those that have *not* been processed by Prolog during its initialization. Note that Prolog's argument processing stops at `--` or the first non-option argument. See also [os_argv](flags.html#flag:os_argv).^(18Prior to version 6.5.2, [argv](flags.html#flag:argv) was defined as [os_argv](flags.html#flag:os_argv) is now. The change was made for compatibility reasons and because the current definition is more practical.)

**associated_file**(`atom`)  
Set if Prolog was started with a prolog file as argument. Used by e.g., [edit/0](edit.html#edit/0) to edit the initial file.

**autoload**(`atom`, changeable)  
This flag controls autoloading predicates based on [autoload/1](module-autoload.html#autoload/1) and [autoload/2](module-autoload.html#autoload/2) as well as predicates from *autoload libraries*. It has the following values:

**false**  
Predicates are never auto-loaded. If predicates have been imported before using [autoload/\[1,2\]](module-autoload.html#autoload/1), load the referenced files immediately using [use_module/\[1,2\]](import.html#use_module/1). Note that most of the development utilities such as [listing/1](listing.html#listing/1) have to be explicitly imported before they can be used at the toplevel.

**explicit**  
Do not autoload from *autoload libraries*, but do use lazy loading for predicates imported using [autoload/\[1,2\]](module-autoload.html#autoload/1).

**user**  
As `false`, but to autoload library predicates into the global `user` module. This makes the development tools and library implicitly available to the toplevel, but not to modules.

**user_or_explicit**  
Combines `explicit` with `user`, providing lazy loading of predicates imported using [autoload/\[1,2\]](module-autoload.html#autoload/1) and implicit access to the whole library for the toplevel.

**true**  
Provide full autoloading everywhere. This is the default.

**back_quotes**(`codes,chars,string,symbol_char`, changeable)  
Defines the term-representation for back-quoted material. The default is `codes`. If **--traditional** is given, the default is `symbol_char`, which allows using `` ` `` in operators composed of symbols.^(19Older versions had a boolean flag `backquoted_strings`, which toggled between `string` and `symbol_char`) See also [section 5.2](string.html#sec:5.2).

**backtrace**(`bool`, changeable)  
If `true` (default), print a backtrace on an uncaught exception.

**backtrace_depth**(`integer`, changeable)  
If backtraces on errors are enabled, this flag defines the maximum number of frames that is printed (default 20).

**backtrace_goal_depth**(`integer`, changeable)  
The frame of a backtrace is printed after making a shallow copy of the goal. This flag determines the depth to which the goal term is copied. Default is‘3’.

**backtrace_show_lines**(`bool`, changeable)  
If `true` (default), try to reconstruct the line number at which the exception happened.

**bounded**(`bool`)  
ISO Prolog flag. If `true`, integer representation is bound by [min_integer](flags.html#flag:min_integer) and [max_integer](flags.html#flag:max_integer). If `false` integers can be arbitrarily large and the [min_integer](flags.html#flag:min_integer) and [max_integer](flags.html#flag:max_integer) are not present. The flag [max_integer_size](flags.html#flag:max_integer_size) may be used to enforce an arbitrary limit rather than exhausting memory. See [section 4.27.2.1](arith.html#sec:4.27.2.1).

**break_level**(`integer`)  
Current break-level. The initial top level (started with **-t**) has value 0. See [break/0](toplevel.html#break/0). This flag is absent from threads that are not running a top-level loop.

**build_type**(`atom`)  
This flag represents the **cmake** `CMAKE_BUILD_TYPE` to build this instance of SWI-Prolog. Possible values depend on the platform. Some common values are `Debug`, `Release` `MinSizeRel`, `RelWithDebInfo`, `Sanitize`, `DEB` or `PGO`.

**bundle**(`bool`)  
True when SWI-Prolog is installed as a stand-alone bundle. This is set for both the Windows and MacOS binary packages as distributed from the SWI-Prolog download page. This is used to adjust the file search configuration.

**c_cc**(`atom`, changeable)  
Name of the C compiler used to compile SWI-Prolog. Normally one of `gcc`, `clang` or `cc`. See [section 12.5](plld.html#sec:12.5).

**c_cflags**(`atom`, changeable)  
CFLAGS used to compile SWI-Prolog. See [section 12.5](plld.html#sec:12.5).

**c_cxx**(`atom`, changeable)  
Name of the C++ compiler used to test the SWI-Prolog C++ binding. This is the default C++ compiler used by **swipl-ld** (see [section 12.5](plld.html#sec:12.5)) as well as compiling packs using the default setup. Note that SWI-Prolog itself does not contain C++ code and the C++ binding is *header only*. This implies that C++ ABI compatibility issues can not occur.

**c_ldflags**(`atom`, changeable)  
LDFLAGS used to link SWI-Prolog. See [section 12.5](plld.html#sec:12.5).

**c_libplso**(`atom`, changeable)  
Libraries needed to link extensions (shared object, DLL) to SWI-Prolog. Typically empty on ELF systems and `-lswipl` on COFF-based systems. See [section 12.5](plld.html#sec:12.5).

**c_libs**(`atom`, changeable)  
Libraries needed to link executables that embed SWI-Prolog. Typically `-lswipl` if the SWI-Prolog kernel is a shared (DLL). If the SWI-Prolog kernel is in a static library, this flag also contains the dependencies.

**char_conversion**(`bool`, changeable)  
Determines whether character conversion takes place while reading terms. See also [char_conversion/2](charconv.html#char_conversion/2).

**character_escapes**(`bool`, changeable)  
If `true` (default), [read/1](termrw.html#read/1) interprets `\` escape sequences in quoted atoms and strings. May be changed. This flag is local to the module in which it is changed. See [section 2.15.1.3](syntax.html#sec:2.15.1.3).

**character_escapes_unicode**(`bool`, changeable)  
If `true` (default), [write/1](termrw.html#write/1) and friends write escaped characters using the `\uXXXX` or `\UXXXXXXXX` syntax rather than the ISO Prolog `\x<hex>\` syntax. SWI-Prolog reads both.

**ci_speedup**(`float`, changeable)  
Consider generating a hash for clause indexing if the hash has a speedup of at least this flag. Default is 1.5.

**ci_max_var_fraction**(`float`, changeable)  
Do not create a clause indexing hash table of the argument is unbound in the clause for more than this fraction of clauses. Default is 0.1.

**ci_min_speedup_ratio**(`float`, changeable)  
Consider adding a multi-argument hash if it is at lease this values as efficient. Default is 3.0.

**ci_max_lookahead**(`integer`, changeable)  
If a clause is found, scan the clause list of a possible alternative match for at max this number of clauses. Default is 100.

**ci_min_clauses**(`integer`, changeable)  
If the primary index argument (first) is instantiated, still consider a hash of the predicate has more than this number of clauses. Default is 10.

**cmake_build_type**(`atom`, changeable)  
Provides the [cmake](https://cmake.org/) *build type* used to build this version of SWI-Prolog.

**colon_sets_calling_context**(`bool`, changeable)  
Using the construct \<`module`\>:\<`goal`\> sets the *calling context* for executing \<`goal`\>. This flag is defined by ISO/IEC 13211-2 (Prolog modules standard). See [section 6](modules.html#sec:6).

**color_term**(`bool`, changeable)  
This flag is managed by library `library(ansi_term)`, which is loaded at startup if the two conditions below are both true. Note that this implies that setting this flag to `false` from the system or personal initialization file (see [section 2.2](initfile.html#sec:2.2) disables colored output. The predicate [message_property/2](printmsg.html#message_property/2) can be used to control the actual color scheme depending in the message type passed to [print_message/2](printmsg.html#print_message/2).

- `stream_property(current_output, tty(true))`
- `\+ current_prolog_flag(color_term, false)`

**compile_meta_arguments**(`atom`, changeable)  
This flag controls compilation of arguments passed to meta-calls marked‘0’or‘`^`’(see [meta_predicate/1](metapred.html#meta_predicate/1)). Supported values are:

**false**  
(default). Meta-arguments are passed verbatim. If the argument is a control structure ((A,B), (A;B), (A-\>B;C), etc.) it is compile to an temporary clause allocated on the environment stack when the meta-predicate is called.

**control**  
Compile meta-arguments that contain control structures to an auxiliary predicate. This generally improves performance as well as the debugging experience.

**always**  
Always create an intermediate clause, even for system predicates.^(20This may be used in the future for replacing the normal head of the generated predicate with a special reference (similar to database references as used by, e.g., [assert/2](db.html#assert/2)) that provides direct access to the executable code, thus avoiding runtime lookup of predicates for meta-calling.)

**compiled_at**(`atom`)  
Describes when the system has been compiled. Only available if the C compiler used to compile SWI-Prolog provides the \_\_DATE\_\_ and \_\_TIME\_\_ macros.

**conda**(`bool`)  
Set to `true` when built in a [Conda](https://docs.conda.io/) environment.

**console_menu**(`bool`)  
Set to `true` when the I/O is bound to an Epilog (**swipl-win**) Prolog console to indicate that the console supports menus. See also [section 4.35.4](system.html#sec:4.35.4).

**cpu_count**(`integer`, changeable)  
Number of physical CPUs or cores in the system. The flag is marked read-write both to allow pretending the system has more or less processors. See also [thread_setconcurrency/2](threadcreate.html#thread_setconcurrency/2) and the library `library(thread)`. This flag is not available on systems where we do not know how to get the number of CPUs. This flag is not included in a saved state (see [qsave_program/1](saved-states.html#qsave_program/1)).

**dde**(`bool`)  
Set to `true` if this instance of Prolog supports DDE as described in [section 4.44](DDE.html#sec:4.44).

**debug**(`bool`, changeable)  
Switch debugging mode on/off. If debug mode is activated the system traps encountered spy points (see [spy/1](debugger.html#spy/1)) and break points. In addition, last-call optimisation is disabled and the system is more conservative in destroying choice points to simplify debugging.

Disabling these optimisations can cause the system to run out of memory on programs that behave correctly if debug mode is off.

**debug_on_error**(`bool`, changeable)  
If `true`, start the tracer after an error is detected. Otherwise just continue execution. The goal that raised the error will normally fail. See also the Prolog flag [report_error](flags.html#flag:report_error). Default is `true`.

**debug_on_interrupt**(`bool`, changeable)  
If `true`, start the debugger on Control-C.^(21More precisely when receiving `SIGINT`). The initial value is `false` and the value is set to `true` when entering the interactive top level. See **--debug-on-interrupt** to start handling interrupts immediately.

**debugger_show_context**(`bool`, changeable)  
If `true`, show the context module while printing a stack-frame in the tracer. Normally controlled using the‘C’option of the tracer.

**debugger_write_options**(`term`, changeable)  
This argument is given as option-list to [write_term/2](termrw.html#write_term/2) for printing goals by the debugger. Modified by the‘w’,‘p’and‘\<`N`\> d’commands of the debugger. Default is `[quoted(true), portray(true), max_depth(10), attributes(portray)]`.

**determinism_error**(`atom`, changeable)  
This flag defines the behaviour when the predicate determinism is not according to its declaration. See [det/1](debug-determinism.html#det/1). Possible values are `error` (default), `warning` and `silent`.

**dialect**(`atom`)  
Fixed to `swi`. The code below is a reliable and portable way to detect SWI-Prolog.

``` code
is_dialect(swi) :-
        catch(current_prolog_flag(dialect, swi), _, fail).
```

**dir_sep**(`atom`)  
Separator for directories in a file name the OS. Normally `/`, but `\` on Windows.

**double_quotes**(`codes,chars,atom,string`, changeable)  
This flag determines how double quoted strings are read by Prolog and is ---like [character_escapes](flags.html#flag:character_escapes) and [back_quotes](flags.html#flag:back_quotes)--- maintained for each module. The default is `string`, which produces a string as described in [section 5.2](string.html#sec:5.2). If **--traditional** is given, the default is `codes`, which produces a list of character codes, integers that represent a Unicode code-point. The value `chars` produces a list of one-character atoms and the value `atom` makes double quotes the same as single quotes, creating a atom. See also [section 5](extensions.html#sec:5).

**editor**(`atom`, changeable)  
Determines the editor used by [edit/1](edit.html#edit/1). See [section 4.4.1](edit.html#sec:4.4.1) for details on selecting the editor used.

**emacs_inferior_process**(`bool`)  
If true, SWI-Prolog is running as an *inferior process* of (GNU/X-)Emacs. SWI-Prolog assumes this is the case if the environment variable `EMACS` is `t` and `INFERIOR` is `yes`.

**encoding**(`atom`, changeable)  
Default encoding used for opening files in `text` mode. The initial value is deduced from the environment. See [section 2.18.1](widechars.html#sec:2.18.1) for details.

**executable**(`atom`)  
Pathname of the running executable. Used by [qsave_program/2](saved-states.html#qsave_program/2) as default emulator.

**executable_format**(`atom`)  
Format of the SWI-Prolog executable, e.g. `elf` for when `swipl` is an ELF binary file.

**engines**(`bool`)  
True if engines are supported. This is always the case on the multi-threaded versions. Engines may be enabled on single threaded versions of SWI-Prolog.

**exit_status**(`integer`)  
Set by [halt/1](toplevel.html#halt/1) to its argument, making the exit status available to hooks registered with [at_halt/1](consulting.html#at_halt/1).

**file_name_case_handling**(`atom`, changeable)  
This flag defines how Prolog handles the case of file names. The flag is used for case normalization and to determine whether two names refer to the same file.^(bugNote that file name case handling is typically a properly of the filesystem, while Prolog only has a global flag to determine its file handling.) It has one of the following values:

**case_sensitive**  
The filesystem is fully case sensitive. Prolog does not perform any case modification or case insensitive matching. This is the default on Unix systems.

**case_preserving**  
The filesystem is case insensitive, but it preserves the case with which the user has created a file. This is the default on Windows systems.

**case_insensitive**  
The filesystem doesn't store or match case. In this scenario Prolog maps all file names to lower case.

**file_name_variables**(`bool`, changeable)  
If `true` (default `false`), expand `$\arg{varname}` and `~` in arguments of built-in predicates that accept a file name ([open/3](IO.html#open/3), [exists_file/1](files.html#exists_file/1), [access_file/2](files.html#access_file/2), etc.). The predicate [expand_file_name/2](files.html#expand_file_name/2) can be used to expand environment variables and wildcard patterns. This Prolog flag is intended for backward compatibility with older versions of SWI-Prolog.

**file_search_cache_time**(`number`, changeable)  
Time in seconds for which search results from [absolute_file_name/3](files.html#absolute_file_name/3) are cached. Within this time limit, the system will first check that the old search result satisfies the conditions. Default is 10 seconds, which typically avoids most repetitive searches for (library) files during compilation. Setting this value to 0 (zero) disables the cache.

**float_max**(`float`)  
The biggest representable floating point number.

**float_max_integer**(`float`)  
The highest integer that can be represented precisely as a floating point number.

**float_min**(`float`)  
The smallest representable floating point number above 0.0. See also [nexttoward/2](arith.html#f-nexttoward/2).

**float_overflow**(`atom`, changeable)  
One of `error` (default) or `infinity`. The first is ISO compliant. Using `infinity`, floating point overflow is mapped to positive or negative `Inf`. See [section 4.27.2.4](arith.html#sec:4.27.2.4). This flag also affects [read_term/3](termrw.html#read_term/3) and friends, causing them to read too large floating point number as infinity.

**float_rounding**(`atom`, changeable)  
Defines how arithmetic rounds to a float. Defined values are `to_nearest` (default), `to_positive`, `to_negative` or `to_zero`. For most scenarios the function [roundtoward/2](arith.html#f-roundtoward/2) provides a safer and faster alternative.

**float_undefined**(`atom`, changeable)  
One of `error` (default) or `nan`. The first is ISO compliant. Using `nan`, undefined operations such as `sqrt(-2.0)` is mapped to `NaN`. See [section 4.27.2.4](arith.html#sec:4.27.2.4).

**float_underflow**(`atom`, changeable)  
One of `error` or `ignore` (default). The second is ISO compliant, binding the result to 0.0.

**float_zero_div**(`atom`, changeable)  
One of `error` (default) or `infinity`. The first is ISO compliant. Using `infinity`, division by 0.0 is mapped to positive or negative `Inf`. See [section 4.27.2.4](arith.html#sec:4.27.2.4).

**gc**(`bool`, changeable)  
If true (default), the garbage collector is active. If false, neither garbage collection, nor stack shifts will take place, even not on explicit request. May be changed.

**gc_thread**(`bool`)  
If `true` (default if threading is enabled), atom and clause garbage collection are executed in a separate thread with the *alias* `gc`. Otherwise the thread that detected sufficient garbage executes the garbage collector. As running these global collectors may take relatively long, using a separate thread improves real time behaviour. The `gc` thread can be controlled using [set_prolog_gc_thread/1](memory.html#set_prolog_gc_thread/1), which either enables the gc thread or kills the gc thread and waits for it to die.

**generate_debug_info**(`bool`, changeable)  
If `true` (default) generate code that can be debugged using [trace/0](debugger.html#trace/0), [spy/1](debugger.html#spy/1), etc. Can be set to `false` using the **--no-debug**. This flag is scoped within a source file. Many of the libraries have `:- set_prolog_flag(generate_debug_info, false)` to hide their details from a normal trace.^(22In the current implementation this only causes a flag to be set on the predicate that causes children to be hidden from the debugger. The name anticipates further changes to the compiler.)

**gmp_version**(`integer`)  
If Prolog is linked with GMP, this flag gives the major version of the GMP library used. See also [section 12.4.11](foreigninclude.html#sec:12.4.11). This flag is not present when linked to [LibBF](https://bellard.org/libbf/). Use non-existence of the Prolog flag [bounded](flags.html#flag:bounded) to test for big integer and rational number support.

**gui**(`bool`)  
Set to `true` if XPCE is around and can be used for graphics.

**halt_grace_time**(`float`, changeable)  
Time [halt/1](toplevel.html#halt/1) waits for other threads to die gracefully. Default is 1 second.

**heartbeat**(`integer`, changeable)  
If not zero, call prolog:heartbeat/0 every `N` inferences. `N` is rounded to a multiple of 16.

**home**(`atom`)  
SWI-Prolog's notion of the home directory. SWI-Prolog uses its home directory to find its startup file as `<``home``>/boot.prc` and to find its library as `<``home``>/library`. Some installations may put architecture independent files in a *shared home* and also define [shared_home](flags.html#flag:shared_home). System files can be found using [absolute_file_name/3](files.html#absolute_file_name/3) as `swi(file)`. See [file_search_path/2](consulting.html#file_search_path/2).

**integer_rounding_function**(`down,toward_zero`)  
ISO Prolog flag describing rounding by `//` and `rem` arithmetic functions. Value depends on the C compiler used.

**iso**(`bool`, changeable)  
Include some weird ISO compatibility that is incompatible with normal SWI-Prolog behaviour. Currently it has the following effect:

- The `/``/2` (float division) *always* returns a float, even if applied to integers that can be divided.
- In the standard order of terms (see [section 4.6.1](compare.html#sec:4.6.1)), all floats are before all integers.
- [atom_length/2](manipatom.html#atom_length/2) yields a type error if the first argument is a number.
- [clause/\[2,3\]](examineprog.html#clause/2) raises a permission error when accessing static predicates.
- [abolish/\[1,2\]](db.html#abolish/1) raises a permission error when accessing static predicates.
- Syntax is closer to the ISO standard:
  - Within functional notation and list notation terms must have priority below 1000. That means that rules and control constructs appearing as arguments need bracketing. A term like `[a :- b, c].` must now be disambiguated to mean `[(a :- b), c].` or `[(a :- b, c)].`
  - Operators appearing as operands must be bracketed. Instead of `X == -, true.` write `X == (-), true.` Currently, this is not entirely enforced.
  - Backslash-escaped newlines are interpreted according to the ISO standard. See [section 2.15.1.3](syntax.html#sec:2.15.1.3).

**large_files**(`bool`)  
If present and `true`, SWI-Prolog has been compiled with *large file support* (LFS) and is capable of accessing files larger than 2GB. This flag is always `true` on 64-bit hardware and true on 32-bit hardware if the configuration detected support for LFS. Note that it may still be the case that the *file system* on which a particular file resides puts limits on the file size.

**last_call_optimisation**(`bool`, changeable)  
Determines whether or not last-call optimisation is enabled. Normally the value of this flag is the negation of the [debug](flags.html#flag:debug) flag. As programs may run out of stack if last-call optimisation is omitted, it is sometimes necessary to enable it during debugging.

**libswipl**(`atom`, changeable)  
Path where the SWI-Prolog shared library `libswipl`, the SWI-Prolog shared object that provides Prolog, resides. On some systems this can be determined reliably from the running system. On these systems the flag is *read-only*. On other systems it is the configured target installation location and thus this value can be wrong if the installation has been relocated. As we do not have a cross-platform reliable way to compute this path the flag is read-write on such platforms.^(23When running from the build environment, this flag is adjusted to reflect the location in the build tree.)

Currently, this flag is reliable on Windows and POSIX systems providing the **dladdr()** function. This function is provided on Linux and MacOS.

**linux**(`bool`)  
If present and `true`, the OS is some form of Linux. See also [unix](flags.html#flag:unix).

**malloc**(`atom`)  
Set after a successful identification of the used **malloc()** implementation. Currently possibly values are `tcmalloc` and `ptmalloc`. See [section 4.43.2](memory.html#sec:4.43.2) for details.

**max_answers_for_subgoal**(`integer`, changeable)  
Limit the number of answers in a table. The atom `infinite` clears the flag. By default this flag is not defined. See [section 7.11](tabling-restraints.html#sec:7.11) for details.

**max_answers_for_subgoal_action**(`atom`, changeable)  
The action taken when a table reaches the number of answers specified in [max_answers_for_subgoal](flags.html#flag:max_answers_for_subgoal). Supported values are `bounded_rationality`, `error` (default) or `suspend`.

**max_arity**(`unbounded`)  
ISO Prolog flag describing there is no maximum arity to compound terms.

**max_char_code**(`integer`)  
Highest (Unicode) code point that is supported. SWI-Prolog supports all Unicode code points from 0 (zero) up to and including the value of this flag. Currently `0xffff` on Windows (UCS-2) and `0x10ffff` on most other platforms.

**max_integer**(`integer`)  
Maximum integer value if integers are *bounded*. See also the flag [bounded](flags.html#flag:bounded) and [section 4.27.2.1](arith.html#sec:4.27.2.1).

**max_integer_size**(`integer`, changeable)  
When this tripwire is set, memory allocation on behalf of big integers and rational numbers is limited to given number of bytes. The minimum value is 1,000. When unset, the allocation limit is determined by the stack limit as we cannot represent larger numbers or **malloc()** failures. Notably services that may process arbitrary arithmetic expressions on behalf of a client may set this limit to avoid resource exhaustion.

**max_procedure_arity**(`integer`)  
Maximum arity for a predicate. An attempt to define or call such a predicate results in a `representation_error(max_procedure_arity)` exception. Currently set to 1024.

**max_rational_size**(`integer`, changeable)  
Limit the size in bytes for rational numbers. This *tripwire* can be used to identify cases where setting the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals) to `true` creates excessively big rational numbers and, if precision is not required, one should use floating point arithmetic. Note that rationals are also implicitly limited by the Prolog flag [max_integer_size](flags.html#flag:max_integer_size).

**max_rational_size_action**(`atom`, changeable)  
Action when the [max_rational_size](flags.html#flag:max_rational_size) tripwire is exceeded. Possible values are `error` (default), which throws a tripwire resource error and `float`, which converts the rational number into a floating point number. Note that rational numbers may exceed the range for floating point numbers.

**max_table_answer_size**(`integer`, changeable)  
Limit the size of an answer substitution for tabling. The atom `infinite` clears the flag. By default this flag is not defined. See [section 7.11](tabling-restraints.html#sec:7.11) for details.

**max_table_answer_size_action**(`atom`, changeable)  
The action taken if an answer substitution larger than [max_table_answer_size](flags.html#flag:max_table_answer_size) is added to a table. Supported values are `error` (default), `bounded_rationality`, `suspend` and `fail`.

**max_table_subgoal_size**(`integer`, changeable)  
Limit the size of a goal term accessing a table. The atom `infinite` clears the flag. By default this flag is not defined. See [section 7.11](tabling-restraints.html#sec:7.11) for details.

**max_table_subgoal_size_action**(`atom`, changeable)  
The action taken if a tabled goal exceeds [max_table_subgoal_size](flags.html#flag:max_table_subgoal_size). Supported values are `error` (default), `abstract` and `suspend`.

**max_tagged_integer**(`integer`)  
Maximum integer value represented as a‘tagged’value. Tagged integers require one word storage. Larger integers are represented as‘indirect data’and require significantly more space.

**message_context**(`list(atom)`, changeable)  
Context information to add to messages of the levels `error` and `warning`. The list may contain the elements `thread` to add the thread that generates the message to the message, `time` or `time(Format)` to add a time stamp. The default time format is `%T.%3f`. The default is `[thread]`. See also [format_time/3](system.html#format_time/3) and [print_message/2](printmsg.html#print_message/2).

**min_integer**(`integer`)  
Minimum integer value if integers are *bounded*. See also the flag [bounded](flags.html#flag:bounded) and [section 4.27.2.1](arith.html#sec:4.27.2.1).

**min_tagged_integer**(`integer`)  
Start of the tagged-integer value range.

**mitigate_spectre**(`bool`, changeable)  
When `true` (default `false`), enforce mitigation against the [Spectre](https://en.wikipedia.org/wiki/Spectre_(security_vulnerability)) timing-based security vulnerability. Spectre based attacks can extract information from memory owned by the process that should remain invisible, such as passwords or the private key of a web server. The attacks work by causing speculative access to sensitive data, and leaking the data via side-channels such as differences in the duration of successive instructions. An example of a potentially vulnerable application is [SWISH](https://swish.swi-prolog.org). SWISH allows users to run Prolog code while the swish server must protect the privacy of other users as well as its HTTPS private keys, cookies and passwords.

Currently, enabling this flag reduces the resolution of [get_time/1](system.html#get_time/1) and [statistics/2](builtin-statistics.html#statistics/2) CPU time to `20μs`.

**WARNING**: Although a coarser timer makes a successful attack of this type harder, it does not reliably prevent such attacks in general. Full mitigation may require compiler support to disable speculative access to sensitive data.

**msys2**(`bool`)  
If present, SWI-Prolog is the MS-Windows version running under a [MSYS2](https://www.msys2.org/) shell.

**occurs_check**(`atom`, changeable)  
This flag controls unification that creates an infinite tree (also called *cyclic term*) and can have three values. Using `false` (default), unification succeeds, creating an infinite tree. Using `true`, unification behaves as [unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2), failing silently. Using `error`, an attempt to create a cyclic term results in an `occurs_check` exception. The latter is intended for debugging unintentional creations of cyclic terms. Note that this flag is a global flag modifying fundamental behaviour of Prolog. Changing the flag from its default may cause libraries to stop functioning properly.

**on_error**(`atom`, changeable)  
Determines how to act on an error printed using [print_message/2](printmsg.html#print_message/2), i.e., an error that is reported to the user. The possible values are `print` (default), `status` and `halt`. Using `halt` the process halts immediately with status 1. Otherwise execution continues. Using `status` [halt/0](toplevel.html#halt/0) exits with status 1 if one or more errors were printed by the process. In *compile* mode (see **-c**) the default is `status`. This flag can be set from the commandline using **--on-error**. See also [section 4.3.2.1](consulting.html#sec:4.3.2.1).

**on_warning**(`atom`, changeable)  
As [on_error](flags.html#flag:on_error), but for warnings. The default is always `print`. The commandline option is **--on-warning**.

**open_shared_object**(`bool`)  
If true, [open_shared_object/2](foreignlink.html#open_shared_object/2) and friends are implemented, providing access to shared libraries (`.so` files) or dynamic link libraries (`.DLL` files).

**optimise**(`bool`, changeable)  
If `true`, compile in optimised mode. The initial value is `true` if Prolog was started with the **-O** command line option. The [optimise](flags.html#flag:optimise) flag is scoped to a source file.

Currently optimised compilation implies compilation of arithmetic, and deletion of redundant [true/0](control.html#true/0) that may result from [expand_goal/2](consulting.html#expand_goal/2).

Later versions might imply various other optimisations such as integrating small predicates into their callers, eliminating constant expressions and other predictable constructs. Source code optimisation is never applied to predicates that are declared dynamic (see [dynamic/1](dynamic.html#dynamic/1)).

**optimise_unify**(`bool`, changeable)  
If `true` (default), allow the compiler to (re)move explicit unification calls ([=/2](compare.html#=/2)). While this behaviour can significantly improve performance, it is not yet handled properly by the source-level debugger. See [section 2.17.3](jitindex.html#sec:2.17.3).

**os_argv**(`list`, changeable)  
List is a list of atoms representing the command line arguments used to invoke SWI-Prolog. Please note that **all** arguments are included in the list returned. See [argv](flags.html#flag:argv) to get the application options.

**packs**(`bool`)  
If `true`, extension packs (add-ons) are attached. Can be set to `false` using the **--no-packs**.

**path_max**(`integer`)  
Maximum length of a file pathname as reported by the OS. This length does typically not directly define the number of characters in the file name. The actual limit may be shorter due to jargonencoding (e.g., on POSIX systems it typically defines the length limit of the (often) UTF-8 encoded name). The underlying file system may impose additional limits.

**path_sep**(`atom`)  
Separator for file search paths such as the environment variable `PATH` for the OS. Normally `:`, but `;` on Windows.

**pid**(`int`)  
Process identifier of the running Prolog process. Existence of this flag is implementation-defined.

**pipe**(`bool`, changeable)  
If true, `open(pipe(command), mode, Stream)`, etc. are supported. Can be changed to disable the use of pipes in applications testing this feature. Not recommended.

**portable_vmi**(`bool`, changeable)  
If `true` (default), generate `.qlf` files and saved states that run both on 32 bit and 64-bit hardware. If `false`, some optimized virtual machine instructions are only used if the integer argument is within the range of a tagged integer for 32-bit machines.

**posix_shell**(`atom`, changeable)  
Path to a POSIX compatible shell. This default is typically `/bin/sh`. This flag is used by [shell/1](system.html#shell/1) and [qsave_program/2](saved-states.html#qsave_program/2).

**prefer_rationals**(`bool`, changeable)  
Only provided if the system is compiled with unbounded and rational arithmetic support (see [bounded](flags.html#flag:bounded)). If `true`, prefer arithmetic to produce rational numbers over floats. This implies:

- Division ([//2](arith.html#f-//2)) of two integers produces a rational number.
- Power ([^/2](arith.html#f-%5E/2)) of two integers produces a rational number, *also* if the second operand is a negative number. For example, `2^(-2)` evaluates to `1/4`.

Using `true` can create excessively large rational numbers. The Prolog flag [max_rational_size](flags.html#flag:max_rational_size) can be used to detect and act on this *tripwire*.

If `false`, rational numbers can only be created using the functions [rational/1](arith.html#f-rational/1), [rationalize/1](arith.html#f-rationalize/1) and [rdiv/2](arith.html#f-rdiv/2) or by reading them. See also [rational_syntax](flags.html#flag:rational_syntax), [section 2.15.1.6](syntax.html#sec:2.15.1.6) and [section 4.27.2.2](arith.html#sec:4.27.2.2).

The current default is `false`. We consider changing this to `true` in the future. Users are strongly encouraged to set this flag to `true` and report issues this may cause.

**print_write_options**(`term`, changeable)  
Specifies the options for [write_term/2](termrw.html#write_term/2) used by [print/1](termrw.html#print/1) and [print/2](termrw.html#print/2).

**prompt_alternatives_on**(`atom`, changeable)  
Determines prompting for alternatives in the Prolog top level. Default is `determinism`, which implies the system prompts for alternatives if the goal succeeded while leaving choice points. Many classical Prolog systems behave as `groundness`: they prompt for alternatives if and only if the query contains variables.

**protect_static_code**(`bool`, changeable)  
If `true` (default `false`), [clause/2](examineprog.html#clause/2) does not operate on static code, providing some basic protection from hackers that wish to list the static code of your Prolog program. Once the flag is `true`, it cannot be changed back to `false`. Protection is default in ISO mode (see Prolog flag [iso](flags.html#flag:iso)). Note that many parts of the development environment require [clause/2](examineprog.html#clause/2) to work on static code, and enabling this flag should thus only be used for production code.

**qcompile**(`atom`, changeable)  
This option provides the default for the `qcompile(+Atom)` option of [load_files/2](consulting.html#load_files/2).

**rational_syntax**(`atom`, changeable)  
Determines the read and write syntax for rational numbers. Possible values are `natural` (e.g., `1/3`) or `compatibility` (e.g., `1r3`). The `compatibility` syntax is always accepted. This flag is module sensitive.

The default for this flag is currently `compatibility`, which reads and writes rational numbers as e.g., `1r3`.^(24There is still some discussion on the separating character. See [section 2.15.1.6](syntax.html#sec:2.15.1.6).) We will consider `natural` as a default in the future. Users are strongly encouraged to set this flag to `natural` and report issues this may cause.

**rationals**(`atom`)  
This flag is present and has the value `true` if the system supports rational numbers. For SWI-Prolog this flag is always set if the flag [bounded](flags.html#flag:bounded) is `false`.

**readline**(`atom`, changeable)  
Specifies which form of command line editing is provided. Possible values are below.

**false**  
No command line editing is available.

**editline**  
The library `library(editline)` is loaded, providing line editing based on the BSD libedit. This is the default if `library(editline)` is available.

**report_error**(`bool`, changeable)  
If `true`, print error messages; otherwise suppress them. May be changed. See also the [debug_on_error](flags.html#flag:debug_on_error) Prolog flag. Default is `true`, except for the runtime version.

**resource_database**(`atom`)  
Set to the absolute filename of the attached state. Typically this is the file `boot32.prc`, the file specified with **-x** or the running executable. See also [resource/3](program-resources.html#resource/3).

**runtime**(`bool`)  
If present and `true`, SWI-Prolog is compiled with -DO_RUNTIME, disabling various useful development features (currently the tracer and profiler).

**sandboxed_load**(`bool`, changeable)  
If `true` (default `false`), [load_files/2](consulting.html#load_files/2) calls hooks to allow library(sandbox) to verify the safety of directives.

**saved_program**(`bool`)  
If present and `true`, Prolog has been started from a state saved with [qsave_program/\[1,2\]](saved-states.html#qsave_program/1).

**shared_home**(`atom`)  
Indicates that part of the SWI-Prolog system files are installed in `<``prefix``>/share/swipl` instead of in the home at the `<``prefix``>/lib/swipl`. This flag indicates the location of this *shared home* and the directory is added to the file search path `swi`. See [file_search_path/2](consulting.html#file_search_path/2) and the flag [home](flags.html#flag:home).

**shared_object_extension**(`atom`)  
Extension used by the operating system for shared objects. `.so` for most Unix systems and `.dll` for Windows. Used for locating files using the `file_type` `executable`. See also [absolute_file_name/3](files.html#absolute_file_name/3).

**shared_object_search_path**(`atom`)  
Name of the environment variable used by the system to search for shared objects.

**shared_table_space**(`integer`, changeable)  
Space reserved for storing shared answer tables. See [section 7.9](tabling-shared.html#sec:7.9) and the Prolog flag [table_space](flags.html#flag:table_space).

**shift_check**(`bool`, changeable)  
When `true` (default `false`), check for suspicious delimited continuations captured by [shift_for_copy/1](delcont.html#shift_for_copy/1).

**signals**(`bool`)  
Determine whether Prolog is handling signals (software interrupts). This flag is `false` if the hosting OS does not support signal handling or the command line option **--no-signals** is active. See [section 12.4.25.1](foreigninclude.html#sec:12.4.25.1) for details.

**source**(`bool`, changeable)  
If `true`, ignore `.qlf` files if there is a corresponding `.pl` file. The provided `.qlf` files in the library are compiled with optimization enabled (see [optimise](flags.html#flag:optimise)), macro expansion enabled (see `library(apply_macros)`) and [debug/3](debug.html#debug/3) and [assertion/1](debug.html#assertion/1) statements removed. Using this flag loads the source, providing better support for debugging. If a debugging session may benefit from better access to the debugging facilities of the libraries, either set this Prolog flag at the start of the load file for your program or start Prolog as

``` code
swipl -Dsource [option ...] myfile.pl ...
```

**source_search_working_directory**(`bool`, changeable)  
If set to `true`, loading a relative file name from source code searches relative to the location of the source file as well as relative to the working directory. Searching relative to the working directory is deprecated and a warning is printed if the file is found this way. Future versions are likely to change the default to `false`.^(25Searching the working directory was supported up to version 9.3.8. Version 9.3.9 disabled this and version 9.3.10 re-enables it with a warning.)

**stack_limit**(`int`, changeable)  
Limits the combined sizes of the Prolog stacks for the current thread. See also **--stack-limit** and [section 2.19.1](limits.html#sec:2.19.1).

**stream_type_check**(`atom`, changeable)  
Defines whether and how strictly the system validates that byte I/O should not be applied to text streams and text I/O should not be applied to binary streams. Values are `false` (no checking), `true` (full checking) and `loose`. Using checking mode `loose` (default), the system accepts byte I/O from text stream that use ISO Latin-1 encoding and accepts writing text to binary streams.

**string_stack_tripwire**(`int`, changeable)  
Maintenance for foreign language string management. Prints a warning if the string stack depth hits the tripwire value. See [section 12.4.14](foreigninclude.html#sec:12.4.14) for details.

**system_thread_id**(`int`)  
Available in multithreaded version (see [section 10](threads.html#sec:10)) where the operating system provides system-wide integer thread identifiers. The integer is the thread identifier used by the operating system for the calling thread. On Linux systems this is the PID of the thread.

**table_incremental**(`bool`, changeable)  
Set the default for whether to use incremental tabling or not. Initially set to `false`. See [table/1](tabling-preds.html#table/1).

**table_shared**(`bool`, changeable)  
Set the default for whether to use shared tabling or not. Initially set to `false`. See [table/1](tabling-preds.html#table/1).

**table_space**(`integer`, changeable)  
Space reserved for storing answer tables for *tabled predicates* (see [table/1](tabling-preds.html#table/1)).^(bugCurrently only counts the space occupied by the nodes in the answer tries.) When exceeded a `resource_error(table_space)` exception is raised.

**table_subsumptive**(`bool`, changeable)  
Set the default choice between *variant* tabling and *subsumptive* tabling. Initially set to `false`. See [table/1](tabling-preds.html#table/1).

**threads**(`bool`, changeable)  
True when threads are supported. If the system is compiled without thread support the value is `false` and read-only. Otherwise the value is `true` unless the system was started with the **--no-threads**. Threading may be disabled only if no threads are running. See also the [gc_thread](flags.html#flag:gc_thread) flag.

**timezone**(`integer`)  
Offset in seconds west of GMT of the current time zone. Set at initialization time from the `timezone` variable associated with the POSIX **tzset()** function. See also [format_time/3](system.html#format_time/3).

**tmp_dir**(`atom`, changeable)  
Path to the temporary directory. initialised from the environment variable `TMP` or `TEMP` in windows. If this variable is not defined a default is used. This default is typically `/tmp` or `c:/temp` in windows.

**toplevel_goal**(`term`, changeable)  
Defines the goal that is executed after running the initialization goals and entry point (see **-g**, [initialization/2](consulting.html#initialization/2) and [section 2.11.1.1](compilation.html#sec:2.11.1.1). The initial value is `default`, starting a normal interactive session. This value may be changed using the command line option **-t**. The explicit value `prolog` is equivalent to `default`. If `initialization(Goal,main)` is used and the toplevel is `default`, the toplevel is set to `halt` (see [halt/0](toplevel.html#halt/0)).

**toplevel_list_wfs_residual_program**(`bool`, changeable)  
If `true` (default) and the answer is *undefined* according to the Well Founded Semantics (see [section 7.6](WFS.html#sec:7.6)), list the *residual program* before the answer. Otherwise the answer terminated with **undefined**. See also [undefined/0](WFS.html#undefined/0).

**toplevel_mode**(`atom`, changeable)  
If `backtracking` (default), the toplevel backtracks after completing a query. If `recursive`, the toplevel is implemented as a recursive loop. This implies that global variables set using [b_setval/2](gvar.html#b_setval/2) are maintained between queries. In *recursive* mode, answers to toplevel variables (see [section 2.9](topvars.html#sec:2.9)) are kept in backtrackable global variables and thus **not copied**. In *backtracking* mode answers to toplevel variables are kept in the recorded database (see [section 4.14.2](db.html#sec:4.14.2)).

The recursive mode has been added for interactive usage of CHR (see [section 9](chr.html#sec:9)),^(26Suggested by Falco Nogatz) which maintains the global constraint store in backtrackable global variables.

**toplevel_name_variables**(`bool`, changeable)  
If `true` (default), give names to variables at the toplevel instead of printing them as `_NNN`. The variables are named `_A`, `_B`, ... Variables that appear only once (singletons) are printed as `_`.

**toplevel_print_anon**(`bool`, changeable)  
If `true`, top-level variables starting with an underscore (`_`) are printed normally. If `false` (default) the binding of such variables are omitted from the answer. This may be used to hide bindings in complex queries from the top level. For example, the binding for `_List` below is not printed.

``` code
?- numlist(1,1 000 000,_List), sum_list(_List, Sum).
Sum = 500000500000.
```

**toplevel_print_factorized**(`bool`, changeable)  
If `true` (default `false`) show the internal sharing of subterms in the answer substitution. The example below reveals internal sharing of leaf nodes in *red-black trees* as implemented by the `library(rbtrees)` predicate rb_new/1 :

``` code
?- set_prolog_flag(toplevel_print_factorized, true).
?- rb_new(X).
X = t(_S1, _S1), % where
    _S1 = black('', _G387, _G388, '').
```

If this flag is `false`, the `% where` notation is still used to indicate cycles as illustrated below. This example also shows that the implementation reveals the internal cycle length, and *not* the minimal cycle length. Cycles of different length are indistinguishable in Prolog (as illustrated by `S == R`).

``` code
?- S = s(S), R = s(s(R)), S == R.
S = s(S),
R = s(s(R)).
```

**toplevel_prompt**(`atom`, changeable)  
Define the prompt that is used by the interactive top level. The following `~` (tilde) sequences are replaced:

|  |  |
|----|----|
| `~`m | *Type in* module if not `user` (see [module/1](mtoplevel.html#module/1)) |
| `~`l | *Break level* if not 0 (see [break/0](toplevel.html#break/0)) |
| `~`d | *Debugging state* if not normal execution (see [debug/0](debugger.html#debug/0), [trace/0](debugger.html#trace/0)) |
| `~`! | *History event* if history is enabled (see flag **history**) |

**toplevel_residue_vars**(`bool`, changeable)  
When `true` (default `false`), print residual variables as detected by [call_residue_vars/2](coroutining.html#call_residue_vars/2) that do not appear in the bindings returned by the goal.

**toplevel_thread**(`bool`, changeable)  
When `true`, this thread is running a toplevel REPL loop. See [prolog/0](toplevel.html#prolog/0).

**toplevel_var_size**(`int`, changeable)  
Maximum size counted in literals of a term returned as a binding for a variable in a top-level query that is saved for re-use using the `$` variable reference. When 0 (zero), the variable recording and reuse is disabled. See [section 2.9](topvars.html#sec:2.9).

**trace_gc**(`bool`, changeable)  
If `true` (default `false`), garbage collections and stack-shifts will be reported on the terminal. May be changed. Values are reported in bytes as `G`+`T`, where `G` is the global stack value and `T` the trail stack value.‘Gained’describes the number of bytes reclaimed.‘used’the number of bytes on the stack after GC and‘free’the number of bytes allocated, but not in use. Below is an example output.

``` code
% GC: gained 236,416+163,424 in 0.00 sec;
      used 13,448+5,808; free 72,568+47,440
```

**traditional**(`bool`)  
Available in SWI-Prolog version 7. If `true`,‘traditional’mode has been selected using **--traditional**. Notice that some SWI7 features, like the functional notation on dicts, do not work in this mode. See also [section 5](extensions.html#sec:5).

**tty_control**(`bool`, changeable)  
Determines whether the terminal is switched to raw mode for [get_single_char/1](chario.html#get_single_char/1), which also reads the user actions for the trace. May be set. If this flag is `false` at startup, command line editing is disabled. See also the **--no-tty** command line option.

**unix**(`bool`)  
If present and `true`, the operating system is some version of Unix. Defined if the C compiler used to compile this version of SWI-Prolog either defines `__unix__` or `unix`. On other systems this flag is not available. See also [linux](flags.html#flag:linux), [apple](flags.html#flag:apple) and [windows](flags.html#flag:windows).

**unknown**(`fail,warning,error`, changeable)  
Determines the behaviour if an undefined procedure is encountered. If `fail`, the predicate fails silently. If `warn`, a warning is printed, and execution continues as if the predicate was not defined, and if `error` (default), an `existence_error` exception is raised. This flag is local to each module and inherited from the module's *import-module*. Using default setup, this implies that normal modules inherit the flag from `user`, which in turn inherit the value `error` from `system`. The user may change the flag for module `user` to change the default for all application modules or for a specific module. It is strongly advised to keep the `error` default and use [dynamic/1](dynamic.html#dynamic/1) and/or [multifile/1](dynamic.html#multifile/1) to specify possible non-existence of a predicate.

**unknown_option**(`ignore,warning,error`, changeable)  
Determines the behaviour if a predicate that processes an option list is passed an option that is not understood by the predicate. The ISO standard dictates raising a `domain_error` exception. This is considered impractical as it makes it hard to write portable code if different Prolog systems support different options and it makes it hard to write predicates that process options and pass some of the options to one predicate and others to some other predicate. For example, a predicate reading a file to a list of terms must distribute options to [open/4](IO.html#open/4) and [read_term/3](termrw.html#read_term/3). SWI-Prolog has always ignored unknown options unless in ISO mode (see the [iso](flags.html#flag:iso) flag). This flag provides full control over how options are processed.

**unload_foreign_libraries**(`bool`, changeable)  
If `true` (default `false`), unload all loaded foreign libraries. Default is `false` because modern OSes reclaim the resources anyway and unloading the foreign code may cause registered hooks to point to no longer existing data or code.

**user_flags**(`Atom`, changeable)  
Define the behaviour of [set_prolog_flag/2](flags.html#set_prolog_flag/2) if the flag is not known. Values are `silent`, `warning` and `error`. The first two create the flag on-the-fly, where `warning` prints a message. The value `error` is consistent with ISO: it raises an existence error and does not create the flag. See also [create_prolog_flag/3](flags.html#create_prolog_flag/3). The default is `silent`, but future versions may change that. Developers are encouraged to use another value and ensure proper use of [create_prolog_flag/3](flags.html#create_prolog_flag/3) to create flags for their library.

**var_prefix**(`bool`, changeable)  
If `true` (default `false`), variables must start with an underscore (`_`). May be changed. This flag is local to the module in which it is changed. See [section 2.15.1.8](syntax.html#sec:2.15.1.8).

**var_tag**(`Atom`, changeable)  
This flag controls the interpretation of `Tag`{...} if `Tag` is unbound. Possible values are:

**dict**  
Read as a dict (see [section 5.4](bidicts.html#sec:5.4)) with unbound tag. This is the current default. The use of unbound tags is deprecated. In due time the default will be changed, ultimately to `attvar`. See [section 5.4.1](bidicts.html#sec:5.4.1).

**attvar**(`R`)  
ead `Var{...}` as an attributed variable. This is likely to become the future default.

**`#`**  
Read `Var{...}` as `#{...}`. This flag allows to quickly assess the impact on your code of using `#{...}` rather than dicts with unbound tags.

**warning**  
As `#`, but prints a warning.

**error**  
Considers `Var{...}` a syntax error.

**verbose**(`atom`, changeable)  
This flag is used by [print_message/2](printmsg.html#print_message/2). If its value is `silent`, messages of type `informational` and `banner` are suppressed. The **-q** switches the value from the initial `normal` to `silent`.

**verbose_autoload**(`bool`, changeable)  
If `true` the normal consult message will be printed if a library is autoloaded. By default this message is suppressed. Intended to be used for debugging purposes.

**verbose_file_search**(`bool`, changeable)  
If `true` (default `false`), print messages indicating the progress of [absolute_file_name/\[2,3\]](files.html#absolute_file_name/2) in locating files. Intended for debugging complicated file-search paths. See also [file_search_path/2](consulting.html#file_search_path/2).

**verbose_load**(`atom`, changeable)  
Determines messages printed for loading (compiling) Prolog files. Current values are `full` (print a message at the start and end of each file loaded), `normal` (print a message at the end of each file loaded), `brief` (print a message at end of loading the toplevel file), and `silent` (no messages are printed, default). The value of this flag is normally controlled by the option `silent(Bool)` provided by [load_files/2](consulting.html#load_files/2).

**version**(`integer`)  
The version identifier is an integer with value:

> `10000 × ``Major`` + 100 × ``Minor`` + ``Patch`

**version_data**(`swi(Major, Minor, Patch, Extra)`)  
Part of the dialect compatibility layer; see also the Prolog flag [dialect](flags.html#flag:dialect) and [section C](dialect.html#sec:C). `Extra` provides platform-specific version information as a list. `Extra` is used for *tagged versions* such as “7.4.0-rc1” , in which case `Extra` contains a term `tag(rc1)`.

**version_git**(`atom`)  
Available if created from a git repository. See **git-describe** for details.

**vmi_builtin**(`bool`, changeable)  
Determines whether well known built-ins such as [true/0](control.html#true/0) or [atom/1](typetest.html#atom/1) are handled by their translation into virtual machine code. The default for this flag is `true`, unless debug mode is enabled. Setting this flag to `false` may improve other runtime instrumentation results. Note that optimized arithmetic (**-O**, see Prolog flag [optimise](flags.html#flag:optimise)) is currently not translated into a normal predicate call.

**warn_autoload**(`bool`, changeable)  
If `true` (default `false`), warn when autoloading predicates from a file that defines global term- or goal-expansion rules. These rules typically enhance performance or provide cleaner semantics and thus autoloading is not recommended. Future versions will enable this flag by default.

**warn_override_implicit_import**(`bool`, changeable)  
If `true` (default), a warning is printed if an implicitly imported predicate is clobbered by a local definition. See [use_module/1](import.html#use_module/1) for details.

**win_file_access_check**(`atom`, changeable)  
Controls the behaviour or [access_file/2](files.html#access_file/2) under Windows. There is no reliable way to check access to files and directories on Windows. This flag allows for switching between three alternative approximations.

**access**  
Use Windows **\_waccess()** function. This ignores ACLs (Access Control List) and thus may indicate that access is allowed while it is not.

**getfilesecurity**  
Use the Windows **GetFileSecurity()** function. This does not work on all file systems, but is probably the best choice on file systems that do support it, notably local NTFS volumes.

**openclose**  
Try to open the file and close it. This works reliable for files, but not for directories. Currently directories are checked using **\_waccess()**. This is the default.

**windows**(`bool`)  
If present and `true`, the operating system is an implementation of Microsoft Windows. This flag is only available on MS-Windows based versions. See also [unix](flags.html#flag:unix).

**wine_version**(`atom`)  
If present, SWI-Prolog is the MS-Windows version running under the [Wine](https://www.winehq.org/) emulator.

**write_attributes**(`atom`, changeable)  
Defines how [write/1](termrw.html#write/1) and friends write attributed variables. The option values are described with the `attributes` option of [write_term/2](termrw.html#write_term/2). Default is `ignore`.

**write_help_with_overstrike**(`bool`)  
Internal flag used by [help/1](online-help.html#help/1) when writing to a terminal. If present and `true` it prints bold and underlined text using *overstrike*.

**xdg**(`bool`, changeable)  
This flag defines whether or not the we follow the Free Desktop standard for application data and configuration files. The flag is `true` and read-only for non-Windows systems. On Windows systems the flag is `true` but read-write when compiled under *Conda* or `MSYS2` and not defined otherwise. On Windows, the search order is

**Flag is not defined**  
First search the Windows directories, then the XDG directories. This is the default for the Windows binaries.

**Flag is `true`**  
Only search the XDG directories.

**Flag is `false`**  
Only search the Windows directories.

**xpce**(`bool`)  
Available and set to `true` if the XPCE graphics system is loaded.

**xpce_version**(`atom`)  
Available and set to the version of the loaded XPCE system.

**xref**(`bool`, changeable)  
If `true`, source code is being read for *analysis* purposes such as cross-referencing. Otherwise (default) it is being read to be compiled. This flag is used at several places by [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2) hooks, notably if these hooks use side effects. See also the libraries `library(prolog_source)` and `library(prolog_xref)`.

\[ISO\]**set_prolog_flag**(`:Key, +Value`)  
Set the value of a Prolog flag. `Key` is an atom. If the flag is a system-defined flag that is not marked *changeable* above, an attempt to modify the flag yields a `permission_error`. If the provided `Value` does not match the type of the flag, a `type_error` is raised.

Some flags (e.g., [unknown](flags.html#flag:unknown)) are maintained on a per-module basis. The addressed module is determined by the `Key` *meta* argument.

In addition to ISO, SWI-Prolog allows for user-defined Prolog flags. New Prolog flags should be created using [create_prolog_flag/3](flags.html#create_prolog_flag/3). For historical reasons [set_prolog_flag/2](flags.html#set_prolog_flag/2) silently creates a Prolog flag if the Prolog flag [user_flags](flags.html#flag:user_flags) is `true` (default), i.e., [set_prolog_flag/2](flags.html#set_prolog_flag/2) behaves as if defined below.

``` code
set_prolog_flag(Key, Value) :-
    current_prolog_flag(Key, _),
    !,
    <set the flag>.
set_prolog_flag(Key, Value) :-
    current_prolog_flag(user_flags, true),
    !,
    create_prolog_flag(Key, Value, []).
set_prolog_flag(Key, _) :-
    existence_error(prolog_flag, Key).
```

\[YAP\]**create_prolog_flag**(`+Key, +Value, +Options`)  
Create a new Prolog flag. The ISO standard does not foresee creation of new flags, but many libraries introduce new flags. Prolog flags are particularly useful for managing mutating global settings in a threaded environment. Were predicates are either local to a thread or shared between all threads, a thread *inherits* its flags from the thread that creates it (see [thread_create/3](threadcreate.html#thread_create/3)) while modifications are local to the calling thread.

SWI-Prolog flags are typed. If the type is not explicitly defined using the `type(Type)` option (see below), the type is determined from the initial value. Defined types are `boolean` (if the initial value is one of `false`, `true`, `on` or `off`), `atom` if the initial value is any other atom, `integer` if the value is an integer that can be expressed as a 64-bit signed value. Any other initial value results in an untyped flag that can represent any valid Prolog term.

By default, the new flag is added to the global flag table, making the value available to all threads that have not set the flag. If the flag was already defined locally in the calling thread, the value is updated both for the calling thread and in the global flag table. See the `local(+Boolean)` option.

`Options` is a list of the options below. See also the Prolog flag [user_flags](flags.html#flag:user_flags).

**access**(`+Access`)  
Define access rights for the flag. Values are `read_write` and `read_only`. The default is `read_write`.

**type**(`+Atom`)  
Define a type restriction. Possible values are `boolean`, `atom`, `oneof(ListOfAtoms)`, `integer`, `float` and `term`. The default is determined from the initial value. Note that `term` restricts the term to be ground.

**keep**(`+Boolean`)  
If `true`, do not modify the flag if it already exists. Otherwise (default), this predicate behaves as [set_prolog_flag/2](flags.html#set_prolog_flag/2) if the flag already exists.

If the flag has a value, but this value is incompatible with the specified type, a warning is printed and the flag gets the value and type specified by this call to [create_prolog_flag/3](flags.html#create_prolog_flag/3).

**local**(`+Boolean`)  
If `true` (default `false`), if the flag does not exist, create it *only* in the calling thread. The flag is only visible in the calling thread and threads that inherit from the calling thread.

**warn_not_accessed**(`+Boolean`)  
If `true` and the flag is never read using [current_prolog_flag/2](flags.html#current_prolog_flag/2), print a warning. This option is used for options set using the commandline option `-D<flag>[=<value>]`.
