
## 4.3 Loading Prolog source files

This section deals with loading Prolog source files. A Prolog source file is a plain text file containing a Prolog program or part thereof. Prolog source files come in three flavours:

**A traditional**  
Prolog source file contains Prolog clauses and directives, but no *module declaration* (see [module/1](mtoplevel.html#module/1)). They are normally loaded using [consult/1](consulting.html#consult/1) or [ensure_loaded/1](consulting.html#ensure_loaded/1). Currently, a non-module file can only be loaded into a single module.^(54This limitation may be lifted in the future. Existing limitations in SWI-Prolog's source code administration make this non-trivial.)

**A module**  
Prolog source file starts with a module declaration. The subsequent Prolog code is loaded into the specified module, and only the *exported* predicates are made available to the context loading the module. Module files are normally loaded with [use_module/\[1,2\]](import.html#use_module/1). See [chapter 6](modules.html#sec:6) for details.

**An include**  
Prolog source file is loaded using the [include/1](consulting.html#include/1) directive, textually including Prolog text into another Prolog source. A file may be included into multiple source files and is typically used to share *declarations* such as multifile or dynamic between source files.

Prolog source files are located using [absolute_file_name/3](files.html#absolute_file_name/3) with the following options:

``` code
locate_prolog_file(Spec, Path) :-
        absolute_file_name(Spec,
                           [ file_type(prolog),
                             access(read)
                           ],
                           Path).
```

The `file_type(prolog)` option is used to determine the extension of the file using [prolog_file_type/2](consulting.html#prolog_file_type/2). The default extension is `.pl`. `Spec` allows for the *path alias* construct defined by [absolute_file_name/3](files.html#absolute_file_name/3). The most commonly used path alias is `library(LibraryFile)`. The example below loads the library file `ordsets.pl` (containing predicates for manipulating ordered sets).

``` code
:- use_module(library(ordsets)).
```

SWI-Prolog recognises grammar rules (DCG) as defined in [Clocksin & Melish, 1987](Bibliography.html#Clocksin:87). The user may define additional compilation of the source file by defining the dynamic multifile predicates [term_expansion/2](consulting.html#term_expansion/2), [term_expansion/4](consulting.html#term_expansion/4), [goal_expansion/2](consulting.html#goal_expansion/2) and [goal_expansion/4](consulting.html#goal_expansion/4). It is not allowed to use [assert/1](db.html#assert/1), [retract/1](db.html#retract/1) or any other database predicate in [term_expansion/2](consulting.html#term_expansion/2) other than for local computational purposes.^(55It does work for normal loading, but not for [qcompile/1](consulting.html#qcompile/1).) Code that needs to create additional clauses must use [compile_aux_clauses/1](consulting.html#compile_aux_clauses/1). See `library(library(apply_macros))` for an example.

A *directive* is an instruction to the compiler. Directives are used to set (predicate) properties (see [section 4.15](dynamic.html#sec:4.15)), set flags (see [set_prolog_flag/2](flags.html#set_prolog_flag/2)) and load files (this section). Directives are terms of the form `:-` \<`term`\>. . Here are some examples:

``` code
:- use_module(library(lists)).
:- dynamic
        store/2.                % Name, Value
```

The directive [initialization/1](consulting.html#initialization/1) can be used to run arbitrary Prolog goals. The specified goal is started *after* loading the file in which it appears has completed.

SWI-Prolog compiles code as it is read from the file, and directives are executed as *goals*. This implies that directives may call any predicate that has been defined before the point where the directive appears. It also accepts `?-` \<`term`\>. as a synonym.

SWI-Prolog does not have a separate reconsult/1 predicate. Reconsulting is implied automatically by the fact that a file is consulted which is already loaded.

Advanced topics are handled in subsequent sections: mutually dependent files ([section 4.3.2.2](consulting.html#sec:4.3.2.2)), multithreaded loading ([section 4.3.2.3](consulting.html#sec:4.3.2.3)) and reloading running code ([section 4.3.2](consulting.html#sec:4.3.2)).

The core of the family of loading predicates is [load_files/2](consulting.html#load_files/2). The predicates [consult/1](consulting.html#consult/1), [ensure_loaded/1](consulting.html#ensure_loaded/1), [use_module/1](import.html#use_module/1), [use_module/2](import.html#use_module/2) and [reexport/1](reexport.html#reexport/1) pass the file argument directly to [load_files/2](consulting.html#load_files/2) and pass additional options as expressed in the [table 4](consulting.html#tab:loadpreds):

|  |  |  |  |
|----|:--:|:--:|:--:|
| **Predicate** | **if** | **must_be_module** | **import** |
| [consult/1](consulting.html#consult/1) | `true` | `false` | all |
| [ensure_loaded/1](consulting.html#ensure_loaded/1) | `not_loaded` | `false` | all |
| [use_module/1](import.html#use_module/1) | `not_loaded` | `true` | all |
| [use_module/2](import.html#use_module/2) | `not_loaded` | `true` | specified |
| [reexport/1](reexport.html#reexport/1) | `not_loaded` | `true` | all |
| [reexport/2](reexport.html#reexport/2) | `not_loaded` | `true` | specified |

**Table 4 :** Properties of the file-loading predicates. The *import* column specifies what is imported if the loaded file is a module file.

**load_files**(`:Files`)  
Equivalent to `load_files(Files,[])`. Same as [consult/1](consulting.html#consult/1), See [load_files/2](consulting.html#load_files/2) for supported options.

**load_files**(`:Files, +Options`)  
The predicate [load_files/2](consulting.html#load_files/2) is the parent of all the other loading predicates except for [include/1](consulting.html#include/1). It currently supports a subset of the options of Quintus [load_files/2](consulting.html#load_files/2). `Files` is either a single source file or a list of source files. The specification for a source file is handed to [absolute_file_name/2](files.html#absolute_file_name/2). See this predicate for the supported expansions. `Options` is a list of options using the format `OptionName`(`OptionValue`).

The following options are currently supported:

**autoload**(`Bool`)  
If `true` (default `false`), indicate that this load is a *demand* load. This implies that, depending on the setting of the Prolog flag [verbose_autoload](flags.html#flag:verbose_autoload), the load action is printed at level `informational` or `silent`. See also [print_message/2](printmsg.html#print_message/2) and [current_prolog_flag/2](flags.html#current_prolog_flag/2).

**check_script**(`Bool`)  
If `false` (default `true`), do not check the first character to be `#` and skip the first line when found.

**derived_from**(`File`)  
Indicate that the loaded file is derived from `File`. Used by [make/0](consulting.html#make/0) to time-check and load the original file rather than the derived file.

**dialect**(`+Dialect`)  
Load `Files` with enhanced compatibility with the target Prolog system identified by `Dialect`. See [expects_dialect/1](dialect.html#expects_dialect/1) and [section C](dialect.html#sec:C) for details.

**encoding**(`Encoding`)  
Specify the way characters are encoded in the file. Default is taken from the Prolog flag [encoding](flags.html#flag:encoding). See [section 2.18.1](widechars.html#sec:2.18.1) for details.

**expand**(`Bool`)  
If `true`, run the filenames through [expand_file_name/2](files.html#expand_file_name/2) and load the returned files. Default is `false`, except for [consult/1](consulting.html#consult/1) which is intended for interactive use. Flexible location of files is defined by [file_search_path/2](consulting.html#file_search_path/2).

**format**(`+Format`)  
Used to specify the file format if data is loaded from a stream using the `stream(Stream)` option. Default is `source`, loading Prolog source text. If `qlf`, load QLF data (see [qcompile/1](consulting.html#qcompile/1)).

**if**(`Condition`)  
Load the file only if the specified condition is satisfied. The value `true` loads the file unconditionally, `changed` loads the file if it was not loaded before or has been modified since it was loaded the last time, `not_loaded` loads the file if it was not loaded before, and `exists` is as `changed`, but the call [load_files/2](consulting.html#load_files/2) silently if the file does not exist.

**imports**(`Import`)  
Specify what to import from the loaded module. The default for [use_module/1](import.html#use_module/1) is `all`. `Import` is passed from the second argument of [use_module/2](import.html#use_module/2). Traditionally it is a list of predicate indicators to import. As part of the SWI-Prolog/YAP integration, we also support “`Pred` `as` `Name`” to import a predicate under another name. Finally, `Import` can be the term `except(Exceptions)`, where `Exceptions` is a list of predicate indicators that specify predicates that are *not* imported or `Pred` as `Name` terms to denote renamed predicates. See also [reexport/2](reexport.html#reexport/2) and [use_module/2](import.html#use_module/2).^(bug`Name`/`Arity` as `NewName` is currently implemented using a *link clause*. This harms efficiency and does not allow for querying the relation through [predicate_property/2](examineprog.html#predicate_property/2).)

If `Import` equals `all`, all operators are imported as well. Otherwise, operators are *not* imported. Operators can be imported selectively by adding terms `op(Pri,Assoc,Name)` to the `Import` list. If such a term is encountered, all exported operators that unify with this term are imported. Typically, this construct will be used with all arguments unbound to import all operators or with only `Name` bound to import a particular operator.

**modified**(`TimeStamp`)  
Claim that the source was loaded at `TimeStamp` without checking the source. This option is intended to be used together with the `stream(Input)` option, for example after extracting the time from an HTTP server or database.

**module**(`+Module`)  
Load the indicated file into the given module, overruling the module name specified in the `:- module(Name, ...)` directive. This currently serves two purposes: (1) allow loading two module files that specify the same module into the same process and force and (2): force loading source code in a specific module, even if the code provides its own module name. Experimental.

**must_be_module**(`Bool`)  
If `true`, raise an error if the file is not a module file. Used by [use_module/\[1,2\]](import.html#use_module/1).

**qcompile**(`Atom`)  
How to deal with quick-load-file compilation by [qcompile/1](consulting.html#qcompile/1). Values are:

**never**  
Default. Do not use qcompile unless called explicitly.

**auto**  
Use qcompile for all writeable files. See comment below.

**large**  
Use qcompile if the file is‘large’. Currently, files larger than 100 Kbytes are considered large.

**part**  
If [load_files/2](consulting.html#load_files/2) appears in a directive of a file that is compiled into Quick Load Format using [qcompile/1](consulting.html#qcompile/1), the contents of the argument files are included in the `.qlf` file instead of the loading directive.

If this option is not present, it uses the value of the Prolog flag [qcompile](flags.html#flag:qcompile) as default.

**optimise**(`+Boolean`)  
Explicitly set the optimization for compiling this module. See [optimise](flags.html#flag:optimise).

**redefine_module**(`+Action`)  
Defines what to do if a file is loaded that provides a module that is already loaded from another file. `Action` is one of `false` (default), which prints an error and refuses to load the file, or `true`, which uses [unload_file/1](consulting.html#unload_file/1) on the old file and then proceeds loading the new file. Finally, there is `ask`, which starts interaction with the user. `ask` is only provided if the stream `user_input` is associated with a terminal.

**reexport**(`Bool`)  
If `true` re-export the imported predicate. Used by [reexport/1](reexport.html#reexport/1) and [reexport/2](reexport.html#reexport/2).

**register**(`Bool`)  
If `false`, do not register the load location and options. This option is used by [make/0](consulting.html#make/0) and load_hotfixes/1 to avoid polluting the load-context database. See [source_file_property/2](consulting.html#source_file_property/2).

**sandboxed**(`Bool`)  
Load the file in *sandboxed* mode. This option controls the flag [sandboxed_load](flags.html#flag:sandboxed_load). The only meaningful value for `Bool` is `true`. Using `false` while the Prolog flag is set to `true` raises a permission error.

**scope_settings**(`Bool`)  
Scope [style_check/1](debugger.html#style_check/1) and [expects_dialect/1](dialect.html#expects_dialect/1) to the file and files loaded from the file after the directive. Default is `true`. The system and user initialization files (see **-f** and **-F**) are loading with `scope_settings(false)`.

**silent**(`Bool`)  
If `true`, load the file without printing a message. The specified value is the default for all files loaded as a result of loading the specified files. This option writes the Prolog flag [verbose_load](flags.html#flag:verbose_load) with the negation of `Bool`.

**stream**(`Input`)  
This SWI-Prolog extension compiles the data from the stream `Input`. If this option is used, `Files` must be a single atom which is used to identify the source location of the loaded clauses as well as to remove all clauses if the data is reconsulted.

This option is added to allow compiling from non-file locations such as databases, the web, the *user* (see [consult/1](consulting.html#consult/1)) or other servers. It can be combined with `format(qlf)` to load QLF data from a stream.

The [load_files/2](consulting.html#load_files/2) predicate can be hooked to load other data or data from objects other than files. See [prolog_load_file/2](loadfilehook.html#prolog_load_file/2) for a description and `library(http/http_load)` for an example. All hooks for [load_files/2](consulting.html#load_files/2) are documented in [section B.10](loadfilehook.html#sec:B.10).

**consult**(`:File`)  
Read `File` as a Prolog source file. Calls to [consult/1](consulting.html#consult/1) may be abbreviated by just typing a number of filenames in a list. Examples:

|                        |                                |
|------------------------|--------------------------------|
| `?- consult(load).`    | % consult `load` or `load.pl`  |
| `?- [library(lists)].` | % load library lists           |
| `?- [user].`           | % Type program on the terminal |

The predicate [consult/1](consulting.html#consult/1) is equivalent to `load_files(File, [])`, except for handling the special file `user`, which reads clauses from the terminal. See also the `stream(Input)` option of [load_files/2](consulting.html#load_files/2). Abbreviation using `?- [file1,file2].` does *not* work for the empty list (`[]`). This facility is implemented by defining the list as a predicate. Applications may only rely on using the list abbreviation at the Prolog toplevel and in directives.

**ensure_loaded**(`:File`)  
If the file is not already loaded, this is equivalent to [consult/1](consulting.html#consult/1). Otherwise, if the file defines a module, import all public predicates. Finally, if the file is already loaded, is not a module file, and the context module is not the global user module, [ensure_loaded/1](consulting.html#ensure_loaded/1) will call [consult/1](consulting.html#consult/1).

With this semantics, we hope to get as close as possible to the clear semantics without the presence of a module system. Applications using modules should consider using [use_module/\[1,2\]](import.html#use_module/1).

Equivalent to `load_files(Files, [if(not_loaded)]).`^(56On older versions the condition used to be `if(changed)`. Poor time management on some machines or copying often caused problems. The [make/0](consulting.html#make/0) predicate deals with updating the running system after changing the source code.)

\[ISO\]**include**(`+File`)  
Textually include the content of `File` at the position where the *directive* `:- include(File).` appears. The include construct is only honoured if it appears as a directive in a source file. *Textual* include (similar to C/C++ \#include) is obviously useful for sharing declarations such as [dynamic/1](dynamic.html#dynamic/1) or [multifile/1](dynamic.html#multifile/1) by including a file with directives from multiple files that use these predicates.

Textually including files that contain *clauses* is less obvious. Normally, in SWI-Prolog, clauses are *owned* by the file in which they are defined. This information is used to *replace* the old definition after the file has been modified and is reloaded by, e.g., [make/0](consulting.html#make/0). As we understand it, [include/1](consulting.html#include/1) is intended to include the same file multiple times. Including a file holding clauses multiple times into the same module is rather meaningless as it just duplicates the same clauses. Including a file holding clauses in multiple modules does not suffer from this problem, but leads to multiple equivalent *copies* of predicates. Using [use_module/1](import.html#use_module/1) can achieve the same result while *sharing* the predicates.

If [include/1](consulting.html#include/1) is used to load files holding clauses, and if these files are loaded only once, then these [include/1](consulting.html#include/1) directives can be replaced by other predicates (such as [consult/1](consulting.html#consult/1)). However, there are several cases where either [include/1](consulting.html#include/1) has no alternative, or using any alternative also requires other changes. An example of the former is using [include/1](consulting.html#include/1) to share directives. An example of the latter are cases where clauses of different predicates are distributed over multiple files: If these files are loaded with [include/1](consulting.html#include/1), the directive [discontiguous/1](dynamic.html#discontiguous/1) is appropriate, whereas if they are consulted, one must use the directive [multifile/1](dynamic.html#multifile/1).

To accommodate included files holding clauses, SWI-Prolog distinguishes between the source location of a clause (in this case the included file) and the *owner* of a clause (the file that includes the file holding the clause). The source location is used by, e.g., [edit/1](edit.html#edit/1), the graphical tracer, etc., while the owner is used to determine which clauses are removed if the file is modified. Relevant information is found with the following predicates:

- [source_file/2](consulting.html#source_file/2) describes the owner relation.
- [predicate_property/2](examineprog.html#predicate_property/2) describes the source location (of the first clause).
- [clause_property/2](examineprog.html#clause_property/2) provides access to both source and ownership.
- [source_file_property/2](consulting.html#source_file_property/2) can be used to query include relationships between files.

**require**(`+Predicates`)  
Declare that this file/module requires the specified predicates to be defined “with their commonly accepted definition” . `Predicates` is either a list of predicate indicators or a *comma-list* of predicate indicators. First, all built-in predicates are removed from the set. The remaining predicates are searched using the library index used for autoloading and mapped to a set of [autoload/2](module-autoload.html#autoload/2) directives. This implies that the targets will be loaded lazily if autoloading is not completely disabled and loaded using [use_module/2](import.html#use_module/2) otherwise. See [autoload](flags.html#flag:autoload).

The [require/1](consulting.html#require/1) directive provides less control over the exact nature and location of the predicate. As [autoload/2](module-autoload.html#autoload/2), it prevents a local definition of this predicate. As SWI-Prolog guarantees that the set of built-in predicates and predicates available for autoloading is unambiguous (i.e., has no duplicates) the specification is unambiguous. It provides four advantages over [autoload/2](module-autoload.html#autoload/2): (1) the user does not have to remember the exact library, (2) the directive can be supported in other Prolog systems^(57SICStus provides it), providing compatibility despite differences in library and built-in predicate organization, (3) it is robust against changes to the SWI-Prolog libraries and (4) it is less typing.

**encoding**(`+Encoding`)  
This directive can appear anywhere in a source file to define how characters are encoded in the remainder of the file. It can be used in files that are encoded with a superset of US-ASCII, currently UTF-8 and ISO Latin-1. See also [section 2.18.1](widechars.html#sec:2.18.1).

**make**  
Consult all source files that have been changed since they were consulted. It checks `all` loaded source files: files loaded into a compiled state using `pl -c ...` and files loaded using [consult/1](consulting.html#consult/1) or one of its derivatives. The predicate [make/0](consulting.html#make/0) is called after [edit/1](edit.html#edit/1), automatically reloading all modified files. If the user uses an external editor (in a separate window), [make/0](consulting.html#make/0) is normally used to update the program after editing. In addition, [make/0](consulting.html#make/0) updates the autoload indices (see [section 2.14](autoload.html#sec:2.14)) and runs [list_undefined/0](check.html#list_undefined/0) from the `library(check)` library to report on undefined predicates.

**library_directory**(`?Atom`)  
Dynamic predicate used to specify library directories. Defaults to `app_config(lib)` (see [file_search_path/2](consulting.html#file_search_path/2)) and the system's library (in this order) are defined. The user may add library directories using [assertz/1](db.html#assertz/1), [asserta/1](db.html#asserta/1) or remove system defaults using [retract/1](db.html#retract/1). Deprecated. New code should use [file_search_path/2](consulting.html#file_search_path/2).

**file_search_path**(`+Alias, -Path`)  
Dynamic multifile hook predicate used to specify‘path aliases’. This hook is called by [absolute_file_name/3](files.html#absolute_file_name/3) to search files specified as `Alias(Name)`, e.g., `library(lists)`. This feature is best described using an example. Given the definition:

``` code
file_search_path(demo, '/usr/lib/prolog/demo').
```

the file specification `demo(myfile)` will be expanded to `/usr/lib/prolog/demo/myfile`. The second argument of [file_search_path/2](consulting.html#file_search_path/2) may be another alias.

Below is the initial definition of the file search path. This path implies `swi(<``Path``>)` and refers to a file in the SWI-Prolog home directory. The alias `foreign(<``Path``>)` is intended for storing shared libraries (`.so` or `.DLL` files). See also [use_foreign_library/1](foreignlink.html#use_foreign_library/1).

``` code
user:(file_search_path(library, Dir) :-
        library_directory(Dir)).
user:file_search_path(swi, Home) :-
    current_prolog_flag(home, Home).
user:file_search_path(swi, Home) :-
    current_prolog_flag(shared_home, Home).
user:file_search_path(library, app_config(lib)).
user:file_search_path(library, swi(library)).
user:file_search_path(library, swi(library/clp)).
user:file_search_path(foreign, swi(ArchLib)) :-
    current_prolog_flag(apple_universal_binary, true),
    ArchLib = 'lib/fat-darwin'.
user:file_search_path(foreign, swi(ArchLib)) :-
    \+ current_prolog_flag(windows, true),
    current_prolog_flag(arch, Arch),
    atom_concat('lib/', Arch, ArchLib).
user:file_search_path(foreign, swi(ArchLib)) :-
    current_prolog_flag(msys2, true),
    current_prolog_flag(arch, Arch),
    atomic_list_concat([lib, Arch], /, ArchLib).
user:file_search_path(foreign, swi(SoLib)) :-
    current_prolog_flag(msys2, true),
    current_prolog_flag(arch, Arch),
    atomic_list_concat([bin, Arch], /, SoLib).
user:file_search_path(foreign, swi(SoLib)) :-
    (   current_prolog_flag(windows, true)
    ->  SoLib = bin
    ;   SoLib = lib
    ).
user:file_search_path(path, Dir) :-
    getenv('PATH', Path),
    (   current_prolog_flag(windows, true)
    ->  atomic_list_concat(Dirs, (;), Path)
    ;   atomic_list_concat(Dirs, :, Path)
    ),
    '$member'(Dir, Dirs).
user:file_search_path(user_app_data, Dir) :-
    '$xdg_prolog_directory'(data, Dir).
user:file_search_path(common_app_data, Dir) :-
    '$xdg_prolog_directory'(common_data, Dir).
user:file_search_path(user_app_config, Dir) :-
    '$xdg_prolog_directory'(config, Dir).
user:file_search_path(common_app_config, Dir) :-
    '$xdg_prolog_directory'(common_config, Dir).
user:file_search_path(app_data, user_app_data('.')).
user:file_search_path(app_data, common_app_data('.')).
user:file_search_path(app_config, user_app_config('.')).
user:file_search_path(app_config, common_app_config('.')).
user:file_search_path(app, swi(app)).
user:file_search_path(app, app_data(app)).
user:file_search_path(working_directory, CWD) :-
    working_directory(CWD, CWD).
```

The '\$xdg_prolog_directory'/2 uses either the [XDG Base Directory](https://wiki.archlinux.org/index.php/XDG_Base_Directory) or [win_folder/2](system.html#win_folder/2) on Windows. On Windows, user config is mapped to roaming appdata (CSIDL_APPDATA), user data to the non-roaming (CSIDL_LOCAL_APPDATA) and common data to (CSIDL_COMMON_APPDATA).

The [file_search_path/2](consulting.html#file_search_path/2) expansion is used by all loading predicates as well as by [absolute_file_name/\[2,3\]](files.html#absolute_file_name/2).

The Prolog flag [verbose_file_search](flags.html#flag:verbose_file_search) can be set to `true` to help debugging Prolog's search for files.

\[nondet\]**expand_file_search_path**(`+Spec, -Path`)  
Unifies `Path` with all possible expansions of the filename specification `Spec`. See also [absolute_file_name/3](files.html#absolute_file_name/3).

**prolog_file_type**(`?Extension, ?Type`)  
This dynamic multifile predicate defined in module `user` determines the extensions considered by [file_search_path/2](consulting.html#file_search_path/2). `Extension` is the filename extension without the leading dot, and `Type` denotes the type as used by the `file_type(Type)` option of [file_search_path/2](consulting.html#file_search_path/2). Here is the initial definition of [prolog_file_type/2](consulting.html#prolog_file_type/2):

``` code
user:prolog_file_type(pl,       prolog).
user:prolog_file_type(Ext,      prolog) :-
        current_prolog_flag(associate, Ext),
        Ext \== pl.
user:prolog_file_type(qlf,      qlf).
user:prolog_file_type(Ext,      executable) :-
        current_prolog_flag(shared_object_extension, Ext).
```

Users can add extensions for Prolog source files to avoid conflicts (for example with **perl**) as well as to be compatible with another Prolog implementation. We suggest using `.pro` for avoiding conflicts with **perl**. Overriding the system definitions can stop the system from finding libraries.

**source_file**(`?File`)  
True if `File` is a loaded Prolog source file. `File` is the absolute and canonical path to the source file.

**source_file**(`:Pred, ?File`)  
True if the predicate specified by `Pred` is owned by file `File`, where `File` is an absolute path name (see [absolute_file_name/2](files.html#absolute_file_name/2)). Can be used with any instantiation pattern, but the database only maintains the source file for each predicate. If `Pred` is a *multifile* predicate this predicate succeeds for all files that contribute clauses to `Pred`.^(58The current implementation performs a linear scan through all clauses to establish this set of files.) See also [clause_property/2](examineprog.html#clause_property/2). Note that the relation between files and predicates is more complicated if [include/1](consulting.html#include/1) is used. The predicate describes the *owner* of the predicate. See [include/1](consulting.html#include/1) for details.

**source_file_property**(`?File, ?Property`)  
True when `Property` is a property of the loaded file `File`. If `File` is non-var, it can be a file specification that is valid for [load_files/2](consulting.html#load_files/2). Defined properties are:

**derived_from**(`Original, OriginalModified`)  
`File` was generated from the file `Original`, which was last modified at time `OriginalModified` at the time it was loaded. This property is available if `File` was loaded using the `derived_from(Original)` option to [load_files/2](consulting.html#load_files/2).

**includes**(`IncludedFile, IncludedFileModified`)  
`File` used [include/1](consulting.html#include/1) to include `IncludedFile`. The last modified time of `IncludedFile` was `IncludedFileModified` at the time it was included.

**included_in**(`MasterFile, Line`)  
`File` was included into `MasterFile` from line `Line`. This is the inverse of the `includes` property.

**load_context**(`Module, Location, Options`)  
`Module` is the module into which the file was loaded. If `File` is a module, this is the module into which the exports are imported. Otherwise it is the module into which the clauses of the non-module file are loaded. `Location` describes the file location from which the file was loaded. It is either a term \<`file`\>:\<`line`\> or the atom `user` if the file was loaded from the terminal or another unknown source. `Options` are the options passed to [load_files/2](consulting.html#load_files/2). Note that all predicates to load files are mapped to [load_files/2](consulting.html#load_files/2), using the option argument to specify the exact behaviour.

**load_count**(`-Count`)  
`Count` is the number of times the file have been loaded, i.e., 1 (one) if the file has been loaded once.

**modified**(`Stamp`)  
File modification time when `File` was loaded. This is used by [make/0](consulting.html#make/0) to find files whose modification time is different from when it was loaded.

**source**(`Source`)  
One of `file` if the source was loaded from a file, `resource` if the source was loaded from a resource or `state` if the file was included in the saved state.

**module**(`Module`)  
`File` is a module file that declares the module `Module`.

**number_of_clauses**(`Count`)  
`Count` is the number of clauses associated with `File`. Note that clauses loaded from included files are counted as part of the main file.

**reloading**  
Present if the file is currently being **re**loaded.

\[semidet\]**exists_source**(`+Source`)  
True if `Source` (a term valid for [load_files/2](consulting.html#load_files/2)) exists. Fails without error if this is not the case. The predicate is intended to be used with *conditional compilation* (see [section 4.3.1.2](consulting.html#sec:4.3.1.2) For example:

``` code
:- if(exists_source(library(error))).
:- use_module_library(error).
:- endif.
```

The implementation uses [absolute_file_name/3](files.html#absolute_file_name/3) using `file_type(prolog)`.

\[semidet\]**exists_source**(`+Source, -File`)  
As [exists_source/1](consulting.html#exists_source/1), binding `File` to an atom describing the full absolute path to the source file.

**unload_file**(`+File`)  
Remove all clauses loaded from `File`. If `File` loaded a module, clear the module's export list and disassociate it from the file. `File` is a canonical filename or a file indicator that is valid for [load_files/2](consulting.html#load_files/2).

This predicate should be used with care. The multithreaded nature of SWI-Prolog makes removing static code unsafe. Attempts to do this should be reserved for development or situations where the application can guarantee that none of the clauses associated to `File` are active.

**prolog_load_context**(`?Key, ?Value`)  
Obtain context information during compilation. This predicate can be used from directives appearing in a source file to get information about the file being loaded as well as by the [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2) hooks. See also [source_location/2](consulting.html#source_location/2) and [if/1](consulting.html#if/1). The following keys are defined:

|  |  |
|----|----|
| **Key** | **Description** |
| `directory` | Directory in which `source` lives (absolute path) |
| `dialect` | Compatibility mode. See [expects_dialect/1](dialect.html#expects_dialect/1). |
| `file` | Similar to `source`, but returns the file being included when called while an include file is being processed (absolute path) |
| `module` | Module into which file is loaded |
| `reload` | `true` if the file is being **re**loaded. Not present on first load |
| `script` | Boolean that indicates whether the file is loaded as a script file (see **-s**) |
| `source` | File being loaded (absolute path). If the system is processing an included file, the value is the *main* file. Returns the original Prolog file when loading a `.qlf` file. |
| `stream` | Stream identifier (see [current_input/1](IO.html#current_input/1)) |
| `term_position` | Start position of last term read. See also [stream_property/2](IO.html#stream_property/2) (`position` property and [stream_position_data/3](IO.html#stream_position_data/3).^(59Up to version 7.1.22, the position term carried fake data except for the `line_count` and had **five** arguments, where the position property of a stream only has **four**.) |
| `term` | Term being expanded by [expand_term/2](consulting.html#expand_term/2). |
| `variable_names` | A list of‘`Name` = `Var`’of the last term read. See [read_term/2](termrw.html#read_term/2) for details. |

The `directory` is commonly used to add rules to [file_search_path/2](consulting.html#file_search_path/2), setting up a search path for finding files with [absolute_file_name/3](files.html#absolute_file_name/3). For example:

``` code
:- dynamic user:file_search_path/2.
:- multifile user:file_search_path/2.

:- prolog_load_context(directory, Dir),
   asserta(user:file_search_path(my_program_home, Dir)).

    ...
    absolute_file_name(my_program_home('README.TXT'), ReadMe,
                       [ access(read) ]),
    ...
```

**source_location**(`-File, -Line`)  
If the last term has been read from a physical file (i.e., not from the file `user` or a string), unify `File` with an absolute path to the file and `Line` with the line number in the file. New code should use [prolog_load_context/2](consulting.html#prolog_load_context/2).

**at_halt**(`:Goal`)  
Register `Goal` to be run from [PL_cleanup()](foreigninclude.html#PL_cleanup()), which is called when the system halts. The hooks are run in the reverse order they were registered (FIFO). Success or failure executing a hook is ignored. If the hook raises an exception this is printed using [print_message/2](printmsg.html#print_message/2). An attempt to call [halt/\[0,1\]](toplevel.html#halt/0) from a hook is ignored. Hooks may call [cancel_halt/1](consulting.html#cancel_halt/1), causing [halt/0](toplevel.html#halt/0) and [PL_halt(0)](foreigninclude.html#PL_halt()) to print a message indicating that halting the system has been cancelled.

**cancel_halt**(`+Reason`)  
If this predicate is called from a hook registered with [at_halt/1](consulting.html#at_halt/1), halting Prolog is cancelled and an informational message is printed that includes `Reason`. This is used by the development tools to cancel halting the system if the editor has unsaved data and the user decides to cancel.

\[ISO\]:- **initialization**(`:Goal`)  
Call `Goal` *after* loading the source file in which this directive appears has been completed. In addition, `Goal` is executed if a saved state created using [qsave_program/1](saved-states.html#qsave_program/1) is restored.

The ISO standard only allows for using `:- Term` if `Term` is a *directive*. This means that arbitrary goals can only be called from a directive by means of the [initialization/1](consulting.html#initialization/1) directive. SWI-Prolog does not enforce this rule.

The [initialization/1](consulting.html#initialization/1) directive must be used to do program initialization in saved states (see [qsave_program/1](saved-states.html#qsave_program/1)). A saved state contains the predicates, Prolog flags and operators present at the moment the state was created. Other resources (records, foreign resources, etc.) must be recreated using [initialization/1](consulting.html#initialization/1) directives or from the entry goal of the saved state.

Up to SWI-Prolog 5.7.11, `Goal` was executed immediately rather than after loading the program text in which the directive appears as dictated by the ISO standard. In many cases the exact moment of execution is irrelevant, but there are exceptions. For example, [load_foreign_library/1](foreignlink.html#load_foreign_library/1) must be executed immediately to make the loaded foreign predicates available for exporting. SWI-Prolog now provides the directive [use_foreign_library/1](foreignlink.html#use_foreign_library/1) to ensure immediate loading as well as loading after restoring a saved state. If the system encounters a directive `:- initialization(load_foreign_library(...))`, it will load the foreign library immediately and issue a warning to update your code. This behaviour can be extended by providing clauses for the multifile hook predicate `prolog:initialize_now(Term, Advice)`, where `Advice` is an atom that gives advice on how to resolve the compatibility issue.

**initialization**(`:Goal, +When`)  
Similar to [initialization/1](consulting.html#initialization/1), but allows for specifying when `Goal` is executed while loading the program text:

**now**  
Execute `Goal` immediately.

**after_load**  
Execute `Goal` after loading the program text in which the directive appears. This is the same as [initialization/1](consulting.html#initialization/1).

**prepare_state**  
Execute `Goal` as part of [qsave_program/2](saved-states.html#qsave_program/2). This hook can be used for example to eagerly execute initialization that is normally done lazily on first usage.

**restore_state**  
Do not execute `Goal` while loading the program, but *only* when restoring a saved state.^(60Used to be called `restore`. `restore` is still accepted for backward compatibility.)

**program**  
Execute `Goal` once after executing the **-g** goals at program startup. Registered goals are executed in the order encountered and a failure or exception causes the Prolog to exit with non-zero exit status. These goals are *not* executed if the **-l** is given to merely *load* files. In that case they may be executed explicitly using [initialize/0](consulting.html#initialize/0). See also [section 2.11.1.1](compilation.html#sec:2.11.1.1).

**main**  
When Prolog starts, the last goal registered using `initialization(Goal, main)` is executed as main goal. If `Goal` fails or raises an exception, the process terminates with non-zero exit code. If not explicitly specified using the **-t** the *toplevel goal* is set to [halt/0](toplevel.html#halt/0), causing the process to exit with status 0. An explicitly specified toplevel is executed normally. This implies that `-t prolog` causes the application to start the normal interactive toplevel after completing `Goal`. See also the Prolog flag [toplevel_goal](flags.html#flag:toplevel_goal) and [section 2.11.1.1](compilation.html#sec:2.11.1.1).

\[det\]**initialize**  
Run all initialization goals registered using `initialization(Goal, program)`. Raises an error `initialization_error(Reason, Goal, File:Line)` if `Goal` fails or raises an exception. `Reason` is `failed` or the exception raised.

**compiling**  
True if the system is compiling source files with the **-c** option or [qcompile/1](consulting.html#qcompile/1) into an intermediate code file. Can be used to perform conditional code optimisations in [term_expansion/2](consulting.html#term_expansion/2) (see also the **-O** option) or to omit execution of directives during compilation.

### 4.3.1 Conditional compilation and program transformation

ISO Prolog defines no way for program transformations such as macro expansion or conditional compilation. Expansion through [term_expansion/2](consulting.html#term_expansion/2) and [expand_term/2](consulting.html#expand_term/2) can be seen as part of the de-facto standard. This mechanism can do arbitrary translation between valid Prolog terms read from the source file to Prolog terms handed to the compiler. As [term_expansion/2](consulting.html#term_expansion/2) can return a list, the transformation does not need to be term-to-term.

Various Prolog dialects provide the analogous [goal_expansion/2](consulting.html#goal_expansion/2) and [expand_goal/2](consulting.html#expand_goal/2) that allow for translation of individual body terms, freeing the user of the task to disassemble each clause.

**term_expansion**(`+Term1, -Term2`)  
Dynamic and multifile predicate, normally not defined. When defined by the user all terms read during consulting are given to this predicate. If the predicate succeeds Prolog will assert `Term2` in the database rather than the read term (`Term1`). `Term2` may be a term of the form `?- Goal.` or `:- Goal`. `Goal` is then treated as a directive. If `Term2` is a list, all terms of the list are stored in the database or called (for directives). If `Term2` is of the form below, the system will assert `Clause` and record the indicated source location with it:

> `’$source_location’(<``File``>, <``Line``>):<``Clause``>`

When compiling a module (see [chapter 6](modules.html#sec:6) and the directive [module/2](defmodule.html#module/2)), [expand_term/2](consulting.html#expand_term/2) will first try [term_expansion/2](consulting.html#term_expansion/2) in the module being compiled to allow for term expansion rules that are local to a module. If there is no local definition, or the local definition fails to translate the term, [expand_term/2](consulting.html#expand_term/2) will try [term_expansion/2](consulting.html#term_expansion/2) in module `user`. For compatibility with SICStus and Quintus Prolog, this feature should not be used. See also [expand_term/2](consulting.html#expand_term/2), [goal_expansion/2](consulting.html#goal_expansion/2) and [expand_goal/2](consulting.html#expand_goal/2).

It is possible to act on the beginning and end of a file by expanding the terms `begin_of_file` and `end_of_file`. The latter is supported by most Prolog systems that support term expansion as [read_term/3](termrw.html#read_term/3) returns `end_of_file` on reaching the end of the input. Expanding `begin_of_file` may be used to initialise the compilation, for example base on the file name extension. It was added in SWI-Prolog 8.1.1.

The current macro-expansion mechanism originates from Prolog systems in the 1980s and 1990s. It has several flaws, (1) the hooks act globally (except for definitions in a module), (2) it is hard to deal with interactions between transformations, (3) macros can not be reused between modules using the normal module export/import protocol and (4) it is hard to make source code aware tools such as the graphical debugger act properly in the context of macro expansion. Several Prolog implementations have tried to implement better expansion mechanisms. None of these solve all problems and all are largely incompatible with our current macro expansion. Future versions may provide a new mechanism to solve these issues.

Controlled interaction is provided between macro expansion defined in a module and the `user` and `system` modules. Here, SWI-Prolog uses a *pipeline* where the result of local module expansion is the input for the expansion in `user`, which is the input for the expansion in `system`. See also [section 6.10](importmodule.html#sec:6.10).

*Scoping*, i.e., make a rule defined in a module only active if this module is imported into the module being compiled, can be emulated by defining the macro globally in the `user` module and using [prolog_load_context/2](consulting.html#prolog_load_context/2) and some logic to verify the macro expansion should apply. If (goal) expansion effectively defined *inlining* it is good practice to also define the predicate and have the macro expansion check that the predicate is in scope. Here is an example.

``` code
:- module(m1, [double/2]).

double(X, D) :- D is X*2.

user:goal_expansion(double(X,D), D is X*2) :-
    prolog_load_context(module, M),
    predicate_property(M:double(_,_), imported_from(m1)).
```

For term expansion that is not related to a specific predicate we can define a sentinel predicate rather than using the goal predicate and check it is imported into the current module to verify that the module that defines the expansion is imported into the current compilation context.

**expand_term**(`+Term1, -Term2`)  
This predicate is normally called by the compiler on terms read from the input to perform preprocessing. It consists of four steps, where each step processes the output of the previous step.

1.  Test conditional compilation directives and translate all input to `[]` if we are in a‘false branch’of the conditional compilation. See [section 4.3.1.2](consulting.html#sec:4.3.1.2).
2.  Call [term_expansion/2](consulting.html#term_expansion/2). This predicate is first tried in the module that is being compiled and then in modules from which this module inherits according to [default_module/2](importmodule.html#default_module/2). The output of the expansion in a module is used as input for the next module. Using the default setup and when compiling a normal application module `M`, this implies expansion is executed in `M`, `user` and finally in `system`. Library modules inherit directly from `system` and can thus not be re-interpreted by term expansion rules in `user`.
3.  Call DCG expansion ([dcg_translate_rule/2](consulting.html#dcg_translate_rule/2)).
4.  Call [expand_goal/2](consulting.html#expand_goal/2) on each body term that appears in the output of the previous steps.

**goal_expansion**(`+Goal1, -Goal2`)  
Like [term_expansion/2](consulting.html#term_expansion/2), [goal_expansion/2](consulting.html#goal_expansion/2) provides for macro expansion of Prolog source code. Between [expand_term/2](consulting.html#expand_term/2) and the actual compilation, the body of clauses analysed and the goals are handed to [expand_goal/2](consulting.html#expand_goal/2), which uses the [goal_expansion/2](consulting.html#goal_expansion/2) hook to do user-defined expansion.

The predicate [goal_expansion/2](consulting.html#goal_expansion/2) is first called in the module that is being compiled, and then follows the module inheritance path as defined by [default_module/2](importmodule.html#default_module/2), i.e., by default `user` and `system`. If `Goal` is of the form `Module`:`Goal` where `Module` is instantiated, [goal_expansion/2](consulting.html#goal_expansion/2) is called on `Goal` using rules from module `Module` followed by default modules for `Module`.

Only goals appearing in the body of clauses when reading a source file are expanded using this mechanism, and only if they appear literally in the clause, or as an argument to a defined meta-predicate that is annotated using‘0’(see [meta_predicate/1](metapred.html#meta_predicate/1)). Other cases need a real predicate definition.

The expansion hook can use [prolog_load_context/2](consulting.html#prolog_load_context/2) to obtain information about the context in which the goal is expanded such as the module, variable names or the encapsulating term.

**expand_goal**(`+Goal1, -Goal2`)  
This predicate is normally called by the compiler to perform preprocessing using [goal_expansion/2](consulting.html#goal_expansion/2). The predicate computes a fixed-point by applying transformations until there are no more changes. If optimisation is enabled (see **-O** and [optimise](flags.html#flag:optimise)), [expand_goal/2](consulting.html#expand_goal/2) simplifies the result by removing unneeded calls to [true/0](control.html#true/0) and [fail/0](control.html#fail/0) as well as trivially unreachable branches.

If [goal_expansion/2](consulting.html#goal_expansion/2) *wraps* a goal as in the example below the system still reaches fixed-point as it prevents re-expanding the expanded term while recursing. It does re-enable expansion on the *arguments* of the expanded goal as illustrated in t2/1 in the example.^(61After discussion with Peter Ludemann and Paulo Moura on the forum.)

``` code
:- meta_predicate run(0).

may_not_fail(test(_)).
may_not_fail(run(_)).

goal_expansion(G, (G *-> true ; error(goal_failed(G),_))) :-
    may_not_fail(G).

t1(X) :- test(X).
t2(X) :- run(run(X)).
```

Is expanded into

``` code
t1(X) :-
    (   test(X)
    *-> true
    ;   error(goal_failed(test(X)), _)
    ).

t2(X) :-
    (   run((run(X)*->true;error(goal_failed(run(X)), _)))
    *-> true
    ;   error(goal_failed(run(run(X))), _)
    ).
```

Note that goal expansion should not bind any variables in the clause. Doing so may impact the semantics of the clause if the variable is also used elsewhere. In the general case this is not verified. It is verified for [\\/1](control.html#\+/1) and [;/2](control.html#;/2), resulting in an exception.

**compile_aux_clauses**(`+Clauses`)  
Compile clauses on behalf of [goal_expansion/2](consulting.html#goal_expansion/2). This predicate compiles the argument clauses into static predicates, associating the predicates with the current file but avoids changing the notion of current predicate and therefore discontiguous warnings.

Note that in some cases multiple expansions of similar goals can share the same compiled auxiliary predicate. In such cases, the implementation of [goal_expansion/2](consulting.html#goal_expansion/2) can use [predicate_property/2](examineprog.html#predicate_property/2) using the property `defined` to test whether the predicate is already defined in the current context.

**dcg_translate_rule**(`+In, -Out`)  
This predicate performs the translation of a term `Head-->Body` into a normal Prolog clause. Normally this functionality should be accessed using [expand_term/2](consulting.html#expand_term/2).

**var_property**(`+Var, ?Property`)  
True when `Property` is a property of `Var`. These properties are available during goal- and term-expansion. Defined properties are below. Future versions are likely to provide more properties, such as whether the variable is referenced in the remainder of the term. See also [goal_expansion/2](consulting.html#goal_expansion/2).

**fresh**(`Bool`)  
`Bool` has the value `true` if the variable is guaranteed to be unbound at entry of the goal, otherwise its value is `false`. This implies that the variable first appears in this goal or a previous appearance was in a negation ([\\/1](control.html#\+/1)) or a different branch of a disjunction.

**singleton**(`Bool`)  
`Bool` has the value `true` if the variable is a *syntactic* singleton in the term it appears in. Note that this tests that the variable appears exactly once in the term being expanded without making any claim on the syntax of the variable. Variables that appear only once in multiple branches are *not* singletons according to this property. Future implementations may improve on that.

**name**(`Name`)  
True when variable appears with the given name in the source.

#### 4.3.1.1 Program transformation with source layout info

This sections documents extended versions of the program transformation predicates that also transform the source layout information. Extended layout information is currently processed, but unused. Future versions will use for the following enhancements:

- More precise locations of warnings and errors
- More reliable setting of breakpoints
- More reliable source layout information in the graphical debugger.

**expand_goal**(`+Goal1, ?Layout1, -Goal2, -Layout2`)  
**goal_expansion**(`+Goal1, ?Layout1, -Goal2, -Layout2`)  
**expand_term**(`+Term1, ?Layout1, -Term2, -Layout2`)  
**term_expansion**(`+Term1, ?Layout1, -Term2, -Layout2`)  
**dcg_translate_rule**(`+In, ?LayoutIn, -Out, -LayoutOut`)  
These versions are called *before* their 2-argument counterparts. The input layout term is either a variable (if no layout information is available) or a term carrying detailed layout information as returned by the `subterm_positions` of [read_term/2](termrw.html#read_term/2). The output layout should be a variable if no layout information can be computed for the expansion; a sub-term can also be a variable to indicate “don't know” .

#### 4.3.1.2 Conditional compilation

Conditional compilation builds on the same principle as [term_expansion/2](consulting.html#term_expansion/2), [goal_expansion/2](consulting.html#goal_expansion/2) and the expansion of grammar rules to compile sections of the source code conditionally. One of the reasons for introducing conditional compilation is to simplify writing portable code. See [section C](dialect.html#sec:C) for more information. Here is a simple example:

``` code
:- if(\+source_exports(library(lists), suffix/2)).

suffix(Suffix, List) :-
        append(_, Suffix, List).

:- endif.
```

Note that these directives can only appear as separate terms in the input. SWI-Prolog accommodates syntax extensions under conditional compilation by silently ignoring syntax errors when in the *false* branch. This allow, for example, for the code below. With rational number support `1r3` denotes the rational number 1/3 while without it is a syntax error. Note that this only works properly if (1) the syntax error still allows to re-synchronize on the full stop of the invalid clause and (2) the subsequent conditional compilation directive is valid.

``` code
:- if(current_prolog_flag(bounded, false)).
one_third(1r3).
:- endif.
```

Typical usage scenarios include:

- Load different libraries on different dialects.
- Define a predicate if it is missing as a system predicate.
- Realise totally different implementations for a particular part of the code due to different capabilities.
- Realise different configuration options for your software.

:- **if**(`:Goal`)  
Compile subsequent code only if `Goal` succeeds. For enhanced portability, `Goal` is processed by [expand_goal/2](consulting.html#expand_goal/2) before execution. If an error occurs, the error is printed and processing proceeds as if `Goal` has failed.

:- **elif**(`:Goal`)  
Equivalent to `:- else. :-if(Goal).` ... `:- endif.` In a sequence as below, the section below the first matching `elif` is processed. If no test succeeds, the else branch is processed.

``` code
:- if(test1).
section_1.
:- elif(test2).
section_2.
:- elif(test3).
section_3.
:- else.
section_else.
:- endif.
```

:- **else**  
Start‘else’branch.

:- **endif**  
End of conditional compilation.

### 4.3.2 Reloading files, active code and threads

Traditionally, Prolog environments allow for reloading files holding currently active code. In particular, the following sequence is a valid use of the development environment:

- Trace a goal
- Find unexpected behaviour of a predicate
- Enter a *break* using the **b** command
- Fix the sources and reload them using [make/0](consulting.html#make/0)
- Exit the break, *retry* executing the now fixed predicate using the **r** command

*Reloading* a previously loaded file is safe, both in the debug scenario above and when the code is being executed by another *thread*. Executing threads switch atomically to the new definition of modified predicates, while clauses that belong to the old definition are (eventually) reclaimed by [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0).^(62As of version 7.3.12. Older versions wipe all clauses originating from the file before loading the new clauses. This causes threads that executes the code to (typically) die with an *undefined predicate* exception.) Below we describe the steps taken for *reloading* a file to help understanding the limitations of the process.

1.  If a file is being reloaded, a *reload context* is associated to the file administration. This context includes a table keeping track of predicates and a table keeping track of the module(s) associated with this source.
2.  If a new predicate is found, an entry is added to the context predicate table. Three options are considered:
    1.  The predicate is new. It is handled the same as if the file was loaded for the first time.
    2.  The predicate is foreign or thread local. These too are treated as if the file was loaded for the first time.
    3.  Normal predicates. Here we initialise a pointer to the *current clause*.
3.  New clauses for‘normal predicates’are considered as follows:
    1.  If the clause's byte-code is the same as the predicates current clause, discard the clause and advance the current clause pointer.
    2.  If the clause's byte-code is the same as some clause further into the clause list of the predicate, discard the new clause, mark all intermediate clauses for future deletion, and advance the current clause pointer to the first clause after the matched one.
    3.  If the clause's byte-code matches no clause, insert it for *future activation* before the current clause and keep the current clause.
4.  *Properties* such as `dynamic` or `meta_predicate` are in part applied immediately and in part during the fixup process after the file completes loading. Currently, `dynamic` and `thread_local` are applied immediately.
5.  New modules are recorded in the reload context. Export declarations (the module's public list and [export/1](altmoduleapi.html#export/1) calls) are both applied and recorded.
6.  When the end-of-file is reached, the following fixup steps are taken
    1.  For each predicate
        1.  The current clause and subsequent clauses are marked for future deletion.
        2.  All clauses marked for future deletion or creation are (in)activated by changing their‘erased’or‘created’ *generation*. Erased clauses are (eventually) reclaimed by the *clause garbage collector*, see [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0).
        3.  Pending predicate property changes are applied.
    2.  For each module
        1.  Exported predicates that are not encountered in the reload context are removed from the export list.

The above generally ensures that changes to the *content* of source files can typically be activated safely using [make/0](consulting.html#make/0). Global changes such as operator changes, changes of module names, changes to multi-file predicates, etc. sometimes require a restart. In almost all cases, the need for restart is indicated by permission or syntax errors during the reload or existence errors while running the program.

In some cases the content of a source file refers‘to itself’. This is notably the case if local rules for [goal_expansion/2](consulting.html#goal_expansion/2) or [term_expansion/2](consulting.html#term_expansion/2) are defined or goals are executed using *directives*.^(63Note that [initialization/1](consulting.html#initialization/1) directives are executed *after* loading the file. SWI-Prolog allows for directives that are executed *while* loading the file using `:- Goal.` or [initialization/2](consulting.html#initialization/2)). Up to version 7.5.12 it was typically needed to reload the file *twice*, once for updating the code that was used for compiling the remainder of the file and once to effectuate this. As of version 7.5.13, conventional *transaction semantics* apply. This implies that for the thread performing the reload the file's content is first wiped and gradually rebuilt, while other threads see an *atomic* update from the old file content to the new.^(64This feature was implemented by Keri Harris.)

#### 4.3.2.1 Errors and warnings during compilation

Errors and warnings reported while compiling a file are reported using [print_message/2](printmsg.html#print_message/2). Typical errors are syntax errors, errors during macro expansion by [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2), compiler errors such as illegal clauses or an attempt to redefine a system predicate and errors caused by executing *directives*, notably using [initialization/1](consulting.html#initialization/1) and [initialization/2](consulting.html#initialization/2).

Merely reporting error messages and warnings is typically desirable for interactive usage. Non-interactive applications often require to be notified of such issues, typically using the *exit code* of the process. We can distinguish two types of errors and warnings: (1) those resulting from loading an invalid program and (2) messages that result from running the program. A typical example is user code that wishes to try something and in case of an error report this and continue.

``` code
    ...,
    E = error(_,_),
    catch(do_something, E,
          print_message(error, E)),
    ...
```

User code may be (and often is) started from directives, while running user code may involve compilation due to autoloading, loading of data files, etc. As a result, it is unclear whether an error message should merely be printed, should result in a non-zero exit status at the end or should immediately terminate the process.

The default behaviour is defined by the Prolog flags [on_error](flags.html#flag:on_error) and [on_warning](flags.html#flag:on_warning). It can be fine tuned by defining the *hook predicate* [message_hook/3](printmsg.html#message_hook/3). The compiler calls [print_message/2](printmsg.html#print_message/2) using the level `silent` and the message below if errors or warnings where printed during the execution of [load_files/2](consulting.html#load_files/2).

**load_file_errors**(`File, Errors, Warnings`)  
Here, `File` is the raw file specification handed to [load_files/2](consulting.html#load_files/2), i.e., `’myfile.pl’` or `library(lists)`, `Errors` is the number of errors printed while loading and `Warnings` is the number of warnings printed while loading. Note that these counts include messages from (initialization) directives.

This allows the user to fine tune the behaviour on errors and, for example, halt the process on a non-zero error count right after loading the file with errors using the code below.

``` code
:- multifile prolog:message_action/2.

prolog:message_action(load_file_errors(_File, Errors, _Warnings),
                      _Level) :-
    Errors > 0,
    halt(1).
```

#### 4.3.2.2 Compilation of mutually dependent code

Large programs are generally split into multiple files. If file `A` accesses predicates from file `B` which accesses predicates from file `A`, we consider this a mutual or circular dependency. If traditional load predicates (e.g., [consult/1](consulting.html#consult/1)) are used to include file `B` from `A` and `A` from `B`, loading either file results in a loop. This is because [consult/1](consulting.html#consult/1) is mapped to [load_files/2](consulting.html#load_files/2) using the option `if(true)(.)` Such programs are typically loaded using a *load file* that consults all required (non-module) files. If modules are used, the dependencies are made explicit using [use_module/1](import.html#use_module/1) statements. The [use_module/1](import.html#use_module/1) predicate, however, maps to [load_files/2](consulting.html#load_files/2) with the option `if(not_loaded)(.)` A [use_module/1](import.html#use_module/1) on an already loaded file merely makes the public predicates of the used module available.

Summarizing, mutual dependency of source files is fully supported with no precautions when using modules. Modules can use each other in an arbitrary dependency graph. When using [consult/1](consulting.html#consult/1), predicate dependencies between loaded files can still be arbitrary, but the consult relations between files must be a proper tree.

#### 4.3.2.3 Compilation with multiple threads

This section discusses compiling files for the first time. For reloading, see [section 4.3.2](consulting.html#sec:4.3.2).

Multiple threads can compile files concurrently. This requires special precautions only if multiple threads wish to load the same file at the same time. Therefore, [load_files/2](consulting.html#load_files/2) checks whether some other thread is already loading the file. If not, it starts loading the file. If a thread detects that another thread is already loading the file the thread blocks until the other thread finishes loading the file. After waiting, and if the file is a module file, it imports the exported predicates and operators from the module.

Note that this schema does not prevent deadlocks under all situations. Consider two mutually dependent (see [section 4.3.2.2](consulting.html#sec:4.3.2.2)) module files `A` and `B`, where thread 1 starts loading `A` and thread 2 starts loading `B` at the same time. Both threads will deadlock when trying to load the used module.

The current implementation does not detect such cases and the involved threads will freeze. This problem can be avoided if a mutually dependent collection of files is always loaded from the same start file.

### 4.3.3 Quick load files

SWI-Prolog supports compilation of individual or multiple Prolog source files into‘Quick Load Files’. A‘Quick Load File’(`.qlf` file) stores the contents of the file in a precompiled format.

These files load considerably faster than source files and are normally more compact. They are machine-independent and may thus be loaded on any implementation of SWI-Prolog. Note, however, that clauses are stored as virtual machine instructions. Changes to the compiler will generally make old compiled files unusable.

Quick Load Files are created using [qcompile/1](consulting.html#qcompile/1). They are loaded using [consult/1](consulting.html#consult/1) or one of the other file-loading predicates described in [section 4.3](consulting.html#sec:4.3). If [consult/1](consulting.html#consult/1) is given an explicit `.pl` file, it will load the Prolog source. When given a `.qlf` file, it will load the file. When no extension is specified, it will load the `.qlf` file when present and the `.pl` file otherwise.

**qcompile**(`:File`)  
Takes a file specification as [consult/1](consulting.html#consult/1), etc., and, in addition to the normal compilation, creates a *Quick Load File* from `File`. The file extension of this file is `.qlf`. The basename of the Quick Load File is the same as the input file.

If the file contains‘`:- consult(``+File``)`’,‘`:- [``+File``]`’or‘`:- load_files(``+File``, [qcompile(part), ...])`’statements, the referred files are compiled into the same `.qlf` file. Other directives will be stored in the `.qlf` file and executed in the same fashion as when loading the `.pl` file.

For [term_expansion/2](consulting.html#term_expansion/2), the same rules as described in [section 2.11](compilation.html#sec:2.11) apply.

Conditional execution or optimisation may test the predicate [compiling/0](consulting.html#compiling/0).

Source references ([source_file/2](consulting.html#source_file/2)) in the Quick Load File refer to the Prolog source file from which the compiled code originates.

**qcompile**(`:File, +Options`)  
As [qcompile/1](consulting.html#qcompile/1), but processes additional options as defined by [load_files/2](consulting.html#load_files/2). `Options` are passed to [load_files/2](consulting.html#load_files/2). In addition the following options are processed:

**include**(`+Include`)  
What to include into the QLF file. Currently accepts only a single value: the atom `user`. When specified, files loaded indirectly from `File` that to not come from the Prolog library are included into the `.qlf` file. This may be used to generate a single file from an application. The result is comparable to a *save state* (see [qsave_program/2](saved-states.html#qsave_program/2)) with the following differences:

- Only your application code is included. The Prolog libraries and boot files are not.
- Only Prolog code is included, `.qlf` files cannot include arbitrary *resources*.
- The file can be loaded into a running Prolog process, while a saved state can only be loaded into a virgin Prolog virtual machine.
