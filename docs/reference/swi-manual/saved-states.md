
## 14.2 Understanding saved states

A SWI-Prolog *saved state* is a *resource archive* that contains the compiled program in a machine-independent format,^(247Although the compiled code is independent from the CPU and operating system, 32-bit compiled code does not run on the 64-bit emulator, nor the other way around. Conditionally compiled code (see [if/1](consulting.html#if/1)) may also reduce platform independence.) startup options, optionally shared objects/DLLs and optionally additional *resource* files. As of version 7.7.13, the resource archive format is ZIP. A resource file is normally **created** using the commandline option **-c**:

``` code
swipl -o mystate option ... -c file.pl ...
```

The above causes SWI-Prolog to load the given Prolog files and call [qsave_program/2](saved-states.html#qsave_program/2) using options created from the *option ...* in the command above.

A saved state may be **executed** in several ways. The basic mechanism is to use the **-x**:

``` code
swipl -x mystate app-arg ...
```

Saved states may have an arbitrary payload at the *start*. This allows combining a (shell) script or the emulator with the state to turn the state into a single file executable. By default a state starts with a shell script (Unix) or the emulator (Windows).^(248As the default emulator is a short program while the true emulator is in a DLL this keeps the state short.) The options `emulator(File)` and `stand_alone(Bool)` control what is added at the start of the state. Finally, C/C++ programs that embed Prolog may use a static C string that embeds the state into the executable. See [PL_set_resource_db_mem()](foreigninclude.html#PL_set_resource_db_mem()).

### 14.2.1 Creating a saved state

The predicates in this section support creating a saved state. Note that states are commonly created from the commandline using the **-c**, for example:

``` code
swipl -o mystate --foreign=save -c load.pl
```

Long (`--`) options are translated into options for [qsave_program/2](saved-states.html#qsave_program/2). This transformation uses the same conventions as used by [argv_options/3](main.html#argv_options/3), except that the transformation is guided by the option type. This implies that integer and callable options need to have valid syntax and boolean options may be abbreviated to simply `--autoload` or `--no-autoload` as shorthands for `--autoload=true` and `--autoload=false`.

**qsave_program**(`+File, +Options`)  
Saves the current state of the program to the file `File`. The result is a resource archive `File` containing expresses all Prolog data from the running program, all user-defined resources (see [resource/2](program-resources.html#resource/2) and [open_resource/2](program-resources.html#open_resource/2)) and optionally all shared objects/DLLs required by the program for the current architecture. Depending on the `stand_alone` option, the resource is headed by the emulator, a Unix shell script or nothing. `Options` is a list of additional options:

**stack_limit**(`+Bytes`)  
Sets default stack limit for the new process. See the command line option **--stack-limit** and the Prolog flag [stack_limit](flags.html#flag:stack_limit).

**goal**(`:Callable`)  
Initialization goal for the new executable (see **-g**). Two values have special meaning: `prolog` starts the Prolog toplevel and `default` runs [halt/0](toplevel.html#halt/0) if there are initialization goals and the [prolog/0](toplevel.html#prolog/0) toplevel otherwise.

**toplevel**(`:Callable`)  
Top-level goal for the new executable (see **-t**). Similar to [initialization/2](consulting.html#initialization/2) using `main`, the default toplevel is to enter the Prolog interactive shell unless a goal has been specified using `goal(Callable)`.

**init_file**(`+Atom`)  
Default initialization file for the new executable. See **-f**.

**class**(`+Class`)  
If `runtime` (default), read resources from the state and disconnect the code loaded into the state from the original source. If `development`, save the predicates in their current state and keep reading resources from their source (if present). See also [open_resource/3](program-resources.html#open_resource/3).

**autoload**(`+Boolean`)  
If `true` (default), run autoload/0 first. If the class is `runtime` and `autoload` is `true`, the state is supposed to be self contained and autoloading is disabled in the restored state.

**map**(`+File`)  
Dump a human-readable trace of what has been saved in `File`.

**op**(`+Action`)  
One of `save` (default) to save the current operator table or `standard` to use the initial table of the emulator.

**stand_alone**(`+Boolean`)  
If `true`, the emulator is the first part of the state. If the emulator is started it tests whether a saved state is attached to itself and load this state. Provided the application has all libraries loaded, the resulting executable is completely independent from the runtime environment or location where it was built. See also [section 2.11.1.4](compilation.html#sec:2.11.1.4).

**emulator**(`+File`)  
File to use for the emulator or executable used by the startup script. Default is the running Prolog image *after* following symbolic links, e.g., `/usr/lib/swipl/lib/x86_64-linux/swipl`. To create a saved state based on the public executable such that it can run on multiple architectures one can use e.g.

``` code
$ swipl -o myexe --emulator=$(which swipl) -c myload.pl
```

**foreign**(`+Action`)  
If `save`, include shared objects (DLLs) for the current architecture into the saved state. See [current_foreign_library/2](foreignlink.html#current_foreign_library/2), and current_prolog_flag(arch, Arch). If the program **strip** is available, this is first used to reduce the size of the shared object. If a state is started, [use_foreign_library/1](foreignlink.html#use_foreign_library/1) first tries to locate the foreign resource in the resource database. When found it copies the content of the resource to a temporary file and loads it. If possible (Unix), the temporary object is deleted immediately after opening.^(249This option is experimental and currently disabled by default. It will become the default if it proves robust.250Creating a temporary file is the most portable way to load a shared object from a zip file but requires write access to the file system. Future versions may provide shortcuts for specific platforms that bypass the file system.)

If `Action` is of the form `arch(ListOfArches)` then the shared objects for the specified architectures are stored in the saved state. On the command line, the list of architectures can be passed as `--foreign=<``CommaSepArchesList``>`. In order to obtain the shared object file for the specified architectures, [qsave_program/2](saved-states.html#qsave_program/2) calls a user defined hook: `qsave:arch_shlib(+Arch, +FileSpec, -SoPath)`. This hook needs to unify `SoPath` with the absolute path to the shared object for the specified architecture. `FileSpec` is of the form `foreign(Name)`.

At runtime, SWI-Prolog will try to load the shared library which is compatible with the current architecture, obtained by calling `current_prolog_flag(arch, Arch)`. An architecture is compatible if one of the two following conditions is true (tried in order):

1.  There is a shared object in the saved state file which matches the current architecture name (from [current_prolog_flag/2](flags.html#current_prolog_flag/2)) exactly.
2.  The user definable `qsave:compat_arch(Arch1, Arch2)` hook succeeds.

This last one is useful when one wants to produce one shared object file that works for multiple architectures, usually compiling for the lowest common denominator of a certain CPU type. For example, it is common to compile for armv7 if even if the code will be running on newer arm CPUs. It is also useful to provide highly-optimized shared objects for particular architectures.

If `Action` is **copy**, the foreign extensions are copied to the installation location. This feature is currently only supported for Windows, where both DLLs required to run **swipl.exe** and DLLs loaded through extensions are are copied to the same directory as where the executable is saved. Required DLLs are found using win_process_modules/2. Thus, we can create an executable using e.g.,

``` code
swipl -o dist/myprog.exe --foreign=copy -c myprog.pl
```

The `--foreign=copy` option is introduced in 9.3.6.

**undefined**(`+Value`)  
Defines what happens if an undefined predicate is found during the code analysis. Values are `ignore` (default) or `error`. In the latter case creating the state is aborted with a message that indicates the undefines predicates and from where they are called.

**obfuscate**(`+Boolean`)  
If `true` (default `false`), replace predicate names with generated symbols to make the code harder to assess for reverse engineering. See [section 14.6.1](protect-code.html#sec:14.6.1).

**verbose**(`+Boolean`)  
If `true` (default `false`), report progress and status, notably regarding auto loading.

**qsave_program**(`+File`)  
Equivalent to `qsave_program(File, [])`.

**autoload_all**  
Check the current Prolog program for predicates that are referred to, are undefined, and have a definition in the Prolog library. Load the appropriate libraries.

This predicate is used by [qsave_program/\[1,2\]](saved-states.html#qsave_program/1) to ensure the saved state does not depend on availability of the libraries. The predicate [autoload_all/0](saved-states.html#autoload_all/0) examines all clauses of the loaded program (obtained with [clause/2](examineprog.html#clause/2)) and analyzes the body for referenced goals. Such an analysis cannot be complete in Prolog, which allows for the creation of arbitrary terms at runtime and the use of them as a goal. The current analysis is limited to the following:

- Direct goals appearing in the body
- Arguments of declared meta-predicates that are marked with an integer (0..9). See [meta_predicate/1](metapred.html#meta_predicate/1).

The analysis of meta-predicate arguments is limited to cases where the argument appears literally in the clause or is assigned using =/2 before the meta-call. That is, the following fragment is processed correctly:

``` code
        ...,
        Goal = prove(Theory),
        forall(current_theory(Theory),
               Goal)),
```

But, the calls to prove_simple/1 and prove_complex/1 in the example below are *not* discovered by the analysis and therefore the modules that define these predicates must be loaded explicitly using [use_module/\[1,2\]](import.html#use_module/1).

``` code
        ...,
        member(Goal, [ prove_simple(Theory),
                       prove_complex(Theory)
                     ]),
        forall(current_theory(Theory),
               Goal)),
```

It is good practice to use [gxref/0](xref.html#gxref/0) to make sure that the program has sufficient declarations such that the analysis tools can verify that all required predicates can be resolved and that all code is called. See [meta_predicate/1](metapred.html#meta_predicate/1), [dynamic/1](dynamic.html#dynamic/1), [public/1](dynamic.html#public/1) and [prolog:called_by/2](prologxref.html#prolog:called_by/2).

**volatile** `+Name/Arity, ...`  
Declare that the clauses of specified predicates should **not** be saved to the program. The volatile declaration is normally used to prevent the clauses of dynamic predicates that represent data for the current session from being saved in the state file.

### 14.2.2 Limitations of qsave_program

There are three areas that require special attention when using [qsave_program/\[1,2\]](saved-states.html#qsave_program/1).

- If the program is an embedded Prolog application or uses the foreign language interface, care has to be taken to restore the appropriate foreign context. See [section 14.2.3](saved-states.html#sec:14.2.3) for details.

- If the program uses directives (`:- goal.` lines) that perform other actions than setting predicate attributes ([dynamic/1](dynamic.html#dynamic/1), [volatile/1](saved-states.html#volatile/1), etc.) or loading files ([use_module/1](import.html#use_module/1), etc.). Goals that need to be executed when the state is started must use [initialization/1](consulting.html#initialization/1) (ISO standard) or [initialization/2](consulting.html#initialization/2) (SWI extension that provides more control over when the goal is executed). For example, [initialization/2](consulting.html#initialization/2) can be used to start the application:

  ``` code
  :- initialization(go, main).
  ```

- *Blobs* used as references to the database (see [clause/3](examineprog.html#clause/3), [recorded/3](db.html#recorded/3)), streams, threads, etc. can not be saved. This implies that (dynamic) clauses may not contain such references at the moment the [qsave_program/2](saved-states.html#qsave_program/2) is called. Note that the required foreign context (stream, etc.) cannot be present in the state anyway, making it pointless to save such references. An attempt to save such objects results in a warning.

  The [volatile/1](saved-states.html#volatile/1) directive may be used to prevent saving the clauses of predicates that hold such references. The saved program must reinitialise such references using the normal program initialization techniques: use [initialization/1](consulting.html#initialization/1),2 directives, explicitly create them by the entry point or make the various components recreate the context lazily when required.

- *Blobs* that properly implement the [save()](foreigninclude.html#save()) and [load()](foreigninclude.html#load()) callbacks can be saved and restored. By default a blob is saved as an array of bytes, of the internal form of the blob. This means that any saved program using such a blob is probably not portable to a different architecture.

### 14.2.3 Runtimes and Foreign Code

Many applications use packages that include foreign language components compiled to shared objects or DLLs. This code is normally loaded using [use_foreign_library/1](foreignlink.html#use_foreign_library/1) and the `foreign` file search path. Below is an example from the `socket` library.

``` code
:- use_foreign_library(foreign(socket)).
```

There are two options to handle shared objects in runtime applications. The first is to use the `foreign(save)` option of [qsave_program/2](saved-states.html#qsave_program/2) or the **--foreign=save** commandline option. This causes the dependent shared objects to be included into the resource archive. The [use_foreign_library/1](foreignlink.html#use_foreign_library/1) directive first attempts to find the foreign file in the resource archive. Alternatively, the shared objects may be placed in a directory that is distributed with the application. In this cases the file search path `foreign` must be setup to point at this directory. For example, we can place the shared objects in the same directory as the executable using the definition below. This may be refined further by adding subdirectories depending on the architecture as available from the Prolog flag [arch](flags.html#flag:arch).

``` code
:- multifile user:file_search_path/2.

user:file_search_path(foreign, Dir) :-
    current_prolog_flag(executable, Exe),
    file_directory_name(Exe, Dir).
```
