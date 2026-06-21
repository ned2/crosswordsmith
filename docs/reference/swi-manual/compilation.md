
## 2.11 Loading and running projects

Most Prolog programs are split over multiple files organized in a directory and optionally multiple subdirectories. Typically all files are Prolog *module* files. See [section 6](modules.html#sec:6). Typically, the directory contains a file, often called `load.pl`, that loads all other files (modules) using [use_module/\[1,2\]](import.html#use_module/1) or, for projects that do not use modules, using [ensure_loaded/1](consulting.html#ensure_loaded/1).

If the project is an application (rather than a library), there are several ways to start it. One option is by using the commandline option **-g** `goal`. The classical Prolog way is by using an [initialization/1](consulting.html#initialization/1) *directive*. Th problem with the latter is that such directives are both used for runtime initialization in modules and starting the application while it is hard to control the order in which they are executed. For this reason, SWI-Prolog introduced [initialization/2](consulting.html#initialization/2), adding an argument that specifies the role and (indirectly) the order of initialization. The application entry point is now declared using

``` code
:- initialization(start, main).

start :-
    ...
```

Using these conventions we may run the application using this command line, where `option ...` are Prolog options to control e.g., memory limits. Typically, none are required. `arg ...` are made available to the program using the Prolog flag [argv](flags.html#flag:argv).

``` code
% swipl [option ...] load.pl [arg ...]
```

To merely load the code without running the application, provided the entry point is started using the [initialization/2](consulting.html#initialization/2) directive described above, we can use the **-l**. After loading we can debug and/or edit the application.

``` code
% swipl [option ...] -l load.pl [arg ...]
```

Rather than just using start/0 as above, applications typically use [main/0](main.html#main/0) from the library `library(main)`. The [main/0](main.html#main/0) predicate prepares for non-development usage and calls main/1 with the application `argv` (command line arguments). These are normally processed into *positional arguments* and *options* using argv_options/2 from the same library.

While the above works fine when using Prolog from the commandline, it is less suitable for scenarios that make it hard to control the SWI-Prolog commandline which as using **swipl-win** or running Prolog under some IDE such as Emacs. Loading a program that uses the above [initialization/2](consulting.html#initialization/2) directive into the toplevel using

``` code
?- [load].
```

does **not** start the entry point. Opening a `..pl` file using **swipl-win** does start the entry point.

### 2.11.1 Running an application

There are various options if you want to make your program ready for real usage. The best choice depends on whether the program is to be used only on machines holding the SWI-Prolog development system, the size of the program, and the operating system (Unix vs. Windows). There are four options

- On Unix-like systems one can use the *shebang* magic sequence to turn a Prolog source into an executable. See [section 2.11.1.1](compilation.html#sec:2.11.1.1).

- On any system you can use a *shell script* (Unix **sh** or Windows **cmd**) script to start the application. See [section 2.11.1.2](compilation.html#sec:2.11.1.2).

- On any system you can create a *saved state* that consists of the virtual machine code and a startup sequence. Saved states can be stand-alone and with some precautions they can work without SWI-Prolog itself installed. They start fast, but they are big and creating a state from a program that uses native code extensions and (file) resources is not trivial while details depend on the OS and required resources. See [section 2.11.1.3](compilation.html#sec:2.11.1.3).

- On any system you can add a Prolog file to a designated directory and allow it to be started using

  ``` code
  swipl name [arg ...]
  ```

  New commands can be added to the Prolog installation, by Prolog *packs*, in a user specific directory or in a system-wide directory. See [section 2.11.1.5](compilation.html#sec:2.11.1.5).

#### 2.11.1.1 Using PrologScript

A Prolog source file can be used directly as a Unix program using the Unix `#!` magic start. The Unix `#!` magic is allowed because if the first letter of a Prolog file is `#`, the first line is treated as a comment.^(12The `#`-sign can be the legal start of a normal Prolog clause. In the unlikely case this is required, leave the first line blank or add a header comment.) To create a Prolog script, use one of the two alternatives below as first line. The first can be used to bind a script to a specific Prolog installation, while the latter uses the default prolog installed in `$PATH`.

``` code
#!/path/to/swipl
#!/usr/bin/env swipl
```

The interpretation of arguments to the executable in the *HashBang* line differs between Unix-derived systems. For portability, the `#!` must be followed immediately with an absolute path to the executable and should have none or one argument. Neither the executable path, nor the argument shall use quotes or spaces. When started this way, the Prolog flag [argv](flags.html#flag:argv) contains the command line arguments that follow the script invocation.

Starting with version 7.5.8, [initialization/2](consulting.html#initialization/2) support the `When` options `program` and `main`, allowing for the following definition of a Prolog script that evaluates an arithmetic expression on the command line. Note that [main/0](main.html#main/0) is defined lib the library `library(main)`. It calls main/1 with the command line arguments after disabling signal handling.

``` code
#!/usr/bin/env swipl

:- initialization(main, main).

main(Argv) :-
    atomic_list_concat(Argv, ' ', SingleArg),
    term_to_atom(Term, SingleArg),
    Val is Term,
    format('~w~n', [Val]).
```

And here are two example runs:

``` code
% ./eval 1+2
3
% ./eval foo
ERROR: is/2: Arithmetic: `foo/0' is not a function
```

Prolog script may be launched for debugging or inspection purposes using the **-l** or **-t**. For example, **-l** merely loads the script, ignoring `main` and `program` initialization.

``` code
swipl -l eval 1+1
<banner>

?- main.
2
true.

?-
```

We can also force the program to enter the interactive toplevel after the application is completed using `-t prolog`:

``` code
swipl -t prolog eval 1+1
2
?-
```

The Windows version simply ignores the `#!` line.^(13Older versions extracted command line arguments from the *HashBang* line. As of version 5.9 all relevant setup can be achieved using *directives*. Due to the compatibility issues around *HashBang* line processing, we decided to remove it completely.)

#### 2.11.1.2 Creating a shell script

With the introduction of *PrologScript* (see [section 2.11.1.1](compilation.html#sec:2.11.1.1)), using shell scripts as explained in this section has become redundant for most applications.

Especially on Unix systems and not-too-large applications, writing a shell script that simply loads your application and calls the entry point is often a good choice. A skeleton for the script is given below, followed by the Prolog code to obtain the program arguments. See library `library(main)` and [argv_options/3](main.html#argv_options/3) for details.

``` code
#!/bin/sh

base=<absolute-path-to-source>
SWIPL=swipl

exec $SWIPL "$base/load.pl" -- "$@"
```

``` code
:- use_module(library(main)).
:- initialization(main,main).

main(Argv) :-
    argv_options(Argv, Positional, Options),
    go(Positional, Options).

go(Positional, Options) :-
    ...
```

On Windows systems, similar behaviour can be achieved by creating a shortcut to Prolog, passing the proper options or writing a `.bat` file.

#### 2.11.1.3 Creating a saved state

For larger programs, as well as for programs that are required to run on systems that do not have the SWI-Prolog development system installed, creating a saved state is the best solution. A saved state is created using [qsave_program/\[1,2\]](saved-states.html#qsave_program/1) or the **-c** command line option. A saved state is a file containing machine-independent^(14The saved state does not depend on the CPU instruction set or endianness. Saved states for 32- and 64-bits are not compatible. Typically, saved states only run on the same version of Prolog on which they have been created.) intermediate code in a format dedicated for fast loading. Optionally, the emulator may be integrated in the saved state, creating a single file, but machine-dependent, executable. This process is described in [chapter 14](runtime.html#sec:14).

#### 2.11.1.4 Compilation using the -c command line option

This mechanism loads a series of Prolog source files and then creates a saved state as [qsave_program/2](saved-states.html#qsave_program/2) does. The command syntax is:

``` code
% swipl [option ...] [-o output] -c file.pl ...
```

The `options` argument are options to [qsave_program/2](saved-states.html#qsave_program/2) written in the format below. The option names and their values are described with [qsave_program/2](saved-states.html#qsave_program/2).

> `--`*option-name*`=`*option-value*

For example, to create a stand-alone executable that starts by executing main/0 and for which the source is loaded through `load.pl`, use the command

``` code
% swipl --goal=main --stand_alone=true -o myprog -c load.pl
```

This performs exactly the same as executing

``` code
% swipl
<banner>

?- [load].
?- qsave_program(myprog,
                 [ goal(main),
                   stand_alone(true)
                 ]).
?- halt.
```

#### 2.11.1.5 SWI-Prolog app scripts

As of version 9.1.18, SWI-Prolog allows starting an application using the command below.

``` code
swipl [option ...] [path:]name [arg ...]
```

This command line first processes Prolog options described in [section 2.4](cmdline.html#sec:2.4). Note that most standard Prolog commandline options are not relevant. The **-f** defaults to `none`, which implies that the user init file is by default not loaded. If an application wishes to load the user init file, it should load `user_app_config(init)` if this file exists (see [exists_source/1](consulting.html#exists_source/1)).

Next, it locates `path(name)` using SWI-Prolog's file search mechanism defined by [absolute_file_name/3](files.html#absolute_file_name/3). After loading this file it finds the last goal registered for `main` using [initialization/2](consulting.html#initialization/2) as described in [section 2.11](compilation.html#sec:2.11) - if there is no initialization directive for `main`, the program terminates with an error. By default, the application terminates after the entry point terminates. The entry point may enable the interactive Prolog REPL loop by calling [cli_enable_development_system/0](main.html#cli_enable_development_system/0). Other forms of the [initialization/2](consulting.html#initialization/2) directive are also allowed, in addition to‘main\`.

All command line options after `[path:]name` are accessible in the Prolog flag [argv](flags.html#flag:argv).

The optional `path` defaults to `app`. By default, apps are searched in the directories below. See [file_search_path/2](consulting.html#file_search_path/2) for details.

1.  The `app` directory of the SWI-Prolog installation
2.  User and site configuration. On POSIX systems using the XDG file name conventions, this is normally ` /.local/share/swi-prolog/app/` and `/usr/share/swi-prolog/app`.
3.  The `app` directory of a Prolog *pack*.

The following apps are provided by the installation

**app**  
Print information on installed apps. For example, to list all available apps, run

``` code
swipl app list
```

**pack**  
Command line driven management of Prolog packs. This is a front-end to the Prolog library `library(prolog_pack)`. For example, to find packages related to *type*, use the command below.

``` code
swipl pack find type
```
