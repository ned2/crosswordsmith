
## 2.4 Command line options

SWI-Prolog can be executed in one of the following modes:

**`swipl --help`**  
**`swipl --version`**  
**`swipl --arch`**  
**`swipl --dump-runtime-variables`**  
These options must appear as only option. They cause Prolog to print an informational message and exit. See [section 2.4.1](cmdline.html#sec:2.4.1).

**`swipl` \[`option` ...\] `script-file` \[`arg` ...\]**  
These arguments are passed on Unix systems if file that starts with `#!/path/to/executable` \[`option` ...\] is executed. Arguments after the script file are made available in the Prolog flag [argv](flags.html#flag:argv).

**`swipl` \[`option` ...\] `prolog-file` ... \[\[`--`\] `arg` ...\]**  
This is the normal way to start Prolog. The options are described in [section 2.4.2](cmdline.html#sec:2.4.2), [section 2.4.3](cmdline.html#sec:2.4.3) and [section 2.4.4](cmdline.html#sec:2.4.4). The Prolog flag [argv](flags.html#flag:argv) provides access to `arg` ... If the `options` are followed by one or more Prolog file names (i.e., names with extension `.pl`, `.prolog` or (on Windows) the user preferred extension registered during installation), these files are loaded. The first file is registered in the Prolog flag [associated_file](flags.html#flag:associated_file). In addition, **pl-win\[.exe\]** switches to the directory in which this primary source file is located using [working_directory/2](files.html#working_directory/2).

**`swipl` -o `output` -c `prolog-file` ...**  
The **-c** option is used to compile a set of Prolog files into an executable. See [section 2.4.5](cmdline.html#sec:2.4.5).

**`swipl` -o `output` -b `bootfile` `prolog-file` ...**  
Bootstrap compilation. See [section 2.4.6](cmdline.html#sec:2.4.6).

### 2.4.1 Informational command line options

**--arch**  
When given as the only option, it prints the architecture identifier (see Prolog flag [arch](flags.html#flag:arch)) and exits. See also **--dump-runtime-variables**.

**--dump-runtime-variables** `[=format]`  
When given as the only option, it prints a sequence of variable settings that can be used in shell scripts to deal with Prolog parameters. This feature is also used by **swipl-ld** (see [section 12.5](plld.html#sec:12.5)). Below is a typical example of using this feature.

``` code
eval `swipl --dump-runtime-variables`
cc -I$PLBASE/include -L$PLBASE/lib/$PLARCH ...
```

The option can be followed by `=sh` to dump in POSIX shell format (default) or `=cmd` to dump in MS-Windows **cmd.exe** compatible format.

**--help**  
When given as the only option, it summarises the most important options.

**--version**  
When given as the only option, it summarises the version and the architecture identifier.

**--abi-version**  
Print a key (string) that represents the binary compatibility on a number of aspects. See [section 2.21](abi-versions.html#sec:2.21).

### 2.4.2 Command line options for running Prolog

Note that *boolean options* may be written as `--name` (true), `--noname` or `--no-name` (false). They are written as `--no-name` below as the default is‘true’.

**-D** `name[=value]`  
Set the Prolog flag `name` to `value`. The flags are set immediately after loading the initial *saved state*. If the flag is already defined, `value` is converted to the type of the flag. If the flag is undefined it is set to a number of `value` represents a number and an atom otherwise. If no `=``value` is given, a Boolean value is used. If `name` is `no-``flag`, `flag` is set to `false`. Otherwise, the flag `name` is set to `true`. The `name[=value]` may follow the `-D` immediately or appear as the next commandline argument.

Note that many of the commandline options are reflected by a Prolog flag. We intend to handle these as synonyms. Currently, some of the commandline flags affect the Prolog initialization before loading the saved state has completed, while other may not be changed after Prolog initialization. For example, future versions will support `-Dhome=dir` to change the notion of the Prolog installation directory.

**--debug-on-interrupt**  
Enable debugging on an interrupt signal (Control-C, `SIGINT`) immediately. Normally debugging on interrupt is enabled when entering the interactive toplevel. This flag can be used to start the debugger on an interrupt while executing goals from **-g** or [initialization/\[1,2\]](consulting.html#initialization/1). See also the Prolog flag [debug_on_interrupt](flags.html#flag:debug_on_interrupt).

**--home\[=DIR\]**  
Use `DIR` as home directory. See [section 12.6](findhome.html#sec:12.6) for details. If `DIR` is omitted, the found location is printed and the process exits. If the location cannot be found an error is printed and the process exits with status 1. If the home directory is set using this option, the environment variable `SWI_HOME_DIR` holding the specified directory is added to the process.

**--quiet**  
Set the Prolog flag [verbose](flags.html#flag:verbose) to `silent`, suppressing informational and banner messages. Also available as **-q**.

**--no-debug**  
Disable debugging. See the [current_prolog_flag/2](flags.html#current_prolog_flag/2) flag [generate_debug_info](flags.html#flag:generate_debug_info) for details.

**--no-signals**  
Inhibit any signal handling by Prolog, a property that is sometimes desirable for embedded applications. This option sets the flag [signals](flags.html#flag:signals) to `false`. See [section 12.4.25.1](foreigninclude.html#sec:12.4.25.1) for details. Note that the handler to unblock system calls is still installed. This can be prevented using `--sigalert=0` additionally. See **--sigalert**.

**--no-threads**  
Disable threading for the multi-threaded version at runtime. See also the flags [threads](flags.html#flag:threads) and [gc_thread](flags.html#flag:gc_thread).

**--no-packs**  
Do *not* attach extension packages (add-ons). See also [attach_packs/0](pack-attach.html#attach_packs/0) and the Prolog flag [packs](flags.html#flag:packs).

**--no-pce**  
Enable/disable the xpce GUI subsystem. The default is to make it available as autoload component if it is installed and the system can access the graphics. Using `--pce` load the xpce system in user space and `--no-pce` makes it unavailable in the session.

**--on-error** `=style`  
How to handle on errors. See the Prolog flag [on_error](flags.html#flag:on_error) for details.

**--on-warning** `=style`  
How to handle on warnings. See the Prolog flag [on_warning](flags.html#flag:on_warning) for details.

**--pldoc** `[=port]`  
Start the PlDoc documentation system on a free network port and launch the user's browser on `http://localhost:``port`. If `port` is specified, the server is started at the given port and the browser is *not* launched.

**--sigalert=NUM**  
Use signal `NUM` (1 ... 31) for alerting a thread. This is needed to make [thread_signal/2](threadcom.html#thread_signal/2), and derived Prolog signal handling act immediately when the target thread is blocked on an interruptible system call (e.g., [sleep/1](miscpreds.html#sleep/1), read/write to most devices). The default is to use `SIGUSR2`. If `NUM` is 0 (zero), this handler is not installed. See [prolog_alert_signal/2](signal.html#prolog_alert_signal/2) to query or modify this value at runtime.

**--no-tty**  
Unix only. Switches controlling the terminal for allowing single-character commands to the tracer and [get_single_char/1](chario.html#get_single_char/1). By default, manipulating the terminal is enabled unless the system detects it is not connected to a terminal or it is running as a GNU-Emacs inferior process. See also [tty_control](flags.html#flag:tty_control).

**--win-app**  
This option is available only in **swipl-win.exe** and is used for the start-menu item. If causes **plwin** to start in the folder `...\My Documents\Prolog` or local equivalent thereof (see [win_folder/2](system.html#win_folder/2)). The `Prolog` subdirectory is created if it does not exist.

**-O**  
Optimised compilation. See [current_prolog_flag/2](flags.html#current_prolog_flag/2) flag [optimise](flags.html#flag:optimise) for details.

**-l** `file`  
Load `file`. This flag provides compatibility with some other Prolog systems.^(10YAP, SICStus) It is used in SWI-Prolog to skip the program initialization specified using [initialization/2](consulting.html#initialization/2) directives. See also [section 2.11.1.1](compilation.html#sec:2.11.1.1), and [initialize/0](consulting.html#initialize/0).

**-s** `file`  
Use `file` as a script file. The script file is loaded after the initialisation file specified with the **-f** `file` option. Unlike **-f** `file`, using **-s** does not stop Prolog from loading the personal initialisation file.

**-f** `file`  
Use `file` as initialisation file instead of the default `init.pl`.‘**-f** `none`’stops SWI-Prolog from searching for a startup file. This option can be used as an alternative to **-s** `file` that stops Prolog from loading the personal initialisation file. See also [section 2.2](initfile.html#sec:2.2).

**-F** `script`  
Select a startup script from the SWI-Prolog home directory. The script file is named `<``script``>.rc`. The default `script` name is deduced from the executable, taking the leading alphanumerical characters (letters, digits and underscore) from the program name. **-F** `none` stops looking for a script. Intended for simple management of slightly different versions. One could, for example, write a script `iso.rc` and then select ISO compatibility mode using `pl -F iso` or make a link from **iso-pl** to **pl**.

**-x** `bootfile`  
Boot from `bootfile` instead of the system's default boot file. A boot file is a file resulting from a Prolog compilation using the **-b** or **-c** option or a program saved using [qsave_program/\[1,2\]](saved-states.html#qsave_program/1).

**-p** `alias=path1[:path2 ...`  
Define a path alias for file_search_path. `alias` is the name of the alias, and arg path1 ... is a list of values for the alias. On Windows the list separator is `;`. On other systems it is `:`. A value is either a term of the form alias(value) or pathname. The computed aliases are added to [file_search_path/2](consulting.html#file_search_path/2) using [asserta/1](db.html#asserta/1), so they precede predefined values for the alias. See [file_search_path/2](consulting.html#file_search_path/2) for details on using this file location mechanism.

**--traditional**  
This flag disables the most important extensions of SWI-Prolog version 7 (see [section 5](extensions.html#sec:5)) that introduce incompatibilities with earlier versions. In particular, lists are represented in the traditional way, double quoted text is represented by a list of character codes and the functional notation on dicts is not supported. Dicts as a syntactic entity, and the predicates that act on them, are still supported if this flag is present.

**--**  
Stops scanning for more arguments, so you can pass arguments for your application after this one. See [current_prolog_flag/2](flags.html#current_prolog_flag/2) using the flag [argv](flags.html#flag:argv) for obtaining the command line arguments.

### 2.4.3 Controlling the stack sizes

As of version 7.7.14 the stacks are no longer limited individually. Instead, only the combined size is limited. Note that 32 bit systems still pose a 128Mb limit. See [section 2.19.1](limits.html#sec:2.19.1). The combined limit is by default 1Gb on 64 bit machines and 512Mb on 32 bit machines.

For example, to limit the stacks to 32Gb use the command below. Note that the stack limits apply *per thread*. Individual threads may be controlled using the `stack_limit(+Bytes)` option of thread_create. Any thread can call `set_prolog_flag(stack_limit, Limit)` (see [stack_limit](flags.html#flag:stack_limit)) to adjust the stack limit. This limit is inherited by threads created from this thread.

``` code
$ swipl --stack-limit=32g
```

**--stack-limit**`=size[bkmg]`  
Limit the combined size of the Prolog stacks to the indicated `size`. The suffix specifies the value as *bytes*, *Kbytes*, *Mbytes* or *Gbytes*.

**--table-space**`=size[bkmg]`  
Limit for the `table space`. This is where tries holding memoized^(11The letter M is used because the T was already in use. It is a mnemonic for **M**emoizing.) answers for *tabling* are stored. The default is 1Gb on 64 bit machines and 512Mb on 32 bit machines. See the Prolog flag [table_space](flags.html#flag:table_space).

**--shared-table-space**`=size[bkmg]`  
Limit for the table space for *shared* tables. See [section 7.9](tabling-shared.html#sec:7.9).

### 2.4.4 Running goals from the command line

**-g** `goal`  
`Goal` is executed just before entering the top level. This option may appear multiple times. See [section 2.3](initgoal.html#sec:2.3) for details. If no initialization goal is present the system calls [version/0](printmsg.html#version/0) to print the welcome message. The welcome message can be suppressed with **--quiet**, but also with **-g** `true`. `goal` can be a complex term. In this case quotes are normally needed to protect it from being expanded by the shell. A safe way to run a goal non-interactively is below. If go/0 succeeds **-g** `halt` causes the process to stop with exit code 0. If it fails, the exit code is 1; and if it raises an exception, the exit code is 2.

``` code
% swipl <options> -g go -g halt
```

**-t** `goal`  
Use `goal` as interactive top level instead of the default goal [prolog/0](toplevel.html#prolog/0). The `goal` can be a complex term. If the top-level goal succeeds SWI-Prolog exits with status 0. If it fails the exit status is 1. If the top level raises an exception, this is printed as an uncaught error and the top level is *restarted*. This flag also determines the goal started by [break/0](toplevel.html#break/0) and [abort/0](toplevel.html#abort/0). If you want to prevent the user from entering interactive mode, start the application with‘**-g** `goal` **-t** `halt`’.

### 2.4.5 Compilation options

**-c** `file ...`  
Compile files into an‘intermediate code file’. See [section 2.11](compilation.html#sec:2.11).

**-o** `output`  
Used in combination with **-c** or **-b** to determine output file for compilation.

### 2.4.6 Maintenance options

The following options are for system maintenance. They are given for reference only.

**-b** `initfile ...`**`-c`**` file ...`  
Boot compilation. `initfile ...` are compiled by the C-written bootstrap compiler, `file ...` by the normal Prolog compiler. System maintenance only.

**-d** `token1,token2,...`  
Print debug messages for DEBUG statements tagged with one of the indicated tokens. Only has effect if the system is compiled with the `-DO_DEBUG` flag. System maintenance only.
