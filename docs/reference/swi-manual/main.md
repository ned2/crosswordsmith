
## A.27 library(main): Provide entry point for scripts

See also  
\- `library(prolog_stack)` to force backtraces in case of an uncaught exception.  
- XPCE users should have a look at `library(pce_main)`, which starts the GUI and processes events until all windows have gone.

This library is intended for supporting PrologScript on Unix using the `#!` magic sequence for scripts using commandline options. The entry point [main/0](main.html#main/0) calls the user-supplied predicate main/1 passing a list of commandline options. Below is a simle `echo` implementation in Prolog.

``` code
#!/usr/bin/env swipl

:- initialization(main, main).

main(Argv) :-
    echo(Argv).

echo([]) :- nl.
echo([Last]) :- !,
    write(Last), nl.
echo([H|T]) :-
    write(H), write(' '),
    echo(T).
```

**main**  
Call main/1 using the passed command-line arguments. Before calling main/1 this predicate installs a signal handler for `SIGINT` (Control-C) that terminates the process with status 1.

When [main/0](main.html#main/0) is called interactively it simply calls main/1 with the arguments. This allows for debugging scripts as follows:

``` code
$ swipl -l script.pl -- arg ...
?- gspy(suspect/1).     % setup debugging
?- main.            % run program
```

\[det\]**argv_options**(`:Argv, -Positional, -Options`)  
Parse command line arguments. This predicate acts in one of two modes.

- If the calling module defines opt_type/3, full featured parsing with long and short options, type conversion and help is provided.
- If opt_type/3 is not defined, only unguided transformation using long options is supported. See argv_untyped_options/3 for details.

When **guided**, three predicates are called in the calling module. opt_type/3 **must** be defined, the others need not. Note that these three predicates *may* be defined as *multifile* to allow multiple modules contributing to the provided commandline options. Defining them as *discontiguous* allows for creating blocks that describe a group of related options.

**opt_type**(`Opt, Name, Type`)  
Defines `Opt` to add an option `Name`(Value), where Value statisfies `Type`. `Opt` does not include the leading `-`. A single character implies a short option, multiple a long option. Long options use `_` as *word separator*, user options may use either `_` or `-`. `Type` is one of:

`A` **`|`** `B`  
Disjunctive type. Disjunction can be used create long options with optional values. For example, using the type `nonneg|boolean`, for an option `http` handles `--http` as `http(true)`, `--no-http` as `http(false)` and `--http=3000` as `http(3000)`. Note that with an optional boolean a option is considered boolean unless it has a value written as `--longopt=value`.

**boolean**(`Default`)  
**boolean**  
Boolean options are special. They do not take a value except for when using the long `--opt=value` notation. This explicit value specification converts `true`, `True`, `TRUE`, `on`, `On`, `ON`, `1` and the obvious false equivalents to Prolog `true` or `false`. If the option is specified, Default is used. If `--no-opt` or `--noopt` is used, the inverse of Default is used.

**integer**  
Argument is converted to an integer

**float**  
Argument is converted to a float. User may specify an integer

**nonneg**  
As `integer`. Requires value `>=` 0.

**natural**  
As `integer`. Requires value `>=` 1.

**number**  
Any number (integer, float, rational).

**between**(`Low, High`)  
If both one of `Low` and `High` is a float, convert as `float`, else convert as `integer`. Then check the range.

**atom**  
No conversion

**oneof**(`List`)  
As `atom`, but requires the value to be a member of `List` (*enum* type).

**string**  
Convert to a SWI-Prolog string

**file**  
Convert to a file name in Prolog canonical notation using [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2).

**directory**  
Convert to a file name in Prolog canonical notation using [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2). No checking is done and thus this type is the same as `file`

**file**(`Access`)  
As `file`, and check access using [access_file/2](files.html#access_file/2). A value `-` is not checked for access, assuming the application handles this as standard input or output.

**directory**(`Access`)  
As `directory`, and check access. `Access` is one of `read` `write` or `create`. In the latter case the parent directory must exist and have write access.

**term**  
Parse option value to a Prolog term.

**term**(`+Options`)  
As `term`, but passes `Options` to [term_string/3](string.html#term_string/3). If the option `variable_names(Bindings)` is given the option value is set to the *pair* `Term-Bindings`.

**opt_help**(`Name, HelpString`)  
Help string used by [argv_usage/1](main.html#argv_usage/1).

**opt_meta**(`Name, Meta`)  
If a typed argument is required this defines the placeholder in the help message. The default is the uppercase version of the type *functor name*. This produces the `FILE` in e.g. `-f FILE`.

By default, `-h`, `-?` and `--help` are bound to help. If `opt_type(Opt, help, boolean)` is true for some `Opt`, the default help binding and help message are disabled and the normal user rules apply. In particular, the user should also provide a rule for `opt_help(help, String)`.

\[det\]**argv_options**(`:Argv, -Positional, -Options, +ParseOptions`)  
As [argv_options/3](main.html#argv_options/3) in **guided** mode, Currently this version allows parsing argument options throwing an exception rather than calling [halt/1](toplevel.html#halt/1) by passing an empty list to `ParseOptions`. `ParseOptions`:

**on_error**(`+Goal`)  
If `Goal` is `halt(Code)`, exit with Code. Other goals are currently not supported.

**options_after_arguments**(`+Boolean`)  
If `false` (default `true`), stop parsing after the first positional argument, returning options that follow this argument as positional arguments. E.g, `-x file -y` results in positional arguments `[file, '-y']`

**unknown_option**(`+Mode`)  
One of `error` (default) or `pass`. Using `pass`, the option is passed in `Positional`. Multi-flag short options may be processed partially. For example, if `-v` is defined and `-iv` is in `Argv`, `Positional` receives `'-i'` and the option defined with `-v` is added to `Options`.

To be done  
When passing unknown options we may wish to process multi-flag options as a whole or not at all rather than passing the unknown flags.

\[det\]**argv_usage**(`:Level`)  
Use [print_message/2](printmsg.html#print_message/2) to print a usage message at `Level`. To print the message as plain text indefault color, use `debug`. Other meaningful options are `informational` or `warning`. The help page consists of four sections, two of which are optional:

1.  The **header** is created from `opt_help(help(header), String)`. It is optional.
2.  The **usage** is added by default. The part behind `Usage: <command>` is by default `[options]` and can be overruled using `opt_help(help(usage), String)`.
3.  The actual option descriptions. The options are presented in the order they are defined in opt_type/3. Subsequent options for the same *destination* (option name) are joined with the first.
4.  The *footer\_* is created from `opt_help(help(footer), String)`. It is optional.

The help provided by `help(header)`, `help(usage)` and `help(footer)` are either a simple string or a list of elements as defined by [print_message_lines/3](printmsg.html#print_message_lines/3). In the latter case, the construct `\Callable` can be used to call a DCG rule in the module from which the user calls [argv_options/3](main.html#argv_options/3). For example, we can add a bold title using

``` code
opt_help(help(header), [ansi(bold, '~w', ['My title'])]).
```

\[det\]**cli_parse_debug_options**(`+OptionsIn, -Options`)  
Parse certain commandline options for debugging and development purposes. `Options` processed are below. Note that the option argument is an atom such that these options may be activated as e.g., `--debug='http(_)'`.

**debug**(`Topic`)  
Call `debug(Topic)`. See [debug/1](debug.html#debug/1) and [debug/3](debug.html#debug/3).

**spy**(`Predicate`)  
Place a spy-point on `Predicate`.

**gspy**(`Predicate`)  
As spy using the graphical debugger. See [tspy/1](threadutil.html#tspy/1).

**interactive**(`true`)  
Start the Prolog toplevel after main/1 completes.

**cli_debug_opt_type**(`-Flag, -Option, -Type`)  
**cli_debug_opt_help**(`-Option, -Message`)  
**cli_debug_opt_meta**(`-Option, -Arg`)  
Implements opt_type/3, [opt_help/2](optparse.html#opt_help/2) and opt_meta/2 for debug arguments. Applications that wish to use these features can call these predicates from their own hook. Fot example:

``` code
opt_type(..., ..., ...).    % application types
opt_type(Flag, Opt, Type) :-
    cli_debug_opt_type(Flag, Opt, Type).
% similar for opt_help/2 and opt_meta/2

main(Argv) :-
    argv_options(Argv, Positional, Options0),
    cli_parse_debug_options(Options0, Options),
    ...
```

**cli_enable_development_system**  
Re-enable the development environment. Currently re-enables xpce if this was loaded, but not initialised and causes the interactive toplevel to be re-enabled.

This predicate may be called from main/1 to enter the Prolog toplevel rather than terminating the application after main/1 completes.
