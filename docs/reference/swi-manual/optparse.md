
## A.33 library(optparse): command line parsing

author  
Marcus Uneson

version  
0.20 (2011-04-27)

To be done  
: validation? e.g, numbers; file path existence; one-out-of-a-set-of-atoms

This module helps in building a command-line interface to an application. In particular, it provides functions that take an option specification and a list of atoms, probably given to the program on the command line, and return a parsed representation (a list of the customary Key(Val) by default; or optionally, a list of Func(Key, Val) terms in the style of [current_prolog_flag/2](flags.html#current_prolog_flag/2)). It can also synthesize a simple help text from the options specification.

The terminology in the following is partly borrowed from python, see [http://docs.python.org/library/optparse.html\\terminology](http://docs.python.org/library/optparse.html\#terminology) . Very briefly, *arguments* is what you provide on the command line and for many prologs show up as a list of atoms `Args` in `current_prolog_flag(argv, Args)`. For a typical prolog incantation, they can be divided into

- *runtime arguments*, which controls the prolog runtime; conventionally, they are ended by’--’;
- *options*, which are key-value pairs (with a boolean value possibly implicit) intended to control your program in one way or another; and
- *positional arguments*, which is what remains after all runtime arguments and options have been removed (with implicit arguments -- true/false for booleans -- filled in).

Positional arguments are in particular used for mandatory arguments without which your program won't work and for which there are no sensible defaults (e.g,, input file names). Options, by contrast, offer flexibility by letting you change a default setting. Options are optional not only by etymology: this library has no notion of mandatory or required options (see the python docs for other rationales than laziness).

The command-line arguments enter your program as a list of atoms, but the programs perhaps expects booleans, integers, floats or even prolog terms. You tell the parser so by providing an *options specification*. This is just a list of individual option specifications. One of those, in turn, is a list of ground prolog terms in the customary Name(Value) format. The following terms are recognized (any others raise error).

**opt**(`Key`)  
`Key` is what the option later will be accessed by, just like for `current_prolog_flag(Key, Value)`. This term is mandatory (an error is thrown if missing).

**shortflags**(`ListOfFlags`)  
`ListOfFlags` denotes any single-dashed, single letter args specifying the current option (`-s , -K`, etc). Uppercase letters must be quoted. Usually `ListOfFlags` will be a singleton list, but sometimes aliased flags may be convenient.

**longflags**(`ListOfFlags`)  
`ListOfFlags` denotes any double-dashed arguments specifying the current option (`--verbose, --no-debug`, etc). They are basically a more readable alternative to short flags, except

1.  long flags can be specified as `--flag value` or `--flag=value` (but not as `--flagvalue`); short flags as `-f val` or `-fval` (but not `-f=val`)

2.  boolean long flags can be specified as `--bool-flag` or `--bool-flag=true` or `--bool-flag true`; and they can be negated as `--no-bool-flag` or `--bool-flag=false` or `--bool-flag false`.

    Except that shortflags must be single characters, the distinction between long and short is in calling convention, not in namespaces. Thus, if you have `shortflags([v])`, you can use it as `-v2` or `-v 2` or `--v=2` or `--v 2` (but not `-v=2` or `--v2`).

    Shortflags and longflags both default to `[]`. It can be useful to have flagless options -- see example below.

**meta**(`Meta`)  
`Meta` is optional and only relevant for the synthesized usage message and is the name (an atom) of the metasyntactic variable (possibly) appearing in it together with type and default value (e.g, `x:integer=3`, `interest:float=0.11`). It may be useful to have named variables (`x`, `interest`) in case you wish to mention them again in the help text. If not given the `Meta:` part is suppressed -- see example below.

**type**(`Type`)  
`Type` is one of `boolean, atom, integer, float, term`. The corresponding argument will be parsed appropriately. This term is optional; if not given, defaults to `term`.

**default**(`Default`)  
`Default` value. This term is optional; if not given, or if given the special value’\_’, an uninstantiated variable is created (and any type declaration is ignored).

**help**(`Help`)  
`Help` is (usually) an atom of text describing the option in the help text. This term is optional (but obviously strongly recommended for all options which have flags).

Long lines are subject to basic word wrapping -- split on white space, reindent, rejoin. However, you can get more control by supplying the line breaking yourself: rather than a single line of text, you can provide a list of lines (as atoms). If you do, they will be joined with the appropriate indent but otherwise left untouched (see the option `mode` in the example below).

Absence of mandatory option specs or the presence of more than one for a particular option throws an error, as do unknown or incompatible types.

As a concrete example from a fictive application, suppose we want the following options to be read from the command line (long `flag(s)`, short `flag(s)`, meta:type=default, help)

``` code
--mode                  -m     atom=SCAN       data gathering mode,
                                               one of
                                                SCAN: do this
                                                READ: do that
                                                MAKE: make numbers
                                                WAIT: do nothing
--rebuild-cache         -r     boolean=true    rebuild cache in
                                               each iteration
--heisenberg-threshold  -t,-h  float=0.1       heisenberg threshold
--depths, --iters       -i,-d  K:integer=3     stop after K
                                               iterations
--distances                    term=[1,2,3,5]  initial prolog term
--output-file           -o     FILE:atom=_     write output to FILE
--label                 -l     atom=REPORT     report label
--verbosity             -v     V:integer=2     verbosity level,
                                               1 <= V <= 3
```

We may also have some configuration parameters which we currently think not needs to be controlled from the command line, say `path('/some/file/path')`.

This interface is described by the following options specification (order between the specifications of a particular option is irrelevant).

``` code
ExampleOptsSpec =
    [ [opt(mode    ), type(atom), default('SCAN'),
        shortflags([m]),   longflags(['mode'] ),
        help([ 'data gathering mode, one of'
             , '  SCAN: do this'
             , '  READ: do that'
             , '  MAKE: fabricate some numbers'
             , '  WAIT: don''t do anything'])]

    , [opt(cache), type(boolean), default(true),
        shortflags([r]),   longflags(['rebuild-cache']),
        help('rebuild cache in each iteration')]

    , [opt(threshold), type(float), default(0.1),
        shortflags([t,h]),  longflags(['heisenberg-threshold']),
        help('heisenberg threshold')]

    , [opt(depth), meta('K'), type(integer), default(3),
        shortflags([i,d]),longflags([depths,iters]),
        help('stop after K iterations')]

    , [opt(distances), default([1,2,3,5]),
        longflags([distances]),
        help('initial prolog term')]

    , [opt(outfile), meta('FILE'), type(atom),
        shortflags([o]),  longflags(['output-file']),
        help('write output to FILE')]

    , [opt(label), type(atom), default('REPORT'),
        shortflags([l]), longflags([label]),
        help('report label')]

    , [opt(verbose),  meta('V'), type(integer), default(2),
        shortflags([v]),  longflags([verbosity]),
        help('verbosity level, 1 <= V <= 3')]

    , [opt(path), default('/some/file/path/')]
    ].
```

The help text above was accessed by `opt_help(ExamplesOptsSpec, HelpText)`. The options appear in the same order as in the OptsSpec.

Given `ExampleOptsSpec`, a command line (somewhat syntactically inconsistent, in order to demonstrate different calling conventions) may look as follows

``` code
ExampleArgs = [ '-d5'
              , '--heisenberg-threshold', '0.14'
              , '--distances=[1,1,2,3,5,8]'
              , '--iters', '7'
              , '-ooutput.txt'
              , '--rebuild-cache', 'true'
              , 'input.txt'
              , '--verbosity=2'
              ].
```

`opt_parse(ExampleOptsSpec, ExampleArgs, Opts, PositionalArgs)` would then succeed with

``` code
Opts =    [ mode('SCAN')
          , label('REPORT')
          , path('/some/file/path')
          , threshold(0.14)
          , distances([1,1,2,3,5,8])
          , depth(7)
          , outfile('output.txt')
          , cache(true)
          , verbose(2)
          ],
PositionalArgs = ['input.txt'].
```

Note that `path('/some/file/path')` showing up in Opts has a default value (of the implicit type’term’), but no corresponding flags in OptsSpec. Thus it can't be set from the command line. The rest of your program doesn't need to know that, of course. This provides an alternative to the common practice of asserting such hard-coded parameters under a single predicate (for instance `setting(path, '/some/file/path')`), with the advantage that you may seamlessly upgrade them to command-line options, should you one day find this a good idea. Just add an appropriate flag or two and a line of help text. Similarly, suppressing an option in a cluttered interface amounts to commenting out the flags.

[opt_parse/5](optparse.html#opt_parse/5) allows more control through an additional argument list as shown in the example below.

``` code
?- opt_parse(ExampleOptsSpec, ExampleArgs,  Opts, PositionalArgs,
             [ output_functor(appl_config)
             ]).

Opts =    [ appl_config(verbose, 2),
          , appl_config(label, 'REPORT')
          ...
          ]
```

This representation may be preferable with the empty-flag configuration parameter style above (perhaps with asserting appl_config/2).

### A.33.1 Notes and tips

- In the example we were mostly explicit about the types. Since the default is `term`, which subsumes `integer, float, atom`, it may be possible to get away cheaper (e.g., by only giving booleans). However, it is recommended practice to always specify types: parsing becomes more reliable and error messages will be easier to interpret.
- Note that `-sbar` is taken to mean `-s bar`, not `-s -b -a -r`, that is, there is no clustering of flags.
- `-s=foo` is disallowed. The rationale is that although some command-line parsers will silently interpret this as `-s =foo`, this is very seldom what you want. To have an option argument start with’=’(very un-recommended), say so explicitly.
- The example specifies the option `depth` twice: once as `-d5` and once as `--iters 7`. The default when encountering duplicated flags is to `keeplast` (this behaviour can be controlled, by ParseOption duplicated_flags).
- The order of the options returned by the parsing functions is the same as given on the command line, with non-overridden defaults prepended and duplicates removed as in previous item. You should not rely on this, however.
- Unknown flags (not appearing in OptsSpec) will throw errors. This is usually a Good Thing. Sometimes, however, you may wish to pass along flags to an external program (say, one called by [shell/2](system.html#shell/2)), and it means duplicated effort and a maintenance headache to have to specify all possible flags for the external program explicitly (if it even can be done). On the other hand, simply taking all unknown flags as valid makes error checking much less efficient and identification of positional arguments uncertain. A better solution is to collect all arguments intended for passing along to an indirectly called program as a single argument, probably as an atom (if you don't need to inspect them first) or as a prolog term (if you do).

\[det\]**opt_arguments**(`+OptsSpec, -Opts, -PositionalArgs`)  
Extract commandline options according to a specification. Convenience predicate, assuming that command-line arguments can be accessed by [current_prolog_flag/2](flags.html#current_prolog_flag/2) (as in swi-prolog). For other access mechanisms and/or more control, get the args and pass them as a list of atoms to [opt_parse/4](optparse.html#opt_parse/4) or [opt_parse/5](optparse.html#opt_parse/5) instead.

`Opts` is a list of parsed options in the form Key(Value). Dashed args not in `OptsSpec` are not permitted and will raise error (see tip on how to pass unknown flags in the module description). `PositionalArgs` are the remaining non-dashed args after each flag has taken its argument (filling in `true` or `false` for booleans). There are no restrictions on non-dashed arguments and they may go anywhere (although it is good practice to put them last). Any leading arguments for the runtime (up to and including’--’) are discarded.

\[det\]**opt_parse**(`+OptsSpec, +ApplArgs, -Opts, -PositionalArgs`)  
Equivalent to `opt_parse(OptsSpec, ApplArgs, Opts, PositionalArgs, [])`.

\[det\]**opt_parse**(`+OptsSpec, +ApplArgs, -Opts, -PositionalArgs, +ParseOptions`)  
Parse the arguments Args (as list of atoms) according to `OptsSpec`. Any runtime arguments (typically terminated by’--’) are assumed to be removed already.

`Opts` is a list of parsed options in the form Key(Value), or (with the option `functor(Func)` given) in the form Func(Key, Value). Dashed args not in `OptsSpec` are not permitted and will raise error (see tip on how to pass unknown flags in the module description). `PositionalArgs` are the remaining non-dashed args after each flag has taken its argument (filling in `true` or `false` for booleans). There are no restrictions on non-dashed arguments and they may go anywhere (although it is good practice to put them last). `ParseOptions` are

**output_functor**(`Func`)  
Set the functor `Func` of the returned options `Func`(Key,Value). Default is the special value’OPTION’(upper-case), which makes the returned options have form Key(Value).

**duplicated_flags**(`Keep`)  
Controls how to handle options given more than once on the commad line. `Keep` is one of `keepfirst, keeplast, keepall` with the obvious meaning. Default is `keeplast`.

**allow_empty_flag_spec**(`Bool`)  
If true (default), a flag specification is not required (it is allowed that both shortflags and longflags be either `[]` or absent). Flagless options cannot be manipulated from the command line and will not show up in the generated help. This is useful when you have (also) general configuration parameters in your `OptsSpec`, especially if you think they one day might need to be controlled externally. See example in the module overview. `allow_empty_flag_spec(false)` gives the more customary behaviour of raising error on empty flags.

\[det\]**opt_help**(`+OptsSpec, -Help:atom`)  
True when `Help` is a help string synthesized from `OptsSpec`.

\[semidet,multifile\]**parse_type**(`+Type, +Codes:list(code), -Result`)  
Hook to parse option text `Codes` to an object of type `Type`.
