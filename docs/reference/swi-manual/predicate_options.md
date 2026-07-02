
## A.39 library(predicate_options): Declare option-processing of predicates

*Discussions with Jeff Schultz helped shaping this library*

### A.39.1 The strength and weakness of predicate options

Many ISO predicates accept options, e.g., [open/4](IO.html#open/4), [write_term/3](termrw.html#write_term/3). Options offer an attractive alternative to proliferation into many predicates and using high-arity predicates. Properly defined and used, they also form a mechanism for extending the API of both system and application predicates without breaking portability. I.e., previously fixed behaviour can be replaced by dynamic behaviour controlled by an option where the default is the previously defined fixed behaviour. The alternative to using options is to add an additional argument and maintain the previous definition. While a series of predicates with increasing arity is adequate for a small number of additional parameters, the untyped positional argument handling of Prolog quickly makes this unmanageable.

The ISO standard uses the extensibility offered by options by allowing implementations to extend the set of accepted options. While options form a perfect solution to maintain backward portability in a linear development model, it is not well equipped to deal with concurrent branches because

1.  There is no API to find which options are supported in a particular implementation.
2.  While the portability problem caused by a missing predicate in Prolog *A* can easily be solved by implementing this predicate, it is much harder to add processing of an additional option to an already existing predicate.

Different Prolog implementations can be seen as concurrent development branches of the Prolog language. Different sets of supported options pose a serious portability issue. Using an option *O* that establishes the desired behaviour on system *A* leads (on most systems) to an error or system *B*. Porting may require several actions:

- Drop *O* (if the option is not vital, such as the layout options to [write_term/3](termrw.html#write_term/3))
- Replace *O* by *O2* (i.e., a differently named option doing the same)
- Something else (cannot be ported; requires a totally different approach, etc.)

Predicates that process options are particularly a problem when writing a compatibility layer to run programs developed for System *A* on System *B* because complete emulation is often hard, may cause a serious slowdown and is often not needed because the application-to-be-ported only uses options that are shared by all target Prolog implementations. Unfortunately, the consequences of a partial emulation cannot be assessed by tools.

### A.39.2 Options as arguments or environment?

We distinguish two views on options. One is to see them as additional parameters that require strict existence, type and domain-checking and the other is to consider them‘locally scoped environment variables’. Most systems adopt the first option. SWI-Prolog adopts the second: it silently ignores options that are not supported but does type and domain checking of option-values. The‘environment’view is commonly used in applications to create predicates supporting more options using the skeleton below. This way of programming requires that *pred1* and *pred2* do not interpret the same option differently. In cases where this is not true, the options must be distributed by *some_pred*. We have been using this programming style for many years and in practice it turns out that the need for active distribution of options is rare. I.e., options either have distinct names or multiple predicates implement the same option but this has the desired effect. An example of the latter is the `encoding` option, which typically needs to be applied consistently.

``` code
some_pred(..., Options) :-
      pred1(..., Options),
      pred2(..., Options).
```

As stated before, options provide a readable alternative to high-arity predicates and offer a robust mechanism to evolve the API, but at the cost of some runtime overhead and weaker consistency checking, both at compiletime and runtime. From our experience, the‘environment’approach is productive, but the consequence is that mistyped options are silently ignored. The option infrastructure described in this section tries to remedy these problems.

### A.39.3 Improving on the current situation

Whether we see options as arguments or locally scoped environment variables, the most obvious way to improve on the current situation is to provide reflective support for options: discover that an argument is an option-list and find what options are supported. Reflective access to options can be used by the compiler and development environment as well as by the runtime system to warn or throw errors.

#### A.39.3.1 Options as types

An obvious approach to deal with options is to define the different possible option values as a type and type the argument that processes the option as list(\<option_type\>), as illustrated below. Considering options as types fully covers the case where we consider options as additional parameters.

``` code
:- type open_option ---> type(stream_type) |
                         alias(atom) | ... .
:- pred open(source_sink, open_mode, stream, list(open_option)).
```

There are three reasons for considering a different approach:

- There is no consensus about types in the Prolog world, neither about what types should look like, nor whether or not they are desirable. It is not likely that this debate will be resolved shortly.
- Considering options as types does not support the‘environment’view, which we consider the most productive.
- Even when using types, we need reflective access to what options are provided in order to be able to write compile or runtime conditional code.

#### A.39.3.2 Reflective access to options

From the above, we conclude that we require reflective access to find out whether an option is supported and valid for a particular predicate. Possible option values must be described by types. Due to lack of a type system, we use `library(error)` to describe allowed option values. Predicate options are declared using [predicate_options/3](predicate_options.html#predicate_options/3):

\[det\]**predicate_options**(`:PI, +Arg, +Options`)  
Declare that the predicate `PI` processes options on `Arg`. `Options` is a list of options processed. Each element is one of:

- Option(ModeAndType) `PI` processes Option. The option-value must comply to ModeAndType. Mode is one of + or - and Type is a type as accepted by [must_be/2](error.html#must_be/2).
- pass_to(:`PI`,`Arg`) The option-list is passed to the indicated predicate.

Below is an example that processes the option `header(boolean)` and passes all options to [open/4](IO.html#open/4):

``` code
:- predicate_options(write_xml_file/3, 3,
                     [ header(boolean),
                       pass_to(open/4, 4)
                     ]).

write_xml_file(File, XMLTerm, Options) :-
    open(File, write, Out, Options),
    (   option(header(true), Options, true)
    ->  write_xml_header(Out)
    ;   true
    ),
    ...
```

This predicate may only be used as a *directive* and is processed by [expand_term/2](consulting.html#expand_term/2). Option processing can be specified at runtime using assert_predicate_options/3, which is intended to support program analysis.

\[semidet\]**assert_predicate_options**(`:PI, +Arg, +Options, ?New`)  
As predicate_options(:`PI`, +`Arg`, +`Options`). `New` is a boolean indicating whether the declarations have changed. If `New` is provided and `false`, the predicate becomes semidet and fails without modifications if modifications are required.

The predicates below realise the support for compile and runtime checking for supported options.

\[nondet\]**current_predicate_option**(`:PI, ?Arg, ?Option`)  
True when `Arg` of `PI` processes `Option`. For example, the following is true:

``` code
?- current_predicate_option(open/4, 4, type(text)).
true.
```

This predicate is intended to support conditional compilation using [if/1](consulting.html#if/1) ... [endif/0](consulting.html#endif/0). The predicate [current_predicate_options/3](predicate_options.html#current_predicate_options/3) can be used to access the full capabilities of a predicate.

\[det\]**check_predicate_option**(`:PI, +Arg, +Option`)  
Verify predicate options at runtime. Similar to [current_predicate_option/3](predicate_options.html#current_predicate_option/3), but intended to support runtime checking.

Errors  
\- `existence_error(option, OptionName)` if the option is not supported by `PI`.  
- `type_error(Type, Value)` if the option is supported but the value does not match the option type. See [must_be/2](error.html#must_be/2).

The predicates below can be used in a development environment to inform the user about supported options. PceEmacs uses this for colouring option names and values.

\[nondet\]**current_option_arg**(`:PI, ?Arg`)  
True when `Arg` of `PI` processes predicate options. Which options are processed can be accessed using [current_predicate_option/3](predicate_options.html#current_predicate_option/3).

\[nondet\]**current_predicate_options**(`:PI, ?Arg, ?Options`)  
True when `Options` is the current active option declaration for `PI` on `Arg`. See [predicate_options/3](predicate_options.html#predicate_options/3) for the argument descriptions. If `PI` is ground and refers to an undefined predicate, the autoloader is used to obtain a definition of the predicate.

The library can execute a complete check of your program using [check_predicate_options/0](predicate_options.html#check_predicate_options/0):

\[det\]**check_predicate_options**  
Analyse loaded program for erroneous options. This predicate decompiles the current program and searches for calls to predicates that process options. For each option list, it validates whether the provided options are supported and validates the argument type. This predicate performs partial dataflow analysis to track option-lists inside a clause.

See also  
[derive_predicate_options/0](predicate_options.html#derive_predicate_options/0) can be used to derive declarations for predicates that pass options. This predicate should normally be called before [check_predicate_options/0](predicate_options.html#check_predicate_options/0).

The library offers predicates that may be used to create declarations for your application. These predicates are designed to cooperate with the module system.

\[det\]**derive_predicate_options**  
Derive new predicate option declarations. This predicate analyses the loaded program to find clauses that process options using one of the predicates from `library(option)` or passes options to other predicates that are known to process options. The process is repeated until no new declarations are retrieved.

See also  
autoload/0 may be used to complete the loaded program.

\[det\]**retractall_predicate_options**  
Remove all dynamically (derived) predicate options.

\[nondet\]**derived_predicate_options**(`:PI, ?Arg, ?Options`)  
Derive option arguments using static analysis. True when `Options` is the current *derived* active option declaration for `PI` on `Arg`.

\[det\]**derived_predicate_options**(`+Module`)  
Derive predicate option declarations for a module. The derived options are printed to the `current_output` stream.
