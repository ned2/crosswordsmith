
## A.32 library(option): Option list processing

See also  
\- `library(record)`  
- Option processing capabilities may be declared using the directive [predicate_options/3](predicate_options.html#predicate_options/3).

The `library(option)` provides some utilities for processing option lists. Option lists are commonly used as an alternative for many arguments. Examples of built-in predicates are [open/4](IO.html#open/4) and [write_term/3](termrw.html#write_term/3). Naming the arguments results in more readable code, and the list nature makes it easy to extend the list of options accepted by a predicate. Option lists come in two styles, both of which are handled by this library.

- Name(Value)  
  This is the preferred style.
- Name = Value  
  This is often used, but deprecated.

SWI-Prolog *dicts* provide a convenient and efficient alternative to option lists. For this reason, both built-in predicates and predicates that use this library support dicts transparantly.

Processing option lists inside time-critical code (loops) can cause serious overhead. The above mentioned *dicts* is the preferred mitigation. A more portable alternative is to define a record using `library(record)` and initialise this using make\_`<`record[\>/2](arith.html#%3E/2). In addition to providing good performance, this also provides type-checking and central declaration of defaults.

Options typically have exactly one argument. The library does support options with 0 or more than one argument with the following restrictions:

- The predicate [option/3](option.html#option/3) and [select_option/4](option.html#select_option/4), involving default are meaningless. They perform an `arg(1, Option, Default)`, causing failure without arguments and filling only the first option-argument otherwise.
- [meta_options/3](option.html#meta_options/3) can only qualify options with exactly one argument.

\[semidet\]**option**(`?Option, +Options`)  
Get an `Option` from `Options`. Fails silently if the option does not appear in `Options`. If `Option` appears multiple times in `Options`, the first value is used.

|           |                                                     |
|-----------|-----------------------------------------------------|
| `Option`  | Term of the form Name(?Value).                      |
| `Options` | is a list of Name(Value) or `Name=Value` or a dict. |

\[det\]**option**(`?Option, +Options, +Default`)  
Get an `Option` from `Options`. If `Option` does not appear in `Options`, unify the value with `Default`. If `Option` appears multiple times in `Options`, the first value is used. For example

``` code
?- option(max_depth(D), [x(a), max_depth(20)], 10).
D = 20.
?- option(max_depth(D), [x(a)], 10).
D = 10.
```

|           |                                                     |
|-----------|-----------------------------------------------------|
| `Option`  | Term of the form Name(?Value).                      |
| `Options` | is a list of Name(Value) or `Name=Value` or a dict. |

\[semidet\]**select_option**(`?Option, +Options, -RestOptions`)  
Get and remove `Option` from `Options`. As [option/2](option.html#option/2), removing the matching option from `Options` and unifying the remaining options with `RestOptions`. If `Option` appears multiple times in `Options`, the first value is used. Note that if `Options` contains multiple terms that are compatible to `Option`, the first is used to set the value of `Option` and the duplicate appear in `RestOptions`.

\[det\]**select_option**(`?Option, +Options, -RestOptions, +Default`)  
Get and remove `Option` with default value. As [select_option/3](option.html#select_option/3), but if `Option` is not in `Options`, its value is unified with `Default` and `RestOptions` with `Options`.

\[det\]**merge_options**(`+New, +Old, -Merged`)  
Merge two option sets. If `Old` is a dict, `Merged` is a dict. Otherwise `Merged` is a sorted list of options using the canonical format Name(Value) holding all options from `New` and `Old`, after removing conflicting options from `Old`.

Multi-values options (e.g., `proxy(Host, Port)`) are allowed, where both option-name and arity define the identity of the option.

\[det\]**meta_options**(`+IsMeta, :Options0, -Options`)  
Perform meta-expansion on options that are module-sensitive. Whether an option name is module-sensitive is determined by calling `call(IsMeta, Name)`. Here is an example:

``` code
    meta_options(is_meta, OptionsIn, Options),
    ...

is_meta(callback).
```

Meta-options must have exactly one argument. This argument will be qualified.

To be done  
Should be integrated with declarations from [predicate_options/3](predicate_options.html#predicate_options/3).

\[det\]**dict_options**(`?Dict, ?Options`)  
Convert between an option list and a dictionary. One of the arguments must be instantiated. If the option list is created, it is created in canonical form, i.e., using Option(Value) with the `Options` sorted in the standard order of terms. Note that the conversion is not always possible due to different constraints and conversion may thus lead to (type) errors.

- `Dict` keys can be integers. This is not allowed in canonical option lists.
- `Options` can hold multiple options with the same key. This is not allowed in dicts. This predicate removes all but the first option on the same key.
- `Options` can have more than one value (`name(V1,V2)`). This is not allowed in dicts.

Also note that most system predicates and predicates using this library for processing the option argument can both work with classical Prolog options and dicts objects.
