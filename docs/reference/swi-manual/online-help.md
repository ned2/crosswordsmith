
## 2.7 Online Help

### 2.7.1 library(help): Text based manual

This module provides [help/1](online-help.html#help/1) and [apropos/1](online-help.html#apropos/1) that give help on a topic or searches the manual for relevant topics.

By default the result of [help/1](online-help.html#help/1) is sent through a *pager* such as `less`. This behaviour is controlled by the following:

- The Prolog flag **help_pager**, which can be set to one of the following values:
  **false**  
  Never use a pager.

  **default**  
  Use default behaviour. This tries to determine whether Prolog is running interactively in an environment that allows for a pager. If so it examines the environment variable `PAGER` or otherwise tries to find the `less` program.

  **`Callable`**  
  A `Callable` term is interpreted as `program_name(Arg, ...)`. For example, `less('-r')` would be the default. Note that the program name can be an absolute path if single quotes are used.

\[det\]**help**  
\[det\]**help**(`+What`)  
Show help for `What`. `What` is a term that describes the `topics(s)` to give help for. Notations for `What` are:

**`Atom`**  
This ambiguous form is most commonly used and shows all matching documents. For example:

``` code
?- help(append).
```

`Name` **`/`** `Arity`  
Give help on predicates with matching `Name`/`Arity`. `Arity` may be unbound.

`Name` **`//`** `Arity`  
Give help on the matching DCG rule (non-terminal)

`Module`**`:`**`Name`  
Give help on predicates with `Name` in `Module` and any arity. Used for loaded code only.

`Module`**`:`**`Name``/``Arity`  
Give help on predicates with `Name` in `Module` and `Arity`. Used for loaded code only.

**f**(`Name/Arity`)  
Give help on the matching Prolog arithmetic functions.

**c**(`Name`)  
Give help on the matching C interface function

**section**(`Label`)  
Show the section from the manual with matching `Label`.

[help/1](online-help.html#help/1) shows documentation from the manual as well as from loaded user code if the code is documented using PlDoc. To show only the documentatoion of the loaded predicate we may prefix predicate indicator with the module in which it is defined.

If an exact match fails this predicates attempts fuzzy matching and, when successful, display the results headed by a warning that the matches are based on fuzzy matching.

If possible, the results are sent through a *pager* such as the `less` program. This behaviour is controlled by the Prolog flag `help_pager`. See section level documentation.

See also  
[apropos/1](online-help.html#apropos/1) for searching the manual names and summaries.

\[semidet,multifile\]**show_html_hook**(`+HTML:string`)  
Hook called to display the extracted `HTML` document. If this hook fails the `HTML` is rendered to the console as plain text using html_text/2.

\[det\]**apropos**(`+Query`)  
Print objects from the manual whose name or summary match with `Query`. `Query` takes one of the following forms:

`Type`**`:`**`Text`  
Find objects matching `Text` and filter the results by `Type`. `Type` matching is a case intensitive *prefix* match. Defined types are `section`, `cfunction`, `function`, `iso_predicate`, `swi_builtin_predicate`, `library_predicate`, `dcg` and aliases `chapter`, `arithmetic`, `c_function`, `predicate`, `nonterminal` and `non_terminal`. For example:

``` code
?- apropos(c:close).
?- apropos(f:min).
```

**`Text`**  
`Text` is broken into tokens. A topic matches if all tokens appear in the name or summary of the topic. Matching is case insensitive. Results are ordered depending on the quality of the match.

\[nondet\]**help_apropos**(`+Query, -Obj, -Summary, -Score`)  
Find matching documented objects in the help database. `Obj` is the formal object identifier, `Summary` its summary description and `Score` is a number indicating the quality of the match.

\[semidet\]**help_text**(`+Predicate:term, -HelpText:string`)  
When `Predicate` is a term of the form `Name/Arity` for which documentation exists, `HelpText` is the documentation in textual format (parsed from the HTML help).

### 2.7.2 library(explain): Describe Prolog Terms

The `library(explain)` describes prolog-terms. The most useful functionality is its cross-referencing function.

``` code
?- explain(subset(_,_)).
"subset(_, _)" is a compound term
    from 2-th clause of lists:subset/2
    Referenced from 46-th clause of prolog_xref:imported/3
    Referenced from 68-th clause of prolog_xref:imported/3
lists:subset/2 is a predicate defined in
    /staff/jan/lib/pl-5.6.17/library/lists.pl:307
    Referenced from 2-th clause of lists:subset/2
    Possibly referenced from 2-th clause of lists:subset/2
```

Note that PceEmacs can jump to definitions and [gxref/0](xref.html#gxref/0) can be used for an overview of dependencies.

\[det\]**explain**(`@Term`)  
Give an explanation on `Term`. `Term` can be any Prolog data object. Some terms have a specific meaning:

- A (partial) reference to a predicate gives the predicates, its main properties and references to the predicates. Partial references are:
  - Module:Name/Arity
  - Module:Head
  - Name/Arity
  - Name`//`Arity
  - Name
  - Module:Name
- Some predicate properties. This lists predicates as above the have this property. The specification can be of the shape `Module:Property` or just `Property`. The qualified version limits the result to predicates defined in Module. Supported properties are:
  - dynamic
  - thread_local
  - multifile
  - tabled

\[nondet\]**explain**(`@Term, -Explanation`)  
True when `Explanation` is an explanation of `Term`. The explaination is a list of elements that is printed using `print_message(information, explain(Explanation))`.
