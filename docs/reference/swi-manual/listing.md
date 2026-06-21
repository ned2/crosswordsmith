
## A.24 library(listing): List programs and pretty print clauses

To be done  
\- More settings, support *Coding Guidelines for Prolog* and make the suggestions there the default.  
- Provide persistent user customization

This module implements listing code from the internal representation in a human readable format.

- [listing/0](listing.html#listing/0) lists a module.
- [listing/1](listing.html#listing/1) lists a predicate or matching clause
- [listing/2](listing.html#listing/2) lists a predicate or matching clause with options
- [portray_clause/2](listing.html#portray_clause/2) pretty-prints a clause-term

Layout can be customized using `library(settings)`. The effective settings can be listed using [list_settings/1](settings.html#list_settings/1) as illustrated below. Settings can be changed using [set_setting/2](settings.html#set_setting/2).

``` code
?- list_settings(listing).
========================================================================
Name                      Value (*=modified) Comment
========================================================================
listing:body_indentation  4              Indentation used goals in the body
listing:tab_distance      0              Distance between tab-stops.
...
```

**listing**  
Lists all predicates defined in the calling module. Imported predicates are not listed. To list the content of the module `mymodule`, use one of the calls below.

``` code
?- mymodule:listing.
?- listing(mymodule:_).
```

\[det\]**listing**(`:What`)  
\[det\]**listing**(`:What, +Options`)  
List matching clauses. `What` is either a plain specification or a list of specifications. Plain specifications are:

- Predicate indicator (Name/Arity or Name`//`Arity) Lists the indicated predicate. This also outputs relevant *declarations*, such as [multifile/1](dynamic.html#multifile/1) or [dynamic/1](dynamic.html#dynamic/1).

- A *Head* term. In this case, only clauses whose head unify with *Head* are listed. This is illustrated in the query below that only lists the first clause of [append/3](lists.html#append/3).

  ``` code
  ?- listing(append([], _, _)).
  lists:append([], L, L).
  ```

- A clause reference as obtained for example from [nth_clause/3](examineprog.html#nth_clause/3).

The following options are defined:

**variable_names**(`+How`)  
One of `source` (default) or `generated`. If `source`, for each clause that is associated to a source location the system tries to restore the original variable names. This may fail if macro expansion is not reversible or the term cannot be read due to different operator declarations. In that case variable names are generated.

**source**(`+Bool`)  
If `true` (default `false`), extract the lines from the source files that produced the clauses, i.e., list the original source text rather than the *decompiled* clauses. Each set of contiguous clauses is preceded by a comment that indicates the file and line of origin. Clauses that cannot be related to source code are decompiled where the comment indicates the decompiled state. This is notably practical for collecting the state of *multifile* predicates. For example:

``` code
?- listing(file_search_path, [source(true)]).
```

**thread**(`+ThreadId`)  
If a predicate is *thread local*, list the clauses as seen by the given `ThreadId`. Ignored if the predicate is not thread local.

\[det\]**portray_clause**(`+Clause`)  
\[det\]**portray_clause**(`+Out:stream, +Clause`)  
\[det\]**portray_clause**(`+Out:stream, +Clause, +Options`)  
Portray‘`Clause`’on the current output stream. Layout of the clause is to our best standards. Deals with control structures and calls via meta-call predicates as determined using the predicate property meta_predicate. If `Clause` contains attributed variables, these are treated as normal variables.

Variable names are by default generated using [numbervars/4](manipterm.html#numbervars/4) using the option `singletons(true)`. This names the variables `A`, `B`, ... and the singletons `_`. Variables can be named explicitly by binding them to a term `'$VAR'(Name)`, where `Name` is an atom denoting a valid variable name (see the option `numbervars(true)` from [write_term/2](termrw.html#write_term/2)) as well as by using the `variable_names(Bindings)` option from [write_term/2](termrw.html#write_term/2).

`Options` processed in addition to [write_term/2](termrw.html#write_term/2) options:

**variable_names**(`+Bindings`)  
See above and [write_term/2](termrw.html#write_term/2).

**indent**(`+Columns`)  
Left margin used for the clause. Default `0`.

**module**(`+Module`)  
`Module` used to determine whether a goal resolves to a meta predicate. Default `user`.
