
## 4.20 Term reading and writing

This section describes the basic term reading and writing predicates. The predicates [format/\[1,2\]](format.html#format/1) and [writef/2](writef.html#writef/2) provide formatted output. Writing to Prolog data structures such as atoms or code-lists is supported by [with_output_to/2](IO.html#with_output_to/2) and [format/3](format.html#format/3).

Reading is sensitive to the Prolog flag [character_escapes](flags.html#flag:character_escapes), which controls the interpretation of the `\` character in quoted atoms and strings.

\[ISO\]**write_term**(`+Term, +Options`)  
The predicate [write_term/2](termrw.html#write_term/2) is the generic form of all Prolog term-write predicates. Valid options are:

**attributes**(`Atom`)  
Define how attributed variables (see [section 8.1](attvar.html#sec:8.1)) are written. The default is determined by the Prolog flag [write_attributes](flags.html#flag:write_attributes). Defined values are `ignore` (ignore the attribute), `dots` (write the attributes as `{...}`), `write` (simply hand the attributes recursively to [write_term/2](termrw.html#write_term/2)) and `portray` (hand the attributes to [attr_portray_hook/2](attvar.html#attr_portray_hook/2)).

**back_quotes**(`Atom`)  
Fulfills the same role as the [back_quotes](flags.html#flag:back_quotes) prolog flag. Notably, the value `string` causes string objects to be printed between back quotes and `symbol_char` causes the backquote to be printed unquoted. In all other cases the backquote is printed as a quoted atom. The default is derived from the Prolog flag using the value associated with the `module(Module)` option or the `user` module.

**brace_terms**(`Bool`)  
If `true` (default), write `{}(X)` as `{X}`. See also `dotlists` and `ignore_ops`.

**blobs**(`Atom`)  
Define how non-text blobs are handled. By default, this is left to the write handler specified with the blob type. Using `portray`, [portray/1](termrw.html#portray/1) is called for each blob encountered. See [section 12.4.10](foreigninclude.html#sec:12.4.10).

**character_escapes**(`Bool`)  
If `true` and `quoted(true)` is active, special characters in quoted atoms and strings are emitted as ISO escape sequences. Default is taken from the reference module (see below).

**character_escapes_unicode**(`Bool`)  
If `true` and `character_escapes(true)` and `quoted(true)` are active escaped characters are written using `\uXXXX` or `\UXXXXXXXX` syntax. The default depends on the Prolog flag [character_escapes_unicode](flags.html#flag:character_escapes_unicode)

**cycles**(`Bool`)  
If `true` (default), cyclic terms are written as `@(Template, Substitutions)`, where `Substitutions` is a list `Var` = `Value`. If `cycles` is `false`, `max_depth` is not given, and `Term` is cyclic, [write_term/2](termrw.html#write_term/2) raises a `domain_error`.^(106The cycles option and the cyclic term representation using the @-term are copied from SICStus Prolog. However, the default in SICStus is set to `false` and SICStus writes an infinite term if not protected by, e.g., the `depth_limit` option.) See also the `cycles` option in [read_term/2](termrw.html#read_term/2).

**dotlists**(`Bool`)  
If `true` (default `false`), write lists using the dotted term notation rather than the list notation.^(107Copied from ECLiPSe.) Note that as of version 7, the list constructor is `'[|]'`. Using `dotlists(true)`, [write_term/2](termrw.html#write_term/2) writes a list using‘.’as constructor. This is intended for communication with programs such as other Prolog systems, that rely on this notation. See also the option `no_lists(true)` to use the actual SWI-Prolog list functor.

**fullstop**(`Bool`)  
If `true` (default `false`), add a fullstop token to the output. The dot is preceded by a space if needed and followed by a space (default) or newline if the `nl(true)` option is also given.^(108Compatible with [ECLiPSe](http://eclipseclp.org/doc/bips/kernel/ioterm/write_term-3.html))

**float_format**(`+Atom`)  
Print floating point numbers using [format/2](format.html#format/2) as `format(Atom, [Float])`. The default is `~h`. This option is compatible with SICStus. See [format/2](format.html#format/2) for valid format specifiers and the `integer_format` option for additional comments.

**ignore_ops**(`Bool`)  
If `true`, the generic term representation (\<`functor`\>(\<`args`\> ... )) will be used for all terms. Otherwise (default), operators will be used where appropriate. The ISO standard also dictates lists to be written as `.(Head,Tail)`. SWI-Prolog writes lists using Prolog list notation unless the option `dotlists(true)` is also given. This is done because SWI-Prolog lists do not use the dot as list cell functor. PIP-0105 added the option `portable(Bool)` to accomplish unambiguous transfer of terms while maximizing readability.

**integer_format**(`+Atom`)  
Print integers using [format/2](format.html#format/2) as `format(Atom, [Int])`. The default is `~d`. This allows to print integers using an alternative *radix*, using e.g. `~16r` or `0x~16r` or to use digit grouping using e.g. `~D`. Note that the user is responsible to provide a format that produces valid Prolog syntax if the term must be readable by Prolog. The format must accept exactly one argument. If that is not satisfied, printing an integer results in an exception. See [format/2](format.html#format/2) for for valid format specifiers.

**max_depth**(`Integer`)  
If the term is nested deeper than `Integer`, print the remainder as ellipses ( ... ). A 0 (zero) value (default) imposes no depth limit. This option also delimits the number of printed items in a list. Example:

``` code
?- write_term(a(s(s(s(s(0)))), [a,b,c,d,e,f]),
              [max_depth(3)]).
a(s(s(...)), [a, b|...])
true.
```

Used by the top level and debugger to limit screen output. See also the options `max_text(Length)` and `truncated(-Bool)` as well as the prolog flags [answer_write_options](flags.html#flag:answer_write_options) and [debugger_write_options](flags.html#flag:debugger_write_options).

**max_text**(`Length`)  
Abbreviate atoms and strings that are longer than `Length` characters. The abbreviated characters are indicated using elipsis (`...`), if possible using the Unicode character U-2026 and using `...` otherwise. On sufficiantly long texts the elipsis is followed by the three last characters. If `Length` is -1, no length limit is imposed. If `Length` is 0, it only prints the elipsis. When `quoted(true)` is active, the quotes are not included in the length restriction. The text between the quotes represents at max `Length` characters, but the lexical representation can be longer due to required escape sequences. This option is defined by PIP-0105.^(109The current draft calls this option `text_max(Max)`. Considering all other options and flags use max\_\*, SWI-Prolog uses `max_text`. This may change depending on how PIP-0105 evolves.)

**module**(`Module`)  
Define the reference module (default `user`). This defines the default value for the [character_escapes](flags.html#flag:character_escapes) option as well as the operator definitions to use. If `Module` does not exist it is *not* created and the `user` module is used. See also [op/3](operators.html#op/3) and [read_term/2](termrw.html#read_term/2), providing the same option.

Note that, currently, [write_term/2](termrw.html#write_term/2) is not a *meta predicate* and the default module is thus not derived from the calling context. This is likely to change in the future.

**nl**(`Bool`)  
Add a newline to the output. See also the `fullstop` option.

**no_lists**(`Bool`)  
Do not use list notation. This is similar to `dotlists(true)`, but uses the SWI-Prolog list functor, which is by default `'[|]'` instead of the ISO Prolog `'.'`. Used by display/1.

**numbervars**(`Bool`)  
If `true`, terms of the format `$VAR(N)`, where `N` is an integer that fits in 64-bit,^(110Larger integers are ignored. As no term that fits into memory can have that many variables, this is not a restriction.) will be written as a variable name. For `N` in 0..25 it emits A..Z. For higher numbers it emits An..Zn, where `n` is `N`//26. For negative numbers it emits S\_`N`, which is used for representing shared sub-terms and cyclic terms.

If `N` is an atom that is syntactically a valid variable it is written without quotes. This extension allows for writing variables with user-provided names.

The default for this flag is `false` unless the `portrayed` option is enabled. See also [numbervars/3](manipterm.html#numbervars/3) and the option `variable_names`.

**partial**(`Bool`)  
If `true` (default `false`), do not reset the logic that inserts extra spaces that separate tokens where needed. This is intended to solve the problems with the code below. Calling `write_value(``.``)` writes `..`, which cannot be read. By adding `partial(true)` to the option list, it correctly emits `. .`. Similar problems appear when emitting operators using multiple calls to [write_term/3](termrw.html#write_term/3).

``` code
write_value(Value) :-
        write_term(Value, [partial(true)]),
        write('.'), nl.
```

In addition, if the priority is not 1200 or 999 this assumes we are printing an operand of an operator. If `Term` is an atom that is also an operator it will always be embraced.^(111If the priority is 1200 it is assumed to be a toplevel term and if the priority is 999 it is assumed to be a list element or argument of a compound term.)

**portable**(`Bool`)  
Similar to `ignore_ops(Bool)`, but when `true`, lists, *brace terms* and the infix operator “`,`” are printed normally. The name *portable* refers to the fact that any term written this way can be read unambigously regardless of differences in the (operator) context while providing maximal readability and avoiding high stack usage implied by the ISO `ignore_ops` list representation as `.(a,.(b,[]))`.^(112The dotted list notation cannot be reduced to list cells before we reach the closing parenthesis.) This option is defined by PIP-0105.

For example:

``` code
?- write_term((a,b+c,[a],{d}), [portable]).
a,+(b,c),[a],{d}
```

**portray**(`Bool`)  
Same as `portrayed(Bool)`. Deprecated.

**portray_goal**(`:Goal`)  
Implies `portray(true)`, but calls `Goal` rather than the predefined hook [portray/1](termrw.html#portray/1). `Goal` is called through call/3, where the first argument is `Goal`, the second is the term to be printed and the 3rd argument is the current write option list. The write option list is copied from the write_term call, but the list is guaranteed to hold an option `priority` that reflects the current priority.

**portrayed**(`Bool`)  
If `true`, the hook [portray/1](termrw.html#portray/1) is called before printing a term that is not a variable. If [portray/1](termrw.html#portray/1) succeeds, the term is considered printed. See also [print/1](termrw.html#print/1). The default is `false`. This option is an extension to the ISO write_term options. If this option is set, the `numbervars` option *defaults* to `true`.

**priority**(`Integer`)  
An integer between 0 and 1200 representing the‘context priority’. Default is 1200. Can be used to write partial terms appearing as the argument to an operator. For example:

``` code
        format('~w = ', [VarName]),
        write_term(Value, [quoted(true), priority(699)])
```

**quoted**(`Bool`)  
If `true`, atoms and strings that need quotes will be quoted. The default is `false`. If [character_escapes](flags.html#flag:character_escapes) is `true` (default) characters in the quoted atom or string are escaped using backslash (`\`) sequences. To the minimum, the quote itself, newlines and backslash characters are escaped to make the output valid for [read/1](termrw.html#read/1). All unassigned unicode characters and characters in the Unicode *separator* (Z\*) and *control* (C\*) classes except for the ASCII space (`\u0020`) are escaped. For those characters for which an ISO Prolog single character escape, e.g., `\t` is defined, this is used. Otherwise the output depends on the option `character_escapes_unicode`. If this flag applies(default) the widely accepted `\uXXXX` or `\UXXXXXXXX` is used. Otherwise the ISO Prolog `\x<hex>\` syntax is used.

**quote_non_ascii**(`Bool`)  
Quote an atom that contains non-ASCII, i.e., larger than 127 code points. The Prolog standard only describes non-quoted atom syntax containing ASCII characters. While SWI-Prolog extends this to Unicode (see [section 2.15.1.9](syntax.html#sec:2.15.1.9)), transferring atoms holding non-ASCII text to other Prolog implementations may cause problems. This flag is used by [write_canonical/1](termrw.html#write_canonical/1).

**spacing**(`+Spacing`)  
Determines whether and where extra white space is added to enhance readability. The default is `standard`, adding only space where needed for proper tokenization by [read_term/3](termrw.html#read_term/3). Currently, the only other value is `next_argument`, adding a space after a comma used to separate arguments in a term or list.

**truncated**(`-Bool`)  
If `max_depth(+Depth)` and/or `max_text(+Length)` is active and the limit has been exceeded, `Bool` is unified to `true`. Else it is unified to `false`. This option is part of PIP-0105.

**variable_names**(`+List`)  
Assign names to variables in `Term`. `List` is a list of terms `Name` = `Var`, where `Name` is an atom that represents a valid Prolog variable name. Terms where `Var` is bound or is a variable that does not appear in `Term` are ignored. Raises an error if `List` is not a list, one of the members is not a term `Name` = `Var`, `Name` is not an atom or `Name` does not represent a valid Prolog variable name.

The implementation binds the variables from `List` to a term `'$VAR'`(`Name`). Like [write_canonical/1](termrw.html#write_canonical/1), terms that where already bound to `'$VAR'`(`X`) before [write_term/2](termrw.html#write_term/2) are printed normally, unless the option `numbervars(true)` is also provided. If the option `numbervars(true)` is used, the user is responsible for avoiding collisions between assigned names and numbered names. See also the `variable_names` option of [read_term/2](termrw.html#read_term/2).

Possible variable attributes (see [section 8.1](attvar.html#sec:8.1)) are ignored. In most cases one should use [copy_term/3](attvar.html#copy_term/3) to obtain a copy that is free of attributed variables and handle the associated constraints as appropriate for the use-case.

\[ISO\]**write_term**(`+Stream, +Term, +Options`)  
As [write_term/2](termrw.html#write_term/2), but output is sent to `Stream` rather than the current output.

\[semidet\]**write_length**(`+Term, -Length, +Options`)  
True when `Length` is the number of characters emitted for `write_term(Term, Options)`. In addition to valid options for [write_term/2](termrw.html#write_term/2), it processes the option:

**max_length**(`+MaxLength`)  
If provided, fail if `Length` would be larger than `MaxLength`. The implementation ensures that the runtime is limited when computing the length of a huge term with a bounded maximum.

\[ISO\]**write_canonical**(`+Term`)  
Write `Term` on the current output stream using standard parenthesised prefix notation (i.e., ignoring operator declarations). Atoms that need quotes are quoted. Terms written with this predicate can always be read back, regardless of current operator declarations. Equivalent to [write_term/2](termrw.html#write_term/2) using the options `ignore_ops`, `quoted`, `quote_non_ascii`, `brace_terms(false)` and `numbervars` after [numbervars/4](manipterm.html#numbervars/4) using the `singletons` option.

Note that due to the use of [numbervars/4](manipterm.html#numbervars/4), non-ground terms must be written using a *single* [write_canonical/1](termrw.html#write_canonical/1) call. This used to be the case anyhow, as garbage collection between multiple calls to one of the write predicates can change the `_`\<`NNN`\> identity of the variables.

\[ISO\]**write_canonical**(`+Stream, +Term`)  
Write `Term` in canonical form on `Stream`.

\[ISO\]**write**(`+Term`)  
Write `Term` to the current output, using brackets and operators where appropriate.

\[ISO\]**write**(`+Stream, +Term`)  
Write `Term` to `Stream`.

\[ISO\]**writeq**(`+Term`)  
Write `Term` to the current output, using brackets and operators where appropriate. Atoms that need quotes are quoted. Terms written with this predicate can be read back with [read/1](termrw.html#read/1) provided the currently active operator declarations are identical and Term. Equivalent to `write_term(Term, [quoted(true), numbervars(true)])`.

\[ISO\]**writeq**(`+Stream, +Term`)  
Write `Term` to `Stream`, inserting quotes.

**writeln**(`+Term`)  
Equivalent to `write(Term), nl.`. The output stream is locked, which implies no output from other threads can appear between the term and newline.

**writeln**(`+Stream, +Term`)  
Equivalent to `write(Stream, Term), nl(Stream).`. The output stream is locked, which implies no output from other threads can appear between the term and newline.

**print**(`+Term`)  
Print a term for debugging purposes. The predicate [print/1](termrw.html#print/1) acts as if defined as below.

``` code
print(Term) :-
    current_prolog_flag(print_write_options, Options), !,
    write_term(Term, Options).
print(Term) :-
    write_term(Term, [ portray(true),
                       numbervars(true),
                       quoted(true)
                     ]).
```

The [print/1](termrw.html#print/1) predicate is used primarily through the `~p` escape sequence of [format/2](format.html#format/2), which is commonly used in the recipes used by [print_message/2](printmsg.html#print_message/2) to emit messages.

The classical definition of this predicate is equivalent to the ISO predicate [write_term/2](termrw.html#write_term/2) using the options `portray(true)` and `numbervars(true)`. The `portray(true)` option allows the user to implement application-specific printing of terms printed during debugging to facilitate easy understanding of the output. See also [portray/1](termrw.html#portray/1) and `library(portray_text)`. SWI-Prolog adds `quoted(true)` to (1) facilitate the copying/pasting of terms that are not affected by [portray/1](termrw.html#portray/1) and to (2) allow numbers, atoms and strings to be more easily distinguished, e.g., `42`, `'42'` and `"42"`.

**print**(`+Stream, +Term`)  
Print `Term` to `Stream`.

**portray**(`+Term`)  
A dynamic predicate, which can be defined by the user to change the behaviour of [print/1](termrw.html#print/1) on (sub)terms. For each subterm encountered that is not a variable [print/1](termrw.html#print/1) first calls [portray/1](termrw.html#portray/1) using the term as argument. For lists, only the list as a whole is given to [portray/1](termrw.html#portray/1). If [portray/1](termrw.html#portray/1) succeeds [print/1](termrw.html#print/1) assumes the term has been written.

\[ISO\]**read**(`-Term`)  
Read the next **Prolog term** from the current input stream and unify it with `Term`. On reaching end-of-file `Term` is unified with the atom `end_of_file`. This is the same as [read_term/2](termrw.html#read_term/2) using an empty option list.

**\[NOTE\]** You might have found this while looking for a predicate to read input from a file or the user. Quite likely this is not what you need in this case. This predicate is for reading a **Prolog term** which may span multiple lines and must end in a *full stop* (dot character followed by a layout character). The predicates for reading and writing Prolog terms are particularly useful for storing Prolog data in a file or transferring them over a network communication channel (socket) to another Prolog process. The libraries provide a wealth of predicates to read data in other formats. See e.g., `library(readutil)`, `library(pure_input)` or libraries from the extension packages to read XML, JSON, YAML, etc.

\[ISO\]**read**(`+Stream, -Term`)  
Read the next **Prolog term** from `Stream`. See [read/1](termrw.html#read/1) and [read_term/2](termrw.html#read_term/2) for details.

**read_clause**(`+Stream, -Term, +Options`)  
Equivalent to [read_term/3](termrw.html#read_term/3), but sets options according to the current compilation context and optionally processes comments. Defined options:

**syntax_errors**(`+Atom`)  
See [read_term/3](termrw.html#read_term/3), but the default is `dec10` (report and restart).

**term_position**(`-TermPos`)  
Same as for [read_term/3](termrw.html#read_term/3).

**subterm_positions**(`-TermPos`)  
Same as for [read_term/3](termrw.html#read_term/3).

**variable_names**(`-Bindings`)  
Same as for [read_term/3](termrw.html#read_term/3).

**process_comment**(`+Boolean`)  
If `true` (default), call `prolog:comment_hook(Comments, TermPos, Term)` if this multifile hook is defined (see [prolog:comment_hook/3](loadfilehook.html#prolog:comment_hook/3)). This is used to drive PlDoc.

**comments**(`-Comments`)  
If provided, unify `Comments` with the comments encountered while reading `Term`. This option implies `process_comment(false)`.

The `singletons` option of [read_term/3](termrw.html#read_term/3) is initialised from the active style-checking mode. The `module` option is initialised to the current compilation module (see [prolog_load_context/2](consulting.html#prolog_load_context/2)).

\[ISO\]**read_term**(`-Term, +Options`)  
Read a term from the current input stream and unify the term with `Term`. The reading is controlled by options from the list of `Options`. If this list is empty, the behaviour is the same as for [read/1](termrw.html#read/1). The options are upward compatible with Quintus Prolog. The argument order is according to the ISO standard. Syntax errors are always reported using exception-handling (see [catch/3](exception.html#catch/3)). Options:

**backquoted_string**(`Bool`)  
If `true`, read `` ` ``...`` ` `` to a string object (see [section 5.2](string.html#sec:5.2)). The default depends on the Prolog flag [back_quotes](flags.html#flag:back_quotes).

**character_escapes**(`Bool`)  
Defines how to read `\` escape sequences in quoted atoms. See the Prolog flag [character_escapes](flags.html#flag:character_escapes) in [current_prolog_flag/2](flags.html#current_prolog_flag/2). (SWI-Prolog).

**comments**(`-Comments`)  
Unify `Comments` with a list of `Position`-`Comment`, where `Position` is a stream position object (see [stream_position_data/3](IO.html#stream_position_data/3)) indicating the start of a comment and `Comment` is a string object containing the text including delimiters of a comment. It returns all comments from where the [read_term/2](termrw.html#read_term/2) call started up to the end of the term read.

**cycles**(`Bool`)  
If `true` (default `false`), re-instantiate templates as produced by the corresponding [write_term/2](termrw.html#write_term/2) option. Note that the default is `false` to avoid misinterpretation of `@(Template, Substitutions)`, while the default of [write_term/2](termrw.html#write_term/2) is `true` because emitting cyclic terms without using the template construct produces an infinitely large term (read: it will generate an error after producing a huge amount of output).

**dotlists**(`Bool`)  
If `true` (default `false`), read `.(a,[])` as a list, even if lists are internally constructed a different functor (`[|](Head,Tail)`). This is primarily intended to read the output from [write_canonical/1](termrw.html#write_canonical/1) from other Prolog systems. See [section 5.1](ext-lists.html#sec:5.1).

**double_quotes**(`Atom`)  
Defines how to read " ... " strings. See the Prolog flag [double_quotes](flags.html#flag:double_quotes). (SWI-Prolog).

**module**(`Module`)  
Specify `Module` for operators, [character_escapes](flags.html#flag:character_escapes) flag and [double_quotes](flags.html#flag:double_quotes) flag. The value of the latter two is overruled if the corresponding [read_term/3](termrw.html#read_term/3) option is provided. If no module is specified, the current‘source module’is used. If the options is provided but the target module does not exist, module `user` is used because new modules by default inherit from `user`

**quasi_quotations**(`-List`)  
If present, unify `List` with the quasi quotations (see [section A.46](quasiquotations.html#sec:A.46)) instead of evaluating quasi quotations. Each quasi quotation is a term `quasi_quotation(+Syntax, +Quotation, +VarDict, -Result)`, where `Syntax` is the term in `{|Syntax||..|}`, `Quotation` is a list of character codes that represent the quotation, `VarDict` is a list of `Name`=`Variable` and `Result` is a variable that shares with the place where the quotation must be inserted. This option is intended to support tools that manipulate Prolog source text.

**singletons**(`Vars`)  
As `variable_names`, but only reports the variables occurring only once in the `Term` read (ISO). If `Vars` is the constant `warning`, singleton variables are reported using [print_message/2](printmsg.html#print_message/2). The variables appear in the order they have been read. The latter option provides backward compatibility and is used to read terms from source files. Not all singleton variables are reported as a warning. See [section 2.15.1.10](syntax.html#sec:2.15.1.10) for the rules that apply for warning about a singleton variable.^(113As of version 7.7.17, *all* variables starting with an underscore except for the truly anonymous variable are returned in `Vars`. Older versions only reported those that would have been reported if `warning` is used.)

**syntax_errors**(`Atom`)  
If `error` (default), throw an exception on a syntax error. Other values are `fail`, which causes a message to be printed using [print_message/2](printmsg.html#print_message/2), after which the predicate fails, `quiet` which causes the predicate to fail silently, and `dec10` which causes syntax errors to be printed, after which [read_term/\[2,3\]](termrw.html#read_term/2) continues reading the next term. Using `dec10`, [read_term/\[2,3\]](termrw.html#read_term/2) never fails. (Quintus, SICStus).

**subterm_positions**(`TermPos`)  
Describes the detailed layout of the term. The formats for the various types of terms are given below. All positions are character positions. If the input is related to a normal stream, these positions are relative to the start of the input; when reading from the terminal, they are relative to the start of the term.

**`From`-`To`**  
Used for primitive types (atoms, numbers, variables).

**string_position**(`From``, ``To`)  
Used to indicate the position of a string enclosed in double quotes (`"`).

**brace_term_position**(`From``, ``To``, ``Arg`)  
Term of the form `{...}`, as used in DCG rules. `Arg` describes the argument.

**list_position**(`From``, ``To``, ``Elms``, ``Tail`)  
A list. `Elms` describes the positions of the elements. If the list specifies the tail as `|`\<`TailTerm`\>, `Tail` is unified with the term position of the tail, otherwise with the atom `none`.

**term_position**(`From``, ``To``, ``FFrom``, ``FTo``, ``SubPos`)  
Used for a compound term not matching one of the above. `FFrom` and `FTo` describe the position of the functor. `SubPos` is a list, each element of which describes the term position of the corresponding subterm.

**dict_position**(`From``, ``To``, ``TagFrom``, ``TagTo``, ``KeyValuePosList`)  
Used for a dict (see [section 5.4](bidicts.html#sec:5.4)). The position of the key-value pairs is described by `KeyValuePosList`, which is a list of `key_value_position/7` terms. The `key_value_position/7` terms appear in the order of the input. Because maps do not preserve ordering, the key is provided in the position description.

**key_value_position**(`From``, ``To``, ``SepFrom``, ``SepTo``, ``Key``, ``KeyPos``, ``ValuePos`)  
Used for key-value pairs in a map (see [section 5.4](bidicts.html#sec:5.4)). It is similar to the `term_position/5` that would be created, except that the key and value positions do not need an intermediate list and the key is provided in `Key` to enable synchronisation of the file position data with the data structure.

**parentheses_term_position**(`From``, ``To``, ``ContentPos`)  
Used for terms between parentheses. This is an extension compared to the original Quintus specification that was considered necessary for secure refactoring of terms.

**quasi_quotation_position**(`From``, ``To``, ``SyntaxTerm``, ``SyntaxPos``, ``ContentPos`)  
Used for quasi quotations. Given the input `{|Syntax||Content|}`, `SyntaxTerm` is the parsed term representation from `Syntax`, e.g., `{|string(X)||Hello {{X}}|}` produces `Syntax` `string(X)` and `SyntaxPos` describes the layout of this term. `ContentPos` is always a term `From`-`To` describing the character range of `Content`.^(114The layout of the term produced by the quasi quotation parser is not available. Future versions may provide an interface that allows contributing a layout term.)

**term_position**(`Pos`)  
Unifies `Pos` with the starting position of the term read. `Pos` is of the same format as used by [stream_property/2](IO.html#stream_property/2).

**var_prefix**(`Bool`)  
If `true`, demand variables to start with an underscore. See [section 2.15.1.8](syntax.html#sec:2.15.1.8).

**variables**(`Vars`)  
Unify `Vars` with a list of variables in the term. The variables appear in the order they have been read. See also [term_variables/2](manipterm.html#term_variables/2). (ISO).

**variable_names**(`Vars`)  
Unify `Vars` with a list of‘`Name` = `Var`’, where `Name` is an atom describing the variable name and `Var` is a variable that shares with the corresponding variable in `Term`. (ISO). The variables appear in the order they have been read.

\[ISO\]**read_term**(`+Stream, -Term, +Options`)  
Read term with options from `Stream`. See [read_term/2](termrw.html#read_term/2).

**read_term_from_atom**(`+Atom, -Term, +Options`)  
Use [read_term/3](termrw.html#read_term/3) to read the next term from `Atom`. `Atom` is either an atom or a string object (see [section 5.2](string.html#sec:5.2)). It is not required for `Atom` to end with a full-stop. If `Atom` only contains white space and/or comments, an `syntax_error(end_of_string)` exception is raised. This predicate supersedes [atom_to_term/3](manipatom.html#atom_to_term/3).

**read_term_with_history**(`-Term, +Options`)  
Read a term while providing history substitutions. [read_term_with_history/2](termrw.html#read_term_with_history/2) is used by the top level to read the user's actions. In addition to the options recognised by [read_term/2](termrw.html#read_term/2), the following options are recognised:

**prompt**(`+Prompt`)  
Define the prompt to use. The default is `~! ?-`. A sequence `~!` is replaced by the current history event number.

**show**(`+Command`)  
Using `Command` lists the saved history events. Default is `!history`.

**help**(`+Command`)  
Using `Command` shows help on the history system. Default is `!help`.

**no_save**(`+Commands`)  
Do not save the command into the history if it appears in the list `Commands`.

**module**(`+Module`)  
Defines the module from which to extract module-specific syntax such as operators and handling of the various quotes. Default is the *typein* module which is set using [module/1](mtoplevel.html#module/1) and is initially set to `user`.

**input**(`+Stream`)  
Stream from which to read `Term`. Default is `user_input`.

Most applications will use the [read_term/2](termrw.html#read_term/2) option `variable_names` to get access to the names of the variables in `Term`. SWI-Prolog calls [read_term_with_history/2](termrw.html#read_term_with_history/2) as follows:

``` code
    read_term_with_history(
        Goal,
        [ show(h),
          help('!h'),
          no_save([trace, end_of_file]),
          prompt('~! ?-'),
          variable_names(Bindings)
        ]).
```

**prompt**(`-Old, +New`)  
Set prompt associated with reading from the `user_input` stream. `Old` is first unified with the current prompt. On success the prompt will be set to `New` (an atom). A prompt is printed if data is read from `user_input`, the cursor is at the left margin and the `user_input` is considered to be connected to a terminal. See the `tty(Bool)` property of [stream_property/2](IO.html#stream_property/2) and [set_stream/2](IO.html#set_stream/2).

The default prompt is `'|: '`. Note that the toplevel loop (see [prolog/0](toplevel.html#prolog/0)) sets the prompt for the first prompt (see [prompt1/1](termrw.html#prompt1/1)) to `'?- '`, possibly decorated by the history event number, *break level* and debug mode. If the first line does not complete the term, subsequent lines are prompted for using the prompt as defined by [prompt/2](termrw.html#prompt/2).

**prompt1**(`?Prompt`)  
Sets the prompt for the next line to be read. Continuation lines will be read using the prompt defined by [prompt/2](termrw.html#prompt/2). If `Prompt` is unbound it is unified with the current first-line prompt.
