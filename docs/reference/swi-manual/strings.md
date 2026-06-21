
## A.55 library(strings): String utilities

See also  
\- [format/3](format.html#format/3) can format to a string as well. The `library(lynx/format)` provides primitive to wrap long strings.  
- The core system provides many additional string processing predicates.

To be done  
There are probably many other high level string predicates that belong in this library. For example, predicates similar to the functions in [https://docs.python.org/3/library/textwrap.html](https://docs.python.org/3/library/textwrap.html)

This module provides string handling utilities, currently notably for dealing with multi-line strings and *interpolation*. The library provides a couple of primitives as well definitions for the `string` *quasi quotation* syntax. The latter allows for constructing both single line and multi-line long strings based on template interpolation. Below is a simple example using the quasi quotation syntax.

``` code
test(To) :-
    write({|string(To)||
           | Dear {To},
           |
           | I'm happy to announce a string interpolation quasi quoter.
           |}.
```

**Warning**

The general purpose string interpolation implemented by this library should **not** be used to create strings for a formal language such as HTML, JavaScript, SQL, etc. because the result will be subject to **injection attacks**, providing a serious **security risc**. The core idea of quasi quotation is to know about the target language and interpolate Prolog data into the template **while respecting the syntax of the target language**, notable to **escape certain characters where needed**. See also `library(http/html_write)` and `library(http/js_write)` which define quasi quotation rules for HTML and JavaScript.

**string**(`+Content, +Args, +Binding, -DOM`)  
Implements the quasi quotation syntax `string`. If the first character of the content is a newline (i.e., there is a newline *immediately* after the `||` token) this first uses [dedent_lines/3](strings.html#dedent_lines/3) to the remove common white space prefix from the lines. This is called with the option `chars("\s\t|")`, i.e., also removing `|` characters and `tab(8)`.

If the quasi quotation syntax carries arguments (e.g., `string(To)`), the string is compiled into a function that produces the result of interpolating the arguments into the template. See user functions on dict objects. If there are no arguments, the result is simply the final string.

See also  
\- [interpolate_string/4](strings.html#interpolate_string/4) for the interpolation syntax.  
- Section for examples and discussion.

To be done  
Specify tab width and allow for {@Goal} templates.

**interpolate_string**(`:In, -Out, +Map, +Options`)  
Establish a string from a template by replacing patterns. Supported patterns are:

{`Name`}  
If `Map` contains `Name=Value`, insert `Value` using [write/1](termrw.html#write/1). If `Name` does not appear in `Map`, raise an existence error. `Name` must satisfy the rules for a Prolog variable.

{`Name,Default`}  
As above, but if `Name` does not appear in `Map`, use `Value`

{`@(Goal)`}  
Insert the output (to `current_output`) of `Goal` here. For safety reasons only accepted if `Options` contains `goals(true)`

\[det\]**string_lines**(`?String, ?Lines`)  
True when `String` represents `Lines`. This follows the normal text convention that a line is defined as a possible empty string followed by a newline character ("`\n`"). E.g.

``` code
?- string_lines("a\nb\n", L).
L = ["a", "b"].
?- string_lines(S, ["a", "b"]).
S = "a\nb\n".
```

This predicate is a true *relation* if both arguments are in canonical form, i.e. all text is represented as strings and the first argument ends with a newline. The implementation tolerates non-canonical input: other types than strings are accepted and `String` does not need to end with a newline.

See also  
[split_string/4](string.html#split_string/4). Using `split_string(String, "\n", "", Lines)` on a string that ends in a newline adds an additional empty string compared to [string_lines/2](strings.html#string_lines/2).

**dedent_lines**(`+In, -Out, +Options`)  
Remove shared indentation for all lines in a string. Lines are separated by "`\n`" -- conversion to and from external forms (such as "`\`r`\n`") are typically done by the I/O predicates. A final "`\n`" is preserved.

`Options`:

**tab**(`N`)  
Assume tabs at columns of with `N`. When omitted, tabs are taken literally and only exact matches are removed.

**chars**(`CodesOrString`)  
Characters to remove. This can notably be used to remove additional characters such as `*` or‘`|`‘. Default is `" \t"`.

\[det\]**indent_lines**(`+Prefix, +In, -Out`)  
Add `Prefix` to the beginning of lines in `In`. Lines are separated by "`\n`" -- conversion to and from external forms (such as "`\`r`\n`") are typically done by the I/O predicates. Lines that consist entirely of whitespace are left as-is.

\[det\]**indent_lines**(`:Filter, +Prefix, +In, -Out`)  
Similar to [indent_lines/3](strings.html#indent_lines/3), but only adds `Prefix` to lines for which `call(Filter, Line)` succeeds.
