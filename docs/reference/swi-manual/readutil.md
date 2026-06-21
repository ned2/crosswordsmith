
## A.49 library(readutil): Read utilities

See also  
\- `library(pure_input)` allows for processing files with DCGs.  
- `library(lazy_lists)` for creating lazy lists from input.

This library provides some commonly used reading predicates. As these predicates have proven to be time-critical in some applications we moved them to C. For compatibility as well as to reduce system dependency, we link the foreign code at runtime and fallback to the Prolog implementation if the shared object cannot be found.

\[det\]**read_line_to_codes**(`+Stream, -Line:codes`)  
Read the next line of input from `Stream`. Unify content of the lines as a list of character codes with `Line` *after* the line has been read. A line is ended by a newline character or end-of-file. Unlike [read_line_to_codes/3](readutil.html#read_line_to_codes/3), this predicate removes a trailing newline character.

\[det\]**read_line_to_codes**(`+Stream, -Line, ?Tail`)  
Difference-list version to read an input line to a list of character codes. Reading stops at the newline or end-of-file character, but unlike [read_line_to_codes/2](readutil.html#read_line_to_codes/2), the newline is retained in the output. This predicate is especially useful for reading a block of lines up to some delimiter. The following example reads an HTTP header ended by a blank line:

``` code
read_header_data(Stream, Header) :-
    read_line_to_codes(Stream, Header, Tail),
    read_header_data(Header, Stream, Tail).

read_header_data("\r\n", _, _) :- !.
read_header_data("\n", _, _) :- !.
read_header_data("", _, _) :- !.
read_header_data(_, Stream, Tail) :-
    read_line_to_codes(Stream, Tail, NewTail),
    read_header_data(Tail, Stream, NewTail).
```

\[det\]**read_line_to_string**(`+Stream, -String`)  
Read the next line from `Stream` into `String`. `String` does not contain the line terminator. `String` is unified with the *atom* `end_of_file` if the end of the file is reached.

See also  
[read_string/5](string.html#read_string/5) can be used to read lines with separated records without creating intermediate strings.

\[det\]**read_stream_to_codes**(`+Stream, -Codes`)  
\[det\]**read_stream_to_codes**(`+Stream, -Codes, ?Tail`)  
Read input from `Stream` to a list of character codes. The version [read_stream_to_codes/3](readutil.html#read_stream_to_codes/3) creates a difference-list.

\[det\]**read_file_to_codes**(`+Spec, -Codes, +Options`)  
Read the file `Spec` into a list of `Codes`. `Options` is split into options for [absolute_file_name/3](files.html#absolute_file_name/3) and [open/4](IO.html#open/4). In addition, the following option is provided:

**tail**(`?Tail`)  
Read the data into a *difference list* `Codes``\``Tail`.

See also  
[phrase_from_file/3](pio.html#phrase_from_file/3) and [read_file_to_string/3](readutil.html#read_file_to_string/3).

\[det\]**read_file_to_string**(`+Spec, -String, +Options`)  
Read the file `Spec` into a the string `String`. `Options` is split into options for [absolute_file_name/3](files.html#absolute_file_name/3) and [open/4](IO.html#open/4).

See also  
[phrase_from_file/3](pio.html#phrase_from_file/3) and [read_file_to_codes/3](readutil.html#read_file_to_codes/3).

\[det\]**read_file_to_terms**(`+Spec, -Terms, +Options`)  
Read the file `Spec` into a list of terms. `Options` is split over [absolute_file_name/3](files.html#absolute_file_name/3), [open/4](IO.html#open/4) and [read_term/3](termrw.html#read_term/3). In addition, the following option is processed:

**tail**(`?Tail`)  
If present, `Terms``\``Tail` forms a *difference list*.

Note that the *output* options of [read_term/3](termrw.html#read_term/3), such as `variable_names` or `subterm_positions` will cause [read_file_to_terms/3](readutil.html#read_file_to_terms/3) to fail if `Spec` contains multiple terms because the values for the different terms will not unify.
