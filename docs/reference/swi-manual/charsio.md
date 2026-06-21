
## A.6 library(charsio): I/O on Lists of Character Codes

Compatibility  
The naming of this library is not in line with the ISO standard. We believe that the SWI-Prolog native predicates form a more elegant alternative for this library.

This module emulates the Quintus/SICStus library `charsio.pl` for reading and writing from/to lists of character codes. Most of these predicates are straight calls into similar SWI-Prolog primitives. Some can even be replaced by ISO standard predicates.

\[det\]**format_to_chars**(`+Format, +Args, -Codes`)  
Use [format/2](format.html#format/2) to write to a list of character codes.

\[det\]**format_to_chars**(`+Format, +Args, -Codes, ?Tail`)  
Use [format/2](format.html#format/2) to write to a difference list of character codes.

**write_to_chars**(`+Term, -Codes`)  
Write a term to a code list. True when `Codes` is a list of character codes written by [write/1](termrw.html#write/1) on `Term`.

**write_to_chars**(`+Term, -Codes, ?Tail`)  
Write a term to a code list. `Codes``\``Tail` is a difference list of character codes produced by [write/1](termrw.html#write/1) on `Term`.

\[det\]**atom_to_chars**(`+Atom, -Codes`)  
Convert `Atom` into a list of character codes.

deprecated  
Use ISO [atom_codes/2](manipatom.html#atom_codes/2).

\[det\]**atom_to_chars**(`+Atom, -Codes, ?Tail`)  
Convert `Atom` into a difference list of character codes.

\[det\]**number_to_chars**(`+Number, -Codes`)  
Convert Atom into a list of character codes.

deprecated  
Use ISO [number_codes/2](manipatom.html#number_codes/2).

\[det\]**number_to_chars**(`+Number, -Codes, ?Tail`)  
Convert `Number` into a difference list of character codes.

\[det\]**read_from_chars**(`+Codes, -Term`)  
Read `Codes` into `Term`.

Compatibility  
The SWI-Prolog version does not require `Codes` to end in a full-stop.

\[det\]**read_term_from_chars**(`+Codes, -Term, +Options`)  
Read `Codes` into `Term`. `Options` are processed by [read_term/3](termrw.html#read_term/3).

Compatibility  
sicstus

\[det\]**open_chars_stream**(`+Codes, -Stream`)  
Open `Codes` as an input stream.

See also  
[open_string/2](string.html#open_string/2).

\[det\]**with_output_to_chars**(`:Goal, -Codes`)  
Run `Goal` as with [once/1](metacall.html#once/1). Output written to `current_output` is collected in `Codes`.

\[det\]**with_output_to_chars**(`:Goal, -Codes, ?Tail`)  
Run `Goal` as with [once/1](metacall.html#once/1). Output written to `current_output` is collected in `Codes``\``Tail`.

\[det\]**with_output_to_chars**(`:Goal, -Stream, -Codes, ?Tail`)  
Same as [with_output_to_chars/3](charsio.html#with_output_to_chars/3) using an explicit stream. The difference list `Codes``\``Tail` contains the character codes that `Goal` has written to `Stream`.
