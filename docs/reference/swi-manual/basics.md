
## A.12 library(dcg/basics): Various general DCG utilities

To be done  
This is just a starting point. We need a comprehensive set of generally useful DCG primitives.

This library provides various commonly used DCG primitives acting on list of character **codes**. Character classification is based on [code_type/2](chartype.html#code_type/2).

This module started its life as `library(http/dcg_basics)` to support the HTTP protocol. Since then, it was increasingly used in code that has no relation to HTTP and therefore this library was moved to the core library.

\[det\]**string_without**(`+EndCodes, -Codes`)`//`  
Take as many codes from the input until the next character code appears in the list `EndCodes`. The terminating code itself is left on the input. Typical use is to read upto a defined delimiter such as a newline or other reserved character. For example:

``` code
    ...,
    string_without("\n", RestOfLine)
```

|            |                               |
|------------|-------------------------------|
| `EndCodes` | is a list of character codes. |

See also  
[string//1](basics.html#string//1).

\[nondet\]**string**(`-Codes`)`//`  
Take as few as possible tokens from the input, taking one more each time on backtracking. This code is normally followed by a test for a delimiter. For example:

``` code
upto_colon(Atom) -->
        string(Codes), ":", !,
        { atom_codes(Atom, Codes) }.
```

See also  
[string_without//2](basics.html#string_without//2).

\[det\]**blanks**`//`  
Skip zero or more white-space characters.

\[semidet\]**blank**`//`  
Take next `space` character from input. Space characters include newline.

See also  
[white//0](basics.html#white//0)

\[det\]**nonblanks**(`-Codes`)`//`  
Take all `graph` characters

\[semidet\]**nonblank**(`-Code`)`//`  
`Code` is the next non-blank (`graph`) character.

\[semidet\]**blanks_to_nl**`//`  
Take a sequence of [blank//0](basics.html#blank//0) codes if blanks are followed by a newline or end of the input.

\[det\]**whites**`//`  
Skip white space *inside* a line.

See also  
[blanks//0](basics.html#blanks//0) also skips newlines.

\[semidet\]**white**`//`  
Take next `white` character from input. White characters do *not* include newline.

\[semidet\]**alpha_to_lower**(`?C`)`//`  
Read a letter (class `alpha`) and return it as a lowercase letter. If `C` is instantiated and the DCG list is already bound, `C` must be `lower` and matches both a lower and uppercase letter. If the output list is unbound, its first element is bound to `C`. For example:

``` code
?- alpha_to_lower(0'a, `AB`, R).
R = [66].
?- alpha_to_lower(C, `AB`, R).
C = 97, R = [66].
?- alpha_to_lower(0'a, L, R).
L = [97|R].
```

\[det\]**digits**(`?Chars`)`//`  
\[det\]**digit**(`?Char`)`//`  
\[det\]**integer**(`?Integer`)`//`  
Number processing. The predicate [digits//1](basics.html#digits//1) matches a possibly empty set of digits, [digit//1](basics.html#digit//1) processes a single digit and integer processes an optional sign followed by a non-empty sequence of digits into an integer.

\[det\]**float**(`?Float`)`//`  
Process a floating point number. The actual conversion is controlled by [number_codes/2](manipatom.html#number_codes/2).

\[det\]**number**(`+Number`)`//`  
\[semidet\]**number**(`-Number`)`//`  
Generate extract a number. Handles both integers and floating point numbers.

\[det\]**xinteger**(`+Integer`)`//`  
\[semidet\]**xinteger**(`-Integer`)`//`  
Generate or extract an integer from a sequence of hexadecimal digits. Hexadecimal characters include both uppercase (A-F) and lowercase (a-f) letters. The value may be preceded by a sign (+/-)

\[semidet\]**xdigit**(`-Weight`)`//`  
True if the next code is a hexdecimal digit with `Weight`. `Weight` is between 0 and 15. Hexadecimal characters include both uppercase (A-F) and lowercase (a-f) letters.

\[det\]**xdigits**(`-WeightList`)`//`  
List of weights of a sequence of hexadecimal codes. `WeightList` may be empty. Hexadecimal characters include both uppercase (A-F) and lowercase (a-f) letters.

**eol**`//`  
Matches end-of-line. Matching `\`r`\`n, `\n` or end of input ([eos//0](basics.html#eos//0)).

**eos**`//`  
Matches end-of-input. The implementation behaves as the following portable implementation:

``` code
eos --> call(eos_).
eos_([], []).
```

To be done  
This is a difficult concept and violates the *context free* property of DCGs. Explain the exact problems.

**remainder**(`-List`)`//`  
Unify `List` with the remainder of the input.

\[semidet\]**prolog_var_name**(`-Name:atom`)`//`  
Matches a Prolog variable name. Primarily intended to deal with quasi quotations that embed Prolog variables.

\[semidet\]**csym**(`?Symbol:atom`)`//`  
Recognise a C symbol according to the `csymf` and `csym` code type classification provided by the C library.

\[det\]**atom**(`++Atom`)`//`  
Generate codes of `Atom`. Current implementation uses [write/1](termrw.html#write/1), dealing with any Prolog term. `Atom` must be ground though.
