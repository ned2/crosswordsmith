
## 4.24 Character properties

SWI-Prolog offers two comprehensive predicates for classifying characters and character codes. These predicates are defined as built-in predicates to exploit the C-character classification's handling of *locale* (handling of local character sets). These predicates are fast, logical and deterministic if applicable.

In addition, there is the library `library(ctypes)` providing compatibility with some other Prolog systems. The predicates of this library are defined in terms of [code_type/2](chartype.html#code_type/2).

**code_type**(`?Code, ?Type`)  
**char_type**(`?Char, ?Type`)  
Tests or generates alternative `Type`s or `Char/Code`s. The character types are inspired by the standard C `<ctype.h>` primitives. The types are sensitive to the active *locale*, see [setlocale/3](system.html#setlocale/3). Most of the `Type`s are mapped to the Unicode classification functions from `<wctype.h>`, e.g., `alnum` uses **iswalnum()**. The types `prolog_var_start`, `prolog_atom_start`, `prolog_identifier_continue` and `prolog_symbol` are based on the locale-independent built-in classification routines that are also used by [read/1](termrw.html#read/1) and friends.

Note that the mode (-,+) is only efficient if the `Type` has a parameter, e.g., `char_type(C, digit(8))`. If `Type` is a atomic, the whole unicode range (0..0x1ffff) is generated and tested against the character classification function.

**alnum**  
`Char` is a letter (upper- or lowercase) or digit.

**alpha**  
`Char` is a letter (upper- or lowercase).

**csym**  
`Char` is a letter (upper- or lowercase), digit or the underscore (`_`). These are valid C and Prolog symbol characters.

**csymf**  
`Char` is a letter (upper- or lowercase) or the underscore (`_`). These are valid first characters for C and Prolog symbols.

**ascii**  
`Char` is a 7-bit ASCII character (0..127).

**white**  
`Char` is a space or tab, i.e. white space inside a line.

**cntrl**  
`Char` is an ASCII control character (0..31), ASCII DEL character (127), or non-ASCII character in the range 128..159 or 8232..8233.

**digit**  
`Char` is a digit, i.e., `Char` is in `0 ... 9`. See also `decimal`.

**digit**(`Weight`)  
`Char` is a digit with value `Weight`. I.e. `char_type(X, digit(6))` yields `X` = `’6’`. Useful for parsing numbers.

**xdigit**(`Weight`)  
`Char` is a hexadecimal digit with value `Weight`. I.e. `char_type(a, xdigit(X))` yields `X` = `’10’`. Useful for parsing numbers.

**decimal**  
`Char` is a decimal digit in any script. This implies it has the Unicode *general category* **Nd**).

**decimal**(`Weight`)  
`Char` is a decimal digit in any script with `Weight` `0 ... 9`.

**print**  
`Char` is printable character.

**graph**  
`Char` produces a visible mark on a page when printed. Note that the space is not included!

**lower**  
`Char` is a lowercase letter.

**lower**(`Upper`)  
`Char` is a lowercase version of `Upper`. Only true if `Char` is lowercase and `Upper` uppercase.

**to_lower**(`Upper`)  
`Char` is a lowercase version of `Upper`. For non-letters, or letter without case, `Char` and `Lower` are the same. See also [upcase_atom/2](chartype.html#upcase_atom/2) and [downcase_atom/2](chartype.html#downcase_atom/2).

**upper**  
`Char` is an uppercase letter.

**upper**(`Lower`)  
`Char` is an uppercase version of `Lower`. Only true if `Char` is uppercase and `Lower` lowercase.

**to_upper**(`Lower`)  
`Char` is an uppercase version of `Lower`. For non-letters, or letter without case, `Char` and `Lower` are the same. See also [upcase_atom/2](chartype.html#upcase_atom/2) and [downcase_atom/2](chartype.html#downcase_atom/2).

**punct**  
`Char` is a punctuation character. This is a `graph` character that is not a letter or digit.

**space**  
`Char` is some form of layout character (tab, vertical tab, newline, etc.).

**end_of_file**  
`Char` is -1.

**end_of_line**  
`Char` is one of the four ISO/POSIX line-ending control codes U+000A LF, U+000B VT, U+000C FF, U+000D CR. This is the original ISO Prolog and C-string-literal definition; `prolog_end_of_line` is the wider set used by the SWI-Prolog reader.

**newline**  
`Char` is a newline character (10).

**period**  
`Char` counts as the end of a sentence (.,!,?).

**quote**  
`Char` is a quote character (`"`, `'`, `` ` ``).

**paren**(`Close`)  
`Char` is an opening bracket and `Close` is its matching close. Covers the three ASCII bracket pairs `()`, `[]` and `{}`, plus every Unicode `Ps`/`Pe` pair (about 60 pairs in Unicode 17, including angle, corner, ceiling, floor, mathematical, ornamental, fullwidth and CJK brackets). The mapping is reversible: with `Close` bound, `Char` unifies with the matching open.

**quote**(`Close`)  
`Char` is an opening quotation mark and `Close` is its matching close. The ASCII quotes `’`, `"`, and `‘` have `Close` = `Char`; Unicode `Pi`/`Pf` quote pairs (the guillemets, the standard left/right curly single and double quotes, and the single/double angle and reversed quotation marks) have `Close` different from `Char`. The mapping is reversible.

**width**(`Width`)  
`Width` is the number of columns for fixed-width usage used by `Char`. True for all printable characters. Most characters require 1 column. *Unicode combining characters* require no space (they follow the *base character*). Many Asian characters and the Emojis require 2 columns. These values are used by [stream_property/2](IO.html#stream_property/2) for the `position(-Pos)` property.

**prolog_layout**  
`Char` is a Prolog *layout* character: a member of the Unicode `Pattern_White_Space` set used by [read_term/2](termrw.html#read_term/2) to separate tokens. The eleven code points are U+0009..U+000D, U+0020, U+0085, U+200E, U+200F, U+2028 and U+2029. Locale- independent; pinned to [unicode_syntax_version](flags.html#flag:unicode_syntax_version). `prolog_end_of_line` is the seven-element line-terminator subset.

**prolog_end_of_line**  
`Char` ends a line of Prolog source text. Covers the seven line-terminator-like `Pattern_White_Space` code points: U+000A (LF), U+000B (VT), U+000C (FF), U+000D (CR), U+0085 (NEL), U+2028 (LINE SEPARATOR), and U+2029 (PARAGRAPH SEPARATOR). The same set terminates `%` comments and increments the source line counter. See [section 2.15.1.9](syntax.html#sec:2.15.1.9).

**prolog_var_start**  
`Char` can start a Prolog variable name.

**prolog_atom_start**  
`Char` can start a unquoted Prolog atom that is not a symbol.

**prolog_identifier_continue**  
`Char` can continue a Prolog variable name or atom.

**prolog_symbol**  
`Char` is a Prolog symbol character. Sequences of Prolog symbol characters glue together to form an unquoted atom. Examples are `=..`, `\=`, etc.

**prolog_solo**  
`Char` is a Prolog *solo* character: a punctuation code point that forms an atom on its own and never combines with neighbouring symbol characters. In ASCII the solo set is `!`, `;` and `%`; the same flag carries over to non-ASCII code points via the Unicode syntax map (see [section 2.15.1.9](syntax.html#sec:2.15.1.9)). Solo characters are written unquoted by [writeq/1](termrw.html#writeq/1) and are accepted as single-character atoms by the reader.

**pattern_syntax**  
`Char` has the Unicode `Pattern_Syntax` property (UAX #31 R3). This is the immutable set of punctuation and symbol code points whose classification is guaranteed not to change across Unicode versions. Used by [write_canonical/1](termrw.html#write_canonical/1) and the `pattern_syntax_solo` option of [write_term/2](termrw.html#write_term/2) to decide which single-character atoms can be printed bare with round-trip safety across Unicode upgrades.

### 4.24.1 Case conversion

There is nothing in the Prolog standard for converting case in textual data. The SWI-Prolog predicates [code_type/2](chartype.html#code_type/2) and [char_type/2](chartype.html#char_type/2) can be used to test and convert individual characters. We have started some additional support:

**downcase_atom**(`+AnyCase, -LowerCase`)  
Converts the characters of `AnyCase` into lowercase as [char_type/2](chartype.html#char_type/2) does (i.e. based on the defined *locale* if Prolog provides locale support on the hosting platform) and unifies the lowercase atom with `LowerCase`.

**upcase_atom**(`+AnyCase, -UpperCase`)  
Converts, similar to [downcase_atom/2](chartype.html#downcase_atom/2), an atom to uppercase.

### 4.24.2 White space normalization

**normalize_space**(`-Out, +In`)  
Normalize white space in `In`. All leading and trailing white space is removed. All non-empty sequences for Unicode white space characters are replaced by a single space (`\u0020`) character. `Out` uses the same conventions as [with_output_to/2](IO.html#with_output_to/2) and [format/3](format.html#format/3).

### 4.24.3 Language-specific comparison

This section deals with predicates for language-specific string comparison operations.

**collation_key**(`+Atom, -Key`)  
Create a `Key` from `Atom` for locale-specific comparison. The key is defined such that if the key of atom `A` precedes the key of atom `B` in the standard order of terms, `A` is alphabetically smaller than `B` using the sort order of the current locale.

The predicate [collation_key/2](chartype.html#collation_key/2) is used by [locale_sort/2](chartype.html#locale_sort/2) from library(sort). Please examine the implementation of [locale_sort/2](chartype.html#locale_sort/2) as an example of using this call.

The `Key` is an implementation-defined and generally unreadable string. On systems that do not support locale handling, `Key` is simply unified with `Atom`.

**locale_sort**(`+List, -Sorted`)  
Sort a list of atoms using the current locale. `List` is a list of atoms or string objects (see [section 5.2](string.html#sec:5.2)). `Sorted` is unified with a list containing all atoms of `List`, sorted to the rules of the current locale. See also [collation_key/2](chartype.html#collation_key/2) and [setlocale/3](system.html#setlocale/3).
