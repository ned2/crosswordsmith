
## 4.24 Character properties

SWI-Prolog offers two comprehensive predicates for classifying characters and character codes. These predicates are defined as built-in predicates to exploit the C-character classification's handling of *locale* (handling of local character sets). These predicates are fast, logical and deterministic if applicable.

In addition, there is the library `library(ctypes)` providing compatibility with some other Prolog systems. The predicates of this library are defined in terms of [code_type/2](chartype.html#code_type/2).

**char_type**(`?Char, ?Type`)  
Tests or generates alternative `Type`s or `Char`s. The character types are inspired by the standard C `<ctype.h>` primitives. The types are sensitive to the active *locale*, see [setlocale/3](system.html#setlocale/3). Most of the `Type`s are mapped to the Unicode classification functions from `<wctype.h>`, e.g., `alnum` uses **iswalnum()**. The types `prolog_var_start`, `prolog_atom_start`, `prolog_identifier_continue` and `prolog_symbol` are based on the locale-independent built-in classification routines that are also used by [read/1](termrw.html#read/1) and friends.

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
`Char` is a digit, i.e., `Char` is in `0 ...`. See also `decimal`.

**digit**(`Weight`)  
`Char` is a digit with value `Weight`. I.e. `char_type(X, digit(6))` yields `X` = `’6’`. Useful for parsing numbers.

**xdigit**(`Weight`)  
`Char` is a hexadecimal digit with value `Weight`. I.e. `char_type(a, xdigit(X))` yields `X` = `’10’`. Useful for parsing numbers.

**decimal**  
`Char` is a decimal digit in any script. This implies it has the Unicode *general category* **Nd**).

**decimal**(`Weight`)  
`Char` is a decimal digit in any script with `Weight` `0 ...`.

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
`Char` ends a line (ASCII: 10..13).

**newline**  
`Char` is a newline character (10).

**period**  
`Char` counts as the end of a sentence (.,!,?).

**quote**  
`Char` is a quote character (`"`, `'`, `` ` ``).

**paren**(`Close`)  
`Char` is an open parenthesis and `Close` is the corresponding close parenthesis.

**prolog_var_start**  
`Char` can start a Prolog variable name.

**prolog_atom_start**  
`Char` can start a unquoted Prolog atom that is not a symbol.

**prolog_identifier_continue**  
`Char` can continue a Prolog variable name or atom.

**prolog_symbol**  
`Char` is a Prolog symbol character. Sequences of Prolog symbol characters glue together to form an unquoted atom. Examples are `=..`, `\=`, etc.

**code_type**(`?Code, ?Type`)  
As [char_type/2](chartype.html#char_type/2), but uses character codes rather than one-character atoms. Please note that both predicates are as flexible as possible. They handle either representation if the argument is instantiated and will instantiate only with an integer code or a one-character atom, depending of the version used. See also the Prolog flag [double_quotes](flags.html#flag:double_quotes), [atom_chars/2](manipatom.html#atom_chars/2) and [atom_codes/2](manipatom.html#atom_codes/2).

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
