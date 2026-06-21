
## 4.2 Character representation

In traditional (Edinburgh) Prolog, characters are represented using *character codes*. Character codes are integer indices into a specific character set. Traditionally the character set was 7-bit US-ASCII. 8-bit character sets have been allowed for a long time, providing support for national character sets, of which iso-latin-1 (ISO 8859-1) is applicable to many Western languages.

ISO Prolog introduces three types, two of which are used for characters and one for accessing binary streams (see [open/4](IO.html#open/4)). These types are:

- *code*  
  A *character code* is an integer representing a single character. As files may use multi-byte encoding for supporting different character sets (utf-8 encoding for example), reading a code from a text file is in general not the same as reading a byte.
- *char*  
  Alternatively, characters may be represented as *one-character atoms*. This is a natural representation, hiding encoding problems from the programmer as well as providing much easier debugging.
- *byte*  
  Bytes are used for accessing binary streams.

In SWI-Prolog, character codes are *always* the Unicode equivalent of the encoding. That is, if [get_code/1](chario.html#get_code/1) reads from a stream encoded as `KOI8-R` (used for the Cyrillic alphabet), it returns the corresponding Unicode code points. Similarly, assembling or disassembling atoms using [atom_codes/2](manipatom.html#atom_codes/2) interprets the codes as Unicode points. See [section 2.18.1](widechars.html#sec:2.18.1) for details.

To ease the pain of the two character representations (code and char), SWI-Prolog's built-in predicates dealing with character data work as flexible as possible: they accept data in any of these formats as long as the interpretation is unambiguous. In addition, for output arguments that are instantiated, the character is extracted before unification. This implies that the following two calls are identical, both testing whether the next input character is an `a`.

``` code
peek_code(Stream, a).
peek_code(Stream, 97).
```

The two character representations are handled by a large number of built-in predicates, all of which are ISO-compatible. For converting between code and character there is [char_code/2](manipatom.html#char_code/2). For breaking atoms and numbers into characters there are [atom_chars/2](manipatom.html#atom_chars/2), [atom_codes/2](manipatom.html#atom_codes/2), [number_chars/2](manipatom.html#number_chars/2) and [number_codes/2](manipatom.html#number_codes/2). For character I/O on streams there are [get_char/\[1,2\]](chario.html#get_char/1), [get_code/\[1,2\]](chario.html#get_code/1), [get_byte/\[1,2\]](chario.html#get_byte/1), [peek_char/\[1,2\]](chario.html#peek_char/1), [peek_code/\[1,2\]](chario.html#peek_code/1), [peek_byte/\[1,2\]](chario.html#peek_byte/1), [put_code/\[1,2\]](chario.html#put_code/1), [put_char/\[1,2\]](chario.html#put_char/1) and [put_byte/\[1,2\]](chario.html#put_byte/1). The Prolog flag [double_quotes](flags.html#flag:double_quotes) controls how text between double quotes is interpreted.
