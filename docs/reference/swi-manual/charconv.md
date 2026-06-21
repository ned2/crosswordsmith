
## 4.26 Character Conversion

Although I wouldn't really know why you would like to use these features, they are provided for ISO compliance.

\[ISO\]**char_conversion**(`+CharIn, +CharOut`)  
Define that term input (see [read_term/3](termrw.html#read_term/3)) maps each character read as `CharIn` to the character `CharOut`. Character conversion is only executed if the Prolog flag [char_conversion](flags.html#flag:char_conversion) is set to `true` and not inside quoted atoms or strings. The initial table maps each character onto itself. See also [current_char_conversion/2](charconv.html#current_char_conversion/2).

\[ISO\]**current_char_conversion**(`?CharIn, ?CharOut`)  
Queries the current character conversion table. See [char_conversion/2](charconv.html#char_conversion/2) for details.
