
## A.29 library(writef): Old-style formatted write

copyright  
Copied from Edinburgh C-Prolog. Original version by Byrd, changed many times since.

This library provides [writef/1](writef.html#writef/1) and friends. These predicates originate from Edinburgh C-Prolog and and provided for compatibility purposes. New code should use [format/1](format.html#format/1), [format/2](format.html#format/2) and friends, which are currently supported by more Prolog implementations.

The writef-family of predicates conflicts with the modern *character-escapes* flag about the interpretation of `\`-sequences. This can be avoided by

1.  Disable character escapes (not recommended unless one wants to run really outdated code unmodified).
2.  Double the `\` for conflicting interpretations
3.  Use ISO compliant alternatives for conflicting interpretations

\[det\]**writef**(`+Format`)  
\[det\]**writef**(`+Format, +Arguments`)  
Formatted write to the `current_output`. `Format` is a format specifier. Some escape sequences require arguments that must be provided in the list `Arguments`. There are two types of escape sequences: special characters start with `\` and include arguments start with `%`. The special character sequences are:

> |          |                                                           |
> |----------|-----------------------------------------------------------|
> | `\n`     | Output a newline character                                |
> | `\l`     | Output a line separator (same as `\n`)                    |
> | `\r`     | Output a carriage-return character (ASCII 13)             |
> | `\``t`   | Output a TAB character (ASCII 9)                          |
> | `\\`     | Output `\`                                                |
> | `\``%`   | Output `%`                                                |
> | `\``nnn` | Output character \<nnn\>. \<nnn\> is a 1-3 decimal number |

Escape sequences to include arguments from `Arguments`. Each time a %-escape sequence is found in `Format` the next argument from `Arguments` is formatted according to the specification.

> |  |  |
> |----|----|
> | `%t` | [print/1](termrw.html#print/1) the next item (mnemonic: term) |
> | `%w` | [write/1](termrw.html#write/1) the next item |
> | `%q` | [writeq/1](termrw.html#writeq/1) the next item |
> | `%d` | display/1 the next item |
> | `%n` | Put the next item as a character |
> | `%r` | Write the next item N times where N is the second item (an integer) |
> | `%s` | Write the next item as a String (so it must be a list of characters) |
> | `%f` | Perform a [ttyflush/0](chario.html#ttyflush/0) (no items used) |
> | `%Nc` | Write the next item Centered in N columns. |
> | `%Nl` | Write the next item Left justified in N columns. |
> | `%Nr` | Write the next item Right justified in N columns. |

deprecated  
New code should use [format/1](format.html#format/1), [format/2](format.html#format/2), etc.

\[det\]**swritef**(`-String, +Format`)  
\[det\]**swritef**(`-String, +Format, +Arguments`)  
Use [writef/1](writef.html#writef/1) or [writef/2](writef.html#writef/2) and write the result to a *string*. Note that this is a string in the sense of [string_codes/2](string.html#string_codes/2), *not* a list of `character(-code)`s.

deprecated  
. See [format/2](format.html#format/2),3 and/or [with_output_to/2](IO.html#with_output_to/2).
