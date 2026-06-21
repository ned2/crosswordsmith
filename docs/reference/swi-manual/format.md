
## 4.32 Formatted Write

The format family of predicates is the most versatile and portable way to produce textual output. These predicates where introduced by Quintus Prolog. The SWI-Prolog implementation is compatible with [PIP-0110](https://prolog-lang.org/ImplementersForum/0110-format.html). SWI-Prolog implements a couple of extensions:

- Support for `~@` for capturing the output of any predicate. This extension is implemented by several Prolog systems.
- Support for `~h` and `~H` for more control over floating point number output. This is compatible with SICStus.
- Support the `:` modifier, notably used to switch between *locale* specific output and Prolog syntax.
- The predicate [format/3](format.html#format/3) allows for an output specification compatible to [with_output_to/2](IO.html#with_output_to/2).

**format**(`+Format`)  
**format**(`+Format, :Arguments`)  
**format**(`+Output, +Format, :Arguments`)  
Write formatted output to `Output` or `current_output`. `Format` is an atom, list of character codes, or a Prolog string. `Arguments` is a list of arguments required by the format specification. For backward compatibility, if `Format` needs exactly one argument and the required argument is not a list, single argument needs not be nested in a list. This feature is deprecated as it easily leads to mistakes and make static analysis by [check/0](check.html#check/0) less accurate.

`Output` is normally a stream. SWI-Prolog's [format/3](format.html#format/3) also accepts any other output accepted by [with_output_to/2](IO.html#with_output_to/2). For example:

``` code
?- format(atom(A), '~D', [1000000]).
A = '1,000,000'
```

Format control sequences start with the tilde (`~`), followed by an optional numeric argument, optionally followed by a colon modifier (:), ^(144The colon modifiers is a SWI-Prolog extension, proposed by Richard O'Keefe.) followed by a character describing the action to be undertaken. A numeric argument is either a sequence of digits, representing a positive decimal number, a sequence `‘<``character``>`, representing the character code value of the character (only useful for `~t`) or an asterisk (`*`), in which case the numeric argument is taken from the next argument of the argument list, which should be a positive integer or an arithmetic expression that evaluates to a possitive integer. The following four examples all pass 46 (`.`) to `~t` and use a tab at 72.

``` code
?- format('~w ~46t ~w~72|~n', ['Title', 'Page']).
?- format('~w ~`.t ~w~72|~n', ['Title', 'Page']).
?- format('~w ~*t ~w~72|~n', ['Title', 46, 'Page']).
?- format('~w ~*t ~w~*|~n', ['Title', 46, 'Page', 80-8]).
```

Some format expressions may call back Prolog, i.e., `~p`, `~W`, `~@` and user defined extensions registered with [format_predicate/2](format.html#format_predicate/2). Output written to the stream `current_output` is merged into the [format/2](format.html#format/2) output. If there is no pending *rubber* (`~t`) and the the position notation aligns, only the output is switched. Otherwise the output is captured in a temporary memory buffer and emitted after the callback finishes. The system attempts to preserve the position and alignment promises. It sets the `tty` property of the temporary stream to reflect the main stream and uses the position information of the temporary stream to update its notion of the position. Notable [ansi_format/3](ansiterm.html#ansi_format/3) cooperates properly in callbacks.^(145As of version 8.3.30).

Numeric conversion (`d`, `D`, `e`, `E`, `f`, `F`, `g`, `G`, `h` and `H`) accept an arithmetic expression as argument. This is introduced to handle rational numbers transparently (see [section 4.27.2.2](arith.html#sec:4.27.2.2)). The floating point conversions allow for unlimited precision for printing rational numbers in decimal form. E.g., the following will write as many 3's as you want by changing the‘50’.

``` code
?- format('~50f', [10r3]).
3.33333333333333333333333333333333333333333333333333
```

- `~`  
  Output the tilde itself.

- `a`  
  Output the next argument, which must be an atom. This option is equivalent to **w**, except that it requires the argument to be an atom. For flexibility, SWI-Prolog also accepts a packed string (see [section 5.2](string.html#sec:5.2)) as valid argument.

- `c`  
  Interpret the next argument as a character code and add it to the output. This argument must be a valid Unicode character code. Note that the actually emitted bytes are defined by the character encoding of the output stream and an exception may be raised if the output stream is not capable of representing the requested Unicode character. See [section 2.18.1](widechars.html#sec:2.18.1) for details.

- `d`  
  Output next argument as a decimal number. It should be an integer. If a numeric argument is specified, a dot is inserted `argument` positions from the right (useful for doing fixed point arithmetic with integers, such as handling amounts of money).

  The colon modifier (e.g., `~:d`) causes the number to be printed according to the locale of the output stream. See [section 4.23](locale.html#sec:4.23).

- `D`  
  Same as **d**, but makes large values easier to read by inserting a comma every three digits left or right of the dot. This is the same as `~:d`, but using the fixed English locale. If **D** is modified using the colon (`~:D`), it uses the Prolog grouping character `_`. See also setup_prolog_integer_grouping/0. Future versions may also support line breaks in big integers.

- `e`  
  Output next argument as a floating point number in exponential notation. The numeric argument specifies the precision. Default is 6 digits. Exact representation depends on the C library function **printf()**. This function is invoked with the format `%.<``precision``>e`.

- `E`  
  Equivalent to **e**, but outputs a capital E to indicate the exponent.

- `f`  
  Floating point in non-exponential notation. The numeric argument defines the number of digits right of the decimal point. If the numeric argument is zero (0), the value is printed as an integer. If the colon modifier (:) is used, the float is formatted using conventions from the current locale, which may define the decimal point as well as grouping of digits left of the decimal point.

- `F`  
  As **f**, but used uppercase letters for non-normal floats, e.g.

  ``` code
  ?- format('~F', [nan]).
  NAN
  ```

- `g`  
  Floating point in **e** or **f** notation, whichever is shorter.

- `G`  
  Floating point in **E** or **f** notation, whichever is shorter.

- `h`  

- `H`  
  Print a floating point number with the minimal number of digits such that [read/1](termrw.html#read/1) produces exactly (as in [==/2](compare.html#==/2)) the same number. The argument specifies whether a number is written using exponential notation (using **e** (h) or **E** (H)) or fixed point notation (as `~f`). If the argument is -1, the number is always written using exponential notation. Otherwise, number is written using exponential notation if the exponent is less than `Arg`-1 or greater than `Arg`+`d`, where `d` is the number of digits emitted to establish the required precision. Using an argument larger than the the maximum exponent such as `~999h` never uses exponential notation. The default argument is 3. The predicate [write/1](termrw.html#write/1) and friends act as if using the format `~3h`. This option is compatible to SICStus Prolog.

- `i`  
  Ignore next argument of the argument list. Produces no output.

- `I`  
  Emit a decimal number using Prolog digit grouping (the underscore, `_`). The argument describes the size of each digit group. The default is 3. See also [section 2.15.1.5](syntax.html#sec:2.15.1.5). For example:

  ``` code
  ?- A is 1<<100, format('~10I', [A]).
  1_2676506002_2822940149_6703205376
  ```

- `k`  
  Give the next argument to [write_canonical/1](termrw.html#write_canonical/1).

- `n`  
  Output a newline character.

- `N`  
  Only output a newline if the last character output on this stream was not a newline. Not properly implemented yet.

- `p`  
  Give the next argument to [print/1](termrw.html#print/1).

- `q`  
  Give the next argument to [writeq/1](termrw.html#writeq/1).

- `r`  
  Print integer in radix numeric argument notation (default 8). Thus `~16r` prints its argument hexadecimal. The argument should be in the range `[2, ... , 36]`. Lowercase letters are used for digits above 9. The colon modifier may be used to form locale-specific digit groups.

- `R`  
  Same as **r**, but uses uppercase letters for digits above 9.

- `s`  
  Output text from a list of character codes, characters, string (see [string/1](typetest.html#string/1) and [section 5.2](string.html#sec:5.2)) or atom from the next argument. If an numeric argument the specified number of characters is emitted. This implies the text is truncated if too long or right-padded with spaces if it is too short. Note that this way of padding is not in the spirit of [format/2](format.html#format/2). Using `~t` for padding is more flexible and in line with the spirit of [format/2](format.html#format/2).

- `@`  
  Interpret the next argument as a goal and execute it. Output written to the `current_output` stream is inserted at this place. Goal is called in the module calling [format/3](format.html#format/3). This option is not present in the original definition by Quintus, but supported by some other Prolog systems. The goal is executed as `\+ \+ Goal`, i.e., bindings created by the goal are discarded.

- `t`  
  All remaining space between two tab stops is distributed equally over `~t` statements between the tab stops. This space is padded with spaces by default. If an argument is supplied, it is taken to be the character code of the character used for padding. This can be used to do left or right alignment, centering, distributing, etc. See also `~|` and `~+` to set tab stops. A tab stop is assumed at the start of each line. If the space between the two tab stops cannot be divided by the number of `~t` statements, one additional padding character is added to each `~t` statement until the tab stop is reached. Additional padding is added from the right.^(146Until 9.3.33, additional padding was added from the middle. Current behaviour is defined by PIP-0110.).

- `|`  
  Set a tab stop on the current position. If an argument is supplied set a tab stop on the position of that argument. This will cause all `~t`’s to be distributed between the previous and this tab stop.

  If the current column is at or past the requested tabstop and the modifier (:) is used, a newline is inserted and the padding character of the last `~t` is used to pad to the requested position.

- `+`  
  Set a tab stop (as `~|`) relative to the last tab stop or the beginning of the line if no tab stops are set before the `~+`. This constructs can be used to fill fields. The partial format sequence below prints an integer right-aligned and padded with zeros in 6 columns. The ... sequences in the example illustrate that the integer is aligned in 6 columns regardless of the remainder of the format specification.

  ``` code
          format('...~|~`0t~d~6+...', [..., Integer, ...])
  ```

- `w`  
  Give the next argument to [write/1](termrw.html#write/1).

- `W`  
  Give the next two arguments to [write_term/2](termrw.html#write_term/2). For example, `format('~W', [Term, [numbervars(true)]])`. This option is SWI-Prolog specific.

Example:

``` code
simple_statistics :-
    <obtain statistics>         % left to the user
    format('~tStatistics~t~72|~n~n'),
    format('Runtime: ~`.t ~2f~34|  Inferences: ~`.t ~D~72|~n',
                                            [RunT, Inf]),
    ....
```

will output

``` code
                             Statistics

Runtime: .................. 3.45  Inferences: .......... 60,345
```

### 4.32.1 Programming Format

**format_predicate**(`+Char, +Head`)  
If a sequence `~c` (tilde, followed by some character) is found, the [format/3](format.html#format/3) and friends first check whether the user has defined a predicate to handle the format. If not, the built-in formatting rules described above are used. `Char` is either a character code or a one-character atom, specifying the letter to be (re)defined. `Head` is a term, whose name and arity are used to determine the predicate to call for the redefined formatting character. The first argument to the predicate is the numeric argument of the format command, or the atom `default` if no argument is specified. The remaining arguments are filled from the argument list. The example below defines `~T` to print a timestamp in ISO8601 format (see [format_time/3](system.html#format_time/3)). The subsequent block illustrates a possible call.

``` code
:- format_predicate('T', format_time(_Arg,_Time)).

format_time(_Arg, Stamp) :-
        must_be(number, Stamp),
        format_time(current_output, '%FT%T%z', Stamp).
```

``` code
?- get_time(Now),
   format('Now, it is ~T~n', [Now]).
Now, it is 2012-06-04T19:02:01+0200
Now = 1338829321.6620328.
```

**current_format_predicate**(`?Code, ?:Head`)  
True when `~``Code` is handled by the user-defined predicate specified by `Head`.
