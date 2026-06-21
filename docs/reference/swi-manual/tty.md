
## 4.34 Terminal Control

The following predicates form a simple access mechanism to the Unix termcap library to provide terminal-independent I/O for screen terminals. These predicates are only available on Unix machines. The SWI-Prolog Windows console accepts the ANSI escape sequences.

**tty_get_capability**(`+Name, +Type, -Result`)  
Get the capability named `Name` from the termcap library. See termcap(5) for the capability names. `Type` specifies the type of the expected result, and is one of `string`, `number` or `bool`. String results are returned as an atom, number results as an integer, and bool results as the atom `on` or `off`. If an option cannot be found, this predicate fails silently. The results are only computed once. Successive queries on the same capability are fast. This predicate can raise several exceptions if the terminal environment is incomplete, notably if the environment variable `TERM` does not exist or there is no matching entry in the *termcap database*.

**tty_goto**(`+X, +Y`)  
Goto position (`X`, `Y`) on the screen. Note that the predicates [line_count/2](streamstat.html#line_count/2) and [line_position/2](streamstat.html#line_position/2) will not have a well-defined behaviour while using this predicate.

**tty_put**(`+Atom, +Lines`)  
Put an atom via the termcap library function **tputs()**. This function decodes padding information in the strings returned by [tty_get_capability/3](tty.html#tty_get_capability/3) and should be used to output these strings. `Lines` is the number of lines affected by the operation, or 1 if not applicable (as in almost all cases).

**tty_size**(`-Rows, -Columns`)  
Determine the size of the terminal. This predicate is provided for POSIX terminals, Windows console and the **swipl-win** Prolog consoles. Portable code should wrap this into [catch/3](exception.html#catch/3) and use some default in case the console window size cannot be determined.
