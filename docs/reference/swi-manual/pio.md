
## A.37 library(pio): Pure I/O

This library provides pure list-based I/O processing for Prolog, where the communication to the actual I/O device is performed transparently through coroutining. This module itself is just an interface to the actual implementation modules.

### A.37.1 library(pure_input): Pure Input from files and streams

To be done  
Provide support for alternative input readers, e.g. reading terms, tokens, etc.

This module is part of `pio.pl`, dealing with *pure* *input*: processing input streams from the outside world using pure predicates, notably grammar rules (DCG). Using pure predicates makes non-deterministic processing of input much simpler.

Pure input uses attributed variables to read input from the external source into a list *on demand*. The overhead of lazy reading is more than compensated for by using block reads based on [read_pending_codes/3](chario.html#read_pending_codes/3).

Ulrich Neumerkel came up with the idea to use coroutining for creating a *lazy list*. His implementation repositioned the file to deal with re-reading that can be necessary on backtracking. The current implementation uses destructive assignment together with more low-level attribute handling to realise pure input on any (buffered) stream.

\[nondet\]**phrase_from_file**(`:Grammar, +File`)  
Process the content of `File` using the DCG rule `Grammar`. The space usage of this mechanism depends on the length of the not committed part of `Grammar`. Committed parts of the temporary list are reclaimed by the garbage collector, while the list is extended on demand due to unification of the attributed tail variable. Below is an example that counts the number of times a string appears in a file. The library dcg/basics provides [string//1](basics.html#string//1) matching an arbitrary string and [remainder//1](basics.html#remainder//1) which matches the remainder of the input without parsing.

``` code
:- use_module(library(dcg/basics)).

file_contains(File, Pattern) :-
        phrase_from_file(match(Pattern), File).

match(Pattern) -->
        string(_),
        string(Pattern),
        remainder(_).

match_count(File, Pattern, Count) :-
        aggregate_all(count, file_contains(File, Pattern), Count).
```

This can be called as (note that the pattern must be a string (code list)):

``` code
?- match_count('pure_input.pl', `file`, Count).
```

\[nondet\]**phrase_from_file**(`:Grammar, +File, +Options`)  
As [phrase_from_file/2](pio.html#phrase_from_file/2), providing additional `Options`. `Options` are passed to [open/4](IO.html#open/4).

**phrase_from_stream**(`:Grammar, +Stream`)  
Run Grammer against the character codes on `Stream`. `Stream` must be buffered.

**syntax_error**(`+Error`)`//`  
Throw the syntax error `Error` at the current location of the input. This predicate is designed to be called from the handler of [phrase_from_file/3](pio.html#phrase_from_file/3).

throws  
`error(syntax_error(Error), Location)`

\[det\]**lazy_list_location**(`-Location`)`//`  
Determine current (error) location in a lazy list. True when `Location` is an (error) location term that represents the current location in the DCG list.

|  |  |
|----|----|
| `Location` | is a term `file(Name, Line, LinePos, CharNo)` or `stream(Stream, Line, LinePos, CharNo)` if no file is associated to the stream RestLazyList. Finally, if the Lazy list is fully materialized (ends in `[]`), `Location` is unified with `end_of_file-CharCount`. |

See also  
[lazy_list_character_count//1](pio.html#lazy_list_character_count//1) only provides the character count.

**lazy_list_character_count**(`-CharCount`)`//`  
True when `CharCount` is the current character count in the Lazy list. The character count is computed by finding the distance to the next frozen tail of the lazy list. `CharCount` is one of:

- An integer
- A term end_of_file-Count

See also  
[lazy_list_location//1](pio.html#lazy_list_location//1) provides full details of the location for error reporting.

\[det\]**stream_to_lazy_list**(`+Stream, -List`)  
Create a lazy list representing the character codes in `Stream`. `List` is a partial list ending in an attributed variable. Unifying this variable reads the next block of data. The block is stored with the attribute value such that there is no need to re-read it.

Compatibility  
Unlike the previous version of this predicate this version does not require a repositionable stream. It does require a buffer size of at least the maximum number of bytes of a multi-byte sequence (6).
