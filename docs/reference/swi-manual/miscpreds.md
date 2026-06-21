
## 4.45 Miscellaneous

**dwim_match**(`+Atom1, +Atom2`)  
True if `Atom1` matches `Atom2` in the‘Do What I Mean’sense. Both `Atom1` and `Atom2` may also be integers or floats. The two atoms match if:

- They are identical
- They differ by one character (spy `==` spu)
- One character is inserted/deleted (debug `==` deug)
- Two characters are transposed (trace `==` tarce)
- ‘Sub-words’are glued differently (existsfile `==` existsFile `==` exists_file)
- Two adjacent sub-words are transposed (existsFile `==` fileExists)

**dwim_match**(`+Atom1, +Atom2, -Difference`)  
Equivalent to [dwim_match/2](miscpreds.html#dwim_match/2), but unifies `Difference` with an atom identifying the difference between `Atom1` and `Atom2`. The return values are (in the same order as above): `equal`, `mismatched_char`, `inserted_char`, `transposed_char`, `separated` and `transposed_word`.

**wildcard_match**(`+Pattern, +String`)  
**wildcard_match**(`+Pattern, +String, +Options`)  
True if `String` matches the wildcard pattern `Pattern`. `Pattern` is very similar to the Unix `csh` pattern matcher. The patterns are given below:

|  |  |
|----|----|
| `?` | Matches one arbitrary character. |
| `*` | Matches any number of arbitrary characters. |
| `[ ...` | Matches one of the characters specified between the brackets. |
|  | `<``char1``>-<``char2``>` indicates a range. |
| `{...}` | Matches any of the patterns of the comma-separated list between the braces. |

Example:

``` code
?- wildcard_match('[a-z]*.{pro,pl}[%~]', 'a_hello.pl%').
true.
```

The [wildcard_match/3](miscpreds.html#wildcard_match/3) version processes the following option:

**case_sensitive**(`+Boolean`)  
When `false` (default `true`), match case insensitively.

**sleep**(`+Time`)  
Suspend execution `Time` seconds. `Time` is either a floating point number or an integer. Granularity is dependent on the system's timer granularity. A negative time causes the timer to return immediately. A zero time yields the CPU if this is supported on the target OS. On most non-realtime operating systems we can only ensure execution is suspended for **at least** `Time` seconds.

On Unix systems the [sleep/1](miscpreds.html#sleep/1) predicate is realised ---in order of preference--- by **nanosleep()**, **usleep()**, **select()** if the time is below 1 minute, or **sleep()**. On Windows systems **Sleep()** is used.
