
## 4.22 Analysing and Constructing Atoms

These predicates convert between certain Prolog atomic values on one hand and lists of *character codes* (or, for [atom_chars/2](manipatom.html#atom_chars/2), *characters*) on the other. The Prolog atomic values can be atoms, *character*s (which are atoms of length 1), SWI-Prolog strings, as well as numbers (integers, floats and non-integer rationals).

The *character codes*, also known as *code values*, are integers. In SWI-Prolog, these integers are Unicode code points.^(bugOn Windows the range is limited to UCS-2, 0..65535.)

To ease the pain of all text representation variations in the Prolog community, all SWI-Prolog predicates behave as *flexible as possible*. This implies the‘list-side’accepts both a character-code-list and a character-list and the‘atom-side’accepts all atomic types (atom, number and string). For example, the predicates [atom_codes/2](manipatom.html#atom_codes/2), [number_codes/2](manipatom.html#number_codes/2) and [name/2](manipatom.html#name/2) behave the same in mode (+,-), i.e.,‘listwards’, from a constant to a list of character codes. When converting the other way around:

- [atom_codes/2](manipatom.html#atom_codes/2) will generate an atom;
- [number_codes/2](manipatom.html#number_codes/2) will generate a number or throw an exception;
- [name/2](manipatom.html#name/2) will generate a number if possible and an atom otherwise.

\[ISO\]**atom_codes**(`?Atom, ?CodeList`)  
Convert between an atom and a list of *character codes* (integers denoting characters).

- If `Atom` is instantiated, it will be translated into a list of character codes, which are unified with `CodeList`.
- If `Atom` is uninstantiated and `CodeList` is a list of character codes, then `Atom` will be unified with an atom constructed from this list.

``` code
?- atom_codes(hello, X).
X = [104, 101, 108, 108, 111].
```

The‘listwards’call to [atom_codes/2](manipatom.html#atom_codes/2) can also be written (functionally) using backquotes instead:

``` code
?- Cs = `hello`.
Cs = [104, 101, 108, 108, 111].
```

Backquoted strings can be mostly found in the body of DCG rules that process lists of character codes.

Note that this is the default interpretation for backquotes. It can be changed on a per-module basis by setting the value of the Prolog flag [back_quotes](flags.html#flag:back_quotes).

\[ISO\]**atom_chars**(`?Atom, ?CharList`)  
Similar to [atom_codes/2](manipatom.html#atom_codes/2), but `CharList` is a list of *character*s (atoms of length 1) rather than a list of *character codes* (integers denoting characters).

``` code
?- atom_chars(hello, X).
X = [h, e, l, l, o]
```

\[ISO\]**char_code**(`?Atom, ?Code`)  
Convert between a single *character* (an atom of length 1), and its *character code* (an integer denoting the corresponding character). The predicate alternatively accepts an SWI-Prolog string of length 1 at `Atom` place.

\[ISO\]**number_chars**(`?Number, ?CharList`)  
Similar to [atom_chars/2](manipatom.html#atom_chars/2), but converts between a number and its representation as a list of *characters* (atoms of length 1).

- If `CharList` is a *proper list*, i.e., not unbound or a *partial list*, `CharList` is parsed according to the Prolog syntax for numbers and the resulting number is unified with `Number`. A `syntax_error` exception is raised if `CharList` is instantiated to a ground, proper list but does not represent a valid Prolog number.
- Otherwise, if `Number` is indeed a number, `Number` is serialized and the result is unified with `CharList`.

Following the ISO standard, the Prolog syntax for number allows for *leading* white space (including newlines) and does not allow for *trailing* white space.^(119ISO also allows for Prolog comments in leading white space. We--and most other implementations--believe this is incorrect. We also believe it would have been better not to allow for white space, or to allow for both leading and trailing white space.)

Prolog syntax-based conversion can also be achieved using [format/3](format.html#format/3) and [read_from_chars/2](charsio.html#read_from_chars/2).

\[ISO\]**number_codes**(`?Number, ?CodeList`)  
As [number_chars/2](manipatom.html#number_chars/2), but converts to a list of character codes rather than characters. In the mode (-,+), both predicates behave identically to improve handling of non-ISO source.

**atom_number**(`?Atom, ?Number`)  
Realises the popular combination of [atom_codes/2](manipatom.html#atom_codes/2) and [number_codes/2](manipatom.html#number_codes/2) to convert between atom and number (integer, float or non-integer rational) in one predicate, avoiding the intermediate list. Unlike the ISO standard [number_codes/2](manipatom.html#number_codes/2) predicates, [atom_number/2](manipatom.html#atom_number/2) fails silently in mode (+,-) if `Atom` does not represent a number.

**name**(`?Atomic, ?CodeList`)  
`CodeList` is a list of character codes representing the same text as `Atomic`. Each of the arguments may be a variable, but not both.

- When `CodeList` describes an integer or floating point number and `Atomic` is a variable, `Atomic` will be unified with the numeric value described by `CodeList` (e.g., `name(N, "300"), 400 is N + 100` succeeds).
- If `CodeList` is not a representation of a number, `Atomic` will be unified with the atom with the name given by the character code list.
- If `Atomic` is an atom or number, the unquoted print representation of it as a character code list is unified with `CodeList`.

This predicate is part of the Edinburgh tradition. It should be considered *deprecated* although, given its long tradition, it is unlikely to be removed from the system. It still has some value for converting input to a number or an atom (depending on the syntax). New code should consider the ISO predicates [atom_codes/2](manipatom.html#atom_codes/2), [number_codes/2](manipatom.html#number_codes/2) or the SWI-Prolog predicate [atom_number/2](manipatom.html#atom_number/2).

**term_to_atom**(`?Term, ?Atom`)  
True if `Atom` describes a term that unifies with `Term`. When `Atom` is instantiated, `Atom` is parsed and the result unified with `Term`. If `Atom` has no valid syntax, a `syntax_error` exception is raised. Otherwise `Term` is “written” on `Atom` using [write_term/2](termrw.html#write_term/2) with the option `quoted(true)`. See also [format/3](format.html#format/3), [with_output_to/2](IO.html#with_output_to/2) and [term_string/2](string.html#term_string/2).

\[deprecated\]**atom_to_term**(`+Atom, -Term, -Bindings`)  
Use `Atom` as input to [read_term/2](termrw.html#read_term/2) using the option `variable_names` and return the read term in `Term` and the variable bindings in `Bindings`. `Bindings` is a list of `Name`` = ``Var` couples, thus providing access to the actual variable names. See also [read_term/2](termrw.html#read_term/2). If `Atom` has no valid syntax, a `syntax_error` exception is raised. If `Atom` only contains white space and/or comments, an `syntax_error(end_of_string)` exception is raised. New code should use [read_term_from_atom/3](termrw.html#read_term_from_atom/3).

\[ISO\]**atom_concat**(`?Atom1, ?Atom2, ?Atom3`)  
`Atom3` forms the concatenation of `Atom1` and `Atom2`. At least two of the arguments must be instantiated to atoms. This predicate also allows for the mode (-,-,+), non-deterministically splitting the 3rd argument into two parts (as [append/3](lists.html#append/3) does for lists). SWI-Prolog allows for atomic arguments. Portable code must use [atomic_concat/3](manipatom.html#atomic_concat/3) if non-atom arguments are involved.

**atomic_concat**(`+Atomic1, +Atomic2, -Atom`)  
`Atom` represents the text after converting `Atomic1` and `Atomic2` to text and concatenating the result:

``` code
?- atomic_concat(name, 42, X).
X = name42.
```

\[commons\]**atomic_list_concat**(`+List, -Atom`)  
`List` is a list of strings, atoms, integers, floating point numbers or non-integer rationals. Succeeds if `Atom` can be unified with the concatenated elements of `List`. Equivalent to `atomic_list_concat(List,’’, Atom)`.

\[commons\]**atomic_list_concat**(`+List, +Separator, -Atom`)  
Creates an atom just like [atomic_list_concat/2](manipatom.html#atomic_list_concat/2), but inserts `Separator` between each pair of inputs. For example:

``` code
?- atomic_list_concat([gnu, gnat], ', ', A).

A = 'gnu, gnat'
```

The‘atomwards\` transformation is usually called a *string join* operation in other programming languages.

The SWI-Prolog version of this predicate can also be used to split atoms by instantiating `Separator` and `Atom` as shown below. We kept this functionality to simplify porting old SWI-Prolog code where this predicate was called concat_atom/3. When used in mode (-,+,+), `Separator` must be a non-empty atom. See also [split_string/4](string.html#split_string/4).

``` code
?- atomic_list_concat(L, -, 'gnu-gnat').

L = [gnu, gnat]
```

\[ISO\]**atom_length**(`+Atom, -Length`)  
True if `Atom` is an atom of `Length` characters. The SWI-Prolog version accepts all atomic types, as well as code-lists and character-lists. New code should avoid this feature and use [write_length/3](termrw.html#write_length/3) to get the number of characters that would be written if the argument was handed to [write_term/3](termrw.html#write_term/3).

\[deprecated\]**atom_prefix**(`+Atom, +Prefix`)  
True if `Atom` starts with the characters from `Prefix`. Its behaviour is equivalent to `?- sub_atom(``Atom``, 0, _, _, ``Prefix``)`. Deprecated.

\[ISO\]**sub_atom**(`+Atom, ?Before, ?Length, ?After, ?SubAtom`)  
ISO predicate for breaking atoms. It maintains the following relation: `SubAtom` is a sub-atom of `Atom` that starts at (0-based index) `Before`, has `Length` characters, and `Atom` contains `After` characters after the match. The implementation minimises non-determinism and creation of atoms. This is a flexible predicate that can do search, prefix- and suffix-matching, etc. Scenarios that use this predicate often generate atoms that with a short lifetime; in such cases [sub_string/5](string.html#sub_string/5) may be a better alternative. Examples:

Pick out a sub-atom of length 3 starting a 0-based index 2:

``` code
?- sub_atom(aaxyzbbb, 2, 3, After, SubAtom).
After = 3,
SubAtom = xyz.
```

The following example splits a string of the form \<`name`\>=\<`value`\> into the name part (an atom) and the value (a string).

``` code
name_value(String, Name, Value) :-
    sub_atom(String, Before, _, After, "="),
    !,
    sub_atom(String, 0, Before, _, Name),
    sub_atom(String, _, After, 0, Value).
```

This example defines a predicate that inserts a value at a position. Note that [sub_string/5](string.html#sub_string/5) is used here instead of [sub_atom/5](manipatom.html#sub_atom/5) to avoid the overhead of creating atoms for the intermediate results.

``` code
atom_insert(Str, Val, At, NewStr) :-
    sub_string(Str, 0, At, A1, S1),
    sub_string(Str, At, A1, _, S2),
    atomic_list_concat([S1,Val,S2], NewStr).
```

On backtracking, matches are delivered in order left-to-right (i.e. `Before` increases monotonically):

``` code
?- sub_atom('xATGATGAxATGAxATGAx', Before, Length, After, 'ATGA').
Before = 1, Length = 4, After = 14 ;
Before = Length, Length = 4, After = 11 ;
Before = 9, Length = 4, After = 6 ;
Before = 14, Length = 4, After = 1 ;
false.
```

See also [sub_string/5](string.html#sub_string/5), the corresponding predicate for SWI-Prolog strings.

\[semidet\]**sub_atom_icasechk**(`+Haystack, ?Start, +Needle`)  
True when `Needle` is a sub atom of `Haystack` starting at `Start`. The match is‘half case insensitive’, i.e., uppercase letters in `Needle` only match themselves, while lowercase letters in `Needle` match case insensitively. `Start` is the first 0-based offset inside `Haystack` where `Needle` matches.^(120This predicate replaces \$apropos_match/2, used by the help system, while extending it with locating the (first) match and performing case insensitive prefix matching. We are still not happy with the name and interface.)
