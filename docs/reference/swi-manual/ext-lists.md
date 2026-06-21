
## 5.1 Lists are special

As of version 7, SWI-Prolog lists can be distinguished unambiguously at runtime from `.``/2` terms and the atom `’[]’`.

``` code
   Traditional list               SWI-Prolog 7 list

       '.'                              '[|]'
      /   \                             /   \
     1    '.'                          1   '[|]'
         /   \                             /   \
        2    '.'                          2   '[|]'
            /   \                             /   \
           3   '[]'                          3     []

           terminated with                   terminated with
           the atom '[]',                    a special constant
           indistinguishable from text       which is printed as []
```

The constant `[]` is special constant that is not an atom. It has the following properties:

``` code
atom([]).        fails
atomic([]).      succeeds
[] == '[]'.      fails
[] == [].        succeeds
```

The‘cons’operator for creating *list cells* has changed from the pretty atom‘`.`’to the ugly atom‘`[|]`’, so we can use the‘`.`’for other purposes, notably functional notation on *dicts*. See [section 5.4.2](bidicts.html#sec:5.4.2).

This modification has minimal impact on typical Prolog code. It does affect foreign code (see [section 12](foreign.html#sec:12)) that uses the normal atom and compound term interface for manipulating lists. In most cases this can be avoided by using the dedicated list functions. For convenience, the macros `ATOM_nil` and `ATOM_dot` are provided by `SWI-Prolog.h`.

Another place that is affected is [write_canonical/1](termrw.html#write_canonical/1). Impact is minimized by using the list syntax for lists. The predicates [read_term/2](termrw.html#read_term/2) and [write_term/2](termrw.html#write_term/2) support the option `dotlists(true)`, which causes [read_term/2](termrw.html#read_term/2) to read `.(a,[])` as `[a]` and [write_term/2](termrw.html#write_term/2) to write `[a]` as `.(a,[])`.

### 5.1.1 Motivating‘`[|]`’and `[]` for lists

Representing lists the conventional way using `.``/2` as list cell and the atom `'[]'` as list terminator both (independently) pose conflicts, while these conflicts are easily avoided.

- Using `.``/2` prevents using this commonly used symbol as an operator because `a.B` cannot be distinguished from `[a|B]`. Freeing `.``/2` provides us with a unique term that we can use for functional notation on dicts as described in [section 5.4.2](bidicts.html#sec:5.4.2).

- Using the atom `'[]'` as list terminator prevents dynamic distinction between atoms and the empty list. As a result, we cannot use type polymorphism that involve both atoms and lists. For example, we cannot use *multi lists* (arbitrary deeply nested lists) of atoms. Multi lists of atoms are in some situations a good representation of a flat list that is assembled from sub sequences. The alternative, using difference lists or DCGs, is often less natural and sometimes requires‘opening’proper lists (i.e., copying the list while replacing the terminating atom `'[]'` with a variable) that have to be added to the sequence. The ambiguity of atom and list is particularly painful when mapping external data representations that do not suffer from this ambiguity.

  At the same time, avoiding atom `'[]'` as a list terminator makes the various text representations unambiguous, which allows us to write predicates that require a textual argument to accept any of atoms, strings, lists of character codes or characters. Traditionally, the empty list, as an atom, is afflicted with an ambiguous interpretation as it can stand for any of the strings `"[]"` and `""`.
