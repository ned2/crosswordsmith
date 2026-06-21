
## 4.13 DCG Grammar rules

Grammar rules form a comfortable interface to *difference lists*. They are designed both to support writing parsers that build a parse tree from a list of characters or tokens and for generating a flat list from a term.

Grammar rules look like ordinary clauses using `-->``/2` for separating the head and body rather than `:-``/2`. Expanding grammar rules is done by [expand_term/2](consulting.html#expand_term/2), which adds two additional arguments to each term for representing the difference list.

The body of a grammar rule can contain three types of terms. A callable term is interpreted as a reference to a grammar rule. Code between `{`...`}` is interpreted as plain Prolog code, and finally, a list is interpreted as a sequence of *literals*. The Prolog control-constructs (`\+``/1`, `->``/2`, `;``//`2, `,``/2` and `!``/0`) can be used in grammar rules.

We illustrate the behaviour by defining a rule set for parsing an integer.

``` code
integer(I) -->
        digit(D0),
        digits(D),
        { number_codes(I, [D0|D])
        }.

digits([D|T]) -->
        digit(D), !,
        digits(T).
digits([]) -->
        [].

digit(D) -->
        [D],
        { code_type(D, digit)
        }.
```

Grammar rule sets are called using the built-in predicates [phrase/2](DCG.html#phrase/2) and [phrase/3](DCG.html#phrase/3):

**phrase**(`:DCGBody, ?List`)  
Equivalent to `phrase(``DCGBody``, ``InputList``, [])`.

**phrase**(`:DCGBody, ?List, ?Rest`)  
True when `DCGBody` applies to the difference `List`/`Rest`. Although `DCGBody` is typically a *callable* term that denotes a grammar rule, it can be any term that is valid as the body of a DCG rule.

The example below calls the rule set integer//1 defined in [section 4.13](DCG.html#sec:4.13) and available from `library(library(dcg/basics))`, binding `Rest` to the remainder of the input after matching the integer.

``` code
?- [library(dcg/basics)].
?- atom_codes('42 times', Codes),
   phrase(integer(X), Codes, Rest).
X = 42
Rest = [32, 116, 105, 109, 101, 115]
```

The next example exploits a complete body. Given the following definition of digit_weight//1 , we can pose the query below.

``` code
digit_weight(W) -->
        [D],
        { code_type(D, digit(W)) }.
```

``` code
?- atom_codes('Version 3.4', Codes),
   phrase(("Version ",
           digit_weight(Major),".",digit_weight(Minor)),
          Codes).
Major = 3,
Minor = 4.
```

The SWI-Prolog implementation of [phrase/3](DCG.html#phrase/3) verifies that the `List` and `Rest` arguments are unbound, bound to the empty list or a list *cons cell*. Other values raise a type error.^(84The ISO standard allows for both raising a type error and accepting any term as input and output. Note the tail of the list is not checked for performance reasons.) The predicate [call_dcg/3](DCG.html#call_dcg/3) is provided to use grammar rules with terms that are not lists.

Note that the syntax for lists of codes changed in SWI-Prolog version 7 (see [section 5.2](string.html#sec:5.2)). If a DCG body is translated, both `"text"` and `` `text` `` is a valid code-list literal in version 7. A version 7 string (`"text"`) is **not** acceptable for the second and third arguments of [phrase/3](DCG.html#phrase/3). This is typically not a problem for applications as the input of a DCG rarely appears in the source code. For testing in the toplevel, one must use double quoted text in versions prior to 7 and back quoted text in version 7 or later.

See also [portray_text/1](portraytext.html#portray_text/1), which can be used to print lists of character codes as a string to the top level and debugger to facilitate debugging DCGs that process character codes. The library `library(apply_macros)` compiles [phrase/3](DCG.html#phrase/3) if the argument is sufficiently instantiated, eliminating the runtime overhead of translating `DCGBody` and meta-calling.

**call_dcg**(`:DCGBody, ?State0, ?State`)  
As [phrase/3](DCG.html#phrase/3), but without type checking `State0` and `State`. This allows for using DCG rules for threading an arbitrary state variable. This predicate was introduced after type checking was added to [phrase/3](DCG.html#phrase/3).^(85After discussion with Samer Abdallah.)

A portable solution for threading state through a DCG can be implemented by wrapping the state in a list and use the DCG semicontext facility. Subsequently, the following predicates may be used to access and modify the state:^(86This solution was proposed by Markus Triska.)

``` code
state(S), [S] --> [S].
state(S0, S), [S] --> [S0].
```

As stated above, grammar rules are a general interface to difference lists. To illustrate, we show a DCG-based implementation of [reverse/2](lists.html#reverse/2):

``` code
reverse(List, Reversed) :-
        phrase(reverse(List), Reversed).

reverse([])    --> [].
reverse([H|T]) --> reverse(T), [H].
```
