
# 5 SWI-Prolog extensions

This chapter describes extensions to the Prolog language introduced with SWI-Prolog version 7 in 2014. The changes bring more modern syntactical conventions to Prolog such as key-value maps, called *dicts*, as primary citizens and a restricted form of *functional notation*. They also extend Prolog basic types with strings, providing a natural notation to textual material as opposed to identifiers (atoms) and lists.

These extensions make the syntax more intuitive to new users, simplify the integration of domain specific languages (DSLs) and facilitate a more natural Prolog representation for popular exchange languages such as XML and JSON.

While many programs run unmodified in SWI-Prolog version 7, some require modifications, especially those that pass double quoted strings to general purpose list processing predicates. See [section 5.2.4](string.html#sec:5.2.4) and [section 5.2.5](string.html#sec:5.2.5) for information and tools on porting. We provide a tool ([list_strings/0](string.html#list_strings/0)) that we used to port a huge code base in half a day.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[5.1 Lists are special](ext-lists.html)

[5.1.1 Motivating‘`[|]`’and `[]` for lists](ext-lists.html#sec:5.1.1)

[5.2 The string type and its double quoted syntax](string.html)

[5.2.1 Representing text: strings, atoms and code lists](string.html#sec:5.2.1)

[5.2.2 Predicates that operate on strings](string.html#sec:5.2.2)

[5.2.3 Why has the representation of double quoted text changed?](string.html#sec:5.2.3)

[5.2.4 Adapting code for double quoted strings](string.html#sec:5.2.4)

[5.2.5 Predicates to support adapting code for double quoted strings](string.html#sec:5.2.5)

[5.3 Syntax changes since SWI-Prolog 7](ext-syntax.html)

[5.3.1 Operators and quoted atoms](ext-syntax.html#sec:5.3.1)

[5.3.2 Compound terms with zero arguments](ext-syntax.html#sec:5.3.2)

[5.3.3 Block operators](ext-syntax.html#sec:5.3.3)

[5.4 Dicts: structures with named arguments](bidicts.html)

[5.4.1 Dict types and tags](bidicts.html#sec:5.4.1)

[5.4.2 Functions on dicts](bidicts.html#sec:5.4.2)

[5.4.2.1 User defined functions on dicts](bidicts.html#sec:5.4.2.1)

[5.4.2.2 Predefined functions on dicts](bidicts.html#sec:5.4.2.2)

[5.4.3 Predicates for managing dicts](bidicts.html#sec:5.4.3)

[5.4.3.1 Destructive assignment in dicts](bidicts.html#sec:5.4.3.1)

[5.4.4 When to use dicts?](bidicts.html#sec:5.4.4)

[5.4.5 A motivation for dicts as primary citizens](bidicts.html#sec:5.4.5)

[5.4.6 Implementation notes about dicts](bidicts.html#sec:5.4.6)

[5.5 Integration of strings and dicts in the libraries](ext-integration.html)

[5.5.1 Dicts and option processing](ext-integration.html#sec:5.5.1)

[5.5.2 Dicts in core data structures](ext-integration.html#sec:5.5.2)

[5.5.3 Dicts, strings and XML](ext-integration.html#sec:5.5.3)

[5.5.4 Dicts, strings and JSON](ext-integration.html#sec:5.5.4)

[5.5.5 Dicts, strings and HTTP](ext-integration.html#sec:5.5.5)

[5.6 Single Sided Unification rules](ssu.html)

[5.6.1 Single Sided Unification Guards](ssu.html#sec:5.6.1)

[5.6.2 Consequences of `=>` single sided unification rules](ssu.html#sec:5.6.2)

[5.6.3 Single sided unification for Definite Clause Grammars](ssu.html#sec:5.6.3)

[5.6.4 SSU: Future considerations](ssu.html#sec:5.6.4)

[5.7 Remaining issues](ext-issues.html)
