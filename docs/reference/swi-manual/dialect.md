
# C Compatibility with other Prolog dialects

This chapter explains issues for writing portable Prolog programs. It was started after discussion with Vitor Santos Costa, the leading developer of YAP Prolog^(256[http://yap.sourceforge.net/](http://yap.sourceforge.net/)) YAP and SWI-Prolog have expressed the ambition to enhance the portability beyond the trivial Prolog examples, including complex libraries involving foreign code.

Although it is our aim to enhance compatibility, we are still faced with many incompatibilities between the dialects. As a first step both YAP and SWI will provide some instruments that help developing portable code. A first release of these tools appeared in SWI-Prolog 5.6.43. Some of the facilities are implemented in the base system, others in the library `library(dialect.pl)`.

- The Prolog flag [dialect](flags.html#flag:dialect) is an unambiguous and fast way to find out which Prolog dialect executes your program. It has the value `swi` for SWI-Prolog and `yap` on YAP.
- The Prolog flag [version_data](flags.html#flag:version_data) is bound to a term `swi(Major, Minor, Patch, Extra)`
- Conditional compilation using `:- if(Condition)` ...`:- endif` is supported. See [section 4.3.1.2](consulting.html#sec:4.3.1.2).
- The predicate [expects_dialect/1](dialect.html#expects_dialect/1) allows for specifying for which Prolog system the code was written.
- The predicates [exists_source/1](consulting.html#exists_source/1) and [source_exports/2](dialect.html#source_exports/2) can be used to query the library content. The [require/1](consulting.html#require/1) directive can be used to get access to predicates without knowing their location.
- The module predicates [use_module/1](import.html#use_module/1), [use_module/2](import.html#use_module/2) have been extended with a notion for‘import-except’and‘import-as’. This is particularly useful together with [reexport/1](reexport.html#reexport/1) and [reexport/2](reexport.html#reexport/2) to compose modules from other modules and mapping names.
- Foreign code can expect `__SWI_PROLOG__` when compiled for SWI-Prolog and `__YAP_PROLOG__` when compiled on YAP.

:- **expects_dialect**(`+Dialect`)  
This directive states that the code following the directive is written for the given Prolog `Dialect`. See also [dialect](flags.html#flag:dialect). The declaration holds until the end of the file in which it appears. The current dialect is available using [prolog_load_context/2](consulting.html#prolog_load_context/2).

The exact behaviour of this predicate is still subject to discussion. Of course, if `Dialect` matches the running dialect the directive has no effect. Otherwise we check for the existence of `library(dialect/Dialect)` and load it if the file is found. Currently, this file has this functionality:

- Define system predicates of the requested dialect we do not have.
- Apply [goal_expansion/2](consulting.html#goal_expansion/2) rules that map conflicting predicates to versions emulating the requested dialect. These expansion rules reside in the dialect compatibility module, but are applied if prolog_load_context(dialect, Dialect) is active.
- Modify the search path for library directories, putting libraries compatible with the target dialect before the native libraries.
- Setup support for the default filename extension of the dialect.

**source_exports**(`+Spec, +Export`)  
Is true if source `Spec` exports `Export`, a predicate indicator. Fails without error otherwise.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[C.1 Some considerations for writing portable code](portabilitystrategies.html)

[C.2 Notes on specific dialects](dialect-notes.html)

[C.2.1 Notes on specific dialects](dialect-notes.html#sec:C.2.1)

[C.2.1.1 Loading XSB source files](dialect-notes.html#sec:C.2.1.1)

[C.2.2 The XSB import directive](dialect-notes.html#sec:C.2.2)
