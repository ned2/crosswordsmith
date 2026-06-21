
## 2.16 Rational trees (cyclic terms)

SWI-Prolog supports rational trees, also known as cyclic terms.‘Supports’is so defined that most relevant built-in predicates terminate when faced with rational trees. Almost all SWI-Prolog's built-in term manipulation predicates process terms in a time that is linear to the amount of memory used to represent the term on the stack. The following set of predicates safely handles rational trees: [=../2](manipterm.html#=../2), [==/2](compare.html#==/2), [=@=/2](compare.html#=@=/2), [=/2](compare.html#=/2), [@\</2](compare.html#@%3C/2), [@=\</2](compare.html#@=%3C/2), [@\>=/2](compare.html#@%3E=/2), [@\>/2](compare.html#@%3E/2), [\\=/2](compare.html#\==/2), [\\@=/2](compare.html#\=@=/2), [\\/2](compare.html#\=/2), [acyclic_term/1](typetest.html#acyclic_term/1), [bagof/3](allsolutions.html#bagof/3), [compare/3](compare.html#compare/3), [copy_term/2](manipterm.html#copy_term/2), [cyclic_term/1](typetest.html#cyclic_term/1), [dif/2](coroutining.html#dif/2), [duplicate_term/2](manipterm.html#duplicate_term/2), [findall/3](allsolutions.html#findall/3), [ground/1](typetest.html#ground/1), [term_hash/2](db.html#term_hash/2), [numbervars/3](manipterm.html#numbervars/3), [numbervars/4](manipterm.html#numbervars/4), [recorda/3](db.html#recorda/3), [recordz/3](db.html#recordz/3), [setof/3](allsolutions.html#setof/3), [subsumes_term/2](compare.html#subsumes_term/2), [term_variables/2](manipterm.html#term_variables/2), [throw/1](exception.html#throw/1), [unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2), [unifiable/3](compare.html#unifiable/3), [when/2](coroutining.html#when/2), [write/1](termrw.html#write/1) (and related predicates) .

In addition, some built-ins recognise rational trees and raise an appropriate exception. Arithmetic evaluation belongs to this group. The compiler ([asserta/1](db.html#asserta/1), etc.) also raises an exception. Future versions may support rational trees. Predicates that could provide meaningful processing of rational trees raise a `representation_error`. Predicates for which rational trees have no meaningful interpretation raise a `type_error`. For example:

``` code
1 ?- A = f(A), asserta(a(A)).
ERROR: asserta/1: Cannot represent due to `cyclic_term'
2 ?- A = 1+A, B is A.
ERROR: is/2: Type error: `expression' expected, found
             `@(S_1,[S_1=1+S_1])' (cyclic term)
```
