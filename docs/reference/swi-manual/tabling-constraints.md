
## 7.10 Tabling and constraints

Starting with version 9.3.24, SWI-Prolog offers some support for tabled execution with constraints (attributed variables, see [section 8.1](attvar.html#sec:8.1)) All relevant data structures support attributed variables, notably tries (see [section 4.14.4](db.html#sec:4.14.4)). The basic attributed variable and tabling with attributed variables tests from XSB have been ported and integrated in SWI-Prolog's test suite.^(192Thanks to Theresa Swift and David Warren.) Some remarks:

- The *variant* is defined by the attributes and their values. Note however that SWI-Prolog represents multiple attributes in a linked list where the ordering depends on the order in which the attributes were added while, ideally, the order is not relevant for the semantics of attributes.

- Solving a goal with attributed variables may *modify* attributes. As a result, enumeration of answers from the completed trie *replaces* attributes rather than unifying with the attributes. The new set of attributes is always a copy of the original set of attributes. For example:

  ``` code
  :- use_module(library(clpfd)).
  :- table p/1.

  p(X) :- X #>= 1, X #=< 6.
  p(20).
  ```

  ``` code
  ?- X #> 0, p(X).
  X in 1..6 ;
  X = 20.
  ```

  Note that this behaviour is unlike [trie_gen/2](db.html#trie_gen/2). If an attributed variable is inserted into a trie, [trie_gen/2](db.html#trie_gen/2) unifies the stored attributed term with the second argument of the call.
