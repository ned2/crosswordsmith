
## 9.5 CHR Examples

Here are two example constraint solvers written in CHR.

- The program below defines a solver with one constraint, leq/2 , which is a less-than-or-equal constraint, also known as a partial order constraint.

  ``` code
  :- module(leq,[leq/2]).
  :- use_module(library(chr)).

  :- chr_constraint leq/2.
  reflexivity  @ leq(X,X) <=> true.
  antisymmetry @ leq(X,Y), leq(Y,X) <=> X = Y.
  idempotence  @ leq(X,Y) \ leq(X,Y) <=> true.
  transitivity @ leq(X,Y), leq(Y,Z) ==> leq(X,Z).
  ```

  When the above program is saved in a file and loaded in SWI-Prolog, you can call the leq/2 constraints in a query, e.g.:

  ``` code
  ?- leq(X,Y), leq(Y,Z).
  leq(_G23837, _G23841)
  leq(_G23838, _G23841)
  leq(_G23837, _G23838)
  true .
  ```

  When the query succeeds, the SWI-Prolog top level prints the content of the CHR constraint store and displays the bindings generated during the query. Some of the query variables may have been bound to attributed variables, as you see in the above example.

- The program below implements a simple finite domain constraint solver.

  ``` code
  :- module(dom,[dom/2]).
  :- use_module(library(chr)).

  :- chr_constraint dom(?int,+list(int)).
  :- chr_type list(T) ---> [] ; [T|list(T)].

  dom(X,[]) <=> fail.
  dom(X,[Y]) <=> X = Y.
  dom(X,L) <=> nonvar(X) | memberchk(X,L).
  dom(X,L1), dom(X,L2) <=> intersection(L1,L2,L3), dom(X,L3).
  ```

  When the above program is saved in a file and loaded in SWI-Prolog, you can call the dom/2 constraints in a query, e.g.:

  ``` code
  ?- dom(A,[1,2,3]), dom(A,[3,4,5]).
  A = 3.
  ```
