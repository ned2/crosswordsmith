
## 7.13 About the tabling implementation

The SWI-Prolog implementation uses *Delimited continuations* (see [section 4.9](delcont.html#sec:4.9) to realise suspension of variant calls. The initial version was written by Benoit Desouter and described in [Desouter *et al.*, 2015](Bibliography.html#DBLP:journals/tplp/DesouterDS15). We moved the main data structures required for tabling, the *answer tables* (see [section 4.14.4](db.html#sec:4.14.4)) and the *worklist* to SWI-Prolog's C core. *Mode directed tabling* ([section 7.3](tabling-mode-directed.html#sec:7.3)) is based on a prototype implementation by Fabrizio Riguzzi.

The implementation of dynamic SCCs, dynamically stratified negation and Well Founded Semantics was initiated by Benjamin Grosof from Kyndi and was realised with a lot of help by Theresa Swift, David Warren and Fabrizio Riguzzi, as well as publications about XSB [Sagonas & Swift, 1998](Bibliography.html#DBLP:journals/toplas/SagonasS98), [Sagonas *et al.*, 2000](Bibliography.html#SAGONAS20001).

The [table/1](tabling-preds.html#table/1) directive causes the creation of a wrapper calling the renamed original predicate. For example, the program in [section 7.2](tabling-non-termination.html#sec:7.2) is translated into the following program. We give this information to improve your understanding of the current tabling implementation. Future versions are likely to use a more low-level translation that is not based on wrappers.

``` code
connection(A, B) :-
        start_tabling(user:connection(A, B),
                      'connection tabled'(A, B)).

'connection tabled'(X, Y) :-
        connection(X, Z),
        connection(Z, Y).
'connection tabled'(X, Y) :-
        connection(Y, X).

'connection tabled'('Amsterdam', 'Schiphol').
'connection tabled'('Amsterdam', 'Haarlem').
'connection tabled'('Schiphol', 'Leiden').
'connection tabled'('Haarlem', 'Leiden').
```

#### 7.13.1 Status of tabling

The current implementation is merely a first prototype. It needs several enhancements before we can consider it a serious competitor to Prolog systems with mature tabling such as XSB, YAP and B-Prolog. In particular,

- The performance needs to be improved.
- Memory usage needs to be reduced.
- Tables must be shared between threads, both to reduce space and avoid recomputation.
- Tables must be invalidated and reclaimed automatically.
- Notably XSB supports incremental tabling and well-founded semantics under negation.
