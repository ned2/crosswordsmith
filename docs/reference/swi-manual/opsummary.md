
## F.4 Operators

|      |     |                    |                           |
|------|-----|--------------------|---------------------------|
| 1    | fx  | \$                 | Bind top-level variable   |
| 200  | xfy | ^                  | Existential qualification |
| 200  | xfy | ^                  | Arithmetic function       |
| 300  | xfx | mod                | Arithmetic function       |
| 400  | yfx | \*                 | Arithmetic function       |
| 400  | yfx | /                  | Arithmetic function       |
| 400  | yfx | //                 | Arithmetic function       |
| 400  | yfx | \<\<               | Arithmetic function       |
| 400  | yfx | \>\>               | Arithmetic function       |
| 400  | yfx | xor                | Arithmetic function       |
| 500  | fx  | \+                 | Arithmetic function       |
| 500  | fx  | \-                 | Arithmetic function       |
| 500  | fx  | ?                  | XPCE: obtainer            |
| 500  | fx  | \\                 | Arithmetic function       |
| 500  | yfx | \+                 | Arithmetic function       |
| 500  | yfx | \-                 | Arithmetic function       |
| 500  | yfx | /\\                | Arithmetic function       |
| 500  | yfx | \\                 | Arithmetic function       |
| 600  | xfy | :                  | module:term separator     |
| 700  | xfx | \<                 | Predicate                 |
| 700  | xfx | =                  | Predicate                 |
| 700  | xfx | =..                | Predicate                 |
| 700  | xfx | =:=                | Predicate                 |
| 700  | xfx | \<                 | Predicate                 |
| 700  | xfx | ==                 | Predicate                 |
| 700  | xfx | =@=                | Predicate                 |
| 700  | xfx | =\\                | Predicate                 |
| 700  | xfx | \>                 | Predicate                 |
| 700  | xfx | \>=                | Predicate                 |
| 700  | xfx | @\<                | Predicate                 |
| 700  | xfx | @=\<               | Predicate                 |
| 700  | xfx | @\>                | Predicate                 |
| 700  | xfx | @\>=               | Predicate                 |
| 700  | xfx | as                 | Predicate                 |
| 700  | xfx | is                 | Predicate                 |
| 700  | xfx | \\                 | Predicate                 |
| 700  | xfx | \\=                | Predicate                 |
| 700  | xfx | =@=                | Predicate                 |
| 900  | fy  | not                | Predicate                 |
| 900  | fy  | \\                 | Predicate                 |
| 1000 | xfy | ,                  | Predicate                 |
| 1050 | xfy | -\>                | Predicate                 |
| 1050 | xfy | \*-\>              | Predicate                 |
| 1100 | xfy | ;                  | Predicate                 |
| 1105 | xfy | \|                 | DCG disjunction           |
| 1150 | fx  | discontiguous      | Directive                 |
| 1150 | fx  | dynamic            | Directive                 |
| 1150 | fx  | module_transparent | Directive                 |
| 1150 | fx  | meta_predicate     | Head                      |
| 1150 | fx  | multifile          | Directive                 |
| 1150 | fx  | thread_local       | Directive                 |
| 1150 | fx  | volatile           | Directive                 |
| 1150 | fx  | initialization     | Directive                 |
| 1200 | fx  | :-                 | Introduces a directive    |
| 1200 | fx  | ?-                 | Introduces a directive    |
| 1200 | xfx | --\>               | DCGrammar: rewrite        |
| 1200 | xfx | ==\>               | DCGrammar: rewrite        |
| 1200 | xfx | :-                 | head `:-` body. separator |
