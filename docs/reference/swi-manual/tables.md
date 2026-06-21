
## A.58 library(tables): XSB interface to tables

This module provides an XSB compatible library to access tables as created by tabling (see [table/1](tabling-preds.html#table/1)). The aim of this library is first of all compatibility with XSB. This library contains some old and internal XSB predicates that are marked deprecated.

**tnot**(`:Goal`)  
Tabled negation.

deprecated  
This is a synonym to [tnot/1](tabling-preds.html#tnot/1).

**tfindall**(`+Template, :Goal, -Answers`)  
This predicate emerged in XSB in an attempt to provide a safer alternative to [findall/3](allsolutions.html#findall/3). This doesn't really work in XSB and the SWI-Prolog emulation is a simple call to [findall/3](allsolutions.html#findall/3). Note that `Goal` may not be a variant of an *incomplete* table.

deprecated  
Use [findall/3](allsolutions.html#findall/3)

**set_pil_on**  
**set_pil_off**  
Dummy predicates for XSB compatibility.

deprecated  
These predicates have no effect.

\[semidet\]**get_call**(`:CallTerm, -Trie, -Return`)  
True when `Trie` is an answer trie for a variant of `CallTerm`. `Return` is a term ret/N with N variables that share with variables in `CallTerm`. The `Trie` contains zero or more instances of the `Return` term. See also [get_calls/3](tables.html#get_calls/3).

\[nondet\]**get_calls**(`:CallTerm, -Trie, -Return`)  
True when `Trie` is an answer trie for a variant that unifies with `CallTerm` and Skeleton is the answer skeleton. See [get_call/3](tables.html#get_call/3) for details.

\[nondet\]**get_returns**(`+ATrie, -Return`)  
True when `Return` is an answer template for the AnswerTrie.

|  |  |
|----|----|
| `Return` | is a term `ret(...)`. See [get_calls/3](tables.html#get_calls/3). |

\[nondet\]**get_returns**(`+AnswerTrie, -Return, -NodeID`)  
True when `Return` is an answer template for the `AnswerTrie` and the answer is represented by the trie node `NodeID`.

|  |  |
|----|----|
| `Return` | is a term `ret(...)`. See [get_calls/3](tables.html#get_calls/3). |

\[nondet\]**get_returns_and_tvs**(`+AnswerTrie, -Return, -TruthValue`)  
Identical to [get_returns/2](tables.html#get_returns/2), but also obtains the truth value of a given answer, setting `TruthValue` to `t` if the answer is unconditional and to `u` if it is conditional. If a conditional answer has multiple delay lists, this predicate will succeed only once, so that using this predicate may be more efficient than [get_residual/2](tables.html#get_residual/2) (although less informative)

\[nondet\]**get_returns_and_dls**(`+AnswerTrie, -Return, :DelayLists`)  
True when `Return` appears in `AnswerTrie` with the given `DelayLists`. `DelayLists` is a list of lists, where the inner lists expresses a conjunctive condition and and outer list a disjunction.

\[nondet\]**get_residual**(`:CallTerm, -DelayList`)  
True if `CallTerm` appears in a table and has `DelayList`. SWI-Prolog's representation for a delay is a body term, more specifically a disjunction of conjunctions. The XSB representation is non-deterministic and uses a list to represent the conjunction.

The delay condition is a disjunction of conjunctions and is represented as such in the native SWI-Prolog interface as a nested term of ;/2 and ,/2, using `true` if the answer is unconditional. This XSB predicate returns the associated conjunctions non-deterministically as a list.

See also [call_residual_program/2](WFS.html#call_residual_program/2) from `library(wfs)`.

\[nondet\]**get_returns_for_call**(`:CallTerm, -AnswerTerm`)  
True if `AnswerTerm` appears in the tables for the *variant* `CallTerm`.

**abolish_table_pred**(`:CallTermOrPI`)  
Invalidates all tabled subgoals for the predicate denoted by the predicate or term indicator Pred.

To be done  
If Pred has a subgoal that contains a conditional answer, the default behavior will be to transitively abolish any tabled predicates with subgoals having answers that depend on any conditional answers of S.

\[det\]**abolish_table_call**(`+Head`)  
\[det\]**abolish_table_call**(`+Head, +Options`)  
Same as [abolish_table_subgoals/1](tabling-preds.html#abolish_table_subgoals/1). See also [abolish_table_pred/1](tables.html#abolish_table_pred/1).

deprecated  
Use abolish_table_subgoals/\[1,2\].

**abolish_table_subgoals**(`:Head, +Options`)  
Behaves as [abolish_table_subgoals/1](tabling-preds.html#abolish_table_subgoals/1), but allows the default `table_gc_action` to be over-ridden with a flag, which can be either `abolish_tables_transitively` or `abolish_tables_singly`.

Compatibility  
`Options` is compatible with XSB, but does not follow the ISO option handling conventions.
