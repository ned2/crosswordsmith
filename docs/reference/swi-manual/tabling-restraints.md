
## 7.11 Tabling restraints: bounded rationality and tripwires

Tabling avoids non-termination due to *self-recursion*. As Prolog allows for infinitely nested *compound terms* (*function symbols* in logic) and arbitrary numbers, the set of possible answers is not finite and thus there is no guaranteed termination.

This section describes *restraints* [Grosof & Swift, 2013](Bibliography.html#DBLP:conf/aaai/GrosofS13) that can be enforced to specific or all tabled predicates. Currently there are three defined restraints, limiting (1) the size of (the arguments to) goals, (2) the size of the answer substitution added to a table and (3) the number of answers allowed in any table. If any of these events occurs we can specify the action taken. We distinguish two classes of actions. First, these events can trap a *tripwire* which can be handled using a hook or a predefined action such as raising an exception, printing a warning or enter a *break level*. This can be used for limiting resources, be notified of suspicious events (debugging) or dynamically adjust the (tabling) strategy of the program. Second, they may continue the computation that results in a partial answer (*bounded rationality*). Unlike just not exploring part of the space though, we use the third truth value of well founded semantics to keep track of answers that have not been affected by the restraints and those that have been affected.

The tripwire actions apply for all restraints. If a tripwire action is triggered, the system takes the steps below.

1.  Call the [prolog:tripwire/2](tabling-restraints.html#prolog:tripwire/2) hook.
2.  If [prolog:tripwire/2](tabling-restraints.html#prolog:tripwire/2) fails, take one of the predefined actions:
    **warning**  
    Print a message indicating the trapped tripwire and continue execution as normal, i.e., the final answer is the same as if no restraint was active.

    **error**  
    Throw an exception `error(resource_error(tripwire(Wire,Context)))`.

    **suspend**  
    Print a message and start a *break level* (see [break/0](toplevel.html#break/0)).

\[multifile\]**prolog:tripwire**(`Wire, Context`)  
Called when tripwire `Wire` is trapped. `Context` provides additional context for interpreting the tripwire. The hook can take one of three actions:

- Succeed. In this case the tripwire is considered handled and execution proceeds as if there was no tripwire. Note that tripwires only trigger at the exact value, which implies that a wire on a count will be triggered only once. The hook can install a new tripwire at a higher count.
- Fail. In this case the default action is taken.
- Throw an exception. Exceptions are propagated normally.

*Radial restraints* limit the sizes of subgoals or answers. Abstraction of a term according to the size limit is implemented by [size_abstract_term/3](tabling-restraints.html#size_abstract_term/3).

\[det\]**size_abstract_term**(`+Size, +Term, -Abstract`)  
The size of a term is defined as the number of compound subterms (*function symbols*) that appear in term. `Abstract` is an abstract copy of `Term` where each argument is abstracted by copying only the first `Size` function symbols and constants. Excess function symbols are replaced by fresh variables.

This predicate is a helper for tabling where `Term` is the `ret/N` *answer skeleton* that is added to the *answer table*. Examples:

|          |                    |                  |
|:--------:|--------------------|------------------|
| **Size** | **Term**           | **Abstract**     |
|    0     | ret(f(x), a)       | ret(\_, a)       |
|    1     | ret(f(x), a)       | ret(f(x), a)     |
|    1     | ret(f(A), a)       | ret(f(A), a)     |
|    1     | ret(f(x), x(y(Z))) | ret(f(x), x(\_)) |

\[undefined\]**radial_restraint**  
This predicate is *undefined* in the sense of well founded semantics (see [section 7.6](WFS.html#sec:7.6) and [undefined/0](WFS.html#undefined/0)). Any answer that depends on this condition is undefined because either the restraint on the subgoal size or answer size was violated.

### 7.11.1 Restraint subgoal size

Using the `subgoal_abstract(Size)` attribute, a tabled subgoal that that is too large is *abstracted* by replacing compound subterms of the goal with variables. In a nutshell, a goal `p(s(s(s(s(s(0))))))` is converted into the semantically equivalent subgoal if the subgoal size is limited to 3.

``` code
    ...,
    p(s(s(s(X)))), X = s(s(0)),
    ...,
```

As a result of this, terms stored in the *variant trie* that maps goal variants into *answer tables* is limited. Note that does not limit the number of answer tables as atomic values are never abstracted and there are, for example, an infinite number of integers. Note that restraining the subgoal size does not affect the semantics, provided more general queries on the predicate include all answers that more specific queries do. See also *call substitution* as described in [section 7.5](tabling-subsumptive.html#sec:7.5). In addition to the tripwire actions, the [max_table_subgoal_size_action](flags.html#flag:max_table_subgoal_size_action) can be set to the value `abstract`:

**abstract**  
Abstract the goal as described above and provide correctness by adding the required unification instructions after the goal.

### 7.11.2 Restraint answer size

Using the `answer_abstract(Size)` attribute, a tabled subgoal that produces answer substitutions (instances of the variables in the goal) whose size exceed `Size` are trapped. In addition to the tripwire actions, answer abstraction defines two additional modes for dealing with too large answers as defines by the Prolog flag [max_table_answer_size_action](flags.html#flag:max_table_answer_size_action):

**fail**  
Ignore the too large answer. Note that this is semantically incorrect.

**bounded_rationality**  
In this mode, the large answer is *abstracted* in the same way as subgoals are abstracted (see [section 7.11.1](tabling-restraints.html#sec:7.11.1)). This is semantically incorrect, but our third truth value *undefined* is used to remedy this problem. In other words, the abstracted value is added to the table as *undefined* and all conclusions that depend on usage of this abstracted value are thus undefined unless they can also be proved some other way.

### 7.11.3 Restraint answer count

Finally, using “as `max_answers(Count)`” or the Prolog flag [max_answers_for_subgoal](flags.html#flag:max_answers_for_subgoal), the number of answers in a table is restrained. In addition to the tripwire actions this restraint supports the action `bounded_rationality`^(193The action `complete_soundly` is supported as a synonym for XSB compatibility). If the restraint is reached in the bounded rationality mode the system takes the following actions:

- Ignore the answer that triggered the restraint.
- Prune the choice points of the tabled goal to avoid more answers.
- Add an new answer to the table that does not bind any variables, i.e., an empty answer substitution. This answer is conditional on [answer_count_restraint/0](tabling-restraints.html#answer_count_restraint/0).

\[undefined\]**answer_count_restraint**  
This predicate is *undefined* in the sense of well founded semantics (see [section 7.6](WFS.html#sec:7.6) and [undefined/0](WFS.html#undefined/0)). Any answer that depends on this condition is undefined because the `max_answers` restraint on some table was violated.

The program and subsequent query below illustrate the behavior.

``` code
:- table p/2 as max_answers(3).

p(M,N) :-
    between(1,M,N).
```

``` code
?- p(1 000 000, X).
X = 3 ;
X = 2 ;
X = 1 ;
% WFS residual program
    p(1000000, X) :-
        answer_count_restraint.
p(1000000, X).
```
