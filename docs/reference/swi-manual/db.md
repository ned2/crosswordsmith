
## 4.14 Database

SWI-Prolog offers several ways to store data in globally accessible memory, i.e., outside the Prolog *stacks*. Data stored this way notably does not change on *backtracking*. Typically it is a bad idea to use any of the predicates in this section for realising global variables that can be assigned to. Typically, first consider representing data processed by your program as terms passed around as predicate arguments. If you need to reason over multiple solutions to a goal, consider [findall/3](allsolutions.html#findall/3), [aggregate/3](aggregate.html#aggregate/3) and related predicates.

Nevertheless, there are scenarios where storing data outside the Prolog stacks is a good option. Below are the main options for storing data:

**Using dynamic predicates**  
Dynamic predicates are predicates for which the list of clauses is modified at runtime using [asserta/1](db.html#asserta/1), [assertz/1](db.html#assertz/1), [retract/1](db.html#retract/1) or [retractall/1](db.html#retractall/1). Following the ISO standard, predicates that are modified this way need to be declared using the [dynamic/1](dynamic.html#dynamic/1) *directive*. These facilities are defined by the ISO standard and widely supported. The mechanism is often considered slow in the literature. Performance depends on the Prolog implementation. In SWI-Prolog, querying dynamic predicates has the same performance as static ones. The manipulation predicates are fast. Using [retract/1](db.html#retract/1) or [retractall/1](db.html#retractall/1) on a predicate registers the predicate as‘dirty’. Dirty predicates are cleaned by [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0), which is normally automatically invoked. Some workloads may result in significant performance reduction due to skipping retracted clauses and/or clause garbage collection.

Dynamic predicates can be wrapped using library `library(persistency)` to maintain a backup of the data on disk. Dynamic predicates come in two flavours, *shared* between threads and *local* to each thread. The latter version is created using the directive [thread_local/1](threadcom.html#thread_local/1).

**The recorded database**  
The‘recorded database’registers a list of terms with a *key*, an atom or compound term. The list is managed using [recorda/3](db.html#recorda/3), [recordz/3](db.html#recordz/3) and [erase/1](db.html#erase/1). It is queried using [recorded/3](db.html#recorded/3). The recorded database is not part of the ISO standard but fairly widely supported, notably in implementations building on the‘Edinburgh tradition’. There are few reasons to use this database in SWI-Prolog due to the good performance of dynamic predicates. Advantages are (1) the handle provides a direct reference to a term, (2) cyclic terms can be stored and (3) attributes ([section 8.1](attvar.html#sec:8.1)) are preserved. Disadvantages are (1) the terms in a list associated with a key are not indexed, (2) the poorly specified *immediate update semantics* (see [section 4.14.1.1](db.html#sec:4.14.1.1) applies to the recorded database and (3) reduced portability.

**The [flag/3](db.html#flag/3) predicate**  
The predicate [flag/3](db.html#flag/3) associates one simple value (number or atom) with a key (atom, integer or compound). It is an old SWI-Prolog specific predicate that should be considered deprecated, although there is no plan to remove it.

**Using global variables**  
The predicates [b_setval/2](gvar.html#b_setval/2) and [nb_setval/2](gvar.html#nb_setval/2) associate a term living on the Prolog stack with a name, either backtrackable or non-backtrackable. Backtrackable and non-backtrackable assignment without using a global name can be realised with [setarg/3](manipterm.html#setarg/3) and [nb_setarg/3](manipterm.html#nb_setarg/3). Notably the latter are used to realise aggregation as e.g., [aggregate_all/3](aggregate.html#aggregate_all/3) performs.

**Tries**  
As of version 7.3.21, SWI-Prolog provides *tries* (prefix trees) to associate a term *variant* with a value. Tries have been introduced to support `tabling` and are described in [section 4.14.4](db.html#sec:4.14.4).

### 4.14.1 Managing (dynamic) predicates

\[ISO\]**abolish**(`:PredicateIndicator`)  
Removes all clauses of a predicate with functor `Functor` and arity `Arity` from the database. All predicate attributes (dynamic, multifile, index, etc.) are reset to their defaults. Abolishing an imported predicate only removes the import link; the predicate will keep its old definition in its definition module.

According to the ISO standard, [abolish/1](db.html#abolish/1) can only be applied to dynamic procedures. This is odd, as for dealing with dynamic procedures there is already [retract/1](db.html#retract/1) and [retractall/1](db.html#retractall/1). The [abolish/1](db.html#abolish/1) predicate was introduced in DEC-10 Prolog precisely for dealing with static procedures. In SWI-Prolog, [abolish/1](db.html#abolish/1) works on static procedures, unless the Prolog flag [iso](flags.html#flag:iso) is set to `true`.

It is advised to use [retractall/1](db.html#retractall/1) for erasing all clauses of a dynamic predicate.

**abolish**(`+Name, +Arity`)  
Same as `abolish(Name/Arity)`. The predicate [abolish/2](db.html#abolish/2) conforms to the Edinburgh standard, while [abolish/1](db.html#abolish/1) is ISO compliant.

**copy_predicate_clauses**(`:From, :To`)  
Copy all clauses of predicate `From` to `To`. The predicate `To` must be dynamic or undefined. If `To` is undefined, it is created as a dynamic predicate holding a copy of the clauses of `From`. If `To` is a dynamic predicate, the clauses of `From` are added (as in [assertz/1](db.html#assertz/1)) to the clauses of `To`. `To` and `From` must have the same arity. Acts as if defined by the program below, but at a much better performance by avoiding decompilation and compilation.

``` code
copy_predicate_clauses(From, To) :-
        head(From, MF:FromHead),
        head(To, MT:ToHead),
        FromHead =.. [_|Args],
        ToHead =.. [_|Args],
        forall(clause(MF:FromHead, Body),
               assertz(MT:ToHead, Body)).

head(From, M:Head) :-
        strip_module(From, M, Name/Arity),
        functor(Head, Name, Arity).
```

**redefine_system_predicate**(`+Head`)  
This directive may be used both in module `user` and in normal modules to redefine any system predicate. If the system definition is redefined in module `user`, the new definition is the default definition for all sub-modules. Otherwise the redefinition is local to the module. The system definition remains in the module `system`.

Redefining system predicate facilitates the definition of compatibility packages. Use in other contexts is discouraged.

\[ISO,nondet\]**retract**(`+Term`)  
When `Term` is an atom or a term it is unified with the first unifying fact or clause in the database. The fact or clause is removed from the database. The [retract/1](db.html#retract/1) predicate respects the *logical update view*. This implies that [retract/1](db.html#retract/1) succeeds for all clauses that match `Term` when the predicate was *called*. The example below illustrates that the first call to [retract/1](db.html#retract/1) succeeds on `bee` on backtracking despite the fact that `bee` is already retracted.^(87Example by Jan Burse)

``` code
:- dynamic insect/1.
insect(ant).
insect(bee).

?- (   retract(insect(I)),
       writeln(I),
       retract(insect(bee)),
       fail
   ;   true
   ).
ant ;
bee.
```

If multiple threads start a retract on the same predicate at the same time their notion of the *entry generation* is adjusted such that they do not retract the same first clause. This implies that, if multiple threads use `once(retract(Term))`, no two threads will retract the same clause. Note that on backtracking over [retract/1](db.html#retract/1), multiple threads may retract the same clause as both threads respect the logical update view.

\[ISO,det\]**retractall**(`+Head`)  
All facts or clauses in the database for which the `head` unifies with `Head` are removed. If all arguments of `Head` are non-sharing variables (see [is_most_general_term/1](manipterm.html#is_most_general_term/1)), all clauses are removed without inspecting the clauses. Cleaning all clauses of a dynamic predicate must use [retractall/1](db.html#retractall/1) rather than [abolish/1](db.html#abolish/1) as the latter completely wipes the predicate, including its properties. If `Head` refers to a predicate that is not defined, it is implicitly created as a dynamic predicate. See also [dynamic/1](dynamic.html#dynamic/1).^(88The ISO standard only allows using [dynamic/1](dynamic.html#dynamic/1) as a *directive*.)

\[ISO\]**asserta**(`+Term`)  
\[ISO\]**assertz**(`+Term`)  
\[deprecated\]**assert**(`+Term`)  
Assert a clause (fact or rule) into the database. The predicate [asserta/1](db.html#asserta/1) asserts the clause as first clause of the predicate while [assertz/1](db.html#assertz/1) assert the clause as last clause. The deprecated [assert/1](db.html#assert/1) is equivalent to [assertz/1](db.html#assertz/1). If the program space for the target module is limited (see [set_module/1](manipmodule.html#set_module/1)), [asserta/1](db.html#asserta/1) can raise a `resource_error(program_space)` exception. The example below adds two facts and a rule. Note the double parentheses around the rule.

``` code
?- assertz(parent('Bob', 'Jane')).
?- assertz(female('Jane')).
?- assertz((mother(Child, Mother) :-
                parent(Child, Mother),
                female(Mother))).
```

**asserta**(`+Term, -Reference`)  
**assertz**(`+Term, -Reference`)  
\[deprecated\]**assert**(`+Term, -Reference`)  
Equivalent to [asserta/1](db.html#asserta/1), [assertz/1](db.html#assertz/1), [assert/1](db.html#assert/1), but in addition unifies `Reference` with a handle to the asserted clauses. The handle can be used to access this clause with [clause/3](examineprog.html#clause/3) and [erase/1](db.html#erase/1).

#### 4.14.1.1 Update view

SWI-Prolog adheres to the *logical update view*, where backtrackable predicates that enter the definition of a predicate will not see any changes (either caused by [assert/1](db.html#assert/1) or [retract/1](db.html#retract/1)) to the predicate. This view is the ISO standard. Logical updates are realised by keeping *generation* information on clauses. Each change to the database causes an increment of the generation of the database. Each goal is tagged with the generation in which it was started. Each clause is flagged with the generation it was created in as well as the generation it was erased. Only clauses with a‘created’ ...‘erased’interval that encloses the generation of the current goal are considered visible. The generation mechanism is also used to implement *transactions* See [section 4.14.1.3](db.html#sec:4.14.1.3).

Erased clauses are (eventually) reclaimed by the clause garbage collector implemented by [garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0). By default, the clause garbage collector runs in a thread named `gc`, together with the atom garbage collector ([garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0)). See also the Prolog flag [gc_thread](flags.html#flag:gc_thread).

#### 4.14.1.2 Indexing databases

The indexing capabilities of SWI-Prolog are described in [section 2.17](jitindex.html#sec:2.17). Summarizing, SWI-Prolog creates indexes for any applicable argument, pairs of arguments and indexes on the arguments of compound terms when applicable. Extended JIT indexing is not widely supported among Prolog implementations. Programs that aim at portability should consider using [term_hash/2](db.html#term_hash/2) and [term_hash/4](db.html#term_hash/4) to design their database such that indexing on constant or functor (name/arity reference) on the first argument is sufficient. In some cases, using the predicates below to add one or more additional columns (arguments) to a database predicate may improve performance. The overall design of code using these predicates is given below. Note that as [term_hash/2](db.html#term_hash/2) leaves the hash unbound if `Term` is not ground. This causes the lookup to be fast if `Term` is ground and correct (but slow) otherwise.

``` code
:- dynamic
    x/2.

assert_x(Term) :-
    term_hash(Term, Hash),
    assertz(x(Hash, Term)).

x(Term) :-
    term_hash(Term, Hash),
    x(Hash, Term).
```

\[det\]**term_hash**(`+Term, -HashKey`)  
If `Term` is a ground term (see [ground/1](typetest.html#ground/1)), `HashKey` is unified with a positive integer value that may be used as a hash key to the value. If `Term` is not ground, the predicate leaves `HashKey` an unbound variable. Hash keys are in the range `0 ... 72,057,594,037,927,935` (`2^56`), the maximal integer that can be stored efficiently (see [max_tagged_integer](flags.html#flag:max_tagged_integer).

This predicate may be used to build hash tables as well as to exploit argument indexing to find complex terms more quickly.

The hash key does not rely on temporary information like addresses of atoms and may be assumed constant over different invocations and versions of SWI-Prolog.^(89Last change: version 9.3.28) The [term_hash/2](db.html#term_hash/2) predicate is cycle-safe.^(bugAll arguments that (indirectly) lead to a cycle have the same hash key.) ^(bugHashes differ between big and little endian machines and hashes for big integers (currently `> 2^56`) differ depending on the big integer *backend* (GMP or LibBF). See [gmp_version](flags.html#flag:gmp_version).)

\[det\]**term_hash**(`+Term, +Depth, +Range, -HashKey`)  
As [term_hash/2](db.html#term_hash/2), but only considers `Term` to the specified `Depth`. The top-level term has depth 1, its arguments have depth 2, etc. That is, `Depth`` = 0` hashes nothing; `Depth`` = 1` hashes atomic values or the functor and arity of a compound term, not its arguments; `Depth`` = 2` also indexes the immediate arguments, etc. If `Depth`` = -1`, no depth limit is applied.

`HashKey` is in the range `[0 ...``Range``-1]`. `Range` must be in the range `[1 ... 2147483647]`.

\[det\]**variant_sha1**(`+Term, -SHA1`)  
Compute a SHA1-hash from `Term`. The hash is represented as a 40-byte hexadecimal atom. Unlike [term_hash/2](db.html#term_hash/2) and friends, this predicate produces a hash key for non-ground terms. The hash is invariant over variable-renaming (see [=@=/2](compare.html#=@=/2)) and constants over different invocations of Prolog.^(bugThe hash depends on word order (big/little-endian) and the wordsize (32/64 bits).)

This predicate raises an exception when trying to compute the hash on a cyclic term or attributed term. Attributed terms are not handled because [subsumes_chk/2](terms.html#subsumes_chk/2) is not considered well defined for attributed terms. Cyclic terms are not supported because this would require establishing a canonical cycle. That is, given A=\[a\|A\] and B=\[a,a\|B\], `A` and `B` should produce the same hash. This is not (yet) implemented.

This hash was developed for lookup of solutions to a goal stored in a table. By using a cryptographic hash, heuristic algorithms can often ignore the possibility of hash collisions and thus avoid storing the goal term itself as well as testing using [=@=/2](compare.html#=@=/2).

\[det\]**variant_hash**(`+Term, -HashKey`)  
Similar to [variant_sha1/2](db.html#variant_sha1/2), but using a non-cryptographic hash and produces an integer result like [term_hash/2](db.html#term_hash/2). This version does deal with attributed variables, processing them as normal variables. This hash is primarily intended to speedup finding variant terms in a set of terms. ^(bugAs [variant_sha1/2](db.html#variant_sha1/2), cyclic terms result in an exception.)

#### 4.14.1.3 Transactions

Traditionally, Prolog database updates add or remove individual clauses. The *Logical Update View* ensures that a goal that is started on a dynamic predicate does not see modifications due to [assert/1](db.html#assert/1) or [retract/1](db.html#retract/1) during its life time. See [section 4.14.1.1](db.html#sec:4.14.1.1). In a multi-threaded context this assumption still holds for individual predicates: concurrent modifications to a dynamic predicate are invisible.

*Transactions* allow running a goal in *isolation*. The goals running inside the transaction‘see’the database as it was when the transaction was started together with database changes done by the transaction goal. Other threads see no changes until the transaction is *committed*. The commit, also if it involved multiple clauses spread over multiple predicates, becomes *atomically* visible to other threads. Transactions have several benefits [Wielemaker, 2013](Bibliography.html#DBLP:journals/corr/abs-1301-7669)

- If a database update requires multiple [assert/1](db.html#assert/1) and/or [retract/1](db.html#retract/1) operations, a transaction ensure either all are executed or the database remains unchanged. Notably unexpected exceptions or failures cannot leave the database in an inconsistent state.

- Other threads do not see the intermediate inconsistent states when a database update that consists of multiple assert and/or retract is performed in a transaction. This notably avoids the need to use locks (see [with_mutex/2](threadsync.html#with_mutex/2)) in threads that read the data. A reading thread may still need to use [snapshot/1](db.html#snapshot/1) if a goal depends on multiple calls to dynamic predicates. Unlike locks, transaction and snapshot based synchronization allows both readers and writers to make progress simultaneously.^(90*Read-write* locks also provide readers and writers to make progress simultaneously, but readers see all intermediate states rather than a consistent state.)

  Transactions on their own **do not guarantee consistency**. For example, when running the code below to update the temperature concurrently from multiple threads it is possible for the global state to have multiple temperature/1 clauses.

  ``` code
  update_temperature(Temp) :-
      transaction(( retractall(temperature(_)),
                    asserta(temperature(Temp)))).
  ```

  Global *consistency* can be achieved by wrapping the above transaction using [with_mutex/2](threadsync.html#with_mutex/2) or by using [transaction/3](db.html#transaction/3) with a *constraint* that demands a single clause for temperature/1

- Transactions allow for “what if” reasoning over the dynamic database. This is particularly useful when combined with the deductive database facilities provided by tabling (see [section 7](tabling.html#sec:7)).

SWI-Prolog transactions only affect the *dynamic* database. Static predicates are globally visible and shared at all times. In particular, transactions do not affect loading source files and thus, source files loaded inside a transaction (e.g., due to *autoloading*) are immediately globally visible. This may pose problems if loading source files provide clauses for dynamic predicates.

**transaction**(`:Goal`)  
**transaction**(`:Goal, +Options`)  
Run `Goal` as [once/1](metacall.html#once/1) in a transaction. This implies that access to dynamic predicates‘sees’the dynamic predicates at the moment the transaction is started, together with the modifications issued by `Goal`. Thus, `Goal` does not see changes to dynamic predicates from other threads and other threads do not see modifications by `Goal` (*isolation*). If `Goal` succeeds, all modifications become *atomically* visible to the other threads. If `Goal` fails or raises an exception all local modifications are discarded and [transaction/1](db.html#transaction/1) fails or passes the exception.

Currently the number of database changes inside a transaction (or snapshot, see [snapshot/1](db.html#snapshot/1)) is limited to `2 ** 32 -1`. If this limit is exceeded a `representation_error(transaction_generations)` exception is raised.

Transactions may be nested. The above mentioned limitation for the number of database changes applies to the combined number in nested transactions.

If `Goal` succeeds, the transaction is *committed*. This implies that (1) any clause that is asserted in the transaction and not retracted in the same transaction is made *globally visible* and (2) and clause the existed before the transaction and is retracted in the transaction becomes *globally invisible*. Multiple transactions may retract the same clause and be committed, i.e., committing a retract that was already performed is a no-op. All modifications become *atomically* visible to other threads. The [transaction/3](db.html#transaction/3) variation allows for verifying *constraints* just before the commit takes place.

**Clause ordering** Inside a transaction clauses can be added using [asserta/1](db.html#asserta/1) and [assertz/1](db.html#assertz/1). If only a single transaction is active at any point in time transactions preserve the usual ordering of clauses. However, if multiple transactions manipulate the same predicate(s) concurrently (typically using [transaction/3](db.html#transaction/3)), the final order of the clauses is the order in which the transactions asserted the clauses and **not** the order in which the transactions are committed.

The [transaction/1](db.html#transaction/1) variant is equivalent to `transaction(Goal,[])`. The [transaction/2](db.html#transaction/2) variant processed the following options:

**bulk**(`+Boolean`)  
When `true`, accumulate events from changes to dynamic predicates (see [prolog_listen/2](prolog-event.html#prolog_listen/2)) and trigger these events as part of the commit phase. This implies that if the transaction is not committed the events are never triggered. Failure to trigger the events causes the transaction to be discarded. Experimental.

**transaction**(`:Goal, :Constraint, +Mutex`)  
Similar to [transaction/1](db.html#transaction/1), but allows verifying `Constraint` during the commit phase. This predicate follows the steps below. Any failure or exception during this process discards the transaction and releases `Mutex` when applicable. `Constraint` may modify the database. Such modifications follow the semantics that apply for `Goal`.

- Call `once(Goal)`
- Lock `Mutex`
- Change the visibility to the *current* global state combined with the changes made by `Goal`
- Call `once(Constraint)`
- Commit the changes
- Unlock `Mutex`.

This predicate is intended to execute multiple transactions with a time consuming `Goal` in part concurrently. For example, it can be used for a *Compare And Swap* (CAS) like design. We illustrate this using a simple counter in the code below. Note that the transaction fails if some other thread concurrently updated the counter. This is why we need the [repeat/0](control.html#repeat/0) and a final [!/0](control.html#!/0). The CAS-style update is in general useful if `Goal` is expensive and conflicts are rare.

``` code
:- dynamic counter/1.

increment_counter(Delta) :-
    repeat,
      transaction(( counter(Value),
                    Value2 is Value+Delta,
                  ),
                  ( retract(counter(Value)),
                    asserta(counter(Value2))
                  ),
                  counter_lock),
    !.
```

**snapshot**(`:Goal`)  
Similar to [transaction/1](db.html#transaction/1), but *always* discards the local modifications. In other words, [snapshot/1](db.html#snapshot/1) allows a thread to examine a frozen state of the dynamic predicates and/or make isolated modifications without affecting other threads and without making permanent changes to the database. Where transactions allow the global state to be updated atomically from one consistent state to the next, a snapshot allows reasoning about a consistent state.

\[nondet\]**current_transaction**(`-Goal`)  
True when called inside a transaction running `Goal`. This predicate generates candidates from the current (nested) transaction outward. `Goal` is a plain goal if the calling context module is the same as matching [transaction/1](db.html#transaction/1) or [snapshot/1](db.html#snapshot/1) and a qualified callable term otherwise. Note that this only enumerates transactions in the current thread.

**transaction_updates**(`-Updates`)  
Unify `Updates` with a list of database updates that would be effectuated if the transaction is going to be committed at this stage. `Updates` is a list of terms defined below. The elements are sorted on the change generation, i.e., the order in which the operations were performed.

**asserta**(`+ClauseRef`)  
**assertz**(`+ClauseRef`)  
The given clause will be asserted at the start or end. Note that due to competing transactions the clause may no longer be the first/last clause of the predicate.

**erased**(`+ClauseRef`)  
The given clause will be removed. This may be due to [erase/1](db.html#erase/1), [retract/1](db.html#retract/1) or [retractall/1](db.html#retractall/1).

#### 4.14.1.4 Impact of transactions

Transactions interact with other facilities that depend on changing dynamic predicates. This section discusses these interactions.

**Last modified generation**  
Using the [predicate_property/2](examineprog.html#predicate_property/2) property `last_modified_generation(Generation)` we can determine whether a predicate was modified. When a predicate is changed inside a transaction this generation is not updated. The generation for dynamic predicates that are modified in the transaction is updated to the *commit generation* when the transaction is committed. Asking for the last modified generation *inside* the transaction examines the log of modified clauses and reports the generation as one of

- The global modified generation if the predicate was not modified in the transaction and not modified outside the transaction to beyond the start generation of the transaction. If the modified generation is higher than the transaction start generation, this generation is reported. ^(bugNote that the above implies that inside a transaction we observe a changing last modified generation for predicates that have only been modified outside the transaction while these changes are not visible.)
- The transaction start generation plus the local generation of the last change if the predicate is modified inside the transaction.

**Wait for database changes**  
The predicate [thread_wait/2](threadcom.html#thread_wait/2) does not wakeup threads for changes inside a transaction. The wakeup is delayed until the transaction is committed. Note that [thread_wait/2](threadcom.html#thread_wait/2) cannot be meaningfully called from inside a transaction because no external entities can cause changes to the dynamic database inside the transaction.

**Incremental tabling**  
Consistency of tables must be restored if the transaction is rolled back. For local tables this is realised as follows:

- Tables are either marked to be *invalidated* on rollback or, for *monotonic* tabling individual answers are marked to be removed on rollback.
- A table is marked to be *invalidated* if, while it is created or reevaluated, at least one dependent dynamic predicate has been modified inside the transaction.
- Answers are marked to be retracted when they result from monotonic reevaluation based on changes *inside* the transaction.

In other words: tables being reevaluated inside a transaction that do not depend on predicates modified inside the transaction remain valid. Monotonic tables that get new answers due to asserts inside the transaction have these answers removed during the rollback while the table remains valid. Monotonic tables that are for some reason invalidated inside the transaction are invalidated during the rollback.

Correct interaction between tabling and transaction currently **only deals with local tables**. *Shared* tables should not be combined with transactions. Future versions may improve on that. A possible route is to make a local copy from a shared table when (re)evaluation is performed inside a transaction.

**Status** SWI-Prolog transaction basics and API are stable. Interaction with other parts of the system that depend on dynamic predicates is still unsettled. Future versions may support non-determinism through transactions and snapshots.

### 4.14.2 The recorded database

**recorda**(`+Key, +Term, -Reference`)  
Assert `Term` in the recorded database under key `Key`. `Key` is a small integer (range [min_tagged_integer](flags.html#flag:min_tagged_integer) ...[max_tagged_integer](flags.html#flag:max_tagged_integer), atom or compound term. If the key is a compound term, only the name and arity define the key. `Reference` is unified with an opaque handle to the record (see [erase/1](db.html#erase/1)).

**recorda**(`+Key, +Term`)  
Equivalent to `recorda(``Key``, ``Term``, _)`.

**recordz**(`+Key, +Term, -Reference`)  
Equivalent to [recorda/3](db.html#recorda/3), but puts the `Term` at the tail of the terms recorded under `Key`.

**recordz**(`+Key, +Term`)  
Equivalent to `recordz(``Key``, ``Term``, _)`.

**recorded**(`?Key, ?Value, ?Reference`)  
True if `Value` is recorded under `Key` and has the given database `Reference`. If `Reference` is given, this predicate is semi-deterministic. Otherwise, it must be considered non-deterministic. If neither `Reference` nor `Key` is given, the triples are generated as in the code snippet below.^(91Note that, without a given `Key`, some implementations return triples in the order defined by [recorda/2](db.html#recorda/2) and [recordz/2](db.html#recordz/2).) See also [current_key/1](examineprog.html#current_key/1).

``` code
        current_key(Key),
        recorded(Key, Value, Reference)
```

**recorded**(`+Key, -Value`)  
Equivalent to `recorded(``Key``, ``Value``, _)`.

**erase**(`+Reference`)  
Erase a record or clause from the database. `Reference` is a db-reference returned by [recorda/3](db.html#recorda/3), [recordz/3](db.html#recordz/3) or [recorded/3](db.html#recorded/3), [clause/3](examineprog.html#clause/3), [assert/2](db.html#assert/2), [asserta/2](db.html#asserta/2) or [assertz/2](db.html#assertz/2). Fail silently if the referenced object no longer exists. Notably, if multiple threads attempt to erase the same clause one will succeed and the others will fail.

**instance**(`+Reference, -Term`)  
Unify `Term` with the referenced clause or database record. Unit clauses are represented as `Head` :- `true`.

### 4.14.3 Flags

The predicate [flag/3](db.html#flag/3) is the oldest way to store global non-backtrackable data in SWI-Prolog. Flags are global and shared by all threads. Their value is limited to atoms, small (64-bit) integers and floating point numbers. Flags are thread-safe. The flags described in this section must not be confused with *Prolog flags* described in [section 2.12](flags.html#sec:2.12).

**get_flag**(`+Key, -Value`)  
True when `Value` is the value currently associated with `Key`. If `Key` does not exist, a new flag with value‘0’(zero) is created.

**set_flag**(`+Key, Value`)  
Set flag `Key` to `Value`. Value must be an atom, small (64-bit) integer or float.

**flag**(`+Key, -Old, +New`)  
True when `Old` is the current value of the flag `Key` and the flag has been set to `New`. `New` can be an arithmetic expression. The update is *atomic*. This predicate can be used to create a *shared* global counter as illustrated in the example below.

``` code
next_id(Id) :-
    flag(my_id, Id, Id+1).
```

### 4.14.4 Tries

Tries (also called *digital tree*, *radix tree* or *prefix tree* maintain a mapping between a variant of a term (see [=@=/2](compare.html#=@=/2)) and a value. They have been introduced in SWI-Prolog 7.3.21 as part of the implementation of *tabling*. The current implementation is rather immature. In particular, the following limitations currently apply:

- Tries only offer partial thread-safety. Multiple threads may concurrently insert values in a trie. Concurrent deletion and concurrent non-deterministic access is not supported.
- Tries should not be modified while non-deterministic predicates such as [trie_gen/3](db.html#trie_gen/3) are running on the trie.
- Terms cannot be *cyclic*. Possibly this will not change because cyclic terms can only be supported after creating a canonical form of the term.
- Starting with SWI-Prolog 9.3.23, tries support attributed variables both in for keys and values. Tries holding attributed variables can also be *compiled* using [trie_gen_compiled/2](db.html#trie_gen_compiled/2) or [trie_gen_compiled/3](db.html#trie_gen_compiled/3).

**We give the definition of these predicates for reference and debugging tabled predicates. Future versions are likely to get a more stable and safer implementation. The API to tries should not be considered stable.**

**trie_new**(`-Trie`)  
Create a new trie and unify `Trie` with a handle to the trie. The trie handle is a *blob*. Tries are subject to atom garbage collection.

**trie_destroy**(`+Trie`)  
Destroy `Trie`. This removes all nodes from the trie and causes further access to `Trie` to raise an existence_error exception. The handle itself is reclaimed by atom garbage collection.

\[semidet\]**is_trie**(`@Trie`)  
True when `Trie` is a trie object. See also [current_trie/1](db.html#current_trie/1).

\[nondet\]**current_trie**(`-Trie`)  
True if `Trie` is a currently existing trie. As this enumerates and then filters all known atoms this predicate is slow and should only be used for debugging purposes. See also [is_trie/1](db.html#is_trie/1).

**trie_insert**(`+Trie, +Key`)  
Insert the term `Key` into `Trie`. If `Key` is already part of `Trie` the predicates *fails* silently. This is the same as [trie_insert/3](db.html#trie_insert/3), but using a fixed reserved `Value`.

**trie_insert**(`+Trie, +Key, +Value`)  
Insert the term `Key` into `Trie` and associate it with `Value`. `Value` can be any term. If `Key`-`Value` is already part of `Trie`, the predicates *fails* silently. If `Key` is in `Trie` associated with a different value, a `permission_error` is raised.

**trie_update**(`+Trie, +Key, +Value`)  
As [trie_insert/3](db.html#trie_insert/3), but if `Key` is in `Trie`, its associated value is *updated*.

**trie_insert**(`+Trie, +Term, +Value, -Handle`)  
As [trie_insert/3](db.html#trie_insert/3), returning a handle to the trie node. This predicate is currently unsafe as `Handle` is an integer used to encode a pointer. It was used to implement a pure Prolog version of the `library(tabling)` library.

**trie_delete**(`+Trie, +Key, ?Value`)  
Delete `Key` from `Trie` if the value associated with `Key` unifies with `Value`.

**trie_lookup**(`+Trie, +Key, -Value`)  
True if the term `Key` is in `Trie` and associated with `Value`.

**trie_term**(`+Handle, -Term`)  
True when `Term` is a copy of the term associated with `Handle`. The result is undefined (including crashes) if `Handle` is not a handle returned by trie_insert_new/3 or the node has been removed afterwards.

\[nondet\]**trie_gen**(`+Trie, ?Key`)  
True when `Key` is a member of `Trie`. See also [trie_gen_compiled/2](db.html#trie_gen_compiled/2).

\[nondet\]**trie_gen**(`+Trie, ?Key, -Value`)  
True when `Key` is associated with `Value` in `Trie`. Backtracking retrieves all pairs. Currently scans the entire trie, even if `Key` is partly known. Currently unsafe if `Trie` is modified while the values are being enumerated. See also [trie_gen_compiled/3](db.html#trie_gen_compiled/3).

\[nondet\]**trie_gen_compiled**(`+Trie, ?Key`)  
\[nondet\]**trie_gen_compiled**(`+Trie, ?Key, -Value`)  
Similar to [trie_gen/3](db.html#trie_gen/3), but uses a *compiled* representation of `Trie`. The compiled representation is created lazily and manipulations of the trie (insert, delete) invalidate the current compiled representation. The compiled representation generates answers faster and, as it runs on a snapshot of the trie, is immune to concurrent modifications of the trie. This predicate is used to generate answers from *answer tries* as used for tabled execution. See [section 7](tabling.html#sec:7).

\[nondet\]**trie_property**(`?Trie, ?Property`)  
True if `Trie` exists with `Property`. Intended for debugging and statistical purposes. Retrieving some of these properties visit all nodes of the trie. Defined properties are

**value_count**(`-Count`)  
Number of key-value pairs in the trie.

**node_count**(`-Count`)  
Number of nodes in the trie.

**size**(`-Bytes`)  
Required storage space of the trie.

**compiled_size**(`-Bytes`)  
Required storage space for the compiled representation as used by [trie_gen_compiled/2](db.html#trie_gen_compiled/2),3.

**hashed**(`-Count`)  
Number of nodes that use a hashed index to its children.

**lookup_count**(`-Count`)  
Number of [trie_lookup/3](db.html#trie_lookup/3) calls (only when compiled with `O_TRIE_STATS`).

**gen_call_count**(`-Count`)  
Number of [trie_gen/3](db.html#trie_gen/3) calls (only when compiled with `O_TRIE_STATS`).

**wait**(`-Count`)  
Number of times a thread waited on this trie for another thread to complete it (shared tabling, only when compiled with `O_TRIE_STATS`).

**deadlock**(`-Count`)  
Number of times this trie was part of a deadlock and its completion was abandoned (shared tabling, only when compiled with `O_TRIE_STATS`).

In addition, a number of additional properties are defined on *answer tries*.

**invalidated**(`-Count`)  
Number of times the trie was invalidated (incremental tabling).

**reevaluated**(`-Count`)  
Number of times the trie was re-evaluated (incremental tabling).

**idg_affected_count**(`-Count`)  
Number of answer tries affected by this one (incremental tabling).

**idg_dependent_count**(`-Count`)  
Number of answer tries this one depends on (incremental tabling).

**idg_size**(`-Bytes`)  
Number of bytes in the IDG node representation.
