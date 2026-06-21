
## 11.1 Examples using engines

We introduce engines by describing application areas and providing simple example programs. The predicates are defined in [section 11.3](engine-predicates.html#sec:11.3). We identify the following application areas for engines.

1.  Aggregating solutions from one or more goals. See [section 11.1.1](engine-examples.html#sec:11.1.1).
2.  Access the terms produced in *forward execution* through backtracking without collecting all of them first. [Section 11.1.1](engine-examples.html#sec:11.1.1) illustrates this as well.
3.  State accumulation and sharing. See [section 11.1.2](engine-examples.html#sec:11.1.2).
4.  Scalable many-agent applications. See [section 11.1.3](engine-examples.html#sec:11.1.3).

### 11.1.1 Aggregation using engines

Engines can be used to reason about solutions produced by a goal through backtracking. In this scenario we create an engine with the goal we wish to backtrack over and we enumerate all its solution using [engine_next/2](engine-predicates.html#engine_next/2). This usage scenario competes with the all solution predicates ([findall/3](allsolutions.html#findall/3), [bagof/3](allsolutions.html#bagof/3), etc.) and the predicates from library `library(aggregate)`. Below we implement [findall/3](allsolutions.html#findall/3) using engines.

``` code
findall(Templ, Goal, List) :-
        setup_call_cleanup(
            engine_create(Templ, Goal, E),
            get_answers(E, List),
            engine_destroy(E)).

get_answers(E, [H|T]) :-
        engine_next(E, H), !,
        get_answers(E, T).
get_answers(_, []).
```

The above is not a particularly attractive alternative for the built-in [findall/3](allsolutions.html#findall/3). It is mostly slower due to time required to create and destroy the engine as well as the (currently^(210The current implementation of engines is built on top of primitives that are not optimal for the engine use case. There is considerable opportunity to reduce the overhead.)) higher overhead of copying terms between engines than the overhead required by the dedicated representation used by [findall/3](allsolutions.html#findall/3).

It gets more interesting if we wish to combine answers from multiple backtracking predicates. Assume we have two predicates that, on backtracking, return ordered solutions and we wish to merge the two answer streams into a single ordered stream of answers. The solution in classical Prolog is below. It collects both answer sets, merges them using ordered set merging and extract the answers. The code is clean and short, but it doesn't produce any answers before both generators are fully enumerated and it uses memory that is proportional to the combined set of answers.

``` code
:- meta_predicate merge(?,0, ?,0, -).

merge_answers(T1,G1, T2,G2, A) :-
        findall(T1, G1, L1),
        findall(T2, G2, L2),
        ord_union(L1, L2, Ordered),
        member(A, Ordered).
```

We can achieve the same using engines. We create two engines to generate the solutions to both our generators. Now, we can ask both for an answer, put the smallest in the answer set and ask the engine that created the smallest for its next answer, etc. This way we can create an ordered list of answers as above, but now without creating intermediate lists. We can avoid creating the intermediate list by introducing a third engine that controls the two generators and *yields* the answers rather than putting them in a list. This is a general example of turning a predicate that builds a set of terms into a non-deterministic predicate that produces the results on backtracking. The full code is below. Merging the answers of two generators, each generating a set of 10,000 integers is currently about 20% slower than the code above, but the engine-based solution runs in constant space and generates the first solution immediately after both our generators have produced their first solution.

``` code
:- meta_predicate merge(?,0, ?,0, -).

merge(T1,G1, T2,G2, A) :-
        engine_create(A, merge(T1,G1, T2,G2), E),
        repeat,
            (   engine_next(E, A)
            ->  true
            ;   !, fail
            ).

merge(T1,G1, T2,G2) :-
        engine_create(T1, G1, E1),
        engine_create(T2, G2, E2),
        (   engine_next(E1, S1)
        ->  (   engine_next(E2, S2)
            ->  order_solutions(S1, S2, E1, E2)
            ;   yield_remaining(S1, E1)
            )
        ;   engine_next(E2, S2),
            yield_remaining(S2, E2)
        ).

order_solutions(S1, S2, E1, E2) :- !,
        (   S1 @=< S2
        ->  engine_yield(S1),
            (   engine_next(E1, S11)
            ->  order_solutions(S11, S2, E1, E2)
            ;   yield_remaining(S2, E2)
            )
        ;   engine_yield(S2),
            (   engine_next(E2, S21)
            ->  order_solutions(S1, S21, E1, E2)
            ;   yield_remaining(S1, E1)
            )
        ).

yield_remaining(S, E) :-
        engine_yield(S),
        engine_next(E, S1),
        yield_remaining(S1, E).
```

### 11.1.2 State accumulation using engines

Applications that need to manage a state can do so by passing the state around in an additional argument, storing it in a global variable or update it in the dynamic database using [assertz/1](db.html#assertz/1) and [retract/1](db.html#retract/1). Both using an additional argument and a global variable (see [b_setval/2](gvar.html#b_setval/2)), make the state subject to backtracking. This may or may not be desirable. If having a state is that subject to backtracking is required, using an additional argument or backtrackable global variable is the right approach. Otherwise, non-backtrackable global variables ([nb_setval/2](gvar.html#nb_setval/2)) and dynamic database come into the picture, where global variables are always local to a thread and the dynamic database may or may not be shared between threads (see [thread_local/1](threadcom.html#thread_local/1)).

Engines bring an alternative that packages a state inside the engine where it is typically represented in a (threaded) Prolog variable. The state may be updated, while controlled backtracking to a previous state belongs to the possibilities. It can be accessed and updated by anyone with access to the engines’handle. Using an engine to represent state has the following advantages:

- The programming style needed inside the engine is much more‘Prolog friendly’, using [engine_fetch/1](engine-predicates.html#engine_fetch/1) to read a request and [engine_yield/1](engine-predicates.html#engine_yield/1) to reply to it.
- The state is packaged and subject to (atom) garbage collection.
- The state may be accessed from multiple threads. Access to the state is synchronized without the need for explicit locks.

The example below implements a shared priority heap based on library `library(heaps)`. The predicate update_heap/1 shows the typical update loop for maintaining state inside an engine: fetch a command, update the state, yield with the reply and call the updater recursively. The update step is guarded against failure. For robustness one may also guard it against exceptions using [catch/3](exception.html#catch/3). Note that heap_get/3 passes the `Priority` and `Key` it wishes to delete from the heap such that if the unification fails, the heap remains unchanged.

The resulting heap is a global object with either a named or anonymous handle that evolves independently from the Prolog thread(s) that access it. If the heap is anonymous, it is subject to (atom) garbage collection.

``` code
:- use_module(library(heaps)).

create_heap(E) :-
        empty_heap(H),
        engine_create(_, update_heap(H), E).

update_heap(H) :-
        engine_fetch(Command),
        (   update_heap(Command, Reply, H, H1)
        ->  true
        ;   H1 = H,
            Reply = false
        ),
        engine_yield(Reply),
        update_heap(H1).

update_heap(add(Priority, Key), true, H0, H) :-
        add_to_heap(H0, Priority, Key, H).
update_heap(get(Priority, Key), Priority-Key, H0, H) :-
        get_from_heap(H0, Priority, Key, H).

heap_add(Priority, Key, E) :-
        engine_post(E, add(Priority, Key), true).

heap_get(Priority, Key, E) :-
        engine_post(E, get(Priority, Key), Priority-Key).
```

### 11.1.3 Scalable many-agent applications

The final application area we touch are agent systems were we wish to capture an agent in a Prolog goal. Such systems can be implemented using threads (see [section 10](threads.html#sec:10)) that use [thread_send_message/2](threadcom.html#thread_send_message/2) and [thread_get_message/1](threadcom.html#thread_get_message/1) to communicate. The main problem is that each thread is associated by an operating system thread. OS threads are, depending on the OS, relatively expensive. Scalability of this design typically ends, depending on OS and hardware, somewhere between 1,000 and 100,000 agents.

Engines provide an alternative. A detached Prolog engine currently requires approximately 20 Kbytes memory on 64 bit hardware, growing with the size of the Prolog stacks. The Prolog stacks may be minimised by calling [garbage_collect/0](memory.html#garbage_collect/0) followed by [trim_stacks/0](memory.html#trim_stacks/0), providing a *deep sleep* mode. The set of agents, each represented by an engine can be controlled by a static or dynamic pool of threads. Scheduling the execution of agents and their communication is completely open and can be optimised to satisfy the requirements of the application.

> This section needs an example. Preferably something that fits on one page and would not scale using threads. Engines might work nice to implement *Antrank: An ant colony algorithm for ranking web pages*.^(211[http://www.ijettcs.org/Volume3Issue2/IJETTCS-2014-04-23-113.pdf](http://www.ijettcs.org/Volume3Issue2/IJETTCS-2014-04-23-113.pdf))
