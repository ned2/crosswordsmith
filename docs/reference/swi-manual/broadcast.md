
## A.5 library(broadcast): Broadcast and receive event notifications

The `library(broadcast)` library was invented to realise GUI applications consisting of stand-alone components that use the Prolog database for storing the application data. [Figure 11](broadcast.html#fig:broadcast) illustrates the flow of information using this design

![](broadcast.png)

**Figure 11 :** Information-flow using broadcasting service

The broadcasting service provides two services. Using the‘shout’service, an unknown number of agents may listen to the message and act. The broadcaster is not (directly) aware of the implications. Using the‘request’service, listening agents are asked for an answer one-by-one and the broadcaster is allowed to reject answers using normal Prolog failure.

Shouting is often used to inform about changes made to a common database. Other messages can be “save yourself” or “show this” .

Requesting is used to get information while the broadcaster is not aware who might be able to answer the question. For example “who is showing `X`?” .

**broadcast**(`+Term`)  
Broadcast `Term`. There are no limitations to `Term`, though being a global service, it is good practice to use a descriptive and unique principal functor. All associated goals are started and regardless of their success or failure, [broadcast/1](broadcast.html#broadcast/1) always succeeds. Exceptions are passed.

**broadcast_request**(`+Term`)  
Unlike [broadcast/1](broadcast.html#broadcast/1), this predicate stops if an associated goal succeeds. Backtracking causes it to try other listeners. A broadcast request is used to fetch information without knowing the identity of the agent providing it. C.f. “Is there someone who knows the age of John?” could be asked using

``` code
        ...,
        broadcast_request(age_of('John', Age)),
```

If there is an agent (*listener*) that registered an‘age-of’service and knows about the age of‘John’this question will be answered.

**listen**(`+Template, :Goal`)  
Register a *listen* channel. Whenever a term unifying `Template` is broadcasted, call `Goal`. The following example traps all broadcasted messages as a variable unifies to any message. It is commonly used to debug usage of the library.

``` code
?- listen(Term, (writeln(Term),fail)).
?- broadcast(hello(world)).
hello(world)
true.
```

**listen**(`+Listener, +Template, :Goal`)  
Declare `Listener` as the owner of the channel. Unlike a channel opened using [listen/2](broadcast.html#listen/2), channels that have an owner can terminate the channel. This is commonly used if an object is listening to broadcast messages. In the example below we define a‘name-item’displaying the name of an identifier represented by the predicate name_of/2.

``` code
:- pce_begin_class(name_item, text_item).

variable(id,    any,    get, "Id visualised").

initialise(NI, Id:any) :->
        name_of(Id, Name),
        send_super(NI, initialise, name, Name,
                   message(NI, set_name, @arg1)),
        send(NI, slot, id, Id),
        listen(NI, name_of(Id, Name),
               send(NI, selection, Name)).

unlink(NI) :->
        unlisten(NI),
        send_super(NI, unlink).

set_name(NI, Name:name) :->
        get(NI, id, Id),
        retractall(name_of(Id, _)),
        assert(name_of(Id, Name)),
        broadcast(name_of(Id, Name)).

:- pce_end_class.
```

**unlisten**(`+Listener`)  
Deregister all entries created with [listen/3](broadcast.html#listen/3) whose `Listener` unify.

**unlisten**(`+Listener, +Template`)  
Deregister all entries created with [listen/3](broadcast.html#listen/3) whose `Listener` and `Template` unify.

**unlisten**(`+Listener, +Template, :Goal`)  
Deregister all entries created with [listen/3](broadcast.html#listen/3) whose `Listener`, `Template` and `Goal` unify.

**listening**(`?Listener, ?Template, ?Goal`)  
Examine the current listeners. This predicate is useful for debugging purposes.
