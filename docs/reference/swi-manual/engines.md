
# 11 Coroutining using Prolog engines

Where the term *coroutine* in Prolog typically refer to hooks triggered by *attributed variables* ([section 8.1](attvar.html#sec:8.1)), SWI-Prolog provides two other forms of coroutines. Delimited continuations (see [section 4.9](delcont.html#sec:4.9)) allow creating coroutines that run in the same Prolog engine by capturing and restarting the *continuation*. This section discusses *engines*, also known as *interactors*. The idea was developed by Paul Tarau [Tarau, 2011](Bibliography.html#DBLP:conf/coordination/Tarau11). The API described in this chapter has been established together with Paul Tarau and Paulo Moura.

Engines are closely related to *threads* ([section 10](threads.html#sec:10)). An engine is a Prolog virtual machine that has its own stacks and (virtual) machine state. Unlike normal Prolog threads though, they are not associated with an operating system thread. Instead, you *ask* an engine for a next answer ([engine_next/2](engine-predicates.html#engine_next/2)). Asking an engine for the next answer attaches the engine to the calling operating system thread and cause it to run until the engine calls [engine_yield/1](engine-predicates.html#engine_yield/1) or its associated goal completes with an answer, failure or an exception. After the engine yields or completes, it is detached from the operating system thread and the answer term is made available to the calling thread. Communicating with an engine is similar to communicating with a Prolog system though the terminal. In this sense engines are related to *Pengines* as provided by library `library(pengines)`, but where Pengines aim primarily at accessing Prolog engines through the network, engines are in-process entities.

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[11.1 Examples using engines](engine-examples.html)

[11.1.1 Aggregation using engines](engine-examples.html#sec:11.1.1)

[11.1.2 State accumulation using engines](engine-examples.html#sec:11.1.2)

[11.1.3 Scalable many-agent applications](engine-examples.html#sec:11.1.3)

[11.2 Engine resource usage](engine-resources.html)

[11.3 Engine predicate reference](engine-predicates.html)
