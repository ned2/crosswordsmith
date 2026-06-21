Pengines: Web Logic Programming Made Easy

Torbjörn Lager  
University of Gothenburg  
Sweden  
E-mail: [lager@ling.gu.se](mailto:lager@ling.gu.se)  
Jan Wielemaker  
VU University Amsterdam  
The Netherlands  
E-mail: [J.Wielemaker@vu.nl](mailto:J.Wielemaker@vu.nl)

Abstract

Pengines is short for Prolog Engines. The pengines package greatly simplifies (1) developing JavaScript based web-applications that must talk to a Prolog server and (2) realise distributed programming in Prolog by providing RPC (*Remote Procedure Calling*) over HTTP.

See also [http://www.swi-prolog.org/pengines](http://www.swi-prolog.org/pengines).

# Table of Contents

[1 An overview of Pengines](#sec:1)

[1.1 Pengine references](#sec:1.1)

[1.2 Pengine by examples](#sec:1.2)

[1.3 Making predicates available to clients](#sec:1.3)

[1.4 Mapping Prolog terms into JSON](#sec:1.4)

[1.5 Pengine settings](#sec:1.5)

[2 Pengine libraries](#sec:2)

[2.1 library(pengines): Pengines: Web Logic Programming Made Easy](#sec:2.1)

[2.2 library(term_to_json)](#sec:2.2)

## 1 An overview of Pengines

This package provides a powerful high-level programming abstraction implemented on top of SWI-Prolog's thread predicates \[1\] and its HTTP client and server libraries \[2\]. The package makes it easy to create and query *Prolog engines* (or *Pengines* for short), over HTTP, from an ordinary Prolog thread, from a pengine, or from JavaScript running in a web client. Querying follows Prolog's default one-tuple-at-a-time generation of solutions. I/O is also supported.

Possible applications abound, but in particular three kinds of applications stick out: 1) The package provides us with a useful point of departure for the design and implementation of more advanced Prolog-based agent programming platforms, 2) it suggests an elegant and very straightforward approach to the building of a Semantic Web which is Prolog-based in a very radical sense, and, 3) it constitutes an ideal way to interface Prolog with JavaScript, the programming language most commonly available in web browsers.

A pengine is comprised of:

- A Prolog thread
- A dynamic clause database, private to the pengine, into which other processes may assert clauses. These clauses reside in the module `pengine_sandbox`.
- A message queue for incoming requests
- A message queue for outgoing responses

Everything needed to work with pengines is included in the package, including a JavaScript library for creating and interacting with pengines from a web client. However, the web server (in the file `examples/server.pl`) should only be regarded as a minimal example.

Underlying the design of the package is a careful analysis of the conversations taking place between Prolog and a user (which could be a human or another piece of software). Such conversations follow a communication protocol that we refer to as the Prolog Transport Protocol (PLTP). The protocol has been modelled by means of so called *communicating finite-state machines* \[3\]. A slight modification of the protocol -- referred to as PLTP(HTTP) -- enables us to synchronize it with HTTP. The diagram below depicts the communicating finite-state machines for PLTP(HTTP) and HTTP. Labels in bold indicate requests, and labels with a slash in front indicate responses.

![](pltpsynch.png)

The diagram below depicts a PLTP run (on the right) corresponding to a user's interaction with Prolog (on the left).‘1234’is the Pengine's identifier, which is a UUID in the actual implementation.

![](pltpruncolour.png)

As for the relations between pengines, and for the time being, we have opted for a simple *master-slave architecture*. Once the master/slave relationships are established, the direction of control is always from the master to the slaves. One or more pengines can be *orchestrated* by a common master which can be an ordinary Prolog thread, another pengine, or a JavaScript process. A slave is always a pengine, running either locally or remotely with respect to its master. Subject to a setting, slaves are also dependent on their masters in the sense that if a master terminates, so do its slaves. (Note that in the source code we often use the term *parent* instead of *master* and *child* instead of *slave*. That is, we treat *parent/child* as synonymous to *master/slave*.)

![](penarch.png)

The transport format is different depending on the nature of the master. If the master is a JavaScript process, it will (by default) formulate its requests using Prolog syntax, and get responses back as Prologs terms encoded in JSON. If the master is a Prolog process (a Prolog thread or a pengine) it will (again only by default) get responses back as Prolog terms.

Most of the pengine predicates are deterministic, yet they can control one or more pengines solving possibly non-deterministic queries. But the package also offers a number of non-deterministic predicates, built on top of the deterministic ones, that can solve queries "the Prolog way", binding query variables in the process, backtracking for more solutions. Of these predicates, [pengine_rpc/3](#pengine_rpc/3) is the most important. By means of [pengine_rpc/3](#pengine_rpc/3) a pengine running in a pengine server A can call and try to solve a query in the context of another pengine server B, taking advantage of the data being offered by B, just as if the data was local to A. Thus, in theory, a Prolog program, be it a pure Horn clause theory or not, can be as big as the Web. This is something that should make us think about a *Semantic Web*, especially when we consider the excellent fit between the Pengine library and SWI-Prolog's Semantic Web Library \[4\]. Adding Pengines functionality to the Cliopatria platform \[5\] is straightforward.

A note about safety: Because PLTP is layered on top of HTTP, it may utilize any standard HTTP security feature, such as HTTP authentication or SSL. Moreover, subject to a setting, the library uses safe_goal/1 \[6\], which determines whether it is safe for a slave pengine to try to solve queries asked by a master.

### 1.1 Pengine references

1.  [http://www.swi-prolog.org/pldoc/man?section=threads](http://www.swi-prolog.org/pldoc/man?section=threads)
2.  [http://www.swi-prolog.org/pldoc/package/http.html](http://www.swi-prolog.org/pldoc/package/http.html)
3.  D. Brand and P. Zafiropulo. On communicating finite-state machines. *Journal of the ACM*, 30(2):323-342, 1983.
4.  [http://www.swi-prolog.org/pldoc/package/semweb.html](http://www.swi-prolog.org/pldoc/package/semweb.html)
5.  [http://cliopatria.swi-prolog.org/home](http://cliopatria.swi-prolog.org/home)
6.  [http://www.swi-prolog.org/pldoc/doc/home/vnc/prolog/lib/swipl/library/sandbox.pl](http://www.swi-prolog.org/pldoc/doc/home/vnc/prolog/lib/swipl/library/sandbox.pl)

### 1.2 Pengine by examples

In this example we load the pengines library and use [pengine_create/1](#pengine_create/1) to create a slave pengine in a remote pengine server, and inject a number of clauses in it. We then use [pengine_event_loop/2](#pengine_event_loop/2) to start an event loop that listens for three kinds of event terms. Running `main/0` will write the terms `q(a)`, `q(b)` and `q(c)` to standard output. Using [pengine_ask/3](#pengine_ask/3) with the option `template(X)` would instead produce the output `a`, `b` and `c`.

``` code
:- use_module(library(pengines)).

main :-
    pengine_create([
        server('https://pengines.swi-prolog.org'),
        src_text("
            q(X) :- p(X).
            p(a). p(b). p(c).
        ")
    ]),
    pengine_event_loop(handle, []).

handle(create(ID, _)) :-
    pengine_ask(ID, q(_X), []).
handle(success(_ID, [X], false)) :-
    writeln(X).
handle(success(ID, [X], true)) :-
    writeln(X),
    pengine_next(ID, []).
```

Here is another example, showing how to create and interact with a pengine from JavaScript in a way that seems ideal for Prolog programmers and JavaScript programmers alike. Loading the page brings up the browser's prompt dialog, waits for the user's input, and writes that input in the browser window. If the input was’stop’, it stops there, else it repeats. Note that I/O works as expected. All we need to do is to use pengine_input/1 instead of read/1 and [pengine_output/1](#pengine_output/1) instead of write/1.

**See Also:**

- [pengines.js documentation](https://pengines.swi-prolog.org/docs/documentation.html)

``` code
<html lang="en">
    <head>
        <title>Pengine Example</title>
    </head>
    <body>
            <h1>Pengine Example</h1>
        <div id="out"></div>

                <script src="https://ajax.googleapis.com/ajax/libs/jquery/2.0.3/jquery.min.js"></script>
        <script src="https://pengines.swi-prolog.org/pengine/pengines.js"></script>

        <script type="text/x-prolog">

            main :-
                repeat,
                pengine_input(X),
                pengine_output(X),
                X == stop.

        </script>
        <script>
            var pengine = new Pengine({
                oncreate: handleCreate,
                onprompt: handlePrompt,
                onoutput: handleOutput
            });
            function handleCreate() {
                pengine.ask('main');
            }
            function handlePrompt() {
                pengine.input(prompt(this.data));
            }
            function handleOutput() {
                $('#out').html(this.data);
            }
        </script>
    </body>
</html>
```

Our third example shows that a non-deterministic predicate can be called remotely by means of [pengine_rpc/2](#pengine_rpc/2), yet behave exactly as if called locally:

``` code
?- use_module(library(pengines)).

?- member(X, [a, b, c, d]),
   pengine_rpc('https://pengines.swi-prolog.org', p(X), [
       src_list([p(b), p(c), p(d), p(e)])
   ]),
   member(X, [c, d, e, f]).
X = c ;
X = d.

?-
```

### 1.3 Making predicates available to clients

The code sent to a pengine is executed in the context of the module `pengine_sandbox` and the safety of goals is validated using safe_goal/1 prior to execution. Any pengine has access to the safe predicates defined in `library(sandbox)`. If a server wishes to extend the set of predicates, it must:

1.  Define one or more modules that export the desired additional predicates.

2.  Makes this code available to the sandbox using the call below, assuming that the additional predicates are defined in the Prolog module file `domain_predicates.pl`

    ``` code
    :- use_module(pengine_sandbox:domain_predicates).
    ```

3.  Register **safe** foreign predicates with `library(sandbox)`, i.e., predicates that do not have side effects such as accessing the file system, load foreign extensions, define other predicates outside the sandbox environment, etc.

    Note that the safety of Prolog predicate can typically be proven by `library(sandbox)`. This may not be the case if untracktable forms of meta-calling are used. In this case it is adviced to avoid such code. If this is not possible, the code must be carefully reviewed by hand and of proven to be safe it may be registered with the sandbox library.

For example, basic RDF access can be granted to pengines using the code below. Please **study the sandboxing code carefully before adding declarations**.

``` code
:- use_module(pengine_sandbox:library(semweb/rdf_db)).
:- use_module(library(sandbox)).

:- multifile sandbox:safe_primitive/1.

sandbox:safe_primitive(rdf_db:rdf(_,_,_)).
```

### 1.4 Mapping Prolog terms into JSON

In Prolog, solutions to queries are given as bindings which map variable names into Prolog terms. A programmer using Pengines in a JavaScript evironment needs to understand how bindings are converted into JSON. For example, the programmer needs to understand that the second solution to `append(Xs, Ys, [a,b,c])` is given by the bindings `['Xs'=[a],'Ys'=[b,c]]` and that these binding can be represented in JSON as `{"Xs":["a"], "Ys":["b","c"]}`.

Pengines defines the following mapping between ground Prolog terms and JSON.

- A Prolog atom is mapped to a JSON string.
- A Prolog number is mapped to a JSON number.
- A Prolog list is mapped to a JSON array.
- The Prolog terms `@(true)` and `@(false)` are mapped to the JSON constants `true` and `false`, respectively.
- The Prolog term `@(null)` is mapped to the JSON constant `null`.
- A Prolog term `json(NameValueList)`, where `NameValueList` is a list of `Name=Value` pairs, is mapped to a JSON object.
- Any other complex Prolog term `T` is mapped to a JSON object of the form `{"functor": F, "args": A}` where `F` is a string representing the functor of `T` and `A` is the list of JSON values representing `T`s arguments.

### 1.5 Pengine settings

Settings currently recognized by the Pengines library:

> |  |  |  |  |
> |----|----|----|----|
> | **Name** | **Type** | **Default** | **Description** |
> | max_session_pengines | integer | 1 | Maximum number of pengines a client can create. -1 is infinite |
> | time_limit | number | 60 | Maximum time between output (in seconds) |
> | allow_from | `list(atom)` | \[\*\] | Specify allowed IP addresses |
> | deny_from | `list(atom)` | `[]` | Specify denied IP addresses. Applied after `allow_from`. |

## 2 Pengine libraries

### 2.1 library(pengines): Pengines: Web Logic Programming Made Easy

author  
Torbjörn Lager and Jan Wielemaker

The `library(pengines)` provides an infrastructure for creating Prolog engines in a (remote) pengine server and accessing these engines either from Prolog or JavaScript.

\[det\]**pengine_create**(`:Options`)  
Creates a new pengine. Valid options are:

**id**(`-ID`)  
`ID` gets instantiated to the id of the created pengine. `ID` is atomic.

**alias**(`+Name`)  
The pengine is named `Name` (an atom). A slave pengine (child) can subsequently be referred to by this name.

**application**(`+Application`)  
`Application` in which the pengine runs. See [pengine_application/1](#pengine_application/1).

**server**(`+URL`)  
The pengine will run in (and in the Prolog context of) the pengine server located at `URL`.

**src_list**(`+List_of_clauses`)  
Inject a list of Prolog clauses into the pengine.

**src_text**(`+Atom_or_string`)  
Inject the clauses specified by a source text into the pengine.

**src_url**(`+URL`)  
Inject the clauses specified in the file located at `URL` into the pengine.

**src_predicates**(`+List`)  
Send the local predicates denoted by `List` to the remote pengine. `List` is a list of predicate indicators.

Remaining options are passed to http_open/3 (meaningful only for non-local pengines) and thread_create/3. Note that for thread_create/3 only options changing the stack-sizes can be used. In particular, do not pass the detached or alias options..

Successful creation of a pengine will return an *event term* of the following form:

**create**(`ID, Term`)  
`ID` is the id of the pengine that was created. `Term` is not used at the moment.

An error will be returned if the pengine could not be created:

**error**(`ID, Term`)  
`ID` is invalid, since no pengine was created. `Term` is the exception's error term.

\[det\]**pengine_ask**(`+NameOrID, @Query, +Options`)  
Asks pengine `NameOrID` a query `Query`.

`Options` is a list of options:

**template**(`+Template`)  
`Template` is a variable (or a term containing variables) shared with the query. By default, the template is identical to the query.

**chunk**(`+IntegerOrFalse`)  
Retrieve solutions in chunks of Integer rather than one by one. 1 means no chunking (default). Other integers indicate the maximum number of solutions to retrieve in one chunk. If `false`, the Pengine goal is not executed using findall/3 and friends and we do not backtrack immediately over the goal. As a result, changes to backtrackable global state are retained. This is similar that using `set_prolog_flag(toplevel_mode, recursive)`.

**bindings**(`+Bindings`)  
Sets the global variable’\$variable_names’to a list of `Name = Var` terms, providing access to the actual variable names.

Any remaining options are passed to pengine_send/3.

Note that the predicate [pengine_ask/3](#pengine_ask/3) is deterministic, even for queries that have more than one solution. Also, the variables in `Query` will not be bound. Instead, results will be returned in the form of *event terms*.

**success**(`ID, Terms, Projection, Time, More`)  
`ID` is the id of the pengine that succeeded in solving the query. `Terms` is a list holding instantiations of `Template`. `Projection` is a list of variable names that should be displayed. `Time` is the CPU time used to produce the results and finally, `More` is either `true` or `false`, indicating whether we can expect the pengine to be able to return more solutions or not, would we call [pengine_next/2](#pengine_next/2).

**failure**(`ID`)  
`ID` is the id of the pengine that failed for lack of a solutions.

**error**(`ID, Term`)  
`ID` is the id of the pengine throwing the exception. `Term` is the exception's error term.

**output**(`ID, Term`)  
`ID` is the id of a pengine running the query that called [pengine_output/1](#pengine_output/1). `Term` is the term that was passed in the first argument of [pengine_output/1](#pengine_output/1) when it was called.

**prompt**(`ID, Term`)  
`ID` is the id of the pengine that called [pengine_input/2](#pengine_input/2) and `Term` is the prompt.

Defined in terms of pengine_send/3, like so:

``` code
pengine_ask(ID, Query, Options) :-
    partition(pengine_ask_option, Options, AskOptions, SendOptions),
    pengine_send(ID, ask(Query, AskOptions), SendOptions).
```

\[det\]**pengine_next**(`+NameOrID, +Options`)  
Asks pengine `NameOrID` for the next solution to a query started by [pengine_ask/3](#pengine_ask/3). Defined options are:

**chunk**(`+Count`)  
Modify the chunk-size to `Count` before asking the next set of solutions. This may not be used if the goal was started with `chunk(false)`.

Remaining options are passed to pengine_send/3. The result of re-executing the current goal is returned to the caller's message queue in the form of *event terms*.

**success**(`ID, Terms, Projection, Time, More`)  
See [pengine_ask/3](#pengine_ask/3).

**failure**(`ID`)  
`ID` is the id of the pengine that failed for lack of more solutions.

**error**(`ID, Term`)  
`ID` is the id of the pengine throwing the exception. `Term` is the exception's error term.

**output**(`ID, Term`)  
`ID` is the id of a pengine running the query that called [pengine_output/1](#pengine_output/1). `Term` is the term that was passed in the first argument of [pengine_output/1](#pengine_output/1) when it was called.

**prompt**(`ID, Term`)  
`ID` is the id of the pengine that called [pengine_input/2](#pengine_input/2) and `Term` is the prompt.

Defined in terms of pengine_send/3, as follows:

``` code
pengine_next(ID, Options) :-
    pengine_send(ID, next, Options).
```

\[det\]**pengine_stop**(`+NameOrID, +Options`)  
Tells pengine `NameOrID` to stop looking for more solutions to a query started by [pengine_ask/3](#pengine_ask/3). `Options` are passed to pengine_send/3.

Defined in terms of pengine_send/3, like so:

``` code
pengine_stop(ID, Options) :-
    pengine_send(ID, stop, Options).
```

\[det\]**pengine_abort**(`+NameOrID`)  
Aborts the running query. The pengine goes back to state‘2’, waiting for new queries.

See also  
[pengine_destroy/1](#pengine_destroy/1).

\[det\]**pengine_destroy**(`+NameOrID`)  
\[det\]**pengine_destroy**(`+NameOrID, +Options`)  
Destroys the pengine `NameOrID`. With the option `force(true)`, the pengine is killed using abort/0 and [pengine_destroy/2](#pengine_destroy/2) succeeds.

\[det\]**pengine_self**(`-Id`)  
True if the current thread is a pengine with `Id`.

\[det\]**pengine_application**(`+Application`)  
Directive that must be used to declare a pengine application module. The module must not be associated to any file. The default application is `pengine_sandbox`. The example below creates a new application `address_book` and imports the API defined in the module file `adress_book_api.pl` into the application.

``` code
:- pengine_application(address_book).
:- use_module(address_book:adress_book_api).
```

\[nondet\]**current_pengine_application**(`?Application`)  
True when `Application` is a currently defined application.

See also  
[pengine_application/1](#pengine_application/1)

\[nondet\]**pengine_property**(`?Pengine, ?Property`)  
True when `Property` is a property of the given `Pengine`. Enumerates all pengines that are known to the calling Prolog process. Defined properties are:

**self**(`ID`)  
Identifier of the pengine. This is the same as the first argument, and can be used to enumerate all known pengines.

**alias**(`Name`)  
`Name` is the alias name of the pengine, as provided through the `alias` option when creating the pengine.

**thread**(`Thread`)  
If the pengine is a local pengine, `Thread` is the Prolog thread identifier of the pengine.

**remote**(`Server`)  
If the pengine is remote, the URL of the server.

**application**(`Application`)  
`Pengine` runs the given application

**module**(`Module`)  
Temporary module used for running the `Pengine`.

**destroy**(`Destroy`)  
`Destroy` is `true` if the pengines is destroyed automatically after completing the query.

**parent**(`Queue`)  
Message queue to which the (local) pengine reports.

**source**(`?SourceID, ?Source`)  
`Source` is the source code with the given `SourceID`. May be present if the setting `debug_info` is present.

**detached**(`?Time`)  
`Pengine` was detached at `Time`.

\[det\]**pengine_output**(`+Term`)  
Sends `Term` to the parent pengine or thread.

\[det\]**pengine_debug**(`+Format, +Args`)  
Create a message using format/3 from `Format` and `Args` and send this to the client. The default JavaScript client will call `console.log(Message)` if there is a console. The predicate [pengine_rpc/3](#pengine_rpc/3) calls `debug(pengine(debug), '~w', [Message])`. The debug topic `pengine(debug)` is enabled by default.

See also  
\- debug/1 and nodebug/1 for controlling the `pengine(debug)` topic  
- format/2 for format specifications

\[det,multifile\]thread_pool:**create_pool**(`+Application`)  
On demand creation of a thread pool for a pengine application.

\[det\]**pengine_done**  
Called from the pengine thread `at_exit` option. Destroys *child* pengines using [pengine_destroy/1](#pengine_destroy/1). Cleaning up the Pengine is synchronised by the `pengine_done` mutex. See read_event/6.

\[semidet,multifile\]**prepare_module**(`+Module, +Application, +Options`)  
Hook, called to initialize the temporary private module that provides the working context of a pengine. This hook is executed by the pengine's thread. Preparing the source consists of three steps:

1.  Add `Application` as (first) default import module for `Module`
2.  Call this hook
3.  Compile the source provided by the the `src_text` and `src_url` options

|  |  |
|----|----|
| `Module` | is a new temporary module (see in_temporary_module/3) that may be (further) prepared by this hook. |
| `Application` | (also a module) associated to the pengine. |
| `Options` | is passed from the environment and should (currently) be ignored. |

\[semidet,multifile\]**prepare_goal**(`+Goal0, -Goal1, +Options`)  
Pre-preparation hook for running `Goal0`. The hook runs in the context of the pengine. Goal is the raw goal given to *ask*. The returned `Goal1` is subject to goal expansion (expand_goal/2) and sandbox validation (safe_goal/1) prior to execution. If this goal fails, `Goal0` is used for further processing.

|           |                                        |
|-----------|----------------------------------------|
| `Options` | provides the options as given to *ask* |

\[semidet,multifile\]**not_sandboxed**(`+User, +Application`)  
This hook is called to see whether the Pengine must be executed in a protected environment. It is only called after [authentication_hook/3](#authentication_hook/3) has confirmed the authentity of the current user. If this hook succeeds, both loading the code and executing the query is executed without enforcing sandbox security. Typically, one should:

1.  Provide a safe user authentication hook.
2.  Enable HTTPS in the server or put it behind an HTTPS proxy and ensure that the network between the proxy and the pengine server can be trusted.

\[det\]**pengine_pull_response**(`+Pengine, +Options`)  
Pulls a response (an event term) from the slave `Pengine` if `Pengine` is a remote process, else does nothing at all.

\[det\]**pengine_input**(`+Prompt, -Term`)  
Sends `Prompt` to the master (parent) pengine and waits for input. Note that `Prompt` may be any term, compound as well as atomic.

\[det\]**pengine_respond**(`+Pengine, +Input, +Options`)  
Sends a response in the form of the term `Input` to a slave (child) pengine that has prompted its master (parent) for input.

Defined in terms of pengine_send/3, as follows:

``` code
pengine_respond(Pengine, Input, Options) :-
    pengine_send(Pengine, input(Input), Options).
```

\[det\]**pengine_event_loop**(`:Closure, +Options`)  
Starts an event loop accepting event terms sent to the current pengine or thread. For each such event E, calls `ignore(call(Closure, E))`. A closure thus acts as a *handler* for the event. Some events are also treated specially:

**create**(`ID, Term`)  
The `ID` is placed in a list of active pengines.

**destroy**(`ID`)  
The `ID` is removed from the list of active pengines. When the last pengine `ID` is removed, the loop terminates.

**output**(`ID, Term`)  
The predicate [pengine_pull_response/2](#pengine_pull_response/2) is called.

Valid options are:

**autoforward**(`+To`)  
Forwards received event terms to slaves. `To` is either `all`, `all_but_sender` or a Prolog list of NameOrIDs. \[not yet implemented\]

\[nondet\]**pengine_rpc**(`+URL, +Query`)  
\[nondet\]**pengine_rpc**(`+URL, +Query, +Options`)  
Semantically equivalent to the sequence below, except that the query is executed in (and in the Prolog context of) the pengine server referred to by `URL`, rather than locally.

``` code
  copy_term_nat(Query, Copy),  % attributes are not copied to the server
  call(Copy),            % executed on server at URL
  Query = Copy.
```

Valid options are:

**chunk**(`+IntegerOrFalse`)  
Can be used to reduce the number of network roundtrips being made. See [pengine_ask/3](#pengine_ask/3).

**timeout**(`+Time`)  
Wait at most `Time` seconds for the next event from the server. The default is defined by the setting `pengines:time_limit`.

Remaining options (except the server option) are passed to [pengine_create/1](#pengine_create/1).

\[semidet,multifile\]**prompt**(`+ID, +Prompt, -Term`)  
Hook to handle [pengine_input/2](#pengine_input/2) from the remote pengine. If the hooks fails, [pengine_rpc/3](#pengine_rpc/3) calls read/1 using the current prompt.

\[semidet,multifile\]**output**(`+ID, +Term`)  
Hook to handle [pengine_output/1](#pengine_output/1) from the remote pengine. If the hook fails, it calls print/1 on `Term`.

\[det\]**portray_blob**(`+Blob, +Options`)  
Portray non-text blobs that may appear in output terms. Not really sure about that. Basically such terms need to be avoided as they are meaningless outside the process. The generated error is hard to debug though, so now we send them as `'$BLOB'(Type)`. Future versions may include more info, depending on `Type`.

\[semidet,multifile\]**write_result**(`+Lang, +Event, +Dict`)  
Hook that allows for different output formats. The core Pengines library supports `prolog` and various JSON dialects. The hook [event_to_json/3](#event_to_json/3) can be used to refine the JSON dialects. This hook must be used if a completely different output format is desired.

**add_error_details**(`+Error, +JSON0, -JSON`)  
Add format error code and location information to an error. Also used by `pengines_io.pl`.

\[semidet,multifile\]**event_to_json**(`+Event, -JSONTerm, +Lang`)  
Hook that translates a Pengine event structure into a term suitable for reply_json_dict/1, according to the language specification `Lang`. This can be used to massage general Prolog terms, notably associated with `success(ID, Bindings, Projection, Time, More)` and `output(ID, Term)` into a format suitable for processing at the client side.

\[semidet,multifile\]**authentication_hook**(`+Request, +Application, -User`)  
This hook is called from the =/pengine/create= HTTP handler to discover whether the server is accessed by an authorized user. It can react in three ways:

- Succeed, binding `User` to a ground term. The authentity of the user is available through [pengine_user/1](#pengine_user/1).
- Fail. The =/create= succeeds, but the pengine is not associated with a user.
- Throw an exception to prevent creation of the pengine. Two meaningful exceptions are:
  - `throw(http_reply(authorise(basic(Realm))))` Start a normal HTTP login challenge (reply 401)
  - `throw(http_reply(forbidden(Path)))`) Reject the request using a 403 repply.

See also  
http_authenticate/3 can be used to implement this hook using default HTTP authentication data.

\[semidet\]**pengine_user**(`-User`)  
True when the pengine was create by an HTTP request that authorized `User`.

See also  
[authentication_hook/3](#authentication_hook/3) can be used to extract authorization from the HTTP header.

\[det\]**pengine_event**(`?EventTerm`)  
\[det\]**pengine_event**(`?EventTerm, +Options`)  
Examines the pengine's event queue and if necessary blocks execution until a term that unifies to Term arrives in the queue. After a term from the queue has been unified to Term, the term is deleted from the queue.

Valid options are:

**timeout**(`+Time`)  
`Time` is a float or integer and specifies the maximum time to wait in seconds. If no event has arrived before the time is up `EventTerm` is bound to the atom `timeout`.

**listen**(`+Id`)  
Only listen to events from the pengine identified by `Id`.

### 2.2 library(term_to_json)

\[det\]**term_to_json**(`+Term, +Bindings, -JsonTerm`)  
\[det\]**term_to_json**(`+Term, -JsonTerm`)  
Convert any general Prolog term into a JSON term. Prolog lists are treated in a special way. Also, JSON terms are not converted. Mapping:

- Variable: `{"type":"var", "name":<string>}`
- Atom: `{"type":"atom", "value":<string>}`
- Integer: `{"type":"integer", "value":<integer>}`
- Float: `{"type":"float", "value":<float>}`
- List: JSON array
- Dict: a JSON object. Values are processed recursively. (the tag is ignored)
- `json([Key=Value, ...])`: a JSON object Values are processed recursively.
- compound: `{"type":"compound", "functor":<string>, "args":<array>}`

|  |  |
|----|----|
| `Bindings` | is a list of Name=Var terms for variables that get their name from the environment. |

# Index

?  
[add_error_details/3](#add_error_details/3)  
[authentication_hook/3](#authentication_hook/3)  
[current_pengine_application/1](#current_pengine_application/1)  
[event_to_json/3](#event_to_json/3)  
[not_sandboxed/2](#not_sandboxed/2)  
[output/2](#output/2)  
[pengine_abort/1](#pengine_abort/1)  
[pengine_application/1](#pengine_application/1)  
[pengine_ask/3](#pengine_ask/3)  
[pengine_create/1](#pengine_create/1)  
[pengine_debug/2](#pengine_debug/2)  
[pengine_destroy/1](#pengine_destroy/1)  
[pengine_destroy/2](#pengine_destroy/2)  
[pengine_done/0](#pengine_done/0)  
[pengine_event/1](#pengine_event/1)  
[pengine_event/2](#pengine_event/2)  
[pengine_event_loop/2](#pengine_event_loop/2)  
[pengine_input/2](#pengine_input/2)  
[pengine_next/2](#pengine_next/2)  
[pengine_output/1](#pengine_output/1)  
[pengine_property/2](#pengine_property/2)  
[pengine_pull_response/2](#pengine_pull_response/2)  
[pengine_respond/3](#pengine_respond/3)  
[pengine_rpc/2](#pengine_rpc/2)  
[pengine_rpc/3](#pengine_rpc/3)  
[pengine_self/1](#pengine_self/1)  
[pengine_stop/2](#pengine_stop/2)  
[pengine_user/1](#pengine_user/1)  
[portray_blob/2](#portray_blob/2)  
[prepare_goal/3](#prepare_goal/3)  
[prepare_module/3](#prepare_module/3)  
[prompt/3](#prompt/3)  
[term_to_json/2](#term_to_json/2)  
[term_to_json/3](#term_to_json/3)  
[thread_pool:create_pool/1](#thread_pool:create_pool/1)  
[write_result/3](#write_result/3)  
