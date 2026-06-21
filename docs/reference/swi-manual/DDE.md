
## 4.44 Windows DDE interface

The predicates in this section deal with MS-Windows‘Dynamic Data Exchange’or DDE protocol.^(168This interface is contributed by Don Dwiggins.) A Windows DDE conversation is a form of interprocess communication based on sending reserved window events between the communicating processes.

Failing DDE operations raise an error of the structure below, where `Operation` is the name of the (partial) operation that failed and `Message` is a translation of the operator error code. For some errors, `Context` provides additional comments.

``` code
        error(dde_error(Operation, Message), Context)
```

### 4.44.1 DDE client interface

The DDE client interface allows Prolog to talk to DDE server programs. We will demonstrate the use of the DDE interface using the Windows PROGMAN (Program Manager) application:

``` code
1 ?- open_dde_conversation(progman, progman, C).

C = 0
2 ?- dde_request(0, groups, X)

--> Unifies X with description of groups

3 ?- dde_execute(0, '[CreateGroup("DDE Demo")]').
true.

4 ?- close_dde_conversation(0).
true.
```

For details on interacting with **progman**, use the SDK online manual section on the Shell DDE interface. See also the Prolog `library(progman)`, which may be used to write simple Windows setup scripts in Prolog.

**open_dde_conversation**(`+Service, +Topic, -Handle`)  
Open a conversation with a server supporting the given service name and topic (atoms). If successful, `Handle` may be used to send transactions to the server. If no willing server is found this predicate fails silently.

**close_dde_conversation**(`+Handle`)  
Close the conversation associated with `Handle`. All opened conversations should be closed when they're no longer needed, although the system will close any that remain open on process termination.

**dde_request**(`+Handle, +Item, -Value`)  
Request a value from the server. `Item` is an atom that identifies the requested data, and `Value` will be a string (`CF_TEXT` data in DDE parlance) representing that data, if the request is successful.

**dde_execute**(`+Handle, +Command`)  
Request the DDE server to execute the given command string. Succeeds if the command could be executed and fails with an error message otherwise.

**dde_poke**(`+Handle, +Item, +Command`)  
Issue a `POKE` command to the server on the specified `Item`. `command` is passed as data of type `CF_TEXT`.

### 4.44.2 DDE server mode

The `library(dde)` defines primitives to realise simple DDE server applications in SWI-Prolog. These features are provided as of version 2.0.6 and should be regarded as prototypes. The C part of the DDE server can handle some more primitives, so if you need features not provided by this interface, please study `library(dde)`.

**dde_register_service**(`+Template, +Goal`)  
Register a server to handle DDE request or DDE `execute` requests from other applications. To register a service for a DDE request, `Template` is of the form:

> \+**Service(+Topic, +Item, +Value)**

`Service` is the name of the DDE service provided (like **progman** in the client example above). `Topic` is either an atom, indicating `Goal` only handles requests on this topic, or a variable that also appears in `Goal`. `Item` and `Value` are variables that also appear in `Goal`. `Item` represents the request data as a Prolog atom.^(169Up to version 3.4.5 this was a list of character codes. As recent versions have atom garbage collection there is no need for this anymore.)

The example below registers the Prolog [current_prolog_flag/2](flags.html#current_prolog_flag/2) predicate to be accessible from other applications. The request may be given from the same Prolog as well as from another application.

``` code
?- dde_register_service(prolog(current_prolog_flag, F, V),
                        current_prolog_flag(F, V)).

?- open_dde_conversation(prolog, current_prolog_flag, Handle),
   dde_request(Handle, home, Home),
   close_dde_conversation(Handle).

Home = '/usr/local/lib/pl-2.0.6/'
```

Handling DDE `execute` requests is very similar. In this case the template is of the form:

> \+**Service(+Topic, +Item)**

Passing a `Value` argument is not needed as `execute` requests either succeed or fail. If `Goal` fails, a‘not processed’is passed back to the caller of the DDE request.

**dde_unregister_service**(`+Service`)  
Stop responding to `Service`. If Prolog is halted, it will automatically call this on all open services.

**dde_current_service**(`-Service, -Topic`)  
Find currently registered services and the topics served on them.

**dde_current_connection**(`-Service, -Topic`)  
Find currently open conversations.
