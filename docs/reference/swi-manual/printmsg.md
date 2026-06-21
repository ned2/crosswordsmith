
## 4.11 Printing messages

The predicate [print_message/2](printmsg.html#print_message/2) is used to print a message term in a human-readable format. The other predicates from this section allow the user to refine and extend the message system. A common usage of [print_message/2](printmsg.html#print_message/2) is to print error messages from exceptions. The code below prints errors encountered during the execution of `Goal`, without further propagating the exception and without starting the debugger.

``` code
        ...,
        catch(Goal, E,
              ( print_message(error, E),
                fail
              )),
        ...
```

Another common use is to define [message_hook/3](printmsg.html#message_hook/3) for printing messages that are normally *silent*, suppressing messages, redirecting messages or make something happen in addition to printing the message.

**print_message**(`+Kind, +Term`)  
The predicate [print_message/2](printmsg.html#print_message/2) is used by the system and libraries to print messages. `Kind` describes the nature of the message, while `Term` is a Prolog term that describes the content. Printing messages through this indirection instead of using [format/3](format.html#format/3) to the stream `user_error` allows displaying the message appropriate to the application (terminal, logfile, graphics), acting on messages based on their content instead of a string (see [message_hook/3](printmsg.html#message_hook/3)) and creating language specific versions of the messages. See also [section 4.11.1](printmsg.html#sec:4.11.1). The following message kinds are known:

**banner**  
The system banner message. Banner messages can be suppressed by setting the Prolog flag [verbose](flags.html#flag:verbose) to `silent`.

**debug**(`Topic`)  
Message from library(debug). See [debug/3](debug.html#debug/3).

**error**  
The message indicates an erroneous situation. This kind is used to print uncaught exceptions of type `error(Formal, Context)`. See section introduction ([section 4.11](printmsg.html#sec:4.11)). An error message causes the process to halt with status 1 if the Prolog flag [on_error](flags.html#flag:on_error) is set to `halt` and the message is not intercepted by [message_hook/3](printmsg.html#message_hook/3). Not intercepted error messages increment the `errors` key for [statistics/2](builtin-statistics.html#statistics/2).

**help**  
User requested help message, for example after entering‘h’or‘?’to a prompt.

**information**  
Information that is requested by the user. An example is [statistics/0](statistics.html#statistics/0).

**informational**  
Typically messages of events and progress that are considered useful to a developer. Such messages can be suppressed by setting the Prolog flag [verbose](flags.html#flag:verbose) to `silent`.

**silent**  
Message that is normally not printed. Applications may define [message_hook/3](printmsg.html#message_hook/3) to act upon such messages.

**trace**  
Messages from the (command line) tracer.

**warning**  
The message indicates something dubious that is not considered fatal. For example, discontiguous predicates (see [discontiguous/1](dynamic.html#discontiguous/1)). A warning message causes the process to halt with status 1 if the Prolog flag [on_warning](flags.html#flag:on_warning) is set to `halt` and the message is not intercepted by [message_hook/3](printmsg.html#message_hook/3). Not intercepted warning messages increment the `warnings` key for [statistics/2](builtin-statistics.html#statistics/2).

The predicate [print_message/2](printmsg.html#print_message/2) first translates the `Term` into a list of‘message lines’(see [print_message_lines/3](printmsg.html#print_message_lines/3) for details). Next, it calls the hook [message_hook/3](printmsg.html#message_hook/3) to allow the user to intercept the message. If [message_hook/3](printmsg.html#message_hook/3) fails it prints the message unless `Kind` is `silent`.

The [print_message/2](printmsg.html#print_message/2) predicate and its rules are in the file `<``plhome``>/boot/messages.pl`, which may be inspected for more information on the error messages and related error terms. If you need to write messages from your own predicates, it is recommended to reuse the existing message terms if applicable. If no existing message term is applicable, invent a fairly unique term that represents the event and define a rule for the multifile predicate prolog:message//1. See [section 4.11.1](printmsg.html#sec:4.11.1) for a deeper discussion and examples.

See also [message_to_string/2](printmsg.html#message_to_string/2).

**print_message_lines**(`+Stream, +Prefix, +Lines`)  
Print a message (see [print_message/2](printmsg.html#print_message/2)) that has been translated to a list of message elements. The elements of this list are:

**\<`Format`\>-\<`Args`\>**  
Where `Format` is an atom and `Args` is a list of format arguments. Handed to [format/3](format.html#format/3).

**flush**  
If this appears as the last element, `Stream` is flushed (see [flush_output/1](chario.html#flush_output/1)) and no final newline is generated. This is combined with a subsequent message that starts with `at_same_line` to complete the line.

**at_same_line**  
If this appears as first element, no prefix is printed for the first line and the line position is not forced to 0 (see [format/1](format.html#format/1), `~N`).

**ansi**(`+Attributes, +Format, +Args`)  
This message may be intercepted by means of the hook [prolog:message_line_element/2](printmsg.html#prolog:message_line_element/2). The library `library(ansi_term)` implements this hook to achieve coloured output. If it is not intercepted it invokes `format(Stream, Format, Args)`.

**url**(`Location`)  
Print a source location. `Location` is one of the terms `File:Line:Column`, `File:Line` or `File`. When using library `library(ansi_term)`, this is translated into a hyperlink for modern terminals.

**url**(`URL, Label`)  
Print `Label`. When using library `library(ansi_term)`, this is translated into a hyperlink for modern terminals.

**nl**  
A new line is started. If the message is not complete, `Prefix` is printed before the remainder of the message.

**begin**(`Kind, Var`)  
**end**(`Var`)  
The entire message is headed by `begin(Kind, Var)` and ended by `end(Var)`. This feature is used by, e.g., library `library(ansi_term)` to colour entire messages.

**\<`Format`\>**  
Handed to [format/3](format.html#format/3) as `format(Stream, Format,[])`. Deprecated because it is ambiguous if `Format` collides with one of the atomic commands.

See also [print_message/2](printmsg.html#print_message/2) and [message_hook/3](printmsg.html#message_hook/3).

**message_hook**(`+Term, +Kind, +Lines`)  
Hook predicate that may be defined in the module `user` to intercept messages from [print_message/2](printmsg.html#print_message/2). `Term` and `Kind` are the same as passed to [print_message/2](printmsg.html#print_message/2). `Lines` is a list of format statements as described with [print_message_lines/3](printmsg.html#print_message_lines/3). See also [message_to_string/2](printmsg.html#message_to_string/2).

If calling this hook succeeds, the message is considered printed and no further action is taken. Version 9.3.34 introduced [prolog:message_action/2](printmsg.html#prolog:message_action/2) to allow modules to act on certain actions. Before [prolog:message_action/2](printmsg.html#prolog:message_action/2) was introduced a program could associate an action with a message by adding a clause for [message_hook/3](printmsg.html#message_hook/3) that realises the desired action and then fails. The problem with this approach for actions is that the order of clauses for a multi-file predicate is poorly defined. If the application defines a hook that succeeds, this hook may cause intended actions not to take place.

This predicate must be defined dynamic and multifile to allow other modules defining clauses for it too.

**thread_message_hook**(`+Term, +Kind, +Lines`)  
As [message_hook/3](printmsg.html#message_hook/3), but this predicate is local to the calling thread (see [thread_local/1](threadcom.html#thread_local/1)). This hook is called *before* [message_hook/3](printmsg.html#message_hook/3). The‘pre-hook’is indented to catch messages they may be produced by calling some goal without affecting other threads.

\[nondet\]**prolog:message_action**(`+Term, +Kind`)  
This hook is called before [message_hook/3](printmsg.html#message_hook/3) as below. This hook appeared in 9.3.34. It allows to associate side-effects with messages in a reliable way by decoupling printing (or reporting some other way) from associated side effects. See also [broadcast/1](broadcast.html#broadcast/1) from library `library(broadcast)`.

``` code
    forall(prolog:message_action(Term, Kind), true).
```

**message_property**(`+Kind, ?Property`)  
This hook can be used to define additional message kinds and the way they are displayed. The following properties are defined:

**color**(`-Attributes`)  
Print message using ANSI terminal attributes. See [ansi_format/3](ansiterm.html#ansi_format/3) for details. Here is an example, printing help messages in blue:

``` code
:- multifile user:message_property/2.

user:message_property(help, color([fg(blue)])).
```

**prefix**(`-Prefix`)  
Prefix printed before each line. This argument is handed to [format/3](format.html#format/3). The default is `'~N'`. For example, messages of kind `warning` use `'~NWarning: '`.

**tag**(`-Tag`)  
Defines the text part for the `prefix` property for error and warning messages.

**location_prefix**(`+Location, -FirstPrefix, -ContinuePrefix`)  
Used for printing messages that are related to a source location. Currently, `Location` is a term `File`:`Line`. `FirstPrefix` is the prefix for the first line and `-ContinuePrefix` is the prefix for continuation lines. For example, the default for errors is

``` code
location_prefix(File:Line,
                '~NERROR: ~w:~d:'-[File,Line], '~N\t')).
```

**stream**(`-Stream`)  
Stream to which to print the message. Default is `user_error`.

**wait**(`-Seconds`)  
Amount of time to wait after printing the message. Default is not to wait.

**prolog:message_line_element**(`+Stream, +Term`)  
This hook is called to print the individual elements of a message from [print_message_lines/3](printmsg.html#print_message_lines/3). This hook is used by e.g., library `library(ansi_term)` to colour messages on ANSI-capable terminals.

**prolog:message_prefix_hook**(`+ContextTerm, -Prefix`)  
This hook is called to add context to the message prefix. `ContextTerm` is a member of the list provided by the [message_context](flags.html#flag:message_context). `Prefix` must be unified with an atomic value that is added to the message prefix.

**message_to_string**(`+Term, -String`)  
Translates a message term into a string object (see [section 5.2](string.html#sec:5.2)).

**version**  
Write the SWI-Prolog banner message as well as additional messages registered using [version/1](printmsg.html#version/1). This is the default *initialization goal* which can be modified using **-g**.

**version**(`+Message`)  
Register additional messages to be printed by [version/0](printmsg.html#version/0). Each registered message is handed to the message translation DCG and can thus be defined using the hook prolog:message//1. If not defined, it is simply printed.

### 4.11.1 Printing from libraries

Libraries should *not* use [format/3](format.html#format/3) or other output predicates directly. Libraries that print informational output directly to the console are hard to use from code that depend on your textual output, such as a CGI script. The predicates in [section 4.11](printmsg.html#sec:4.11) define the API for dealing with messages. The idea behind this is that a library that wants to provide information about its status, progress, events or problems calls [print_message/2](printmsg.html#print_message/2). The first argument is the *level*. The supported levels are described with [print_message/2](printmsg.html#print_message/2). Libraries typically use `informational` and `warning`, while libraries should use exceptions for errors (see [throw/1](exception.html#throw/1), [type_error/2](error.html#type_error/2), etc.).

The second argument is an arbitrary Prolog term that carries the information of the message, but *not* the precise text. The text is defined by the grammar rule prolog:message//1. This distinction is made to allow for translations and to allow hooks processing the information in a different way (e.g., to translate progress messages into a progress bar).

For example, suppose we have a library that must download data from the Internet (e.g., based on http_open/3). The library wants to print the progress after each downloaded file. The code below is a good skeleton:

``` code
download_urls(List) :-
        length(List, Total),
        forall(nth1(I, List, URL),
               (   download_url(URL),
                   print_message(informational,
                                 download_url(URL, I, Total)))).
```

The programmer can now specify the default textual output using the rule below. Note that this rule may be in the same file or anywhere else. Notably, the application may come with several rule sets for different languages. This, and the user-hook example below are the reason to represent the message as a compound term rather than a string. This is similar to using message numbers in non-symbolic languages. The documentation of [print_message_lines/3](printmsg.html#print_message_lines/3) describes the elements that may appear in the output list.

``` code
:- multifile
        prolog:message//1.

prolog:message(download_url(URL, I, Total)) -->
        { Perc is round(I*100/Total) },
        [ 'Downloaded ~w; ~D from ~D (~d%)'-[URL, I, Total, Perc] ].
```

A *user* of the library may define rules for [message_hook/3](printmsg.html#message_hook/3). The rule below acts on the message content. Other applications can act on the message level and, for example, popup a message box for warnings and errors.

``` code
:- multifile user:message_hook/3.

message_hook(download_url(URL, I, Total), _Kind, _Lines) :-
        <send this information to a GUI component>
```

In addition, using the command line option **-q**, the user can disable all *informational* messages.
