
## A.14 library(debug): Print debug messages and test assertions

This library is a replacement for [format/3](format.html#format/3) for printing debug messages. Messages are assigned a *topic*. By dynamically enabling or disabling topics the user can select desired messages. Calls to [debug/3](debug.html#debug/3) and [assertion/1](debug.html#assertion/1) are removed when the code is compiled for optimization unless the Prolog flag **optimise_debug** is set to `true`.

Using the predicate [assertion/1](debug.html#assertion/1) you can make assumptions about your program explicit, trapping the debugger if the condition does not hold.

Output and actions by these predicates can be configured using *hooks* to fit your environment. With XPCE, you can use the call below to start a graphical monitoring tool.

``` code
?- prolog_ide(debug_monitor).
```

\[semidet\]**debugging**(`+Topic`)  
\[nondet\]**debugging**(`-Topic`)  
\[nondet\]**debugging**(`?Topic, ?Bool`)  
Examine debug topics. The form `debugging(+Topic)` may be used to perform more complex debugging tasks. A typical usage skeleton is:

``` code
      (   debugging(mytopic)
      ->  <perform debugging actions>
      ;   true
      ),
      ...
```

The other two calls are intended to examine existing and enabled debugging tokens and are typically not used in user programs.

\[det\]**debug**(`+Topic`)  
\[det\]**nodebug**(`+Topic`)  
Add/remove a topic from being printed. `nodebug(_)` removes all topics. Gives a warning if the topic is not defined unless it is used from a directive. The latter allows placing debug topics at the start of a (load-)file without warnings.

For [debug/1](debug.html#debug/1), `Topic` can be a term `Topic > Out`, where `Out` is either a stream or stream-alias or a filename (an atom). This redirects debug information on this topic to the given output. On Linux systems redirection can be used to make the message appear, even if the `user_error` stream is redefined using

``` code
?- debug(Topic > '/proc/self/fd/2').
```

A platform independent way to get debug messages in the current console (for example, a `swipl-win` window, or login using `ssh` to Prolog running an SSH server from the `libssh` pack) is to use:

``` code
?- stream_property(S, alias(user_error)),
   debug(Topic > S).
```

Do not forget to disable the debugging using [nodebug/1](debug.html#nodebug/1) before quitting the console if Prolog must remain running.

\[det\]**list_debug_topics**  
\[det\]**list_debug_topics**(`+Options`)  
List currently known topics for [debug/3](debug.html#debug/3) and their setting. `Options` is either an atom or string, which is a shorthand for `[search(String)]` or a normal option list. Defined options are:

**search**(`String`)  
Only show topics that match `String`. Match is case insensitive on the printed representation of the term.

**active**(`+Boolean`)  
Only print topics that are active (`true`) or inactive (`false`).

**output**(`+To`)  
Only print topics whose target location matches `To`. This option implicitly restricts the output to active topics.

\[det\]**debug_message_context**(`+What`)  
Specify additional context for debug messages.

deprecated  
New code should use the Prolog flag [message_context](flags.html#flag:message_context). This predicates adds or deletes topics from this list.

\[det\]**debug**(`+Topic, +Format, :Args`)  
`Format` a message if debug topic is enabled. Similar to [format/3](format.html#format/3) to `user_error`, but only prints if `Topic` is activated through [debug/1](debug.html#debug/1). `Args` is a meta-argument to deal with goal for the @-command. Output is first handed to the hook [prolog:debug_print_hook/3](debug.html#prolog:debug_print_hook/3). If this fails, `Format`+`Args` is translated to text using the message-translation (see [print_message/2](printmsg.html#print_message/2)) for the term `debug(Format, Args)` and then printed to every matching destination (controlled by [debug/1](debug.html#debug/1)) using [print_message_lines/3](printmsg.html#print_message_lines/3).

The message is preceded by’% ’and terminated with a newline.

See also  
[format/3](format.html#format/3).

\[semidet,multifile\]prolog:**debug_print_hook**(`+Topic, +Format, +Args`)  
Hook called by [debug/3](debug.html#debug/3). This hook is used by the graphical frontend that can be activated using [prolog_ide/1](idepreds.html#prolog_ide/1):

``` code
?- prolog_ide(debug_monitor).
```

\[det\]**assertion**(`:Goal`)  
Acts similar to C `assert()` macro. It has no effect if `Goal` succeeds. If `Goal` fails or throws an exception, the following steps are taken:

- call [prolog:assertion_failed/2](debug.html#prolog:assertion_failed/2). If [prolog:assertion_failed/2](debug.html#prolog:assertion_failed/2) fails, then:
  - If this is an interactive toplevel thread, print a message, the stack-trace, and finally trap the debugger.
  - Otherwise, throw `error(assertion_error(Reason, G),_)` where Reason is one of `fail` or the exception raised.

\[semidet,multifile\]prolog:**assertion_failed**(`+Reason, +Goal`)  
This hook is called if the `Goal` of [assertion/1](debug.html#assertion/1) fails. `Reason` is unified with either `fail` if `Goal` simply failed or an exception call otherwise. If this hook fails, the default behaviour is activated. If the hooks throws an exception it will be propagated into the caller of [assertion/1](debug.html#assertion/1).
