
## B.9 Hooks for integrating libraries

Some libraries realise an entirely new programming paradigm on top of Prolog. An example is XPCE which adds an object system to Prolog as well as an extensive set of graphical primitives. SWI-Prolog provides several hooks to improve the integration of such libraries. See also [section A.24](listing.html#sec:A.24) for editing hooks and [section 4.11](printmsg.html#sec:4.11) for hooking into the message system.

**prolog_list_goal**(`:Goal`)  
Hook, normally not defined. This hook is called by the’L’command of the tracer in the module `user` to list the currently called predicate. This hook may be defined to list only relevant clauses of the indicated `Goal` and/or show the actual source code in an editor. See also [portray/1](termrw.html#portray/1) and [multifile/1](dynamic.html#multifile/1).

**prolog:debug_control_hook**(`:Action`)  
Hook for the debugger control predicates that allows the creator of more high-level programming languages to use the common front-end predicates to control the debugger. For example, XPCE uses these hooks to allow for spying methods rather than predicates. `Action` is one of:

**spy**(`Spec`)  
Hook in [spy/1](debugger.html#spy/1). If the hook succeeds [spy/1](debugger.html#spy/1) takes no further action.

**nospy**(`Spec`)  
Hook in [nospy/1](debugger.html#nospy/1). If the hook succeeds [nospy/1](debugger.html#nospy/1) takes no further action. If [spy/1](debugger.html#spy/1) is hooked, it is advised to place a complementary hook for [nospy/1](debugger.html#nospy/1).

**nospyall**  
Hook in [nospyall/0](debugger.html#nospyall/0). Should remove all spy points. This hook is called in a failure-driven loop.

**debugging**(`DebugMode`)  
Hook in [debugging/0](debugger.html#debugging/0). `DebugMode` holds the current value of the [debug](flags.html#flag:debug) flag. The hook can be used in two ways. It can report the status of the additional debug points controlled by the above hooks and fail to let the system report the others, or it succeeds, overruling the entire behaviour of [debugging/0](debugger.html#debugging/0).

**prolog:help_hook**(`+Action`)  
Hook into [help/0](online-help.html#help/0) and [help/1](online-help.html#help/1). If the hook succeeds, the built-in actions are not executed. For example, `?- help(picture).` is caught by the XPCE help hook to give help on the class *picture*. Defined actions are:

**help**  
User entered plain [help/0](online-help.html#help/0) to give default help. The default performs `help(`[`help/1`](online-help.html#help/1)`)`, giving help on help.

**help**(`What`)  
Hook in [help/1](online-help.html#help/1) on the topic `What`.

**apropos**(`What`)  
Hook in [apropos/1](online-help.html#apropos/1) on the topic `What`.
