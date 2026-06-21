SWI-Prolog binding to BSD libedit

Jan Wielemaker  
SWI-Prolog Solutions b.v.  
The Netherlands  
E-mail: [jan@swi-prolog.org](mailto:jan@swi-prolog.org)

Abstract

This package enables editing commands in the Prolog toplevel using the BSD libedit library.

# Table of Contents

[1 library(editline): BSD libedit based command line editing](#sec:1)

## 1 library(editline): BSD libedit based command line editing

This library wraps the BSD libedit command line editor. The binding provides a high level API to enable command line editing on the Prolog user streams and low level predicates to apply the library on other streams and program the library.

\[det\]**el_wrap**  
\[det\]**el_wrap**(`+Options`)  
Enable using editline on the standard user streams if `user_input` is connected to a terminal. This is the high level predicate used for most purposes. The remainder of the library interface deals with low level predicates that allows for applying and programming libedit in non-standard situations.

The library is registered with *ProgName* set to `swipl` (see [el_wrap/4](#el_wrap/4)).

`Options` processed:

**pipes**(`+Boolean`)  
Used by Epilog windows to indicate we are reading from a Windows named pipe in *overlapped* mode. Ignored on other platforms.

**history**(`+Size`)  
`Size` of the history. Default is defined by the Prolog flag `history` or `100` if this flag is not defined.

**alert_signo**(`+Integer`)  
Signal used for making thread_signal/2 work while the thread is in a blocking system call.

\[det\]**el_wrap**(`+ProgName:atom, +In:stream, +Out:stream, +Error:stream`)  
\[det\]**el_wrap**(`+ProgName:atom, +In:stream, +Out:stream, +Error:stream, +Options`)  
Enable editline on the stream-triple `<``In`,`Out`,`Error``>`. From this moment on `In` is a handle to the command line editor. `Options`:

**pipes**(`true`)  
Windows only. Assume the I/O is using pipes rather than a console. This is used for the Epilog terminal.

|  |  |
|----|----|
| `ProgName` | is the name of the invoking program, used when reading the `editrc(5)` file to determine which settings to use. |

\[nondet,multifile\]**el_setup**(`+In:stream`)  
This hooks is called as `forall(el_setup(Input), true)` *after* the input stream has been wrapped, the default Prolog commands have been added and the default user setup file has been sourced using [el_source/2](#el_source/2). It can be used to define and bind additional commands.

\[semidet\]**el_wrapped**(`+In:stream`)  
True if `In` is a stream wrapped by el_wrap/3.

\[det\]**el_unwrap**(`+In:stream`)  
Remove the libedit wrapper for `In` and the related output and error streams.

bug  
The wrapper creates `FILE*` handles that cannot be closed and thus wrapping and unwrapping implies a (modest) memory leak.

\[det\]**el_source**(`+In:stream, +File`)  
Initialise editline by reading the contents of `File`. If `File` is unbound try `$HOME/.editrc`

\[det\]**el_bind**(`+In:stream, +Args`)  
Invoke the libedit `bind` command with the given arguments. The example below lists the current key bindings.

``` code
?- el_bind(user_input, ['-a']).
```

The predicate [el_bind/2](#el_bind/2) is typically used to bind commands defined using [el_addfn/4](#el_addfn/4). Note that the C proxy function has only the last character of the command as context to find the Prolog binding. This implies we cannot both bind e.g., "`^`\[?" *and* "?" to a Prolog function.

See also  
`editrc(5)` for more information.

\[det\]**el_addfn**(`+Input:stream, +Command, +Help, :Goal`)  
Add a new command to the command line editor associated with `Input`. `Command` is the name of the command, `Help` is the help string printed with e.g. `bind -a` (see [el_bind/2](#el_bind/2)) and `Goal` is called of the associated key-binding is activated. `Goal` is called as

``` code
call(:Goal, +Input, +Char, -Continue)
```

where `Input` is the input stream providing access to the editor, Char the activating character and Continue must be instantated with one of the known continuation codes as defined by libedit: `norm`, `newline`, `eof`, `arghack`, `refresh`, `refresh_beep`, `cursor`, `redisplay`, `error` or `fatal`. In addition, the following Continue code is provided.

**electric**(`Move, TimeOut, Continue`)  
Show *electric caret* at `Move` positions to the left of the normal cursor positions for the given `TimeOut`. `Continue` as defined by the `Continue` value.

The registered `Goal` typically used [el_line/2](#el_line/2) to fetch the input line and [el_cursor/2](#el_cursor/2), [el_insertstr/2](#el_insertstr/2) and/or [el_deletestr/2](#el_deletestr/2) to manipulate the input line.

Normally [el_bind/2](#el_bind/2) is used to associate the defined command with a keyboard sequence.

See also  
`el_set(3)` `EL_ADDFN` for details.

\[semidet\]**el_set**(`+Input:stream, +Action`)  
Interface to `el_set()` and `el_wset()`. Currently provided values for `Action` are:

**wordchars**(`+Text`)  
Set the characters considered part of a *word*. This feature depends on `el_wsey()` `EL_WORDCHARS`, which is only provided in some recent versions of `libedit`.

This predicate fails silently of `Action` is not implemented. Illegal input raises in an exception.

\[det\]**el_line**(`+Input:stream, -Line`)  
Fetch the currently buffered input line. `Line` is a term `line(Before, After)`, where `Before` is a string holding the text before the cursor and `After` is a string holding the text after the cursor.

\[det\]**el_cursor**(`+Input:stream, +Move:integer`)  
`Move` the cursor `Move` character forwards (positive) or backwards (negative).

\[det\]**el_insertstr**(`+Input:stream, +Text`)  
Insert `Text` at the cursor.

\[det\]**el_deletestr**(`+Input:stream, +Count`)  
Delete `Count` characters before the cursor.

\[det\]**el_history**(`+In:stream, ?Action`)  
Perform a generic action on the history. This provides an incomplete interface to `history()` from libedit. Supported actions are:

**clear**  
Clear the history.

**setsize**(`+Integer`)  
Set size of history to size elements.

**getsize**(`-Integer`)  
Unify `Integer` with the maximum size of the history. Note that this is *not* the same as `el_history()` using `H_GETSIZE`, which returns the number of currently saved events. The number of saved events may be computed from `first` and `last` or using [el_history_events/2](#el_history_events/2).

**setunique**(`+Boolean`)  
Set flag that adjacent identical event strings should not be entered into the history.

**first**(`-Num, -String`)  
**last**(`-Num, -String`)  
**curr**(`-Num, -String`)  
**prev**(`-Num, -String`)  
**next**(`-Num, -String`)  
Retrieve an event. `Num` is the event number and `String` is the event string. Note that `first` is the most recent event and `last` the oldest.

**set**(`Num`)  
Set the notion of *current* to `Num`.

**prev_str**(`+Search, -Num, -String`)  
**next_str**(`+Search, -Num, -String`)  
Retrieve the previous or next event whose `String` starts with `Search`.

**event**(`+Num, -String`)  
True when `String` represents event `Num`. This is an extension to the `history()` API, retrieving a numbered event without changing the current notion.

\[det\]**el_history_events**(`+In:stream, -Events:list(pair)`)  
Unify `Events` with a list of pairs of the form `Num-String`, where `Num` is the event number and `String` is the associated string without terminating newline.

\[det\]**el_add_history**(`+In:stream, +Line:text`)  
Add a line to the command line history.

\[det\]**el_read_history**(`+In:stream, +File:file`)  
Read the history saved using [el_write_history/2](#el_write_history/2).

|        |                                                   |
|--------|---------------------------------------------------|
| `File` | is a file specification for absolute_file_name/3. |

\[det\]**el_write_history**(`+In:stream, +File:file`)  
Save editline history to `File`. The history may be reloaded using [el_read_history/2](#el_read_history/2).

|        |                                                   |
|--------|---------------------------------------------------|
| `File` | is a file specification for absolute_file_name/3. |

**el_version**(`-Version`)  
True when `Version` is `LIBEDIT_MAJOR*10000 + LIBEDIT_MINOR*100`. The version is generated from the include file `histedit.h`, which implies that the actual version of the shared library may be different.

\[semidet,multifile\]prolog:**history**(`+Input, ?Action`)  
Provide the plugable interface into the system command line management.

# Index

?  
[el_add_history/2](#el_add_history/2)  
[el_addfn/4](#el_addfn/4)  
[el_bind/2](#el_bind/2)  
[el_cursor/2](#el_cursor/2)  
[el_deletestr/2](#el_deletestr/2)  
[el_history/2](#el_history/2)  
[el_history_events/2](#el_history_events/2)  
[el_insertstr/2](#el_insertstr/2)  
[el_line/2](#el_line/2)  
[el_read_history/2](#el_read_history/2)  
[el_set/2](#el_set/2)  
[el_setup/1](#el_setup/1)  
[el_source/2](#el_source/2)  
[el_unwrap/1](#el_unwrap/1)  
[el_version/1](#el_version/1)  
[el_wrap/0](#el_wrap/0)  
[el_wrap/1](#el_wrap/1)  
[el_wrap/4](#el_wrap/4)  
[el_wrap/5](#el_wrap/5)  
[el_wrapped/1](#el_wrapped/1)  
[el_write_history/2](#el_write_history/2)  
[prolog:history/2](#prolog:history/2)  
