
## 4.17 Input and output

SWI-Prolog provides two different packages for input and output. The native I/O system is based on the ISO standard predicates [open/3](IO.html#open/3), [close/1](IO.html#close/1) and friends.^(97Actually based on Quintus Prolog, providing this interface before the ISO standard existed.) Being more widely portable and equipped with a clearer and more robust specification, new code is encouraged to use these predicates for manipulation of I/O streams.

[Section 4.17.3](IO.html#sec:4.17.3) describes [tell/1](IO.html#tell/1), [see/1](IO.html#see/1) and friends, providing I/O in the spirit of the traditional Edinburgh standard. These predicates are layered on top of the ISO predicates. Both packages are fully integrated; the user may switch freely between them.

### 4.17.1 Predefined stream aliases

Each thread has five stream aliases: `user_input`, `user_output`, `user_error`, `current_input`, and `current_output`. Newly created threads inherit these stream aliases from their parent. The `user_input`, `user_output` and `user_error` aliases of the `main` thread are initially bound to the standard operating system I/O streams (*stdin*, *stdout* and *stderr*, normally bound to the POSIX file handles 0, 1 and 2). These aliases may be re-bound, for example if standard I/O refers to a window such as in the **swipl-win.exe** GUI executable for Windows. They can be re-bound by the user using [set_prolog_IO/3](IO.html#set_prolog_IO/3) and [set_stream/2](IO.html#set_stream/2) by setting the alias of a stream (e.g, `set_stream(S, alias(user_output))`). An example of rebinding can be found in library `library(prolog_server)`, providing a **telnet** service. The aliases `current_input` and `current_output` define the source and destination for predicates that do not take a stream argument (e.g., [read/1](termrw.html#read/1), [write/1](termrw.html#write/1), [get_code/1](chario.html#get_code/1), ... ). Initially, these are bound to the same stream as `user_input` and `user_error`. They are re-bound by [see/1](IO.html#see/1), [tell/1](IO.html#tell/1), [set_input/1](IO.html#set_input/1) and [set_output/1](IO.html#set_output/1). The `current_output` stream is also temporary re-bound by [with_output_to/2](IO.html#with_output_to/2) or [format/3](format.html#format/3) using e.g., `format(atom(A), ...`. Note that code which explicitly writes to the streams `user_output` and `user_error` will not be redirected by [with_output_to/2](IO.html#with_output_to/2).

**Compatibility**

Note that the ISO standard only defines the `user_*` streams. The‘current’streams can be accessed using [current_input/1](IO.html#current_input/1) and [current_output/1](IO.html#current_output/1). For example, an ISO compatible implementation of [write/1](termrw.html#write/1) is

``` code
write(Term) :- current_output(Out), write_term(Out, Term).
```

while SWI-Prolog additionally allows for

``` code
write(Term) :- write(current_output, Term).
```

### 4.17.2 ISO Input and Output Streams

The predicates described in this section provide ISO compliant I/O, where streams are explicitly created using the predicate [open/3](IO.html#open/3). The resulting stream identifier is then passed as a parameter to the reading and writing predicates to specify the source or destination of the data.

This schema is not vulnerable to filename and stream ambiguities as well as changes to the working directory. On the other hand, using the notion of current-I/O simplifies reusability of code without the need to pass arguments around. E.g., see [with_output_to/2](IO.html#with_output_to/2).

SWI-Prolog streams are, compatible with the ISO standard, either input or output streams. To accommodate portability to other systems, a pair of streams can be packed into a *stream-pair*. See [stream_pair/3](IO.html#stream_pair/3) for details.

SWI-Prolog stream handles are unique symbols that have no syntactical representation. They are written as `<stream>(hex-number)`, which is not valid input for [read/1](termrw.html#read/1). They are realised using a *blob* of type `stream` (see [blob/2](typetest.html#blob/2) and [section 12.4.10](foreigninclude.html#sec:12.4.10)).

\[ISO\]**open**(`+SrcDest, +Mode, --Stream, +Options`)  
True when `SrcDest` can be opened in `Mode` and `Stream` is an I/O stream to/from the object. `SrcDest` is normally the name of a file, represented as an atom or string. `Mode` is one of `read`, `write`, `append` or `update`. Mode `append` opens the file for writing, positioning the file pointer at the end. Mode `update` opens the file for writing, positioning the file pointer at the beginning of the file without truncating the file. `Stream` is either a variable, in which case it is bound to an integer identifying the stream, or an atom, in which case this atom will be the stream identifier.^(98New code should use the `alias(Alias)` option for compatibility with the ISO standard.)

SWI-Prolog also allows `SrcDest` to be a term `pipe(Command)`. In this form, `Command` is started as a child process and if `Mode` is `write`, output written to `Stream` is sent to the standard input of `Command`. Vice versa, if `Mode` is `read`, data written by `Command` to the standard output can be read from `Stream`. On Unix systems, `Command` is handed to **popen()** which hands it to the Unix shell. On Windows, `Command` is executed directly and therefore shell syntax such as redirecting (using e.g., `>` `file`) does not work. Use of the `pipe(Command)` feature is deprecated. The predicate process_create/3 from `library(process)` provides a richer and more portable alternative for interacting with processes including handling all three standard streams.

If `SrcDest` is an *IRI*, i.e., starts with \<`scheme`\>`://`, where \<`scheme`\> is a non-empty sequence of lowercase ASCII letters [open/3](IO.html#open/3),4 calls hooks registered by [register_iri_scheme/3](IO.html#register_iri_scheme/3). Currently the only predefined IRI scheme is `res`, providing access to the *resource database*. See [section 14.4](program-resources.html#sec:14.4).

The following `Options` are recognised by [open/4](IO.html#open/4):

**alias**(`Atom`)  
Gives the stream a name and unifies Stream with Atom. Below is an example. Be careful with this option as stream names are global. See also [set_stream/2](IO.html#set_stream/2).

``` code
?- open(data, read, Fd, [alias(input)]).

        ...,
        read(input, Term),
        ...
```

**bom**(`Bool`)  
Check for a BOM (*Byte Order Marker*) or write one. If omitted, the default is `true` for mode `read` and `false` for mode `write`. See also [stream_property/2](IO.html#stream_property/2) and especially [section 2.18.1.1](widechars.html#sec:2.18.1.1) for a discussion of this feature.

**buffer**(`Buffering`)  
Defines output buffering. The atom `full` (default) defines full buffering, `line` buffering by line, and `false` implies the stream is fully unbuffered. Smaller buffering is useful if another process or the user is waiting for the output as it is being produced. See also [flush_output/\[0,1\]](chario.html#flush_output/0). This option is not an ISO option.

**close_on_abort**(`Bool`)  
If `true` (default), the stream is closed on an abort (see [abort/0](toplevel.html#abort/0)). If `false`, the stream is not closed. If it is an output stream, however, it will be flushed. Useful for logfiles and if the stream is associated to a process (using the `pipe/1` construct).

**create**(`+List`)  
Specifies how a new file is created when opening in `write`, `append` or `update` mode. Currently, `List` is a list of atoms that describe the permissions of the created file.^(99Added after feedback from Joachim Shimpf and Per Mildner.) Defined values are below. Not recognised values are silently ignored, allowing for adding platform specific extensions to this set.

**read**  
Allow read access to the file.

**write**  
Allow write access to the file.

**execute**  
Allow execution access to the file.

**default**  
Allow read and write access to the file.

**all**  
Allow any access provided by the OS.

Note that if `List` is empty, the created file has no associated access permissions. The create options map to the POSIX `mode` option of **open()**, where `read` map to 0444, `write` to 0222 and `execute` to 0111. On POSIX systems, the final permission is defined as (mode & `~`umask).

**encoding**(`Encoding`)  
Define the encoding used for reading and writing text to this stream. The default encoding for type `text` is derived from the Prolog flag [encoding](flags.html#flag:encoding). For `binary` streams the default encoding is `octet`. For details on encoding issues, see [section 2.18.1](widechars.html#sec:2.18.1).

**eof_action**(`Action`)  
Defines what happens if the end of the input stream is reached. The default value for Action is `eof_code`, which makes [get0/1](chario.html#get0/1) and friends return -1, and [read/1](termrw.html#read/1) and friends return the atom `end_of_file`. Repetitive reading keeps yielding the same result. Action `error` is like `eof_code`, but repetitive reading will raise an error. With action `reset`, Prolog will examine the file again and return more data if the file has grown.

**locale**(`+Locale`)  
Set the locale that is used by notably [format/2](format.html#format/2) for output on this stream. See [section 4.23](locale.html#sec:4.23).

**lock**(`LockingMode`)  
Try to obtain a lock on the open file. Default is `none`, which does not lock the file. The value `read` or `shared` means other processes may read the file, but not write it. The value `write` or `exclusive` means no other process may read or write the file.

Locks are acquired through the POSIX function **fcntl()** using the command `F_SETLKW`, which makes a blocked call wait for the lock to be released. Please note that **fcntl()** locks are *advisory* and therefore only other applications using the same advisory locks honour your lock. As there are many issues around locking in Unix, especially related to NFS (network file system), please study the **fcntl()** manual page before trusting your locks!

The `lock` option is a SWI-Prolog extension.

**newline**(`Mode`)  
Set end-of-line processing for the stream. `Mode` is one of `posix`, `dos` or `detect`. This option is ignored for binary streams. Using `detect` on an *output stream* raises an exception. See also [set_stream/2](IO.html#set_stream/2).

**reposition**(`+Bool`)  
If `false` (default `true`), drop the position tracking logic from the stream. This disables the use of stream_position/3 on this stream.

**type**(`Type`)  
Using type `text` (default), Prolog will write a text file in an operating system compatible way. Using type `binary` the bytes will be read or written without any translation. See also the option `encoding`.

**wait**(`Bool`)  
This option can be combined with the `lock` option. If `false` (default `true`), the open call returns immediately with an exception if the file is locked. The exception has the format `permission_error(lock, source_sink, SrcDest)`.

**unicode_atoms**(`Mode`)  
Set the per-stream atom-content policy used by the term reader on this stream. See [read_term/2](termrw.html#read_term/2),3 for the meaning of each value (`accept`, `nfc`, `error`, `reject`).

\[ISO\]**open**(`+SrcDest, +Mode, --Stream`)  
Equivalent to [open/4](IO.html#open/4) with an empty option list.

**open_null_stream**(`--Stream`)  
Open an output stream that produces no output. All counting functions are enabled on such a stream. It can be used to discard output (like Unix `/dev/null`) or exploit the counting properties. The initial encoding of `Stream` is `utf8`, enabling arbitrary Unicode output. The encoding can be changed to determine byte counts of the output in a particular encoding or validate if output is possible in a particular encoding. For example, the code below determines the number of characters emitted when writing `Term`.

``` code
write_length(Term, Len) :-
        open_null_stream(Out),
        write(Out, Term),
        character_count(Out, Len0),
        close(Out),
        Len = Len0.
```

\[ISO\]**close**(`+Stream`)  
Close the specified stream. If `Stream` is not open, an existence error is raised. See [stream_pair/3](IO.html#stream_pair/3) for the implications of closing a *stream pair*.

If the closed stream is the current input, output or error stream, the stream alias is bound to the initial standard I/O streams of the process. Calling [close/1](IO.html#close/1) on the initial standard I/O streams of the process is a no-op for an input stream and flushes an output stream without closing it.^(100This behaviour was defined with purely interactive usage of Prolog in mind. Applications should not count on this behaviour. Future versions may allow for closing the initial standard I/O streams.)

\[ISO\]**close**(`+Stream, +Options`)  
Provides `close(Stream, [force(true)])` as the only option. Called this way, any resource errors (such as write errors while flushing the output buffer) are ignored.

\[ISO\]**stream_property**(`?Stream, ?StreamProperty`)  
True when `StreamProperty` is a property of `Stream`. If enumeration of streams or properties is demanded because either `Stream` or `StreamProperty` are unbound, the implementation enumerates all candidate streams and properties while locking the stream database. Properties are fetched without locking the stream and may be outdated before this predicate returns due to asynchronous activity.

**alias**(`Atom`)  
If `Atom` is bound, test if the stream has the specified alias. Otherwise unify `Atom` with the first alias of the stream.^(bugBacktracking does not give other aliases.)

**buffer**(`Buffering`)  
SWI-Prolog extension to query the buffering mode of this stream. `Buffering` is one of `full`, `line` or `false`. See also [open/4](IO.html#open/4).

**buffer_size**(`Integer`)  
SWI-Prolog extension to query the size of the I/O buffer associated to a stream in bytes. Fails if the stream is not buffered.

**bom**(`Bool`)  
If present and `true`, a BOM (*Byte Order Mark*) was detected while opening the file for reading, or a BOM was written while opening the stream. See [section 2.18.1.1](widechars.html#sec:2.18.1.1) for details.

**close_on_abort**(`Bool`)  
Determine whether or not [abort/0](toplevel.html#abort/0) closes the stream. By default streams are closed.

**close_on_exec**(`Bool`)  
Determine whether or not the stream is closed when executing a new process (**exec()** in Unix, **CreateProcess()** in Windows). Default is to close streams. This maps to **fcntl()** `F_SETFD` using the flag `FD_CLOEXEC` on Unix and (negated) `HANDLE_FLAG_INHERIT` on Windows.

**encoding**(`Encoding`)  
Query the encoding used for text. See [section 2.18.1](widechars.html#sec:2.18.1) for an overview of wide character and encoding issues in SWI-Prolog.

**end_of_stream**(`E`)  
If `Stream` is an input stream, unify `E` with one of the atoms `not`, `at` or `past`. See also [at_end_of_stream/\[0,1\]](chario.html#at_end_of_stream/0).

**eof_action**(`A`)  
Unify `A` with one of `eof_code`, `reset` or `error`. See [open/4](IO.html#open/4) for details.

**error**(`Bool`)  
When `true`, the stream is in an error state. Applies to both input and output streams.

**file_name**(`Atom`)  
If `Stream` is associated to a file, unify `Atom` to the name of this file.

**file_no**(`Integer`)  
If the stream is associated with a POSIX file descriptor, unify `Integer` with the descriptor number. SWI-Prolog extension used primarily for integration with foreign code. See also [Sfileno()](foreign-streams.html#Sfileno()) from `SWI-Stream.h`.

**input**  
True if `Stream` has mode `read`.

**locale**(`Locale`)  
True when `Locale` is the current locale associated with the stream. See [section 4.23](locale.html#sec:4.23).

**mode**(`IOMode`)  
Unify `IOMode` to the mode given to [open/4](IO.html#open/4) for opening the stream. Values are: `read`, `write`, `append` and the SWI-Prolog extension `update`.

**newline**(`NewlineMode`)  
One of `posix` or `dos`. If `dos`, text streams will emit `\r\n` for `\n` and discard `\r` from input streams. Default depends on the operating system.

**nlink**(`-Count`)  
Number of hard links to the file. This expresses the number of‘names’the file has. Not supported on all operating systems and the value might be bogus. See the documentation of **fstat()** for your OS and the value `st_nlink`.

**output**  
True if `Stream` has mode `write`, `append` or `update`.

**position**(`Pos`)  
Unify `Pos` with the current stream position. A stream position is an opaque term whose fields can be extracted using [stream_position_data/3](IO.html#stream_position_data/3). See also [set_stream_position/2](IO.html#set_stream_position/2).

**reposition**(`Bool`)  
Unify `Bool` with `true` if the position of the stream can be set (see [seek/4](IO.html#seek/4)). It is assumed the position can be set if the stream has a *seek-function* and is not based on a POSIX file descriptor that is not associated to a regular file.

**representation_errors**(`Mode`)  
Determines behaviour of character output if the stream cannot represent a character. For example, an ISO Latin-1 stream cannot represent Cyrillic characters. The behaviour is one of `error` (throw an I/O error exception), `prolog` (write `\x<hex>\`), `unicode` (write `\uXXXX` or `\UXXXXXXXX` escape sequences) or `xml` (write `&#...;` XML character entity). The initial mode is `unicode` for the user streams and `error` for all other streams. See also [section 2.18.1](widechars.html#sec:2.18.1) and [set_stream/2](IO.html#set_stream/2).

**timeout**(`-Time`)  
`Time` is the timeout currently associated with the stream. See [set_stream/2](IO.html#set_stream/2) with the same option. If no timeout is specified, `Time` is unified to the atom `infinite`.

**type**(`Type`)  
Unify `Type` with `text` or `binary`.

**tty**(`Bool`)  
This property is reported with `Bool` equal to `true` if the stream is associated with a terminal. See also [set_stream/2](IO.html#set_stream/2).

**unicode_atoms**(`Mode`)  
`Mode` is the per-stream atom-content policy in effect for the term reader: one of `accept`, `nfc`, `error` or `reject`. See [read_term/2](termrw.html#read_term/2),3 for the meaning of each value.

**write_errors**(`Atom`)  
`Atom` is one of `error` (default) or `ignore`. The latter is intended to deal with service processes for which the standard output handles are not connected to valid streams. In these cases write errors may be ignored on `user_error`.

\[deprecated\]**current_stream**(`?Object, ?Mode, ?Stream`)  
The predicate [current_stream/3](IO.html#current_stream/3) is used to access the status of a stream as well as to generate all open streams. `Object` is the name of the file opened if the stream refers to an open file, an integer file descriptor if the stream encapsulates an operating system stream, or the atom `[]` if the stream refers to some other object. `Mode` is one of `read` or `write`.

This predicate is deprecated. New code should use the ISO predicate [stream_property/2](IO.html#stream_property/2).

**is_stream**(`+Term`)  
True if `Term` is a stream name or valid stream handle. This predicate realises a safe test for the existence of a stream alias or handle.

**stream_pair**(`?StreamPair, ?Read, ?Write`)  
This predicate can be used in mode (-,+,+) to create a *stream-pair* from an input stream and an output stream. Mode (+,-,-) can be used to get access to the underlying streams. If a stream has already been closed, the corresponding argument is left unbound. If mode (+,-,-) is used on a single stream, either `Read` or `Write` is unified with the stream while the other argument is left unbound. This behaviour simplifies writing code that must operate both on streams and stream pairs.

Stream-pairs can be used by all I/O operations on streams, where the operation selects the appropriate member of the pair. The predicate [close/1](IO.html#close/1) closes the still open streams of the pair.^(101As of version 7.1.19, it is allowed to close one of the members of the stream directly and close the pair later.) The output stream is closed before the input stream. If closing the output stream results in an error, the input stream is still closed. Success is only returned if both streams were closed successfully.

\[ISO\]**set_stream_position**(`+Stream, +Pos`)  
Set the current position of `Stream` to `Pos`. `Pos` is a term as returned by [stream_property/2](IO.html#stream_property/2) using the `position(Pos)` property. See also [seek/4](IO.html#seek/4).

**stream_position_data**(`?Field, +Pos, -Data`)  
Extracts information from the opaque stream position term as returned by [stream_property/2](IO.html#stream_property/2) requesting the `position(Pos)` property. `Field` is one of `line_count`, `line_position`, `char_count` or `byte_count`. See also [line_count/2](streamstat.html#line_count/2), [line_position/2](streamstat.html#line_position/2), [character_count/2](streamstat.html#character_count/2) and [byte_count/2](streamstat.html#byte_count/2).^(102Introduced in version 5.6.4 after extending the position term with a byte count. Compatible with SICStus Prolog.)

**seek**(`+Stream, +Offset, +Method, -NewLocation`)  
Reposition the current point of the given `Stream`. `Method` is one of `bof`, `current` or `eof`, indicating positioning relative to the start, current point or end of the underlying object. `NewLocation` is unified with the new offset, relative to the start of the stream.

Positions are counted in‘units’. A unit is 1 byte, except for text files using 2-byte Unicode encoding (2 bytes) or *wchar* encoding (sizeof(wchar_t)). The latter guarantees comfortable interaction with wide-character text objects. Otherwise, the use of [seek/4](IO.html#seek/4) on non-binary files (see [open/4](IO.html#open/4)) is of limited use, especially when using multi-byte text encodings (e.g. UTF-8) or multi-byte newline files (e.g. DOS/Windows). On text files, SWI-Prolog offers reliable backup to an old position using [stream_property/2](IO.html#stream_property/2) and [set_stream_position/2](IO.html#set_stream_position/2). Skipping `N` character codes is achieved calling [get_code/2](chario.html#get_code/2) `N` times or using [copy_stream_data/3](chario.html#copy_stream_data/3), directing the output to a null stream (see [open_null_stream/1](IO.html#open_null_stream/1)). If the seek modifies the current location, the line number and character position in the line are set to 0.

If the stream cannot be repositioned, a `permission_error` is raised. If applying the offset would result in a file position less than zero, a `domain_error` is raised. Behaviour when seeking to positions beyond the size of the underlying object depend on the object and possibly the operating system. The predicate [seek/4](IO.html#seek/4) is compatible with Quintus Prolog, though the error conditions and signalling is ISO compliant. See also [stream_property/2](IO.html#stream_property/2) and [set_stream_position/2](IO.html#set_stream_position/2).

**set_stream**(`+Stream, +Attribute`)  
Modify an attribute of an existing stream. `Attribute` specifies the stream property to set. If stream is a *pair* (see [stream_pair/3](IO.html#stream_pair/3)) both streams are modified, unless the property is only meaningful on one of the streams or setting both is not meaningful. In particular, `eof_action` only applies to the *read* stream, `representation_errors` only applies to the *write* stream and trying to set `alias` or `line_position` on a pair results in a `permission_error` exception. See also [stream_property/2](IO.html#stream_property/2) and [open/4](IO.html#open/4).

**alias**(`AliasName`)  
Set the alias of an already created stream. If `AliasName` is the name of one of the standard streams, this stream is rebound. Thus, `set_stream(S, current_input)` is the same as [set_input/1](IO.html#set_input/1), and by setting the alias of a stream to `user_input`, etc., all user terminal input is read from this stream. See also interactor/0.

**buffer**(`Buffering`)  
Set the buffering mode of an already created stream. Buffering is one of `full`, `line` or `false`.

**buffer_size**(`+Size`)  
Set the size of the I/O buffer of the underlying stream to `Size` bytes.

**close_on_abort**(`Bool`)  
Determine whether or not the stream is closed by [abort/0](toplevel.html#abort/0). By default, streams are closed.

**close_on_exec**(`Bool`)  
Set the `close_on_exec` property. See [stream_property/2](IO.html#stream_property/2).

**encoding**(`Atom`)  
Defines the mapping between bytes and character codes used for the stream. See [section 2.18.1](widechars.html#sec:2.18.1) for supported encodings. The value `bom` causes the stream to check whether the current character is a Unicode BOM marker. If a BOM marker is found, the encoding is set accordingly and the call succeeds. Otherwise the call fails.

**eof_action**(`Action`)  
Set end-of-file handling to one of `eof_code`, `reset` or `error`.

**file_name**(`FileName`)  
Set the filename associated to this stream. This call can be used to set the file for error locations if `Stream` corresponds to `FileName` and is not obtained by opening the file directly but, for example, through a network service.

**line_position**(`LinePos`)  
Set the line position attribute of the stream. This feature is intended to correct position management of the stream after sending a terminal escape sequence (e.g., setting ANSI character attributes). Setting this attribute raises a permission error if the stream does not record positions. See [line_position/2](streamstat.html#line_position/2) and [stream_property/2](IO.html#stream_property/2) (property `position`).

**locale**(`+Locale`)  
Change the locale of the stream. See [section 4.23](locale.html#sec:4.23).

**newline**(`NewlineMode`)  
Set input or output translation for newlines. See corresponding [stream_property/2](IO.html#stream_property/2) for details. In addition to the detected modes, an input stream can be set in mode `detect`. It will be set to `dos` if a `\r` character was removed.

**timeout**(`Seconds`)  
This option can be used to make streams generate an exception if it takes longer than `Seconds` before any new data arrives at the stream. The value `infinite` (default) makes the stream block indefinitely. Like [wait_for_input/3](streamstat.html#wait_for_input/3), this call only applies to streams that support the **select()** system call. For further information about timeout handling, see [wait_for_input/3](streamstat.html#wait_for_input/3). The exception is of the form

> `error(``timeout_error(read, Stream)``, _)`

**type**(`Type`)  
Set the type of the stream to one of `text` or `binary`. See also [open/4](IO.html#open/4) and the `encoding` property of streams. Switching to `binary` sets the encoding to `octet`. Switching to `text` sets the encoding to the default text encoding.

**record_position**(`Bool`)  
Do/do not record the line count and line position (see [line_count/2](streamstat.html#line_count/2) and [line_position/2](streamstat.html#line_position/2)). Calling `set_stream(S, record_position(true))` resets the position the start of line 1.

**representation_errors**(`Mode`)  
Change the behaviour when writing characters to the stream that cannot be represented by the encoding. See also [stream_property/2](IO.html#stream_property/2) and [section 2.18.1](widechars.html#sec:2.18.1).

**tty**(`Bool`)  
Modify whether Prolog thinks there is a terminal (i.e. human interaction) connected to this stream. On Unix systems the initial value comes from **isatty()**. On Windows, the initial user streams are supposed to be associated to a terminal. See also [stream_property/2](IO.html#stream_property/2).

**unicode_atoms**(`Mode`)  
Set the per-stream atom-content policy used by the term reader. `Mode` is one of `accept`, `nfc`, `error` or `reject`; see [read_term/2](termrw.html#read_term/2),3 for the meaning of each value. The Prolog flag [unicode_atoms](flags.html#flag:unicode_atoms) provides the default seeded into newly opened streams.

**set_prolog_IO**(`+In, +Out, +Error`)  
Prepare the given streams for interactive behaviour normally associated to the terminal. `In` becomes the `user_input` and `current_input` of the calling thread. `Out` becomes `user_output` and `current_output`. If `Error` equals `Out` an unbuffered stream is associated to the same destination and linked to `user_error`. Otherwise `Error` is used for `user_error`. Output buffering for `Out` is set to `line` and buffering on `Error` is disabled. See also [prolog/0](toplevel.html#prolog/0) and [set_stream/2](IO.html#set_stream/2). The *clib* package provides the library `library(prolog_server)`, creating a TCP/IP server for creating an interactive session to Prolog.

**set_system_IO**(`+In, +Out, +Error`)  
Bind the given streams to the operating system I/O streams 0-2 using POSIX **dup2()** API. `In` becomes `stdin`. `Out` becomes `stdout`. If `Error` equals `Out` an unbuffered stream is associated to the same destination and linked to `stderr`. Otherwise `Error` is used for `stderr`. Output buffering for `Out` is set to line and buffering on `Error` is disabled. The operating system I/O streams are shared across all threads. The three streams must be related to a *file descriptor* or a `domain_error` `file_stream` is raised. See also [stream_property/2](IO.html#stream_property/2), property `file_no(Fd)`.

Where [set_prolog_IO/3](IO.html#set_prolog_IO/3) rebinds the Prolog streams `user_input`, `user_output` and `user_error` for a specific thread providing a private interactive session, [set_system_IO/3](IO.html#set_system_IO/3) rebinds the shared console I/O and also captures Prolog kernel events (e.g., low-level debug messages, unexpected events) as well as messages from foreign libraries that are directly written to `stdout` or `stderr`.

This predicate is intended to capture all output in situations where standard I/O is normally lost, such as when Prolog is running as a service on Windows.

### 4.17.3 Edinburgh-style I/O

The package for implicit input and output destinations is (almost) compatible with Edinburgh DEC-10 and C-Prolog. The reading and writing predicates refer to, resp., the *current* input and output streams. Initially these streams are connected to the terminal. The current output stream is changed using [tell/1](IO.html#tell/1) or [append/1](IO.html#append/1). The current input stream is changed using [see/1](IO.html#see/1). The stream's current value can be obtained using [telling/1](IO.html#telling/1) for output and [seeing/1](IO.html#seeing/1) for input.

Source and destination are either a file, `user`, or a term‘pipe(`Command`)’. The reserved stream name `user` refers to the terminal.^(103The ISO I/O layer uses `user_input`, `user_output` and `user_error`.) In the predicate descriptions below we will call the source/destination argument‘`SrcDest`’. Below are some examples of source/destination specifications.

|                       |                                  |
|-----------------------|----------------------------------|
| `?- see(data).`       | % Start reading from file‘data’. |
| `?- tell(user).`      | % Start writing to the terminal. |
| `?- tell(pipe(lpr)).` | % Start writing to the printer.  |

Another example of using the `pipe/1` construct is shown below. Note that the `pipe/1` construct is not part of Prolog's standard I/O repertoire. See also process_create/3.

``` code
getwd(Wd) :-
        seeing(Old), see(pipe(pwd)),
        collect_wd(String),
        seen, see(Old),
        atom_codes(Wd, String).

collect_wd([C|R]) :-
        get0(C), C \== -1, !,
        collect_wd(R).
collect_wd([]).
```

The effect of [tell/1](IO.html#tell/1) is not undone on backtracking, and since the stream handle is not specified explicitly in further I/O operations when using Edinburgh-style I/O, you may write to unintended streams more easily than when using ISO compliant I/O. For example, the following query writes both "a" and "b" into the file‘out’:

``` code
?- (tell(out), write(a), false ; write(b)), told.
```

#### Compatibility notes

Unlike Edinburgh Prolog systems, [telling/1](IO.html#telling/1) and [seeing/1](IO.html#seeing/1) do not return the filename of the current input/output but rather the stream identifier, to ensure the design pattern below works under all circumstances:^(104Filenames can be ambiguous and SWI-Prolog streams can refer to much more than just files.)

``` code
        ...,
        telling(Old), tell(x),
        ...,
        told, tell(Old),
        ...,
```

The predicates [tell/1](IO.html#tell/1) and [see/1](IO.html#see/1) first check for `user`, the `pipe(command)` and a stream handle. Otherwise, if the argument is an atom it is first compared to open streams associated to a file with *exactly* the same name. If such a stream exists, created using [tell/1](IO.html#tell/1) or [see/1](IO.html#see/1), output (input) is switched to the open stream. Otherwise a file with the specified name is opened.

The behaviour is compatible with Edinburgh Prolog. This is not without problems. Changing directory, non-file streams, and multiple names referring to the same file easily lead to unexpected behaviour. New code, especially when managing multiple I/O channels, should consider using the ISO I/O predicates defined in [section 4.17.2](IO.html#sec:4.17.2).

**see**(`+SrcDest`)  
Open `SrcDest` for reading and make it the current input (see [set_input/1](IO.html#set_input/1)). If `SrcDest` is a stream handle, just make this stream the current input. See the introduction of [section 4.17.3](IO.html#sec:4.17.3) for details.

**tell**(`+SrcDest`)  
Open `SrcDest` for writing and make it the current output (see [set_output/1](IO.html#set_output/1)). If `SrcDest` is a stream handle, just make this stream the current output. See the introduction of [section 4.17.3](IO.html#sec:4.17.3) for details.

**append**(`+File`)  
Similar to [tell/1](IO.html#tell/1), but positions the file pointer at the end of `File` rather than truncating an existing file. The pipe construct is not accepted by this predicate.

**seeing**(`?SrcDest`)  
Same as [current_input/1](IO.html#current_input/1), except that `user` is returned if the current input is the stream `user_input` to improve compatibility with traditional Edinburgh I/O. See the introduction of [section 4.17.3](IO.html#sec:4.17.3) for details.

**telling**(`?SrcDest`)  
Same as [current_output/1](IO.html#current_output/1), except that `user` is returned if the current output is the stream `user_output` to improve compatibility with traditional Edinburgh I/O. See the introduction of [section 4.17.3](IO.html#sec:4.17.3) for details.

**seen**  
Close the current input stream. The new input stream becomes `user_input`.

**told**  
Close the current output stream. The new output stream becomes `user_output`.

### 4.17.4 Switching between Edinburgh and ISO I/O

The predicates below can be used for switching between the implicit and the explicit stream-based I/O predicates.

\[ISO\]**set_input**(`+Stream`)  
Set the current input stream to become `Stream`. Thus, `open(file, read, Stream), set_input(Stream)` is equivalent to `see(file)`.

\[ISO\]**set_output**(`+Stream`)  
Set the current output stream to become `Stream`. See also [with_output_to/2](IO.html#with_output_to/2).

\[ISO\]**current_input**(`-Stream`)  
Get the current input stream. Useful for getting access to the status predicates associated with streams.

\[ISO\]**current_output**(`-Stream`)  
Get the current output stream.

### 4.17.5 Adding IRI schemas

The file handling predicates may be *hooked* to deal with *IRIs*. An IRI starts with \<`scheme`\>`://`, where \<`scheme`\> is a non-empty sequence of lowercase ASCII letters. After detecting the scheme the file manipulation predicates call a hook that is registered using [register_iri_scheme/3](IO.html#register_iri_scheme/3).

Hooking the file operations using extensible IRI schemas allows us to place any resource that is accessed through Prolog I/O predicates on arbitrary devices such as web servers or the ZIP archive used to store program resources (see [section 14.2](saved-states.html#sec:14.2)). This is typically combined with [file_search_path/2](consulting.html#file_search_path/2) declarations to switch between accessing a set of resources from local files, from the program resource database, from a web-server, etc.

**register_iri_scheme**(`+Scheme, :Hook, +Options`)  
Register `Hook` to be called by all file handling predicates if a name that starts with `Scheme`:// is encountered. The `Hook` is called by call/4 using the *operation*, the `IRI` and a term that receives the *result* of the operation. The following operations are defined:

**open**(`Mode,Options`)  
Called by [open/3](IO.html#open/3),4. The result argument must be unified with a stream.

**access**(`Mode`)  
Called by [access_file/2](files.html#access_file/2), [exists_file/1](files.html#exists_file/1) (`Mode` is `file`) and [exists_directory/1](files.html#exists_directory/1) (`Mode` is `directory`). The result argument must be unified with a boolean.

**time**  
Called by [time_file/2](files.html#time_file/2). The result must be unified with a time stamp.

**size**  
Called by [size_file/2](files.html#size_file/2). The result must be unified with an integer representing the size in bytes.

### 4.17.6 Write onto atoms, code-lists, etc.

**with_output_to**(`+Output, :Goal`)  
Run `Goal` as [once/1](metacall.html#once/1), while characters written to the current output are sent to `Output`. The predicate is SWI-Prolog-specific, inspired by various posts to the mailinglist. It provides a flexible replacement for predicates such as sformat/3 , [swritef/3](writef.html#swritef/3), [term_to_atom/2](manipatom.html#term_to_atom/2), [atom_number/2](manipatom.html#atom_number/2) converting numbers to atoms, etc. The predicate [format/3](format.html#format/3) accepts the same terms as output argument.

For capturing other streams, see with_output_to/3.

Applications should generally avoid creating atoms by breaking and concatenating other atoms, as the creation of large numbers of intermediate atoms generally leads to poor performance, even more so in multithreaded applications. This predicate supports creating difference lists from character data efficiently. The example below defines the DCG rule term//1 to insert a term in the output:

``` code
term(Term, In, Tail) :-
        with_output_to(codes(In, Tail), write(Term)).

?- phrase(term(hello), X).

X = [104, 101, 108, 108, 111]
```

`Output` takes one of the shapes below. Except for the first, the system creates a temporary stream using the `wchar_t` internal encoding that points at a memory buffer. The encoding cannot be changed and an attempt to call [set_stream/2](IO.html#set_stream/2) using `encoding(Encoding)` results in a `permission_error` exception.

**A Stream handle or alias**  
Temporarily switch current output to the given stream. Redirection using [with_output_to/2](IO.html#with_output_to/2) guarantees the original output is restored, also if `Goal` fails or raises an exception. See also [call_cleanup/2](metacall.html#call_cleanup/2).

**atom**(`-Atom`)  
Create an atom from the emitted characters. Please note the remark above.

**string**(`-String`)  
Create a string object as defined in [section 5.2](string.html#sec:5.2).

**codes**(`-Codes`)  
Create a list of character codes from the emitted characters, similar to [atom_codes/2](manipatom.html#atom_codes/2).

**codes**(`-Codes, -Tail`)  
Create a list of character codes as a difference list.

**chars**(`-Chars`)  
Create a list of one-character atoms from the emitted characters, similar to [atom_chars/2](manipatom.html#atom_chars/2).

**chars**(`-Chars, -Tail`)  
Create a list of one-character atoms as a difference list.

### 4.17.7 Fast binary term I/O

The predicates in this section provide fast binary I/O of arbitrary Prolog terms, including cyclic terms and terms holding attributed variables. Library `library(fastrw)` is a SICSTus/Ciao compatible library that extends the core primitives described below.

The binary representation the same as used by [PL_record_external()](foreigninclude.html#PL_record_external()). The use of these primitives instead of using [write_canonical/2](termrw.html#write_canonical/2) has advantages and disadvantages. Below are the main considerations:

- Using [write_canonical/2](termrw.html#write_canonical/2) allows or exchange of terms with other Prolog systems. The format is stable and, as it is text based, it can be inspected and corrected.
- Using the binary format improves the performance roughly 3 times.
- The size of both representations is comparable.
- The binary format can deal with cycles, sharing and attributes. Special precautions are needed to transfer such terms using [write_canonical/2](termrw.html#write_canonical/2). See [term_factorized/3](terms.html#term_factorized/3) and [copy_term/3](attvar.html#copy_term/3).
- In the current version, reading the binary format has only incomplete consistency checks. This implies a user must be able to **trust the source** as crafted messages may compromise the reading Prolog system.

**fast_term_serialized**(`?Term, ?String`)  
(De-)serialize `Term` to/from `String`.

**fast_write**(`+Output, +Term`)  
Write `Term` using the fast serialization format to the `Output` stream. `Output` *must* be a binary stream.

**fast_read**(`+Input, -Term`)  
Read `Term` using the fast serialization format from the `Input` stream. `Input` *must* be a binary stream.^(bugThe predicate [fast_read/2](IO.html#fast_read/2) may crash on arbitrary input.)
