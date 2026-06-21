
## 4.18 Status of streams

\[det\]**wait_for_input**(`+ListOfStreams, -ReadyList, +TimeOut`)  
Wait for input on one of the streams in `ListOfStreams` and return a list of streams on which input is available in `ReadyList`. Each element of `ListOfStreams` is either a stream or an integer. Integers are consider waitable OS handles. This can be used to (also) wait for handles that are not associated with Prolog streams such as UDP sockets. See tcp_setopt/2.

This predicate waits for at most `TimeOut` seconds. `TimeOut` may be specified as a floating point number to specify fractions of a second. If `TimeOut` equals `infinite`, [wait_for_input/3](streamstat.html#wait_for_input/3) waits indefinitely. If `Timeout` is 0 or 0.0 this predicate returns without waiting.^(105Prior to 7.3.23, the integer value‘0’was the same as `infinite`.)

This predicate can be used to implement timeout while reading and to handle input from multiple sources and is typically used to wait for multiple (network) sockets. On Unix systems it may be used on any stream that is associated with a system file descriptor. On Windows it can only be used on sockets. If `ListOfStreams` contains a stream that is not associated with a supported device, a `domain_error(waitable_stream, Stream)` is raised.

The example below waits for input from the user and an explicitly opened secondary terminal stream. On return, `Inputs` may hold `user_input` or `P4` or both.

``` code
?- open('/dev/ttyp4', read, P4),
   wait_for_input([user_input, P4], Inputs, 0).
```

When available, the implementation is based on the **poll()** system call. The **poll()** puts no additional restriction on the number of open files the process may have. It does limit the time to `2^31-1` milliseconds (a bit less than 25 days). Specifying a too large timeout raises a `representation_error(timeout)` exception. If **poll()** is not supported by the OS, **select()** is used. The **select()** call can only handle file descriptors up to `FD_SETSIZE`. If the set contains a descriptor that exceeds this limit a `representation_error(’FD_SETSIZE’)` is raised.

Note that [wait_for_input/3](streamstat.html#wait_for_input/3) returns streams that have data waiting. This does not mean you can, for example, call [read/2](termrw.html#read/2) on the stream without blocking as the stream might hold an incomplete term. The predicate [set_stream/2](IO.html#set_stream/2) using the option `timeout(Seconds)` can be used to make the stream generate an exception if no new data arrives within the timeout period. Suppose two processes communicate by exchanging Prolog terms. The following code makes the server immune for clients that write an incomplete term:

``` code
    ...,
    tcp_accept(Server, Socket, _Peer),
    tcp_open(Socket, In, Out),
    set_stream(In, timeout(10)),
    catch(read(In, Term), _, (close(Out), close(In), fail)),
    ...,
```

**byte_count**(`+Stream, -Count`)  
Byte position in `Stream`. For binary streams this is the same as [character_count/2](streamstat.html#character_count/2). For text files the number may be different due to multi-byte encodings or additional record separators (such as Control-M in Windows).

**character_count**(`+Stream, -Count`)  
Unify `Count` with the current character index. For input streams this is the number of characters read since the open; for output streams this is the number of characters written. Counting starts at 0.

**line_count**(`+Stream, -Count`)  
Unify `Count` with the number of lines read or written. Counting starts at 1.

**line_position**(`+Stream, -Count`)  
Unify `Count` with the position on the current line. Note that this assumes the position is 0 after the open. Tabs are assumed to be defined on each 8-th character, and backspaces are assumed to reduce the count by one, provided it is positive.
