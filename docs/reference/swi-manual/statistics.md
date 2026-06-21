
## A.54 library(statistics): Get information about resource usage

This library provides predicates to obtain information about resource usage by your program. The predicates of this library are for human use at the toplevel: information is *printed*. All predicates obtain their information using public low-level primitives. These primitives can be use to obtain selective statistics during execution.

\[det\]**statistics**  
Print information about resource usage using [print_message/2](printmsg.html#print_message/2).

See also  
All statistics printed are obtained through [statistics/2](builtin-statistics.html#statistics/2).

\[det\]**statistics**(`-Stats:dict`)  
`Stats` is a dict representing the same information as [statistics/0](statistics.html#statistics/0). This convience function is primarily intended to pass statistical information to e.g., a web client. Time critical code that wishes to collect statistics typically only need a small subset and should use [statistics/2](builtin-statistics.html#statistics/2) to obtain exactly the data they need.

\[nondet\]**thread_statistics**(`?Thread, -Stats:dict`)  
Obtain statistical information about a single thread. Fails silently of the `Thread` is no longer alive.

|  |  |
|----|----|
| `Stats` | is a dict containing status, time and stack-size information about `Thread`. |

\[nondet\]**time**(`:Goal`)  
Execute `Goal`, reporting statistics to the user. If `Goal` succeeds non-deterministically, retrying reports the statistics for providing the next answer.

Note that is no portable way to get thread-specific CPU time. SWI-Prolog has implementations for Linux, Windows and MacOS. The automatic detection may work on some other operating systems.

See also  
\- [statistics/2](builtin-statistics.html#statistics/2) for obtaining statistics in your program and understanding the reported values.  
- [call_time/2](statistics.html#call_time/2), [call_time/3](statistics.html#call_time/3) to obtain the timing in a dict.

bug  
Inference statistics are often a few off.

**call_time**(`:Goal, -Time:dict`)  
**call_time**(`:Goal, -Time:dict, -Result`)  
Call `Goal` as [call/1](metacall.html#call/1), unifying `Time` with a dict that provides information on the resource usage. If `Goal` succeeds with a choice point, backtracking reports the time used to find the *next answer*, failure or exception. If `Goal` succeeds deterministically no choice point is left open. Currently `Time` contains the keys below. Future versions may provide additional keys.

- wall:Seconds
- cpu:Seconds
- inferences:Count

[call_time/2](statistics.html#call_time/2) is defined as below. Note that for [call_time/2](statistics.html#call_time/2) the time is only available if `Goal` succeeds.

``` code
call_time(Goal, Time) :-
    call_time(Goal, Time, Result),
    call(Result).
```

|  |  |
|----|----|
| `Result` | is one of `true`, `false` or `throw(E)`, depending on whether or not the goal succeeded or raised an exception. Note that `Result` may be called using [call/1](metacall.html#call/1) to propagate the failure or exception. |
