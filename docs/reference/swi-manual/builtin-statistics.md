
## 4.41 Obtaining Runtime Statistics

The predicate [statistics/2](builtin-statistics.html#statistics/2) is built-in. More high level predicates are available from library `library(statistics)`. See [section A.54](statistics.html#sec:A.54).

**statistics**(`+Key, -Value`)  
Unify system statistics determined by `Key` with `Value`. The possible keys are given in the [table 6](builtin-statistics.html#tab:statistics). This predicate supports additional keys for compatibility reasons. These keys are described in [table 7](builtin-statistics.html#tab:qpstatistics). CPU time results are based on **clock_gettime()**, **times()** or wall time since the process was started (in that order of preference). On Windows **GetProcessTimes()** is used. Both **clock_gettime()** and **GetProcessTimes()** provide a nanosecond resolution interface. The actual resolution depends on the platform.

Starting with version 9.1.9, the `cputime` and `inferences` keys include the final value for threads that have been created by the calling thread *and* has been *joined* by the calling thread. The new keys `self_cputime` and `self_inferences` may be used to get statistics for the calling thread only. Both keys also exist in the single threaded version, where the “self” key always returns the same value as the one without “self” .

Native keys (times as float in seconds)

agc
