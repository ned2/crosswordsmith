
## 4.42 Execution profiling

This section describes the hierarchical execution profiler. This profiler is based on ideas from **gprof** described in [Graham *et al.*, 1982](Bibliography.html#graham82gprof). The profiler consists of two parts: the information-gathering component built into the kernel,^(163There are two implementations; one based on **setitimer()** using the `SIGPROF` signal and one using Windows Multi Media (MM) timers. On other systems the profiler is not provided.) and a presentation component which is defined in the `library(statistics)` library. The latter can be hooked, which is used by the XPCE module `library(swi/pce_profile)` to provide an interactive graphical frontend for the results.

### 4.42.1 library(prolog_profile): Execution profiler

This module provides a simple frontend on the execution profiler with a hook to the GUI visualiser for profiling results defined in `library(swi/pce_profile)`.

**profile**(`:Goal`)  
**profile**(`:Goal, +Options`)  
Run `once(Goal)` under the execution profiler. If the (xpce) GUI is enabled this predicate is hooked by `library(swi/pce_profile)` and results are presented in a gui that enables navigating the call tree and jump to predicate implementations. Without the GUI, a simple textual report is generated. Defined options are:

**time**(`Which`)  
Profile `cpu` or `wall` time. The default is CPU time.

**sample_rate**(`Rate`)  
Samples per second, any numeric value between 1 and 1000. Default is defined by the Prolog flag **profile_sample_rate**, which defaults to 200.

**ports**(`Bool`)  
Specifies ports counted - `true` (all ports), `false` (call port only) or `classic` (all with some errors). Accomodates space/accuracy tradeoff building call tree. Default is defined by the Prolog flag **profile_ports**, which defaults to `true`.

**top**(`N`)  
When generating a textual report, show the top `N` predicates.

**cumulative**(`Bool`)  
If `true` (default `false`), show cumulative output in a textual report.

See also  
show_coverage/2 from `library(test_cover)`.

To be done  
The textual input reflects only part of the information.

**show_profile**(`+Options`)  
Display last collected profiling data. `Options` are

**top**(`N`)  
When generating a textual report, show the top `N` predicates.

**cumulative**(`Bool`)  
If `true` (default `false`), show cumulative output in a textual report.

\[det\]**profile_data**(`-Data`)  
Gather all relevant data from profiler. This predicate may be called while profiling is active in which case it is suspended while collecting the data. `Data` is a dict providing the following fields:

`summary`**`:`**`Dict`  
Overall statistics providing

- samples:Count: Times the statistical profiler was called
- ticks:Count Virtual ticks during profiling
- accounting:Count Tick spent on accounting
- time:Seconds Total time sampled
- nodes:Count Nodes in the call graph.
- sample_period: MicroSeconds Same interval timer period in micro seconds
- ports: Ports One of `true`, `false` or `classic`

**nodes**  
List of nodes. Each node provides:

- predicate:PredicateIndicator
- ticks_self:Count
- ticks_siblings:Count
- call:Count
- redo:Count
- exit:Count
- callers:`list_of(Relative)`
- callees:`list_of(Relative)`

*Relative* is a term of the shape below that represents a caller or callee. Future versions are likely to use a dict instead.

``` code
node(PredicateIndicator, CycleID, Ticks, TicksSiblings,
     Calls, Redos, Exits)
```

\[nondet\]**profile_procedure_data**(`?Pred, -Data:dict`)  
Collect data for `Pred`. If `Pred` is unbound data for each predicate that has profile data available is returned. `Data` is described in [profile_data/1](profile.html#profile_data/1) as an element of the `nodes` key.

### 4.42.2 Visualizing profiling data

Browsing the annotated call-tree as described in [section 4.42.3](profile.html#sec:4.42.3) itself is not very attractive. Therefore, the results are combined per predicate, collecting all *callers* and *callees* as well as the propagation of time and activations in both directions. [Figure 4](profile.html#fig:profnode) illustrates this. The central yellowish line is the‘current’predicate with counts for time spent in the predicate (\`Self’), time spent in its children (\`Siblings’), activations through the call and redo ports. Above that are the *callers*. Here, the two time fields indicate how much time is spent serving each of the callers. The columns sum to the time in the yellowish line. The caller *`<`recursive`>`* is the number of recursive calls. Below the yellowish lines are the callees, with the time spent in the callee itself for serving the current predicate and the time spent in the callees of the callee ('Siblings’), so the whole time-block adds up to the‘Siblings’field of the current predicate. The‘Access’fields show how many times the current predicate accesses each of the callees.

The predicates have a menu that allows changing the view of the detail window to the given caller or callee, showing the documentation (if it is a built-in) and/or jumping to the source.

![](profnode.png)

**Figure 4 :** Execution profiler showing the activity of the predicate chat:inv_map_list/5.

The statistics shown in the report field of [figure 4](profile.html#fig:profnode) show the following information:

- *samples*  
  Number of times the call-tree was sampled for collecting time statistics. On most hardware, the resolution of `SIGPROF` is 1/100 second. This number must be sufficiently large to get reliable timing figures. The **Time** menu allows viewing time as samples, relative time or absolute time.
- *sec*  
  Total user CPU time with the profiler active.
- *predicates*  
  Total count of predicates that have been called at least one time during the profile.
- *nodes*  
  Number of nodes in the call-tree.
- *distortion*  
  How much of the time is spent building the call-tree as a percentage of the total execution time. Timing samples while the profiler is building the call-tree are not added to the call-tree.

### 4.42.3 Information gathering

While the program executes under the profiler, the system builds a *dynamic* call-tree. It does this using three hooks from the kernel: one that starts a new goal (*profCall*), one that tells the system which goal is resumed after an *exit* (*profExit*) and one that tells the system which goal is resumed after a *fail* (i.e., which goal is used to *retry* (*profRedo*)). The **profCall()** function finds or creates the subnode for the argument predicate below the current node, increments the call-count of this link and returns the sub-node which is recorded in the Prolog stack-frame. Choice-points are marked with the current profiling node. **profExit()** and **profRedo()** pass the profiling node where execution resumes.

Just using the above algorithm would create a much too big tree due to recursion. For this reason the system performs detection of recursion. In the simplest case, recursive procedures increment the‘recursive’count on the current node. Mutual recursion, however, is not easily detected. For example, [call/1](metacall.html#call/1) can call a predicate that uses [call/1](metacall.html#call/1) itself. This can be viewed as a recursive invocation, but this is generally not desirable. Recursion is currently assumed if the same predicate *with the same parent* appears higher in the call-graph. Early experience with some non-trivial programs are promising.

The last part of the profiler collects statistics on the CPU time used in each node. On systems providing **setitimer()** with `SIGPROF`, it‘ticks’the current node of the call-tree each time the timer fires. On Windows, a MM-timer in a separate thread checks 100 times per second how much time is spent in the profiled thread and adds this to the current node. See [section 4.42.3.1](profile.html#sec:4.42.3.1) for details.

#### 4.42.3.1 Profiling in the Windows Implementation

Profiling in the Windows version is similar, but as profiling is a statistical process it is good to be aware of the implementation^(164We hereby acknowledge Lionel Fourquaux, who suggested the design described here after a newsnet enquiry.) for proper interpretation of the results.

Windows does not provide timers that fire asynchronously, frequent and proportional to the CPU time used by the process. Windows does provide multi-media timers that can run at high frequency. Such timers, however, run in a separate thread of execution and they are fired on the wall clock rather than the amount of CPU time used. The profiler installs such a timer running, for saving CPU time, rather inaccurately at about 100 Hz. Each time it is fired, it determines the CPU time in milliseconds used by Prolog since the last time it was fired. If this value is non-zero, active predicates are incremented with this value.
