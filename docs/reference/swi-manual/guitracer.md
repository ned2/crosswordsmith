
## 3.5 The Graphical Debugger

SWI-Prolog offers two debuggers. One is the traditional text console-based 4-port Prolog tracer and the other is a window-based source level debugger. The window-based debugger requires XPCE installed. It operates based on the [prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4) hook and other low-level functionality described in [chapter B](hack.html#sec:B).

Window-based tracing provides a much better overview due to the eminent relation to your source code, a clear list of named variables and their bindings as well as a graphical overview of the call and choice point stack. There are some drawbacks though. Using a textual trace on the console, one can scroll back and examine the past, while the graphical debugger just presents a (much better) overview of the current state.

### 3.5.1 Invoking the window-based debugger

Whether the text-based or window-based debugger is used is controlled using the predicates [guitracer/0](guitracer.html#guitracer/0) and [noguitracer/0](guitracer.html#noguitracer/0). Entering debug mode is controlled using the normal predicates for this: [trace/0](debugger.html#trace/0) and [spy/1](debugger.html#spy/1). In addition, PceEmacs prolog mode provides the command **Prolog/Break at** (`Control-c b`) to insert a break-point at a specific location in the source code.

The graphical tracer is particularly useful for debugging threads. The tracer must be loaded from the `main` thread before it can be used from a background thread.

**guitracer**  
This predicate installs the above-mentioned hooks that redirect tracing to the window-based environment. No window appears. The debugger window appears as actual tracing is started through [trace/0](debugger.html#trace/0), by hitting a spy point defined by [spy/1](debugger.html#spy/1) or a break point defined using the PceEmacs command **Prolog/Break at** (`Control-c b`).

**noguitracer**  
Disable the hooks installed by [guitracer/0](guitracer.html#guitracer/0), reverting to normal text console-based tracing.

**gtrace**  
Utility defined as `guitracer,trace`.

**gdebug**  
Utility defined as `guitracer,debug`.

**gspy**(`+Predicate`)  
Utility defined as `guitracer,spy(Predicate)`.
