
## 4.1 Notation of Predicate Descriptions

We have tried to keep the predicate descriptions clear and concise. First, the predicate name is printed in **bold face**, followed by the arguments in *italics*. Arguments are preceded by a *mode indicator*.

### 4.1.1 The argument mode indicator

An *argument mode indicator* gives information about the intended direction in which information carried by a predicate argument is supposed to flow. Mode indicators (and types) are not a formal part of the Prolog language but help in explaining intended semantics to the programmer. There is no complete agreement on argument mode indicators in the Prolog community. We use the following definitions:^(51These definitions are taken from the *PlDoc* markup language description. *PldDoc* markup is used for source code markup (as well as for the commenting tool). The current manual has only one mode declaration per predicate and therefore predicates with mode (`+`,`-`) and (`-`,`+`) are described as (`?`,`?`). The `@`-mode is often replaced by  
chr`+`.)

|  |  |
|----|----|
| `++` | At call time, the argument must be *ground*, i.e., the argument may not contain any variables that are still unbound. |
| `+` | At call time, the argument must be instantiated to a term satisfying some (informal) type specification. The argument need not necessarily be ground. For example, the term `[_]` is a list, although its only member is the anonymous variable, which is always unbound (and thus nonground). |
| `-` | Argument is an *output* argument. It may or may not be bound at call-time. If the argument is bound at call time, the goal behaves as if the argument were unbound, and then unified with that term after the goal succeeds. This is what is called being *steadfast*: instantiation of output arguments at call-time does not change the semantics of the predicate, although optimizations may be performed. For example, the goal `findall(X, Goal, [T])` is good style and equivalent to `findall(X, Goal, Xs), Xs = [T]`^(52The ISO standard dictates that `findall(X, Goal, 1)` raise a `type_error` exception, breaking steadfastness. SWI-Prolog does not follow the standard here.) Note that any *determinism* specification, e.g., `det`, only applies if the argument is unbound. For the case where the argument is bound or involved in constraints, `det` effectively becomes `semidet`, and `multi` effectively becomes `nondet`. |
| `--` | At call time, the argument must be unbound. This is typically used by predicates that create‘something’and return a handle to the created object, such as [open/3](IO.html#open/3), which creates a *stream*. |
| `?` | At call time, the argument must be bound to a *partial term* (a term which may or may not be ground) satisfying some (informal) type specification. Note that an unbound variable *is* a partial term. Think of the argument as either providing input or accepting output or being used for both input and output. For example, in `stream_property(S, reposition(Bool))`, the `reposition` part of the term provides input and the unbound-at-call-time `Bool` variable accepts output. |
| `:` | Argument is a *meta-argument*, for example a term that can be called as goal. The predicate is thus a *meta-predicate*. This flag implies `+`. |
| `@` | Argument will not be further instantiated than it is at call-time. Typically used for type tests. |
| `!` | Argument contains a mutable structure that may be modified using [setarg/3](manipterm.html#setarg/3) or [nb_setarg/3](manipterm.html#nb_setarg/3). |

See also [section 4.8](metacall.html#sec:4.8) for examples of meta-predicates, and [section 6.5](metapred.html#sec:6.5) for mode flags to label meta-predicate arguments in module export declarations.

### 4.1.2 Predicate indicators

Referring to a predicate in running text is done using a *predicate indicator*. The canonical and most generic form of a predicate indicator is a term `[<``module``>:]<``name``>/<``arity``>`. The module is generally omitted if it is irrelevant (case of a built-in predicate) or if it can be inferred from context.

#### 4.1.2.1 Non-terminal indicators

Compliant to the ISO standard draft on Definite Clause Grammars (see [section 4.13](DCG.html#sec:4.13)), SWI-Prolog also allows for the *non-terminal indicator* to refer to a *DCG grammar rule*. The non-terminal indicator is written as `[<``module``>]:<``name``>//<``arity``>`.

A non-terminal indicator `<``name``>//<``arity``>` is understood to be equivalent to `<``name``>/<``arity``>+2`, regardless of whether or not the referenced predicate is defined or can be used as a grammar rule.^(53This, however, makes a specific assumption about the implementation of DCG rules, namely that DCG rules are preprocessed into standard Prolog rules taking two additional arguments, the input list and the output list, in accumulator style. This *need* not be true in all implementations.) The `//`-notation can be used in all places that traditionally allow for a predicate indicator, e.g., the module declaration, [spy/1](debugger.html#spy/1), and [dynamic/1](dynamic.html#dynamic/1).

### 4.1.3 Predicate behaviour and determinism

To describe the general behaviour of a predicate, the following vocabulary is employed. In source code, structured comments contain the corresponding keywords:

|  |  |
|----|----|
| `det` | A *deterministic* predicate always succeeds exactly once and does not leave a choicepoint. |
| `semidet` | A *semi-deterministic* predicate succeeds at most once. If it succeeds it does not leave a choicepoint. |
| `nondet` | A *non-deterministic* predicate is the most general case and no claims are made on the number of solutions (which may be zero, i.e., the predicate may *fail*) and whether or not the predicate leaves an choicepoint on the last solution. |
| `multi` | As `nondet`, but succeeds at least once. |
| `undefined` | Well founded semantics third value. See [undefined/0](WFS.html#undefined/0). |
