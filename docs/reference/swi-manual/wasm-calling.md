
## 13.2 Calling Prolog from JavaScript

The `Prolog` class provides several methods for calling Prolog from JavaScript.

`Boolean` **Prolog.call**(`Goal`)  
Processes a Prolog goal represented as a `String` and returns `true` or `false`. This simple calling pattern is intended for trivial goals such as setting a Prolog flag. For example, the call below limits the Prolog stacks to 10Mb.

``` code
Prolog.call("set_prolog_flag(stack_limit, 10 000 000)");
```

`Query` **Prolog.query**(`Goal`)  
`Query` **Prolog.query**(`Goal, Input`)  
`Query` **Prolog.query**(`Goal, Input, Options`)  
Create a Prolog query from a `String`, optionally binding Prolog variables embedded in `Goal` from properties of the `Object` `Input`. The returned object is an instance of class `Query`. This instance can be used as a JavaScript *iterator*. The value returned in each iteration is an `Object` with properties for each variable in `Goal` that is not in `Input` and does not start with an underscore. For example, we can iterate over the members of a list like below. Further details on class `Query` are provided in [section 13.2.1](wasm-calling.html#sec:13.2.1). The translation of data between Prolog and JavaScript is described in [section 13.2.3](wasm-calling.html#sec:13.2.3).

``` code
for(let r of Prolog.query("member(Elem,List)",
                          {List: ["aap", "noot", "mies"]}))
{ console.log(r.Elem);
}
```

This interface is also practical for calling (fast) Prolog predicates to compute a single answer from an input using the [Query.once()](wasm-calling.html#Query.once()) method. Assuming a Prolog predicate fib/2 that computes the nth *Fibonacci* number, we can call this using the code below. Note that if the fib/2 fails or raises an exception the object returned by [Query.once()](wasm-calling.html#Query.once()) does not contain the `Out` key and thus our function returns `undefined`.

``` code
function fib(in, out)
{ return Prolog.query("fib(In,Out)", {In:in}).once().Out;
}
```

The `.`[`query()`](wasm-calling.html#query()) method is indented for fast queries that do not require the *yield* mechanism, i.e., the execution should not require asynchronous operations and the browser is not responsive during the execution. Use [Prolog.forEach()](wasm-calling.html#Prolog.forEach()) for asynchronous queries.

The optional `Options` parameter defines the following options

**engine**: `Boolean`  
If `true`, run the query in a temporary engine. Note that JavaScript can only interact with the *innermost query* of an engine. By using a new engine we can interact with multiple queries, using them as *coroutines*. [Prolog.Engine()](wasm-calling.html#Prolog.Engine()) for details.

**string**: `Type`  
Prolog type to use for JavaScript strings when converting the input. Default is `string`. Alternatively one may use `atom`.

**nodebug**: `Boolean`  
If set to `true`, the execution cannot be seen in the debugger.

`Promise` **Prolog.forEach**(`Goal, [Input], [OnAnswer], [Options]`)  
This method executes `Goal` asynchronously. This implies that `Goal` may execute asynchronous operations and the browser remains responsive while executing `Goal`. `Goal` and `Input` are processed as with [Prolog.query()](wasm-calling.html#Prolog.query()). If `OnAnswer` is provided, this `Function` is called with a `Object` that holds the bindings for the output arguments of `Goal` for each answer. `Options` supports the following options

**engine**: `Boolean`  
If `true`, create an engine that will be destroyed when the query completes.

**heartbeat**: `Integer`  
Sets the *heartbeat rate*. This is the number of Prolog inferences executed before yielding. The default is 10,000. Lower values improve interactive behaviour but lower Prolog performance.

All default parameters may be omitted. A JavaScript type check is used to discover whether the `Input` or `OnAnswer` parameter is omitted. Use an empty `Input` to specify `Options` without inputs, e.g. `forEach(goal, {}, {engine:true})`

The returned `Promise` is resolved when the query completes. The value passed to the `.`**`then()`** method of the `Promise` is the number of answers if `OnAnswer` is provided or an `Array` of answers if `OnAnswer` is omitted. If `Goal` raises an exception the `Promise` is rejected.

Multiple calls to Prolog can be activated on the same engine at any time. Prolog processes such queries in *LIFO* *(Last In, First Out)* mode. If queries need to be processed sequentially use JavaScript `await` or the `Promise.`**`finally()`** method to wait for completion. Multiple [Prolog.forEach()](wasm-calling.html#Prolog.forEach()) queries that run in separate engines act as *cooperative threads*, i.e., all make progress. For example, a goal can be started to run as cooperative thread using:

``` code
  setTimeout(async () => {
    await Prolog.forEach(goal, input, onanswer,
                         {engine:true});
  });
```

### 13.2.1 The JavaScript class Query

The method [Prolog.query()](wasm-calling.html#Prolog.query()) returns an instance of the JavaScript class `Query` that may be used to explore the solutions of the query. The `Query` class implements the JavaScript *iterator* protocol.

`Object` **Query.next**()  
Implements the *iterator* protocol. This returns an object holding the keys `done` and `value`. If exception handling is enabled it returns an object {`done`:`true`, `error`:`true`, `message`:`String`}.

`Object` **Query.once**()  
Close the query after the first answer. Returns the `.value` of the object returned by `.`**`next()`** on success and the complete object on failure or error. In addition, on a logical result (no error), a field `success` is added with a boolean value. Thus, the return value may contain these keys:

**{Bindings}**  
Query succeeded. Objects holds a key for each output variable. In addition the `success` key is set to `true`.

**{`success`:`false`}**  
Query failed. Note that the binding keys all start with an uppercase letter and this is thus not ambiguous.

**{`error`:true, `message`:`String`}**  
Query raised an exception.

`Object` **Query.close**()  
Closes the query. This can be used inside the iterator to stop further enumeration.

### 13.2.2 Using engines

The WASM version of SWI-Prolog supports *engines*. The initial engine is called `main`. Additional engines can be used to enumerate answers of multiple open queries as well as for implementing *coroutines*. Combined with JavaScript *async* functions, engines can provide *cooperative multi-threading*. For example, to enumerate the answers of two queries we may use the following code

``` code
  const e1 = new Prolog.Engine({auto_close:true});
  const e2 = new Prolog.Engine({auto_close:true});

  const q1 = e1.query(Query1);  // see Prolog.query()
  const q2 = e2.query(Query2);

  try
  { for(;;)
    { const n1 = q1.next();
      const n2 = q2.next();
      if ( n1.done && n2.done )
        break;
      // Handle answers
    }
  } finally
  { q1.close(); // also closes e1
    q2.close();
  }
  ...
```

Engines can also be used to create a cooperative thread. The example below creates a Prolog task that prints the numbers 1..20 to the console, one number every second.

``` code
  ...
  setTimeout(async () => {
    const e = new Prolog.Engine({auto_close:true});
    await e.forEach("between(1,20,X),sleep(1)",
                    (a) => console.out(a.X));
  });
```

`Engine` **Prolog.Engine**(`[name], [options]`)  
Create a new engine. The `name` argument names the engine. When omitted, engines are named `engine``N`, where `N` is defined by a counter. The `options` argument provides additional configuration. Both arguments are optional.

**Boolean `auto_close`**  
When `true` (default `false`, closing the last query associated with the engine also closes the engine.

`undefined` **close**()  
Terminate the engine. This may be called safely multiple times on the same instance.

`Boolean` **Engine.call**(`Goal`)  
`Object` **query**(`...args`)  
`Object` **forEach**(`goal, ...args`)  
`Any` **with_frame**(`function, persist`)  
Same as the corresponding methods on class **Prolog**, using the specified engine for running the Prolog goals.

`Any` **with**(`function`)  
Run `function` using the specified **Engine** instance.

### 13.2.3 Translating data between JavaScript and Prolog

JavaScript and Prolog are both dynamically typed languages. The WASM module defines a faithful translation between JavaScript data and Prolog data that aims at completeness as well as keeping the data representation clean in the common cases. We describe the translation in two descriptions because round tripping does not always result in the original object.

#### 13.2.3.1 Translating JavaScript data to Prolog

This section describes how data from JavaScript is translated into Prolog. The interface is primarily designed for passing JavaScript data as typically used to a natural Prolog representation. In addition a number of classes are provided to create Prolog specific data structures such as strings (as opposed to atoms), variables, compound terms, etc.

**Number**  
Translate to a Prolog integer or floating point number.

**BigInt**  
Translate to a Prolog integer.

**String**  
Translate to a Prolog atom. Use `new Prolog.`**`String(text)`** to create a Prolog string. See below.

**Boolean**  
Translate to one of the Prolog atoms `true` or `false`.

**undefined**  
Translate the Prolog atom `undefined`.

**null**  
Translate the Prolog atom `null`.

**Array**  
Translate to a Prolog list.

**Objects holding the key `$t`:`Type`**  
Such objects are converted depending on the value for this key. The interface defines classes to simplify creating such objects.

**s**  
Represent a Prolog string. The key `v` holds the text. May be created using `new Prolog.string(text)`. May be created using `new Prolog.`**`String(text)`**.

**r**  
Represent a Prolog *rational number*. The keys `n` and `d` represent the *numerator* and *denominator*. For example, to represent `1r3`, use {`$t`:"r", `n`:1, `d`:3}. May be created using `new Prolog.Rational(n, d)`, where `n` and `d` can be JavaScript numbers or big integers.

**t**  
Represent a Prolog *compound term*. The object should hold exactly one key whose value is an array that holds the argument values. For example a term `point(1,2)` is constructed using {`$t`:"t", `point`:\[1,2\]}. May be created using `new Prolog.Compound(functor, args)`

**v**  
Represent a variable. If the key `v` is present this identifies the variable. Two variables processed in the same translation with the same identifier represent the same Prolog variable. If the `v` key is omitted the variable will be unique. May be created using `new Prolog.Var(id)`.

**l**  
Represent a Prolog list. As a JavaScript `Array` we only need this typed object to create a *partial list*. The `v` key contains the “normal” elements and the key `tail` contains the tail of the list. May be created using `new Prolog.List(array, tail)`.

**Object of class `Object`**  
Plain JavaScript objects are translated into a Prolog `dict`. Note that JavaScript object keys are always strings and (thus) all dict keys are atoms. This, {1:"one"} is translated into `_{'1': one}`.

**ArrayBuffer**  
Instances of `ArrayBuffer` are translated into a Prolog string that consists of characters in the range `0 ... 255`.

**Objects of a one class not being `Object`**  
Instances of non-plain JavaScript objects are translated into a Prolog *blob*. Such objects are written as `<js_Class(id)>`. The Prolog interface allows for passing the objects back and calling methods on them. See [section 13.3](wasm-js-call.html#sec:13.3).

#### 13.2.3.2 Translating Prolog data to JavaScript

Most of the translation from Prolog data to JavaScript is the reverse of the translation described in [section 13.2.3.1](wasm-calling.html#sec:13.2.3.1). In some cases however reverse translation is ambiguous. For example, both `42` and `42n` (a JavaScript `BigInt`) translate to a simple Prolog integer. The other way around, as JavaScript `Number` is a float, both Prolog `42` and `42.0` translate to `42` in JavaScript.

**Variable**  
Translate to a JavaScript `Prolog.Variable` instance where the identifier is a unique number of each unique variable.

**Integer**  
Translate to a JavaScript `Number` when possible or `BigInt` otherwise. Currently JavaScript `Number` can represent integers up to `2^53` precisely.

**Rational**  
Translate to a JavaScript `Prolog.Rational` instance.

**Float**  
Translate to a JavaScript `Number`.

**Atom**  
Translate to a JavaScript `String`.

**String**  
Translate to a JavaScript `Prolog.String` instance.

**List**  
When a *proper list* create a JavaScript `Array`, otherwise create a JavaScript `Prolog.List` instance.

**Compound term**  
Create a JavaScript `Prolog.Compound` instance.

**Dict**  
Create a plain JavaScript `Object` with the same keys. If the dict has a non-var *tag*, add a `$tag` property.
