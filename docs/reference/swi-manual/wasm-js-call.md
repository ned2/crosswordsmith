
## 13.3 Accessing JavaScript from Prolog

This section describes how we can interact with JavaScript from Prolog. The interface is captured in a single predicate :=/2/.

`Left` **:=** `Right`  
Depending on `Left`, this predicate implements two different actions. If `Left` is a Prolog variable, it evaluates the expression `Right` in JavaScript and unifies the result to `Left`. If `Left` is a term `Obj`\[`Key`\], where `Key` is an atom, it accesses a JavaScript *setter*. The general form of an expression is `Expression`\[`Callable`\] or simply `Callable`. If `Callable` is compound it expresses a function (or method) call. Otherwise we call JavaScript **eval()**, except for these special values:

**window**  
The main browser window itself (`undefined` when not in a browser).

**prolog**  
The `Prolog` instance.

Prolog values are translated according to the rules in [section 13.2.3.2](wasm-calling.html#sec:13.2.3.2) and the result is translated back to Prolog according to the rules in [section 13.2.3.1](wasm-calling.html#sec:13.2.3.1). Because callables are translated to function calls, object properties or global variables we need an escape to pass them as data. This is achieved using the prefix operator `#`. Note that lists are passed as JavaScript arrays rather than calls to the list functor. For convenience Prolog strings are by default translated to JavaScript `String` objects rather than `Prolog.String` instances. Below are some examples:

``` code
?- Res := myfunc([1,2,3]).
?- Max := 'Math'.max(10, 20).
?- Out := document.getElementById('output').
?- Par := document.createElement(p),
   Par.textContent := #Text.
?- Par.textContent := "aap" + " " + "noot".
```

Some JavaScript expressions are not implemented as functions. The following “functions” are handled directly by the implementation.

**instanceof**  
Returns the name of the class to which the object belongs. Same as `Obj.constructor.name`.

**instanceof**(`ClassName`)  
Returns a `Boolean` indicating whether the object is an instance of `ClassName`. Note that the class name must be an atom and as JavaScript class names normally start with a capital, the names typically need to be quoted using *single* quotes. For example:

``` code
?- W := window, T := W.instanceof('Window').
W = <js_Window>(1),
T = true.
```

**`-`**(`Any`)  
Numerical negation

**`!`**(`Any`)  
Logical negation.

**`+`**(`Any, Any`)  
**`-`**(`Any, Any`)  
**`*`**(`Any, Any`)  
**`/`**(`Any, Any`)  
**`&`**(`Any, Any`)  
**`|`**(`Any, Any`)  
**`&&`**(`Any, Any`)  
**`||`**(`Any, Any`)  
Binary operators. Note that some are not defined as Prolog operators and thus one must write e.g. `A := &&(true,false)`. `||` is not a Prolog atom, so logical disjunction gets `A := '||'(false,false)`.

\[semidet\]**is_object**(`@Term`)  
True if `Term` is a reference to a JavaScript object.

\[semidet\]**is_object**(`@Term, ?Class`)  
True when `Term` is an instance of `Class`. If `Class` is unbound it is unified with the name of the *constructor*, otherwise a JavaScript `Term instanceof Class` is executed.

**js_script**(`+String, +Options`)  
Evaluate `String` as JavaScript. This is designed to cooperate with string *quasi quotations*, so we can write e.g.,

``` code
:- use_module(library(strings)).
:- js_script({|string||
function myfunc(a)
...
|}).
```

The implementation uses =:/2, calling the JavaScript function **eval()**.

**fetch**(`+URL, +Type, -Data`)  
Wrapper around JavaScript **fetch()**, conversion of the `Response` object and waiting for the `Promise`. Type is an atom or string that is used as method on the `Response` object. Examples are `text`, `json`, `html` or `blob`. The `blob` type returns the `Data` as a string of *bytes*, i.e., character codes in the range `0 ... 255`.

#### 13.3.1 Asynchronous access to JavaScript from Prolog

While [section 13.3](wasm-js-call.html#sec:13.3) describes synchronous calls from Prolog to JavaScript, we also need asynchronous calling to implement [sleep/1](miscpreds.html#sleep/1), wait for user input, downloading documents from the web, etc. Asynchronous calling is achieved by *yielding* from the Prolog virtual machine. This can only be done when Prolog is told to expect that the VM may yield. This is implemented by `Prolog.`[`forEach()`](wasm-calling.html#forEach()) as described in [section 13.2](wasm-calling.html#sec:13.2).

\[det\]**await**(`+Promise, -Result`)  
Yield the Prolog VM, returning control back to JavaScript. When this is called from Prolog invoked using `Prolog.`[`forEach()`](wasm-calling.html#forEach()), execution of [await/2](wasm-js-call.html#await/2) completes when the `Promise` resolves and `Result` is unified with the value passed to the `Promise.`**`then()`** method. As an exception to the normal conversion rules, if the result is a single `String`, it is returned as a Prolog string rather than an atom. When the `Promise` is rejected [await/2](wasm-js-call.html#await/2) throws an exception. Note that [await/2](wasm-js-call.html#await/2) allows, for example, downloading a URL from Prolog:

``` code
?- FP := fetch("test.pl"), await(FP, Response),
   TP := Response.text(), await(TP, T).
FP = <js_Promise>(4),
Response = <js_Response>(5),
TP = <js_Promise>(6),
T = "% :- debug(js) ...".
```

Calls to [await/2](wasm-js-call.html#await/2) may be asynchronously aborted by calling `Prolog.`**`abort()`** if `Promise` implements `.`**`abort()`**. See [section 13.3.2](wasm-js-call.html#sec:13.3.2) for implementing such a promise.

\[semidet\]**is_async**  
True when we can call [await/2](wasm-js-call.html#await/2) in the current state. This implies Prolog has been called from JavaScript code that is prepared to deal with Prolog yielding and Prolog is not inside a callback from C (WASM).

#### 13.3.2 JavaScript Promise that can be aborted

A `Promise` resolves or is rejected. As Prolog waits for a specific promise on a call to [await/2](wasm-js-call.html#await/2) we may want to abort long running operations. This may be achieved using the class `Prolog.Promise` which extends `Promise`. To make the promise abortable the *executor* function must have an `abort` property. Below is the code for `Prolog.`**`promise_sleep()`** that implements this schema. First we create the *executor* and use properties on the function itself to represent the necessary state information (here, the running timer). Next, we add an `abort` property the clears the timer and runs the `reject` callback of the `Promise`. Finally we return an instance of `Prolog.Promise` which implements `.`**`abort()`**.

``` code
promise_sleep(time)
{ const f = function(resolve, reject)
  { f.reject = reject;
    f.timer = setTimeout(() =>
      { f.timer = undefined;
        resolve(true);
      }, time*1000);
  };

  f.abort = function()
  { if ( f.timer )
    { clearTimeout(f.timer);
      f.timer = undefined;
      f.reject("abort");
    }
  }

  return new Prolog.Promise(f);
}
```
