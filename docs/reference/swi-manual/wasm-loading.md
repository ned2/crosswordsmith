
## 13.1 Loading and initializing Prolog

The WASM SWI-Prolog distribution consists of three files:

**`swipl-web.js`**  
This is the main file that must be loaded using a `<script>` element. It defines a global function `SWIPL` that loads the other components and gives access to the SWI-Prolog system.

**`swipl-web.wasm`**  
This is the actual WASM binary containing the compiled C code of the core, foreign packages and required libraries.

**`swipl-web.data`**  
This data contains a *file system* that is *mounted* on `/swipl`. It contains the Prolog startup code and libraries.

Below is the skeleton for getting access to the Prolog system. We first define the global `Prolog` and `Module` objects. The `options` object provides module configuration options. **SWIPL()** returns a `Promise` that resolves when the WASM system is loaded and initialized. It returns the WASM *module*, which contains `Module.prolog` to an instance of the class `Prolog` that provides a high level interface from JavaScript.

``` code
let Prolog;
let Module;
const options = {
  // Provide options for customization
};

SWIPL(options).then((module) =>
{ Module = module;
  Prolog = Module.prolog;

  // Start using Prolog
};
```

The `options` object defines customization properties for the Emscripten module as well as for Prolog. We highlight the important properties below.

**arguments**  
An `Array` of `String` objects that define the commandline arguments for initializing Prolog. `argv[0]` is *not* part of this array. Few arguments are useful in this context. The **-q** may be used to suppress the Prolog banner.

**locateFile**  
A `Function` that is used to translate `swipl-web.wasm` and `swipl-web.data` into a (relative) URL. Default is to find these resources in the same directory of the server. For example, to load `swipl-web.wasm` and `swipl-web.data` from `/wasm` use

``` code
var Module = {
  ...,
  locateFile: (file) => '/wasm/' + file
}
```

**on_output**  
A `Function` that is called when Prolog writes to `user_output` or `user_error`. It is passed two arguments: a `String` containing the text to emit and one of the constant strings `"stdout"` or `"stderr"` to indicate the output stream. This uses the Emscripten `Module.FS.init` option to rebind the output and error streams, providing behaviour that is similar to the Emscripten properties `print` and `printErr`. However, our passed string contains the newline character and the handler is called when Prolog *flushes* the output. Normally the callback should insert a `<span>` element that has (at least) the following style:

``` code
.stderr, .stdout {
  white-space: pre-wrap;
  font-family: monospace;
  overflow-wrap: anywhere;
}
```

Here is a simple implementation of **print()**, assuming the document has a `<div id="output">` node.

``` code
function print(line, cls)
{ const output = document.getElementById("output");
  const node   = document.createElement('span');

  node.className = cls;
  node.textContent = line;
  output.appendChild(node);
};
```

### 13.1.1 Loading Prolog files

The WASM build ships with the Prolog library and thus Prolog libraries can be loaded as normal using [use_module/1](import.html#use_module/1), etc., for example, we can include the `lists` library using this directive. Note that the normal *autoloading* of library code works in the WASM version.

``` code
:- use_module(library(lists)).
```

When Prolog is in *asynchronous* mode, i.e., called through [Prolog.forEach()](wasm-calling.html#Prolog.forEach()), we can also load code from a *URL*. For example, we can load the `CHAT80` demo program directly from [GitHub](https://github.com/JanWielemaker/chat80) using^(246The `\`c continues the quoted atom from the next line after removing leading white space.)

``` code
?- consult('https://raw.githubusercontent.com/JanWielemaker/\c
            chat80/master/prolog/chat80.pl').
```

Larger files can be loaded as `.qlf` files. See [section 4.3.3](consulting.html#sec:4.3.3) and [qcompile/2](consulting.html#qcompile/2). Notably we can create a single qlf file from an application using the `include(user)` option. Below we create a `.qlf` file from [CHAT80](https://github.com/JanWielemaker/chat80). The resulting `chat80.qlf` can be loaded from a URL using [consult/1](consulting.html#consult/1) as above.

``` code
?- qcompile('chat80.pl', [include(user)]).
```

There are three ways to load Prolog code from JavaScript: (1) loading from a string, (2) loading from `<script>` elements and (3) loading from URL. Note that all the loading methods return a `Promise` that is resolved when loading the data is completed.

`Promise` **Prolog.load_string**(`String, Id`)  
Load Prolog code from `String`, pretending it was loaded from the file `Id`. The `Id` is optional. When omitted it generates `/string/1`, `/string/2`, ... .

`Promise` **Prolog.load_scripts**()  
Load all scripts from the current document that have their `type` set to `text/prolog`. The file reference for the loaded script is `/script/Id`, where `Id` is derived from (1) the `id` of the script, (2) the `name` of the script or (3) being the nth Prolog script in the document. When resolved, the promise returns an array with the names of the loaded scripts.

`Promise` **Prolog.consult**(`...Sources[, Options]`)  
Load the given `Sources`. Each source is either a file from the local file system, e.g., `library(lists)` or a URL. The sources are downloaded and processed sequentially. This uses [Prolog.forEach()](wasm-calling.html#Prolog.forEach()) calling [load_files/1](consulting.html#load_files/1). The returned `Promise` returns 1 on success. If the last argument is an object, it is treated as options. Options processed:

**module**: `string`  
Defines the module into which a non-module file is loaded or into which the exports of a module file are imported. Default is `"user"`.

**engine**: `Boolean`  
If `true`, run the compilation in a temporary engine. This keeps the main Prolog engine available for other tasks while the file is being loaded.
