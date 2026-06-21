
# 13 Using SWI-Prolog in your browser (WASM)

The SWI-Prolog WebAssembly (WASM) port lets you run SWI-Prolog directly in your browser. This is a fairly comprehensive version of SWI-Prolog that supports the core system as well as a good selection of *packages*, including many of the *foreign packages*.

The WASM port uses [Emscripten](https://emscripten.org/) to compile the SWI-Prolog source code to WASM, which runs on a virtual machine that is provided by almost all modern browsers.

To build the SWI-Prolog WASM port, see the building instructions on [the wiki page](https://swi-prolog.discourse.group/t/swi-prolog-in-the-browser-using-wasm/5650)

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[13.1 Loading and initializing Prolog](wasm-loading.html)

[13.1.1 Loading Prolog files](wasm-loading.html#sec:13.1.1)

[13.2 Calling Prolog from JavaScript](wasm-calling.html)

[13.2.1 The JavaScript class Query](wasm-calling.html#sec:13.2.1)

[13.2.2 Using engines](wasm-calling.html#sec:13.2.2)

[13.2.3 Translating data between JavaScript and Prolog](wasm-calling.html#sec:13.2.3)

[13.2.3.1 Translating JavaScript data to Prolog](wasm-calling.html#sec:13.2.3.1)

[13.2.3.2 Translating Prolog data to JavaScript](wasm-calling.html#sec:13.2.3.2)

[13.3 Accessing JavaScript from Prolog](wasm-js-call.html)

[13.3.1 Asynchronous access to JavaScript from Prolog](wasm-js-call.html#sec:13.3.1)

[13.3.2 JavaScript Promise that can be aborted](wasm-js-call.html#sec:13.3.2)
