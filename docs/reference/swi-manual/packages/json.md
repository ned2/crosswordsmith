SWI-Prolog JSON library

Jan Wielemaker  
SWI-Prolog Solutions b.v.  
The Netherlands  
E-mail: [jan@swi-prolog.org](mailto:jan@swi-prolog.org)

Abstract

This package reads and writes JSON documents from and to SWI-Prolog streams, files and strings.

# Table of Contents

[1 Supporting JSON](#sec:1)

[1.1 library(json): Reading and writing JSON serialization](#sec:1.1)

[1.2 library(json_schema): JSON Schema reader and validator](#sec:1.2)

[1.2.1 Status](#sec:1.2.1)

[1.2.2 Predicates](#sec:1.2.2)

[1.3 library(json_rpc_client): JSON RPC client](#sec:1.3)

[1.4 library(json_rpc_server): JSON RPC Server](#sec:1.4)

[1.5 library(json_convert): Convert between JSON terms and Prolog application terms](#sec:1.5)

[1.6 library(http/http_json): HTTP JSON Plugin module](#sec:1.6)

## 1 Supporting JSON

From [http://json.org](http://json.org), " JSON (JavaScript Object Notation) is a lightweight data-interchange format. It is easy for humans to read and write. It is easy for machines to parse and generate. It is based on a subset of the JavaScript Programming Language, Standard ECMA-262 3rd Edition - December 1999. JSON is a text format that is completely language independent but uses conventions that are familiar to programmers of the C-family of languages, including C, C++, C#, Java, JavaScript, Perl, Python, and many others. These properties make JSON an ideal data-interchange language."

Although JSON is nowadays used a lot outside the context of web applications, SWI-Prolog's support for JSON started life as part of the HTTP package. SWI-Prolog supports two Prolog representations for JSON terms. The first and oldest map JSON objects to a term `json(PropertyList)` and use the `@` functor to disambiguate e.g. `null` from the string `"null"`, leading to `@(null)`. As of SWI-Prolog version 7, JSON objects may be represented using *dict* objects and JSON strings using Prolog strings. Predicates following this convention are suffixed with `_dict`, e.g. [json_read_dict/2](#json_read_dict/2). For example, given the JSON document

``` code
{ "name": "Bob", "children": ["Mary", "John"], "age":42, "married": true }
```

we get either (using [json_read/2](#json_read/2)):

``` code
json([name='Bob', children=['Mary', 'John'], age=42, married= @(true)]).
```

or (using [json_read_dict/2](#json_read_dict/2)):

``` code
#{age:42, children:["Mary", "John"], married:true, name:"Bob"}
```

The SWI-Prolog JSON interface consists of the following libraries:

- `library(json)` provides support for the core JSON object serialization and parsing.

- `library(json_schema)` implements JSON schema validation.

- `library(json_convert)` converts between the primary representation of JSON terms in Prolog and more application oriented Prolog terms. E.g. `point(X,Y)` vs. `object([x=X,y=Y])`.

- `library(http/http_json)` hooks the conversion libraries into the HTTP client and server libraries.

  ### 1.1 library(json): Reading and writing JSON serialization

  See also  
  \- `http_json.pl` links JSON to the HTTP client and server modules.  
  - `json_convert.pl` converts JSON Prolog terms to more comfortable terms.

  This module supports reading and writing JSON objects. This library supports two Prolog representations (the *new* representation is only supported in SWI-Prolog version 7 and later):

  - The **classical** representation is provided by [json_read/3](#json_read/3) and [json_write/3](#json_write/3). This represents a JSON object as `json(NameValueList)`, a JSON string as an atom and the JSON constants `null`, `true` and `false` as @(null), @(true) and @false.
  - The **new** representation is provided by [json_read_dict/3](#json_read_dict/3) and [json_write_dict/3](#json_write_dict/3). This represents a JSON object as a dict, a JSON string as a Prolog string and the JSON constants using the Prolog atoms `null`, `true` and `false`.

  This module provides the `json` *Quasi Quotation* syntax that allows for embedding JSON documents in Prolog.

  \[det\]**atom_json_term**(`?Atom, ?JSONTerm, +Options`)  
  Convert between textual representation and a JSON term. In *write* mode (`JSONTerm` to `Atom`), the option

  **as**(`Type`)  
  defines the output type, which is one of `atom` (default), `string`, `codes` or `chars`.

  \[det\]**json_read**(`+Stream, -Term`)  
  \[det\]**json_read**(`+Stream, -Term, +Options`)  
  Read next JSON value from `Stream` into a Prolog term. The canonical representation for `Term` is:

  - A JSON object is mapped to a term `json(NameValueList)`, where NameValueList is a list of Name=Value. Name is an atom created from the JSON string.
  - A JSON array is mapped to a Prolog list of JSON values.
  - A JSON string is mapped to a Prolog atom
  - A JSON number is mapped to a Prolog number
  - The JSON constants `true` and `false` are mapped -like JPL- to @(true) and @(false).
  - The JSON constant `null` is mapped to the Prolog term @(null)

  Here is a complete example in JSON and its corresponding Prolog term.

  ``` code
  { "name":"Demo term",
    "created": {
      "day":null,
      "month":"December",
      "year":2007
    },
    "confirmed":true,
    "members":[1,2,3]
  }
  ```

  ``` code
  json([ name='Demo term',
         created=json([day= @null, month='December', year=2007]),
         confirmed= @true,
         members=[1, 2, 3]
       ])
  ```

  The following options are processed:

  **null**(`+NullTerm`)  
  `Term` used to represent JSON `null`. Default @(null)

  **true**(`+TrueTerm`)  
  `Term` used to represent JSON `true`. Default @(true)

  **false**(`+FalseTerm`)  
  `Term` used to represent JSON `false`. Default @(false)

  **end_of_file**(`+ErrorOrTerm`)  
  If end of file is reached after skipping white space but before any input is processed take the following action (default `error`):

  - If `ErrorOrTerm` `==` `error`, throw an unexpected end of file syntax error
  - Otherwise return `ErrorOrTerm`.

  Returning an status term is required to process [Concatenated JSON](https://en.wikipedia.org/wiki/JSON_streaming\#Concatenated_JSON). Suggested values are `@(eof)` or `end_of_file`.

  **value_string_as**(`+Type`)  
  Prolog type used for strings used as value. Default is `atom`. The alternative is `string`, producing a packed string object. Please note that `codes` or `chars` would produce ambiguous output and are therefore not supported.

  See also  
  [json_read_dict/3](#json_read_dict/3) to read a JSON term using the version 7 extended data types.

  \[det\]**json_write**(`+Stream, +Term`)  
  \[det\]**json_write**(`+Stream, +Term, +Options`)  
  Write a JSON term to `Stream`. The JSON object is of the same format as produced by [json_read/2](#json_read/2), though we allow for some more flexibility with regard to pairs in objects. All of Name=Value, Name-Value and Name(Value) produce the same output.

  Values can be of the form \#(`Term`), which causes `Term` to be *stringified* if it is not an atom or string. Stringification is based on term_string/2.

  Rational numbers are emitted as floating point numbers. The hook [json_write_hook/4](#json_write_hook/4) can be used to realize domain specific alternatives.

  The version 7 *dict* type is supported as well. Optionally, if the dict has a *tag*, a property "type":"tag" can be added to the object. This behaviour can be controlled using the `tag` option (see below). For example:

  ``` code
  ?- json_write(current_output, point{x:1,y:2}).
  {
    "x":1,
    "y":2
  }
  ```

  ``` code
  ?- json_write(current_output, point{x:1,y:2}, [tag(type)]).
  {
    "type":"point",
    "x":1,
    "y":2
  }
  ```

  In addition to the options recognised by [json_read/3](#json_read/3), we process the following options are recognised:

  **width**(`+Width`)  
  `Width` in which we try to format the result. Too long lines switch from *horizontal* to *vertical* layout for better readability. If performance is critical and human readability is not an issue use `Width` = 0, which causes a single-line output.

  **step**(`+Step`)  
  Indentation increment for next level. Default is 2.

  **tab**(`+TabDistance`)  
  Distance between tab-stops. If equal to Step, layout is generated with one tab per level.

  **serialize_unknown**(`+Boolean`)  
  If `true` (default `false`), serialize unknown terms and print them as a JSON string. The default raises a type error. Note that this option only makes sense if you can guarantee that the passed value is not an otherwise valid Prolog representation of a Prolog term.

  If a string is emitted, the sequence `</` is emitted as `<\/`. This is valid JSON syntax which ensures that JSON objects can be safely embedded into an HTML `<script>` element.

  \[semidet,multifile\]**json_write_hook**(`+Term, +Stream, +State, +Options`)  
  Hook that can be used to emit a JSON representation for `Term` to `Stream`. If the predicate succeeds it **must** have written a **valid** JSON data element and if it fails it may not have produced any output. This facility may be used to map arbitrary Prolog terms to JSON. It was added to manage the precision with which floating point numbers are emitted.

  Note that this hook is shared by all users of this library. It is generally advised to map a unique compound term to avoid interference with normal output.

  |  |  |
  |----|----|
  | `State` | and `Options` are opaque handles to the current output state and settings. Future versions may provide documented access to these terms. Currently it is advised to ignore these arguments. |

  \[semidet,multifile\]**json_dict_pairs**(`+Dict, -Pairs`)  
  This hook may be used to order the keys of an object. If it fails, dict_pairs/3 is used which produces an ordered list of keys.

  \[semidet\]**is_json_term**(`@Term`)  
  \[semidet\]**is_json_term**(`@Term, +Options`)  
  True if `Term` is a json term. `Options` are the same as for [json_read/2](#json_read/2), defining the Prolog representation for the JSON `true`, `false` and `null` constants.

  \[det\]**json_read_dict**(`+Stream, -Dict`)  
  \[det\]**json_read_dict**(`+Stream, -Dict, +Options`)  
  Read a JSON object, returning objects as a dicts. The representation depends on the options, where the default is:

  - String values are mapped to Prolog strings
  - JSON `true`, `false` and `null` are represented using these Prolog atoms.
  - JSON objects are mapped to dicts.
  - Optionally, a `type` field in an object assigns a tag for the dict.

  The predicate [json_read_dict/3](#json_read_dict/3) processes the same options as [json_read/3](#json_read/3), but with different defaults. In addition, it processes the `tag` option. See [json_read/3](#json_read/3) for details about the shared options.

  **tag**(`+Name`)  
  When converting to/from a dict, map the indicated JSON attribute to the dict *tag*. No mapping is performed if `Name` is the empty atom (” , default). See [json_read_dict/2](#json_read_dict/2) and [json_write_dict/2](#json_write_dict/2).

  **default_tag**(`+Tag`)  
  Provide the default tag if the above `tag` option does not apply.

  **null**(`+NullTerm`)  
  Default the atom `null`.

  **true**(`+TrueTerm`)  
  Default the atom `true`.

  **false**(`+FalseTerm`)  
  Default the atom `false`

  **end_of_file**(`+ErrorOrTerm`)  
  Action on reading end-of-file. See [json_read/3](#json_read/3) for details.

  **value_string_as**(`+Type`)  
  Prolog type used for strings used as value. Default is `string`. The alternative is `atom`, producing a packed string object.

  \[det\]**json_write_dict**(`+Stream, +Dict`)  
  \[det\]**json_write_dict**(`+Stream, +Dict, +Options`)  
  Write a JSON term, represented using dicts. This is the same as [json_write/3](#json_write/3), but assuming the default representation of JSON objects as dicts.

  \[det\]**atom_json_dict**(`+Atom, -JSONDict, +Options`)  
  \[det\]**atom_json_dict**(`-Text, +JSONDict, +Options`)  
  Convert between textual representation and a JSON term represented as a dict. `Options` are as for [json_read/3](#json_read/3). In *write* mode, the additional option

  **as**(`Type`)  
  defines the output type, which is one of `atom`, `string` or `codes`.

  \[det\]**json**(`+Content, +Vars, +VarDict, -JSON`)  
  The predicate [json/4](#json/4) implements `JSON` quasi quotations. These quotations produce a `JSON` dict that is suitable for [json_write_dict/2](#json_write_dict/2). The quasi quoter only accepts valid, but possibly partial `JSON` documents. The quoter replaces content whose value is a Prolog variable that appears in the argument list of the `json` indicator. Notably, you can't use a Prolog variable in place of an object key. Here is an example.

  ``` code
    {|json(Name)||
        { "name": Name,
          "created": {
            "day":null,
            "month":"December",
            "year":2007
          },
          "confirmed":true,
          "members":[1,2,3]
        }
    |}.
  ```

  ### 1.2 library(json_schema): JSON Schema reader and validator

  This module provides a JSON Schema reader and validator. This module is based on the 2020-12 draft of the specification.

  The API consists of two primitives and a simple high level predicate ([json_validate/3](#json_validate/3)):

  - [json_compile_schema/3](#json_compile_schema/3) translates a file or parsed JSON schema data into a Prolog term that represents the schema in a way that is more suitable for checking a JSON document.
  - [json_check/3](#json_check/3) validates a document against the compiled schema.

  #### 1.2.1 Status

  The implementation is validated against [https://github.com/json-schema-org/JSON-Schema-Test-Suite.git](https://github.com/json-schema-org/JSON-Schema-Test-Suite.git). It fails 4 out of 1,261 tests. Issues:

  - `$dynamicRef` fails 3 tests, probably mixing up physical context and resolution context that defines the *dynamic scopes*.
  - There is no support for non-default vocabulary selection.

  The current implementation is the result of an incremental process. It should be refactored to make it a cleaner translation of the specification.

  #### 1.2.2 Predicates

  \[det\]**json_validate**(`+SchemaFile, +DataDict, +Options`)  
  Given a file holding a JSON Schema and a Prolog dict holding JSON data, validate the data against the schema. `Options` are passed to [json_compile_schema/3](#json_compile_schema/3) and [json_check/3](#json_check/3).

  throws  
  `error(Formal, json_path(Path))`, where `Path` is a list of properties from the root element to the culprit element. Formal is typically a type, domain or existence error. This file contains the message hooks to generate a human readable error from these exceptions using print_message/2.

  \[det\]**json_compile_schema**(`+Input, -Type, +Options`)  
  Load and translated a JSON Schema. `Input` is either a file name, a specification for absolute_file_name/3 or the output of [json_read_dict/2](#json_read_dict/2).

  If `Input` is a file name, the loaded and compiled schema is cached. Reusing the cache validates the modification file of the schema file and reloads it if the file's time stamp has changed. Note that `true` and `false` are valid schemas and cannot be used as file names.

  \[semidet\]**json_check**(`+Spec, ?JSON, +Options`)  
  Validate a `JSON` object. `Spec` is a Prolog representation of the schema that is optimized for validation. This representation is derived from `JSON` data using [json_compile_schema/3](#json_compile_schema/3). `Options`:

  **on_error**(`Mode`)  
  What to do if an error is found. Defined modes are

  **error**  
  Raise an exception. This is the default. Note that only the first error is reported this way.

  **warning**  
  Print a message

  **silent**  
  Fail

  **value_string_as**(`Type`)  
  Same as for [json_read/3](#json_read/3).

  This predicate is often used through validate_json_dict/3, which mantains a cached mapping from the `JSON` Schema to `Spec`.

  ### 1.3 library(json_rpc_client): JSON RPC client

  This module implements a JSON RPC compliant client. The three predicates require a *stream pair* (see stream_pair/2) that connects us to a JSON RPC server.

  \[det\]**json_call**(`+Stream, +Goal, -Result, +Options`)  
  Run `Goal` on a JSON RPC service identified by `Stream` and wait for `Result`. This predicate may be called from multiple threads. As replies come in in arbitrary order, this predicate starts a thread the reads the replies from `Stream` and informs the calling thread using a Prolog message queue.

  If `Stream` is closed this library terminates the thread and related message queue.

  `Options` are passed to [json_write_dict/3](#json_write_dict/3) and thread_get_message/3. Additional options:

  **async**(`:Closure`)  
  Do not wait for the request to complete. Instead, call `call(Closure, Data)` from the client reading thread when the request is completed. If `Closure` is `true`, ignore the reply. As we cannot inject errors as exceptions in the calling thread, possible errors are printed.

  **thread_alias**(`+Atom`)  
  Alias name to use for the thread that deals with incomming replies and requests. Defaults to `json_rpc_client:<N>`, where *N* is a unique number.

  [TABLE]

  \[det\]**json_notify**(`+Stream, +Goal, +Options`)  
  Run `Goal` on a JSON RPC service identified by `Stream` without waiting for the result.

  \[det\]**json_batch**(`+Stream, +Notifications:list, +Calls:list, -Results:list, +Options`)  
  Run a batch of notifications and normal calls on the JSON server at the other end of `Stream`. On success, Result is unified to a list with the same length as `Calls`. Each element either contains a value, similar to [json_call/4](#json_call/4) or a term `error(Dict)`, where `Dict` holds the `code`, `message` and optional `data` field. Note that `error(Dict)` is not a valid JSON type and this is thus unambiguous. While the JSON RPC standard allows the server to process the messages in any order and allows for concurrent processing, all results are sent in one message and this client ensures the elements of the `Results` list are in the same order as the `Calls` list. If the `Calls` list is empty this predicate does not wait for a reply.

  \[det\]**json_full_duplex**(`+Stream, :Options`)  
  Start the thread for incomming data and on requests, dispatch them using `library(jso_rpc_server)` in the module derived from the `Options` list.

  ### 1.4 library(json_rpc_server): JSON RPC Server

  See also  
  [JSON-RPC](https://www.jsonrpc.org/specification)

  This module implements an JSON RPC server. It provides declarations that bind Prolog predicates to JSON RPC methods and a dispatch loop that acts on a bi-directional stream. This module assumes a two-directional stream and provides [json_rpc_dispatch/2](#json_rpc_dispatch/2) that receiveds JSON messages on the input side of this stream and sends the replies through the output. This module does not implement obtaining such a stream. Obvious candidates for obtaining a stream are:

  - Using standard I/O to a child process. See process_create/3.
  - Using sockets. See `library(socket)`. Using the SSL package this also provides secure sockets.
  - Using the HTTP package to extablish a *web socket*.

  This library defines [json_method/1](#json_method/1) for declaring predicates to act as a JSON method. The declaration accepts a JSON Schema specification, represented as a SWI-Prolog dict to specify the input parameters as well as the output.

  **json_method**(`+Methods`)  
  `Methods` is a comma-list of JSON RPC method declarations. Each declaration takes one of the forms below:

  `Callable`**`:`**`Reply`  
  Here, `Callable` is a Prolog callable term whose name and number of argument match a predicate in this module. The arguments are JSON Schema types and `Reply` is a JSON Schema type.

  **`Callable`**  
  `Callable` is as above, but there is no return value. This implements JSON RPC *notifications*, i.e., asynchronously processed messages for which we do not wait for a reply.

  For example:

  ``` code
  :- json_method
      subtract(#{type:number}, #{type:number}): #{type:number}.

  subtract(A, B, R) :- R is A-B.
  ```

  `Methods` with *named arguments* can be implemented using a single argument that is an object with specified properties. For example, the program below implements a depositing to a bank account. The method takes an `account` and `amount` parameter and returns the new balance. The [json_rpc_error/2](#json_rpc_error/2) throws a JSON RPC *application error*.

  ``` code
  :- json_method
      deposit(#{ properties:
                 #{ account: #{type:string},
                    amount:  #{type:number}
                  }}): #{type:number},

  deposit(Request, Reply),
      #{account: Account, amount: Amount} :< Request =>
      transaction((   retract(account(Account, Old))
                  ->  New is Old+Amount,
                      asserta(account(Account, New))
                  ;   json_rpc_error(2, "Account does not exist")
                  )),
      Reply = New.
  ```

  \[det\]**json_rpc_dispatch**(`:Stream, +Options`)  
  Run the JSON RPC dispatch loop until end of file is reached on `Stream`.

  |  |  |
  |----|----|
  | `Stream` | is stream pair (see stream_pair/2). Normally, the stream should use `utf8` encoding. If the stream is a binary stream, it will be processed as if `utf8` encoding is enabled. If it is a text stream the encoding of the stream is respected. |

  **json_rpc_dispatch_request**(`+Module, +Stream, +Request, +Options`)  
  Handle a request that has been read from `Stream`, possibly sending a reply to `Stream`.

  **json_rpc_error**(`+Code, +Message`)  
  **json_rpc_error**(`+Code, +Message, +Data`)  
  Normally called from a method implementation to raise an *application error*.

  |  |  |
  |----|----|
  | `Code` | is an integer. The range -32768 to -32000 is reserved for JSON RPC server errors. |
  | `Message` | is a short string decribing the error |
  | `Data` | is optional JSON data that provides context for the error. |

  Errors  
  `json_rpc_error(Dict)`, where `Dict` contains the JSON RPC defined fields `code`, `message` and optionally `data`.

  ### 1.5 library(json_convert): Convert between JSON terms and Prolog application terms

  To be done  
  \- Ignore extra fields. Using a partial list of *extra*?  
  - Consider a sensible default for handling JSON `null`. Conversion to Prolog could translate @null into a variable if the desired type is not `any`. Conversion to JSON could map variables to `null`, though this may be unsafe. If the Prolog term is known to be non-ground and JSON @null is a sensible mapping, we can also use this simple snippet to deal with that fact.

  ``` code
          term_variables(Term, Vars),
          maplist(=(@null), Vars).
  ```

  The idea behind this module is to provide a flexible high-level mapping between Prolog terms as you would like to see them in your application and the standard representation of a JSON object as a Prolog term. For example, an X-Y point may be represented in JSON as `{"x":25, "y":50}`. Represented in Prolog this becomes `json([x=25,y=50])`, but this is a pretty non-natural representation from the Prolog point of view.

  This module allows for defining records (just like `library(record)`) that provide transparent two-way transformation between the two representations.

  ``` code
  :- json_object
          point(x:integer, y:integer).
  ```

  This declaration causes [prolog_to_json/2](#prolog_to_json/2) to translate the native Prolog representation into a JSON Term:

  ``` code
  ?- prolog_to_json(point(25,50), X).

  X = json([x=25, y=50])
  ```

  A [json_object/1](#json_object/1) declaration can define multiple objects separated by a comma (,), similar to the dynamic/1 directive. Optionally, a declaration can be qualified using a module. The conversion predicates [prolog_to_json/2](#prolog_to_json/2) and [json_to_prolog/2](#json_to_prolog/2) first try a conversion associated with the calling module. If not successful, they try conversions associated with the module `user`.

  JSON objects have no *type*. This can be solved by adding an extra field to the JSON object, e.g. `{"type":"point", "x":25, "y":50}`. As Prolog records are typed by their functor we need some notation to handle this gracefully. This is achieved by adding +Fields to the declaration. I.e.

  ``` code
  :- json_object
          point(x:integer, y:integer) + [type=point].
  ```

  Using this declaration, the conversion becomes:

  ``` code
  ?- prolog_to_json(point(25,50), X).

  X = json([x=25, y=50, type=point])
  ```

  The predicate [json_to_prolog/2](#json_to_prolog/2) is often used after [http_read_json/2](#http_read_json/2) and [prolog_to_json/2](#prolog_to_json/2) before [reply_json/1](#reply_json/1). For now we consider them separate predicates because the transformation may be too general, too slow or not needed for dedicated applications. Using a separate step also simplifies debugging this rather complicated process.

  \[multifile\]**current_json_object**(`Term, Module, Fields`)  
  Multifile predicate computed from the [json_object/1](#json_object/1) declarations. `Term` is the most general Prolog term representing the object. `Module` is the module in which the object is defined and `Fields` is a list of `f(Name, Type, Default, Var)`, ordered by Name. Var is the corresponding variable in `Term`.

  **json_object**(`+Declaration`)  
  Declare a JSON object. The declaration takes the same format as using in record/1 from `library(record)`. E.g.

  ``` code
  ?- json_object
        point(x:int, y:int, z:int=0).
  ```

  The type arguments are either types as know to `library(error)` or functor names of other JSON objects. The constant `any` indicates an untyped argument. If this is a JSON term, it becomes subject to [json_to_prolog/2](#json_to_prolog/2). I.e., using the type `list(any)` causes the conversion to be executed on each element of the list.

  If a field has a default, the default is used if the field is not specified in the JSON object. Extending the record type definition, types can be of the form (Type1`|`Type2). The type `null` means that the field may *not* be present.

  Conversion of JSON to Prolog applies if all non-defaulted arguments can be found in the JSON object. If multiple rules match, the term with the highest arity gets preference.

  \[semidet\]**prolog_bool_to_json**(`+Prolog, -JSON`)  
  `JSON` is the `JSON` boolean for `Prolog`. It is a flexible the `Prolog` notation for truth-value, accepting one of `true`, `on` or `1` for @true and one of `false`, `fail`, `off` or `0` for @false.

  Errors  
  instantiation_error if `Prolog` is unbound.

  \[det\]**prolog_to_json**(`:Term, -JSONObject`)  
  Translate a Prolog application `Term` into a JSON object term. This transformation is based on `:-` [json_object/1](#json_object/1) declarations. If a [json_object/1](#json_object/1) declaration declares a field of type `boolean`, commonly used truth-values in Prolog are converted to JSON booleans. Boolean translation accepts one of `true`, `on`, `1`, @true, `false`, `fail`, `off` or `0`, @false.

  Errors  
  \- `type_error(json_term, X)`  
  - instantiation_error

  \[det\]**json_to_prolog**(`+JSON, -Term`)  
  Translate a `JSON` term into an application term. This transformation is based on `:-` [json_object/1](#json_object/1) declarations. An efficient transformation is non-trivial, but we rely on the assumption that, although the order of fields in `JSON` terms is irrelevant and can therefore vary a lot, practical applications will normally generate the `JSON` objects in a consistent order.

  If a field in a json_object is declared of type `boolean`, @true and @false are translated to `true` or `false`, the most commonly used Prolog representation for truth-values.

  ### 1.6 library(http/http_json): HTTP JSON Plugin module

  See also  
  \- JSON Requests are discussed in [http://json.org/JSONRequest.html](http://json.org/JSONRequest.html)  
  - `json.pl` describes how JSON objects are represented in Prolog terms.  
  - `json_convert.pl` converts between more natural Prolog terms and json terms.

  Most code doesn't need to use this directly; instead use `library(http/http_server)`, which combines this library with the typical HTTP libraries that most servers need.

  This module adds hooks to several parts of the HTTP libraries, making them JSON-aware. Notably:

  - Make http_read_data/3 convert `application/json` and `application/jsonrequest` content to a JSON term.
  - Cause http_open/3 to accept `post(json(Term))` to issue a POST request with JSON content.
  - Provide HTTP server and client utility predicates for reading and replying JSON:
    - [http_read_json/2](#http_read_json/2)
    - [http_read_json/3](#http_read_json/3)
    - [http_read_json_dict/2](#http_read_json_dict/2)
    - [http_read_json_dict/3](#http_read_json_dict/3)
    - [reply_json/1](#reply_json/1)
    - [reply_json/2](#reply_json/2)
    - [reply_json_dict/1](#reply_json_dict/1)
    - [reply_json_dict/2](#reply_json_dict/2)
  - Reply to exceptions in the server using an JSON document rather then HTML if the `Accept` header prefers application/json over text/html.

  Typically JSON is used by Prolog HTTP servers. This module supports two JSON representations: the classical representation and the new representation supported by the SWI-Prolog version 7 extended data types. Below is a skeleton for handling a JSON request, answering in JSON using the classical interface.

  ``` code
  handle(Request) :-
        http_read_json(Request, JSONIn),
        json_to_prolog(JSONIn, PrologIn),
        <compute>(PrologIn, PrologOut),         % application body
        prolog_to_json(PrologOut, JSONOut),
        reply_json(JSONOut).
  ```

  When using dicts, the conversion step is generally not needed and the code becomes:

  ``` code
  handle(Request) :-
        http_read_json_dict(Request, DictIn),
        <compute>(DictIn, DictOut),
        reply_json(DictOut).
  ```

  This module also integrates JSON support into the http client provided by `http_client.pl`. Posting a JSON query and processing the JSON reply (or any other reply understood by http_read_data/3) is as simple as below, where Term is a JSON term as described in `json.pl` and reply is of the same format if the server replies with JSON.

  ``` code
        ...,
        http_post(URL, json(Term), Reply, [])
  ```

  \[multifile\]http_client:**http_convert_data**(`+In, +Fields, -Data, +Options`)  
  Hook implementation that supports reading JSON documents. It processes the following option:

  **json_object**(`+As`)  
  Where `As` is one of `term` or `dict`. If the value is `dict`, [json_read_dict/3](#json_read_dict/3) is used.

  \[semidet\]**is_json_content_type**(`+ContentType`)  
  True if `ContentType` is a header value (either parsed or as atom/string) that denotes a JSON value.

  \[semidet,multifile\]**json_type**(`?MediaType`)  
  True if `MediaType` is a JSON media type. http_json:json_type/1 is a multifile predicate and may be extended to facilitate non-conforming clients.

  |  |  |
  |----|----|
  | `MediaType` | is a term `Type`/`SubType`, where both `Type` and `SubType` are atoms. |

  \[semidet,multifile\]http:**post_data_hook**(`+Data, +Out:stream, +HdrExtra`)  
  Hook implementation that allows http_post_data/3 posting JSON objects using one of the forms below.

  ``` code
  http_post(URL, json(Term), Reply, Options)
  http_post(URL, json(Term, Options), Reply, Options)
  ```

  If Options are passed, these are handed to [json_write/3](#json_write/3). In addition, this option is processed:

  **json_object**(`As`)  
  If `As` is `dict`, [json_write_dict/3](#json_write_dict/3) is used to write the output. This is default if `json(Dict)` is passed.

  To be done  
  avoid creation of intermediate data using chunked output.

  \[det\]**http_read_json**(`+Request, -JSON`)  
  \[det\]**http_read_json**(`+Request, -JSON, +Options`)  
  Extract `JSON` data posted to this HTTP request. `Options` are passed to [json_read/3](#json_read/3). In addition, this option is processed:

  **json_object**(`+As`)  
  One of `term` (default) to generate a classical Prolog term or `dict` to exploit the SWI-Prolog version 7 data type extensions. See [json_read_dict/3](#json_read_dict/3).

  Errors  
  \- `domain_error(mimetype, Found)` if the mimetype is not known (see [json_type/1](#json_type/1)).  
  - `domain_error(method, Method)` if the request method is not a `POST`, `PUT` or `PATCH`.

  \[det\]**http_read_json_dict**(`+Request, -Dict`)  
  \[det\]**http_read_json_dict**(`+Request, -Dict, +Options`)  
  Similar to [http_read_json/2](#http_read_json/2),3, but by default uses the version 7 extended datatypes.

  \[det\]**reply_json**(`+JSONTerm`)  
  \[det\]**reply_json**(`+JSONTerm, +Options`)  
  Formulate a JSON HTTP reply. See [json_write/2](#json_write/2) for details. The processed options are listed below. Remaining options are forwarded to [json_write/3](#json_write/3).

  **content_type**(`+Type`)  
  The default `Content-type` is `application/json; charset=UTF8`. `charset=UTF8` should not be required because JSON is defined to be UTF-8 encoded, but some clients insist on it.

  **status**(`+Code`)  
  The default status is 200. REST API functions may use other values from the 2XX range, such as 201 (created).

  **json_object**(`+As`)  
  One of `term` (classical json representation) or `dict` to use the new dict representation. If omitted and Term is a dict, `dict` is assumed. SWI-Prolog Version 7.

  \[det\]**reply_json_dict**(`+JSONTerm`)  
  \[det\]**reply_json_dict**(`+JSONTerm, +Options`)  
  As [reply_json/1](#reply_json/1) and [reply_json/2](#reply_json/2), but assumes the new dict based data representation. Note that this is the default if the outer object is a dict. This predicate is needed to serialize a list of objects correctly and provides consistency with [http_read_json_dict/2](#http_read_json_dict/2) and friends.

# Index

?  
[atom_json_dict/3](#atom_json_dict/3)  
[atom_json_term/3](#atom_json_term/3)  
[current_json_object/3](#current_json_object/3)  
[http:post_data_hook/3](#http:post_data_hook/3)  
[http_client:http_convert_data/4](#http_client:http_convert_data/4)  
[http_read_json/2](#http_read_json/2)  
[http_read_json/3](#http_read_json/3)  
[http_read_json_dict/2](#http_read_json_dict/2)  
[http_read_json_dict/3](#http_read_json_dict/3)  
[is_json_content_type/1](#is_json_content_type/1)  
[is_json_term/1](#is_json_term/1)  
[is_json_term/2](#is_json_term/2)  
[json/4](#json/4)  
[json_batch/5](#json_batch/5)  
[json_call/4](#json_call/4)  
[json_check/3](#json_check/3)  
[json_compile_schema/3](#json_compile_schema/3)  
[json_dict_pairs/2](#json_dict_pairs/2)  
[json_full_duplex/2](#json_full_duplex/2)  
[json_method/1](#json_method/1)  
[json_notify/3](#json_notify/3)  
[json_object/1](#json_object/1)  
[json_read/2](#json_read/2)  
[json_read/3](#json_read/3)  
[json_read_dict/2](#json_read_dict/2)  
[json_read_dict/3](#json_read_dict/3)  
[json_rpc_dispatch/2](#json_rpc_dispatch/2)  
[json_rpc_dispatch_request/4](#json_rpc_dispatch_request/4)  
[json_rpc_error/2](#json_rpc_error/2)  
[json_rpc_error/3](#json_rpc_error/3)  
[json_to_prolog/2](#json_to_prolog/2)  
[json_type/1](#json_type/1)  
[json_validate/3](#json_validate/3)  
[json_write/2](#json_write/2)  
[json_write/3](#json_write/3)  
[json_write_dict/2](#json_write_dict/2)  
[json_write_dict/3](#json_write_dict/3)  
[json_write_hook/4](#json_write_hook/4)  
[prolog_bool_to_json/2](#prolog_bool_to_json/2)  
[prolog_to_json/2](#prolog_to_json/2)  
[reply_json/1](#reply_json/1)  
[reply_json/2](#reply_json/2)  
[reply_json_dict/1](#reply_json_dict/1)  
[reply_json_dict/2](#reply_json_dict/2)  
