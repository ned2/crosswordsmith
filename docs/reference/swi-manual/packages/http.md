SWI-Prolog HTTP support

Jan Wielemaker  
VU University Amsterdam  
University of Amsterdam  
The Netherlands  
E-mail: [J.Wielemaker@vu.nl](mailto:J.Wielemaker@vu.nl)

Abstract

This article documents the package HTTP, a series of libraries for accessing data on HTTP servers as well as providing HTTP server capabilities from SWI-Prolog (up to HTTP 1.1). Both server and client are modular libraries. Further reading material is available from the locations below.

- [HOWTO collection](http://www.swi-prolog.org/howto/http/)
- [Tutorial by Anne Ogborn](https://github.com/Anniepoo/swiplwebtut/blob/master/web.adoc)
- [Wikipedia entry on HTTP](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol)
- [Wikipedia entry on REST](https://en.wikipedia.org/wiki/Representational_state_transfer)

# Table of Contents

[1 Introduction](#sec:1)

[2 The HTTP client libraries](#sec:2)

[2.1 library(http/http_open): HTTP client library](#sec:2.1)

[2.2 library(http/http_client): HTTP client library](#sec:2.2)

[3 The HTTP server libraries](#sec:3)

[3.1 Creating an HTTP reply](#sec:3.1)

[3.1.1 Returning special status codes](#sec:3.1.1)

[3.2 library(http/http_dispatch): Dispatch requests in the HTTP server](#sec:3.2)

[3.3 library(http/http_dirindex): HTTP directory listings](#sec:3.3)

[3.4 library(http/http_files): Serve plain files from a hierarchy](#sec:3.4)

[3.5 library(http/http_session): HTTP Session management](#sec:3.5)

[3.6 library(http/http_cors): Enable CORS: Cross-Origin Resource Sharing](#sec:3.6)

[3.7 library(http/http_authenticate): Authenticate HTTP connections using 401 headers](#sec:3.7)

[3.8 library(http/http_digest): HTTP Digest authentication](#sec:3.8)

[3.9 library(http/http_dyn_workers): Dynamically schedule HTTP workers.](#sec:3.9)

[3.9.1 Providing Server-Sent Events (sse)](#sec:3.9.1)

[3.10 Custom Error Pages](#sec:3.10)

[3.11 library(http/http_openid): OpenID consumer and server library](#sec:3.11)

[3.12 Get parameters from HTML forms](#sec:3.12)

[3.13 Request format](#sec:3.13)

[3.13.1 Handling POST requests](#sec:3.13.1)

[3.14 Running the server](#sec:3.14)

[3.14.1 Common server interface options](#sec:3.14.1)

[3.14.2 Multi-threaded Prolog](#sec:3.14.2)

[3.14.3 library(http/http_unix_daemon): Run SWI-Prolog HTTP server as a Unix system daemon](#sec:3.14.3)

[3.14.4 From (Unix) inetd](#sec:3.14.4)

[3.14.5 MS-Windows](#sec:3.14.5)

[3.14.6 As CGI script](#sec:3.14.6)

[3.14.7 Using a reverse proxy](#sec:3.14.7)

[3.15 The wrapper library](#sec:3.15)

[3.16 library(http/http_host): Obtain public server location](#sec:3.16)

[3.17 library(http/http_log): HTTP Logging module](#sec:3.17)

[3.18 library(http/http_server_health): HTTP Server health statistics](#sec:3.18)

[3.19 Debugging HTTP servers](#sec:3.19)

[3.20 library(http/http_header): Handling HTTP headers](#sec:3.20)

[3.21 The `library(http/html_write)` library](#sec:3.21)

[3.21.1 Emitting HTML documents](#sec:3.21.1)

[3.21.2 Repositioning HTML for CSS and javascript links](#sec:3.21.2)

[3.21.3 Adding rules for html//1](#sec:3.21.3)

[3.21.4 Generating layout](#sec:3.21.4)

[3.21.5 Examples for using the HTML write library](#sec:3.21.5)

[3.21.6 Remarks on the `library(http/html_write)` library](#sec:3.21.6)

[3.22 library(http/js_write): Utilities for including JavaScript](#sec:3.22)

[3.23 library(http/http_path): Abstract specification of HTTP server locations](#sec:3.23)

[3.24 library(http/html_head): Automatic inclusion of CSS and scripts links](#sec:3.24)

[3.24.1 About resource ordering](#sec:3.24.1)

[3.24.2 Debugging dependencies](#sec:3.24.2)

[3.24.3 Predicates](#sec:3.24.3)

[3.25 library(http/http_pwp): Serve PWP pages through the HTTP server](#sec:3.25)

[3.26 library(http/htmx): Support htmx.org](#sec:3.26)

[4 HTTP and IPv6](#sec:4)

[5 Transfer encodings](#sec:5)

[5.1 The `library(http/http_chunked)` library](#sec:5.1)

[6 library(http/websocket): WebSocket support](#sec:6)

[7 library(http/hub): Manage a hub for websockets](#sec:7)

[8 MIME support](#sec:8)

[8.1 library(http/mimepack): Create a MIME message](#sec:8.1)

[9 Security](#sec:9)

[10 Tips and tricks](#sec:10)

[11 Status](#sec:11)

## 1 Introduction

HTTP (Hypertext Transfer Protocol) is the W3C standard protocol for transferring information between a web-client (e.g., a browser) and a web-server. The protocol is a simple *envelope* protocol where standard name/value pairs in the header are used to split the stream into messages and communicate about the connection-status. Many languages have client and server libraries to deal with the HTTP protocol, making this protocol an excellent candidate for building client-server applications. In particular, HTTP is a natural fit for networked systems built according to the principles of “Representational State Transfer” (*REST*).

In this document we describe a modular infrastructure to access web-servers from SWI-Prolog and turn Prolog into a web-server.

### Acknowledgements

This work has been carried out under the following projects: [GARP](http://hcs.science.uva.nl/projects/GARP/), [MIA](http://www.ins.cwi.nl/projects/MIA/) (dead link), [IBROW](http://hcs.science.uva.nl/projects/ibrow/home.html) (dead link), [KITS](http://kits.edte.utwente.nl/) (dead link) and [MultiMediaN](http://e-culture.multimedian.nl/) (dead link).

The following people have pioneered parts of this library and contributed with bug reports and suggestions for improvements: Anjo Anjewierden, Bert Bredeweg, Wouter Jansweijer, Bob Wielinga, Jacco van Ossenbruggen, Michiel Hildebrandt, Matt Lilley and Keri Harris.

*Path wildcards* (see [http_handler/3](#http_handler/3)) have been modelled after the “arouter” add-on pack by Raivo Laanemets. *Request rewriting* has been added after discussion with Raivo Laanemets and Anne Ogborn on the SWI-Prolog mailinglist.

## 2 The HTTP client libraries

This package provides two client libraries for accessing HTTP servers.

**`library(http/http_open)`**  
This library provides [http_open/3](#http_open/3) and friends. It is a library for opening an endpoint identified by an HTTP URL as a Prolog stream. The general skeleton for using this library is given below, where process/1 processes the data from the HTTP server.^(1One may opt to use cleanup/2 instead of setup_call_cleanup/3 to allow for aborting while [http_open/3](#http_open/3) is waiting for the connection.)

``` code
    setup_call_cleanup(
        http_open(URL, In, []),
        process(In),
        close(In)).
```

**`library(http/http_client)`**  
This library provides [http_get/3](#http_get/3) and [http_post/4](#http_post/4) and friends. These predicates process the reply using plugins to convert the data based on the `Content-Type` of the reply. This library supports a plugin infrastructure that can register hooks for converting additional document types.

### 2.1 library(http/http_open): HTTP client library

See also  
\- load_html/3 and xpath/3 can be used to parse and navigate HTML documents.  
- [http_get/3](#http_get/3) and [http_post/4](#http_post/4) provide an alternative interface that convert the reply depending on the `Content-Type` header.

This library defines [http_open/3](#http_open/3), which opens an URL as a Prolog stream. The functionality of the library can be extended by loading two additional modules that act as plugins:

**library**(`http/http_ssl_plugin`)  
Loading this library causes [http_open/3](#http_open/3) to handle HTTPS connections. Relevant options for SSL certificate handling are handed to ssl_context/3. This plugin is loaded automatically if the scheme `https` is requested using a default SSL context. See the plugin for additional information regarding security.

**library**(`zlib`)  
Loading this library supports the `gzip` transfer encoding. This plugin is lazily loaded if a connection is opened that claims this transfer encoding.

**library**(`http/http_cookie`)  
Loading this library adds tracking cookies to [http_open/3](#http_open/3). Returned cookies are collected in the Prolog database and supplied for subsequent requests.

**library**(`http/http_stream`)  
This library adds support for *chunked* encoding. It is lazily loaded if the server sends a `Transfer-encoding: chunked` header.

Here is a simple example to fetch a web-page:

``` code
?- http_open('http://www.google.com/search?q=prolog', In, []),
   copy_stream_data(In, user_output),
   close(In).
<!doctype html><head><title>prolog - Google Search</title><script>
...
```

The example below fetches the modification time of a web-page. Note that `Modified` is `''` (the empty atom) if the web-server does not provide a time-stamp for the resource. See also parse_time/2.

``` code
modified(URL, Stamp) :-
       http_open(URL, In,
                 [ method(head),
                   header(last_modified, Modified)
                 ]),
       close(In),
       Modified \== '',
       parse_time(Modified, Stamp).
```

Then next example uses Google search. It exploits `library(uri)` to manage URIs, `library(sgml)` to load an HTML document and `library(xpath)` to navigate the parsed HTML. Note that you may need to adjust the XPath queries if the data returned by Google changes (this example indeed no longer works and currently fails at the first xpath/3 call)

``` code
:- use_module(library(http/http_open)).
:- use_module(library(xpath)).
:- use_module(library(sgml)).
:- use_module(library(uri)).

google(For, Title, HREF) :-
        uri_encoded(query_value, For, Encoded),
        atom_concat('http://www.google.com/search?q=', Encoded, URL),
        http_open(URL, In, []),
        call_cleanup(
            load_html(In, DOM, []),
            close(In)),
        xpath(DOM, //h3(@class=r), Result),
        xpath(Result, //a(@href=HREF0, text), Title),
        uri_components(HREF0, Components),
        uri_data(search, Components, Query),
        uri_query_components(Query, Parts),
        memberchk(q=HREF, Parts).
```

An example query is below:

``` code
?- google(prolog, Title, HREF).
Title = 'SWI-Prolog',
HREF = 'http://www.swi-prolog.org/' ;
Title = 'Prolog - Wikipedia',
HREF = 'https://nl.wikipedia.org/wiki/Prolog' ;
Title = 'Prolog - Wikipedia, the free encyclopedia',
HREF = 'https://en.wikipedia.org/wiki/Prolog' ;
Title = 'Pro-Log is logistiek dienstverlener m.b.t. vervoer over water.',
HREF = 'http://www.pro-log.nl/' ;
Title = 'Learn Prolog Now!',
HREF = 'http://www.learnprolognow.org/' ;
Title = 'Free Online Version - Learn Prolog
...
```

\[det\]**http_open**(`+URL, -Stream, +Options`)  
Open the data at the HTTP server as a Prolog stream. `URL` is either an atom specifying a `URL` or a list representing a broken-down `URL` as specified below. After this predicate succeeds the data can be read from `Stream`. After completion this stream must be closed using the built-in Prolog predicate close/1. `Options` provides additional options:

**authenticate**(`+Boolean`)  
If `false` (default `true`), do *not* try to automatically authenticate the client if a 401 (Unauthorized) status code is received.

**authorization**(`+Term`)  
Send authorization. See also [http_set_authorization/2](#http_set_authorization/2). Supported schemes:

**basic**(`+User, +Password`)  
HTTP Basic authentication.

**bearer**(`+Token`)  
HTTP Bearer authentication.

**digest**(`+User, +Password`)  
HTTP Digest authentication. This option is only provided if the plugin `library(http/http_digest)` is also loaded.

**unix_socket**(`+Path`)  
Connect to the given Unix domain socket. In this scenario the host name and port or ignored. If the server replies with a *redirect* message and the host differs from the original host as normal TCP connection is used to handle the redirect. This option is inspired by `curl(1)`’s option‘--unix-socket\`.

**connection**(`+Connection`)  
Specify the `Connection` header. Default is `close`. The alternative is `Keep-alive`. This maintains a pool of available connections as determined by keep_connection/1. The `library(http/websockets)` uses `Keep-alive, Upgrade`. Keep-alive connections can be closed explicitly using [http_close_keep_alive/1](#http_close_keep_alive/1). Keep-alive connections may significantly improve repetitive requests on the same server, especially if the IP route is long, HTTPS is used or the connection uses a proxy.

**final_url**(`-FinalURL`)  
Unify `FinalURL` with the final destination. This differs from the original `URL` if the returned head of the original indicates an HTTP redirect (codes 301, 302 or 303). Without a redirect, `FinalURL` is the same as `URL` if `URL` is an atom, or a `URL` constructed from the parts.

**header**(`Name, -AtomValue`)  
If provided, `AtomValue` is unified with the value of the indicated field in the reply header. `Name` is matched case-insensitive and the underscore (\_) matches the hyphen (-). Multiple of these options may be provided to extract multiple header fields. If the header is not available `AtomValue` is unified to the empty atom (” ).

**headers**(`-List`)  
If provided, `List` is unified with a list of Name(Value) pairs corresponding to fields in the reply header. Name and Value follow the same conventions used by the `header(Name,Value)` option. A pseudo header `status_code(Code)` is added to provide the HTTP status as an integer. See also `raw_headers(-List)` which provides the entire HTTP reply header in unparsed representation.

**method**(`+Method`)  
One of `get` (default), `head`, `delete`, `post`, `put` or `patch`. The `head` message can be used in combination with the `header(Name, Value)` option to access information on the resource without actually fetching the resource itself. The returned stream must be closed immediately.

If `post(Data)` is provided, the default is `post`.

**size**(`-Size`)  
`Size` is unified with the integer value of `Content-Length` in the reply header.

**version**(`-Version`)  
`Version` is a *pair* `Major-Minor`, where `Major` and `Minor` are integers representing the HTTP version in the reply header.

**range**(`+Range`)  
Ask for partial content. `Range` is a term *Unit(From,To)*, where `From` is an integer and `To` is either an integer or the atom `end`. HTTP 1.1 only supports Unit = `bytes`. E.g., to ask for bytes 1000-1999, use the option `range(bytes(1000,1999))`

**raw_encoding**(`+Encoding`)  
Do not install a decoding filter for `Encoding`. For example, using `raw_encoding('applocation/gzip')` the system will not decompress the stream if it is compressed using `gzip`.

**raw_headers**(`-Lines`)  
Unify `Lines` with a list of strings that represents the complete reply header returned by the server. See also `headers(-List)`.

**redirect**(`+Boolean`)  
If `false` (default `true`), do *not* automatically redirect if a 3XX code is received. Must be combined with `status_code(Code)` and one of the header options to read the redirect reply. In particular, without `status_code(Code)` a redirect is mapped to an exception.

**status_code**(`-Code`)  
If this option is present and `Code` unifies with the HTTP status code, do **not** translate errors (4xx, 5xx) into an exception. Instead, [http_open/3](#http_open/3) behaves as if 2xx (success) is returned, providing the application to read the error document from the returned stream.

**output**(`-Out`)  
Unify the output stream with `Out` and do not close it. This can be used to upgrade a connection.

**timeout**(`+Timeout`)  
If provided, set a timeout on the stream using set_stream/2. With this option if no new data arrives within `Timeout` seconds the stream raises an exception. Default is to wait forever (`infinite`).

**post**(`+Data`)  
Issue a `POST` request on the HTTP server. `Data` is handed to [http_post_data/3](#http_post_data/3).

**proxy**(`+Host:Port`)  
Use an HTTP proxy to connect to the outside world. See also socket:proxy_for_url/3. This option overrules the proxy specification defined by socket:proxy_for_url/3.

**proxy**(`+Host, +Port`)  
Synonym for `proxy(+Host:Port)`. Deprecated.

**proxy_authorization**(`+Authorization`)  
Send authorization to the proxy. Otherwise the same as the `authorization` option.

**bypass_proxy**(`+Boolean`)  
If `true`, bypass proxy hooks. Default is `false`.

**request_header**(`Name=Value`)  
Additional name-value parts are added in the order of appearance to the HTTP request header. No interpretation is done.

**max_redirect**(`+Max`)  
Sets the maximum length of a redirection chain. This is needed for some IRIs that redirect indefinitely to other IRIs without looping (e.g., redirecting to IRIs with a random element in them). `Max` must be either a non-negative integer or the atom `infinite`. The default value is `10`.

**user_agent**(`+Agent`)  
Defines the value of the `User-Agent` field of the HTTP header. Default is `SWI-Prolog`.

The hook [http:open_options/2](#http:open_options/2) can be used to provide default options based on the broken-down `URL`. The option `status_code(-Code)` is particularly useful to query **REST** interfaces that commonly return status codes other than `200` that need to be be processed by the client code.

[TABLE]

throws  
`error(existence_error(url, Id),Context)` is raised if the HTTP result code is not in the range 200..299. Context has the shape `context(Message, status(Code, TextCode))`, where `Code` is the numeric HTTP code and `TextCode` is the textual description thereof provided by the server. `Message` may provide additional details or may be unbound.

See also  
ssl_context/3 for SSL related options if `library(http/http_ssl_plugin)` is loaded.

\[multifile\]**map_method**(`+MethodID, -Method`)  
Support additional `METHOD` keywords. Default are the official HTTP methods as defined by the various RFCs.

\[semidet,multifile\]http:**disable_encoding_filter**(`+ContentType`)  
Do not use the `Content-encoding` as `Transfer-encoding` encoding for specific values of `ContentType`. This predicate is multifile and can thus be extended by the user.

\[det\]**http_set_authorization**(`+URL, +Authorization`)  
Set user/password to supply with URLs that have `URL` as prefix. If `Authorization` is the atom `-`, possibly defined authorization is cleared. For example:

``` code
?- http_set_authorization('http://www.example.com/private/',
                          basic('John', 'Secret'))
```

To be done  
Move to a separate module, so [http_get/3](#http_get/3), etc. can use this too.

\[semidet,multifile\]iostream:**open_hook**(`+Spec, +Mode, -Stream, -Close, +Options0, -Options`)  
Hook implementation that makes open_any/5 support `http` and `https` URLs for `Mode == read`.

\[det\]**keep_alive**(`+StreamPair, +Host, +In, -Left`)  
Callback when closing the range stream used to process the content of the reply. This callback makes the stream available for future keep-alive connections or closes the stream. The stream is closed if

- There are too many bytes left unprocessed in the range stream.
- There are too many pooled connections.

\[det\]**http_close_keep_alive**(`+Address`)  
Close all keep-alive connections matching `Address`. `Address` is of the form Host:Port. In particular, `http_close_keep_alive(_)` closes all currently known keep-alive connections.

\[nondet,multifile\]http:**open_options**(`+Parts, -Options`)  
This hook is used by the HTTP client library to define default options based on the the broken-down request-URL. The following example redirects all traffic, except for localhost over a proxy:

``` code
:- multifile
    http:open_options/2.

http:open_options(Parts, Options) :-
    option(host(Host), Parts),
    Host \== localhost,
    Options = [proxy('proxy.local', 3128)].
```

This hook may return multiple solutions. The returned options are combined using merge_options/3 where earlier solutions overrule later solutions.

\[semidet,multifile\]http:**write_cookies**(`+Out, +Parts, +Options`)  
Emit a `Cookie:` header for the current connection. `Out` is an open stream to the HTTP server, `Parts` is the broken-down request (see uri_components/2) and `Options` is the list of options passed to http_open. The predicate is called as if using ignore/1.

See also  
\- complements [http:update_cookies/3](#http:update_cookies/3).  
- `library(http/http_cookie)` implements cookie handling on top of these hooks.

\[semidet,multifile\]http:**update_cookies**(`+CookieData, +Parts, +Options`)  
Update the cookie database. `CookieData` is the value of the `Set-Cookie` field, `Parts` is the broken-down request (see uri_components/2) and `Options` is the list of options passed to http_open.

See also  
\- complements http:write_cookies  
- `library(http/http_cookies)` implements cookie handling on top of these hooks.

### 2.2 library(http/http_client): HTTP client library

This library provides the four basic HTTP client actions: `GET`, `DELETE`, `POST` and `PUT`. In addition, it provides [http_read_data/3](#http_read_data/3), which is used by `library(http/http_parameters)` to decode `POST` data in server applications.

This library is based on [http_open/3](#http_open/3), which opens a URL as a Prolog stream. The reply is processed by [http_read_data/3](#http_read_data/3). The following content-types are supported. Options passed to [http_get/3](#http_get/3) and friends are passed to [http_read_data/3](#http_read_data/3), which in turn passes them to the conversion predicates. Support for additional content types can be added by extending the multifile predicate http_client:http_convert_data/4.

**application/x-www-form-urlencoded**  
Built in. Converts form-data into a list of `Name=Value` terms.

**application/x-prolog**  
Built in. Reads a single Prolog term.

**multipart/form-data**  
Processed if `library(http/http_multipart_plugin)` is loaded. This format should be used to handle web forms that upload a file.

`text/html` **`|`** `text/xml`  
Processed if `library(http/http_sgml_plugin)` is loaded. See load_html/3 for details and load_xml/3 for details. The output is often processed using xpath/3.

`application/json` **`|`** `application/jsonrequest`  
Processed if `library(http/http_json)` is loaded. The option `json_object(As)` can be used to return a term `json(Attributes)` (`As` is `term`) or a dict (`As` is `dict`).

&nbsp;

\[det\]**http_get**(`+URL, -Data, +Options`)  
Get data from a `URL` server and convert it to a suitable Prolog representation based on the `Content-Type` header and plugins. This predicate is the common implementation of the HTTP client operations. The predicates [http_delete/3](#http_delete/3), [http_post/4](#http_post/4) and [http_put/4](#http_put/4) call this predicate with an appropriate `method(+Method)` option and ---for [http_post/4](#http_post/4) and [http_put/4](#http_put/4)--- a `post(+Data)` option.

`Options` are passed to [http_open/3](#http_open/3) and [http_read_data/3](#http_read_data/3). Other options:

**reply_header**(`-Fields`)  
Synonym for `headers(Fields)` from [http_open/3](#http_open/3). Provided for backward compatibility. Note that `http_version(Major-Minor)` is missing in the new version.

\[det\]**http_delete**(`+URL, -Data, +Options`)  
Execute a `DELETE` method on the server. Arguments are the same as for [http_get/3](#http_get/3). Typically one should pass the option `status_code(-Code)` to assess and evaluate the returned status code. Without, codes other than 200 are interpreted as an error.

See also  
Implemented on top of [http_get/3](#http_get/3).

To be done  
Properly map the 201, 202 and 204 replies.

\[det\]**http_post**(`+URL, +Data, -Reply, +Options`)  
Issue an HTTP `POST` request. `Data` is posted using [http_post_data/3](#http_post_data/3). The HTTP server reply is returned in `Reply`, using the same rules as for [http_get/3](#http_get/3).

See also  
Implemented on top of [http_get/3](#http_get/3).

**http_put**(`+URL, +Data, -Reply, +Options`)  
Issue an HTTP `PUT` request. Arguments are the same as for [http_post/4](#http_post/4).

See also  
Implemented on top of [http_post/4](#http_post/4).

**http_patch**(`+URL, +Data, -Reply, +Options`)  
Issue an HTTP `PATCH` request. Arguments are the same as for [http_post/4](#http_post/4).

See also  
Implemented on top of [http_post/4](#http_post/4).

\[det\]**http_read_data**(`+Request, -Data, +Options`)  
Read data from an HTTP connection and convert it according to the supplied `to(Format)` option or based on the `Content-type` in the `Request`. The following options are supported:

**to**(`Format`)  
Convert data into `Format`. Values are:

- `stream(+WriteStream)`) Append the content of the message to Stream
- atom Return the reply as an atom
- string Return the reply as a string
- codes Return the reply as a list of codes

**form_data**(`AsForm`)  
**input_encoding**(`+Encoding`)  
**on_filename**(`:CallBack`)  
These options are implemented by the plugin `library(http/http_multipart_plugin)` and apply to processing `multipart/form-data` content.

**content_type**(`+Type`)  
Overrule the content-type that is part of `Request` as a work-around for wrongly configured servers.

Without plugins, this predicate handles

**application/x-www-form-urlencoded**  
Converts form-data into a list of `Name=Value` terms.

**application/x-prolog**  
Converts data into a Prolog term.

|  |  |
|----|----|
| `Request` | is a parsed HTTP request as returned by [http_read_request/2](#http_read_request/2) or available from the HTTP server's request dispatcher. `Request` must contain a term `input(In)` that provides the input stream from the HTTP server. |

\[semidet,multifile\]**http_convert_data**(`+In, +Fields, -Data, +Options`)  
Multi-file hook to convert a HTTP payload according to the *Content-Type* header. The default implementation deals with application/x-prolog. The HTTP framework provides implementations for JSON (`library(http/http_json)`), HTML/XML (`library(http/http_sgml_plugin)`)

\[det\]**http_disconnect**(`+Connections`)  
Close down some connections. Currently `Connections` must have the value `all`, closing all connections.

deprecated  
New code should use [http_close_keep_alive/1](#http_close_keep_alive/1) from `library(http/http_open)`.

\[semidet,multifile\]http:**post_data_hook**(`+Term, +Out, +Options`)  
Hook to extend the datatypes supported by the `post(Data)` option of [http_open/3](#http_open/3). The default implementation supports `prolog(Term)`, sending a Prolog term as `application/x-prolog`.

## 3 The HTTP server libraries

The HTTP server infra structure consists of a number of small modular libraries that are combined into `library(http/http_server)`. These modules are:

**`library(http/thread_httpd)`**  
This library is responsible for accepting and managing connections.^(2In older versions there were two alternative libraries for managing connections based on XPCE and Unix inetd.)

**`library(http/http_dyn_workers)`**  
This library dynamically adds and removes workers based on the workload of the server.

**`library(http/http_wrapper)`**  
This library takes a connection, parses the HTTP request header and runs a goal that produces a CGI document based on the parsed request. It watches for exceptions and turns these into (error) status pages. The status page generation may be hooked to provide custom pages.

**`library(http/http_dispatch)`**  
This library associates the *path* of the HTTP request with a *handler* that services this particular path. It also manages timeouts and may pass the execution of a request to a dedicated thread with specified resource limits using [http_spawn/2](#http_spawn/2). The module supports plugable request rewrite handlers that may be used to implement identification, authorization, input argument processing, etc.

**`library(http/http_parameters)`**  
This library parses HTTP request parameters, both dealing with GET and POST style parameter passing.

**`library(http/html_write)`**  
This library translates a Prolog term into an HTML document using Prolog *grammar rules* (DCG). It provides a modular infrastructure to build pages that are guaranteed to be valid HTML. The HTTP server libraries provide several alternatives for generating HTML ranging from simple printing to `current_output` to XML-based templates (PWP).

**`library(http/http_json)`**  
This library parses a POSTed HTTP document into a Prolog dict and formulates an HTTP JSON reply from a Prolog dict and is typically used to implement REST services.

Most server implementation simply load the `library(http/http_server)` library, which loads the above modules and reexports all predicates except for those used for internal communication and older deprecated predicates. Specific use cases may load a subset of the individual libraries and may decide to replace one or more of them.

A typical skeleton for building a server is given below. If this file is loaded as main file (using e.g., `swipl server.pl`) it creates a simple server that listens on port 8080. If the root is accessed it redirects to the home page and shows **Hello world!**.

``` code
:- use_module(library(http/http_server)).

:- initialization
    http_server([port(8080)]).

:- http_handler(root(.),
                http_redirect(moved, location_by_id(home_page)),
                []).
:- http_handler(root(home), home_page, []).

home_page(_Request) :-
    reply_html_page(
        title('Demo server'),
        [ h1('Hello world!')
        ]).
```

### 3.1 Creating an HTTP reply

The *handler* (e.g., home_page/1 above) is called with the parsed request (see [section 3.13](#sec:3.13)) as argument and `current_output` set to a temporary buffer. Its task is closely related to the task of a CGI script; it must write a header declaring at least the `Content-type` field and a body. Below is a simple body writing the request as an HTML table.^(3Note that writing an HTML reply this way is deprecated. In fact, the code is subject to *injection attacks* as the HTTP request field values are literally injected in the output while HTML reserved characters should be properly escaped.)

``` code
reply(Request) :-
        format('Content-type: text/html~n~n', []),
        format('<html>~n', []),
        format('<table border=1>~n'),
        print_request(Request),
        format('~n</table>~n'),
        format('</html>~n', []).

print_request([]).
print_request([H|T]) :-
        H =.. [Name, Value],
        format('<tr><td>~w<td>~w~n', [Name, Value]),
        print_request(T).
```

The infrastructure recognises the header fields described below. Other header lines are passed verbatim to the client. Typical examples are `Set-Cookie` and authentication headers (see [section 3.7](#sec:3.7)).

**Content-type: `Type`**  
This field is passed to the client and used by the infrastructure to determine the *encoding* to use for the stream. If `type` matches `text/*` or the type matches with `UTF-8` (case insensitive), the server uses UTF-8 encoding. The user may force UTF-8 encoding for arbitrary content types by adding `; charset=UTF-8` to the end of the `Content-type` header.

**Transfer-encoding: chunked**  
Causes the server to use *chunked* encoding if the client allows for it. See also [section 5](#sec:5) and the `chunked` option in [http_handler/3](#http_handler/3).

**Connection: close**  
Causes the connection to be closed after the transfer. The default is to keep it open‘Keep-Alive’if possible.

**Location: `URL`**  
This header may be combined with the `Status` header to force a *redirect* response to the given `URL`. The message body must be empty. Handling this header is primarily intended for compatibility with the CGI conventions. Prolog code should use [http_redirect/3](#http_redirect/3).

**Status: `Status`**  
This header can be combined with `Location`, where `Status` must be one of 301 (moved), 302 (moved temporary, default) or 303 (see other). Using the status field also allows for formulating replies such as 201 (created).

Note that the handler may send any type of document instead of HTML. After the header has been written, the *encoding* of the `current_output` stream encoding is established as follows:

1.  If the content type is `text/*` the stream is switched to UTF-8 encoding. If the content type does not provide attributes, `; charset=UTF-8` is added.
2.  The content type contains `UTF-8` the stream is switched to UTF-8 encoding.
3.  [http:mime_type_encoding/2](#http:mime_type_encoding/2) succeeds the returned encoding is used. The returned encoding must be valid for set_stream/2.
4.  If the content type matches a list of known encodings, this is used. See mime_type_encoding/2 is `http_header`. The current list deals with JSON, Turtle and SPARQL.
5.  Otherwise the stream uses octed (binary) encoding.

#### 3.1.1 Returning special status codes

Besides returning a page by writing it to the current output stream, the server goal can raise an exception using throw/1 to generate special pages such as `not_found`, `moved`, etc. The defined exceptions are:

**http_reply**(`+Reply, +HdrExtra, +Context`)  
Return a result page using [http_reply/3](#http_reply/3). See [http_reply/3](#http_reply/3) for supported values for Reply and [section 3.10](#sec:3.10) for providing a custom error page.

**http_reply**(`+Reply, +HdrExtra`)  
Return a result page using [http_reply/3](#http_reply/3). Equivalent to `http_reply(Reply, HdrExtra,[])`.

**http_reply**(`+Reply`)  
Equivalent to `http_reply(Reply, [],[])`.

**http**(`not_modified`)  
Equivalent to `http_reply(not_modified,[])`. This exception is for backward compatibility and can be used by the server to indicate the referenced resource has not been modified since it was requested last time.

In addition, the normal `"200 OK"` reply status may be overruled by writing a CGI `Status` header prior to the remainder of the message. This is particularly useful for defining REST APIs. The following handler replies with a `"201 Created"` header:

``` code
handle_request(Request) :-
        process_data(Request, Id),      % application predicate
        format('Status: 201~n'),
        format('Content-type: text/plain~n~n'),
        format('Created object as ~q~n', [Id]).
```

### 3.2 library(http/http_dispatch): Dispatch requests in the HTTP server

Most code doesn't need to use this directly; instead use `library(http/http_server)`, which combines this library with the typical HTTP libraries that most servers need.

This module can be placed between `http_wrapper.pl` and the application code to associate HTTP *locations* to predicates that serve the pages. In addition, it associates parameters with locations that deal with timeout handling and user authentication. The typical setup is:

``` code
server(Port, Options) :-
        http_server(http_dispatch,
                    [ port(Port)
                    | Options
                    ]).

:- http_handler('/index.html', write_index, []).

write_index(Request) :-
        ...
```

\[det\]**http_handler**(`+Path, :Closure, +Options`)  
Register `Closure` as a handler for HTTP requests. `Path` is either an absolute path such as `'/home.html'` or a term Alias(Relative). Where Alias is associated with a concrete path using [http:location/3](#http:location/3) and resolved using [http_absolute_location/3](#http_absolute_location/3). `Relative` can be a single atom or a term‘Segment1/Segment2/...\`, where each element is either an atom or a variable. If a segment is a variable it matches any segment and the binding may be passed to the closure. If the last segment is a variable it may match multiple segments. This allows registering REST paths, for example:

``` code
:- http_handler(root(user/User), user(Method, User),
                [ method(Method),
                  methods([get,post,put])
                ]).

user(get, User, Request) :-
    ...
user(post, User, Request) :-
    ...
```

If an HTTP request arrives at the server that matches `Path`, `Closure` is called as below, where `Request` is the parsed HTTP request.

``` code
call(Closure, Request)
```

`Options` is a list containing the following options:

**authentication**(`+Type`)  
Demand authentication. Authentication methods are plugable. The library `http_authenticate.pl` provides a plugin for user/password based `Basic` HTTP authentication.

**chunked**  
Use `Transfer-encoding: chunked` if the client allows for it.

**condition**(`:Goal`)  
If present, the handler is ignored if `Goal` does not succeed.

**content_type**(`+Term`)  
Specifies the content-type of the reply. This value is currently not used by this library. It enhances the reflexive capabilities of this library through [http_current_handler/3](#http_current_handler/3).

**id**(`+Atom`)  
Identifier of the handler. The default identifier is the predicate name. Used by [http_location_by_id/2](#http_location_by_id/2) and [http_link_to_id/3](#http_link_to_id/3).

**hide_children**(`+Bool`)  
If `true` on a prefix-handler (see prefix), possible children are masked. This can be used to (temporary) overrule part of the tree.

**method**(`+Method`)  
Declare that the handler processes `Method`. This is equivalent to `methods([Method])`. Using `method(*)` allows for all methods.

**methods**(`+ListOfMethods`)  
Declare that the handler processes all of the given methods. If this option appears multiple times, the methods are combined.

**prefix**  
Call Pred on any location that is a specialisation of `Path`. If multiple handlers match, the one with the longest path is used. `Options` defined with a prefix handler are the default options for paths that start with this prefix. Note that the handler acts as a fallback handler for the tree below it:

``` code
:- http_handler(/, http_404([index('index.html')]),
                [spawn(my_pool),prefix]).
```

**priority**(`+Integer`)  
If two handlers handle the same path, the one with the highest priority is used. If equal, the last registered is used. Please be aware that the order of clauses in multifile predicates can change due to reloading files. The default priority is 0 (zero).

**spawn**(`+SpawnOptions`)  
Run the handler in a separate thread. If `SpawnOptions` is an atom, it is interpreted as a thread pool name (see create_thread_pool/3). Otherwise the options are passed to [http_spawn/2](#http_spawn/2) and from there to thread_create/3. These options are typically used to set the stack limits.

**time_limit**(`+Spec`)  
One of `infinite`, `default` or a positive number (seconds). If `default`, the value from the setting `http:time_limit` is taken. The default of this setting is 300 (5 minutes). See setting/2.

Note that [http_handler/3](#http_handler/3) is normally invoked as a directive and processed using term-expansion. Using term-expansion ensures proper update through make/0 when the specification is modified.

Errors  
\- `existence_error(http_location, Location)`  
- `permission_error(http_method, Method, Location)`

See also  
[http_reply_file/3](#http_reply_file/3) and [http_redirect/3](#http_redirect/3) are generic handlers to serve files and achieve redirects.

\[det\]**http_delete_handler**(`+Spec`)  
Delete handler for `Spec`. Typically, this should only be used for handlers that are registered dynamically. `Spec` is one of:

**id**(`Id`)  
Delete a handler with the given id. The default id is the handler-predicate-name.

**path**(`Path`)  
Delete handler that serves the given path.

\[det\]**http_dispatch**(`Request`)  
Dispatch a `Request` using [http_handler/3](#http_handler/3) registrations. It performs the following steps:

1.  Find a matching handler based on the `path` member of `Request`. If multiple handlers match due to the `prefix` option or variables in path segments (see [http_handler/3](#http_handler/3)), the longest specification is used. If multiple specifications of equal length match the one with the highest priority is used.
2.  Check that the handler matches the `method` member of the `Request` or throw `permission_error(http_method, Method, Location)`
3.  Expand the request using expansion hooks registered by http_request_expansion/3. This may add fields to the request, such the authenticated user, parsed parameters, etc. The hooks may also throw exceptions, notably using [http_redirect/3](#http_redirect/3) or by throwing `http_reply(Term, ExtraHeader, Context)` exceptions.
4.  Extract possible fields from the `Request` using e.g. `method(Method)` as one of the options.
5.  Call the registered *closure*, optionally spawning the request to a new thread or enforcing a time limit.

**http_request_expansion**(`:Goal, +Rank:number`)  
Register `Goal` for expanding the HTTP request handler. `Goal` is called as below. If `Goal` fail the request is passed to the next expansion unmodified.

``` code
call(Goal, Request0, Request, Options)
```

If multiple goals are registered they expand the request in a pipeline starting with the expansion hook with the lowest rank.

Besides rewriting the request, for example by validating the user identity based on HTTP authentication or cookies and adding this to the request, the hook may raise HTTP exceptions to indicate a bad request, permission error, etc. See [http_status_reply/4](#http_status_reply/4).

Initially, auth_expansion/3 is registered with rank `100` to deal with the older [http:authenticate/3](#http:authenticate/3) hook.

\[semidet\]**http_current_handler**(`+Location, :Closure`)  
\[nondet\]**http_current_handler**(`-Location, :Closure`)  
True if `Location` is handled by `Closure`.

\[semidet\]**http_current_handler**(`+Location, :Closure, -Options`)  
\[nondet\]**http_current_handler**(`?Location, :Closure, ?Options`)  
Resolve the current handler and options to execute it.

\[det\]**http_location_by_id**(`+ID, -Location`)  
True when `Location` represents the HTTP path to which the handler with identifier `ID` is bound. Handler identifiers are deduced from the [http_handler/3](#http_handler/3) declaration as follows:

**Explicit id**  
If a term `id(ID)` appears in the option list of the handler, `ID` it is used and takes preference over using the predicate.

**Using the handler predicate**  
`ID` matches a handler if the predicate name matches `ID`. The `ID` may have a module qualification, e.g., `Module:Pred`

If the handler is declared with a pattern, e.g., `root(user/User)`, the location to access a particular *user* may be accessed using e.g., `user('Bob')`. The number of arguments to the compound term must match the number of variables in the path pattern.

A plain atom `ID` can be used to find a handler with a pattern. The returned location is the path up to the first variable, e.g., `/user/` in the example above.

User code is advised to use [http_link_to_id/3](#http_link_to_id/3) which can also add query parameters to the URL. This predicate is a helper for [http_link_to_id/3](#http_link_to_id/3).

Errors  
`existence_error(http_handler_id, Id)`.

See also  
[http_link_to_id/3](#http_link_to_id/3) and the `library(http/html_write)` construct `location_by_id(ID)` or its abbreviation `#(ID)`

**http_link_to_id**(`+HandleID, +Parameters, -HREF`)  
`HREF` is a link on the local server to a handler with given ID, passing the given `Parameters`. This predicate is typically used to formulate a `HREF` that resolves to a handler implementing a particular predicate. The code below provides a typical example. The predicate user_details/1 returns a page with details about a user from a given id. This predicate is registered as a handler. The DCG user_link//1 renders a link to a user, displaying the name and calling user_details/1 when clicked. Note that the location (`root(user_details)`) is irrelevant in this equation and HTTP locations can thus be moved freely without breaking this code fragment.

``` code
:- http_handler(root(user_details), user_details, []).

user_details(Request) :-
    http_parameters(Request,
                    [ user_id(ID)
                    ]),
    ...

user_link(ID) -->
    { user_name(ID, Name),
      http_link_to_id(user_details, [id(ID)], HREF)
    },
    html(a([class(user), href(HREF)], Name)).
```

[TABLE]

See also  
[http_location_by_id/2](#http_location_by_id/2) and [http_handler/3](#http_handler/3) for defining and specifying handler IDs.

\[det\]**http_reload_with_parameters**(`+Request, +Parameters, -HREF`)  
Create a request on the current handler with replaced search parameters.

\[det\]**http_reply_file**(`+FileSpec, +Options, +Request`)  
`Options` is a list of

**cache**(`+Boolean`)  
If `true` (default), handle If-modified-since and send modification time.

**mime_type**(`+Type`)  
Overrule mime-type guessing from the filename as provided by file_mime_type/2.

**static_gzip**(`+Boolean`)  
If `true` (default `false`) and, in addition to the plain file, there is a `.gz` file that is not older than the plain file and the client accepts `gzip` encoding, send the compressed file with `Transfer-encoding: gzip`.

**cached_gzip**(`+Boolean`)  
If `true` (default `false`) the system maintains cached gzipped files in a directory accessible using the file search path `http_gzip_cache` and serves these similar to the `static_gzip(true)` option. If the gzip file does not exist or is older than the input the file is recreated.

**unsafe**(`+Boolean`)  
If `false` (default), validate that `FileSpec` does not contain references to parent directories. E.g., specifications such as `www('../../etc/passwd')` are not allowed.

**headers**(`+List`)  
Provides additional reply-header fields, encoded as a list of *Field(Value)*.

If caching is not disabled, it processes the request headers `If-modified-since` and `Range`.

throws  
\- `http_reply(not_modified)`  
- `http_reply(file(MimeType, Path))`

\[det\]**http_safe_file**(`+FileSpec, +Options`)  
True if `FileSpec` is considered *safe*. If it is an atom, it cannot be absolute and cannot have references to parent directories. If it is of the form `alias(Sub)`, than Sub cannot have references to parent directories.

Errors  
\- instantiation_error  
- `permission_error(read, file, FileSpec)`

\[det\]**http_redirect**(`+How, +To, +Request`)  
Redirect to a new location. The argument order, using the `Request` as last argument, allows for calling this directly from the handler declaration:

``` code
:- http_handler(root(.),
                http_redirect(moved, myapp('index.html')),
                []).
```

|  |  |
|----|----|
| `How` | is one of `moved`, `moved_temporary` or `see_other` |
| `To` | is an atom, a aliased path as defined by [http_absolute_location/3](#http_absolute_location/3). or a term `location_by_id(Id)` or its abbreviations `#(Id)` or `#(Id)+Parameters`. If `To` is not absolute, it is resolved relative to the current location. |

\[det\]**http_404**(`+Options, +Request`)  
Reply using an "HTTP 404 not found" page. This handler is intended as fallback handler for *prefix* handlers. `Options` processed are:

**index**(`Location`)  
If there is no path-info, redirect the request to `Location` using [http_redirect/3](#http_redirect/3).

Errors  
`http_reply(not_found(Path))`

**http_switch_protocol**(`:Goal, +Options`)  
Send an `"HTTP 101 Switching Protocols"` reply. After sending the reply, the HTTP library calls `call(Goal, InStream, OutStream)`, where InStream and OutStream are the raw streams to the HTTP client. This allows the communication to continue using an an alternative protocol.

If `Goal` fails or throws an exception, the streams are closed by the server. Otherwise `Goal` is responsible for closing the streams. Note that `Goal` runs in the HTTP handler thread. Typically, the handler should be registered using the `spawn` option if [http_handler/3](#http_handler/3) or `Goal` must call thread_create/3 to allow the HTTP worker to return to the worker pool.

The streams use binary (octet) encoding and have their I/O timeout set to the server timeout (default 60 seconds). The predicate set_stream/2 can be used to change the encoding, change or cancel the timeout.

This predicate interacts with the server library by throwing an exception.

The following options are supported:

**header**(`+Headers`)  
Backward compatible. Use `headers(+Headers)`.

**headers**(`+Headers`)  
Additional headers send with the reply. Each header takes the form Name(Value).

### 3.3 library(http/http_dirindex): HTTP directory listings

To be done  
Provide more options (sorting, selecting columns, hiding files)

This module provides a simple API to generate an index for a physical directory. The index can be customised by overruling the dirindex.css CSS file and by defining additional rules for icons using the hook http:file_extension_icon/2.

\[det\]**http_reply_dirindex**(`+DirSpec, :Options, +Request`)  
Provide a directory listing for `Request`, assuming it is an index for the physical directrory Dir. If the request-path does not end with /, first return a moved (301 Moved Permanently) reply.

The calling conventions allows for direct calling from [http_handler/3](#http_handler/3).

\[det\]**directory_index**(`+Dir, :Options`)`//`  
Show index for a directory. `Options` processed:

**order_by**(`+Field`)  
Sort the files in the directory listing by `Field`. `Field` is one of `name` (default), `size` or `time`.

**order**(`+AscentDescent`)  
Sorting order. Default is `ascending`. The alternative is `descending`

**name**(`:RenderName`)  
DCG used to render a name in the table. The File is passed.

\[nondet,multifile\]http:**mime_type_icon**(`+MimeType, -IconName`)  
Multi-file hook predicate that can be used to associate icons to files listed by [http_reply_dirindex/3](#http_reply_dirindex/3). The actual icon file is located by `absolute_file_name(icons(IconName), Path, [])`.

See also  
serve_files_in_directory/2 serves the images.

### 3.4 library(http/http_files): Serve plain files from a hierarchy

See also  
[pwp_handler/2](#pwp_handler/2) provides similar facilities, where .pwp files can be used to add dynamic behaviour.

Although the SWI-Prolog Web Server is intended to serve documents that are computed dynamically, serving plain files is sometimes necessary. This small module combines the functionality of [http_reply_file/3](#http_reply_file/3) and [http_reply_dirindex/3](#http_reply_dirindex/3) to act as a simple web-server. Such a server can be created using the following code sample, which starts a server at port 8080 that serves files from the current directory (’.’). Note that the handler needs a `prefix` option to specify that it must handle all paths that begin with the registered location of the handler.

``` code
:- use_module(library(http/http_server)).
:- use_module(library(http/http_files)).

:- http_handler(root(.), http_reply_from_files('.', []), [prefix]).

:- initialization(http_server([port(8080)]), main).
```

**http_reply_from_files**(`+Dir, +Options, +Request`)  
HTTP handler that serves files from the directory `Dir`. This handler uses [http_reply_file/3](#http_reply_file/3) to reply plain files. If the request resolves to a directory, it uses the option `indexes` to locate an index file (see below) or uses [http_reply_dirindex/3](#http_reply_dirindex/3) to create a listing of the directory.

`Options`:

**indexes**(`+List`)  
`List` of files tried to find an index for a directory. The default is `['index.html']`.

**not_found**(`+Action`)  
`Action` defines what happens if the target file was not found. `Options`: `fail` makes the handler fail silently. `404` make the handler call [http_404/2](#http_404/2). Default is `fail`.

Note that this handler must be tagged as a `prefix` handler (see [http_handler/3](#http_handler/3) and module introduction). This also implies that it is possible to override more specific locations in the hierarchy using [http_handler/3](#http_handler/3) with a longer path-specifier.

When using [http_handler/3](#http_handler/3) to bind this predicate to an HTTP location, make sure it is bound to a location that ends in a `/`. When using [http:location/3](#http:location/3) to define symbolic names to HTTP locations this is written as

`:-` `http_handler(aliasname(.), http_reply_from_files(srcdir, []), [prefix])`.

|  |  |
|----|----|
| `Dir` | is either a directory or an path-specification as used by absolute_file_name/3. This option provides great flexibility in (re-)locating the physical files and allows merging the files of multiple physical locations into one web-hierarchy by using multiple user:file_search_path/2 clauses that define the same alias. |

See also  
The hookable predicate file_mime_type/2 is used to determine the `Content-type` from the file name.

### 3.5 library(http/http_session): HTTP Session management

This library defines session management based on HTTP cookies. Session management is enabled simply by loading this module. Details can be modified using [http_set_session_options/1](#http_set_session_options/1). By default, this module creates a session whenever a request is processes that is inside the hierarchy defined for session handling (see path option in [http_set_session_options/1](#http_set_session_options/1)). Automatic creation of a session can be stopped using the option `create(noauto)`. The predicate [http_open_session/2](#http_open_session/2) must be used to create a session if `noauto` is enabled. Sessions can be closed using [http_close_session/1](#http_close_session/1).

If a session is active, [http_in_session/1](#http_in_session/1) returns the current session and [http_session_assert/1](#http_session_assert/1) and friends maintain data about the session. If the session is reclaimed, all associated data is reclaimed too.

Begin and end of sessions can be monitored using `library(broadcast)`. The broadcasted messages are:

**http_session**(`begin(SessionID,Peer)`)  
Broadcasted if a session is started

**http_session**(`end(SessionId,Peer)`)  
Broadcasted if a session is ended. See [http_close_session/1](#http_close_session/1).

For example, the following calls `end_session(SessionId)` whenever a session terminates. Please note that sessions ends are not scheduled to happen at the actual timeout moment of the session. Instead, creating a new session scans the active list for timed-out sessions. This may change in future versions of this library.

``` code
:- listen(http_session(end(SessionId, Peer)),
          end_session(SessionId)).
```

\[det\]**http_set_session_options**(`+Options`)  
Set options for the session library. Provided options are:

**timeout**(`+Seconds`)  
Session timeout in seconds. Default is 600 (10 min). A timeout of `0` (zero) disables timeout.

**cookie**(`+Cookiekname`)  
Name to use for the cookie to identify the session. Default `swipl_session`.

**path**(`+Path`)  
`Path` to which the cookie is associated. Default is `/`. Cookies are only sent if the HTTP request path is a refinement of `Path`.

**route**(`+Route`)  
Set the route name. Default is the unqualified hostname. To cancel adding a route, use the empty atom. See route/1.

**enabled**(`+Boolean`)  
Enable/disable session management. Session management is enabled by default after loading this file.

**create**(`+Atom`)  
Defines when a session is created. This is one of `auto` (default), which creates a session if there is a request whose path matches the defined session path or `noauto`, in which cases sessions are only created by calling [http_open_session/2](#http_open_session/2) explicitly.

**proxy_enabled**(`+Boolean`)  
Enable/disable proxy session management. Proxy session management associates the *originating* IP address of the client to the session rather than the *proxy* IP address. Default is false.

**gc**(`+When`)  
`When` is one of `active`, which starts a thread that performs session cleanup at close to the moment of the timeout or `passive`, which runs session GC when a new session is created.

**samesite**(`+Restriction`)  
One of `none`, `lax` (default), or `strict` - The SameSite attribute prevents the CSRF vulnerability. strict has best security, but prevents links from external sites from operating properly. lax stops most CSRF attacks against REST endpoints but rarely interferes with legit image operations. `none` removes the samesite attribute entirely. **Caution: The value `none` exposes the entire site to CSRF attacks**.

**http_only**(`+Boolean`)  
If `true` (default `false`), add the `HttpOnly` property to the session cookie. This causes the browser to deny access from JavaScript.

**secure**(`+Boolean`)  
If `true`, (default `false`), add the `Secure` property to the session cookie. This causes the browser to report the cookie only over HTTPS connections.

**granularity**(`+Integer`)  
Granularity for updating that the session is active. Default is 60 (seconds). Smaller values lead to more precise session timeout at the cost of more database updates. This may notably a problem when using Redis.

In addition, extension libraries can define session_option/2 to make this predicate support more options. In particular, `library(http/http_redis_plugin)` defines the following additional options:

**redis_db**(`+DB`)  
Alias name of the redis database to access. See redis_server/3.

**redis_ro**(`+DB`)  
Alias name of the redis database for read-only access. See redis_server/3.

**redis_prefix**(`+Atom`)  
Prefix to use for all HTTP session related keys. Default is `'swipl:http:session'`

\[nondet\]**http_session_option**(`?Option`)  
True if `Option` is a current option of the session system.

\[semidet\]**session_setting**(`+SessionID, ?Setting`)  
Find setting for `SessionID`. It is possible to overrule some session settings using `http_session_set(Setting)`.

\[det\]**http_set_session**(`Setting`)  
\[det\]**http_set_session**(`SessionId, Setting`)  
Overrule a setting for the current or specified session. Currently, the only setting that can be overruled is `timeout`.

Errors  
`permission_error(set, http_session, Setting)` if setting a setting that is not supported on per-session basis.

\[det\]**http_session_id**(`-SessionId`)  
True if `SessionId` is an identifier for the current session.

|             |             |
|-------------|-------------|
| `SessionId` | is an atom. |

Errors  
`existence_error(http_session, _)`

See also  
[http_in_session/1](#http_in_session/1) for a version that fails if there is no session.

\[semidet\]**http_in_session**(`-SessionId`)  
True if `SessionId` is an identifier for the current session. The current session is extracted from `session(ID)` from the current HTTP request (see [http_current_request/1](#http_current_request/1)). The value is cached in a backtrackable global variable `http_session_id`. Using a backtrackable global variable is safe because continuous worker threads use a failure driven loop and spawned threads start without any global variables. This variable can be set from the commandline to fake running a goal from the commandline in the context of a session.

See also  
[http_session_id/1](#http_session_id/1)

\[det\]**http_open_session**(`-SessionID, +Options`)  
Establish a new session. This is normally used if the create option is set to `noauto`. `Options`:

**renew**(`+Boolean`)  
If `true` (default `false`) and the current request is part of a session, generate a new session-id. By default, this predicate returns the current session as obtained with [http_in_session/1](#http_in_session/1).

Errors  
`permission_error(open, http_session, CGI)` if this call is used after closing the CGI header.

See also  
\- [http_set_session_options/1](#http_set_session_options/1) to control the `create` option.  
- [http_close_session/1](#http_close_session/1) for closing the session.

\[det\]**http_session_asserta**(`+Data`)  
\[det\]**http_session_assert**(`+Data`)  
\[nondet\]**http_session_retract**(`?Data`)  
\[det\]**http_session_retractall**(`?Data`)  
Versions of assert/1, retract/1 and retractall/1 that associate data with the current HTTP session.

\[nondet\]**http_session_data**(`?Data`)  
True if `Data` is associated using [http_session_assert/1](#http_session_assert/1) to the current HTTP session.

Errors  
`existence_error(http_session,_)`

\[det\]**http_session_asserta**(`+Data, +SessionID`)  
\[det\]**http_session_assert**(`+Data, +SessionID`)  
\[nondet\]**http_session_retract**(`?Data, +SessionID`)  
\[det\]**http_session_retractall**(`@Data, +SessionID`)  
\[det\]**http_session_data**(`?Data, +SessionID`)  
Versions of assert/1, retract/1 and retractall/1 that associate data with an explicit HTTP session.

See also  
[http_current_session/2](#http_current_session/2).

\[nondet\]**http_current_session**(`?SessionID, ?Data`)  
Enumerate the current sessions and associated data. There are two *pseudo* data elements:

**idle**(`Seconds`)  
Session has been idle for `Seconds`.

**peer**(`Peer`)  
`Peer` of the connection.

\[det\]**http_close_session**(`+SessionID`)  
Closes an HTTP session. This predicate can be called from any thread to terminate a session. It uses the broadcast/1 service with the message below.

``` code
http_session(end(SessionId, Peer))
```

The broadcast is done **before** the session data is destroyed and the listen-handlers are executed in context of the session that is being closed. Here is an example that destroys a Prolog thread that is associated to a thread:

``` code
:- listen(http_session(end(SessionId, _Peer)),
          kill_session_thread(SessionID)).

kill_session_thread(SessionID) :-
        http_session_data(thread(ThreadID)),
        thread_signal(ThreadID, throw(session_closed)).
```

Succeed without any effect if `SessionID` does not refer to an active session.

If [http_close_session/1](#http_close_session/1) is called from a handler operating in the current session and the CGI stream is still in state `header`, this predicate emits a `Set-Cookie` to expire the cookie.

Errors  
`type_error(atom, SessionID)`

See also  
listen/2 for acting upon closed sessions

\[det\]**http_session_cookie**(`-Cookie`)  
Generate a random cookie that can be used by a browser to identify the current session. The cookie has the format XXXX-XXXX-XXXX-XXXX\[.\<route\>\], where XXXX are random hexadecimal numbers and \[.\<route\>\] is the optionally added routing information.

\[semidet,multifile\]**hooked**  
\[multifile\]**hook**(`+Goal`)  
These multifile predicates may be used to hook the data storage of this library. An example is implemented by `library(http/http_redis_plugin)`, storing all session data in a redis database.

### 3.6 library(http/http_cors): Enable CORS: Cross-Origin Resource Sharing

See also  
\- [http://en.wikipedia.org/wiki/Cross-site_scripting](http://en.wikipedia.org/wiki/Cross-site_scripting) for understanding Cross-site scripting.  
- [http://www.w3.org/TR/cors/](http://www.w3.org/TR/cors/) for understanding CORS

This small module allows for enabling Cross-Origin Resource Sharing (CORS) for a specific request. Typically, CORS is enabled for API services that you want to have usable from browser client code that is loaded from another domain. An example are the LOD and SPARQL services in ClioPatria.

Because CORS is a security risk (see references), it is disabled by default. It is enabled through the setting http:cors. The value of this setting is a list of domains that are allowed to access the service. Because \* is used as a wildcard match, the value \[\*\] allows access from anywhere.

Services for which CORS is relevant must call [cors_enable/0](#cors_enable/0) as part of the HTTP response, as shown below. Note that [cors_enable/0](#cors_enable/0) is a no-op if the setting http:cors is set to the empty list (`[]`).

``` code
my_handler(Request) :-
      ....,
      cors_enable,
      reply_json(Response, []).
```

If a site uses a *Preflight* `OPTIONS` request to find the server's capabilities and access politics, [cors_enable/2](#cors_enable/2) can be used to formulate an appropriate reply. For example:

``` code
my_handler(Request) :-
      option(method(options), Request), !,
      cors_enable(Request,
                  [ methods([get,post,delete])
                  ]),
      format('~n').                           % 200 with empty body
```

\[det\]**cors_enable**  
Emit the HTTP header `Access-Control-Allow-Origin` using domains from the setting http:cors. This this setting is `[]` (default), nothing is written. This predicate is typically used for replying to API HTTP-request (e.g., replies to an AJAX request that typically serve JSON or XML).

\[det\]**cors_enable**(`+Request, +Options`)  
CORS reply to a *Preflight* `OPTIONS` request. `Request` is the HTTP request. `Options` provides:

**methods**(`+List`)  
`List` of supported HTTP methods. The default is `GET`, only allowing for read requests.

**headers**(`+List`)  
`List` of headers the client asks for and we allow. The default is to simply echo what has been requested for.

Both methods and headers may use Prolog friendly syntax, e.g., `get` for a method and `content_type` for a header.

See also  
[http://www.html5rocks.com/en/tutorials/cors/](http://www.html5rocks.com/en/tutorials/cors/)

### 3.7 library(http/http_authenticate): Authenticate HTTP connections using 401 headers

This module provides the basics to validate an HTTP `Authorization` header. User and password information are read from a Unix/Apache compatible password file.

This library provides, in addition to the HTTP authentication, predicates to read and write password files.

**http_authenticate**(`+Type, +Request, -Fields`)  
True if `Request` contains the information to continue according to `Type`. `Type` identifies the required authentication technique:

**basic**(`+PasswordFile`)  
Use HTTP `Basic` authentication and verify the password from `PasswordFile`. `PasswordFile` is a file holding usernames and passwords in a format compatible to Unix and Apache. Each line is record with `:` separated fields. The first field is the username and the second the password *hash*. Password hashes are validated using crypt/2.

Successful authorization is cached for 60 seconds to avoid overhead of decoding and lookup of the user and password data.

[http_authenticate/3](#http_authenticate/3) just validates the header. If authorization is not provided the browser must be challenged, in response to which it normally opens a user-password dialogue. Example code realising this is below. The exception causes the HTTP wrapper code to generate an HTTP 401 reply.

``` code
(   http_authenticate(basic(passwd), Request, Fields)
->  true
;   throw(http_reply(authorise(basic, Realm)))
).
```

|  |  |
|----|----|
| `Fields` | is a list of fields from the password-file entry. The first element is the user. The hash is skipped. |

To be done  
Should we also cache failures to reduce the risk of DoS attacks?

\[semidet\]**http_authorization_data**(`+AuthorizeText, ?Data`)  
Decode the HTTP `Authorization` header. `Data` is a term

``` code
Method(User, Password)
```

where Method is the (downcased) authorization method (typically `basic`), User is an atom holding the user name and Password is a list of codes holding the password

\[nondet\]**http_current_user**(`+File, ?User, ?Fields`)  
True when `User` is present in the htpasswd file `File` and `Fields` provides the additional fields.

|  |  |
|----|----|
| `Fields` | are the fields from the password file `File`, converted using name/2, which means that numeric values are passed as numbers and other fields as atoms. The password hash is the first element of `Fields` and is a string. |

\[det\]**http_read_passwd_file**(`+Path, -Data`)  
Read a password file. `Data` is a list of terms of the format below, where User is an atom identifying the user, Hash is a string containing the salted password hash and Fields contain additional fields. The string value of each field is converted using name/2 to either a number or an atom.

``` code
passwd(User, Hash, Fields)
```

\[det\]**http_write_passwd_file**(`+File, +Data:list`)  
Write password data `Data` to `File`. `Data` is a list of entries as below. See [http_read_passwd_file/2](#http_read_passwd_file/2) for details.

``` code
passwd(User, Hash, Fields)
```

To be done  
Write to a new file and atomically replace the old one.

\[multifile\]http:**authenticate**(`+AuthData, +Request, -Fields`)  
Plugin for `library(http_dispatch)` to perform basic HTTP authentication.

This predicate throws `http_reply(authorise(basic, Realm))`.

|  |  |
|----|----|
| `AuthData` | must be a term `basic(File, Realm)` |
| `Request` | is the HTTP request |
| `Fields` | describes the authenticated user with the option `user(User)` and with the option `user_details(Fields)` if the password file contains additional fields after the user and password. |

### 3.8 library(http/http_digest): HTTP Digest authentication

See also  
[https://tools.ietf.org/html/rfc2617](https://tools.ietf.org/html/rfc2617)

This library implements HTTP *Digest Authentication* as per RFC2617. Unlike *Basic Authentication*, digest authentication is based on challenge-response and therefore does not need to send the password over the (insecure) connection. In addition, it provides a count mechanism that ensure that old credentials cannot be reused, which prevents attackers from using old credentials with a new request. Digest authentication have the following advantages and disadvantages:

- Advantages
  - Authentication without exchanging the password
  - No re-use of authentication data
- Disadvantages
  - An extra round trip is needed for the first authentication
  - Server-side storage of the password is the MD5 hash of the user, *realm* and password. As MD5 hashes are quick to compute, one needs strong passwords. This fixed algorithm also allows for *rainbow table* attacks, although their value is limited because you need to precompute the rainbow table for every server (*realm*) and user.
  - The connection is sensitive to man-in-the-middle attack, where the attacker can both change the request and response.
  - Both client and server need to keep an administration of issued *nonce* values and associated *nonce count* values.

And, of course, the connection itself remains insecure. Digest based authentication is a viable alternative if HTTPS is not a good option and security of the data itself is not an issue.

This library acts as plugin for `library(http/http_dispatch)`, where the registered handler ([http_handler/3](#http_handler/3)) can be given the option below to initiate digest authentication.

- `authentication(digest(PasswdFile, Realm))`

Above, `PasswdFile` is a file containing lines of the from below, where PasswordHash is computed using [http_digest_password_hash/4](#http_digest_password_hash/4). See also `library(http/http_authenticate)`, [http_read_passwd_file/2](#http_read_passwd_file/2) and [http_write_passwd_file/2](#http_write_passwd_file/2).

``` code
User ":" PasswordHash (":" Extra)*
```

This library also hooks into `library(http/http_open)` if the option `authorization(digest(User, Password))` is given.

**http_digest_challenge**(`+Realm, +Options`)`//`  
Generate the content for a 401 `WWW-Authenticate: Digest` header field.

\[det\]**http_parse_digest_challenge**(`+Challenge, -Fields`)  
Parse the value of an HTTP `WWW-Authenticate` header into a list of Name(Value) terms.

**http_digest_response**(`+Challenge, +User, +Password, -Reply, +Options`)  
Formulate a reply to a digest authentication request. `Options`:

**path**(`+Path`)  
The request URI send along with the authentication. Defaults to `/`

**method**(`+Method`)  
The HTTP method. Defaults to `'GET'`

**nc**(`+Integer`)  
The nonce-count as an integer. This is formatted as an 8 hex-digit string.

|  |  |
|----|----|
| `Challenge` | is a list Name(Value), normally from [http_parse_digest_challenge/2](#http_parse_digest_challenge/2). Must contain `realm` and `nonce`. Optionally contains `opaque`. |
| `User` | is the user we want to authenticated |
| `Password` | is the user's password |
| `Options` | provides additional options |

\[det\]**http_digest_password_hash**(`+User, +Realm, +Password, -Hash`)  
Compute the password hash for the HTTP password file. Note that the HTTP digest mechanism does allow us to use a seeded expensive arbitrary hash function. Instead, the hash is defined as the MD5 of the following components:

``` code
<user>:<realm>:<password>.
```

The inexpensive MD5 algorithm makes the hash sensitive to brute force attacks while the lack of seeding make the hashes sensitive for *rainbow table* attacks, although the value is somewhat limited because the *realm* and *user* are part of the hash.

\[multifile\]http:**authenticate**(`+Digest, +Request, -Fields`)  
Plugin for `library(http_dispatch)` to perform basic HTTP authentication. Note that we keep the authentication details cached to avoid a‘nonce-replay’error in the case that the application tries to verify multiple times.

This predicate throws `http_reply(authorise(digest(Digest)))`

|  |  |
|----|----|
| `Digest` | is a term `digest(File, Realm, Options)` |
| `Request` | is the HTTP request |
| `Fields` | describes the authenticated user with the option `user(User)` and with the option `user_details(Fields)` if the password file contains additional fields after the user and password. |

\[semidet,multifile\]http:**authenticate_client**(`+URL, +Action`)  
This hooks is called by [http_open/3](#http_open/3) with the following `Action` value:

**send_auth_header**(`+AuthData, +Out, +Options`)  
Called when sending the initial request. `AuthData` contains the value for the [http_open/3](#http_open/3) option `authorization(AuthData)` and `Out` is a stream on which to write additional HTTP headers.

**auth_reponse**(`+Headers, +OptionsIn, -Options`)  
Called if the server replies with a 401 code, challenging the client. Our implementation adds a `request_header(authorization=Digest)` header to `Options`, causing [http_open/3](#http_open/3) to retry the request with the additional option.

### 3.9 library(http/http_dyn_workers): Dynamically schedule HTTP workers.

Most code doesn't need to use this directly; instead use `library(http/http_server)`, which combines this library with the typical HTTP libraries that most servers need.

This module defines hooks into the HTTP framework to dynamically schedule worker threads. Dynamic scheduling relieves us from finding a good value for the size of the HTTP worker pool.

The decision to add a worker follows these rules:

- If the load average caused by the worker threads exceeds http:max_load, no worker is added.
- Wait for some time, depending on how close we are to the http:max_workers limit.
  - If the worker is still needed, add it.

The policy depends on three settings:

`http`**`:`**`max_workers`  
The maximum number of workers that will be created. Default is 100.

`http`**`:`**`worker_idle_limit`  
The number of seconds a dynamic worker waits for a new job. If no job arrives in time it terminates. Default is 10 seconds.

`http`**`:`**`max_load`  
Max load average created by **the HTTP server**, i.e. the amount of CPU time consumed per second. Default is 10.

&nbsp;

\[multifile\]http:**schedule_workers**(`+Dict`)  
Called if there is no immediately free worker to handle the incoming request. The request is forwarded to the thread `__http_scheduler` as the hook is called in time critical code.

#### 3.9.1 Providing Server-Sent Events (sse)

[Server-Sent Events](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events) allows for setting up a simple event stream from the server to the client. It can serve roles similar to *long polling* and *web sockets*, enabling the server to notify its clients on some event. *Long polling* uses a normal HTTP (usually) GET request that blocks for a long time on the server. The server finishes the request when it wants to notify the client or after some time (e.g., a minute) to avoid a timeout on the client or some proxy. After receiving an event or timeout, the client repeats the request. *Web sockets* *upgrade* a the socket used for a normal HTTP request to create a bi-directional open communication channel that exchanges encapsulated messages in both directions. *Server-Sent Events* open a normal HTTP channel over which the server can sent simple text messages using a format similar to the HTTP header: a sequence of *Name: Value* lines followed by two newlines. Unlike long polling, the request does not complete after a message.

Following the MDN documentation above, an SSE request can be served using the simple example below, which generates an event counting every minute. The handler is declared to process the request on a new thread and to disable the request timeout, as the response is intended to live for as long as the client stays connected. Note that this design uses one thread per client and therefore does not scale well to large numbers of subscribers.

``` code
:- use_module(library(http/sse)).

:- http_handler(root(events), events,
                [ spawn([]),
                  time_limit(infinite)
                ]).

events(_Request) :-
    sse_open,
    between(1, infinite, Min),
        sse_send(_{event: minute, data: Min}),
        sleep(60),
        fail.
```

The library `library(http/sse)` hides the wire format and the response headers needed to defeat HTTP intermediaries that buffer small responses (see sse_open/0 and sse_send/1). Of course, rather than sleep/1 to decide when to fire the next event this thread typically has to wait for events in the application. This can be achieved using thread_wait/2 or message queues.

### 3.10 Custom Error Pages

It is possible to create arbitrary error pages for responses generated when a http_reply term is thrown. Currently this is only supported for status 403 (*authentication required*). To do this, instead of throwing `http_reply(authorise(Term))` throw `http_reply(authorise(Term), [], Key)`, where `Key` is an arbitrary term relating to the page you want to generate. You must then also define a clause of the multifile predicate [http:status_page_hook/3](#http:status_page_hook/3):

**http:status_page_hook**(`+TermOrCode, +Context, -CustomHTML`)  
TermOrCode is either the first argument of the `http_reply` exception or the HTTP status code, i.e., the hook is called twice. New code should using the `Term`. Context is the third argument of the http_reply exception which was thrown, and CustomHTML is a list of HTML tokens. A page equivalent to the default page for 401 is generated by the example below.

``` code
:- multifile http:status_page_hook/3.

http:status_page_hook(authorise(Term), _Context, HTML) :-
    phrase(page([ title('401 Authorization Required')
                ],
                [ h1('Authorization Required'),
                  p(['This server could not verify that you ',
                     'are authorized to access the document ',
                     'requested.  Either you supplied the wrong ',
                     'credentials (e.g., bad password), or your ',
                     'browser doesn\'t understand how to supply ',
                     'the credentials required.'
                     ]),
                  \address
                ]),
           HTML).
```

### 3.11 library(http/http_openid): OpenID consumer and server library

This library implements the OpenID protocol ([http://openid.net/)](http://openid.net/)). OpenID is a protocol to share identities on the network. The protocol itself uses simple basic HTTP, adding reliability using digitally signed messages.

Steps, as seen from the *consumer* (or *relying partner*).

1.  Show login form, asking for `openid_identifier`
2.  Get HTML page from `openid_identifier` and lookup `<link rel="openid.server" href="server">`
3.  Associate to *server*
4.  Redirect browser (302) to server using mode `checkid_setup`, asking to validate the given OpenID.
5.  OpenID server redirects back, providing digitally signed conformation of the claimed identity.
6.  Validate signature and redirect to the target location.

A **consumer** (an application that allows OpenID login) typically uses this library through [openid_user/3](#openid_user/3). In addition, it must implement the hook http_openid:`openid_hook(trusted(OpenId, Server))` to define accepted OpenID servers. Typically, this hook is used to provide a white-list of acceptable servers. Note that accepting any OpenID server is possible, but anyone on the internet can setup a dummy OpenID server that simply grants and signs every request. Here is an example:

``` code
:- multifile http_openid:openid_hook/1.

http_openid:openid_hook(trusted(_, OpenIdServer)) :-
    (   trusted_server(OpenIdServer)
    ->  true
    ;   throw(http_reply(moved_temporary('/openid/trustedservers')))
    ).

trusted_server('http://www.myopenid.com/server').
```

By default, information who is logged on is maintained with the session using [http_session_assert/1](#http_session_assert/1) with the term `openid(Identity)`. The hooks login/logout/logged_in can be used to provide alternative administration of logged-in users (e.g., based on client-IP, using cookies, etc.).

To create a **server**, you must do four things: bind the handlers [openid_server/2](#openid_server/2) and [openid_grant/1](#openid_grant/1) to HTTP locations, provide a user-page for registered users and define the `grant(Request, Options)` hook to verify your users. An example server is provided in in \<plbase\>/`doc/packages/examples/demo_openid.pl`

\[multifile\]**openid_hook**(`+Action`)  
Call hook on the OpenID management library. Defined hooks are:

**login**(`+OpenID`)  
Consider `OpenID` logged in.

**logout**(`+OpenID`)  
Logout `OpenID`

**logged_in**(`?OpenID`)  
True if `OpenID` is logged in

**grant**(`+Request, +Options`)  
Server: Reply positive on OpenID

**trusted**(`+OpenID, +Server`)  
True if `Server` is a trusted `OpenID` server

**ax**(`Values`)  
Called if the server provided AX attributes

**x_parameter**(`+Server, -Name, -Value`)  
Called to find additional HTTP parameters to send with the OpenID verify request.

\[det\]**openid_login**(`+OpenID`)  
Associate the current HTTP session with `OpenID`. If another `OpenID` is already associated, this association is first removed.

\[det\]**openid_logout**(`+OpenID`)  
Remove the association of the current session with any `OpenID`

\[semidet\]**openid_logged_in**(`-OpenID`)  
True if session is associated with `OpenID`.

\[det\]**openid_user**(`+Request:http_request, -OpenID:url, +Options`)  
True if `OpenID` is a validated `OpenID` associated with the current session. The scenario for which this predicate is designed is to allow an HTTP handler that requires a valid login to use the transparent code below.

``` code
handler(Request) :-
      openid_user(Request, OpenID, []),
      ...
```

If the user is not yet logged on a sequence of redirects will follow:

1.  Show a page for login (default: page /openid/login), predicate reply_openid_login/1)
2.  By default, the `OpenID` login page is a form that is submitted to the `verify`, which calls [openid_verify/2](#openid_verify/2).
3.  [openid_verify/2](#openid_verify/2) does the following:
    - Find the `OpenID` claimed identity and server
    - Associate to the `OpenID` server
    - redirects to the `OpenID` server for validation
4.  The `OpenID` server will redirect here with the authentication information. This is handled by [openid_authenticate/4](#openid_authenticate/4).

`Options`:

**login_url**(`Login`)  
(Local) URL of page to enter `OpenID` information. Default is the handler for openid_login_page/1

See also  
[openid_authenticate/4](#openid_authenticate/4) produces errors if login is invalid or cancelled.

\[det\]**openid_login_form**(`+ReturnTo, +Options`)`//`  
Create the OpenID form. This exported as a separate DCG, allowing applications to redefine /openid/login and reuse this part of the page. `Options` processed:

**action**(`Action`)  
URL of action to call. Default is the handler calling openid_verify/1.

**buttons**(`+Buttons`)  
`Buttons` is a list of `img` structures where the `href` points to an OpenID 2.0 endpoint. These buttons are displayed below the OpenID URL field. Clicking the button sets the URL field and submits the form. Requires Javascript support.

If the `href` is *relative*, clicking it opens the given location after adding’openid.return_to’and‘stay’.

**show_stay**(`+Boolean`)  
If `true`, show a checkbox that allows the user to stay logged on.

**openid_verify**(`+Options, +Request`)  
Handle the initial login form presented to the user by the relying party (consumer). This predicate discovers the OpenID server, associates itself with this server and redirects the user's browser to the OpenID server, providing the extra openid.X name-value pairs. `Options` is, against the conventions, placed in front of the `Request` to allow for smooth cooperation with `http_dispatch.pl`. `Options` processes:

**return_to**(`+URL`)  
Specifies where the OpenID provider should return to. Normally, that is the current location.

**trust_root**(`+URL`)  
Specifies the `openid.trust_root` attribute. Defaults to the root of the current server (i.e., `http://host[.port]/`).

**realm**(`+URL`)  
Specifies the `openid.realm` attribute. Default is the `trust_root`.

**ax**(`+Spec`)  
`Request` the exchange of additional attributes from the identity provider. See http_ax_attributes/2 for details.

The OpenId server will redirect to the `openid.return_to` URL.

throws  
`http_reply(moved_temporary(Redirect))`

\[nondet\]**openid_server**(`?OpenIDLogin, ?OpenID, ?Server`)  
True if `OpenIDLogin` is the typed id for `OpenID` verified by `Server`.

|               |                                 |
|---------------|---------------------------------|
| `OpenIDLogin` | ID as typed by user (canonized) |
| `OpenID`      | ID as verified by server        |
| `Server`      | URL of the `OpenID` server      |

\[det\]**openid_current_url**(`+Request, -URL`)  
Find the public `URL` for `Request` that we can make available to our identity provider. This must be an absolute `URL` where we can be contacted. Before trying a configured version through [http_public_url/2](#http_public_url/2), we try to see whether the login message contains a referrer parameter or whether the browser provided one.

**openid_current_host**(`Request, Host, Port`)  
Find current location of the server.

deprecated  
New code should use [http_current_host/4](#http_current_host/4) with the option `global(true)`.

**ssl_verify**(`+SSL, +ProblemCert, +AllCerts, +FirstCert, +Error`)  
Accept all certificates. We do not care too much. Only the user cares s/he is not entering her credentials with a spoofed side. As we redirect, the browser will take care of this.

\[semidet\]**openid_authenticate**(`+Request, -Server:url, -OpenID:url, -ReturnTo:url`)  
Succeeds if `Request` comes from the `OpenID` server and confirms that User is a verified `OpenID` user. `ReturnTo` provides the URL to return to.

After [openid_verify/2](#openid_verify/2) has redirected the browser to the `OpenID` server, and the `OpenID` server did its magic, it redirects the browser back to this address. The work is fairly trivial. If `mode` is `cancel`, the OpenId server denied. If `id_res`, the OpenId server replied positive, but we must verify what the server told us by checking the HMAC-SHA signature.

This call fails silently if their is no `openid.mode` field in the request.

throws  
\- `openid(cancel)` if request was cancelled by the OpenId server  
- `openid(signature_mismatch)` if the HMAC signature check failed

**openid_server**(`+Options, +Request`)  
Realise the OpenID server. The protocol demands a POST request here.

**openid_grant**(`+Request`)  
Handle the reply from checkid_setup_server/3. If the reply is `yes`, check the authority (typically the password) and if all looks good redirect the browser to ReturnTo, adding the OpenID properties needed by the Relying Party to verify the login.

\[det\]**openid_associate**(`?URL, ?Handle, ?Assoc`)  
Calls [openid_associate/4](#openid_associate/4) as

``` code
openid_associate(URL, Handle, Assoc, []).
```

\[det\]**openid_associate**(`+URL, -Handle, -Assoc, +Options`)  
\[semidet\]**openid_associate**(`?URL, +Handle, -Assoc, +Options`)  
Associate with an open-id server. We first check for a still valid old association. If there is none or it is expired, we establish one and remember it. `Options`:

**ns**(`URL`)  
One of `http://specs.openid.net/auth/2.0` (default) or `http://openid.net/signon/1.1`.

To be done  
Should we store known associations permanently? Where?

### 3.12 Get parameters from HTML forms

The library `library(http/http_parameters)` provides two predicates to fetch HTTP request parameters as a type-checked list easily. The library transparently handles both GET and POST requests. It builds on top of the low-level request representation described in [section 3.13](#sec:3.13).

**http_parameters**(`+Request, ?Parameters`)  
The predicate is passes the `Request` as provided to the handler goal by [http_wrapper/5](#http_wrapper/5) as well as a partially instantiated lists describing the requested parameters and their types. Each parameter specification in `Parameters` is a term of the format `Name`(`-Value`, `+Options`) . `Options` is a list of option terms describing the type, default, etc. If no options are specified the parameter must be present and its value is returned in `Value` as an atom.

If a parameter is missing the exception `error(``existence_error(http_parameter, Name)``, _)` is thrown which. If the argument cannot be converted to the requested type, a `error(``existence_error(Type, Value)``, _)` is raised, where the error context indicates the HTTP parameter. If not caught, the server translates both errors into a `400 Bad request` HTTP message.

Options fall into three categories: those that handle presence of the parameter, those that guide conversion and restrict types and those that support automatic generation of documentation. First, the presence-options:

**default**(`Default`)  
If the named parameter is missing, `Value` is unified to `Default`.

**optional**(`true`)  
If the named parameter is missing, `Value` is left unbound and no error is generated.

**list**(`Type`)  
The same parameter may not appear or appear multiple times. If this option is present, `default` and `optional` are ignored and the value is returned as a list. Type checking options are processed on each value.

**zero_or_more**  
Deprecated. Use `list(Type)`.

The type and conversion options are given below. The type-language can be extended by providing clauses for the multifile hook http:convert_parameter/3.

**`;`**(`Type1, Type2`)  
Succeed if either `Type1` or `Type2` applies. It allows for checks such as `(nonneg;oneof([infinite]))` to specify an integer or a symbolic value.

**oneof**(`List`)  
Succeeds if the value is member of the given list.

**length `> N`**  
Succeeds if value is an atom of more than `N` characters.

**length `>= N`**  
Succeeds if value is an atom of more than or equal to `N` characters.

**length `< N`**  
Succeeds if value is an atom of less than `N` characters.

**length `=< N`**  
Succeeds if value is an atom of length less than or equal to `N` characters.

**atom**  
No-op. Allowed for consistency.

**string**  
Convert value to a string.

**between**(`+Low, +High`)  
Convert value to a number and if either `Low` or `High` is a float, force value to be a float. Then check that the value is in the given range, which includes the boundaries.

**boolean**  
Translate =true=, =yes=, =on= and’1’into =true=; =false=, =no=, =off= and’0’into =false= and raises an error otherwise.

**float**  
Convert value to a float. Integers are transformed into float. Throws a type-error otherwise.

**integer**  
Convert value to an integer. Throws a type-error otherwise.

**nonneg**  
Convert value to a non-negative integer. Throws a type-error of the value cannot be converted to an integer and a domain-error otherwise.

**number**  
Convert value to a number. Throws a type-error otherwise.

The last set of options is to support automatic generation of HTTP API documentation from the sources.^(4This facility is under development in ClioPatria; see `http_help.pl`).

**description**(`+Atom`)  
Description of the parameter in plain text.

**group**(`+Parameters, +Options`)  
Define a logical group of parameters. `Parameters` are processed as normal. `Options` may include a description of the group. Groups can be nested.

Below is an example

``` code
reply(Request) :-
        http_parameters(Request,
                        [ title(Title, [ optional(true) ]),
                          name(Name,   [ length >= 2 ]),
                          age(Age,     [ between(0, 150) ])
                        ]),
        ...
```

Same as `http_parameters(Request, Parameters,[])`

**http_parameters**(`+Request, ?Parameters, +Options`)  
In addition to [http_parameters/2](#http_parameters/2), the following options are defined.

**form_data**(`-Data`)  
Return the entire set of provided `Name`=`Value` pairs from the GET or POST request. All values are returned as atoms.

**attribute_declarations**(`:Goal`)  
If a parameter specification lacks the parameter options, call `call(Goal, +ParamName, -Options)` to find the options. Intended to share declarations over many calls to [http_parameters/3](#http_parameters/3). Using this construct the above can be written as below.

``` code
reply(Request) :-
        http_parameters(Request,
                        [ title(Title),
                          name(Name),
                          age(Age)
                        ],
                        [ attribute_declarations(param)
                        ]),
        ...

param(title, [optional(true)]).
param(name,  [length >= 2 ]).
param(age,   [integer]).
```

### 3.13 Request format

The body-code (see [section 3.1](#sec:3.1)) is driven by a `Request`. This request is generated from [http_read_request/2](#http_read_request/2) defined in `library(http/http_header)`.

**http_read_request**(`+Stream, -Request`)  
Reads an HTTP request from `Stream` and unify `Request` with the parsed request. `Request` is a list of `Name``(Value)` elements. It provides a number of predefined elements for the result of parsing the first line of the request, followed by the additional request parameters. The predefined fields are:

**host**(`Host`)  
If the request contains `Host: ``Host`, Host is unified with the host-name. If `Host` is of the format \<`host`\>:\<`port`\> `Host` only describes \<`host`\> and a field `port(Port)` where `Port` is an integer is added.

**input**(`Stream`)  
The `Stream` is passed along, allowing to read more data or requests from the same stream. This field is always present.

**method**(`Method`)  
`Method` is the HTTP *method* represented as a lower-case atom (i.e., `delete`, `get`, `head`, `options`, `patch`, `post`, `put`, `trace`). This field is present if the header has been parsed successfully.

**path**(`Path`)  
Path associated to the request. This field is always present.

**peer**(`Peer`)  
`Peer` is a term `ip(A,B,C,D)` containing the IP address of the contacting host.

**port**(`Port`)  
Port requested. See `host` for details.

**request_uri**(`RequestURI`)  
This is the untranslated string that follows the method in the request header. It is used to construct the path and search fields of the `Request`. It is provided because reconstructing this string from the path and search fields may yield a different value due to different usage of percent encoding.

**search**(`ListOfNameValue`)  
Search-specification of URI. This is the part after the `?`, normally used to transfer data from HTML forms that use the HTTP GET method. In the URL it consists of a www-form-encoded list of `Name`=`Value` pairs. This is mapped to a list of Prolog `Name`=`Value` terms with decoded names and values. This field is only present if the location contains a search-specification.

The URL specification does not *demand* the query part to be of the form *name=value*. If the field is syntactically incorrect, ListOfNameValue is bound the the empty list (\[\]).

**http_version**(`Major-Minor`)  
If the first line contains the `HTTP/``Major`.`Minor` version indicator this element indicate the HTTP version of the peer. Otherwise this field is not present.

**cookie**(`ListOfNameValue`)  
If the header contains a `Cookie` line, the value of the cookie is broken down in `Name`=`Value` pairs, where the `Name` is the lowercase version of the cookie name as used for the HTTP fields.

**set_cookie**(`set_cookie(Name, Value, Options)`)  
If the header contains a `SetCookie` line, the cookie field is broken down into the `Name` of the cookie, the `Value` and a list of `Name`=`Value` pairs for additional options such as `expire`, `path`, `domain` or `secure`.

If the first line of the request is tagged with `HTTP/``Major`.`Minor`, [http_read_request/2](#http_read_request/2) reads all input upto the first blank line. This header consists of `Name`:`Value` fields. Each such field appears as a term `Name``(Value)` in the `Request`, where `Name` is canonicalised for use with Prolog. Canonisation implies that the `Name` is converted to lower case and all occurrences of the `-` are replaced by `_`. The value for the `Content-length` fields is translated into an integer.

Here is an example:

``` code
?- http_read_request(user_input, X).
|: GET /mydb?class=person HTTP/1.0
|: Host: gollem
|:
X = [ input(user),
      method(get),
      search([ class = person
             ]),
      path('/mydb'),
      http_version(1-0),
      host(gollem)
    ].
```

#### 3.13.1 Handling POST requests

Where the HTTP `GET` operation is intended to get a document, using a `path` and possibly some additional search information, the `POST` operation is intended to hand potentially large amounts of data to the server for processing.

The `Request` parameter above contains the term `method(post)`. The data posted is left on the input stream that is available through the term `input(Stream)` from the `Request` header. This data can be read using [http_read_data/3](#http_read_data/3) from the HTTP client library. Here is a demo implementation simply returning the parsed posted data as plain text (assuming pp/1 pretty-prints the data).

``` code
reply(Request) :-
        member(method(post), Request), !,
        http_read_data(Request, Data, []),
        format('Content-type: text/plain~n~n', []),
        pp(Data).
```

If the POST is initiated from a browser, content-type is generally either `application/x-www-form-urlencoded` or `multipart/form-data`.

### 3.14 Running the server

The functionality of the server should be defined in one Prolog file (of course this file is allowed to load other files). Depending on the wanted server setup this‘body’is wrapped into a small Prolog file combining the body with the appropriate server interface. There are three supported server-setups. For most applications we advise the multi-threaded server. Examples of this server architecture are the [PlDoc](http://www.swi-prolog.org/packages/pldoc.html) documentation system and the [SeRQL](http://www.swi-prolog.org/packages/SeRQL/) Semantic Web server infrastructure.

All the server setups may be wrapped in a *reverse proxy* to make them available from the public web-server as described in [section 3.14.7](#sec:3.14.7).

- *Using `library(thread_httpd)` for a multi-threaded server*  
  This server exploits the multi-threaded version of SWI-Prolog, running the users body code parallel from a pool of worker threads. As it avoids the state engine and copying required in the event-driven server it is generally faster and capable to handle multiple requests concurrently.

  This server is harder to debug due to the involved threading, although the GUI tracer provides reasonable support for multi-threaded applications using the tspy/1 command. It can provide fast communication to multiple clients and can be used for more demanding servers.

- *Using `library(inetd_httpd)` for server-per-client*  
  In this setup the Unix **inetd** user-daemon is used to initialise a server for each connection. This approach is especially suitable for servers that have a limited startup-time. In this setup a crashing client does not influence other requests.

  This server is very hard to debug as the server is not connected to the user environment. It provides a robust implementation for servers that can be started quickly.

#### 3.14.1 Common server interface options

All the server interfaces provide `http_server(:Goal, +Options)` to create the server. The list of options differ, but the servers share common options:

**port**(`?Port`)  
Specify the port to listen to for stand-alone servers. `Port` is either an integer or unbound. If unbound, it is unified to the selected free port.

#### 3.14.2 Multi-threaded Prolog

The `library(http/thread_httpd.pl)` provides the infrastructure to manage multiple clients using a pool of *worker-threads*. This realises a popular server design, also seen in Java Tomcat and Microsoft .NET. As a single persistent server process maintains communication to all clients startup time is not an important issue and the server can easily maintain state-information for all clients.

In addition to the functionality provided by the inetd server, the threaded server can also be used to realise an HTTPS server exploiting the `library(ssl)` library. See option `ssl(+SSLOptions)` below.

**http_server**(`:Goal, +Options`)  
Create the server. `Options` must provide the `port(?Port)` option to specify the port the server should listen to. If `Port` is unbound an arbitrary free port is selected and `Port` is unified to this port-number. The server consists of a small Prolog thread accepting new connection on `Port` and dispatching these to a pool of workers. Defined `Options` are:

**port**(`?Address`)  
Address to bind to. `Address` is either a port (integer) or a term `Host`:`Port`. The port may be a variable, causing the system to select a free port and unify the variable with the selected port. See also tcp_bind/2.

**workers**(`+N`)  
Defines the number of worker threads in the pool. Default is to use `five` workers. Choosing the optimal value for best performance is a difficult task depending on the number of CPUs in your system and how much resources are required for processing a request. Too high numbers makes your system switch too often between threads or even swap if there is not enough memory to keep all threads in memory, while a too low number causes clients to wait unnecessary for other clients to complete. See also [http_workers/2](#http_workers/2).

**timeout**(`+SecondsOrInfinite`)  
Determines the maximum period of inactivity handling a request. If no data arrives within the specified time since the last data arrived, the connection raises an exception, and the worker discards the client and returns to the pool-queue for a new client. If it is `infinite`, a worker may wait forever on a client that doesn't complete its request. Default is 60 seconds.

**keep_alive_timeout**(`+SecondsOrInfinite`)  
Maximum time to wait for new activity on *Keep-Alive* connections. Choosing the correct value for this parameter is hard. Disabling Keep-Alive is bad for performance if the clients request multiple documents for a single page. This may ---for example-- be caused by HTML frames, HTML pages with images, associated CSS files, etc. Keeping a connection open in the threaded model however prevents the thread servicing the client servicing other clients. The default is 2 seconds.

**local**(`+KBytes`)  
Size of the local-stack for the workers. Default is taken from the commandline option.

**global**(`+KBytes`)  
Size of the global-stack for the workers. Default is taken from the commandline option.

**trail**(`+KBytes`)  
Size of the trail-stack for the workers. Default is taken from the commandline option.

**ssl**(`+SSLOptions`)  
Use SSL (Secure Socket Layer) rather than plain TCP/IP. A server created this way is accessed using the `https://` protocol. SSL allows for encrypted communication to avoid others from tapping the wire as well as improved authentication of client and server. The `SSLOptions` option list is passed to ssl_context/3. The port option of the main option list is forwarded to the SSL layer. See the `library(ssl)` library for details.

**http_server_property**(`?Port, ?Property`)  
True if `Property` is a property of the HTTP server running at `Port`. Defined properties are:

**goal**(`:Goal`)  
Goal used to start the server. This is often [http_dispatch/1](#http_dispatch/1).

**scheme**(`-Scheme`)  
Scheme is one of `http` or `https`.

**start_time**(`-Time`)  
Time-stamp when the server was created. See format_time/3 for creating a human-readable representation.

**http_workers**(`+Port, ?Workers`)  
Query or manipulate the number of workers of the server identified by `Port`. If `Workers` is unbound it is unified with the number of running servers. If it is an integer greater than the current size of the worker pool new workers are created with the same specification as the running workers. If the number is less than the current size of the worker pool, this predicate inserts a number of‘quit’requests in the queue, discarding the excess workers as they finish their jobs (i.e. no worker is abandoned while serving a client).

This can be used to tune the number of workers for performance. Another possible application is to reduce the pool to one worker to facilitate easier debugging.

**http_add_worker**(`+Port, +Options`)  
Add a new worker to the HTTP server for port `Port`. `Options` overrule the default queue options. The following additional options are processed:

**max_idle_time**(`+Seconds`)  
The created worker will automatically terminate if there is no new work within Seconds.

**http_stop_server**(`+Port, +Options`)  
Stop the HTTP server at Port. Halting a server is done *gracefully*, which means that requests being processed are not abandoned. The `Options` list is for future refinements of this predicate such as a forced immediate abort of the server, but is currently ignored.

**http_current_worker**(`?Port, ?ThreadID`)  
True if `ThreadID` is the identifier of a Prolog thread serving `Port`. This predicate is motivated to allow for the use of arbitrary interaction with the worker thread for development and statistics.

**http_spawn**(`:Goal, +Spec`)  
Continue handling this request in a new thread running `Goal`. After [http_spawn/2](#http_spawn/2), the worker returns to the pool to process new requests. In its simplest form, `Spec` is the name of a thread pool as defined by thread_pool_create/3. Alternatively it is an option list, whose options are passed to thread_create_in_pool/4 if `Spec` contains `pool(Pool)` or to thread_create/3 of the pool option is not present. If the dispatch module is used (see [section 3.2](#sec:3.2)), spawning is normally specified as an option to the [http_handler/3](#http_handler/3) registration.

We recommend the use of thread pools. They allow registration of a set of threads using common characteristics, specify how many can be active and what to do if all threads are active. A typical application may define a small pool of threads with large stacks for computation intensive tasks, and a large pool of threads with small stacks to serve media. The declaration could be the one below, allowing for max 3 concurrent solvers and a maximum backlog of 5 and 30 tasks creating image thumbnails.

``` code
:- use_module(library(thread_pool)).

:- thread_pool_create(compute, 3,
                      [ local(20000), global(100000), trail(50000),
                        backlog(5)
                      ]).
:- thread_pool_create(media, 30,
                      [ local(100), global(100), trail(100),
                        backlog(100)
                      ]).

:- http_handler('/solve',     solve,     [spawn(compute)]).
:- http_handler('/thumbnail', thumbnail, [spawn(media)]).
```

#### 3.14.3 library(http/http_unix_daemon): Run SWI-Prolog HTTP server as a Unix system daemon

See also  
The file \<swi-home\>/doc/packages/examples/http/linux-init-script provides a /etc/init.d script for controlling a server as a normal Unix service.

To be done  
Cleanup issues wrt. loading and initialization of xpce.

This module provides the logic that is needed to integrate a process into the Unix service (daemon) architecture. It deals with the following aspects, all of which may be used/ignored and configured using commandline options:

- Select the `port(s)` to be used by the server
- Run the startup of the process as root to perform privileged tasks and the server itself as unprivileged user, for example to open ports below 1000.
- Fork and detach from the controlling terminal
- Handle console and debug output using a file and/or the syslog daemon.
- Manage a *pid file*

The typical use scenario is to write a file that loads the following components:

1.  The application code, including http handlers (see [http_handler/3](#http_handler/3)).
2.  This library

In the code below, `?- [load].` loads the remainder of the webserver code. This is often a sequence of use_module/1 directives.

``` code
:- use_module(library(http/http_unix_daemon)).

:- [load].
```

The program entry point is [http_daemon/0](#http_daemon/0), declared using initialization/2. This may be overruled using a new declaration after loading this library. The new entry point will typically call [http_daemon/1](#http_daemon/1) to start the server in a preconfigured way.

``` code
:- use_module(library(http/http_unix_daemon)).
:- initialization(run, main).

run :-
    ...
    http_daemon(Options).
```

Now, the server may be started using the command below. See [http_daemon/0](#http_daemon/0) for supported options.

``` code
% [sudo] swipl mainfile.pl [option ...]
```

Below are some examples. Our first example is completely silent, running on port 80 as user `www`.

``` code
% swipl mainfile.pl --user=www --pidfile=/var/run/http.pid
```

Our second example logs HTTP interaction with the syslog daemon for debugging purposes. Note that the argument to `--debug`= is a Prolog term and must often be escaped to avoid misinterpretation by the Unix shell. The debug option can be repeated to log multiple debug topics.

``` code
% swipl mainfile.pl --user=www --pidfile=/var/run/http.pid \
        --debug='http(request)' --syslog=http
```

**Broadcasting** The library uses broadcast/1 to allow hooking certain events:

**http**(`pre_server_start`)  
Run *after* *fork*, just before starting the HTTP server. Can be used to load additional files or perform additional initialisation, such as starting additional threads. Recall that it is not possible to start threads *before* forking.

**http**(`post_server_start`)  
Run *after* starting the HTTP server.

&nbsp;

**http_daemon**  
Start the HTTP server as a daemon process. This predicate processes the commandline arguments below. Commandline arguments that specify servers are processed in the order they appear using the following schema:

1.  Arguments that act as default for all servers.

2.  `--http=Spec` or `--https=Spec` is followed by arguments for that server until the next `--http=Spec` or `--https=Spec` or the end of the options.

3.  If no `--http=Spec` or `--https=Spec` appears, one HTTP server is created from the specified parameters.

    Examples:

    ``` code
    --workers=10 --http --https
    --http=8080 --https=8443
    --http=localhost:8080 --workers=1 --https=8443 --workers=25
    ```

**--port=Port**  
Start HTTP server at Port. It requires root permission and the option `--user=User` to open ports below 1000. The default port is 80. If `--https` is used, the default port is 443.

**--ip=IP**  
Only listen to the given IP address. Typically used as `--ip=localhost` to restrict access to connections from *localhost* if the server itself is behind an (Apache) proxy server running on the same host.

**--debug=Topic**  
Enable debugging Topic. See debug/3.

**--syslog=Ident**  
Write debug messages to the syslog daemon using Ident

**--user=User**  
When started as root to open a port below 1000, this option must be provided to switch to the target user for operating the server. The following actions are performed as root, i.e., *before* switching to User:

- open the `socket(s)`
- write the pidfile
- setup syslog interaction
- Read the certificate, key and password file (`--pwfile=File`)

**--group=Group**  
May be used in addition to `--user`. If omitted, the login group of the target user is used.

**--pidfile=File**  
Write the PID of the daemon process to File.

**--output=File**  
Send output of the process to File. By default, all Prolog console output is discarded.

**--fork\[=Bool\]**  
If given as `--no-fork` or `--fork=false`, the process runs in the foreground.

**--http\[=(Bool`|`Port`|`BindTo:Port)\]**  
Create a plain HTTP server. If the argument is missing or `true`, create at the specified or default address. Else use the given port and interface. Thus, `--http` creates a server at port 80, `--http=8080` creates one at port 8080 and `--http=localhost:8080` creates one at port 8080 that is only accessible from `localhost`.

**--https\[=(Bool`|`Port`|`BindTo:Port)\]**  
As `--http`, but creates an HTTPS server. Use `--certfile`, `--keyfile`, `-pwfile`, `--password` and `--cipherlist` to configure SSL for this server.

**--certfile=File**  
The server certificate for HTTPS.

**--keyfile=File**  
The server private key for HTTPS.

**--pwfile=File**  
File holding the password for accessing the private key. This is preferred over using `--password=PW` as it allows using file protection to avoid leaking the password. The file is read *before* the server drops privileges when started with the `--user` option.

**--password=PW**  
The password for accessing the private key. See also‘--pwfile\`.

**--cipherlist=Ciphers**  
One or more cipher strings separated by colons. See the OpenSSL documentation for more information. Starting with SWI-Prolog 7.5.11, the default value is always a set of ciphers that was considered secure enough to prevent all critical attacks at the time of the SWI-Prolog release.

**--interactive\[=Bool\]**  
If `true` (default `false`) implies `--no-fork` and presents the Prolog toplevel after starting the server.

**--gtrace=\[Bool\]**  
Use the debugger to trace [http_daemon/1](#http_daemon/1).

**--sighup=Action**  
Action to perform on `kill -HUP <pid>`. Default is `reload` (running make/0). Alternative is `quit`, stopping the server.

Other options are converted by argv_options/3 and passed to [http_server/1](#http_server/1). For example, this allows for:

**--workers=Count**  
Set the number of workers for the multi-threaded server.

[http_daemon/0](#http_daemon/0) is defined as below. The start code for a specific server can use this as a starting point, for example for specifying defaults or additional options. This uses *guided* options processing from argv_options/3 from `library(main)`. The option definitions are available as [http_opt_type/3](#http_opt_type/3), [http_opt_help/2](#http_opt_help/2) and [http_opt_meta/2](#http_opt_meta/2)

``` code
http_daemon :-
    current_prolog_flag(argv, Argv),
    argv_options(Argv, _RestArgv, Options),
    http_daemon(Options).
```

See also  
[http_daemon/1](#http_daemon/1)

**http_opt_type**(`?Flag, ?Option, ?Type`)  
**http_opt_help**(`?Option, ?Help`)  
**http_opt_meta**(`?Option, ?Meta`)  
Allow reusing http option processing

**http_daemon**(`+Options`)  
Start the HTTP server as a daemon process. This predicate processes a Prolog option list. It is normally called from [http_daemon/0](#http_daemon/0), which derives the option list from the command line arguments.

Error handling depends on whether or not `interactive(true)` is in effect. If so, the error is printed before entering the toplevel. In non-interactive mode this predicate calls `halt(1)`.

\[semidet,multifile\]**http_certificate_hook**(`+CertFile, +KeyFile, -Password`)  
Hook called before starting the server if the --https option is used. This hook may be used to create or refresh the certificate. If the hook binds `Password` to a string, this string will be used to decrypt the server private key as if the --password=`Password` option was given.

\[semidet,multifile\]**http_server_hook**(`+Options`)  
Hook that is called to start the HTTP server. This hook must be compatible to `http_server(Handler, Options)`. The default is provided by start_server/1.

\[multi,multifile\]http:**sni_options**(`-HostName, -SSLOptions`)  
Hook to provide Server Name Indication (SNI) for TLS servers. When starting an HTTPS server, all solutions of this predicate are collected and a suitable sni_hook/1 is defined for ssl_context/3 to use different contexts depending on the host name of the client request. This hook is executed *before* privileges are dropped.

#### 3.14.4 From (Unix) inetd

All modern Unix systems handle a large number of the services they run through the super-server *inetd* or one of its descendants (*xinetd*, *systemd* etc.) Such a program reads a configuration file (for example `/etc/inetd.conf`) and opens server-sockets on all ports defined in this file. As a request comes in it accepts it and starts the associated server such that standard I/O is performed through the socket. This approach has several advantages:

- *Simplification of servers*  
  Servers don't have to know about sockets and -operations.
- *Centralised authorisation*  
  Using *tcpwrappers* and similar tools, simple and effective firewalling of all services can be realised.
- *Automatic start and monitor*  
  The inetd automatically starts the server‘just-in-time’and starts additional servers or restarts a crashed server according to its configuration.

The very small generic script for handling inetd based connections is in `inetd_httpd`, defining [http_server/1](#http_server/1):

**http_server**(`:Goal`)  
Initialises and runs [http_wrapper/5](#http_wrapper/5) in a loop until failure or end-of-file. This server does not support the `Port` option as the port is specified with the **inetd** configuration. The only supported option is `After`.

Here is the example from `demo_inetd`

``` code
#!/usr/bin/pl -t main -q -f
:- use_module(demo_body).
:- use_module(inetd_httpd).

main :-
        http_server(reply).
```

With the above file installed in `/home/jan/plhttp/demo_inetd`, the following line in `/etc/inetd` enables the server at port 4001 guarded by *tcpwrappers*. After modifying inetd, send the daemon the `HUP` signal to make it reload its configuration. For more information, please check **inetd.conf**(5).

``` code
4001 stream tcp nowait nobody /usr/sbin/tcpd /home/jan/plhttp/demo_inetd
```

#### 3.14.5 MS-Windows

There are rumours that *inetd* has been ported to Windows.

#### 3.14.6 As CGI script

To be done.

#### 3.14.7 Using a reverse proxy

There are several options for public deployment of a web service. The main decision is whether to run it on a standard port (port 80 for HTTP, port 443 for HTTPS) or a non-standard port such as for example 8000 or 8080. Using a standard port below 1000 requires root access to the machine, and prevents other web services from using the same port. On the other hand, using a non-standard port may cause problems with intermediate proxy- and/or firewall policies that may block the port when you try to access the service from some networks. In both cases, you can either use a physical or a virtual machine running ---for example--- under [VMWARE](http://www.vmware.com) or [XEN](http://www.cl.cam.ac.uk/research/srg/netos/xen/) to host the service. Using a dedicated (physical or virtual) machine to host a service isolates security threats. Isolation can also be achieved using a Unix *chroot* environment, which is however not a security feature.

To make several different web services reachable on the same (either standard or non-standard) port, you can use a so-called *reverse proxy*. A reverse proxy uses rules to relay requests to other web services that use their own dedicated ports. This approach has several advantages:

- We can run the service on a non-standard port, but still access it (via the proxy) on a standard port, just as for a dedicated machine. We do not need a separate machine though: We only need to configure the reverse proxy to relay requests to the intended target servers.
- As the main web server is doing the front-line service, the Prolog server is normally protected from malformed HTTP requests that could result in denial of service or otherwise compromise the server. In addition, the main web server can transparently provide encodings such as compression to the outside world.

Proxy technology can be combined with isolation methods such as dedicated machines, virtual machines and chroot jails. The proxy can also provide load balancing.

**Setting up an Apache reverse proxy**

The Apache reverse proxy setup is really simple. Ensure the modules `proxy` and `proxy_http` are loaded. Then add two simple rules to the server configuration. Below is an example that makes a PlDoc server on port 4000 available from the main Apache server at port 80.

``` code
ProxyPass        /pldoc/ http://localhost:4000/pldoc/
ProxyPassReverse /pldoc/ http://localhost:4000/pldoc/
```

Apache rewrites the HTTP headers passing by, but using the above rules it does not examine the content. This implies that URLs embedded in the (HTML) content must use relative addressing. If the locations on the public and Prolog server are the same (as in the example above) it is allowed to use absolute locations. I.e. `/pldoc/search` is ok, but `http://myhost.com:4000/pldoc/search` is *not*. If the locations on the server differ, locations must be relative (i.e. not start with `/`.

This problem can also be solved using the contributed Apache module `proxy_html` that can be instructed to rewrite URLs embedded in HTML documents. In our experience, this is not trouble-free as URLs can appear in many places in generated documents. JavaScript can create URLs on the fly, which makes rewriting virtually impossible.

### 3.15 The wrapper library

The body is called by the module `library(http/http_wrapper.pl)`. This module realises the communication between the I/O streams and the body described in [section 3.1](#sec:3.1). The interface is realised by [http_wrapper/5](#http_wrapper/5):

**http_wrapper**(`:Goal, +In, +Out, -Connection, +Options`)  
Handle an HTTP request where `In` is an input stream from the client, `Out` is an output stream to the client and `Goal` defines the goal realising the body. `Connection` is unified to `’Keep-alive’` if both ends of the connection want to continue the connection or `close` if either side wishes to close the connection.

This predicate reads an HTTP request-header from `In`, redirects current output to a memory file and then runs `call(Goal, Request)`, watching for exceptions and failure. If `Goal` executes successfully it generates a complete reply from the created output. Otherwise it generates an HTTP server error with additional context information derived from the exception.

[http_wrapper/5](#http_wrapper/5) supports the following options:

**request**(`-Request`)  
Return the executed request to the caller.

**peer**(`+Peer`)  
Add peer(Peer) to the request header handed to `Goal`. The format of `Peer` is defined by tcp_accept/3 from the clib package.

**http:request_expansion**(`+RequestIn, -RequestOut`)  
This *multifile* hook predicate is called just before the goal that produces the body, while the output is already redirected to collect the reply. If it succeeds it must return a valid modified request. It is allowed to throw exceptions as defined in [section 3.1.1](#sec:3.1.1). It is intended for operations such as mapping paths, deny access for certain requests or manage cookies. If it writes output, these must be HTTP header fields that are added *before* header fields written by the body. The example below is from the session management library (see [section 3.5](#sec:3.5)) sets a cookie.

``` code
        ...,
        format('Set-Cookie: ~w=~w; path=~w~n', [Cookie, SessionID, Path]),
        ...,
```

**http_current_request**(`-Request`)  
Get access to the currently executing request. `Request` is the same as handed to `Goal` of [http_wrapper/5](#http_wrapper/5) *after* applying rewrite rules as defined by [http:request_expansion/2](#http:request_expansion/2). Raises an existence error if there is no request in progress.

**http_relative_path**(`+AbsPath, -RelPath`)  
Convert an absolute path (without host, fragment or search) into a path relative to the current page, defined as the path component from the current request (see [http_current_request/1](#http_current_request/1)). This call is intended to create reusable components returning relative paths for easier support of reverse proxies.

If ---for whatever reason--- the conversion is not possible it simply unifies `RelPath` to `AbsPath`.

### 3.16 library(http/http_host): Obtain public server location

This library finds the public address of the running server. This can be used to construct URLs that are visible from anywhere on the internet. This module was introduced to deal with OpenID, where a request is redirected to the OpenID server, which in turn redirects to our server (see `http_openid.pl`).

The address is established from the settings `http:public_host` and `http:public_port` if provided. Otherwise it is deduced from the request.

\[det\]**http_public_url**(`+Request, -URL`)  
True when `URL` is an absolute `URL` for the current request. Typically, the login page should redirect to this `URL` to avoid losing the session.

\[det\]**http_public_host_url**(`+Request, -URL`)  
True when `URL` is the public `URL` at which this server can be contacted. This value is not easy to obtain. See [http_public_host/4](#http_public_host/4) for the hardest part: find the host and port.

\[det\]**http_public_host**(`?Request, -Hostname, -Port, +Options`)  
Current global host and port of the HTTP server. This is the basis to form absolute address, which we need for redirection based interaction such as the OpenID protocol. `Options` are:

**global**(`+Bool`)  
If `true` (default `false`), try to replace a local hostname by a world-wide accessible name.

This predicate performs the following steps to find the host and port:

1.  Use the settings `http:public_host` and `http:public_port`
2.  Use `X-Forwarded-Host` header, which applies if this server runs behind a proxy.
3.  Use the `Host` header, which applies for HTTP 1.1 if we are contacted directly.
4.  Use gethostname/1 to find the host and http_current_server/2 to find the port.

|  |  |
|----|----|
| `Request` | is the current request. If it is left unbound, and the request is needed, it is obtained with [http_current_request/1](#http_current_request/1). |

\[det\]**http_current_host**(`?Request, -Hostname, -Port, +Options`)  
deprecated  
Use [http_public_host/4](#http_public_host/4) (same semantics)

### 3.17 library(http/http_log): HTTP Logging module

Simple module for logging HTTP requests to a file. Logging is enabled by loading this file and ensure the setting http:logfile is not the empty atom. The default file for writing the log is `httpd.log`. See `library(settings)` for details.

The level of logging can be modified using the multifile predicate http_log:nolog/1 to hide HTTP request fields from the logfile and http_log:password_field/1 to hide passwords from HTTP search specifications (e.g. `/topsecret?password=secret`).

\[semidet\]**http_log_stream**(`-Stream`)  
True when `Stream` is a stream to the opened HTTP log file. Opens the log file in `append` mode if the file is not yet open. The log file is determined from the setting `http:logfile`. If this setting is set to the empty atom (” ), this predicate fails.

If a file error is encountered, this is reported using print_message/2, after which this predicate silently fails. Opening is retried every minute when a new message arrives.

Before opening the log file, the message `http_log_open(Term)` is broadcasted. This message allows for creating the directory, renaming, deleting or truncating an existing log file.

\[det\]**http_log_close**(`+Reason`)  
If there is a currently open HTTP logfile, close it after adding a term `server(Reason, Time)`. to the logfile. This call is intended for cooperation with the Unix logrotate facility using the following schema:

- Move logfile (the HTTP server keeps writing to the moved file)
- Inform the server using an HTTP request that calls [http_log_close/1](#http_log_close/1)
- Compress the moved logfile

author  
Suggested by Jacco van Ossenbruggen

\[det\]**http_log**(`+Format, +Args`)  
Write message from `Format` and `Args` to log-stream. See format/2 for details. Succeed without side effects if logging is not enabled.

\[semidet,multifile\]**password_field**(`+Field`)  
Multifile predicate that can be defined to hide passwords from the logfile.

\[multifile\]**nolog**(`+HTTPField`)  
Multifile predicate that can be defined to hide request parameters from the request logfile.

\[semidet,multifile\]**nolog_post_content_type**(`+Type`)  
Multifile hook called with the `Content-type` header. If the hook succeeds, the POST data is not logged. For example, to stop logging anything but application/json messages:

``` code
:- multifile http_log:nolog_post_content_type/1.

http_log:nolog_post_content_type(Type) :-
   Type \= (application/json).
```

|        |                            |
|--------|----------------------------|
| `Type` | is a term MainType/SubType |

\[det\]**post_data_encoded**(`?Bytes:string, ?Encoded:string`)  
Encode the POST body for inclusion into the HTTP log file. The POST data is (in/de)flated using zopen/3 and base64 encoded using base64//1. The encoding makes long text messages shorter and keeps readable logfiles if binary data is posted.

\[det\]**http_logrotate**(`+Options`)  
Rotate the available log files. Note that there are two ways to deal with the rotation of log files:

1.  Use the OS log rotation facility. In that case the OS must (1) move the logfile and (2) have something calling [http_log_close/1](#http_log_close/1) to close the (moved) file and make this server create a new one on the next log message. If `library(http/http_unix_daemon)` is used, closing is achieved by sending SIGHUP or SIGUSR1 to the process.
2.  Call this predicate at scheduled intervals. This can be achieved by calling [http_schedule_logrotate/2](#http_schedule_logrotate/2) in the context of `library(http/http_unix_daemon)` which schedules the maintenance actions.

`Options`:

**min_size**(`+Bytes`)  
Do not rotate if the log file is smaller than `Bytes`. The default is 1Mbytes.

**keep_logs**(`+Count`)  
Number of rotated log files to keep (default 10)

**compress_logs**(`+Format`)  
Compress the log files to the given format. Default `gzip` and this is currently the only supported compressor.

**background**(`+Boolean`)  
If `true`, rotate the log files in the background.

**http_schedule_logrotate**(`When, Options`)  
Schedule log rotation based on maintenance broadcasts. `When` is one of:

**daily**(`Hour:Min`)  
Run each day at `Hour`:`Min`. `Min` is rounded to a multitude of 5.

**weekly**(`Day, Hour:Min`)  
Run at the given `Day` and Time each week. `Day` is either a number 1..7 (1 is Monday) or a weekday name or abbreviation.

**monthly**(`DayOfTheMonth, Hour:Min`)  
Run each month at the given Day (1..31). Note that not all months have all days.

This must be used with a timer that broadcasts a `maintenance(_,_)` message (see broadcast/1). Such a timer is part of `library(http/http_unix_daemon)`.

### 3.18 library(http/http_server_health): HTTP Server health statistics

This module defines an HTTP handler for `/health`. The handler returns a JSON document with elementary health statistics on the running instance. The location can be changed using [http_handler/3](#http_handler/3). Keys may be added using additional clauses for [health/2](#health/2) or hidden using [hide/1](#hide/1).

This library defines an HTTP handler and defines two multifile predicates ([health/2](#health/2) and [hide/1](#hide/1)) to control the information presented.

**server_health**(`+Request`)  
HTTP handler that replies with the overall health of the server. Returns a JSON object from all solutions of [health/2](#health/2).

Processes an optional parameter `fields` to specify the fields that should be returned. The fields content is "," or white space delimited.

\[nondet,multifile\]**health**(`-Key, -Value`)  
Multifile extensible. True when `Key`/`Value` can be reported as a health statistics. Keys may be added by adding clauses to this multifile predicate. Keys may be filtered using [hide/1](#hide/1). Predefined `Key` values are:

**up**  
Defined to be `true`.

**epoch**  
Starting time of the server in seconds after Jan 1, 1970 UTC.

**cpu_time**  
Total process CPU usage in seconds.

**threads**  
Number of active threads

**workers**  
Number of HTTP *worker* threads.

**requests**  
Number of HTTP requests processed.

**bytes_sent**  
Number of bytes send in reply to HTTP requests.

**open_files**  
Number of open file streams. This includes physical files as well as sockets (except for Windows). On Linux we count the file handles in `/proc/self/fd`. Otherwise we use stream_property/2 with the `file_no(Fd)` property.

**loadavg**  
An array holding the load average over the last \[1,5,15\] minutes. This key is only supported on Linux machines. It is based on `/proc/loadavg`

**heap**  
When compiled with TCMalloc, this provides two properties:

`inuse`**`:`**`Bytes`  
Total amount of in-use memory in bytes

`size`**`:`**`Bytes`  
Same as `inuse`, but including the TCMalloc overhead and (thus) memory that has been freed and is not (yet) reused.

|         |                                                       |
|---------|-------------------------------------------------------|
| `Key`   | is the name of the JSON key. Must be an atom          |
| `Value` | is the Prolog representation for a JSON (dict) value. |

\[nondet,multifile\]**hide**(`?Key`)  
Multifile hook. If true for a specific `Key`, hide this statistics from the output. This may be used to hide keys that are considered a security risk.

### 3.19 Debugging HTTP servers

The library `library(http/http_error)` defines a hook that decorates uncaught exceptions with a stack-trace. This will generate a *500 internal server error* document with a stack-trace. To enable this feature, simply load this library. Please do note that providing error information to the user simplifies the job of a hacker trying to compromise your server. It is therefore not recommended to load this file by default.

The example program `calc.pl` has the error handler loaded which can be triggered by forcing a divide-by-zero in the calculator.

### 3.20 library(http/http_header): Handling HTTP headers

The library `library(http/http_header)` provides primitives for parsing and composing HTTP headers. Its functionality is normally hidden by the other parts of the HTTP server and client libraries.

\[det\]**http_read_request**(`+FdIn:stream, -Request`)  
Read an HTTP request-header from `FdIn` and return the broken-down request fields as +Name(+Value) pairs in a list. `Request` is unified to `end_of_file` if `FdIn` is at the end of input.

**http_read_reply_header**(`+FdIn, -Reply`)  
Read the HTTP reply header. Throws an exception if the current input does not contain a valid reply header.

\[det\]**http_reply**(`+Data, +Out:stream`)  
\[det\]**http_reply**(`+Data, +Out:stream, +HdrExtra`)  
\[det\]**http_reply**(`+Data, +Out:stream, +HdrExtra, -Code`)  
\[det\]**http_reply**(`+Data, +Out:stream, +HdrExtra, +Context, -Code`)  
\[det\]**http_reply**(`+Data, +Out:stream, +HdrExtra, +Context, +Request, -Code`)  
Compose a complete HTTP reply from the term `Data` using additional headers from `HdrExtra` to the output stream `Out`. ExtraHeader is a list of Field(Value). `Data` is one of:

**html**(`HTML`)  
`HTML` tokens as produced by [html//1](#html//1) from `html_write.pl`

**file**(`+MimeType, +FileName`)  
Reply content of `FileName` using `MimeType`

**file**(`+MimeType, +FileName, +Range`)  
Reply partial content of `FileName` with given `MimeType`

**tmp_file**(`+MimeType, +FileName`)  
Same as `file`, but do not include modification time

**bytes**(`+MimeType, +Bytes`)  
Send a sequence of `Bytes` with the indicated `MimeType`. `Bytes` is either a string of character codes 0..255 or list of integers in the range 0..255. `Out`-of-bound codes result in a representation error exception.

**stream**(`+In, +Len`)  
Reply content of stream.

**cgi_stream**(`+In, +Len`)  
Reply content of stream, which should start with an HTTP header, followed by a blank line. This is the typical output from a CGI script.

**`Status`**  
HTTP status report as defined by [http_status_reply/4](#http_status_reply/4).

|  |  |
|----|----|
| `HdrExtra` | provides additional reply-header fields, encoded as Name(Value). It can also contain a field `content_length(-Len)` to *retrieve* the value of the Content-length header that is replied. |
| `Code` | is the numeric HTTP status code sent |

To be done  
Complete documentation

\[det\]**http_status_reply**(`+Status, +Out, +HdrExtra, -Code`)  
\[det\]**http_status_reply**(`+Status, +Out, +HdrExtra, +Context, -Code`)  
\[det\]**http_status_reply**(`+Status, +Out, +HdrExtra, +Context, +Request, -Code`)  
Emit HTML non-200 status reports. Such requests are always sent as UTF-8 documents.

`Status` can be one of the following:

**authorise**(`Method`)  
Challenge authorization. `Method` is one of

- `basic(Realm)`
- `digest(Digest)`

**authorise**(`basic, Realm`)  
Same as `authorise(basic(Realm))`. Deprecated.

**bad_request**(`ErrorTerm`)  
**busy**  
**created**(`Location`)  
**forbidden**(`Url`)  
**moved**(`To`)  
**moved_temporary**(`To`)  
**no_content**  
**not_acceptable**(`WhyHtml`)  
**not_found**(`Path`)  
**method_not_allowed**(`Method, Path`)  
**not_modified**  
**resource_error**(`ErrorTerm`)  
**see_other**(`To`)  
**switching_protocols**(`Goal, Options`)  
**server_error**(`ErrorTerm`)  
**unavailable**(`WhyHtml`)  

\[semidet,multifile\]http:**serialize_reply**(`+Reply, -Body`)  
Multifile hook to serialize the result of http:status_reply/3 into a term

**body**(`Type, Encoding, Content`)  
In this term, `Type` is the media type, `Encoding` is the required wire encoding and `Content` a string representing the content.

**http_join_headers**(`+Default, +Header, -Out`)  
Append headers from `Default` to `Header` if they are not already part of it.

**http_update_encoding**(`+HeaderIn, -Encoding, -HeaderOut`)  
Allow for rewrite of the header, adjusting the encoding. We distinguish three options. If the user announces‘text’, we always use UTF-8 encoding. If the user announces charset=utf-8 we use UTF-8 and otherwise we use octet (raw) encoding. Alternatively we could dynamically choose for ASCII, ISO-Latin-1 or UTF-8.

\[semidet,multifile\]http:**mime_type_encoding**(`+MimeType, -Encoding`)  
`Encoding` is the (default) character encoding for `MimeType`. This is used for setting the encoding for HTTP replies after the user calls `format('Content-type: <MIME type>~n')`. This hook is called before mime_type_encoding/2. This default defines `utf8` for JSON and Turtle derived `application/` MIME types.

**http_update_connection**(`+CGIHeader, +Request, -Connection, -Header`)  
Merge keep-alive information from `Request` and `CGIHeader` into `Header`.

**http_update_transfer**(`+Request, +CGIHeader, -Transfer, -Header`)  
Decide on the transfer encoding from the `Request` and the CGI header. The behaviour depends on the setting http:chunked_transfer. If `never`, even explicit requests are ignored. If `on_request`, chunked encoding is used if requested through the CGI header and allowed by the client. If `if_possible`, chunked encoding is used whenever the client allows for it, which is interpreted as the client supporting HTTP 1.1 or higher.

Chunked encoding is more space efficient and allows the client to start processing partial results. The drawback is that errors lead to incomplete pages instead of a nicely formatted complete page.

\[det\]**http_post_data**(`+Data, +Out:stream, +HdrExtra`)  
Send data on behalf on an HTTP POST request. This predicate is normally called by [http_post/4](#http_post/4) from `http_client.pl` to send the POST data to the server. `Data` is one of:

- `html(+Tokens)` Result of [html//1](#html//1) from `html_write.pl`

- `json(+Term)` Posting a JSON query and processing the JSON reply (or any other reply understood by [http_read_data/3](#http_read_data/3)) is simple as `http_post(URL, json(Term), Reply, [])`, where Term is a JSON term as described in `json.pl` and reply is of the same format if the server replies with JSON, when using module `:- use_module(library(http/http_json))`. Note that the module is used in both http server and http client, see `library(http/http_json)`.

- `xml(+Term)` Post the result of xml_write/3 using the Mime-type `text/xml`

- `xml(+Type, +Term)` Post the result of xml_write/3 using the given Mime-type and an empty option list to xml_write/3.

- `xml(+Type, +Term, +Options)` Post the result of xml_write/3 using the given Mime-type and option list for xml_write/3.

- `file(+File)` Send contents of a file. Mime-type is determined by file_mime_type/2.

- `file(+Type, +File)` Send file with content of indicated mime-type.

- `memory_file(+Type, +Handle)` Similar to `file(+Type, +File)`, but using a memory file instead of a real file. See new_memory_file/1.

- `codes(+Codes)` As `codes(text/plain, Codes)`.

- `codes(+Type, +Codes)` Send Codes using the indicated MIME-type.

- `bytes(+Type, +Bytes)` Send Bytes using the indicated MIME-type. Bytes is either a string of character codes 0..255 or list of integers in the range 0..255. `Out`-of-bound codes result in a representation error exception.

- `atom(+Atom)` As `atom(text/plain, Atom)`.

- `atom(+Type, +Atom)` Send Atom using the indicated MIME-type.

- `string(+String)`

- `string(+Type, +String)` Similar to `atom(Atom)` and `atom(Type,Atom)`, accepting a SWI-Prolog string.

- `cgi_stream(+Stream, +Len)` Read the input from Stream which, like CGI data starts with a partial HTTP header. The fields of this header are merged with the provided `HdrExtra` fields. The first Len characters of Stream are used.

- `form(+ListOfParameter)` Send data of the MIME type application/x-www-form-urlencoded as produced by browsers issuing a POST request from an HTML form. ListOfParameter is a list of Name=Value or Name(Value).

- `form_data(+ListOfData)` Send data of the MIME type `multipart/form-data` as produced by browsers issuing a POST request from an HTML form using enctype `multipart/form-data`. ListOfData is the same as for the List alternative described below. Below is an example. Repository, etc. are atoms providing the value, while the last argument provides a value from a file.

  ``` code
  ...,
  http_post([ protocol(http),
              host(Host),
              port(Port),
              path(ActionPath)
            ],
            form_data([ repository = Repository,
                        dataFormat = DataFormat,
                        baseURI    = BaseURI,
                        verifyData = Verify,
                        data       = file(File)
                      ]),
            _Reply,
            []),
  ...,
  ```

- List If the argument is a plain list, it is sent using the MIME type multipart/mixed and packed using [mime_pack/3](#mime_pack/3). See [mime_pack/3](#mime_pack/3) for details on the argument format.

\[det\]**http_reply_header**(`+Out:stream, +What, +HdrExtra`)  
Create a reply header using reply_header//3 and send it to Stream.

\[semidet\]**http_parse_header_value**(`+Field, +Value, -Prolog`)  
Translate `Value` in a meaningful `Prolog` term. `Field` denotes the HTTP request field for which we do the translation. Supported fields are:

**content_length**  
Converted into an integer

**status**  
Converted into an integer

**cookie**  
Converted into a list with Name=`Value` by cookies//1.

**set_cookie**  
Converted into a term `set_cookie(Name, Value, Options)`. Options is a list consisting of Name=`Value` or a single atom (e.g., `secure`)

**host**  
Converted to HostName:Port if applicable.

**range**  
Converted into `bytes(From, To)`, where From is an integer and To is either an integer or the atom `end`.

**accept**  
Parsed to a list of media descriptions. Each media is a term `media(Type, TypeParams, Quality, AcceptExts)`. The list is sorted according to preference.

**content_disposition**  
Parsed into `disposition(Name, Attributes)`, where Attributes is a list of Name=`Value` pairs.

**content_type**  
Parsed into `media(Type/SubType, Attributes)`, where Attributes is a list of Name=`Value` pairs.

**expires**  
Parsed into a time stamp using [http_timestamp/2](#http_timestamp/2).

As some fields are already parsed in the `Request`, this predicate is a no-op when called on an already parsed field.

|  |  |
|----|----|
| `Value` | is either an atom, a list of codes or an already parsed header value. |

\[det\]**http_timestamp**(`?Time:timestamp, ?Text:atom`)  
Convert between a SWI-Prolog time stamp and a string in HTTP format (RFC1123). When parsing, it accepts RFC1123, RFC1036 and ASCTIME formats. See parse_time/3.

Errors  
`syntax_error(http_timestamp(Text))` if the string cannot be parsed.

\[det\]**http_read_header**(`+Fd, -Header`)  
Read Name: Value lines from FD until an empty line is encountered. Field-name are converted to Prolog conventions (all lower, \_ instead of -): Content-Type: text/html `-->` `content_type(text/html)`

\[det\]**http_parse_header**(`+Text:codes, -Header:list`)  
`Header` is a list of Name(Value)-terms representing the structure of the HTTP header in `Text`.

Errors  
`domain_error(http_request_line, Line)`

\[det,multifile\]http:**//**(`http_address`)  
HTML-rule that emits the location of the HTTP server. This hook is called from address//0 to customise the server address. The server address is emitted on non-200-ok replies.

\[semidet,multifile\]http:**status_page**(`+Status, +Context, -HTMLTokens`)  
Hook called by [http_status_reply/4](#http_status_reply/4) and [http_status_reply/5](#http_status_reply/5) that allows for emitting custom error pages for the following HTTP page types:

- 201 - `created(Location)`
- 301 - `moved(To)`
- 302 - `moved_temporary(To)`
- 303 - `see_other(To)`
- 400 - `bad_request(ErrorTerm)`
- 401 - `authorise(AuthMethod)`
- 403 - `forbidden(URL)`
- 404 - `not_found(URL)`
- 405 - `method_not_allowed(Method,URL)`
- 406 - `not_acceptable(Why)`
- 500 - `server_error(ErrorTerm)`
- 503 - `unavailable(Why)`

The hook is tried twice, first using the status term, e.g., `not_found(URL)` and than with the code, e.g. `404`. The second call is deprecated and only exists for compatibility.

|  |  |
|----|----|
| `Context` | is the 4th argument of [http_status_reply/5](#http_status_reply/5), which is invoked after raising an exception of the format `http_reply(Status, HeaderExtra, Context)`. The default context is `[]` (the empty list). |
| `HTMLTokens` | is a list of tokens as produced by [html//1](#html//1). It is passed to [print_html/2](#print_html/2). |

### 3.21 The `library(http/html_write)` library

Producing output for the web in the form of an HTML document is a requirement for many Prolog programs. Just using format/2 is not satisfactory as it leads to poorly readable programs generating poor HTML. This library is based on using DCG rules.

The `library(http/html_write)` structures the generation of HTML from a program. It is an extensible library, providing a *DCG* framework for generating legal HTML under (Prolog) program control. It is especially useful for the generation of structured pages (e.g. tables) from Prolog data structures.

The normal way to use this library is through the DCG html//1. This non-terminal provides the central translation from a structured term with embedded calls to additional translation rules to a list of atoms that can then be printed using [print_html/\[1,2\]](#print_html/1).

**html**(`:Spec`)`//`  
The DCG non-terminal html//1 is the main predicate of this library. It translates the specification for an HTML page into a list of atoms that can be written to a stream using [print_html/\[1,2\]](#print_html/1). The expansion rules of this predicate may be extended by defining the multifile DCG html_write:expand//1. `Spec` is either a single specification or a list of single specifications. Using nested lists is not allowed to avoid ambiguity caused by the atom `[]`

- *Atomic data*  
  Atomic data is quoted using html_quoted//1.

- *`Fmt` - `Args`*  
  `Fmt` and `Args` are used as format-specification and argument list to format/3. The result is quoted and added to the output list.

- *`\``List`*  
  Escape sequence to add atoms directly to the output list. This can be used to embed external HTML code or emit script output. `List` is a list of the following terms:
  - *`Fmt` - `Args`*  
    `Fmt` and `Args` are used as format-specification and argument list to format/3. The result is added to the output list.
  - *`Atomic`*  
    Atomic values are added directly to the output list.

- *`\``Term`*  
  Invoke the non-terminal `Term` in the calling module. This is the common mechanism to realise abstraction and modularisation in generating HTML.

- *`Module`:`Term`*  
  Invoke the non-terminal \<`Module`\>:\<`Term`\>. This is similar to `\``Term` but allows for invoking grammar rules in external packages.

- *&(Entity)*  
  Emit `&<``Entity``>;` or `&#<``Entity``>;` if `Entity` is an integer. SWI-Prolog atoms and strings are represented as Unicode. Explicit use of this construct is rarely needed because code-points that are not supported by the output encoding are automatically converted into character-entities.

- *`Tag(Content)`*  
  Emit HTML element `Tag` using `Content` and no attributes. `Content` is handed to html//1. See [section 3.21.4](#sec:3.21.4) for details on the automatically generated layout.

- *`Tag(Attributes, Content)`*  
  Emit HTML element `Tag` using `Attributes` and `Content`. `Attributes` is either a single attribute of a list of attributes. Each attributes is of the format `Name(Value)` or `Name`=`Value`. `Value` is the atomic attribute value but allows for a limited functional notation:

  - *`A` + `B`*  
    Concatenation of `A` and `B`

  - *`Format`-`Arguments`*  
    Use format/3 and emit the result as quoted value.

  - *`encode(Atom)`*  
    Use uri_encoded/3 to create a valid URL query component.

  - *`location_by_id(ID)`*  
    HTTP location of the HTTP handler with given ID. See [http_location_by_id/2](#http_location_by_id/2).

  - *`#``(ID)`*  
    Abbreviated for for `location_by_id(ID)`.

  - *`A` + `List`*  
    `List` is handled as a URL‘search’component. The list members are terms of the format `Name` = `Value` or `Name(Value)`. Values are encoded as in the encode option described above.

  - *`List`*  
    Emit SGML multi-valued attributes (e.g., `NAMES`). Each value in list is separated by a space. This is particularly useful for setting multiple `class` attributes on an element. For example:

    ``` code
            ...
            span(class([c1,c2]), ...),
    ```

  The example below generates a URL that references the predicate set_lang/1 in the application with given parameters. The [http_handler/3](#http_handler/3) declaration binds `/setlang` to the predicate set_lang/1 for which we provide a very simple implementation. The code between `...` is part of an HTML page showing the English flag which, when pressed, calls `set_lang(Request)` where `Request` contains the search parameter `lang` = `en`. Note that the HTTP location (path) `/setlang` can be moved without affecting this code.

  ``` code
  :- http_handler('/setlang', set_lang, []).

  set_lang(Request) :-
          http_parameters(Request,
                          [ lang(Lang, [])
                          ]),
          http_session_retractall(lang(_)),
          http_session_assert(lang(Lang)),
          reply_html_page(title('Switched language'),
                          p(['Switch language to ', Lang])).

          ...
          html(a(href(location_by_id(set_lang) + [lang(en)]),
                 img(src('/www/images/flags/en.png')))),
          ...
  ```

**page**(`:HeadContent, :BodyContent`)`//`  
The DCG non-terminal page//2 generated a complete page, including the SGML `DOCTYPE` declaration. `HeadContent` are elements to be placed in the `head` element and `BodyContent` are elements to be placed in the `body` element.

To achieve common style (background, page header and footer), it is possible to define DCG non-terminals head//1 and/or body//1. Non-terminal page//1 checks for the definition of these non-terminals in the module it is called from as well as in the `user` module. If no definition is found, it creates a head with only the `HeadContent` (note that the `title` is obligatory) and a `body` with `bgcolor` set to `white` and the provided `BodyContent`.

Note that further customisation is easily achieved using html//1 directly as page//2 is (besides handling the hooks) defined as:

``` code
page(Head, Body) -->
        html([ \['<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 4.0//EN">\n'],
               html([ head(Head),
                      body(bgcolor(white), Body)
                    ])
             ]).
```

**page**(`:Contents`)`//`  
This version of the page/\[1,2\] only gives you the SGML `DOCTYPE` and the `HTML` element. `Contents` is used to generate both the head and body of the page.

**html_begin**(`+Begin`)`//`  
Just open the given element. `Begin` is either an atom or a compound term, In the latter case the arguments are used as arguments to the begin-tag. Some examples:

``` code
        html_begin(table)
        html_begin(table(border(2), align(center)))
```

This predicate provides an alternative to using the `\``Command` syntax in the html//1 specification. The following two fragments are the same. The preferred solution depends on your preferences as well as whether the specification is generated or entered by the programmer.

``` code
table(Rows) -->
        html(table([border(1), align(center), width('80%')],
                   [ \table_header,
                     \table_rows(Rows)
                   ])).

% or

table(Rows) -->
        html_begin(table(border(1), align(center), width('80%'))),
        table_header,
        table_rows,
        html_end(table).
```

**html_end**(`+End`)`//`  
End an element. See html_begin/1 for details.

#### 3.21.1 Emitting HTML documents

The non-terminal html//1 translates a specification into a list of atoms and layout instructions. Currently the layout instructions are terms of the format `nl(N)`, requesting at least `N` newlines. Multiple consecutive `nl(1)` terms are combined to an atom containing the maximum of the requested number of newline characters.

To simplify handing the data to a client or storing it into a file, the following predicates are available from this library:

**reply_html_page**(`:Head, :Body`)  
Same as `reply_html_page(default, Head, Body)`.

**reply_html_page**(`+Style, :Head, :Body`)  
Writes an HTML page preceded by an HTTP header as required by `library(http_wrapper)` (CGI-style). Here is a simple typical example:

``` code
reply(Request) :-
        reply_html_page(title('Welcome'),
                        [ h1('Welcome'),
                          p('Welcome to our ...')
                        ]).
```

The header and footer of the page can be hooked using the grammar-rules user:head//2 and user:body//2. The first argument passed to these hooks is the `Style` argument of [reply_html_page/3](#reply_html_page/3) and the second is the 2nd (for head//2) or 3rd (for body//2) argument of [reply_html_page/3](#reply_html_page/3). These hooks can be used to restyle the page, typically by embedding the real body content in a `div`. E.g., the following code provides a menu on top of each page of that is identified using the style *myapp*.

``` code
:- multifile
        user:body//2.

user:body(myapp, Body) -->
        html(body([ div(id(top), \application_menu),
                    div(id(content), Body)
                  ])).
```

Redefining the `head` can be used to pull in scripts, but typically html_requires//1 provides a more modular approach for pulling scripts and CSS-files.

**reply_html_partial**(`+HTML`)  
Reply with partial HTML document. The reply only contains the element from HTML, i.e., this predicate does not add a `DOCTYPE` header, `html`, `head` or `body`. It is intended for JavaScript handlers that request a partial document and insert that somewhere into the existing page DOM. See [reply_html_page/3](#reply_html_page/3) to reply with a complete (valid) HTML page.

**print_html**(`+List`)  
Print the token list to the Prolog current output stream.

**print_html**(`+Stream, +List`)  
Print the token list to the specified output stream

**html_print_length**(`+List, -Length`)  
When calling html_print/\[1,2\] on `List`, `Length` characters will be produced. Knowing the length is needed to provide the `Content-length` field of an HTTP reply-header.

#### 3.21.2 Repositioning HTML for CSS and javascript links

Modern HTML commonly uses CSS and Javascript. This requires \<link\> elements in the HTML \<head\> element or \<script\> elements in the \<body\>. Unfortunately this seriously harms re-using HTML DCG rules as components as each of these components may rely on their own style sheets or JavaScript code. We added a‘mailing’system to reposition and collect fragments of HTML. This is implemented by [html_post//2](#html_post//2), [html_receive//1](#html_receive//1) and [html_receive//2](#html_receive//2).

\[det\]**html_post**(`+Id, :HTML`)`//`  
Reposition `HTML` to the receiving `Id`. The [html_post//2](#html_post//2) call processes `HTML` using [html//1](#html//1). Embedded `\`-commands are executed by mailman/1 from [print_html/1](#print_html/1) or [html_print_length/2](#html_print_length/2). These commands are called in the calling context of the [html_post//2](#html_post//2) call.

A typical usage scenario is to get required CSS links in the document head in a reusable fashion. First, we define css//1 as:

``` code
css(URL) -->
        html_post(css,
                  link([ type('text/css'),
                         rel('stylesheet'),
                         href(URL)
                       ])).
```

Next we insert the *unique* CSS links, in the pagehead using the following call to [reply_html_page/2](#reply_html_page/2):

``` code
        reply_html_page([ title(...),
                          \html_receive(css)
                        ],
                        ...)
```

\[det\]**html_receive**(`+Id`)`//`  
Receive posted HTML tokens. Unique sequences of tokens posted with [html_post//2](#html_post//2) are inserted at the location where [html_receive//1](#html_receive//1) appears.

See also  
\- The local predicate sorted_html//1 handles the output of [html_receive//1](#html_receive//1).  
- [html_receive//2](#html_receive//2) allows for post-processing the posted material.

\[det\]**html_receive**(`+Id, :Handler`)`//`  
This extended version of [html_receive//1](#html_receive//1) causes `Handler` to be called to process all messages posted to the channel at the time output is generated. `Handler` is called as below, where `PostedTerms` is a list of Module:Term created from calls to [html_post//2](#html_post//2). Module is the context module of html_post and Term is the unmodified term. Members in `PostedTerms` are in the order posted and may contain duplicates.

``` code
  phrase(Handler, PostedTerms, HtmlTerms, Rest)
```

Typically, `Handler` collects the posted terms, creating a term suitable for [html//1](#html//1) and finally calls [html//1](#html//1).

The library predefines the receiver channel `head` at the end of the `head` element for all pages that write the html `head` through this library. The following code can be used anywhere inside an HTML generating rule to demand a javascript in the header:

``` code
js_script(URL) -->
        html_post(head, script([ src(URL),
                                 type('text/javascript')
                               ], [])).
```

This mechanism is also exploited to add XML namespace (`xmlns`) declarations to the (outer) `html` element using xhml_ns//2:

**xhtml_ns**(`+Id, +Value`)`//`  
Demand an xmlns:id=`Value` in the outer html tag. This uses the html_post/2 mechanism to post to the `xmlns` channel. Rdfa ([http://www.w3.org/2006/07/SWD/RDFa/syntax/)](http://www.w3.org/2006/07/SWD/RDFa/syntax/)), embedding RDF in (x)html provides a typical usage scenario where we want to publish the required namespaces in the header. We can define:

``` code
rdf_ns(Id) -->
        { rdf_global_id(Id:'', Value) },
        xhtml_ns(Id, Value).
```

After which we can use rdf_ns//1 as a normal rule in [html//1](#html//1) to publish namespaces from `library(semweb/rdf_db)`. Note that this macro only has effect if the dialect is set to `xhtml`. In `html` mode it is silently ignored.

The required `xmlns` receiver is installed by [html_begin//1](#html_begin//1) using the `html` tag and thus is present in any document that opens the outer `html` environment through this library.

#### 3.21.3 Adding rules for html//1

In some cases it is practical to extend the translations imposed by html//1. We used this technique to define translation rules for the output of the SWI-Prolog `library(sgml)` package.

The html//1 non-terminal first calls the multifile ruleset html_write:expand//1.

**html_write:expand**(`+Spec`)`//`  
Hook to add additional translation rules for html//1.

**html_quoted**(`+Atom`)`//`  
Emit the text in `Atom`, inserting entity-references for the SGML special characters `<&>`.

**html_quoted_attribute**(`+Atom`)`//`  
Emit the text in `Atom` suitable for use as an SGML attribute, inserting entity-references for the SGML special characters `<&>"`.

#### 3.21.4 Generating layout

Though not strictly necessary, the library attempts to generate reasonable layout in SGML output. It does this only by inserting newlines before and after tags. It does this on the basis of the multifile predicate [html_write:layout/3](#html_write:layout/3)

**html_write:layout**(`+Tag, -Open, -Close`)  
Specify the layout conventions for the element `Tag`, which is a lowercase atom. `Open` is a term `Pre`-`Post`. It defines that the element should have at least `Pre` newline characters before and `Post` after the tag. The `Close` specification is similar, but in addition allows for the atom `-`, requesting the output generator to omit the close-tag altogether or `empty`, telling the library that the element has declared empty content. In this case the close-tag is not emitted either, but in addition html//1 interprets `Arg` in `Tag(Arg)` as a list of attributes rather than the content.

A tag that does not appear in this table is emitted without additional layout. See also [print_html/\[1,2\]](#print_html/1). Please consult the library source for examples.

#### 3.21.5 Examples for using the HTML write library

In the following example we will generate a table of Prolog predicates we find from the SWI-Prolog help system based on a keyword. The primary database is defined by the predicate predicate/5 We will make hyperlinks for the predicates pointing to their documentation.

``` code
:- use_module(library(http/html_write)).
:- use_module(library(pldoc/man_index)).
:- use_module(library(uri)).

html_apropos(Kwd) :-
    findall(Pred, apropos_predicate(Kwd, Pred), Matches),
    phrase(apropos_page(Kwd, Matches), Tokens),
    print_html(Tokens).

%       emit page with title, header and table of matches

apropos_page(Kwd, Matches) -->
    page([ title(['Predicates for ', Kwd])
         ],
         [ h2(align(center),
              ['Predicates for ', Kwd]),
           table([ align(center),
                   border(1),
                   width('80%')
                 ],
                 [ tr([ th('Predicate'),
                        th('Summary')
                      ])
                 | \apropos_rows(Matches)
                 ])
         ]).

%       emit the rows for the body of the table.

apropos_rows([]) -->
    [].
apropos_rows([pred(Name, Arity, Summary)|T]) -->
    html([ tr([ td(\predref(Name/Arity)),
                td(em(Summary))
              ])
         ]),
    apropos_rows(T).

%!  predref(Name/Arity)//
%
%   Emit Name/Arity as a hyperlink to
%
%           /cgi-bin/plman?name=Name&arity=Arity

predref(Name/Arity) -->
    { uri_edit([search([name=Name,arity=Arity])],
               '/cgi-bin/plman', Href)
    },
    html(a(href(Href), [Name, /, Arity])).

%       Find predicates from a keyword.

apropos_predicate(Pattern, pred(Name, Arity, Summary)) :-
    man_object_property(Name/Arity, summary(Summary)),
    (   sub_atom_icasechk(Name, _, Pattern)
    ->  true
    ;   sub_atom_icasechk(Summary, _, Pattern)
    ).
```

#### 3.21.6 Remarks on the `library(http/html_write)` library

This library is the result of various attempts to reach at a more satisfactory and Prolog-minded way to produce HTML text from a program. We have been using Prolog for the generation of web pages in a number of projects. Just using format/2 never was not a real option, generating error-prone HTML from clumsy syntax. We started with a layer on top of format/2, keeping track of the current nesting and thus always capable of properly closing the environment.

DCG based translation however, naturally exploits Prolog's term-rewriting primitives. If generation fails for whatever reason it is easy to produce an alternative document (for example holding an error message).

In a future version we will probably define a goal_expansion/2 to do compile-time optimisation of the library. Quotation of known text and invocation of sub-rules using the `\``RuleSet` and \<`Module`\>:\<`RuleSet`\> operators are costly operations in the analysis that can be done at compile-time.

### 3.22 library(http/js_write): Utilities for including JavaScript

This library is a supplement to `library(http/html_write)` for producing JavaScript fragments. Its main role is to be able to call JavaScript functions with valid arguments constructed from Prolog data. For example, suppose you want to call a JavaScript functions to process a list of names represented as Prolog atoms. This can be done using the call below, while without this library you would have to be careful to properly escape special characters.

``` code
numbers_script(Names) -->
    html(script(type('text/javascript'),
         [ \js_call('ProcessNumbers'(Names)
         ]),
```

The accepted arguments are described with [js_expression//1](#js_expression//1).

\[det\]**js_script**(`+Content`)`//`  
Generate a JavaScript `script` element with the given content.

\[det\]**javascript**(`+Content, +Vars, +VarDict, -DOM`)  
Quasi quotation parser for JavaScript that allows for embedding Prolog variables to substitude *identifiers* in the JavaScript snippet. Parameterizing a JavaScript string is achieved using the JavaScript `+` operator, which results in concatenation at the client side.

``` code
    ...,
    js_script({|javascript(Id, Config)||
                $(document).ready(function() {
                   $("#"+Id).tagit(Config);
                 });
               |}),
    ...
```

The current implementation tokenizes the JavaScript input and yields syntax errors on unterminated comments, strings, etc. No further parsing is implemented, which makes it possible to produce syntactically incorrect and partial JavaScript. Future versions are likely to include a full parser, generating syntax errors.

The parser produces a term `\List`, which is suitable for [js_script//1](#js_script//1) and [html//1](#html//1). Embedded variables are mapped to `\js_expression(Var)`, while the remaining text is mapped to atoms.

To be done  
Implement a full JavaScript parser. Users should *not* rely on the ability to generate partial JavaScript snippets.

\[det\]**js_call**(`+Term`)`//`  
Emit a call to a Javascript function. The Prolog functor is the name of the function. The arguments are converted from Prolog to JavaScript using [js_arg_list//1](#js_arg_list//1). Please not that Prolog functors can be quoted atom and thus the following is legal:

``` code
    ...
    html(script(type('text/javascript'),
         [ \js_call('x.y.z'(hello, 42))
         ]),
```

\[det\]**js_new**(`+Id, +Term`)`//`  
Emit a call to a Javascript object declaration. This is the same as:

``` code
['var ', Id, ' = new ', \js_call(Term)]
```

\[det\]**js_arg_list**(`+Expressions:list`)`//`  
Write javascript (function) arguments. This writes "(", Arg, ..., ")". See [js_expression//1](#js_expression//1) for valid argument values.

\[det\]**js_expression**(`+Expression`)`//`  
Emit a single JSON argument. `Expression` is one of:

**Variable**  
Emitted as Javascript `null`

**List**  
Produces a Javascript list, where each element is processed by this library.

**`object(Attributes)`**  
Where Attributes is a Key-Value list where each pair can be written as Key-Value, Key=Value or Key(Value), accommodating all common constructs for this used in Prolog.\< \$ { K:V, ... } Same as `object(Attributes)`, providing a more JavaScript-like syntax. This may be useful if the object appears literally in the source-code, but is generally less friendly to produce as a result from a computation.

**Dict**  
Emit a dict as a JSON object using json_write_dict/3.

**`json(Term)`**  
Emits a term using json_write/3.

**@(Atom)**  
Emits these constants without quotes. Normally used for the symbols `true`, `false` and `null`, but can also be use for emitting JavaScript symbols (i.e. function- or variable names).

**Number**  
Emitted literally

**`symbol(Atom)`**  
Synonym for @(Atom). Deprecated.

**Atom or String**  
Emitted as quoted JavaScript string.

\[semidet\]**js_arg**(`+Expression`)`//`  
Same as [js_expression//1](#js_expression//1), but fails if `Expression` is invalid, where [js_expression//1](#js_expression//1) raises an error.

deprecated  
New code should use [js_expression//1](#js_expression//1).

### 3.23 library(http/http_path): Abstract specification of HTTP server locations

This module provides an abstract specification of HTTP server locations that is inspired on absolute_file_name/3. The specification is done by adding rules to the dynamic multifile predicate [http:location/3](#http:location/3). The specification is very similar to user:file_search_path/2, but takes an additional argument with options. Currently only one option is defined:

**priority**(`+Integer`)  
If two rules match, take the one with highest priority. Using priorities is needed because we want to be able to overrule paths, but we do not want to become dependent on clause ordering.

The default priority is 0. Note however that notably libraries may decide to provide a fall-back using a negative priority. We suggest -100 for such cases.

This library predefines a single location at priority -100:

**root**  
The root of the server. Default is /, but this may be overruled using the setting (see setting/2) `http:prefix`

To serve additional resource files such as CSS, JavaScript and icons, see `library(http/http_server_files)`.

Here is an example that binds `/login` to login/1. The user can reuse this application while moving all locations using a new rule for the admin location with the option `[priority(10)]`.

``` code
:- multifile http:location/3.
:- dynamic   http:location/3.

http:location(admin, /, []).

:- http_handler(admin(login), login, []).

login(Request) :-
        ...
```

\[nondet,multifile\]http:**location**(`+Alias, -Expansion, -Options`)  
Multifile hook used to specify new HTTP locations. `Alias` is the name of the abstract path. `Expansion` is either a term Alias2(Relative), telling [http_absolute_location/3](#http_absolute_location/3) to translate `Alias` by first translating Alias2 and then applying the relative path Relative or, `Expansion` is an absolute location, i.e., one that starts with a `/`. `Options` currently only supports the priority of the path. If [http:location/3](#http:location/3) returns multiple solutions the one with the highest priority is selected. The default priority is 0.

This library provides a default for the abstract location `root`. This defaults to the setting http:prefix or, when not available to the path `/`. It is advised to define all locations (ultimately) relative to `root`. For example, use `root('home.html')` rather than `'/home.html'`.

\[det\]**http_absolute_uri**(`+Spec, -URI`)  
`URI` is the absolute (i.e., starting with `http://`) `URI` for the abstract specification `Spec`. Use [http_absolute_location/3](#http_absolute_location/3) to create references to locations on the same server.

\[det\]**http_absolute_location**(`+Spec, -Path, +Options`)  
`Path` is the HTTP location for the abstract specification `Spec`. `Options`:

**relative_to**(`Base`)  
`Path` is made relative to `Base`. Default is to generate absolute URLs.

See also  
[http_absolute_uri/2](#http_absolute_uri/2) to create a reference that can be used on another server.

**http_clean_location_cache**  
HTTP locations resolved through [http_absolute_location/3](#http_absolute_location/3) are cached. This predicate wipes the cache. The cache is automatically wiped by make/0 and if the setting http:prefix is changed.

### 3.24 library(http/html_head): Automatic inclusion of CSS and scripts links

To be done  
\- Possibly we should add img//2 to include images from symbolic path notation.  
- It would be nice if the HTTP file server could use our location declarations.

This library allows for abstract declaration of available CSS and Javascript resources and their dependencies using [html_resource/2](#html_resource/2). Based on these declarations, html generating code can declare that it depends on specific CSS or Javascript functionality, after which this library ensures that the proper links appear in the HTML head. The implementation is based on mail system implemented by html_post/2 of library `html_write.pl`.

Declarations come in two forms. First of all http locations are declared using the `http_path.pl` library. Second, [html_resource/2](#html_resource/2) specifies HTML resources to be used in the `head` and their dependencies. Resources are currently limited to Javascript files (.js) and style sheets (.css). It is trivial to add support for other material in the head. See html_include//1.

For usage in HTML generation, there is the DCG rule [html_requires//1](#html_requires//1) that demands named resources in the HTML head.

#### 3.24.1 About resource ordering

All calls to [html_requires//1](#html_requires//1) for the page are collected and duplicates are removed. Next, the following steps are taken:

1.  Add all dependencies to the set
2.  Replace multiple members by‘aggregate’scripts or css files. see use_agregates/4.
3.  Order all resources by demanding that their dependencies precede the resource itself. Note that the ordering of resources in the dependency list is **ignored**. This implies that if the order matters the dependency list must be split and only the primary dependency must be added.

#### 3.24.2 Debugging dependencies

Use `?-` `debug(html(script))`. to see the requested and final set of resources. All declared resources are in html_resource/3. The edit/1 command recognises the names of HTML resources.

#### 3.24.3 Predicates

\[det\]**html_resource**(`+About, +Properties`)  
Register an HTML head resource. `About` is either an atom that specifies an HTTP location or a term Alias(Sub). This works similar to absolute_file_name/2. See http:location_path/2 for details. Recognised properties are:

**requires**(`+Requirements`)  
Other required script and css files. If this is a plain file name, it is interpreted relative to the declared resource. `Requirements` can be a list, which is equivalent to multiple requires properties.

**virtual**(`+Bool`)  
If `true` (default `false`), do not include `About` itself, but only its dependencies. This allows for defining an alias for one or more resources.

**ordered**(`+Bool`)  
Defines that the list of requirements is ordered, which means that each requirement in the list depends on its predecessor.

**aggregate**(`+List`)  
States that `About` is an aggregate of the resources in `List`. This means that if both `About` and one of the elements of `List` appears in the dependencies, `About` is kept and the smaller one is dropped. If there are a number of dependencies on the small members, these are replaced with dependency on the big (aggregate) one, for example, to specify that a big javascript is actually the composition of a number of smaller ones.

**mime_type**(`-Mime`)  
May be specified for non-virtual resources to specify the mime-type of the resource. By default, the mime type is derived from the file name using file_mime_type/2.

Registering the same `About` multiple times extends the properties defined for `About`. In particular, this allows for adding additional dependencies to a (virtual) resource.

\[nondet\]**html_current_resource**(`?About`)  
True when `About` is a currently known resource.

\[det\]**html_requires**(`+ResourceOrList`)`//`  
Include `ResourceOrList` and all dependencies derived from it and add them to the HTML `head` using html_post/2. The actual dependencies are computed during the HTML output phase by [html_insert_resource//1](#html_insert_resource//1).

\[det\]**html_insert_resource**(`+ResourceOrList`)`//`  
Actually include HTML head resources. Called through [html_post//2](#html_post//2) from [html_requires//1](#html_requires//1) after rewrite by html_head_expansion/2. We are guaranteed we will only get one call that is passed a flat list of requested requirements. We have three jobs:

1.  Figure out all indirect requirements
2.  See whether we can use any‘aggregate’resources
3.  Put required resources before their requiree.

\[semidet,multifile\]**mime_include**(`+Mime, +Path`)`//`  
Hook called to include a link to an HTML resource of type `Mime` into the HTML head. The `Mime` type is computed from `Path` using file_mime_type/2. If the hook fails, two built-in rules for `text/css` and `text/javascript` are tried. For example, to include a =.pl= files as a Prolog script, use:

``` code
:- multifile
    html_head:mime_include//2.

html_head:mime_include(text/'x-prolog', Path) --> !,
    html(script([ type('text/x-prolog'),
                  src(Path)
                ],  [])).
```

### 3.25 library(http/http_pwp): Serve PWP pages through the HTTP server

To be done  
\- Support elements in the HTML header that allow controlling the page, such as setting the CGI-header, authorization, etc.  
- Allow external styling. Pass through [reply_html_page/2](#reply_html_page/2)? Allow filtering the DOM before/after PWP?

This module provides convenience predicates to include PWP (Prolog Well-formed Pages) in a Prolog web-server. It provides the following predicates:

`pwp_handler` **`/`** `2`  
This is a complete web-server aimed at serving static pages, some of which include PWP. This API is intended to allow for programming the web-server from a hierarchy of pwp files, prolog files and static web-pages.

`reply_pwp_page` **`/`** `3`  
Return a single PWP page that is executed in the context of the calling module. This API is intended for individual pages that include so much text that generating from Prolog is undesirable.

&nbsp;

**pwp_handler**(`+Options, +Request`)  
Handle PWP files. This predicate is defined to create a simple HTTP server from a hierarchy of PWP, HTML and other files. The interface is kept compatible with the `library(http/http_dispatch)`. In the typical usage scenario, one needs to define an http location and a file-search path that is used as the root of the server. E.g., the following declarations create a self-contained web-server for files in `/web/pwp/`.

``` code
user:file_search_path(pwp, '/web/pwp').

:- http_handler(root(.), pwp_handler([path_alias(pwp)]), [prefix]).
```

`Options` include:

**path_alias**(`+Alias`)  
Search for PWP files as `Alias`(Path). See absolute_file_name/3.

**index**(`+Index`)  
Name of the directory index (pwp) file. This option may appear multiple times. If no such option is provided, [pwp_handler/2](#pwp_handler/2) looks for `index.pwp`.

**view**(`+Boolean`)  
If `true` (default is `false`), allow for ?view=source to serve PWP file as source.

**index_hook**(`:Hook`)  
If a directory has no index-file, [pwp_handler/2](#pwp_handler/2) calls `Hook`(PhysicalDir, `Options`, `Request`). If this semidet predicate succeeds, the request is considered handled.

**hide_extensions**(`+List`)  
Hide files of the given extensions. The default is to hide .pl files.

**dtd**(`?DTD`)  
`DTD` to parse the input file with. If unbound, the generated `DTD` is returned

Errors  
`permission_error(index, http_location, Location)` is raised if the handler resolves to a directory that has no index.

See also  
[reply_pwp_page/3](#reply_pwp_page/3)

**reply_pwp_page**(`:File, +Options, +Request`)  
Reply a PWP file. This interface is provided to server individual locations from PWP files. Using a PWP file rather than generating the page from Prolog may be desirable because the page contains a lot of text (which is cumbersome to generate from Prolog) or because the maintainer is not familiar with Prolog.

`Options` supported are:

**mime_type**(`+Type`)  
Serve the file using the given mime-type. Default is text/html.

**unsafe**(`+Boolean`)  
Passed to [http_safe_file/2](#http_safe_file/2) to check for unsafe paths.

**pwp_module**(`+Boolean`)  
If `true`, (default `false`), process the PWP file in a module constructed from its canonical absolute path. Otherwise, the PWP file is processed in the calling module.

Initial context:

**`SCRIPT_NAME`**  
Virtual path of the script.

**`SCRIPT_DIRECTORY`**  
Physical directory where the script lives

**`QUERY`**  
Var=Value list representing the query-parameters

**`REMOTE_USER`**  
If access has been authenticated, this is the authenticated user.

**`REQUEST_METHOD`**  
One of `get`, `post`, `put` or `head`

**`CONTENT_TYPE`**  
Content-type provided with HTTP POST and PUT requests

**`CONTENT_LENGTH`**  
Content-length provided with HTTP POST and PUT requests

While processing the script, the file-search-path pwp includes the current location of the script. I.e., the following will find myprolog in the same directory as where the PWP file resides.

``` code
pwp:ask="ensure_loaded(pwp(myprolog))"
```

See also  
[pwp_handler/2](#pwp_handler/2).

To be done  
complete the initial context, as far as possible from CGI variables. See [http://hoohoo.ncsa.illinois.edu/docs/cgi/env.html](http://hoohoo.ncsa.illinois.edu/docs/cgi/env.html)

### 3.26 library(http/htmx): Support htmx.org

Quoted from htmx.org:

> [htmx](https://htmx.org) gives you access to AJAX, CSS Transitions, WebSockets and Server Sent Events directly in HTML, using attributes, so you can build modern user interfaces with the simplicity and power of hypertext

The idea behind htmx is to allow adding attributes to any HTML element that cause an HTTP request. The HTTP response is typically a (short) HTML fragment that extends or replaces an element on the page. This allows us to program a most functionality interactive seen in modern web applications using the powerful SWI-Prolog HTML generation framework rather than having to write a JSON backend and accompanying JavaScript frontend that runs in the browser.

Below is a minimalistic, yet fully functional application that replaces a button after a click in two ways, using either a direct hx-swap our an out-of-band hx-swap command.

``` code
:- use_module(library(http/http_server)).
:- use_module(library(http/htmx)).
:- use_module(library(main)).

:- initialization(main, main).

main(_Argv) :-
    http_server([port(8080)]),
    thread_get_message(quit).

http:location(htmx, root(htmx), []).

:- http_handler(root(.), home, []).

home(_Request) :-
    reply_html_page(
        [ title('HTMX demo'),
          script(src('https://unpkg.com/htmx.org'), [])
        ],
        [ button([ id(button1),
                   'hx-post'('/htmx/clicked1'),
                   'hx-swap'('outerHTML')
                 ],
                 'Click me (1)'),
          button([ id(button2),
                   'hx-post'('/htmx/clicked2')
                 ],
                 'Click me (2)')
        ]).

:- http_handler(htmx(clicked1), reply_htmx(\clicked1), []).
:- http_handler(htmx(clicked2), reply_htmx(\clicked2), []).

clicked1 -->
    html('Thanks for clicking me! (1)').

clicked2 -->
    htmx_oob(button2, html('Thanks for clicking me! (2)')).
```

HTMX requires no dedicated support from the server. This library provides [reply_htmx/1](#reply_htmx/1),2 to reply with a single HTML element rather than an entire page. Future versions of this library may provide some additional utility predicates.

\[det\]**reply_htmx**(`+HTML`)  
\[det\]**reply_htmx**(`+HTML, +Request`)  
Reply a plain `HTML` element as opposed to a complete `HTML` page as created using [reply_html_page/2](#reply_html_page/2),3. While [reply_htmx/1](#reply_htmx/1) is to be used in a normal HTTP handler (route), [reply_htmx/2](#reply_htmx/2) may be registered directly in the [http_handler/3](#http_handler/3) declaration to deal with simple cases where we do not need the `Request` data.

\[det\]**htmx_oob**(`++Target, :HTML`)`//`  
Emit an htmx out-of-band element. `HTML` is used to swap the content of the DOM element with id `Target`.

## 4 HTTP and IPv6

As of version 9.1.5, SWI-Prolog supports IPv6. This has few implications for the HTTP package because most aspects are handled by `library(socket)` and `library(uri)`. This sections highlights a few aspects.

The *client libraries* use [http_open/3](#http_open/3), which in turn uses tcp_connect/3. This causes the client to use addresses returned by host_address/3, which is based on the C API **getaddrinfo()**, in the order provided by **getaddrinfo()**. The URL is parsed using `library(uri)`, which allows enclosing IPv6 addresses in `[]`. The query below accesses an IPv6 server on localhost at port 8080

``` code
?- http_open('http://[::1]:8080', Stream, []).
```

The predicate [http_server/2](#http_server/2) can be used to create an IPv6 server using one of the queries below. The first binds to all interfaces. The second only binds to the IPv6 equivalent of `localhost`. Note that the IPv6 address needs to be quoted to create the desired `Host`:`Port` term.

``` code
?- http_server(Goal,[port('::':8080)]).
?- http_server(Goal,[port('::1':8080)]).
```

## 5 Transfer encodings

The HTTP protocol provides for *transfer encodings*. These define filters applied to the data described by the `Content-type`. The two most popular transfer encodings are `chunked` and `deflate`. The `chunked` encoding avoids the need for a `Content-length` header, sending the data in chunks, each of which is preceded by a length. The `deflate` encoding provides compression.

Transfer-encodings are supported by filters defined as foreign libraries that realise an encoding/decoding stream on top of another stream. Currently there are two such libraries: `library(http/http_chunked.pl)` and `library(zlib.pl)`.

There is an emerging hook interface dealing with transfer encodings. The `library(http/http_chunked.pl)` provides a hook used by `library(http/http_open.pl)` to support chunked encoding in [http_open/3](#http_open/3). Note that both `http_open.pl` *and* `http_chunked.pl` must be loaded for [http_open/3](#http_open/3) to support chunked encoding.

### 5.1 The `library(http/http_chunked)` library

**http_chunked_open**(`+RawStream, -DataStream, +Options`)  
Create a stream to realise HTTP chunked encoding or decoding. The technique is similar to library(zlib), using a Prolog stream as a filter on another stream. See online documentation at [http://www.swi-prolog.org/](http://www.swi-prolog.org/) for details.

## 6 library(http/websocket): WebSocket support

See also  
RFC 6455, [http://tools.ietf.org/html/rfc6455](http://tools.ietf.org/html/rfc6455)

To be done  
Deal with protocol extensions.

WebSocket is a lightweight message oriented protocol on top of TCP/IP streams. It is typically used as an *upgrade* of an HTTP connection to provide bi-directional communication, but can also be used in isolation over arbitrary (Prolog) streams.

The SWI-Prolog interface is based on *streams* and provides [ws_open/3](#ws_open/3) to create a *websocket stream* from any Prolog stream. Typically, both an input and output stream are wrapped and then combined into a single object using stream_pair/3.

The high-level interface provides [http_upgrade_to_websocket/3](#http_upgrade_to_websocket/3) to realise a websocket inside the HTTP server infrastructure and [http_open_websocket/3](#http_open_websocket/3) as a layer over [http_open/3](#http_open/3) to realise a client connection. After establishing a connection, [ws_send/2](#ws_send/2) and [ws_receive/2](#ws_receive/2) can be used to send and receive messages. The predicate [ws_close/3](#ws_close/3) is provided to perform the closing handshake and dispose of the stream objects.

\[det\]**http_open_websocket**(`+URL, -WebSocket, +Options`)  
Establish a client websocket connection. This predicate calls [http_open/3](#http_open/3) with additional headers to negotiate a websocket connection. In addition to the options processed by [http_open/3](#http_open/3), the following options are recognised:

**subprotocols**(`+List`)  
`List` of subprotocols that are acceptable. The selected protocol is available as ws_property(`WebSocket`, `subprotocol(Protocol)`.

Note that clients often provide an `Origin` header and some servers require this field. See RFC 6455 for details. By default this predicate does not set `Origin`. It may be set using the `request_header` option of [http_open/3](#http_open/3), e.g. by passing this in the `Options` list:

``` code
request_header('Origin' = 'https://www.swi-prolog.org')
```

The following example exchanges a message with the html5rocks.websocket.org echo service:

``` code
?- URL = 'ws://html5rocks.websocket.org/echo',
   http_open_websocket(URL, WS, []),
   ws_send(WS, text('Hello World!')),
   ws_receive(WS, Reply),
   ws_close(WS, 1000, "Goodbye").
URL = 'ws://html5rocks.websocket.org/echo',
WS = <stream>(0xe4a440,0xe4a610),
Reply = websocket{data:"Hello World!", opcode:text}.
```

|             |                                      |
|-------------|--------------------------------------|
| `WebSocket` | is a stream pair (see stream_pair/3) |

**http_upgrade_to_websocket**(`:Goal, +Options, +Request`)  
Create a websocket connection running `call(Goal, WebSocket)`, where WebSocket is a socket-pair. `Options`:

**guarded**(`+Boolean`)  
If `true` (default), guard the execution of `Goal` and close the websocket on both normal and abnormal termination of `Goal`. If `false`, `Goal` itself is responsible for the created websocket if `Goal` succeeds. The websocket is closed if `Goal` fails or raises an exception. This can be used to create a single thread that manages multiple websockets using I/O multiplexing. See `library(http/hub)`.

**subprotocols**(`+List`)  
`List` of acceptable subprotocols.

**timeout**(`+TimeOut`)  
Timeout to apply to the input stream. Default is `infinite`.

Note that the `Request` argument is the last for cooperation with [http_handler/3](#http_handler/3). A simple *echo* server that can be accessed at =/ws/= can be implemented as:

``` code
:- use_module(library(http/websocket)).
:- use_module(library(http/thread_httpd)).
:- use_module(library(http/http_dispatch)).

:- http_handler(root(ws),
                http_upgrade_to_websocket(echo, []),
                [spawn([])]).

echo(WebSocket) :-
    ws_receive(WebSocket, Message),
    (   Message.opcode == close
    ->  true
    ;   ws_send(WebSocket, Message),
        echo(WebSocket)
    ).
```

throws  
`switching_protocols(Goal, Options)`. The recovery from this exception causes the HTTP infrastructure to call `call(Goal, WebSocket)`.

See also  
[http_switch_protocol/2](#http_switch_protocol/2).

\[det\]**ws_send**(`+WebSocket, +Message`)  
Send a message over a websocket. The following terms are allowed for `Message`:

**text**(`+Text`)  
Send a text message. `Text` is serialized using write/1.

**binary**(`+Content`)  
As `text(+Text)`, but all character codes produced by `Content` must be in the range \[0..255\]. Typically, `Content` will be an atom or string holding binary data.

**prolog**(`+Term`)  
Send a Prolog term as a text message. Text is serialized using write_canonical/1.

**json**(`+JSON`)  
Send the Prolog representation of a `JSON` term using json_write_dict/2.

**string**(`+Text`)  
Same as `text(+Text)`, provided for consistency.

**close**(`+Code, +Text`)  
Send a close message. `Code` is 1000 for normal close. See websocket documentation for other values.

**`Dict`**  
A dict that minimally contains an `opcode` key. Other keys used are:

`format`**`:`**`Format`  
Serialization format used for `Message`.data. `Format` is one of `string`, `prolog` or `json`. See [ws_receive/3](#ws_receive/3).

`data`**`:`**`Term`  
If this key is present, it is serialized according to `Message`.format. Otherwise it is serialized using write/1, which implies that string and atoms are just sent verbatim.

Note that ws_start_message/3 does not unlock the stream. This is done by ws_send/1. This implies that multiple threads can use [ws_send/2](#ws_send/2) and the messages are properly serialized.

To be done  
Provide serialization details using options.

\[det\]**ws_receive**(`+WebSocket, -Message:dict`)  
\[det\]**ws_receive**(`+WebSocket, -Message:dict, +Options`)  
Receive the next message from `WebSocket`. `Message` is a dict containing the following keys:

`opcode`**`:`**`OpCode`  
`OpCode` of the message. This is an atom for known opcodes and an integer for unknown ones. If the peer closed the stream, `OpCode` is bound to `close` and data to the atom `end_of_file`.

`data`**`:`**`String`  
The data, represented as a string. This field is always present. `String` is the empty string if there is no data in the message.

`rsv`**`:`**`RSV`  
Present if the `WebSocket` `RSV` header is not 0. `RSV` is an integer in the range \[1..7\].

If `ping` message is received and `WebSocket` is a stream pair, ws_receive/1 replies with a `pong` and waits for the next message.

The predicate [ws_receive/3](#ws_receive/3) processes the following options:

**format**(`+Format`)  
Defines how *text* messages are parsed. `Format` is one of

**string**  
Data is returned as a Prolog string (default)

**json**  
Data is parsed using json_read_dict/3, which also receives `Options`.

**prolog**  
Data is parsed using read_term/3, which also receives `Options`.

To be done  
Add a hook to allow for more data formats?

\[det\]**ws_close**(`+WebSocket:stream_pair, +Code, +Data`)  
Close a `WebSocket` connection by sending a `close` message if this was not already sent and wait for the close reply.

|  |  |
|----|----|
| `Code` | is the numerical code indicating the close status. This is 16-bit integer. The codes are defined in section *7.4.1. Defined Status Codes* of RFC6455. Notably, 1000 indicates a normal closure. |
| `Data` | is currently interpreted as text. |

Errors  
`websocket_error(unexpected_message, Reply)` if the other side did not send a close message in reply.

\[det\]**ws_open**(`+Stream, -WSStream, +Options`)  
Turn a raw TCP/IP (or any other binary stream) into a websocket stream. `Stream` can be an input stream, output stream or a stream pair. `Options` includes

**mode**(`+Mode`)  
One of `server` or `client`. If `client`, messages are sent as *masked*.

**buffer_size**(`+Count`)  
Send partial messages for each `Count` bytes or when flushing the output. The default is to buffer the entire message before it is sent.

**close_parent**(`+Boolean`)  
If `true` (default), closing `WSStream` also closes `Stream`.

**subprotocol**(`+Protocol`)  
Set the subprotocol property of WsStream. This value can be retrieved using [ws_property/2](#ws_property/2). `Protocol` is an atom. See also the `subprotocols` option of [http_open_websocket/3](#http_open_websocket/3) and [http_upgrade_to_websocket/3](#http_upgrade_to_websocket/3).

A typical sequence to turn a pair of streams into a WebSocket is here:

``` code
    ...,
    Options = [mode(server), subprotocol(chat)],
    ws_open(Input, WsInput, Options),
    ws_open(Output, WsOutput, Options),
    stream_pair(WebSocket, WsInput, WsOutput).
```

\[nondet\]**ws_property**(`+WebSocket, ?Property`)  
True if `Property` is a property `WebSocket`. Defined properties are:

**subprotocol**(`Protocol`)  
`Protocol` is the negotiated subprotocol. This is typically set as a property of the websocket by [ws_open/3](#ws_open/3).

**ws_mask**(`-Mask`)  
Produce a good random number of the mask of a client message.

## 7 library(http/hub): Manage a hub for websockets

To be done  
The current design does not use threads to perform tasks for multiple hubs. This implies that the design scales rather poorly for hosting many hubs with few users.

This library manages a hub that consists of clients that are connected using a websocket. Messages arriving at any of the websockets are sent to the *event* queue of the hub. In addition, the hub provides a *broadcast* interface. A typical usage scenario for a hub is a *chat server* A scenario for realizing an chat server is:

1.  Create a new hub using [hub_create/3](#hub_create/3).

2.  Create one or more threads that listen to Hub.queues.event from the created hub. These threads can update the shared view of the world. A message is a dict as returned by [ws_receive/2](#ws_receive/2) or a hub control message. Currently, the following control messages are defined:

    **hub**{`error``:``Error``, ``left``:``ClientId``, ``reason``:``Reason`}  
    A client left us because of an I/O error. `Reason` is `read` or `write` and `Error` is the Prolog I/O exception.

    **hub**{`joined``:``ClientId`}  
    A new client has joined the chatroom.

    The `thread(s)` can talk to clients using two predicates:

    - [hub_send/2](#hub_send/2) sends a message to a specific client
    - [hub_broadcast/2](#hub_broadcast/2) sends a message to all clients of the hub.

A hub consists of (currently) four message queues and a simple dynamic fact. Threads that are needed for the communication tasks are created on demand and die if no more work needs to be done.

\[det\]**hub_create**(`+Name, -Hub, +Options`)  
Create a new hub. `Hub` is a dict containing the following public information:

`Hub` **`.`** `name`  
The name of the hub (the `Name` argument)

`queues` **`.`** `event`  
Message queue to which the hub `thread(s)` can listen.

After creating a hub, the application normally creates a thread that listens to `Hub`.queues.event and exposes some mechanisms to establish websockets and add them to the hub using [hub_add/3](#hub_add/3).

See also  
[http_upgrade_to_websocket/3](#http_upgrade_to_websocket/3) establishes a websocket from the SWI-Prolog webserver.

\[nondet\]**current_hub**(`?Name, ?Hub`)  
True when there exists a hub `Hub` with `Name`.

\[det\]**hub_add**(`+Hub, +WebSocket, ?Id`)  
Add a `WebSocket` to the hub. `Id` is used to identify this user. It may be provided (as a ground term) or is generated as a UUID.

\[nondet\]**hub_member**(`?HubName, ?Id`)  
True when `Id` is a member of the hub `HubName`.

\[semidet\]**hub_send**(`+ClientId, +Message`)  
Send message to the indicated `ClientId`. Fails silently if `ClientId` does not exist.

|  |  |
|----|----|
| `Message` | is either a single message (as accepted by [ws_send/2](#ws_send/2)) or a list of such messages. |

\[det\]**hub_broadcast**(`+Hub, +Message`)  
\[det\]**hub_broadcast**(`+Hub, +Message, :Condition`)  
Send `Message` to all websockets associated with `Hub` for which `call(Condition, Id)` succeeds. Note that this process is *asynchronous*: this predicate returns immediately after putting all requests in a broadcast queue. If a message cannot be delivered due to a network error, the hub is informed through io_error/3.

## 8 MIME support

### 8.1 library(http/mimepack): Create a MIME message

Simple and partial implementation of MIME encoding. MIME is covered by RFC 2045. This library is used by e.g., [http_post_data/3](#http_post_data/3) when using the `form_data(+ListOfData)` input specification.

MIME decoding is now arranged through `library(mime)` from the clib package, based on the external librfc2045 library. Most likely the functionality of this package will be moved to the same library someday. Packing however is a lot simpler then parsing.

\[det\]**mime_pack**(`+Inputs, +Out:stream, ?Boundary`)  
Pack a number of inputs into a MIME package using a specified or generated boundary. The generated boundary consists of the current time in milliseconds since the epoch and 10 random hexadecimal numbers. `Inputs` is a list of *documents* that is added to the mime message. Each element is one of:

`Name` **`=`** `Value`  
`Name` the document. This emits a header of the form below. The `filename` is present if `Value` is of the form `file(File)`. `Value` may be any of remaining value specifications.

``` code
Content-Disposition: form-data; name="Name"[; filename="<File>"
```

**html**(`Tokens`)  
`Tokens` is a list of HTML tokens as produced by [html//1](#html//1). The token list is emitted using [print_html/1](#print_html/1).

**file**(`File`)  
Emit the contents of `File`. The `Content-type` is derived from the `File` using file_mime_type/2. If the content-type is `text/_`, the file data is copied in text mode, which implies that it is read in the default encoding of the system and written using the encoding of the `Out` stream. Otherwise the file data is copied binary.

**stream**(`In, Len`)  
Content is the next `Len` units from `In`. Data is copied using copy_stream_data/3. Units is bytes for binary streams and characters codes for text streams.

**stream**(`In`)  
Content of the stream `In`, copied using copy_stream_data/2. This is often used with memory files (see new_memory_file/1).

**mime**(`Attributes, Value,[]`)  
Create a MIME header from `Attributes` and add `Value`, which can be any of remaining values of this list. `Attributes` may contain `type(ContentType)` and/or `character_set(CharSet)`. This can be used to give a content-type to values that otherwise do not have a content-type. For example:

``` code
mime([type(text/html)], '<b>Hello World</b>', [])
```

**mime**(`[], , Parts`)  
Creates a nested multipart MIME message. `Parts` is passed as `Inputs` to a recursive call to mime_pack/2.

**`Atomic`**  
`Atomic` values are passed to write/1. This embeds simple atoms and numbers.

|  |  |
|----|----|
| `Out` | is a stream opened for writing. Typically, it should be opened in text mode using UTF-8 encoding. |

bug  
Does not validate that the boundary does not appear in any of the input documents.

## 9 Security

Writing servers is an inherently dangerous job that should be carried out with some considerations. You have basically started a program on a public terminal and invited strangers to use it. When using the interactive server or inetd based server the server runs under your privileges. Using CGI scripted it runs with the privileges of your web-server. Though it should not be possible to fatally compromise a Unix machine using user privileges, getting unconstrained access to the system is highly undesirable.

Symbolic languages have an additional handicap in their inherent possibilities to modify the running program and dynamically create goals (this also applies to the popular Perl and PHP scripting languages). Here are some guidelines.

- *Check your input*  
  Hardly anything can go wrong if you check the validity of query-arguments before formulating an answer.

- *Check filenames*  
  If part of the query consists of filenames or directories, check them. This also applies to files you only read. Passing names as `/etc/passwd`, but also `../../../../../etc/passwd` are tried by hackers to learn about the system they want to attack. So, expand provided names using absolute_file_name/\[2,3\] and verify they are inside a folder reserved for the server. Avoid symbolic links from this subtree to the outside world. The example below checks validity of filenames. The first call ensures proper canonisation of the paths to avoid an mismatch due to symbolic links or other filesystem ambiguities.

  ``` code
  check_file(File) :-
          absolute_file_name('/path/to/reserved/area', Reserved),
          absolute_file_name(File, Tried),
          sub_atom(Tried, 0, _, _, Reserved).
  ```

- *Check scripts*  
  Should input in any way activate external scripts using shell/1 or `open(pipe(Command), ...)`, verify the argument once more. Use process_create/3 in preference over shell/1 as this function avoids stringification of arguments (Unix) or ensures proper quoting of arguments (Windows).

- *Check meta-calling*  
  *The* attractive situation for you and your attacker is below:

  ``` code
  reply(Query) :-
          member(search(Args), Query),
          member(action=Action, Query),
          member(arg=Arg, Query),
          call(Action, Arg).              % NEVER EVER DO THIS!
  ```

  All your attacker has to do is specify `Action` as `shell` and `Arg` as `/bin/sh` and he has an uncontrolled shell!

## 10 Tips and tricks

- *URL Locations*  
  With an application in mind, it is tempting to make all URL locations short and directly connected to the root (`/`). This is *not* a good idea. It is advised to have all locations in a server below a directory with an informative name. Consider to make the root location something that can be changed using a global setting.
  - Page generating code can easily be reused. Using locations directly below the root however increases the likelihood of conflicts.
  - Multiple servers can be placed behind the same public server as explained in [section 3.14.7](#sec:3.14.7). Using a common and fairly unique root, redirection is much easier and less likely to lead to conflicts.
- *Debugging*  
  Debugging multi-threaded applications is possible using the graphical debugger. This implies requires that the xpce extension package must be installed. Spy-points may be placed using tspy/1.
- *URL/URI manipulation*  
  Sometimes an [http_handler/3](#http_handler/3) needs to inspect or normalize the URL. There are various utility predicates in `library(uri)`, such as uri_components/2, uri_data/4, uri_edit/3, uri_nomalized/2, etc.

## 11 Status

The SWI-Prolog HTTP library is in active use in a large number of projects. It is considered one of the SWI-Prolog core libraries that is actively maintained and regularly extended with new features. This is particularly true for the multi-threaded server. The inetd based server may be applicable for infrequent requests where the startup time is less relevant. The XPCE based server is considered obsolete.

This library is by no means complete and you are free to extend it.

# Index

?  
absolute_file_name/\[2,3\]  
[9](#idx:absolutefilename23:79)

chunked,encoding  
[5](#idx:chunkedencoding:75)

cleanup/2  
[2](#idx:cleanup2:3)

[cors_enable/0](#cors_enable/0)  
[cors_enable/2](#cors_enable/2)  
[current_hub/2](#current_hub/2)  
deflate,encoding  
[5](#idx:deflateencoding:76)

[directory_index//2](#directory_index//2)  
format/2  
[3.21](#idx:format2:48) [3.21.6](#idx:format2:68) [3.21.6](#idx:format2:69)

format/3  
[3.21](#idx:format3:51) [3.21](#idx:format3:52) [3.21](#idx:format3:53)

format_time/3  
[3.14.2](#idx:formattime3:34)

goal_expansion/2  
[3.21.6](#idx:goalexpansion2:70)

[health/2](#health/2)  
[hide/1](#hide/1)  
[hook/1](#hook/1)  
[hooked/0](#hooked/0)  
host_address/3  
[4](#idx:hostaddress3:73)

[html//1](#html//1)  
[html_begin//1](#html_begin//1)  
html_begin/1  
[3.21](#idx:htmlbegin1:60)

[html_current_resource/1](#html_current_resource/1)  
[html_end//1](#html_end//1)  
[html_insert_resource//1](#html_insert_resource//1)  
[html_post//2](#html_post//2)  
html_print/\[1,2\]  
[3.21.1](#idx:htmlprint12:64)

[html_print_length/2](#html_print_length/2)  
[html_quoted//1](#html_quoted//1)  
[html_quoted_attribute//1](#html_quoted_attribute//1)  
[html_receive//1](#html_receive//1)  
[html_receive//2](#html_receive//2)  
[html_requires//1](#html_requires//1)  
[html_resource/2](#html_resource/2)  
[html_write:expand//1](#html_write:expand//1)  
[html_write:layout/3](#html_write:layout/3)  
[3.21.4](#idx:htmlwritelayout3:65)

[htmx_oob//2](#htmx_oob//2)  
[http:///1](#http:///1)  
[http:authenticate/3](#http:authenticate/3)  
[http:authenticate_client/2](#http:authenticate_client/2)  
http:convert_parameter/3  
[3.12](#idx:httpconvertparameter3:22)

[http:disable_encoding_filter/1](#http:disable_encoding_filter/1)  
[http:location/3](#http:location/3)  
[http:mime_type_encoding/2](#http:mime_type_encoding/2)  
[3.1](#idx:httpmimetypeencoding2:11)

[http:mime_type_icon/2](#http:mime_type_icon/2)  
[http:open_options/2](#http:open_options/2)  
[http:post_data_hook/3](#http:post_data_hook/3)  
[http:request_expansion/2](#http:request_expansion/2)  
[3.15](#idx:httprequestexpansion2:46)

[http:schedule_workers/1](#http:schedule_workers/1)  
[http:serialize_reply/2](#http:serialize_reply/2)  
[http:sni_options/2](#http:sni_options/2)  
[http:status_page/3](#http:status_page/3)  
[http:status_page_hook/3](#http:status_page_hook/3)  
[3.10](#idx:httpstatuspagehook3:20)

[http:update_cookies/3](#http:update_cookies/3)  
[http:write_cookies/3](#http:write_cookies/3)  
[http_404/2](#http_404/2)  
[http_absolute_location/3](#http_absolute_location/3)  
[http_absolute_uri/2](#http_absolute_uri/2)  
[http_add_worker/2](#http_add_worker/2)  
[http_authenticate/3](#http_authenticate/3)  
[http_authorization_data/2](#http_authorization_data/2)  
[http_certificate_hook/3](#http_certificate_hook/3)  
[http_chunked_open/3](#http_chunked_open/3)  
[http_clean_location_cache/0](#http_clean_location_cache/0)  
[http_close_keep_alive/1](#http_close_keep_alive/1)  
[http_close_session/1](#http_close_session/1)  
[http_convert_data/4](#http_convert_data/4)  
[http_current_handler/2](#http_current_handler/2)  
[http_current_handler/3](#http_current_handler/3)  
[http_current_host/4](#http_current_host/4)  
[http_current_request/1](#http_current_request/1)  
[3.15](#idx:httpcurrentrequest1:47)

[http_current_session/2](#http_current_session/2)  
[http_current_user/3](#http_current_user/3)  
[http_current_worker/2](#http_current_worker/2)  
[http_daemon/0](#http_daemon/0)  
[http_daemon/1](#http_daemon/1)  
[http_delete/3](#http_delete/3)  
[http_delete_handler/1](#http_delete_handler/1)  
[http_digest_challenge//2](#http_digest_challenge//2)  
[http_digest_password_hash/4](#http_digest_password_hash/4)  
[http_digest_response/5](#http_digest_response/5)  
[http_disconnect/1](#http_disconnect/1)  
[http_dispatch/1](#http_dispatch/1)  
[3.14.2](#idx:httpdispatch1:33)

[http_get/3](#http_get/3)  
[2](#idx:httpget3:6)

[http_handler/3](#http_handler/3)  
[1](#idx:httphandler3:1) [3.1](#idx:httphandler3:9) [3.14.2](#idx:httphandler3:39) [3.21](#idx:httphandler3:57) [10](#idx:httphandler3:84)

[http_in_session/1](#http_in_session/1)  
[http_join_headers/3](#http_join_headers/3)  
[http_link_to_id/3](#http_link_to_id/3)  
[http_location_by_id/2](#http_location_by_id/2)  
[3.21](#idx:httplocationbyid2:55)

[http_log/2](#http_log/2)  
[http_log_close/1](#http_log_close/1)  
[http_log_stream/1](#http_log_stream/1)  
[http_logrotate/1](#http_logrotate/1)  
[http_open/3](#http_open/3)  
[2](#idx:httpopen3:2) [2](#idx:httpopen3:5) [4](#idx:httpopen3:71) [5](#idx:httpopen3:77) [5](#idx:httpopen3:78)

[http_open_session/2](#http_open_session/2)  
[http_open_websocket/3](#http_open_websocket/3)  
[http_opt_help/2](#http_opt_help/2)  
[http_opt_meta/2](#http_opt_meta/2)  
[http_opt_type/3](#http_opt_type/3)  
[http_parameters/2](#http_parameters/2)  
[3.12](#idx:httpparameters2:23)

[http_parameters/3](#http_parameters/3)  
[3.12](#idx:httpparameters3:24)

[http_parse_digest_challenge/2](#http_parse_digest_challenge/2)  
[http_parse_header/2](#http_parse_header/2)  
[http_parse_header_value/3](#http_parse_header_value/3)  
[http_patch/4](#http_patch/4)  
[http_post/4](#http_post/4)  
[2](#idx:httppost4:7)

[http_post_data/3](#http_post_data/3)  
[http_public_host/4](#http_public_host/4)  
[http_public_host_url/2](#http_public_host_url/2)  
[http_public_url/2](#http_public_url/2)  
[http_put/4](#http_put/4)  
[http_read_data/3](#http_read_data/3)  
[3.13.1](#idx:httpreaddata3:27)

[http_read_header/2](#http_read_header/2)  
[http_read_passwd_file/2](#http_read_passwd_file/2)  
[http_read_reply_header/2](#http_read_reply_header/2)  
[http_read_request/2](#http_read_request/2)  
[3.13](#idx:httpreadrequest2:25) [3.13](#idx:httpreadrequest2:26)

[http_redirect/3](#http_redirect/3)  
[3.1](#idx:httpredirect3:10)

[http_relative_path/2](#http_relative_path/2)  
[http_reload_with_parameters/3](#http_reload_with_parameters/3)  
[http_reply/2](#http_reply/2)  
[http_reply/3](#http_reply/3)  
[3.1.1](#idx:httpreply3:15) [3.1.1](#idx:httpreply3:16) [3.1.1](#idx:httpreply3:17)

[http_reply/4](#http_reply/4)  
[http_reply/5](#http_reply/5)  
[http_reply/6](#http_reply/6)  
[http_reply_dirindex/3](#http_reply_dirindex/3)  
[http_reply_file/3](#http_reply_file/3)  
[http_reply_from_files/3](#http_reply_from_files/3)  
[http_reply_header/3](#http_reply_header/3)  
[http_request_expansion/2](#http_request_expansion/2)  
[http_safe_file/2](#http_safe_file/2)  
[http_schedule_logrotate/2](#http_schedule_logrotate/2)  
[http_server/1](#http_server/1)  
[3.14.4](#idx:httpserver1:40)

[http_server/2](#http_server/2)  
[4](#idx:httpserver2:74)

[http_server_hook/1](#http_server_hook/1)  
[http_server_property/2](#http_server_property/2)  
[http_session_assert/1](#http_session_assert/1)  
[http_session_assert/2](#http_session_assert/2)  
[http_session_asserta/1](#http_session_asserta/1)  
[http_session_asserta/2](#http_session_asserta/2)  
[http_session_cookie/1](#http_session_cookie/1)  
[http_session_data/1](#http_session_data/1)  
[http_session_data/2](#http_session_data/2)  
[http_session_id/1](#http_session_id/1)  
[http_session_option/1](#http_session_option/1)  
[http_session_retract/1](#http_session_retract/1)  
[http_session_retract/2](#http_session_retract/2)  
[http_session_retractall/1](#http_session_retractall/1)  
[http_session_retractall/2](#http_session_retractall/2)  
[http_set_authorization/2](#http_set_authorization/2)  
[http_set_session/1](#http_set_session/1)  
[http_set_session/2](#http_set_session/2)  
[http_set_session_options/1](#http_set_session_options/1)  
[http_spawn/2](#http_spawn/2)  
[3](#idx:httpspawn2:8) [3.14.2](#idx:httpspawn2:35)

[http_status_reply/4](#http_status_reply/4)  
[http_status_reply/5](#http_status_reply/5)  
[http_status_reply/6](#http_status_reply/6)  
[http_stop_server/2](#http_stop_server/2)  
[http_switch_protocol/2](#http_switch_protocol/2)  
[http_timestamp/2](#http_timestamp/2)  
[http_update_connection/4](#http_update_connection/4)  
[http_update_encoding/3](#http_update_encoding/3)  
[http_update_transfer/4](#http_update_transfer/4)  
[http_upgrade_to_websocket/3](#http_upgrade_to_websocket/3)  
[http_workers/2](#http_workers/2)  
[3.14.2](#idx:httpworkers2:31)

[http_wrapper/5](#http_wrapper/5)  
[3.12](#idx:httpwrapper5:21) [3.14.4](#idx:httpwrapper5:41) [3.15](#idx:httpwrapper5:42) [3.15](#idx:httpwrapper5:43) [3.15](#idx:httpwrapper5:45)

[http_write_passwd_file/2](#http_write_passwd_file/2)  
[hub_add/3](#hub_add/3)  
[hub_broadcast/2](#hub_broadcast/2)  
[hub_broadcast/3](#hub_broadcast/3)  
[hub_create/3](#hub_create/3)  
[hub_member/2](#hub_member/2)  
[hub_send/2](#hub_send/2)  
[iostream:open_hook/6](#iostream:open_hook/6)  
[javascript/4](#javascript/4)  
[js_arg//1](#js_arg//1)  
[js_arg_list//1](#js_arg_list//1)  
[js_call//1](#js_call//1)  
[js_expression//1](#js_expression//1)  
[js_new//2](#js_new//2)  
[js_script//1](#js_script//1)  
[keep_alive/4](#keep_alive/4)  
[map_method/2](#map_method/2)  
[mime_include//2](#mime_include//2)  
[mime_pack/3](#mime_pack/3)  
mime_type_encoding/2  
[3.1](#idx:mimetypeencoding2:13)

[nolog/1](#nolog/1)  
[nolog_post_content_type/1](#nolog_post_content_type/1)  
[openid_associate/3](#openid_associate/3)  
[openid_associate/4](#openid_associate/4)  
[openid_authenticate/4](#openid_authenticate/4)  
[openid_current_host/3](#openid_current_host/3)  
[openid_current_url/2](#openid_current_url/2)  
[openid_grant/1](#openid_grant/1)  
[openid_hook/1](#openid_hook/1)  
[openid_logged_in/1](#openid_logged_in/1)  
[openid_login/1](#openid_login/1)  
[openid_login_form//2](#openid_login_form//2)  
[openid_logout/1](#openid_logout/1)  
[openid_server/2](#openid_server/2)  
[openid_server/3](#openid_server/3)  
[openid_user/3](#openid_user/3)  
[openid_verify/2](#openid_verify/2)  
[page//1](#page//1)  
[page//2](#page//2)  
page/\[1,2\]  
[3.21](#idx:page12:59)

[password_field/1](#password_field/1)  
[post_data_encoded/2](#post_data_encoded/2)  
pp/1  
[3.13.1](#idx:pp1:28)

predicate/5  
[3.21.5](#idx:predicate5:67)

[print_html/1](#print_html/1)  
[print_html/2](#print_html/2)  
print_html/\[1,2\]  
[3.21](#idx:printhtml12:49) [3.21](#idx:printhtml12:50) [3.21.4](#idx:printhtml12:66)

process_create/3  
[9](#idx:processcreate3:81)

[pwp_handler/2](#pwp_handler/2)  
[reply_html_page/2](#reply_html_page/2)  
[reply_html_page/3](#reply_html_page/3)  
[3.21.1](#idx:replyhtmlpage3:61) [3.21.1](#idx:replyhtmlpage3:62) [3.21.1](#idx:replyhtmlpage3:63)

[reply_html_partial/1](#reply_html_partial/1)  
[reply_htmx/1](#reply_htmx/1)  
[reply_htmx/2](#reply_htmx/2)  
[reply_pwp_page/3](#reply_pwp_page/3)  
[server_health/1](#server_health/1)  
[session_setting/2](#session_setting/2)  
set_lang/1  
[3.21](#idx:setlang1:56) [3.21](#idx:setlang1:58)

set_stream/2  
[3.1](#idx:setstream2:12)

setup_call_cleanup/3  
[2](#idx:setupcallcleanup3:4)

shell/1  
[9](#idx:shell1:80) [9](#idx:shell1:82)

sleep/1  
[3.9.1](#idx:sleep1:18)

ssl_context/3  
[3.14.2](#idx:sslcontext3:32)

[ssl_verify/5](#ssl_verify/5)  
tcp_accept/3  
[3.15](#idx:tcpaccept3:44)

tcp_bind/2  
[3.14.2](#idx:tcpbind2:30)

tcp_connect/3  
[4](#idx:tcpconnect3:72)

thread_create/3  
[3.14.2](#idx:threadcreate3:38)

thread_create_in_pool/4  
[3.14.2](#idx:threadcreateinpool4:37)

thread_pool_create/3  
[3.14.2](#idx:threadpoolcreate3:36)

thread_wait/2  
[3.9.1](#idx:threadwait2:19)

throw/1  
[3.1.1](#idx:throw1:14)

tspy/1  
[3.14](#idx:tspy1:29) [10](#idx:tspy1:83)

uri_components/2  
[10](#idx:uricomponents2:85)

uri_data/4  
[10](#idx:uridata4:86)

uri_edit/3  
[10](#idx:uriedit3:87)

uri_encoded/3  
[3.21](#idx:uriencoded3:54)

uri_nomalized/2  
[10](#idx:urinomalized2:88)

[ws_close/3](#ws_close/3)  
[ws_mask/1](#ws_mask/1)  
[ws_open/3](#ws_open/3)  
[ws_property/2](#ws_property/2)  
[ws_receive/2](#ws_receive/2)  
[ws_receive/3](#ws_receive/3)  
[ws_send/2](#ws_send/2)  
[xhtml_ns//2](#xhtml_ns//2)  
