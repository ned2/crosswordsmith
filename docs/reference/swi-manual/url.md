
## A.64 library(url): Analysing and constructing URL

author  
\- Jan Wielemaker  
- Lukas Faulstich

deprecated  
New code should use `library(uri)`, provided by the `clib` package.

This library deals with the analysis and construction of a URL, Universal Resource Locator. URL is the basis for communicating locations of resources (data) on the web. A URL consists of a protocol identifier (e.g. HTTP, FTP, and a protocol-specific syntax further defining the location. URLs are standardized in RFC-1738.

The implementation in this library covers only a small portion of the defined protocols. Though the initial implementation followed RFC-1738 strictly, the current is more relaxed to deal with frequent violations of the standard encountered in practical use.

\[det\]**global_url**(`+URL, +Base, -Global`)  
Translate a possibly relative `URL` into an absolute one.

Errors  
`syntax_error(illegal_url)` if `URL` is not legal.

**is_absolute_url**(`+URL`)  
True if `URL` is an absolute `URL`. That is, a `URL` that starts with a protocol identifier.

**http_location**(`?Parts, ?Location`)  
Construct or analyze an HTTP location. This is similar to [parse_url/2](url.html#parse_url/2), but only deals with the location part of an HTTP URL. That is, the path, search and fragment specifiers. In the HTTP protocol, the first line of a message is

``` code
<Action> <Location> HTTP/<version>
```

|            |                                  |
|------------|----------------------------------|
| `Location` | Atom or list of character codes. |

\[det\]**parse_url**(`?URL, ?Attributes`)  
Construct or analyse a `URL`. `URL` is an atom holding a `URL` or a variable. `Attributes` is a list of components. Each component is of the format Name(Value). Defined components are:

**protocol**(`Protocol`)  
The used protocol. This is, after the optional `url:`, an identifier separated from the remainder of the `URL` using :. [parse_url/2](url.html#parse_url/2) assumes the `http` protocol if no protocol is specified and the `URL` can be parsed as a valid HTTP url. In addition to the RFC-1738 specified protocols, the `file` protocol is supported as well.

**host**(`Host`)  
`Host`-name or IP-address on which the resource is located. Supported by all network-based protocols.

**port**(`Port`)  
Integer port-number to access on the `\`arg{Host}. This only appears if the port is explicitly specified in the `URL`. Implicit default ports (e.g., 80 for HTTP) do *not* appear in the part-list.

**path**(`Path`)  
(File-) path addressed by the `URL`. This is supported for the `ftp`, `http` and `file` protocols. If no path appears, the library generates the path `/`.

**search**(`ListOfNameValue`)  
Search-specification of HTTP `URL`. This is the part after the `?`, normally used to transfer data from HTML forms that use the HTTP GET method. In the `URL` it consists of a www-form-encoded list of Name=Value pairs. This is mapped to a list of Prolog Name=Value terms with decoded names and values.

**fragment**(`Fragment`)  
`Fragment` specification of HTTP `URL`. This is the part after the `#` character.

The example below illustrates all of this for an HTTP `URL`.

``` code
?- parse_url('http://www.xyz.org/hello?msg=Hello+World%21#x',
       P).

P = [ protocol(http),
      host('www.xyz.org'),
      fragment(x),
      search([ msg = 'Hello World!'
             ]),
      path('/hello')
    ]
```

By instantiating the parts-list this predicate can be used to create a `URL`.

\[det\]**parse_url**(`+URL, +BaseURL, -Attributes`)  
Similar to [parse_url/2](url.html#parse_url/2) for relative URLs. If `URL` is relative, it is resolved using the absolute `URL` `BaseURL`.

\[det\]**www_form_encode**(`+Value, -XWWWFormEncoded`)  
\[det\]**www_form_encode**(`-Value, +XWWWFormEncoded`)  
En/decode to/from application/x-www-form-encoded. Encoding encodes all characters except RFC 3986 *unreserved* (ASCII `alnum` (see [code_type/2](chartype.html#code_type/2))), and one of "-.\_`~`" using percent encoding. Newline is mapped to `%OD%OA`. When decoding, newlines appear as a single newline (10) character.

Note that a space is encoded as `%20` instead of `+`. Decoding decodes both to a space.

deprecated  
Use uri_encoded/3 for new code.

\[semidet\]**set_url_encoding**(`?Old, +New`)  
Query and set the encoding for URLs. The default is `utf8`. The only other defined value is `iso_latin_1`.

To be done  
Having a global flag is highly inconvenient, but a work-around for old sites using ISO Latin 1 encoding.

\[det\]**url_iri**(`+Encoded, -Decoded`)  
\[det\]**url_iri**(`-Encoded, +Decoded`)  
Convert between a URL, encoding in US-ASCII and an IRI. An IRI is a fully expanded Unicode string. Unicode strings are first encoded into UTF-8, after which %-encoding takes place.

\[det\]**parse_url_search**(`?Spec, ?Fields:list(Name=Value)`)  
Construct or analyze an HTTP search specification. This deals with form data using the MIME-type `application/x-www-form-urlencoded` as used in HTTP GET requests.

\[det\]**file_name_to_url**(`+File, -URL`)  
\[semidet\]**file_name_to_url**(`-File, +URL`)  
Translate between a filename and a file:`//` `URL`.

To be done  
Current implementation does not deal with paths that need special encoding.
