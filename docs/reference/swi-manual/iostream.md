
## A.23 library(iostream): Utilities to deal with streams

See also  
`library(archive)`, `library(process)`, `library(zlib)`, `library(http/http_stream)`

This library contains utilities that deal with streams, notably originating from non-built-in sources such as URLs, archives, windows, processes, etc.

The predicate [open_any/5](iostream.html#open_any/5) acts as a *broker* between applications that can process data from a stream and libraries that can create streams from diverse sources. Without this predicate, processing data inevitally follows the pattern below. As *call_some_open_variation* can be anything, this blocks us from writing predicates such as `load_xml(From, DOM)` that can operate on arbitrary input sources.

``` code
setup_call_cleanup(
    call_some_open_variation(Spec, In),
    process(In),
    close(In)).
```

Libraries that can open streams can install the hook iostream:open_hook/6 to make their functionality available through [open_any/5](iostream.html#open_any/5).

**open_any**(`+Specification, +Mode, -Stream, -Close, +Options`)  
Establish a stream from `Specification` that should be closed using `Close`, which can either be called or passed to [close_any/1](iostream.html#close_any/1). `Options` processed:

**encoding**(`Enc`)  
Set stream to encoding `Enc`.

Without loaded plugins, the [open_any/5](iostream.html#open_any/5) processes the following values for `Specification`. If no rule matches, [open_any/5](iostream.html#open_any/5) processes `Specification` as `file(Specification)`.

**`Stream`**  
A plain stream handle. Possisible post-processing options such as encoding are applied. `Close` does *not* close the stream, but resets other side-effects such as the encoding.

**stream**(`Stream`)  
Same as a plain `Stream`.

**`FileURL`**  
If `Specification` is of the form =file://...=, the pointed to file is opened using [open/4](IO.html#open/4). Requires `library(uri)` to be installed.

**file**(`Path`)  
Explicitly open the file `Path`. `Path` can be an `Path`(File) term as accepted by [absolute_file_name/3](files.html#absolute_file_name/3).

**string**(`String`)  
Open a Prolog string, atom, list of characters or codes as an *input* stream.

The typical usage scenario is given in the code below, where `<`process`>` processes the input.

``` code
setup_call_cleanup(
    open_any(Spec, read, In, Close, Options),
    <process>(In),
    Close).
```

Currently, the following libraries extend this predicate:

**library**(`http/http_open`)  
Adds support for URLs using the `http` and `https` schemes.

**close_any**(`+Goal`)  
Execute the `Close` closure returned by [open_any/5](iostream.html#open_any/5). The closure can also be called directly. Using [close_any/1](iostream.html#close_any/1) can be considered better style and enhances tractability of the source code.

\[semidet,multifile\]**open_hook**(`+Spec, +Mode, -Stream, -Close, +Options0, -Options`)  
Open `Spec` in `Mode`, producing `Stream`.

|  |  |
|----|----|
| `Close` | is unified to a goal that must be called to undo the side-effects of the action, e.g., typically the term `close(Stream)` |
| `Options0` | are the options passed to [open_any/5](iostream.html#open_any/5) |
| `Options` | are passed to the post processing filters that may be installed by [open_any/5](iostream.html#open_any/5). |
