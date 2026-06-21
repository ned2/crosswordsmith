
## B.10 Hooks for loading files

All loading of source files is achieved by [load_files/2](consulting.html#load_files/2). The hook [prolog_load_file/2](loadfilehook.html#prolog_load_file/2) can be used to load Prolog code from non-files or even load entirely different information, such as foreign files.

**prolog_load_file**(`+Spec, +Options`)  
Load a single object. If this call succeeds, [load_files/2](consulting.html#load_files/2) assumes the action has been taken care of. This hook is only called if `Options` does not contain the `stream(Input)` option. The hook must be defined in the module `user`.

This can be used to load from unusual places as well as dealing with Prolog code that is not represented as a Prolog source text (for example some binary representation). For example, library `library(http/http_load)` loads Prolog directly from an HTTP server. See also [prolog:open_source_hook/3](loadfilehook.html#prolog:open_source_hook/3), which merely allows for changing how a physical file is opened.

**prolog:open_source_hook**(`+Path, -Stream, +Options`)  
This hooks is called by the compiler to overrule the default [open/3](IO.html#open/3) call `open(Path, read, Stream)`. `Options` provide the options as provided to [load_files/2](consulting.html#load_files/2). If the hook succeeds compilation continues by loading from the returned (input) stream. This hook is particularly suited to support running the code to a preprocessor. See also [prolog_load_file/2](loadfilehook.html#prolog_load_file/2).

**prolog:comment_hook**(`+Comments, +Pos, +Term`)  
This hook allows for processing comments encountered by the compiler. If this hook is defined, the compiler calls [read_term/2](termrw.html#read_term/2) with the option `comments(Comments)`. If the list of comments returned by [read_term/2](termrw.html#read_term/2) is not empty it calls this comment hook with the following arguments.

- `Comments` is the non-empty list of comments. Each comment is a pair `Position`-`String`, where `String` is a string object (see [section 5.2](string.html#sec:5.2)) that contains the comment *including* delimiters. Consecutive line comments are returned as a single comment.
- `Pos` is a stream-position term that describes the starting position of `Term`
- `Term` is the term read.

This hook is exploited by the documentation system. See [stream_position_data/3](IO.html#stream_position_data/3). See also [read_term/3](termrw.html#read_term/3).
