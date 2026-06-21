
## 12.9 Foreign access to Prolog IO streams

The SWI-Prolog foreign language interface provides access to Prolog IO streams. This interface may be used to get hold of Prolog streams for reading and writing. In addition, this interface allows to define new stream types. For example, the Windows **swipl-win.exe** executable that runs Prolog in a Windows GUI redefines the Prolog standard IO streams (`user_input`, `user_output` and `user_error` to read from and write to the GUI window.

The interface is built around the `IOSTREAM` type which plays a role similar to the POSIX `FILE` type. Most of the functions are modeled after their `FILE` counterpart, prefixed by **S**, e.g. [Sfwrite()](foreign-streams.html#Sfwrite()). The `IOSTREAM` type has considerably more features though. The `IOSTREAM` type is practically disconnected from the rest of the Prolog system. Prolog refers to streams either by *alias* (`user_input`, etc. or created using the `alias(Name)` option of [open/4](IO.html#open/4)) or using a *stream handle* which is represented as a *blob* (see [section 12.4.10](foreigninclude.html#sec:12.4.10)). Foreign extensions that wish to access or define streams should include `SWI-Stream.h` in addition to `SWI-Prolog.h` as below. Both headers may be used with C as well as C++.

The interface also defines `Sinput`, `Suser`, `Serror` for direct access to the operating system's input and output streams, bypassing Prolog's control - for example, these will not be affected by with_output_to/3. There is also a convenience function for debugging, which goes directly to `stderr`: [Sdprintf()](foreign-streams.html#Sdprintf()).^(243On Windows the output is also emitted using **OutputDebugString()**.)

``` code
#include <SWI-Stream.h>
#include <SWI-Prolog.h>
```

### 12.9.1 Get IO stream handles

There are several ways to get access to an IO Stream handle, basically get them from Prolog, get access to the standard streams and create a new stream. The *standard streams* are available as `Sinput`, `Soutput` and `Serror`. Note that these are thread specific. Creating a new stream is discussed with [Snew()](foreign-streams.html#Snew()). Below are the functions to obtain a stream handle from a Prolog term, obtain and release ownership.

`int` **PL_get_stream**(`term_t t, IOSTREAM **s, int flags`)  
Get a stream handle from the Prolog term `t`. Returns `TRUE` on success and `FALSE` on failure, by default generating an exception. The `flags` argument is a bitwise disjunction of these flags:

**`SIO_INPUT`**  
Get an *input stream*. If `t` is a stream pair (see [stream_pair/3](IO.html#stream_pair/3)), return the input channel. If `t` is an output stream the function fails.

**`SIO_OUTPUT`**  
Get an *output stream*. See `SIO_INPUT` for details. If neither `SIO_OUTPUT` nor `SIO_INPUT` is given `t` may not be a *pair*.

**`SIO_TRYLOCK`**  
Return `FALSE` if the stream cannot be locked immediately. No error is generated.

**`SIO_NOERROR`**  
If the function fails no exception is produced.

The returned stream is owned by the calling thread using [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()).

`int` **PL_get_stream_from_blob**(`atom_t b, IOSTREAM **s, int flags`)  
Same as [PL_get_stream()](foreign-streams.html#PL_get_stream()), but operates directly on the blob `b`. This allows for foreign code that wishes long term access to a stream to maintain a handle to the stream as a (registered) `atom_t` object rather than a `IOSTREAM*`.

`IOSTREAM *` **PL_acquire_stream**(`IOSTREAM *s`)  
Obtain ownership of `s` and return `s`. The application must call [PL_release_stream()](foreign-streams.html#PL_release_stream()) when done. Only one thread can own a stream and this call blocks if some other thread owns the stream. This function may be called multiple times by the same thread (*recursive lock*). Note that [PL_get_stream()](foreign-streams.html#PL_get_stream()) also acquires ownership.

`int` **PL_release_stream**(`IOSTREAM *s`)  
Give up ownership acquired using [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()) or [PL_get_stream()](foreign-streams.html#PL_get_stream()). If the stream is an an error state, return `FALSE` with an exception. Otherwise return `TRUE`.

In general, stream functions do not set any Prolog error state; that is done by [PL_release_stream()](foreign-streams.html#PL_release_stream()). Once a stream is in an error state, all subsequent functions act as no-ops (returning -1) unless [Sclearerr()](foreign-streams.html#Sclearerr()) is called. [Sferror()](foreign-streams.html#Sferror()) may be used to check whether a stream is in an error condition. This error may be turned into a Prolog exception by calling [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()) followed by [PL_release_stream()](foreign-streams.html#PL_release_stream()). In this case, [PL_release_stream()](foreign-streams.html#PL_release_stream()) will set the Prolog exception and return `FALSE`.

Below is an example that writes “Hello World” to a stream provided by Prolog. Note that [PL_release_stream()](foreign-streams.html#PL_release_stream()) raises an exception if the [Sfprintf()](foreign-streams.html#Sfprintf()) failed and (thus) left the stream in an error state.

``` code
static foreign_t
hello_world(term_t to)
{ IOSTREAM *s;

  if ( PL_get_stream(to, &s, SIO_OUTPUT) )
  { Sfprintf(s, "Hello World!\n");
    return PL_release_stream(s);
  }

  return FALSE;
}

  ... // fragment from install function
  PL_register_foreign("hello world", 1, hello_world, 0);
```

### 12.9.2 Creating an IO stream

A new stream is created using [Snew()](foreign-streams.html#Snew()). Before we can create a stream we must create a function block of type `IOFUNCTIONS` that provide function pointers for the basic operations on the stream. This type is defined as follows:

``` code
typedef struct io_functions
{ Sread_function        read;           /* fill the buffer */
  Swrite_function       write;          /* empty the buffer */
  Sseek_function        seek;           /* seek to position */
  Sclose_function       close;          /* close stream */
  Scontrol_function     control;        /* Info/control */
  Sseek64_function      seek64;         /* seek to position (large files) */
} IOFUNCTIONS;
```

`ssize_t` **(\*Sread_function)**(`void *handle, char *buf, size_t bufsize`)  
Read new data into `buf` that has size `bufsize`, return the number of bytes read or -1. Note that this is the same interface as the POSIX **read()** API. See [section 12.9.4](foreign-streams.html#sec:12.9.4) for raising errors.

`ssize_t` **(\*Swrite_function)**(`void *handle, char *buf, size_t bufsize`)  
Write the bytes from `buf` with contains `bufsize` bytes and return the number of bytes written or -1. The number of bytes written may be less than `bufsize`. Bytes that were not written remain in the stream's output buffer. Note that this is the same interface as the POSIX [write()](foreigninclude.html#write()) API. See [section 12.9.4](foreign-streams.html#sec:12.9.4) for raising errors.

`long` **(\*Sseek_function)**(`void *handle, long pos, int whence`)  
`int64_t` **(\*Sseek64_function)**(`void *handle, int64_t pos, int whence`)  
Reposition the file pointer. These functions may be `NULL` if repositioning is not possible on this type or they may return -1 and set `errno` to `EPIPE` if the pointer cannot be repositioned on this instance. The function returns the new file position. See [Sseek()](foreign-streams.html#Sseek()) for details on how repositioning is implemented. See [section 12.9.4](foreign-streams.html#sec:12.9.4) for raising errors.

`int` **(\*Sclose_function)**(`void *handle`)  
Close the stream. This is used by [Sclose()](foreign-streams.html#Sclose()). Note that buffered output is first written using the **Swrite_function()**. See [section 12.9.4](foreign-streams.html#sec:12.9.4) for raising errors.

`int` **(\*Scontrol_function)**(`void *handle, int action, void *arg`)  
Obtain information about the stream or modify the stream. The function should return 0 on success and -1 on failure. If some action is not implemented the function should return -1;

**`SIO_GETPENDING`, `size_t*`**  
Return the number of bytes that may be written without blocking. Used by [Spending()](foreign-streams.html#Spending()).

**`SIO_LASTERROR`, `char*`**  
Called after an error is raised on a stream. May return a C string that sets error details using [Sseterr()](foreign-streams.html#Sseterr()).

**`SIO_SETENCODING`, `IOENC*`**  
Called by [Ssetenc()](foreign-streams.html#Ssetenc()) to change the encoding of the stream. If the call does not return 0 the encoding is not changed.

**`SIO_FLUSHOUTPUT`, `NULL`**  
Called by [Sflush()](foreign-streams.html#Sflush()) after flushing the stream's output buffer. Note that this is only called on an *explicit* flush using [Sflush()](foreign-streams.html#Sflush()) or [flush_output/1](chario.html#flush_output/1). An implicit flush because the output buffer is full does *not* call this hook.

**`SIO_GETSIZE`, `int64_t*`**  
Get the size of the underlying object in bytes. Used by [Ssize()](foreign-streams.html#Ssize()).

**`SIO_GETFILENO`, `int*`**  
If the stream is associated with an OS file handle, return this handle. Used by [Sfileno()](foreign-streams.html#Sfileno()).

**`SIO_GETWINSOCK`, `SOCKET*`**  
Windows only. If the stream is associated to a Windows socket return this handle. Used by [Swinsock()](foreign-streams.html#Swinsock()).

Given an `IOFUNCTIONS` block we can create a new stream from a `handle` using [Snew()](foreign-streams.html#Snew()):

`IOSTREAM*` **Snew**(`void *handle, int flags, IOFUNCTIONS *functions`)  
Create an `IOSTREAM*` from a handle, flags and a block of callback functions. The `flags` argument is a bitwise or of SIO\_\* flags. Flags that control the creation are:

**`SIO_INPUT`**  
**`SIO_OUTPUT`**  
One of these flags mut be present to indicate whether this is an input or output stream.

**`SIO_NBUF`**  
**`SIO_LBUF`**  
**`SIO_FBUF`**  
One of these flags must be present to select the buffering as one of unbuffered (`SIO_NBUF`), line buffered (`SIO_LBUF`) or fully buffered (`SIO_FBUF`)

**`SIO_TEXT`**  
If given, this is a text stream and the encoding is set to the default encoding (see the Prolog flag [encoding](flags.html#flag:encoding)). Otherwise this is a binary stream and the encoding is set to `ENC_OCTET`.

**`SIO_RECORDPOS`**  
If given, enable position maintenance on the stream. This is used by [Stell()](foreign-streams.html#Stell()), [Sseek()](foreign-streams.html#Sseek()), [stream_property/2](IO.html#stream_property/2) using the `position` property and related predicates.

**`SIO_NOMUTEX`**  
Used internally to create a stream that cannot be owned or locked.

If the stream is associated with an OS file handle the system initializes the `SIO_ISATTY` flag (on POSIX systems) and if possible tells the OS not to inherit this stream to child processes.

The symbol `Sfilefunctions` is a `IOFUNCTIONS` struct that contains the callbacks for accessing a regular file. After opening an file using the POSIX **open()** API we can create a stream to this file using [Snew()](foreign-streams.html#Snew()):

``` code
  int fno = open(path, O_RDONLY);
  IOSTREAM *s;

  if ( fno >= 0 )
    s = Snew((void*)fno,
             SIO_INPUT|SIO_FBUF|SIO_RECORDPOS|SIO_TEXT,
             &Sfilefunctions);
  ...
```

[Snew()](foreign-streams.html#Snew()) can only fail if there is not enough memory. In that case the return value is `NULL` and `errno` is set to `ENOMEM`.

`IOSTREAM*` **Sopen_pipe**(`const char *command, const char *type`)  
Start a process from `command` and connect the input or output to the returned stream. This wraps the POSIX **popen()** API. The `type` string starts with `r` or `w` and may be followed by `b` to create a *binary stream*. The default is to create a text stream using the platform conventions and locale.

`IOSTREAM*` **Sopenmem**(`char **buffer, size_t *sizep, const char *mode`)  
Open a memory area as a stream. Output streams are automatically resized using **realloc()** if \*`size` = 0 or the stream is opened with mode `"wa"`. If the buffer is allocated or enlarged, this is achieved using **malloc()** or **realloc()**. In this case the returned buffer should be freed by the caller when done. Example:

``` code
    { char buf[1024];             // don't allocate for small stuff
      char *s = buf;
      IOSTREAM *fd;
      size_t size = sizeof(buf);

      fd = Sopenmem(&s, &size, "w");
      ...
      Sclose(fd);
      ...
      if ( s != buf )             // apparently moved
        Sfree(s);
    }
```

The `mode` is `"r"` or `"w"`. The mode `"rF"` calls [`PL_free(buffer)`](foreignnotes.html#PL_free()) when closed.

**Note:** Its is *not* allowed to access streams created with this call from multiple threads. This is ok for all usage inside Prolog itself. This call is intended to use [Sfprintf()](foreign-streams.html#Sfprintf()) and other output functions to create strings.

`void` **Sfree**(`void *ptr`)  
This function must be used to free objects that are allocated by the stream interface. Currently this only applies to strings allocated by [Sopenmem()](foreign-streams.html#Sopenmem()).

A stream can be made accessible from Prolog using [PL_unify_stream()](foreign-streams.html#PL_unify_stream()):

`int` **PL_unify_stream**(`term_t t, IOSTREAM *s`)  
Unify `t` with a *blob* that points at `s`. Note that a blob provides a unique and reliable reference to a stream. Blobs are subject to *atom garbage collection*. If an open stream is garbage collected the behaviour depends on the Prolog flag [agc_close_streams](flags.html#flag:agc_close_streams). See also [Sgcclose()](foreign-streams.html#Sgcclose()).

### 12.9.3 Interacting with foreign streams

`int` **Sset_timeout**(`IOSTREAM *s, int milliseconds`)  
Set the timeout on an input stream to `milliseconds`. If this value is non-negative the the **poll()** or **select()** API is used to wait until input is available. If no input is available within the specified time an error is raised on the stream.

`int` **Sunit_size**()  
Returns the size of a code unit in bytes depending on the stream's encoding. This returns 2 for the encodings `ENC_UNICODE_BE` and `ENC_UNICODE_LE`, `sizeof(wchar_t)` for `ENC_WCHAR` and 1 for all other encodings (including multibyte encodings such as `ENC_UTF8`.

`int` **Sputc**(`int c, IOSTREAM *s`)  
Emit a byte to `s`. Flushes the buffer on `\n` when in `SIO_LBUF` buffering mode and updates the stream position information if enabled (`SIO_RECORDPOS`). Returns 0 on success, -1 on error.

`int` **Sgetc**(`IOSTREAM *s`)  
Read a byte from `s`. Fills the input buffer if buffering is enabled and the buffer is empty. Updates the stream position information if enabled (`SIO_RECORDPOS`). Returns -1 on end of file or error. Use [Sferror()](foreign-streams.html#Sferror()) or [Sfeof()](foreign-streams.html#Sfeof()) to distinguish end of file from an error. This is a C macro.

`int` **Sfgetc**(`IOSTREAM *s`)  
Function equivalent to [Sgetc()](foreign-streams.html#Sgetc()).

`int` **Sungetc**(`int c, IOSTREAM *s`)  
Put a byte back into the input buffer. Returns -1 if this is not possible. Deprecated. New code should use [Speekcode()](foreign-streams.html#Speekcode()) because that reliably maintains the position information on the stream.

`int` **Sputcode**(`int c, IOSTREAM *s`)  
Emit a Unicode code point to `s`. This function also performs newline encoding (see [section 12.9.6](foreign-streams.html#sec:12.9.6)). If the encoding of `s` cannot represent `c`, the behaviour depends on the the following flags. Only one of these flags may be enabled. If none of these flags is enabled an error is raised and the function returns -1.

**`SIO_REPXML`**  
Emit as XML character entity, e.g. `&#4242;`

**`SIO_REPPL`**  
Emit as ISO escape, e.g., `\x4242\`

**`SIO_REPPLU`**  
Emit as Unicode escape, e.g., `\u4242` or `\U42424242`

Updates the stream position information if enabled (`SIO_RECORDPOS`)

`int` **Sgetcode**(`IOSTREAM *s`)  
Read a Unicode code point from `s`. If it detects an invalid multibyte character a warning is emitted and the code point `0xfffd` is returned. Other errors and end-of-file return -1; Use [Sferror()](foreign-streams.html#Sferror()) or [Sfeof()](foreign-streams.html#Sfeof()) to distinguish end of file from an error.

`int` **Speekcode**(`IOSTREAM *s`)  
As [Sgetcode()](foreign-streams.html#Sgetcode()), but leaves the character in the input buffer and does not update the stream position. Returns -1 if the stream is not buffered (`SIO_NBUF`).

`int` **Sputw**(`int w, IOSTREAM *s`)  
`int` **Sgetw**(`IOSTREAM *s`)  
Reads/writes an integer in native byte order. Deprecated.

`size_t` **Sfread**(`void *data, size_t size, size_t elems, IOSTREAM *s`)  
`size_t` **Sfwrite**(`const void *data, size_t size, size_t elems, IOSTREAM *s`)  
Emulations of the POSIX **fread()** and **fwrite()** calls for Prolog streams. These functions read or write `elems` objects of size `size` and return the number of objects successfully read or written. Data exchange is binary (even if the stream is in text mode) and unlike **read()** and [write()](foreigninclude.html#write()), these functions keep reading or writing until end-of-file (for [Sfread()](foreign-streams.html#Sfread())) or an error.

`int` **Sfeof**(`IOSTREAM *s`)  
Returns non-zero if the stream is at the end. It performs the following checks: (1) test the `SIO_FEOF` flag, (2) test whether the buffer is non-empty, (3) fill the buffer and return non-zero if the **Sread_function()** returned 0 (zero).

`int` **Sfpasteof**(`IOSTREAM *s`)  
Returns non-zero when a read operation was performed after signalling end-of-file. On other words, reaching end-of-file first triggers [Sfeof()](foreign-streams.html#Sfeof()) and after another read triggers [Sfpasteof()](foreign-streams.html#Sfpasteof()).

`int` **Ssetlocale**(`IOSTREAM *s, struct PL_locale *new_loc, struct PL_locale **old_loc`)  
Change the locale associated with a stream. The current system does not provide a public C API for dealing with Prolog locale objects. See [section 4.23](locale.html#sec:4.23).

`int` **Sflush**(`IOSTREAM *s`)  
Flush buffered output, returning 0 on success and -1 after a (write) error occurred. Calls **Scontrol_function()** using the action `SIO_FLUSHOUTPUT` after the buffer was successfully written.

`int64_t` **Ssize**(`IOSTREAM *s`)  
Returns the size in bytes of the object associated to the stream or -1 if this is not known.

`int` **Sseek**(`IOSTREAM *s, long pos, int whence`)  
Deprecated - use [Sseek64()](foreign-streams.html#Sseek64()) instead because some platforms define `long` as 32-bits.

`int` **Sseek64**(`IOSTREAM *s, int64_t pos, int whence`)  
Reposition the file pointer in the object associated to `s`, returning 0 on success and -1 otherwise. If the stream is buffered and position information is maintained these functions readjust the buffer information if possible. Otherwise they call **Sseek64_function()** or **Sseek_function()** as a fallback iff `pos` can be represented as a C `long`. `Whence` is one of `SIO_SEEK_SET`, `SIO_SEEK_CUR` or `SIO_SEEK_END`, seeking relative to the start, current position or end.

`long` **Stell**(`IOSTREAM *s`)  
Deprecated - use [Stell64()](foreign-streams.html#Stell64()) instead because some platforms define `long` as 32-bits.

`int64_t` **Stell64**(`IOSTREAM *s`)  
Return the current position in the stream. This is obtained from the recorded position or based on information from the seek handlers, adjusted with the buffer information.

`int` **Sclose**(`IOSTREAM *s`)  
Close the stream. This first locks the stream (see [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream())). When successful it flushes pending output and calls the **Sclose_function()** hook. Finally, the stream is unlocked and all memory associated to the stream is released. On success, the function returns 0. On failure a Prolog exception is raised and the return value is -1. Regardless of the return value, `s` becomes invalid after completion of [Sclose()](foreign-streams.html#Sclose()). See also [Sgcclose()](foreign-streams.html#Sgcclose()).

`int` **Sgcclose**(`IOSTREAM *s, int flags`)  
As [Sclose()](foreign-streams.html#Sclose()), but intended to be used from the atom garbage collector if a stream is closed because it is garbage. The SWI-Prolog atom garbage collector normally runs in a separate thread and thus may be unable to obtain a lock on `s` if some thread lost access to the stream while it is locked. For this situation `flags` may be `SIO_CLOSE_TRYLOCK` which causes [Sgcclose()](foreign-streams.html#Sgcclose()) to return -1 with `errno` set to `EDEADLK` if the stream is locked. Alternatively, using `SIO_CLOSE_FORCE` the stream is closed and released without gaining a lock. This should be safe because the stream is garbage and thus no thread can use the lock.

In addition, [Sgcclose()](foreign-streams.html#Sgcclose()) never raises a Prolog exception because Prolog interaction is not allowed from the blob release hook and there is no meaningful way to raise a Prolog exception from this context.

`char*` **Sfgets**(`char *buf, int n, IOSTREAM *s`)  
Read a line of input as a sequence of *bytes*. The `buf` is `n` bytes long. On success, `buf` is returned and contains a 0-terminated C string that ends with a `\n` character. On end-of-file or an error, `NULL` is returned. If the input line is longer that `n` bytes `buf` is **not** 0-terminated.

`int` **Sgets**(`char *buf`)  
Shorthand for [`Sfgets(buf, Slinesize, Sinput)`](foreign-streams.html#Sfgets()). Deletes the terminating `\n` character. `Slinesize` is a global variable that defines the length of the input buffer. Deprecated.

`int` **Sread_pending**(`IOSTREAM *s, char *buf, size_t limit, int flags`)  
Return the data buffered on an input stream. If `flags` includes `SIO_RP_BLOCK`, fill the buffer (possibly blocking) if the buffer is empty. Update the stream position information unless `flags` include `SIO_RP_NOPOS`. This function effectively provides functionality similar to POSIX **read()** on a stream. This function is used by [read_pending_codes/3](chario.html#read_pending_codes/3).

`size_t` **Spending**(`IOSTREAM *s`)  
Return the number of bytes that can be read from `s` without blocking. If there is buffered input, this is the number of bytes buffered. Otherwise it is the result of the **Scontrol_function()** using the action `SIO_GETPENDING`.

`int` **Sfputs**(`const char *q, IOSTREAM *s`)  
Emit a 0-terminated C string. The input string `q` is handled as a sequence of unsigned characters (code points `1 ... 255`.

`int` **Sputs**(`const char *q`)  
Equivalent to [`Sfputs(q, Soutput)`](foreign-streams.html#Sfputs()).

`int` **Sfprintf**(`IOSTREAM *s, const char *fm, ...`)  
Similar to POSIX **fprintf()**. This function largely accepts the same `%` escape sequences. The `%` character is followed by numeric arguments and modifier characters. The generic format of this is described by the regular expression `[+-0 #]*(\d*|\*)(.(\d*|\*))?`. Here, `+` implies right alignment, `-` left alignment, `0` 0-padding and, a space white-space padding and `#` *modified* output. The two optional numerical arguments are separated by a full stop and may be `*` to get them from the argument list. The first numerical argument specifies the field width and the second the precision for floating point numbers.

This sequence is followed by optional type information. For integers this is one of `l` (`long`), `ll` (`long long`) or `z` (`size_t`). For strings this is one of `L` (ISO Latin 1), `U` (UTF-8) or `W` (`wchar_t*`).

Finally we come to the format specifier. This is one of

`%`  
Emit the `%` character itself.

`c`  
Emit a Unicode code point.

`p`  
Emit a pointer.

`d`  

`i`  
Emit a a signed integer as decimal. The `l` (`long`), `ll` (`long long`) or `z` (`size_t`) denote the size.

`o`  

`u`  

`x`  

`X`  
Emit a a unsigned integer as octal, decimal or hexadecimal.

`f`  

`e`  

`E`  

`g`  

`G`  
Emit a `double`.

`s`  
Emit a 0-terminated string.

Unlike the POSIX **fprintf()**, this function, and the related functions ([Svprintf()](foreign-streams.html#Svprintf()), etc.) returns the number of characters written. Due to multibyte encodings the number of bytes written can be more. On error, it returns a negative value; in some cases there is extra information (e.g., in `errno`) but it cannot be relied on.

Each call to [Sfprintf()](foreign-streams.html#Sfprintf()) is atomic in the sense that another thread that calls [Sfprintf()](foreign-streams.html#Sfprintf()) on the same stream will block. If you wish to do a series of print statements without any other thread interleaving, you should call [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()) and use its returned `IOSTREAM*` value, then call [PL_release_stream()](foreign-streams.html#PL_release_stream()) at the end of the print statements.

`int` **SfprintfX**(`IOSTREAM *s, const char *fm, ...`)

Same as [Sfprintf()](foreign-streams.html#Sfprintf()) but doesn't have the format-checking attribute, which can trigger compiler warnings if the format does not match the arguments. This is intended for formats that include extended format specifiers such as `"%Ws"` or `"%Us"`.

`int` **Sprintf**(`const char *fm, ...`)

Similar to [Sfprintf()](foreign-streams.html#Sfprintf()), printing to `Soutput`

`int` **Svprintf**(`IOSTREAM *s, const char *fm, va_list args`)

Variadic argument list version of [Sfprintf()](foreign-streams.html#Sfprintf()).

`int` **Ssprintf**(`char *buf, const char *fm, ...`)

Print to a C string. Deprecated. Use [Ssnprintf()](foreign-streams.html#Ssnprintf()) instead.

`int` **Ssnprintf**(`char *buf, size_t size, const char *fm, ...`)

Print to a C string, emitting a maximum of `size` bytes while ensuring `buf` is 0-terminated. The `buf` is written using UTF-8 encoding. Unlike **snprintf()**, the return value is the number of logical code points written rather than the number of bytes and if the buffer is too small, `-1` is returned rather than the number of bytes that would be written. Future versions may improve compatibility with the POSIX functions.

`int` **SsnprintfX**(`char *buf, size_t size, const char *fm, ...`)

Same as [Ssnprintf()](foreign-streams.html#Ssnprintf()) but doesn't have the format-checking attribute. This is intended for formats that include extended format specifiers such as `"%Ws"` or `"%Us"`.

`int` **Svsprintf**(`char *buf, const char *fm, va_list args`)

Variadic argument list version of [Ssprintf()](foreign-streams.html#Ssprintf()). Deprecated. Use [Svsnprintf()](foreign-streams.html#Svsnprintf()) instead.

`int` **Svsnprintf**(`char *buf, size_t size, const char *fm, va_list args`)

Variadic argument list version of [Ssnprintf()](foreign-streams.html#Ssnprintf()).

`int` **Sdprintf**(`const char *fm, ...`)

Print to `Serror`. This function should be used for printing debug output from foreign code.

`int` **SdprintfX**(`const char *fm, ...`)

Same as [Sdprintf()](foreign-streams.html#Sdprintf()) but doesn't have the format-checking attribute. This is intended for formats that include extended format specifiers such as `"%Ws"` and `"%Us"`.

`int` **Svdprintf**(`const char *fm, va_list args`)

Variadic argument list version of [Sdprintf()](foreign-streams.html#Sdprintf()).

`int` **Slock**(`IOSTREAM *s`)

`int` **StryLock**(`IOSTREAM *s`)

`int` **Sunlock**(`IOSTREAM *s`)

Low level versions that perform only the (un)locking part of [PL_acquire_stream()](foreign-streams.html#PL_acquire_stream()) and [PL_release_stream()](foreign-streams.html#PL_release_stream()).

`int` **Sfileno**(`IOSTREAM *s`)

If the stream is associated to a POSIX file handle, return this handle. Returns -1 otherwise.

`SOCKET` **Swinsock**(`IOSTREAM *s`)

Windows only. If the stream is associated to a Windows socket handle, returns this handle. Otherwise return `INVALID_SOCKET`

`int` **Sclosehook**(`void (*hook)(IOSTREAM *s)`)

Register a hook function to be called by [Sclose()](foreign-streams.html#Sclose()) just before the stream is deallocated. This is used internally to update the Prolog administration of open streams on [Sclose()](foreign-streams.html#Sclose()).

`int` **Sset_filter**(`IOSTREAM *parent, IOSTREAM *filter`)

Register `filter` as a stream that reads from or writes to the stream `parent`.

`void` **Ssetbuffer**(`IOSTREAM *s, char *buf, size_t size`)

Set the input or output buffer for `s` to `size`. The `buf` argument is either `NULL`, asking the system to allocate a buffer or points at a buffer of (at least) the indicated size long. The default buffer size is defined by the C macro `SIO_BUFSIZE`

#### 12.9.3.1 Writing Prolog terms to foreign streams

`int` **PL_write_term**(`IOSTREAM *s, term_t term, int precedence, int flags`)  
Write `term` to `s`. `precedence` is the initial operator precedence, typically 1200. `flags` is a bitwise or of the constants below. These flags map to options for [write_term/2](termrw.html#write_term/2).

**`PL_WRT_QUOTED`**  
**`PL_WRT_IGNOREOPS`**  
**`PL_WRT_NUMBERVARS`**  
**`PL_WRT_PORTRAY`**  
**`PL_WRT_CHARESCAPES`**  
**`PL_WRT_NO_CHARESCAPES`**  
The `PL_WRT_NO_CHARESCAPES` does not map to a [write_term/2](termrw.html#write_term/2) option. If one of `PL_WRT_CHARESCAPES` or `PL_WRT_NO_CHARESCAPES` is specified, character escapes are (not) applied. If neither is specified the default depends, like for [write/1](termrw.html#write/1), on the [character_escapes](flags.html#flag:character_escapes) flag on the module `user`.^(244Prior to version 9.1.6 the default (no flag) was to escape the quotes and the backslash (`\`).)

**`PL_WRT_BACKQUOTED_STRING`**  
**`PL_WRT_ATTVAR_IGNORE`**  
**`PL_WRT_ATTVAR_DOTS`**  
**`PL_WRT_ATTVAR_WRITE`**  
**`PL_WRT_ATTVAR_PORTRAY`**  
**`PL_WRT_BLOB_PORTRAY`**  
**`PL_WRT_NO_CYCLES`**  
**`PL_WRT_NEWLINE`**  
**`PL_WRT_VARNAMES`**  
**`PL_WRT_BACKQUOTE_IS_SYMBOL`**  
**`PL_WRT_DOTLISTS`**  
**`PL_WRT_BRACETERMS`**  
**`PL_WRT_NODICT`**  
**`PL_WRT_NODOTINATOM`**  
**`PL_WRT_NO_LISTS`**  
**`PL_WRT_RAT_NATURAL`**  
**`PL_WRT_CHARESCAPES_UNICODE`**  
**`PL_WRT_QUOTE_NON_ASCII`**  
**`PL_WRT_PARTIAL`**  

For example, to print a term to `user_error` as the toplevel does, use

``` code
    PL_write_term(Suser_error, t, 1200,
                  PL_WRT_QUOTED|PL_WRT_PORTRAY|
                  PL_WRT_VARNAMES|PL_WRT_NEWLINE)
```

### 12.9.4 Foreign stream error handling

`int` **Sferror**(`IOSTREAM *s`)  
Returns `TRUE` if the stream is in an error condition, `FALSE` if the stream is valid and in normal condition and -1 if the stream is invalid.

`void` **Sclearerr**(`IOSTREAM *s`)  
Clear the error state of a stream. This includes the end-of-file state, pending warnings and errors and timeout.

`int` **Sseterr**(`IOSTREAM *s, int which, const char *message`)  
Set an error or warning state on the stream. The `which` argument is one of `SIO_WARN` or `SIO_FERR`. This causes [PL_release_stream()](foreign-streams.html#PL_release_stream()) to print a message (`SIO_WARN`) or raise an exception (`SIO_FERR`).

`int` **Sset_exception**(`IOSTREAM *s, term_t ex`)  
Associate a Prolog exception term with the stream or clear the associated exception if `ex` is 0 and set/clear the `SIO_FERR` condition on the stream. If an exception is assocated [PL_release_stream()](foreign-streams.html#PL_release_stream()) raises this exception.

### 12.9.5 Foreign stream encoding

`IOSTREAM` has a field `encoding` that is managed at initialization from `SIO_TEXT`. The available encodings are defined as a C *enum* as below.

``` code
typedef enum
{ ENC_UNKNOWN = 0,                      /* invalid/unknown */
  ENC_OCTET,                            /* raw 8 bit input */
  ENC_ASCII,                            /* US-ASCII (0..127) */
  ENC_ISO_LATIN_1,                      /* ISO Latin-1 (0..256) */
  ENC_ANSI,                             /* default (multibyte) codepage */
  ENC_UTF8,
  ENC_UNICODE_BE,                       /* big endian unicode file */
  ENC_UNICODE_LE,                       /* little endian unicode file */
  ENC_WCHAR                             /* wchar_t */
} IOENC;
```

*Binary* streams always have the encoding `ENC_OCTET`. The default encoding of a text stream depends on the Prolog flag [encoding](flags.html#flag:encoding). The encoding is used by all functions that perform text I/O on a stream. The encoding can be changed at any moment using [Ssetenc()](foreign-streams.html#Ssetenc()) which is available from Prolog using the [set_stream/2](IO.html#set_stream/2) `encoding(Encoding)` property. Functions that explicitly manage the encoding are:

`int` **Ssetenc**(`IOSTREAM *s, IOENC new_enc, IOENC *old_enc`)  
Set the encoding for `s` to `new_enc` and, if `old_enc` is not `NULL`, return the old encoding. This function may fail, returning -1 if the **Scontrol_function()** of the stream returns -1 on the `SIO_SETENCODING` request. On success it returns 0. If `new_enc` is `ENC_OCTET` the stream is switched to binary mode. Otherwise text mode is enabled.

`int` **ScheckBOM**(`IOSTREAM *s`)  
This function may be called on a buffered input stream immediately after opening the stream. If the stream starts with a known *Byte Order Mark* (BOM) the encoding is set accordingly and the flag `SIO_BOM` is set on the stream. Possibly resulting encodings are `ENC_UTF8`, `ENC_UNICODE_BE` and `ENC_UNICODE_LE`.

`int` **SwriteBOM**(`IOSTREAM *s`)  
This function writes a *Byte Order Mark* (BOM) to `s` and should be called immediately after opening a stream for writing. If the encoding is one of `ENC_UTF8`, `ENC_UNICODE_BE` or `ENC_UNICODE_LE` it writes the code point `\ufeff` (a zero-width white space) to the stream in the current encoding and sets the `SIO_BOM` flag on the stream.

`int` **Scanrepresent**(`int c, IOSTREAM *s`)  
Returns 0 if the encoding of `s` can represent the code point `c` and -1 otherwise.

### 12.9.6 Foreign stream line endings

Text streams have a field `newline` that controls the handling of the newline convention. Note that inside Prolog all lines end with a single newline (`\u000a`, `\n`) code point. The values are described below. The default depends on the OS and can be manipulated using the `newline(Mode)` property of [set_stream/2](IO.html#set_stream/2).

**`SIO_NL_DETECT`**  
This mode may be enabled on an input stream. It causes the stream to read up to the first newline to set the newline mode accordingly.

**`SIO_NL_POSIX`**  
Do not do any translation on input or output.

**`SIO_NL_DOS`**  
Emit a newline (`\n`) as `\r\n`. Discard `\r` from the input.^(245The current implementation does not check that the character is followed by `\`n.)

### 12.9.7 Foreign stream position information

The `IOSTREAM` has a field `position` that points at a structure of type `IOPOS`. This structure is defined as below.

``` code
typedef struct io_position
{ int64_t               byteno;         /* byte-position in file */
  int64_t               charno;         /* character position in file */
  int                   lineno;         /* lineno in file */
  int                   linepos;        /* position in line */
  intptr_t              reserved[2];    /* future extensions */
} IOPOS;
```

If a stream is created using the flag `SIO_RECORDPOS` the IO functions maintain the position information. Note that, supporting the ISO stream position data (see [stream_property/2](IO.html#stream_property/2)), both the byte and character position is maintained. These may differ if the stream uses a multibyte encoding.

The `linepos` is updated as follows: `\n` and `\r` reset the position to 0 (zero). The backspace (`\b`) decrements the position if it is positive. The tab (`\t`) tabs to the next multiple of 8. Any other character increments the line position by one. Future versions may change that, notably the tab width might no longer be hard coded.

### 12.9.8 Support functions for blob save/load

The functions in this sections are intended to support *blobs* to define [save()](foreigninclude.html#save()) and [load()](foreigninclude.html#load()) functions so they can be part of a saved state or `.qlf` file. The matching pair of functions is guaranteed to give the same result, regardless of byte ordering (big or little endian). The user must not make any assumptions on the exact data format used for storing the data. The atom read/write functions can only be used from the blob callback functions.

For saving an uninterpreted array of bytes, it is suggested that the length is output as a `size_t` value using [PL_qlf_put_uint32()](foreign-streams.html#PL_qlf_put_uint32()) followed by the bytes using [Sfwrite()](foreign-streams.html#Sfwrite()); and for loading, the length is read using [PL_qlf_get_uint32()](foreign-streams.html#PL_qlf_get_uint32()), a buffer is allocated, and the bytes are read using [Sfread()](foreign-streams.html#Sfread()).

`int` **PL_qlf_put_int64**(`int64_t i, IOSTREAM *s`)  
`int` **PL_qlf_put_int32**(`int32_t i, IOSTREAM *s`)  
`int` **PL_qlf_put_uint32**(`uint32 i, IOSTREAM *s`)  
Write integers of several sizes. Signed integers are written in *zigzag* encoding. For unsigned integers we only write the non-zero bytes. The result is compact and the same for big or little endian.

`int` **PL_qlf_put_double**(`double f, IOSTREAM *s`)  
Write double as 8 byte big endian.

`int` **PL_qlf_put_atom**(`atom_t a, IOSTREAM *s`)  
Write an atom. The atom may be a *blob*. Note that this function may *only* be used from a blob [save()](foreigninclude.html#save()) function. Calling from another context results in a fatal error.

`int` **PL_qlf_get_int64**(`IOSTREAM *s, int64_t *ip`)  
`int` **PL_qlf_get_int32**(`IOSTREAM *s, int32_t *ip`)  
`int` **PL_qlf_get_uint32**(`IOSTREAM *s, uint32_t *ip`)  
`int` **PL_qlf_get_double**(`IOSTREAM *s, double *fp`)  
Counterparts of corresponding PL_qlf_put\_\*() functions.

`int` **PL_qlf_get_atom**(`IOSTREAM *s, atom_t *ap`)  
Counterpart of [PL_qlf_put_atom()](foreign-streams.html#PL_qlf_put_atom()). Again, this may *only* be used in the context of a blob [load()](foreigninclude.html#load()) function.
