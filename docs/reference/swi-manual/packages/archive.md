SWI-Prolog binding to libarchive

Jan Wielemaker  
VU University Amsterdam  
The Netherlands  
E-mail: [J.Wielemaker@vu.nl](mailto:J.Wielemaker@vu.nl)

Abstract

The library [libarchive](https://github.com/libarchive/libarchive) provides a portable way to access archive files as well as encoded (typically compressed) data. This package is a Prolog wrapper around this library. The motivation to introduce this library is twofold. In the first place, it provides a minimal platform independent API to access archives. In the second place, it allows accessing archives through Prolog streams, which often eliminates the need for temporary files and all related consequences for performance, security and platform dependency.

# Table of Contents

[1 Motivation](#sec:1)

[2 library(archive): Access several archive formats](#sec:2)

[3 Status](#sec:3)

## 1 Motivation

Archives play two roles: they combine multiple documents into a single one and they typically provide compression and sometimes encryption or other services. Bundling multiple resources into a single archive may greatly simplify distribution and guarantee that the individual resources are consistent. SWI-Prolog provides archiving using its (rather arcane) saved-state format. See resource/3 and open_resource/3. It also provides compression by means of library(zlib).

External archives may be accessed through the process interface provided by process_create/3, but this has disadvantages. The one that motivated this library was that using external processes provide no decent platform independent access to archives. Most likely zip files come closest to platform independent access, but there are many different programs for accessing zip files that provide slightly different sets of options and the existence of any of these programs cannot be guaranteed without distributing our own bundled version. Similar arguments hold for Unix tar archives, where just about any Unix-derives system has a tar program but except for very basic commands, the command line options are not compatible and tar is not part of Windows. The only format granted on Windows is .cab, but a program to create them is not part of Windows and the .cab format is rare outside the Windows context.

Discarding availability of archive programs, each archive program comes with its own set of command line options and its own features and limitations. Fortunately, [libarchive](https://github.com/libarchive/libarchive) provides a consistent interface to a wealth of compression and archiving formats. The library `library(archive)` wraps this library, providing access to archives using Prolog streams both for the archive as a whole and the archive entries. E.g., archives may be read from Prolog streams and each member in turn may be processed using Prolog streams without materialising data using temporary files.

## 2 library(archive): Access several archive formats

See also  
[https://github.com/libarchive/libarchive/](https://github.com/libarchive/libarchive/)

This library uses *libarchive* to access a variety of archive formats. The following example lists the entries in an archive:

``` code
list_archive(File) :-
    setup_call_cleanup(
        archive_open(File, Archive, []),
        (   repeat,
            (   archive_next_header(Archive, Path)
            ->  format('~w~n', [Path]),
                fail
            ;   !
            )
        ),
        archive_close(Archive)).
```

Here is an alternative way of doing this, using [archive_foldl/4](#archive_foldl/4), a higher level predicate.

``` code
list_archive2(File) :-
    list_archive(File, Headers),
    maplist(writeln, Headers).

list_archive2(File, Headers) :-
    archive_foldl(add_header, File, Headers, []).

add_header(Path, _, [Path|Paths], Paths).
```

Here is another example which counts the files in the archive and prints file type information, also using [archive_foldl/4](#archive_foldl/4):

``` code
print_entry(Path, Handle, Cnt0, Cnt1) :-
    archive_header_property(Handle, filetype(Type)),
    format('File ~w is of type ~w~n', [Path, Type]),
    Cnt1 is Cnt0 + 1.

list_archive_headers(File) :-
    archive_foldl(print_entry, File, 0, FileCount),
    format('We have ~w files', [FileCount]).
```

\[det\]**archive_open**(`+Data, -Archive, +Options`)  
Wrapper around [archive_open/4](#archive_open/4) that opens the archive in read mode.

\[det\]**archive_open**(`+Data, +Mode, -Archive, +Options`)  
Open the archive in `Data` and unify `Archive` with a handle to the opened archive. `Data` is either a file name (as accepted by open/4) or a stream that has been opened with the option `type(binary)`. If `Data` is an already open stream, the caller is responsible for closing it (but see option `close_parent(true)`) and must not close the stream until after [archive_close/1](#archive_close/1) is called. `Mode` is either `read` or `write`. Details are controlled by `Options`. Typically, the option `close_parent(true)` is used to also close the `Data` stream if the archive is closed using [archive_close/1](#archive_close/1). For other options when reading, the defaults are typically fine - for writing, a valid format and optional filters must be specified. The option `format(raw)` must be used to process compressed streams that do not contain explicit entries (e.g., gzip'ed data) unambibuously. The `raw` format creates a *pseudo archive* holding a single member named `data`.

**close_parent**(`+Boolean`)  
If this option is `true` (default `false`), `Data` stream is closed when [archive_close/1](#archive_close/1) is called on `Archive`. If `Data` is a file name, the default is `true`.

**compression**(`+Compression`)  
Synomym for `filter(Compression)`. Deprecated.

**filter**(`+Filter`)  
Support the indicated filter. This option may be used multiple times to support multiple filters. In read mode, If no filter options are provided, `all` is assumed. In write mode, `none` is assumed. Supported values are `all`, `bzip2`, `compress`, `gzip`, `grzip`, `lrzip`, `lzip`, `lzma`, `lzop`, `none`, `rpm`, `uu` and `xz`. The value `all` is default for read, `none` for write.

**format**(`+Format`)  
Support the indicated format. This option may be used multiple times to support multiple formats in read mode. In write mode, you must supply a single format. If no format options are provided, `all` is assumed for read mode. Note that `all` does **not** include `raw` and `mtree`. To open both archive and non-archive files, *both* `format(all)` and `format(raw)` and/or `format(mtree)` must be specified. Supported values are: `all`, `7zip`, `ar`, `cab`, `cpio`, `empty`, `gnutar`, `iso9660`, `lha`, `mtree`, `rar`, `raw`, `tar`, `xar` and `zip`. The value `all` is default for read.

Note that the actually supported compression types and formats may vary depending on the version and installation options of the underlying libarchive library. This predicate raises a domain or permission error if the (explicitly) requested format or filter is not supported.

Errors  
\- `domain_error(filter, Filter)` if the requested filter is invalid (e.g., `all` for writing).  
- `domain_error(format, Format)` if the requested format type is not supported.  
- `permission_error(set, filter, Filter)` if the requested filter is not supported.

\[det\]**archive_close**(`+Archive`)  
Close the archive. If `close_parent(true)` was specified in [archive_open/4](#archive_open/4), the underlying entry stream is closed too. If there is an entry opened with [archive_open_entry/2](#archive_open_entry/2), actually closing the archive is delayed until the stream associated with the entry is closed. This can be used to open a stream to an archive entry without having to worry about closing the archive:

``` code
archive_open_named(ArchiveFile, EntryName, Stream) :-
    archive_open(ArchiveFile, Archive, []),
    archive_next_header(Archive, EntryName),
    archive_open_entry(Archive, Stream),
    archive_close(Archive).
```

\[nondet\]**archive_property**(`+Handle, ?Property`)  
True when `Property` is a property of the archive `Handle`. Defined properties are:

**filters**(`List`)  
True when the indicated filters are applied before reaching the archive format.

\[semidet\]**archive_next_header**(`+Handle, -Name`)  
Forward to the next entry of the archive for which `Name` unifies with the pathname of the entry. Fails silently if the end of the archive is reached before success. `Name` is typically specified if a single entry must be accessed and unbound otherwise. The following example opens a Prolog stream to a given archive entry. Note that *Stream* must be closed using close/1 and the archive must be closed using [archive_close/1](#archive_close/1) after the data has been used. See also setup_call_cleanup/3.

``` code
open_archive_entry(ArchiveFile, EntryName, Stream) :-
    open(ArchiveFile, read, In, [type(binary)]),
    archive_open(In, Archive, [close_parent(true)]),
    archive_next_header(Archive, EntryName),
    archive_open_entry(Archive, Stream).
```

Errors  
`permission_error(next_header, archive, Handle)` if a previously opened entry is not closed.

\[det\]**archive_open_entry**(`+Archive, -Stream`)  
Open the current entry as a stream. `Stream` must be closed. If the stream is not closed before the next call to [archive_next_header/2](#archive_next_header/2), a permission error is raised.

**archive_set_header_property**(`+Archive, +Property`)  
Set `Property` of the current header. Write-mode only. Defined properties are:

**filetype**(`-Type`)  
`Type` is one of `file`, `link`, `socket`, `character_device`, `block_device`, `directory` or `fifo`. It appears that this library can also return other values. These are returned as an integer.

**mtime**(`-Time`)  
True when entry was last modified at time.

**size**(`-Bytes`)  
True when entry is `Bytes` long.

**link_target**(`-Target`)  
`Target` for a link. Currently only supported for symbolic links.

**archive_header_property**(`+Archive, ?Property`)  
True when `Property` is a property of the current header. Defined properties are:

**filetype**(`-Type`)  
`Type` is one of `file`, `link`, `socket`, `character_device`, `block_device`, `directory` or `fifo`. It appears that this library can also return other values. These are returned as an integer.

**mtime**(`-Time`)  
True when entry was last modified at time.

**size**(`-Bytes`)  
True when entry is `Bytes` long.

**link_target**(`-Target`)  
`Target` for a link. Currently only supported for symbolic links.

**format**(`-Format`)  
Provides the name of the archive format applicable to the current entry. The returned value is the lowercase version of the output of `archive_format_name()`.

**permissions**(`-Integer`)  
True when entry has the indicated permission mask.

**archive_extract**(`+ArchiveFile, +Dir, +Options`)  
Extract files from the given archive into `Dir`. Supported options:

**remove_prefix**(`+Prefix`)  
Strip `Prefix` from all entries before extracting. If `Prefix` is a list, then each prefix is tried in order, succeding at the first one that matches. If no prefixes match, an error is reported. If `Prefix` is an atom, then that prefix is removed.

**exclude**(`+ListOfPatterns`)  
Ignore members that match one of the given patterns. Patterns are handed to wildcard_match/2.

**include**(`+ListOfPatterns`)  
Include members that match one of the given patterns. Patterns are handed to wildcard_match/2. The `exclude` options takes preference if a member matches both the `include` and the `exclude` option.

Errors  
\- `existence_error(directory, Dir)` if `Dir` does not exist or is not a directory.  
- `domain_error(path_prefix(Prefix), Path)` if a path in the archive does not start with Prefix

To be done  
Add options

\[det\]**archive_entries**(`+Archive, -Paths`)  
True when `Paths` is a list of pathnames appearing in `Archive`.

\[nondet\]**archive_data_stream**(`+Archive, -DataStream, +Options`)  
True when `DataStream` is a stream to a data object inside `Archive`. This predicate transparently unpacks data inside *possibly nested* archives, e.g., a *tar* file inside a *zip* file. It applies the appropriate decompression filters and thus ensures that Prolog reads the plain data from `DataStream`. `DataStream` must be closed after the content has been processed. Backtracking opens the next member of the (nested) archive. This predicate processes the following options:

**meta_data**(`-Data:list(dict)`)  
If provided, `Data` is unified with a list of filters applied to the (nested) archive to open the current `DataStream`. The first element describes the outermost archive. Each `Data` dict contains the header properties ([archive_header_property/2](#archive_header_property/2)) as well as the keys:

**filters**(`Filters:list(atom)`)  
Filter list as obtained from [archive_property/2](#archive_property/2)

**name**(`Atom`)  
Name of the entry.

Non-archive files are handled as pseudo-archives that hold a single stream. This is implemented by using [archive_open/3](#archive_open/3) with the options `[format(all),format(raw)]`.

\[det\]**archive_create**(`+OutputFile, +InputFiles, +Options`)  
Convenience predicate to create an archive in `OutputFile` with data from a list of `InputFiles` and the given `Options`.

Besides options supported by [archive_open/4](#archive_open/4), the following options are supported:

**directory**(`+Directory`)  
Changes the directory before adding input files. If this is specified, paths of input files must be relative to `Directory` and archived files will not have `Directory` as leading path. This is to simulate `-C` option of the `tar` program.

**format**(`+Format`)  
Write mode supports the following formats:‘7zip\`, `cpio`, `gnutar`, `iso9660`, `xar` and `zip`. Note that a particular installation may support only a subset of these, depending on the configuration of `libarchive`.

**archive_foldl**(`:Goal, +Archive, +State0, -State`)  
Operates like foldl/4 but for the entries in the archive. For each member of the archive, `Goal` called as‘call(:`Goal`, +Path, +Handle, +S0, -S1). Here, `S0` is current state of the *accumulator* (starting with `State0`) and `S1` is the next state of the accumulator, producing `State` after the last member of the archive.

|           |                                                          |
|-----------|----------------------------------------------------------|
| `Archive` | File name or stream to be given to archive_open/\[3,4\]. |

See also  
[archive_header_property/2](#archive_header_property/2), [archive_open/4](#archive_open/4).

## 3 Status

The current version is merely a proof-of-concept. It lacks writing archives and does not support many of the options of the underlying library. The main motivation for starting this library was to achieve portability of the upcomming SWI-Prolog package distribution system. Other functionality will be added on‘as needed’basis.

# Index

?  
[archive_close/1](#archive_close/1)  
[archive_create/3](#archive_create/3)  
[archive_data_stream/3](#archive_data_stream/3)  
[archive_entries/2](#archive_entries/2)  
[archive_extract/3](#archive_extract/3)  
[archive_foldl/4](#archive_foldl/4)  
[archive_header_property/2](#archive_header_property/2)  
[archive_next_header/2](#archive_next_header/2)  
[archive_open/3](#archive_open/3)  
[archive_open/4](#archive_open/4)  
[archive_open_entry/2](#archive_open_entry/2)  
[archive_property/2](#archive_property/2)  
[archive_set_header_property/2](#archive_set_header_property/2)  
open_resource/3  
[1](#idx:openresource3:2)

process_create/3  
[1](#idx:processcreate3:3)

resource/3  
[1](#idx:resource3:1)
