
## 4.36 File System Interaction

The predicates in this section provide interaction with the file system and (syntactic) operations on file names. SWI-Prolog file system interaction is based on the POSIX standard. On Windows we use the Unicode (wide character) versions of the C runtime library functions.

The file operations define a large set of error conditions. Errors are mapped to Prolog exceptions using a generic function that receives the `action` (e.g., `make_directory`), `type` (e.g., `directory`), the term that describes the object (name) of the file system and the `errno` value. Unfortunately, the resulting exceptions are often misleading. For example, calling [make_directory/1](files.html#make_directory/1) such that it must create multiple directories (e.g., `d1/d2/d3`) returns an existence error on the directory `d1/d2/d3` rather than the missing component. On Windows the situation is even worse because many of its runtime functions distinguish only a few error codes. For example, **\_wmkdir()** only produces `EEXIST` or `ENOENT` and all failures except for an already existing target result in an existence_error exception. We can only improve on this situation by implementing the translation of `errno` to a Prolog exception specifically for each file operation and perform additional tests to distinguish the different error conditions that are represented by the same `errno` value. This is hard, in particular because this translation needs to depend on the OS and specific file system limitations.

SWI-Prolog uses the Windows Unicode functions to access the file system. Internally all Prolog's file handling is based on C `char*` strings. On POSIX systems these strings use *multibyte encodings* according to the current *locale* (nowadays often UTF-8). On Windows the encoding is fixed to UTF-8 and a wrapper around the low-level file functions translates this to UTF-16^(152Before 8.5.16 to UCS-2, allowing only for Unicode code points up to 0xffff.) before calling the Win32 \***W()** or C runtime \_w\*() functions.

Windows absolute file paths are traditionally limited to 260 characters (`PATH_MAX`). More recent versions of Windows support long files. Many of the Unicode file functions support longer paths. For others, support for long paths can be forced by prefixing the path with “`\\?\`” (\`\``\\?UNC\`” for UNC paths (//server/path)). This syntax does *not* allow for relative paths. Thus, SWI-Prolog file functions use **GetFullPathName()** to arrive at absolute canonical paths and apply the appropriate prefix before calling Windows low-level functions.

Unfortunately the support is not very robust. Some functions still apply length limits while others do not work with the above mentioned prefix. The result depends on the Windows version, several Windows registry entries and the file system. Exceeding the length limit is often reported as non-existence of the target file or directory.

**access_file**(`+File, +Mode`)  
True if `File` exists and can be accessed by this Prolog process under mode `Mode`. `Mode` is one of the atoms `read`, `write`, `append`, `execute`, `search`, `exist`, or `none`. Fails silently otherwise. `File` may also be the name of a directory. `access_file(File, none)` simply succeeds without testing anything.

If `Mode` is `write` or `append`, this predicate also succeeds if the file does not exist and the user has write access to the directory of the specified location.

The mode `execute` is only intended for use with regular files and the mode `search` only with directories. However, the two modes are currently equivalent and both can be used with either files or directories. This may change in the future, so the results of checking `execute` access on directories or `search` access on regular files should not be relied on.

The behaviour is backed up by the POSIX **access()** API. The Windows replacement (**\_waccess()**) returns incorrect results because it does not consider ACLs (Access Control Lists). The Prolog flag [win_file_access_check](flags.html#flag:win_file_access_check) may be used to control the level of checking performed by Prolog. Please note that checking access never provides a guarantee that a subsequent open succeeds without errors due to inherent concurrency in file operations. It is generally more robust to try and open the file and handle possible exceptions. See [open/4](IO.html#open/4) and [catch/3](exception.html#catch/3).

**exists_file**(`+File`)  
True if `File` exists and is a *regular* file. This does not imply the user has read or write access to the file. See also [exists_directory/1](files.html#exists_directory/1) and [access_file/2](files.html#access_file/2). The current implementation fails silently, also on error error conditions such as `File` being too long.

**file_directory_name**(`+File, -Directory`)  
Extracts the directory part of `File`. This predicate removes the longest match for the regular expression `/*[^/]*/*$`. If the result is empty it binds `Directory` to `/` if the first character of `File` is `/` and `.` otherwise. The behaviour is consistent with the POSIX **dirname** program.^(153Before SWI-Prolog 7.7.13 trailing `/` where *not* removed, translation `/a/b/` into `/a/b`. Volker Wysk pointed at this incorrect behaviour.)

See also directory_file_path/3 from `library(filesex)`. The system ensures that for every valid `Path` using the Prolog (POSIX) directory separators, following is true on systems with a sound implementation of [same_file/2](files.html#same_file/2).^(154On some systems, `Path` and `Path2` refer to the same entry in the file system, but [same_file/2](files.html#same_file/2) may fail.) See also [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2).

``` code
        ...,
        file_directory_name(FilePath, Dir),
        file_base_name(FilePath, File),
        directory_file_path(Dir, File, Path2),
        same_file(FilePath, Path2).
```

**file_base_name**(`+File, -Name`)  
Extracts the file name part from name that may include directories. Similar to [file_directory_name/2](files.html#file_directory_name/2) the extraction is based on the regex `/*([^/]*)/*$`, now capturing the non-`/` segment. If the segment is empty it unifies `-Name` with `/` if `File` starts with `/` and the empty atom (`''`) otherwise. The behaviour is consistent with the POSIX **basename** program.^(155Before SWI-Prolog 7.7.13, if argPath ended with a `/` `-Name` was unified with the empty atom.)

**same_file**(`+File1, +File2`)  
True if both filenames refer to the same physical file. That is, if `File1` and `File2` are the same string or both names exist and point to the same file (due to hard or symbolic links and/or relative vs. absolute paths). On systems that provide **stat()** with meaningful values for `st_dev` and `st_inode`, [same_file/2](files.html#same_file/2) is implemented by comparing the device and inode identifiers. On Windows, [same_file/2](files.html#same_file/2) uses **GetFileInformationByHandle()** and compares the volume serial number and file index.^(156As of version 8.5.16. Earlier versions only compare the canonical name obtained using **GetFullPathName()**.)

**exists_directory**(`+Directory`)  
True if `Directory` exists and is a directory. This does not imply the user has read, search or write permission for the directory. The current implementation fails silently, also on error error conditions such as `Directory` being too long.

**delete_file**(`+File`)  
Remove `File` from the file system. Note that on POSIX systems the **remove()** call works on read-only files as long as the containing directory has write access. The Windows **remove()** call raises a permission error if the file is read-only. SWI-Prolog removes the read-only attribute if **remove()** fails and tries again. If the file still cannot be removed it restores the read-only attribute and [delete_file/1](files.html#delete_file/1) raises a permission error. As a consequence a read-only file that cannot be removed is briefly read-write. Also note that while an open file can be removed on POSIX systems (where it is actually deleted when closed), deleting an open file on Windows is not possible.

**rename_file**(`+File1, +File2`)  
Rename `File1` as `File2`. The semantics is compatible to the semantics of the POSIX **rename()** system call as far as the operating system allows. Notably, if `File2` exists, the operation succeeds (except for possible permission errors) and is *atomic* (meaning there is no window where `File2` does not exist). Note that `File2` cannot be an existing directory.^(157The POSIX semantics describe one exception: a directory can be moved to an existing *empty* directory.) To move a file to another directory one must create `File2` from the target directory and the base name of `File1`. See [file_base_name/2](files.html#file_base_name/2).

The **rename()** system call has a large number of error conditions. Errors are mapped to Prolog exceptions using a generic conversion based on the `File1` argument. As a result, the errors may be confusing. Future versions may improve on this.

**size_file**(`+File, -Size`)  
Unify `Size` with the size of `File` in bytes.

**time_file**(`+File, -Time`)  
Unify the last modification time of `File` with `Time`. `Time` is a floating point number expressing the seconds elapsed since Jan 1, 1970. See also convert_time/\[2,8\] and [get_time/1](system.html#get_time/1).

**absolute_file_name**(`+File, -Absolute`)  
Expand a local filename into an absolute path. The absolute path is canonicalised: references to `.`, `..` and repeated directory separators (`/`) are deleted. This predicate ensures that expanding a filename returns the same absolute path regardless of how the file is addressed. Notably, if a file appears in multiple directories due to symbolic or hard links [absolute_file_name/2](files.html#absolute_file_name/2) returns the same absolute filename. SWI-Prolog uses absolute filenames to register source files independent of the current working directory. The directory separators are always `/`; [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2) can be used to obtain the the operating system's preferred form.

This predicate has a different history than [absolute_file_name/3](files.html#absolute_file_name/3) and should primarily be used to get an absolute canonical name from a relative name. If `File` is a term Alias(Relative) is behaviour is defined as below, i.e., if an accessible file can be found using the provided search path this is returned. Otherwise it returns the the expansion of the alias path.^(158The SICStus implementation behaves as [absolute_file_name/3](files.html#absolute_file_name/3) with an empty option list.) Users are advised to use [absolute_file_name/3](files.html#absolute_file_name/3) with appropriate options for resolving an Alias(Relative) term.

``` code
absolute_file_name(Spec, AbsFile) :-
    absolute_file_name(Spec, File, [access(read), file_errors(fail)]),
    !,
    AbsFile = File.
absolute_file_name(Spec, AbsFile) :-
    absolute_file_name(Spec, AbsFile, []).
```

See also [absolute_file_name/3](files.html#absolute_file_name/3), [file_search_path/2](consulting.html#file_search_path/2), and [expand_file_name/2](files.html#expand_file_name/2).

**absolute_file_name**(`+Spec, -Absolute, +Options`)  
Convert the given file specification into an absolute path. `Spec` is a term Alias(Relative) (e.g., `(library(lists)`), a relative filename or an absolute filename. The primary intention of this predicate is to resolve files specified as Alias(Relative), which use [file_search_path/2](consulting.html#file_search_path/2) to look up the possibilities for Alias. With the default options no file system check is performed and the alias-resolved path is returned as is; in particular, it may name a directory. To require a non-directory pass `file_type(regular)`, `file_type(prolog)`, or an `access(Mode)` other than `none`. Pass `file_type(directory)` to look up directories. For example, when the library holds both a `chr.pl` file and a `chr` directory, `absolute_file_name(library(chr), F)` resolves to the directory, while `absolute_file_name(library(chr), F, [file_type(prolog)])` resolves to `chr.pl`. The result always uses the directory separator `/`; if the operating system uses something different, SWI-Prolog converts the file name before it makes an OS call. If you need the filename in the OS's preferred form, use [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2). Supported `Options` are:

**extensions**(`ListOfExtensions`)  
List of file extensions to try. Default is `[’’]`. For each extension, [absolute_file_name/3](files.html#absolute_file_name/3) will first add the extension and then verify the conditions imposed by the other options. If the condition fails, the next extension on the list is tried. Extensions may be specified both as `.ext` or plain `ext`.

**relative_to**(`+FileOrDir`)  
Resolve the path relative to the given directory or the directory holding the given file. Without this option, paths are resolved relative to the working directory (see [working_directory/2](files.html#working_directory/2)) or, if `Spec` is atomic and [absolute_file_name/\[2,3\]](files.html#absolute_file_name/2) is executed in a directive, it uses the current source file as reference. Note that a directive using [initialization/1](consulting.html#initialization/1) is executed *after* loading the file. This implies that such paths are resolved relative to the working directory for the toplevel file and relative to the file loading the file holding the [initialization/1](consulting.html#initialization/1) directive otherwise.

Up to version 9.3.9, the system tried *both* the directory holding the current source and the current working directory.

**access**(`Mode`)  
Imposes the condition access_file(`File`, `Mode`). `Mode` is one of `read`, `write`, `append`, `execute`, `search`, `exist` or `none`. See also [access_file/2](files.html#access_file/2). The default is `none` which, if `file_type` is **not** specified as `directory` or `regular`, returns absolute file names that result from expanding aliases without inspecting the actual file system.

**file_type**(`Type`)  
Defines extensions. Current mapping: `txt` implies `[’’]`, `prolog` implies `[’.pl’,’’]`, `executable` implies `[’.so’,’’]` and `qlf` implies `[’.qlf’,’’]`. The `Type` `directory` implies `[’’]` and causes this predicate to generate (only) directories. The `Type` `regular` is the opposite of `directory` and is the default if no file type is specified and the requested access is not `none`. Since the default access is `none`, `file_type(regular)` is *not* implied by default; see the introduction above.

The file type `source` is an alias for `prolog` for compatibility with SICStus Prolog. See also [prolog_file_type/2](consulting.html#prolog_file_type/2).

**file_errors**(`fail/error`)  
If `error` (default), throw an `existence_error` exception if the file cannot be found. If `fail`, stay silent.^(159Silent operation was the default up to version 3.2.6.)

**solutions**(`first/all`)  
If `first` (default), the predicate leaves no choice point. Otherwise a choice point will be left and backtracking may yield more solutions.

**expand**(`Boolean`)  
If `true` (default is `false`) and `Spec` is atomic, call [expand_file_name/2](files.html#expand_file_name/2) followed by [member/2](lists.html#member/2) on `Spec` before proceeding. This is a SWI-Prolog extension intended to minimise porting effort after SWI-Prolog stopped expanding environment variables and the `~` by default. This option should be considered deprecated. In particular the use of *wildcard* patterns such as `*` should be avoided.

The Prolog flag [verbose_file_search](flags.html#flag:verbose_file_search) can be set to `true` to help debugging Prolog's search for files. See also [file_search_path/2](consulting.html#file_search_path/2).

This predicate is derived from Quintus Prolog. In Quintus Prolog, the argument order was `absolute_file_name(+Spec, +Options, -Path)`. The argument order has been changed for compatibility with ISO and SICStus. The Quintus argument order is still accepted.

**is_absolute_file_name**(`+File`)  
True if `File` specifies an absolute path name. On POSIX systems, this implies the path starts with a‘/’. For Microsoft-based systems this implies the path starts with `<``letter``>:` or `//<``host``>/`. This predicate is intended to provide platform-independent checking for absolute paths. See also [absolute_file_name/2](files.html#absolute_file_name/2) and [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2).

**file_name_extension**(`?Base, ?Extension, ?Name`)  
This predicate is used to add, remove or test filename extensions. The main reason for its introduction is to deal with different filename properties in a portable manner. If the file system is case-insensitive, testing for an extension will also be done case-insensitive. `Extension` may be specified with or without a leading dot (`.`). If an `Extension` is generated, it will not have a leading dot.

**directory_files**(`+Directory, -Entries`)  
Unify `Entries` with a list of entries in `Directory`. Each member of `Entries` is an atom denoting an entry relative to `Directory`. `Entries` contains all entries, including hidden files and, if supplied by the OS, the special entries `.` and `..`. See also [expand_file_name/2](files.html#expand_file_name/2).^(160This predicate should be considered a misnomer because it returns entries rather than files. We stick to this name for compatibility with, e.g., SICStus, Ciao and YAP.)

**expand_file_name**(`+WildCard, -List`)  
Unify `List` with a sorted list of files or directories matching `WildCard`. The normal Unix wildcard constructs‘`?`’,‘`*`’,‘`[ ... ]`’and‘`{...}`’are recognised. The interpretation of‘`{...}`’is slightly different from the C shell (csh(1)). The comma-separated argument can be arbitrary patterns, including‘`{...}`’patterns. The empty pattern is legal as well:‘`{.pl,}`’matches either‘`.pl`’or the empty string.

If the pattern contains wildcard characters, only existing files and directories are returned. Expanding a‘pattern’without wildcard characters returns the argument, regardless of whether or not it exists.

Before expanding wildcards, the construct `$\arg{var}` is expanded to the value of the environment variable `var`, and a possible leading `~` character is expanded to the user's home directory.^(161On Windows, the home directory is determined as follows: if the environment variable `HOME` exists, this is used. If the variables `HOMEDRIVE` and `HOMEPATH` exist (Windows-NT), these are used. At initialisation, the system will set the environment variable `HOME` to point to the SWI-Prolog home directory if neither `HOME` nor `HOMEPATH` and `HOMEDRIVE` are defined.)

**prolog_to_os_filename**(`?PrologPath, ?OsPath`)  
Convert between the internal Prolog path name conventions and the operating system path name conventions. The internal conventions follow the POSIX standard, which implies that this predicate is equivalent to =/2 (unify) on POSIX (e.g., Unix) systems. On Windows systems it changes the directory separator from `\` into `/`.

**read_link**(`+File, -Link, -Target`)  
If `File` points to a symbolic link, unify `Link` with the value of the link and `Target` to the file the link is pointing to. `Target` points to a file, directory or non-existing entry in the file system, but never to a link. Fails if `File` is not a link. Fails always on systems that do not support symbolic links.

\[deprecated\]**tmp_file**(`+Base, -TmpName`)  
Create a name for a temporary file. `Base` is an identifier for the category of file. The `TmpName` is guaranteed to be unique. If the system halts, it will automatically remove all created temporary files. `Base` is used as part of the final filename. Portable applications should limit themselves to alphanumeric characters. The directory for temporary files is defined by the Prolog flag [tmp_dir](flags.html#flag:tmp_dir).

Because it is possible to guess the generated filename, attackers may create the filesystem entry as a link and possibly create a security issue. New code should use [tmp_file_stream/3](files.html#tmp_file_stream/3).

**tmp_file_stream**(`+Encoding, -FileName, -Stream`)  
**tmp_file_stream**(`-FileName, -Stream, +Options`)  
Create a temporary filename `FileName`, open it for writing and unify `Stream` with the output stream. If the OS supports it, the created file is only accessible to the current user and the file is created using the **open()**-flag `O_EXCL`, which guarantees that the file did not exist before this call. The directory for temporary files is defined by the Prolog flag [tmp_dir](flags.html#flag:tmp_dir). The following options are processed:

**encoding**(`+Encoding`)  
Encoding of `Stream`. Default is the value of the Prolog flag [encoding](flags.html#flag:encoding). The value `binary` opens the file in binary mode.

**extension**(`+Ext`)  
Ensure the created file has the given extension. Default is no extension. Using an extension may be necessary to run external programs on the file.

This predicate is a safe replacement of [tmp_file/2](files.html#tmp_file/2). Note that in those cases where the temporary file is needed to store output from an external command, the file must be closed first. E.g., the following downloads a file from a URL to a temporary file and opens the file for reading (on Unix systems you can delete the file for cleanup after opening it for reading):

``` code
open_url(URL, In) :-
        tmp_file_stream(text, File, Stream),
        close(Stream),
        process_create(curl, ['-o', File, URL], []),
        open(File, read, In),
        delete_file(File).              % Unix-only
```

Temporary files created using this call are removed if the Prolog process terminates *gracefully*. Calling [delete_file/1](files.html#delete_file/1) using `FileName` removes the file and removes the entry from the administration of files-to-be-deleted.

**make_directory**(`+Directory`)  
Create a new directory (folder) on the filesystem. Raises an exception on failure. On Unix systems, the directory is created with default permissions (defined by the process *umask* setting).

**delete_directory**(`+Directory`)  
Delete directory (folder) from the filesystem. Raises an exception on failure. Please note that in general it will not be possible to delete a non-empty directory.

**working_directory**(`-Old, +New`)  
Unify `Old` with an absolute path to the current working directory and change working directory to `New`. Use the pattern `working_directory(CWD, CWD)` to get the current directory. See also [absolute_file_name/2](files.html#absolute_file_name/2) and [chdir/1](files.html#chdir/1).^(bugSome of the file I/O predicates use local filenames. Changing directory while file-bound streams are open causes wrong results on [telling/1](IO.html#telling/1), [seeing/1](IO.html#seeing/1) and [current_stream/3](IO.html#current_stream/3).) Note that the working directory is shared between all threads. Applications are strongly encouraged not to change the working directory or change the working directory once during the initialization.

**chdir**(`+Path`)  
Compatibility predicate. New code should use [working_directory/2](files.html#working_directory/2).
