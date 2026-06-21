
## 2.14 Automatic loading of libraries

If ---at runtime--- an undefined predicate is trapped, the system will first try to import the predicate from the module's default module (see [section 6.10](importmodule.html#sec:6.10). If this fails the *auto loader* is activated.^(27Actually, the hook user:exception/3 is called; only if this hook fails it calls the autoloader.) On first activation an index to all library files in all library directories is loaded in core (see [library_directory/1](consulting.html#library_directory/1), [file_search_path/2](consulting.html#file_search_path/2) and [reload_library_index/0](autoload.html#reload_library_index/0)). If the undefined predicate can be located in one of the libraries, that library file is automatically loaded and the call to the (previously undefined) predicate is restarted. By default this mechanism loads the file silently. The [current_prolog_flag/2](flags.html#current_prolog_flag/2) key [verbose_autoload](flags.html#flag:verbose_autoload) is provided to get verbose loading. The Prolog flag [autoload](flags.html#flag:autoload) can be used to enable/disable the autoload system. A more controlled form of autoloading as well as lazy loading application modules is provided by [autoload/\[1,2\]](module-autoload.html#autoload/1).

Autoloading only handles (library) source files that use the module mechanism described in [chapter 6](modules.html#sec:6). The files are loaded with [use_module/2](import.html#use_module/2) and only the trapped undefined predicate is imported into the module where the undefined predicate was called. Each library directory must hold a file `INDEX.pl` that contains an index to all library files in the directory. This file consists of lines of the following format:

``` code
index(Name, Arity, Module, File).
```

The predicate [make/0](consulting.html#make/0) updates the autoload index. It searches for all library directories (see [library_directory/1](consulting.html#library_directory/1) and [file_search_path/2](consulting.html#file_search_path/2)) holding the file `MKINDEX.pl` or `INDEX.pl`. If the current user can write or create the file `INDEX.pl` and it does not exist or is older than the directory or one of its files, the index for this directory is updated. If the file `MKINDEX.pl` exists, updating is achieved by loading this file, normally containing a directive calling [make_library_index/2](autoload.html#make_library_index/2). Otherwise [make_library_index/1](autoload.html#make_library_index/1) is called, creating an index for all `*.pl` files containing a module.

Below is an example creating an indexed library directory.

``` code
% mkdir ~/${XDG_DATA_HOME-.config}/swi-prolog/lib
% cd ~/${XDG_DATA_HOME-.config}/swi-prolog/lib
% swipl -g 'make_library_index(.)' -t halt
```

If there is more than one library file containing the desired predicate, the following search schema is followed:

1.  If there is a library file that defines the module in which the undefined predicate is trapped, this file is used.
2.  Otherwise library files are considered in the order they appear in the [library_directory/1](consulting.html#library_directory/1) predicate and within the directory alphabetically.

**autoload_path**(`+DirAlias`)  
Add `DirAlias` to the libraries that are used by the autoloader. This extends the search path `autoload` and reloads the library index. For example:

``` code
:- autoload_path(library(http)).
```

If this call appears as a directive, it is term-expanded into a clause for user:file_search_path/2 and a directive calling [reload_library_index/0](autoload.html#reload_library_index/0). This keeps source information and allows for removing this directive.

**make_library_index**(`+Directory`)  
Create an index for this directory. The index is written to the file’INDEX.pl’in the specified directory. Fails with a warning if the directory does not exist or is write protected.

**make_library_index**(`+Directory, +ListOfPatterns`)  
Normally used in `MKINDEX.pl`, this predicate creates `INDEX.pl` for `Directory`, indexing all files that match one of the file patterns in `ListOfPatterns`.

Sometimes library packages consist of one public load file and a number of files used by this load file, exporting predicates that should not be used directly by the end user. Such a library can be placed in a sub-directory of the library and the files containing public functionality can be added to the index of the library. As an example we give the XPCE library's `MKINDEX.pl`, including the public functionality of `trace/browse.pl` to the autoloadable predicates for the XPCE package.

``` code
:- prolog_load_context(directory, Dir),
   make_library_index(Dir,
                      [ '*.pl',
                        'trace/browse.pl',
                        'swi/*.pl'
                      ]).
```

**reload_library_index**  
Force reloading the index after modifying the set of library directories by changing the rules for [library_directory/1](consulting.html#library_directory/1), [file_search_path/2](consulting.html#file_search_path/2), adding or deleting `INDEX.pl` files. This predicate does *not* update the `INDEX.pl` files. Check [make_library_index/\[1,2\]](autoload.html#make_library_index/1) and [make/0](consulting.html#make/0) for updating the index files.

Normally, the index is reloaded automatically if a predicate cannot be found in the index and the set of library directories has changed. Using [reload_library_index/0](autoload.html#reload_library_index/0) is necessary if directories are removed or the order of the library directories is changed.

When creating an executable using either [qsave_program/2](saved-states.html#qsave_program/2) or the **-c** command line options, it is necessary to load all predicates that would normally be autoloaded explicitly. This is discussed in [section 14](runtime.html#sec:14). See [autoload_all/0](saved-states.html#autoload_all/0).
