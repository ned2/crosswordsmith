
## 15.2 Built-in predicates for attaching packs

This section documents the built-in predicates to attach packs. Predicates for creating, registering and installing packs are provided by the library `library(prolog_pack)`.

**attach_packs**  
Attaches all packs in subdirectories of directories that are accessible through the *file search path* (see [absolute_file_name/3](files.html#absolute_file_name/3)) `pack`. The default for this search path is given below. See [file_search_path/2](consulting.html#file_search_path/2) for the `app_data` search path.

``` code
user:file_search_path(pack, app_data(pack)).
```

The default path may be overruled with the environment variable `SWIPL_PACK_PATH`. This variable must contain a list of directories separated by the OS-specific [path_sep](flags.html#flag:path_sep).

The predicate [attach_packs/0](pack-attach.html#attach_packs/0) is called on startup of SWI-Prolog.

**attach_packs**(`+Directory`)  
**attach_packs**(`+Directory, +Options`)  
Attach all packs that are subdirectories of `Directory`. `Directory` is translated into a physical directory using [absolute_file_name/3](files.html#absolute_file_name/3). This implies it can be a term `Alias(SubDir)` and the search is relative to the current source file if `Directory` is not an absolute path and these predicates are used as a directive. Defined options are:

**search**(`+Where`)  
Determines the order in which pack library directories are searched. Default is to add new packages at the end (`last`). Using `first`, new packages are added at the start.

**duplicate**(`+Action`)  
Determines what happens if a pack with the same name is already attached. Default is `warning`, which prints a warning and ignores the new pack. Other options are `keep`, which is like `warning` but operates silently and `replace`, which detaches the old pack and attaches the new.

**replace**(`+Boolean`)  
If `true`, unregister all packs before registering the new packs.

The predicate [attach_packs/2](pack-attach.html#attach_packs/2) can be used to attach packages that are bundled with an application. With the option `replace(true)`, [attach_packs/2](pack-attach.html#attach_packs/2) ensures that the application only relies on bundled packs.

**pack_attach**(`+PackDir, +Options`)  
Attach a single package in `PackDir`. `PackDir` is expected to contain a file‘pack.pl\` with the pack metadata and a‘prolog\` directory. Options processed:

**duplicate**(`+Action`)  
What to do if the same package is already installed in a different directory. `Action` is one of

**warning**  
Warn and ignore the package.

**keep**  
Silently ignore the package.

**replace**  
Unregister the existing and insert the new package

**search**(`+Where`)  
Determines the order of searching package library directories. Default is `last`, alternative is `first`.
