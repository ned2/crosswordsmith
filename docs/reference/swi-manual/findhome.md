
## 12.6 The Prolog‘home’directory

Executables embedding SWI-Prolog should be able to find the‘home’directory of the development environment unless a self-contained saved state has been added to the executable (see [qsave_program/\[1,2\]](saved-states.html#qsave_program/1) and [section 12.5](plld.html#sec:12.5)).

If Prolog starts up, it will try to locate the development environment. To do so, it will try the following steps until one succeeds:

1.  If the **--home=DIR** is provided, use this.
2.  If the environment variable `SWI_HOME_DIR` is defined and points to an existing directory, use this.
3.  If the environment variable `SWIPL` is defined and points to an existing directory, use this.
4.  Locate the primary executable or (Windows only) a component (*module*) thereof and check whether the parent directory of the directory holding this file contains the file `swipl`. If so, this file contains the (relative) path to the home directory. If this directory exists, use this. This is the normal mechanism used by the binary distribution.
5.  If the precompiled path exists, use it. This is only useful for a source installation.

If all fails and there is no state attached to the executable or provided Windows module (see [PL_initialise()](foreigninclude.html#PL_initialise())), SWI-Prolog gives up. If a state is attached, the current working directory is used.

The [file_search_path/2](consulting.html#file_search_path/2) alias `swi` is set to point to the home directory located.
