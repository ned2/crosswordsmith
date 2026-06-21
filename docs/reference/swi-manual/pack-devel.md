
## 15.5 Developing a pack

We recommend using GIT for developing SWI-Prolog packages. To start a new package, invent a name and verify that the name is not yet in use at The [pack landing page](https://www.swi-prolog.org/pack/list). Create a directory with this name, a sub-directory `prolog` and a the metadata file `pack.pl` that contains at least the name and version of the pack. Below is a simple example. See [section 15.5.1](pack-devel.html#sec:15.5.1) for all possible metadata fields.

``` code
name(hello).
version('1.0.0').
title('Hello world').
keywords([demo]).
author( 'Bob Programmer, 'bob123@programmer.me. ).
download('https://github/bob123/hello.git').
```

Now, add the Prolog libraries provided to the `prolog` directory. While doing so, please pay attention to the points below. If your are looking for examples of well structured libraries, please look at the system libraries.

- Only add *module* files to the `prolog` directory.
- Be aware that the modules your pack provides are globally accessible as `library(File)`. Thus, make sure the name is fairly unique and the module name is typically the same as the *base name* of the file.
- Modules that need not be immediately visible to the user should be placed in a subdirectory. Typically one uses the pack name to name the subdirectory.^(251Using e.g., `private` is not a good idea because the `private` directory of each pack using this would be available as `library(Pack/private)`.) Use e.g., pack\_\<`name`\> for the module names of the private files.
- Consider documenting the files using PlDoc.

Once the pack is ready for a very first test, we can make it accessible using the command below. On non-Windows systems, this makes the pack accessible using a *symbolic link* from your personal pack directory to this directory.

``` code
swipl pack install .
```

After this command the new libraries should be available when you start a new SWI-Prolog process. Another way to make the pack accessible is by using the `pack` search path (see [file_search_path/2](consulting.html#file_search_path/2)). The command (from the pack directory) is

``` code
swipl -p pack=..
```

### 15.5.1 The pack meta data

A pack must have a file `pack.pl` in its root directory. The file contains Prolog terms. Defined terms are below. The argument types are types registered with [must_be/2](error.html#must_be/2) and described in the running text.

**name**(`atom`)  
Name of the pack. This should be the same as the directory name. Names can be constructed from the ASCII letters, underscore and digits, e.g., `[a-zA-Z9-0_]+`

**title**(`atom`)  
Short summary of the package. Do not use line breaks and limit respect at maximum length of about 40 characters.

**keywords**(`list(atom)`)  
List of keywords that help finding your pack. There is no fixed set of keywords to choose from.

**description**(`list(atom)`)  
Longer description as a list of lines.

**version**(`version`)  
Current version of the pack. This is a list of integers separated by dots. There is no limit to the number of sub revisions.

**author**(`atom, email_or_url_or_empty`)  
Original author of the code. If the contact address is unknown it may be omitted (empty atom). Repeat this term for multiple authors.

**maintainer**(`atom, email_or_url`)  
**packager**(`atom, email_or_url`)  
As `author`, but the contact cannot be empty. May be repeated.

**pack_version**(`nonneg`)  
Package convention number. Currently 1 (default) or 2. Version 2 provides better support for building foreign extensions.

**home**(`atom`)  
Location of th home page. This is typically a URL.

**download**(`atom`)  
Location for downloading. This is either the URL of the GIT repository or a wildcard URL for downloading the archive, e.g., [https://me.com/packs/mypack-\*.zip](https://me.com/packs/mypack-*.zip). An upgrade request fetches the [https://me.com/packs/](https://me.com/packs/), expecting an HTML page with links to the available versions. It then selects the latest version.

**provides**(`atom`)  
Announce that the pack provides facilities identified by the given token. Optionally, the token may be given a version using `@(Token,Version)`. A pack implicitly provides `@(PackName,PackVersion)`. The supplied tokens operate in the same *name space* as packages and thus the same care must be taken to select a name. Multiple of these claims may be present.

**requires**(`dependency`)  
The pack depends on the availability of `Dependency`. The `Dependency` is a token, normally the name of another package. See `provides`. The dependency may be further refined by writing `Token Cmp Version`, where `Cmp` is one of Prolog's standard numerical comparison operators. See [cmp_versions/3](prologversions.html#cmp_versions/3). This metadata is also used to state requirements on Prolog. See [section 15.5.1.1](pack-devel.html#sec:15.5.1.1). Multiple requirements are expressed with multiple claims.

**conflicts**(`dependency`)  
The pack cannot be use together with the indicated `Dependency`. This is the negation of `requires`.

**replaces**(`atom`)  
This pack replaces some other pack.

**autoload**(`boolean`)  
If `true`, add the library for the package as *autoload* library. This implies that the exported predicates may be used without explicitly importing the library. Use with care.

#### 15.5.1.1 Pack requirements on Prolog

The file `pack.pl` may contain `requires(Requirement)` statements. Normally, `Requirement` is a pack or token, optionally with a version requirement. The requirement `prolog` is reserved for requirements on the Prolog version while `prolog:``Feature` may be used to demand specific features. Feature matching is described with [require_prolog_version/2](prologversions.html#require_prolog_version/2). Multiple requirements on Prolog must all be true. Below are some examples

``` code
requires(prolog >= '9.2').        % 9.2.0 or later
requires(prolog:threads).         % flag threads = true
requires(prolog:library(socket)). % library(socket) exists
requires(prolog:bounded(false)).  % flag bounded = false
```

### 15.5.2 Packs with foreign code

Many packs include C or C++ resources. Such packs include the C or C++ resources in a subdirectory of the pack. There are no restrictions for naming this subdirectory or structuring the source files in this directory. The build process must create native *modules* in the directory `lib/<``arch``>`, where \<`arch`\> is the architecture as obtained by the Prolog flag [arch](flags.html#flag:arch).

The build process identifies control files that tell the package manager which build tool to use. The package manager populates the process environment with variables that provide details about the running Prolog instance. This environment is saved in a file `buildenv.sh` in the pack root or build directory. By *sourcing* this file, the user may run the build tools by hand for debugging purposes.

The build process consists of five steps that are described below

**dependencies**  
This step currently only supports `conan`. It is executed if either `conanfile.txt` or `conanfile.py` is found in the root directory of the pack.

**configure**  
This preparation step is executed if one of `CMakeLists.txt` (**cmake**), `configure`, `configure.in` (**autoconf**), `configure.ac` or `Makefile.am` (**automake**) are found. The program to manage them is in parenthesis.

**build**  
Build the process. When configured using (**cmake**) this will use (**cmake**). Otherwise either `Makefile` or `makefile` is expected and Unix **make** is used to build the process.

**test**  
Test the project. Either uses **cmake** or the GNU convention `make check`.

**install**  
Install the project. Either uses **cmake** or the GNU convention `make install`.

While running the above tools, the environment is populated. The names of the variables provided depends on the `pack_version(Version)` metadata. We give the names for version 2, with the names for version 1 in parenthesis if this differs from the version 2 name.

**`PATH`**  
Contains the environment path with the directory holding the currently running SWI-Prolog instance prepended in front of it. As a result, **swipl** is always present and runs the same SWI-Prolog instance as the current Prolog process.

**`SWIPL`**  
Contains the absolute file name of the running executable.

**`SWIPL_PACK_VERSION`**  
Version of the pack system (1 or 2). If not present we must assume‘1’.

**`SWIPL_VERSION` (`SWIPLVERSION`)**  
Contains the numeric SWI-Prolog version defined as `Major × 10000 + Minor × 100 + Patch`

**`SWIPL_HOME_DIR` (`SWIHOME`)**  
Contains the directory holding the SWI-Prolog home.

**`SWIPL_ARCH` (`SWIARCH`)**  
contains the machine architecture identifier.

**`SWIPL_MODULE_DIR` (`PACKSODIR`)**  
contains the destination directory for shared objects/DLLs relative to a Prolog pack, i.e., `lib/$SWIARCH`.

**`SWIPL_MODULE_LIB` (`SWISOLIB`)**  
The SWI-Prolog library or an empty string when it is not required to link modules against this library (e.g., ELF systems)

**`SWIPL_LIB` (`SWILIB`)**  
The SWI-Prolog library we need to link to for programs that *embed* SWI-Prolog (normally `-lswipl`)

**`SWIPL_INCLUDE_DIRS`**  
CMake style variable that contains the directory holding `SWI-Prolog.h`, `SWI-Stream.h` and `SWI-cpp2.h`.

**`SWIPL_LIBRARIES_DIR`**  
CMake style variable that contains the directory holding `libswipl`

**`SWIPL_CC` (`CC`)**  
C compiler used to build SWI-Prolog.

**`SWIPL_CXX` (`CXX`)**  
C++ compiler used to build SWI-Prolog.

**`SWIPL_LD` (`LD`)**  
Linker used to link SWI-Prolog.

**`SWIPL_CFLAGS` (`CFLAGS`)**  
C-Flags for building extensions. Always contains `-ISWIPL-INCLUDE-DIR`.

**`SWIPL_MODULE_LDFLAGS` (`LDSOFLAGS`)**  
Link flags for linking modules.

**`SWIPL_MODULE_EXT` (`SOEXT`)**  
File name extension for modules (e.g., `.so` or `.dll`)

**`SWIPL_PREFIX` (`PREFIX`)**  
Install prefix for global binaries, libraries and include files.

#### 15.5.2.1 Compiling a foreign extension using a simple Makefile

If the package requires some C code to be compiled that has no dependencies and needs no configuration it is probably easiest to use a simple Unix make file. We assume `pack_version(2)`. Here is a simple `Makefile`. We assume the pack contains a file `c/environ.c` that contains the C source. Following the GNU guidelines, the `Makefile` must define the following targets:

**all (default)**  
Build the foreign extension. In this very simple case we build the resulting module directly in the target directory.

**check**  
Test the package. This is executed after the default build target.

**install**  
Install the package. In this case this does nothing.

**clean**  
Clean the package. This target disposes intermediate build products.

**distclean**  
Restore the package to its fully clean state. This implies that all built products and intermediate build products are removed. The `distclean` target is used by [pack_rebuild/1](prologpack.html#pack_rebuild/1).

``` code
MODULE= $(SWIPL_MODULE_DIR)/environ.$(SOEXT)
CFLAGS= $(SWIPL_CFLAGS)

all:    $(MODULE)

OBJ=c/environ.o

$(MODULE): $(OBJ)
        mkdir -p $(SWIPL_MODULE_DIR)
        $(SWIPL_LD) $(SWIPL_MODULE_LDFLAGS) -o $@ $(OBJ) $(SWIPL_MODULE_LIB)

check::
        $(SWIPL) -g run_tests -t halt test/test_environ.pl
install::
clean:
        rm -f $(OBJ)
distclean: clean
        rm -f $(MODULE)
```

#### 15.5.2.2 Publishing a pack

As described in [section 15.4](pack-structure.html#sec:15.4), a pack is distributed either as an archive file or as a GIT repository. We strongly encourage using a GIT repository as that gives good version and provenance support. Packs may be published by hand by making the archive or git repository available from a globally accessible place on the internet and installing the pack from this location. This process is streamlined, notably for GIT packs using [pack_publish/2](prologpack.html#pack_publish/2) and the *app* `pack`. To publish a pack a local GIT repository that has publicly accessible *origin*,

1.  Update `version(Version)` in `pack.pl`

2.  Commit all changes, make sure the the repository is clean.

3.  Run

    ``` code
    swipl pack publish .
    ```

This will

1.  Verify the repository is clean and on the default branch.
2.  *Tag* the repository with V\<`version`\>. By default, the tag will be *signed*. Please setup signing for GIT or use the “--no-sign\`\` option.
3.  Push the repository and release tag.
4.  Figure out the download location, either from the `download(URL)` metadata or the GIT remote information.
5.  Install the package and its dependencies in a temporary isolated pack environment.
6.  On success, register the pack with the server.
7.  Delete the isolated pack environment.

Similarly, a pack can be published from a public archive using the command below. When using an archive, **never** change the content of the archive but, instead, create a new archive with a new version.

``` code
swipl pack publish URL
```

#### 15.5.2.3 Compiling a foreign extension using CMake

If the package is more complicated, a simple Makefile typically does not suffice. In this case we have two options. One is to use the GNU **autoconf** or **automake**. However, **cmake** is getting more popular and provides much better support for non-POSIX platforms, e.g., Windows. This section discusses building the same package as [section 15.5.2.1](pack-devel.html#sec:15.5.2.1) using **cmake**.

To use **cmake**, add the content below as the file `CMakeLists.txt` to the root directory of the pack. SWI-Prolog ships with a **cmake** *include* file named `swipl.cmake` that deals with most of the configuration issues. Comments in the file below explain the various steps of the process.

``` code
cmake_minimum_required(VERSION 3.10)
project(swipl-pack-environ)

# Include swipl.cmake from the running SWI-Prolog's home
list(INSERT CMAKE_MODULE_PATH 0 $ENV{SWIPL_HOME_DIR}/cmake)
include(swipl)

# Create the library as a CMake module
add_library(environ MODULE c/environ.c)

# Link the library to SWI-Prolog.  This also removes the `lib` prefix
# from the target on systems that define a common library file prefix
target_link_swipl(environ)

# Install the foreign target. `${swipl_module_dir}` contains the
# directory for installing modules for this architecture.

install(TARGETS environ
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/${swipl_module_dir})

# Run  tests.  This  is  executed   before    the   pack  is  installed.
# swipl_test(name) runs Prolog with the command line below.
#
#    swipl -p foreign=${CMAKE_CURRENT_SOURCE_DIR}/${swipl_module_dir} \
#          -p library=${CMAKE_CURRENT_SOURCE_DIR}/prolog \
#          --on-error=status \
#          -g test_${name} \
#          -t halt \
#          ${CMAKE_CURRENT_SOURCE_DIR}/test/test_${name}.pl
#
# This  implies  that  a  test  `name`  must    be  defined  in  a  file
# `test/test_${name}.pl`, which exports a  predicate `test_${name}`. The
# test succeeds if this predicate  succeeds   and  no error messages are
# printed.

enable_testing()
swipl_add_test(environ)
```

### 15.5.3 Updating a package

If a package needs a revision to fix bugs or add functionality it needs to be updated. First, we create a development environment using

1.  Clone the git repository that provides the pack.

2.  Install the pack *as a link* using the command below. If the pack contains foreign build scripts, this creates a file `buildenv.sh` that contains the environment variables for building the pack.

    ``` code
    ?- pack_install(.).
    ```

Next, we can edit the pack sources and rebuild it the chosen build tools after running `source buildenv.sh` to set the appropriate environment variables. After validating that the pack works as expected follow the instructions in [section 15.5.2.2](pack-devel.html#sec:15.5.2.2) to publish the new version.
