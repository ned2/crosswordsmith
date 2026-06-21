
## 15.4 Structure of a pack

A *pack* is a directory that has two obligatory components:

1.  A directory named `prolog`. When the pack is attached, this directory is added to the `library` file search path. This implies that any `.pl` file that appears in this directory can be loaded into Prolog using `:- use_module(library(file)).` Alternatively, a file from a specific package can be loaded using e.g., `:- use_module(pack(environ/prolog/environ)).`
2.  A file `pack.pl`. This file provides the *meta data* for the pack. See [section 15.5.1](pack-devel.html#sec:15.5.1) for details.

In addition, a pack may, and often does, include *foreign code*. The current system provides support for classical Unix make files, GNU autoconf/automake and CMake. See [section 15.5.2](pack-devel.html#sec:15.5.2) for details. This build infrastructure is also used to test the package.

A pack can be made accessible in two ways

1.  As an archive file. This file must be named as below, where version is a dotted version number and \<`ext`\> is either `.tgz` (gzipped tar archive) or `.zip`.

    ``` code
    <pack>-<version>.<ext>
    ```

    The pack contains the contents of the package. The root of the archive is identified by locating the file `pack.pl`. Extraction ignores the path leading to this file. Typically, the archive contains a single directory named after the package name without version.

    Installing packs from archives requires that SWI-Prolog has the `archive` extension installed. When a package is registered with the central package server the server identifies it by the SHA1 hash of the archive. It is therefore important that the archive is never modified after registration. If *any* modification is required (including comments, documentation, etc,) the user *must* create a new version.

2.  A git repository. This is now the preferred option because it provides a persistent location and easy version management.
