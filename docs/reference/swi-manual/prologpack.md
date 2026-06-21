
## 15.3 library(prolog_pack): A package manager for Prolog

The `library(prolog_pack)` provides the SWI-Prolog package manager. This library lets you inspect installed packages, install packages, remove packages, etc. This library complemented by the built-in predicates such as [attach_packs/2](pack-attach.html#attach_packs/2) that makes installed packages available as libraries.

The important functionality of this library is encapsulated in the *app* `pack`. For help, run

``` code
swipl pack help
```

\[det\]**pack_list_installed**  
List currently installed packages and report possible dependency issues.

**pack_info**(`+Pack`)  
Print more detailed information about `Pack`.

\[det\]**pack_list**(`+Query`)  
\[det\]**pack_list**(`+Query, +Options`)  
\[det\]**pack_search**(`+Query`)  
`Query` package server and installed packages and display results. `Query` is matches case-insensitively against the name and title of known and installed packages. For each matching package, a single line is displayed that provides:

- Installation status
  - **p**: package, not installed
  - **i**: installed package; up-to-date with public version
  - **a**: as **i**, but installed only as dependency
  - **U**: installed package; can be upgraded
  - **A**: installed package; newer than publically available
  - **l**: installed package; not on server
- Name@Version
- Name@Version(ServerVersion)
- Title

`Options` processed:

**installed**(`true`)  
Only list packages that are locally installed. Contacts the server to compare our local version to the latest available version.

**outdated**(`true`)  
Only list packages that need to be updated. This option implies `installed(true)`.

**server**(`Server``|``false`)  
If `false`, do not contact the server. This implies `installed(true)`. Otherwise, use the given pack server.

Hint: `?- pack_list('').` lists all known packages.

The predicates [pack_list/1](prologpack.html#pack_list/1) and [pack_search/1](prologpack.html#pack_search/1) are synonyms. Both contact the package server at [https://www.swi-prolog.org](https://www.swi-prolog.org) to find available packages. Contacting the server can be avoided using the `server(false)` option.

\[det\]**pack_install**(`+Spec:atom`)  
\[det\]**pack_install**(`+SpecOrList, +Options`)  
Install one or more packs from `SpecOrList`. `SpecOrList` is a single specification or a list of specifications. A specification is one of

- A pack name. This queries the pack repository at [https://www.swi-prolog.org](https://www.swi-prolog.org)
- Archive file name
- A `http(s)` URL of an archive file name. This URL may contain a star (\*) for the version. In this case [pack_install/1](prologpack.html#pack_install/1) asks for the directory content and selects the latest version.
- An https GIT URL
- A local directory name given as `file://` URL
- `'.'`, in which case a relative symlink is created to the current directory (all other options for `Spec` make a copy of the files). Installation using a symlink is normally used during development of a pack.

Processes the options below. Default options as would be used by [pack_install/1](prologpack.html#pack_install/1) are used to complete the provided `Options`. Note that [pack_install/2](prologpack.html#pack_install/2) can be used through the SWI-Prolog command line app `pack` as below. Most of the options of this predicate are available as command line options.

swipl pack install `<`name`>`

`Options`:

**url**(`+URL`)  
Source for downloading the package

**pack_directory**(`+Dir`)  
Directory into which to install the package.

**global**(`+Boolean`)  
If `true`, install in the XDG common application data path, making the pack accessible to everyone. If `false`, install in the XDG user application data path, making the pack accessible for the current user only. If the option is absent, use the first existing and writable directory. If that doesn't exist find locations where it can be created and prompt the user to do so.

**insecure**(`+Boolean`)  
When `true` (default `false`), do not perform any checks on SSL certificates when downloading using `https`.

**interactive**(`+Boolean`)  
Use default answer without asking the user if there is a default action.

**silent**(`+Boolean`)  
If `true` (default false), suppress informational progress messages.

**upgrade**(`+Boolean`)  
If `true` (default `false`), upgrade package if it is already installed.

**rebuild**(`Condition`)  
Rebuild the foreign components. `Condition` is one of `if_absent` (default, do nothing if the directory with foreign resources exists), `make` (run `make`) or `true` (run‘make distclean\` followed by the default configure and build steps).

**test**(`Boolean`)  
If `true` (default), run the pack tests.

**git**(`+Boolean`)  
If `true` (default `false` unless `URL` ends with `.git`), assume the URL is a GIT repository.

**link**(`+Boolean`)  
Can be used if the installation source is a local directory and the file system supports symbolic links. In this case the system adds the current directory to the pack registration using a symbolic link and performs the local installation steps.

**version**(`+Version`)  
Demand the pack to satisfy some version requirement. `Version` is as defined by [require_version/3](prologversions.html#require_version/3). For example `'1.5'` is the same as `>=('1.5')`.

**branch**(`+Branch`)  
When installing from a git repository, clone this branch.

**commit**(`+Commit`)  
When installing from a git repository, checkout this commit. `Commit` is either a hash, a tag, a branch or `'HEAD'`.

**build_type**(`+Type`)  
When building using CMake, use `-DCMAKE_BUILD_TYPE=Type`. Default is the build type of Prolog or `Release`.

**register**(`+Boolean`)  
If `true` (default), register packages as downloaded after performing the download. This contacts the server with the meta-data of each pack that was downloaded. The server will either register the location as a new version or increment the download count. The server stores the IP address of the client. Subsequent downloads of the same version from the same IP address are ignored.

**server**(`+URL`)  
Pack server to contact. Default is the setting `prolog_pack:server`, by default set to `https://www.swi-prolog.org/pack/`

Non-interactive installation can be established using the option `interactive(false)`. It is adviced to install from a particular *trusted* URL instead of the plain pack name for unattented operation.

\[det\]**pack_install_local**(`:Spec, +Dir, +Options`)  
Install a number of packages in a local directory. This predicate supports installing packages local to an application rather than globally.

\[det\]**pack_url_file**(`+URL, -File`)  
True if `File` is a unique id for the referenced pack and version. Normally, that is simply the base name, but GitHub archives destroy this picture. Needed by the pack manager in the web server.

\[det\]**pack_rebuild**  
\[det\]**pack_rebuild**(`+Pack`)  
Rebuild possible foreign components of `Pack`. The predicate [pack_rebuild/0](prologpack.html#pack_rebuild/0) rebuilds all registered packs.

\[semidet\]**pack_upgrade**(`+Pack`)  
Upgrade `Pack`. Shorthand for `pack_install(Pack, [upgrade(true)])`.

\[det\]**pack_remove**(`+Name`)  
\[det\]**pack_remove**(`+Name, +Options`)  
Remove the indicated package. If packages depend (indirectly) on this pack, ask to remove these as well. `Options`:

**interactive**(`false`)  
Do not prompt the user.

**dependencies**(`Boolean`)  
If `true` delete dependencies without asking.

\[det\]**pack_publish**(`+Spec, +Options`)  
Publish a package. There are two ways typical ways to call this. We recommend developing a pack in a GIT repository. In this scenario the pack can be published using

``` code
?- pack_publish('.', []).
```

Alternatively, an archive file has been uploaded to a public location. In this scenario we can publish the pack using

``` code
?- pack_publish(URL, [])
```

In both scenarios, [pack_publish/2](prologpack.html#pack_publish/2) by default creates an isolated environment and installs the package in this directory from the public URL. On success it triggers the pack server to register the URL as a new pack or a new release of a pack.

Packs may also be published using the *app* `pack`, e.g.

``` code
swipl pack publish .
```

`Options`:

**git**(`Boolean`)  
If `true`, and `Spec` is a git managed directory, install using the remote repo.

**sign**(`Boolean`)  
Sign the repository with the current version. This runs `git tag -s <tag>`.

**force**(`Boolean`)  
Force the git tag. This runs `git tag -f <tag>`.

**branch**(`+Branch`)  
`Branch` used for releases. Defined by git_default_branch/2 if not specified.

**register**(`+Boolean`)  
If `false` (default `true`), perform the installation, but do not upload to the server. This can be used for testing.

**isolated**(`+Boolean`)  
If `true` (default), install and build all packages in an isolated package directory. If `false`, use other packages installed for the environment. The latter may be used to speedup debugging.

**pack_directory**(`+Dir`)  
Install the temporary packages in `Dir`. If omitted [pack_publish/2](prologpack.html#pack_publish/2) creates a temporary directory and deletes this directory after completion. An explict target `Dir` is created if it does not exist and is not deleted on completion.

**clean**(`+Boolean`)  
If `true` (default), clean the destination directory first

\[nondet\]**pack_property**(`?Pack, ?Property`)  
True when `Property` is a property of an installed `Pack`. This interface is intended for programs that wish to interact with the package manager. Defined properties are:

**directory**(`Directory`)  
`Directory` into which the package is installed

**version**(`Version`)  
Installed version

**title**(`Title`)  
Full title of the package

**author**(`Author`)  
Registered author

**download**(`URL`)  
Official download `URL`

**readme**(`File`)  
Package `README` file (if present)

**todo**(`File`)  
Package `TODO` file (if present)
