
## 15.1 Installing packs

As of version 9.1.22, SWI-Prolog supports three models for managing packs: *shared packages* are added to the user or system environment, while *project specific packages* are added to a particular project only. Finally, project specific packages can be managed as *git submodules*. These three approaches are discussed in more detail below.

Using [pack_install/2](prologpack.html#pack_install/2) we can install a package either for the current user or globally for all users.

**Shared packages**

System-wide installations is satisfactory as long as all projects can use the same version of a pack, the packs required by all projects do not conflict, and redistribution of the projects is not a primary concern. For example, if you frequently require RocksDB for several projects you are working on, installing the `rocksdb` pack as user is appropriate.

The shared model is similar to e.g., Python's **pip** installer. Python resolves dealing with packages for a specific project using *virtual environments*, where each virtual environment provides a selection of packages. A Python virtual environment may be *activated* for the current shell, which modifies the shell's environment variables.

**Project specific packages**

Alternatively, SWI-Prolog allows packs to be installed as part of a project. This approach is also found with **npm**, the Node.js package manager. Using project-specific packs with SWI-Prolog requires calling [attach_packs/2](pack-attach.html#attach_packs/2) before loading any library from a pack. To use (only) packs from the local sub directory `packs`, add this directive to the code that uses it:

``` code
:- attach_packs(packs, [replace(true)]).
```

Packs can be installed into the `packs` directory directly using [pack_install/2](prologpack.html#pack_install/2) with the `pack_directory(Dir)` option or using the `pack` *app* as

``` code
swipl pack install --dir=packs <pack>
```

The preferred way is to use [pack_install_local/3](prologpack.html#pack_install_local/3). This predicate takes a *closure* to collect the desired packages, creates an installation plan and executes this. This ensures a set of compatible packs at their latest available version or explicitly specified versions. Typically, one would create a file `packs.pl` according to the example below to install the packages required by a project. By using such a file it is easy to replicate a suitable set of installed packs for anyone who wishes to use your application.

``` code
:- module(packs, []).
:- use_module(library(prolog_pack)).
:- attach_packs(packs, [replace(true)]).

:- initialization(install, main).

pack(scasp, [commit('HEAD')]).
pack(environ, []).
pack(date_time, []).

install :-
    pack_install_local(pack, packs, []).
```

Here, the [attach_packs/2](pack-attach.html#attach_packs/2) must be the same as used by the project. The first argument of pack_install_local/2 refers to pack/2 , generating a list of target packages and options for each package. The options for each pack are defined by [pack_install/2](prologpack.html#pack_install/2). They typically refer to the download location and required version. Given the above, we can install the packages we need for a project using

``` code
swipl packs.pl
```

**Using GIT submodules**

Alternative to the above, if the desired packs are all available as git repository, we can add packs to our git managed projects by adding the packs as git submodules to our project. For example, we add a pack to the `packs` directory as

``` code
mkdir packs
git submodule add https://github.com/SWI-Prolog/sCASP.git packs/scasp
git submodule add https://github.com/fnogatz/tap.git tap
```

As above, we can must make our project use the local packs by calling [pack_attach/2](pack-attach.html#pack_attach/2). After fetching all submodules we can build the foreign components and/or run the tests of the attached packs using the steps below

``` code
?- attach_packs(packs, [replace(true)]).
?- pack_rebuild.
```

Using git submodules gives full control of the pack versions you are using. It also makes you responsible of adding dependencies and taking care of version dependencies between packs. Finally, it limits you to using git based packages.
