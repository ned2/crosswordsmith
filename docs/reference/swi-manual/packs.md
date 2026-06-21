
# 15 Packs: community add-ons

SWI-Prolog has a mechanism for incorporating community extensions called *packs*. See the [pack landing page](https://www.swi-prolog.org/pack/list) for details and available packs. This chapter discusses how packages can be attached to the current Prolog process, how they can be installed as well as developing packages.

Packs are installed as self-containing directories that provide additional Prolog libraries and *foreign modules*, compiled native code plugins. In addition, a pack can define *apps*, command line tools that can be started using `swipl app [args]` (see [section 2.11.1.5](compilation.html#sec:2.11.1.5)). Packs are searched as sub-directories of the Prolog search path `pack`. Initially, this search path is the user's *App data*, followed by the system's *App data*. The searched directories can be found using

``` code
?- absolute_file_name(pack(.), Path, [solutions(all)]).
```

The search path can be managed using the environment variable `SWIPL_PACK_PATH`, the **-p** command line option or using [attach_packs/2](pack-attach.html#attach_packs/2).

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[15.1 Installing packs](pack-install.html)

[15.2 Built-in predicates for attaching packs](pack-attach.html)

[15.3 library(prolog_pack): A package manager for Prolog](prologpack.html)

[15.4 Structure of a pack](pack-structure.html)

[15.5 Developing a pack](pack-devel.html)

[15.5.1 The pack meta data](pack-devel.html#sec:15.5.1)

[15.5.1.1 Pack requirements on Prolog](pack-devel.html#sec:15.5.1.1)

[15.5.2 Packs with foreign code](pack-devel.html#sec:15.5.2)

[15.5.2.1 Compiling a foreign extension using a simple Makefile](pack-devel.html#sec:15.5.2.1)

[15.5.2.2 Publishing a pack](pack-devel.html#sec:15.5.2.2)

[15.5.2.3 Compiling a foreign extension using CMake](pack-devel.html#sec:15.5.2.3)

[15.5.3 Updating a package](pack-devel.html#sec:15.5.3)
