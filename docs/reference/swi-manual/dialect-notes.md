
## C.2 Notes on specific dialects

The level of maturity of the various dialect emulation implementations varies enormously. All of them have been developed to realise portability for one or more, often large, programs. This section provides some notes on emulating a particular dialect.

### C.2.1 Notes on specific dialects

[XSB](http://xsb.sourceforge.net/) Prolog compatibility emerged from a project to integrate XSB's advanced tabling support in SWI-Prolog (see [section 7](tabling.html#sec:7)). This project has been made possible by [Kyndi](https://kyndi.com/).^(257This project was initiated by Benjamin Grosof and carried out in cooperation with Theresa Swift, David S. Warren and Fabrizio Riguzzi.) The XSB dialect implementation has been created to share as much as possible of the XSB test suite as well as some larger programs to evaluate both tabling implementations. The dialect emulation was extended to support [Pharos](https://github.com/cmu-sei/pharos).^(258Pharos was used to evaluate *incremental tabling* ([section 7.7](tabling-incremental.html#sec:7.7)), a protect with Edward Schwatz and Cory Cohen from CMU).

Emulating XSB is relatively complicated due to the large distance from the Quintus descendant Prolog systems. Notably XSB's name based module system is hard to map on SWI-Prolog's predicate based module system. As a result, only non-modular projects or projects with basic usage of modules are supported. For the development of new projects that require modules more advanced module support we suggest using [Logtalk](https://logtalk.org/).

#### C.2.1.1 Loading XSB source files

SWI-Prolog's emulation of XSB depends on the XSB preferred file name extension `.P`. This extension is used by `library(dialect/xsb/source)` to initiate a two phase loading process based on [term_expansion/2](consulting.html#term_expansion/2) of the virtual term `begin_of_file`.

1.  In the first phase the file is read with XSB compatible operator declarations and all directives (:- Term) are extracted. The directives are used to determine that the file defines a module (iff the file contains an [export/1](altmoduleapi.html#export/1) directive) and construct a SWI-Prolog compatible module declaration. As XSB has a two phase compiler where SWI has a single phase compiler, this is also used to move some directives to the start of the file.
2.  The second phase loads the file as normal.

To load a project in both XSB and SWI-Prolog it is advised to make sure all source files use the `.P` file name extension. Next, write a SWI-Prolog loader in a `.pl` file that contains e.g.,

``` code
:- use_module(library(dialect/xsb/source)).

:- [main_file].
```

It is also possible to put the able [use_module/1](import.html#use_module/1) directive in your personal initialization file (see [section 2.2](initfile.html#sec:2.2)), after which XSB files can be loaded as normal SWI-Prolog files using

``` code
% swipl file.P
```

XSB code may depend on the **gpp** preprocessor. We do not provide **gpp**. It is however possible to send XSB source files through **gpp** by loading `library(library/dialect/xsb/gpp)`. This require **gpp** to be accessible through the environment variable `PATH` or the [file_search_path/2](consulting.html#file_search_path/2) alias `path`. We refer to the `gpp` library for details.

### C.2.2 The XSB import directive

The XSB import directive takes the form as below.

``` code
:- import p/1, q/2, ... from <lib>.
```

This import directive is resolved as follows:

- If the referenced library is found as a local file, it is loaded and the requested predicates are imported.
- Otherwise, the referenced library is searched for in the `dialect/xsb` directory of the SWI-Prolog library. If found, the predicates are imported from this library.
- The referenced predicates are searched for in SWI-Prolog built-in predicates and the SWI-Prolog library. If found, they are made available if necessary.
