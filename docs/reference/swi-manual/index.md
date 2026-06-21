
SWI-Prolog 10.0.2 Reference Manual

SWI-Prolog development team

Abstract

SWI-Prolog is a comprehensive and portable implementation of the Prolog programming language. SWI-Prolog aims to be a robust and scalable implementation supporting a wide range of applications. In particular, it ships with a wide range of interface libraries, providing interfaces to other languages, databases, graphics and networking. It provides extensive support for managing HTML/SGML/XML, JSON, YAML and RDF documents. The system is particularly suited for server applications due to robust support for multithreading and HTTP server libraries.

SWI-Prolog extends Prolog with *tabling* (SGL resolution). Tabling provides better termination properties and avoids repetitive recomputation. Following XSB, SWI-Prolog's tabling supports sound negation using the *Well Founded Semantics*. *Incremental tabling* supports usage as a *Deductive database*.

SWI-Prolog is designed in the‘Edinburgh tradition’. In addition to the ISO Prolog standard it is largely compatible to Quintus, SICStus and YAP Prolog. SWI-Prolog provides a compatibility framework developed in cooperation with YAP and instantiated for YAP, SICStus, IF/Prolog and XSB.

SWI-Prolog aims at providing a rich development environment, including extensive editor support, graphical source-level debugger, autoloading, a‘make’facility to reload edited files and much more. GNU-Emacs, SWI-Prolog editor for Windows, the PDT plugin for Eclipse or a Visual Studio Code plugin provide alternative environments. [SWISH](https://swish.swi-prolog.org) provides a web based environment.

This document gives an overview of the features, system limits and built-in predicates.

### About this document

This manual is written and maintained using LaTeX . The LaTeX source is included in the source distribution of SWI-Prolog. The manual is converted into HTML using a converter distributed with the SWI-Prolog sources. From the same source we generate the PDF version. Sources, binaries and documentation can be downloaded from the [SWI-Prolog download page](https://www.swi-prolog.org/Download.html).

The SWI-Prolog project **home page** is [https://www.swi-prolog.org](https://www.swi-prolog.org)

![](by-sa.png)

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License. To view a copy of this license, visit [http://creativecommons.org/licenses/by-sa/3.0/](http://creativecommons.org/licenses/by-sa/3.0/) or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
