
## 4.4 Editor Interface

SWI-Prolog offers an extensible interface which allows the user to edit objects of the program: predicates, modules, files, etc. The editor interface is implemented by [edit/1](edit.html#edit/1) and consists of three parts: *locating*, *selecting* and *starting* the editor. Any of these parts may be customized. See [section 4.4.1](edit.html#sec:4.4.1).

The built-in edit specifications for [edit/1](edit.html#edit/1) (see [prolog_edit:locate/3](edit.html#prolog_edit:locate/3)) are described in the table below:

**Fully specified objects**

\<`Module`\>:\<`Name`\>/\<`Arity`\>

Refers to a predicate

module(\<`Module`\>)

Refers to a module

file(\<`Path`\>)

Refers to a file

source_file(\<`Path`\>)

Refers to a loaded source file

**Ambiguous specifications**

\<`Name`\>/\<`Arity`\>

Refers to this predicate in any module

\<`Name`\>

Refers to (1) the named predicate in any module with any arity, (2) a (source) file, or (3) a module.

**edit**(`+Specification`)  
First, exploit [prolog_edit:locate/3](edit.html#prolog_edit:locate/3) to translate `Specification` into a list of *Locations*. If there is more than one‘hit’, the user is asked to select from the locations found. Finally, [prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1) is used to invoke the user's preferred editor. Typically, [edit/1](edit.html#edit/1) can be handed the name of a predicate, module, basename of a file, XPCE class, XPCE method, etc.

**edit**  
Edit the‘default’file using [edit/1](edit.html#edit/1). The default file is either the first `.pl` file from the commandline (the *associated* file, see the Prolog flag [associated_file](flags.html#flag:associated_file) or the first script file specified using the **-s** or **-l** command line option. When using the Windows shell while SWI-Prolog is associated with the `.pl` extension this is the file loaded by double-clicking a `.pl` file. See also [section 2.11.1.1](compilation.html#sec:2.11.1.1).

### 4.4.1 Customizing the editor interface

The predicates described in this section are *hooks* that can be defined to disambiguate specifications given to [edit/1](edit.html#edit/1), find the related source, and open an editor at the given source location.

**prolog_edit:locate**(`+Spec, -FullSpec, -Location`)  
Where `Spec` is the specification provided through [edit/1](edit.html#edit/1). This multifile predicate is used to enumerate locations where an object satisfying the given `Spec` can be found. `FullSpec` is unified with the complete specification for the object. This distinction is used to allow for ambiguous specifications. For example, if `Spec` is an atom, which appears as the basename of a loaded file and as the name of a predicate, `FullSpec` will be bound to `file(Path)` or `Name`/`Arity`.

`Location` is a list of attributes of the location. Normally, this list will contain the term `file(File)` and, if available, the term `line(Line)`.

**prolog_edit:locate**(`+Spec, -Location`)  
Same as [prolog_edit:locate/3](edit.html#prolog_edit:locate/3), but only deals with fully specified objects.

**prolog_edit:edit_source**(`+Location`)  
Start editor on `Location`. See [prolog_edit:locate/3](edit.html#prolog_edit:locate/3) for the format of a location term. This multifile predicate is normally not defined. If it succeeds, [edit/1](edit.html#edit/1) assumes the editor is started.

If it fails, [edit/1](edit.html#edit/1) uses its internal defaults, which are defined by the Prolog flag [editor](flags.html#flag:editor) and/or the environment variable `EDITOR`. The following rules apply. If the Prolog flag [editor](flags.html#flag:editor) is of the format `$`\<`name`\>, the editor is determined by the environment variable \<`name`\>. Else, if this flag is `pce_emacs` or `built_in` *and* XPCE is loaded or can be loaded, the built-in Emacs clone is used. Else, if the environment `EDITOR` is set, this editor is used. Finally, **vi** is used as default on Unix systems and **notepad** on Windows.

See the default user preferences file `customize/init.pl` for examples.

**prolog_edit:edit_command**(`+Editor, -Command`)  
Determines how `Editor` is to be invoked using [shell/1](system.html#shell/1). `Editor` is the determined editor (see [prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1)), without the full path specification, and without a possible (`.exe`) extension. `Command` is an atom describing the command. The following %-sequences are replaced in `Command` before the result is handed to [shell/1](system.html#shell/1):

|     |                                                 |
|-----|-------------------------------------------------|
| %e  | Replaced by the (OS) command name of the editor |
| %f  | Replaced by the (OS) full path name of the file |
| %d  | Replaced by the line number                     |

If the editor can deal with starting at a specified line, two clauses should be provided. The first pattern invokes the editor with a line number, while the second is used if the line number is unknown.

The default contains definitions for **vi**, **emacs**, **emacsclient**, **vim**, **notepad**`^*` and **wordpad**`^*`. Starred editors do not provide starting at a given line number.

Please contribute your specifications to [bugs@swi-prolog.org](mailto:bugs@swi-prolog.org).

**prolog_edit:load**  
Normally an undefined multifile predicate. This predicate may be defined to provide loading hooks for user extensions to the edit module. For example, XPCE provides the code below to load `library(swi_edit)`, containing definitions to locate classes and methods as well as to bind this package to the PceEmacs built-in editor.

``` code
:- multifile prolog_edit:load/0.

prolog_edit:load :-
        ensure_loaded(library(swi_edit)).
```
