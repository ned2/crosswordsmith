
## A.45 library(prolog_xref): Prolog cross-referencer data collection

See also  
Where this library analyses *source text*, `library(prolog_codewalk)` may be used to analyse *loaded code*. The `library(check)` exploits `library(prolog_codewalk)` to report on e.g., undefined predicates.

bug  
[meta_predicate/1](metapred.html#meta_predicate/1) declarations take the module into consideration. Predicates that are both available as meta-predicate and normal (in different modules) are handled as meta-predicate in all places.

This library collects information on defined and used objects in Prolog source files. Typically these are predicates, but we expect the library to deal with other types of objects in the future. The library is a building block for tools doing dependency tracking in applications. Dependency tracking is useful to reveal the structure of an unknown program or detect missing components at compile time, but also for program transformation or minimising a program saved state by only saving the reachable objects.

The library is exploited by two graphical tools in the SWI-Prolog environment: the XPCE front-end started by [gxref/0](xref.html#gxref/0), and `library(prolog_colour)`, which exploits this library for its syntax highlighting.

For all predicates described below, `Source` is the source that is processed. This is normally a filename in any notation acceptable to the file loading predicates (see [load_files/2](consulting.html#load_files/2)). Input handling is done by the `library(prolog_source)`, which may be hooked to process any source that can be translated into a Prolog stream holding Prolog source text. `Callable` is a callable term (see [callable/1](typetest.html#callable/1)). Callables do not carry a module qualifier unless the referred predicate is not in the module defined by `Source`.

\[semidet,multifile\]prolog:**called_by**(`+Goal, +Module, +Context, -Called`)  
True when `Called` is a list of callable terms called from `Goal`, handled by the predicate `Module`:`Goal` and executed in the context of the module `Context`. Elements of `Called` may be qualified. If not, they are called in the context of the module `Context`.

\[multifile\]prolog:**called_by**(`+Goal, -ListOfCalled`)  
If this succeeds, the cross-referencer assumes `Goal` may call any of the goals in `ListOfCalled`. If this call fails, default meta-goal analysis is used to determine additional called goals.

deprecated  
New code should use [prolog:called_by/4](prologxref.html#prolog:called_by/4)

\[multifile\]prolog:**meta_goal**(`+Goal, -Pattern`)  
Define meta-predicates. See the examples in this file for details.

\[multifile\]prolog:**hook**(`Goal`)  
True if `Goal` is a hook that is called spontaneously (e.g., from foreign code).

\[det\]**xref_source**(`+Source`)  
\[det\]**xref_source**(`+Source, +Options`)  
Generate the cross-reference data for `Source` if not already done and the source is not modified. Checking for modifications is only done for files. `Options` processed:

**silent**(`+Boolean`)  
If `true` (default `false`), emit warning messages.

**module**(`+Module`)  
Define the initial context module to work in.

**register_called**(`+Which`)  
Determines which calls are registerd. `Which` is one of `all`, `non_iso` or `non_built_in` (default).

**comments**(`+CommentHandling`)  
How to handle comments. If `store`, comments are stored into the database as if the file was compiled. If `collect`, comments are entered to the xref database and made available through xref_mode/2 and [xref_comment/4](prologxref.html#xref_comment/4). If `ignore`, comments are simply ignored. Default is to `collect` comments.

**process_include**(`+Boolean`)  
Process the content of included files (default is `true`).

**stream**(`+Stream`)  
Process the input from `Stream` rather than opening `Source`.

|          |                                   |
|----------|-----------------------------------|
| `Source` | File specification or XPCE buffer |

\[det\]**xref_clean**(`+Source`)  
Reset the database for the given source.

**xref_current_source**(`?Source`)  
Check what sources have been analysed.

\[det\]**xref_done**(`+Source, -Time`)  
Cross-reference executed at `Time`

\[nondet\]**xref_called**(`?Source, ?Called, ?By`)  
\[nondet\]**xref_called**(`?Source, ?Called, ?By, ?Cond`)  
\[nondet\]**xref_called**(`?Source, ?Called, ?By, ?Cond, ?Line`)  
True when `By` is called from `Called` in `Source`. Note that [xref_called/3](prologxref.html#xref_called/3) and [xref_called/4](prologxref.html#xref_called/4) use [distinct/2](solutionsequences.html#distinct/2) to return only distinct `Called-By` pairs. The [xref_called/5](prologxref.html#xref_called/5) version may return duplicate `Called-By` if `Called` is called from multiple clauses in `By`, but at most one call per clause.

|  |  |
|----|----|
| `By` | is a head term or one of the reserved terms `'<directive>'(Line)` or `'<public>'(Line)`, indicating the call is from an (often [initialization/1](consulting.html#initialization/1)) directive or there is a [public/1](dynamic.html#public/1) directive that claims the predicate is called from in some untractable way. |
| `Cond` | is the (accumulated) condition as defined by `:- if(Cond)` under which the calling code is compiled. |
| `Line` | is the *start line* of the calling clause. |

\[nondet\]**xref_defined**(`?Source, +Goal, ?How`)  
Test if `Goal` is accessible in `Source`. If this is the case, `How` specifies the reason why the predicate is accessible. Note that this predicate does not deal with built-in or global predicates, just locally defined and imported ones. `How` is one of of the terms below. Location is one of Line (an integer) or File:Line if the definition comes from an included (using `:-` `include(File)`) directive.

- `dynamic(Location)`
- `thread_local(Location)`
- `multifile(Location)`
- `public(Location)`
- `local(Location)`
- `foreign(Location)`
- `constraint(Location)`
- `imported(From)`
- dcg

**xref_definition_line**(`+How, -Line`)  
If the 3th argument of xref_defined contains line info, return this in `Line`.

\[nondet\]**xref_exported**(`?Source, ?Head`)  
True when `Source` exports `Head`.

\[nondet\]**xref_module**(`?Source, ?Module`)  
True if `Module` is defined in `Source`.

\[nondet\]**xref_uses_file**(`?Source, ?Spec, ?Path`)  
True when `Source` tries to load a file using `Spec`.

|  |  |
|----|----|
| `Spec` | is a specification for [absolute_file_name/3](files.html#absolute_file_name/3) |
| `Path` | is either an absolute file name of the target file or the atom `<not_found>`. |

\[nondet\]**xref_op**(`?Source, Op`)  
Give the operators active inside the module. This is intended to setup the environment for incremental parsing of a term from the source-file.

|      |                                             |
|------|---------------------------------------------|
| `Op` | Term of the form `op(Priority, Type, Name)` |

\[nondet\]**xref_prolog_flag**(`?Source, ?Flag, ?Value, ?Line`)  
True when `Flag` is set to `Value` at `Line` in `Source`. This is intended to support incremental parsing of a term from the source-file.

\[nondet\]**xref_comment**(`?Source, ?Title, ?Comment`)  
Is true when `Source` has a section comment with `Title` and `Comment`

\[nondet\]**xref_comment**(`?Source, ?Head, ?Summary, ?Comment`)  
Is true when `Head` in `Source` has the given PlDoc comment.

\[nondet\]**xref_mode**(`?Source, ?Mode, ?Det`)  
Is true when `Source` provides a predicate with `Mode` and determinism.

\[nondet\]**xref_option**(`?Source, ?Option`)  
True when `Source` was processed using `Option`. Options are defined with [xref_source/2](prologxref.html#xref_source/2).

\[semidet\]**xref_meta**(`+Source, +Head, -Called`)  
True when `Head` calls `Called` in `Source`.

|  |  |
|----|----|
| `Called` | is a list of called terms, terms of the form Term+Extra or terms of the form `//`(Term). |

\[semidet\]**xref_meta**(`+Head, -Called`)  
\[semidet\]**xref_meta_src**(`+Head, -Called, +Src`)  
True when `Called` is a list of terms called from `Head`. Each element in `Called` can be of the form Term+Int, which means that Term must be extended with Int additional arguments. The variant [xref_meta/3](prologxref.html#xref_meta/3) first queries the local context.

deprecated  
New code should use [xref_meta/3](prologxref.html#xref_meta/3).

To be done  
\- Split predifined in several categories. E.g., the ISO predicates cannot be redefined.  
- Rely on the meta_predicate property for many predicates.

**xref_hook**(`?Callable`)  
Definition of known hooks. Hooks that can be called in any module are unqualified. Other hooks are qualified with the module where they are called.

\[semidet\]**xref_public_list**(`+Spec, +Source, +Options`)  
Find meta-information about File. If `Spec` resolves to a Prolog source file, this predicate reads all terms upto the first term that is not a directive. If `Spec` resolves to a SWI-Prolog‘.qlf\` file, it extracts part of the information from the QLF file. It uses the module and meta_predicate directives to assemble the information in `Options`. `Options` processed:

**path**(`-Path`)  
`Path` is the full path name of the referenced file. If `Spec` resolves to a .qlf file, `Path` is the name of the embedded Prolog file.

**module**(`-Module`)  
`Module` is the module defines in `Spec`.

**exports**(`-Exports`)  
`Exports` is a list of predicate indicators and operators collected from the [module/2](defmodule.html#module/2) term and reexport declarations.

**public** `Public``-`  
`Public` declarations of the file. Currently always `[]` for .qlf files.

**meta**(`-Meta`)  
`Meta` is a list of heads as they appear in [meta_predicate/1](metapred.html#meta_predicate/1) declarations. Currently always `[]` for .qlf files.

**silent**(`+Boolean`)  
Do not print any messages or raise exceptions on errors.

The information collected by this predicate is cached. The cached data is considered valid as long as the modification time of the file does not change.

|          |                                              |
|----------|----------------------------------------------|
| `Source` | is the file from which `Spec` is referenced. |

\[semidet\]**xref_public_list**(`+File, -Path, -Export, +Src`)  
\[semidet\]**xref_public_list**(`+File, -Path, -Module, -Export, -Meta, +Src`)  
\[semidet\]**xref_public_list**(`+File, -Path, -Module, -Export, -Public, -Meta, +Src`)  
Find meta-information about `File`. This predicate reads all terms upto the first term that is not a directive. It uses the module and meta_predicate directives to assemble the information described below.

These predicates fail if `File` is not a module-file.

|  |  |
|----|----|
| `Path` | is the canonical path to `File` |
| `Module` | is the module defined in `Path` |
| `Export` | is a list of predicate indicators. |
| `Meta` | is a list of heads as they appear in [meta_predicate/1](metapred.html#meta_predicate/1) declarations. |
| `Src` | is the place from which `File` is referenced. |

deprecated  
New code should use [xref_public_list/3](prologxref.html#xref_public_list/3), which unifies all variations using an option list.

\[semidet\]**xref_source_file**(`+Spec, -File, +Src`)  
\[semidet\]**xref_source_file**(`+Spec, -File, +Src, +Options`)  
Find named source file from `Spec`, relative to `Src`.
