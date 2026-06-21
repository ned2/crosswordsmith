
## 3.1 The project source files

Organisation of source files depends largely on the size of your project. If you are doing exercises for a Prolog course you'll normally use one file for each exercise. If you have a small project you'll work with one directory holding a couple of files and some files to link it all together. Even bigger projects will be organised in sub-projects, each using its own directory.

### 3.1.1 File Names and Locations

#### 3.1.1.1 File Name Extensions

The first consideration is what extension to use for the source files. Tradition calls for `.pl`, but conflicts with Perl force the use of another extension on systems where extensions have global meaning, such as MS-Windows. On such systems `.pro` is the common alternative. On MS-Windows, the alternative extension is stored in the registry key `HKEY_CURRENT_USER/Software/SWI/Prolog/fileExtension` or `HKEY_LOCAL_MACHINE/Software/SWI/Prolog/fileExtension`. All versions of SWI-Prolog load files with the extension `.pl` as well as with the registered alternative extension without explicitly specifying the extension. For portability reasons we propose the following convention:

**If there is no conflict**  
because you do not use a conflicting application or the system does not force a unique relation between extension and application, use `.pl`.

**With a conflict**  
choose `.pro` and use this extension for the files you want to load through your file manager. Use `.pl` for all other files for maximal portability.

#### 3.1.1.2 Project Directories

Large projects are generally composed of sub-projects, each using its own directory or directory structure. If nobody else will ever touch your files and you use only one computer, there is little to worry about, but this is rarely the case with a large project.

To improve portability, SWI-Prolog uses the POSIX notation for filenames, which uses the forward slash (`/`) to separate directories. Just before reaching the file system, SWI-Prolog uses [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2) to convert the filename to the conventions used by the hosting operating system. It is *strongly* advised to write paths using the `/`, especially on systems using the `\` for this purpose (MS-Windows). Using `\` violates the portability rules and requires you to *double* the `\` due to the Prolog quoted-atom escape rules.

Portable code should use [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2) to convert computed paths into system paths when constructing commands for [shell/1](system.html#shell/1) and friends.

#### 3.1.1.3 Sub-projects using search paths

Thanks to Quintus, Prolog adapted an extensible mechanism for searching files using [file_search_path/2](consulting.html#file_search_path/2). This mechanism allows for comfortable and readable specifications.

Suppose you have extensive library packages on graph algorithms, set operations and GUI primitives. These sub-projects are likely candidates for re-use in future projects. A good choice is to create a directory with sub-directories for each of these sub-projects.

Next, there are three options. One is to add the sub-projects to the directory hierarchy of the current project. Another is to use a completely dislocated directory. Third, the sub-project can be added to the SWI-Prolog hierarchy. Using local installation, a typical [file_search_path/2](consulting.html#file_search_path/2) is:

``` code
:- prolog_load_context(directory, Dir),
   asserta(user:file_search_path(myapp, Dir)).

user:file_search_path(graph, myapp(graph)).
user:file_search_path(ui,    myapp(ui)).
```

When using sub-projects in the SWI-Prolog hierarchy, one should use the path alias `swi` as basis. For a system-wide installation, use an absolute path.

Extensive sub-projects with a small well-defined API should define a load file with calls to [use_module/1](import.html#use_module/1) to import the various library components and export the API.

### 3.1.2 Project Special Files

There are a number of tasks you typically carry out on your project, such as loading it, creating a saved state, debugging it, etc. Good practice on large projects is to define small files that hold the commands to execute such a task, name this file after the task and give it a file extension that makes starting easy (see [section 3.1.1.1](projectfiles.html#sec:3.1.1.1)). The task *load* is generally central to these tasks. Here is a tentative list:

- *`load.pl`*  
  Use this file to set up the environment (Prolog flags and file search paths) and load the sources. Quite commonly this file also provides convenient predicates to parse command line options and start the application.
- *`run.pl`*  
  Use this file to start the application. Normally it loads `load.pl` in silent-mode, and calls one of the starting predicates from `load.pl`.
- *`save.pl`*  
  Use this file to create a saved state of the application by loading `load.pl` and calling [qsave_program/2](saved-states.html#qsave_program/2) to generate a saved state with the proper options.
- *`debug.pl`*  
  Loads the program for debugging. In addition to loading `load.pl` this file defines rules for [portray/1](termrw.html#portray/1) to modify printing rules for complex terms and customisation rules for the debugger and editing environment. It may start some of these tools.

### 3.1.3 International source files

As discussed in [section 2.18](widechars.html#sec:2.18), SWI-Prolog supports international character handling. Its internal encoding is UNICODE. I/O streams convert to/from this internal format. This section discusses the options for source files not in US-ASCII.

SWI-Prolog can read files in any of the encodings described in [section 2.18](widechars.html#sec:2.18). Two encodings are of particular interest. The `text` encoding deals with the current *locale*, the default used by this computer for representing text files. The encodings `utf8`, `unicode_le` and `unicode_be` are *UNICODE* encodings: they can represent---in the same file---characters of virtually any known language. In addition, they do so unambiguously.

If one wants to represent non US-ASCII text as Prolog terms in a source file, there are several options:

- *Use escape sequences*  
  This approach describes NON-ASCII as sequences of the form `\`*octal*`\`. The numerical argument is interpreted as a UNICODE character.^(45To my knowledge, the ISO escape sequence is limited to 3 octal digits, which means most characters cannot be represented.) The resulting Prolog file is strict 7-bit US-ASCII, but if there are many NON-ASCII characters it becomes very unreadable.
- *Use local conventions*  
  Alternatively the file may be specified using local conventions, such as the EUC encoding for Japanese text. The disadvantage is portability. If the file is moved to another machine, this machine must use the same *locale* or the file is unreadable. There is no elegant way if files from multiple locales must be united in one application using this technique. In other words, it is fine for local projects in countries with uniform locale conventions.
- *Using UTF-8 files*  
  The best way to specify source files with many NON-ASCII characters is definitely the use of UTF-8 encoding. Prolog can be notified of this encoding in two ways, using a UTF-8 *BOM* (see [section 2.18.1.1](widechars.html#sec:2.18.1.1)) or using the directive `:- encoding(utf8).` Many of today's text editors, including PceEmacs, are capable of editing UTF-8 files. Projects that were started using local conventions can be re-coded using the Unix **iconv** tool or often using commands offered by the editor.
