
## 14.4 Using program resources

A *resource* is similar to a file. Resources, however, can be represented in two different formats: on files, as well as part of the resource *archive* of a saved state (see [qsave_program/2](saved-states.html#qsave_program/2)) that acts as a *virtual file system* for the SWI-Prolog I/O predicates (see [open/4](IO.html#open/4), [register_iri_scheme/3](IO.html#register_iri_scheme/3)).

A resource has a *name*. The *source* data of a resource is a file. Resources are declared by adding clauses to the predicate [resource/2](program-resources.html#resource/2) or [resource/3](program-resources.html#resource/3). Resources can be accessed from Prolog as files that start with `res://` or they can be opened using [open_resource/3](program-resources.html#open_resource/3).

### 14.4.1 Resources as files

As of SWI-Prolog 7.7.13, resources that are compiled into the program can be accessed using the normal file handling predicates. Currently the following predicates transparently handle resources as read-only files:

- [open/3](IO.html#open/3), [open/4](IO.html#open/4)
- [access_file/2](files.html#access_file/2)
- [exists_file/1](files.html#exists_file/1)
- [exists_directory/1](files.html#exists_directory/1)
- [time_file/2](files.html#time_file/2)
- [size_file/2](files.html#size_file/2)

In addition, [open_shared_object/3](foreignlink.html#open_shared_object/3), underlying [use_foreign_library/1](foreignlink.html#use_foreign_library/1) handles *shared objects* or DLLs by copying them to a temporary file and opening this file. If the OS allows for it, the copied file is deleted immediately, otherwise it is deleted on program termination.

With the ability to open resources as if they were files we can use them for many tasks without changing the source code as required when using [open_resource/2](program-resources.html#open_resource/2). Below we describe a typical scenario.

- Related resources are placed in one or more directories. Consider a web application where we have several directories holding icons. Add clauses to [file_search_path/2](consulting.html#file_search_path/2) that makes all icons accessible using the term `icon(file)`.

- Add a clause as below before creating the state. This causes all icons to be become available as `res://app/icon/``file`.

  ``` code
  resource(app/icon, icon(.)).
  ```

- Add a clause to [file_search_path/2](consulting.html#file_search_path/2) that make the icons available from the resource data. For example using the code below.

  ``` code
  :- asserta(user:file_search_path(icon, 'res://app/icon').
  ```

### 14.4.2 Access resources using open_resource

Before the system had the ability to open resources as files, resources were opened using the predicates [open_resource/2](program-resources.html#open_resource/2) or [open_resource/3](program-resources.html#open_resource/3). These predicates provide somewhat better dynamic control over resources depending on whether the code is running from files or from a saved state. The main disadvantage is that having a separate open call requires rewriting code to make it work with resources rather than files.

**open_resource**(`+Name, -Stream`)  
**open_resource**(`+Name, -Stream, +Options`)  
Opens the resource specified by `Name`. If successful, `Stream` is unified with an input stream that provides access to the resource. The stream can be tuned using the `Options`, which is a subset of the options provided by [open/4](IO.html#open/4).

**type**(`Type`)  
**encoding**(`Encoding`)  
**bom**(`Bool`)  
Options that determine the binary/text type, encoding for text streams and whether or not the content should be checked for a BOM marker. The options have the same meaning as the corresponding options for [open/4](IO.html#open/4).

The predicate [open_resource/3](program-resources.html#open_resource/3) first checks [resource/2](program-resources.html#resource/2). When successful it will open the returned resource source file. Otherwise it will look in the program's resource database. When creating a saved state, the system normally saves the resource contents into the resource archive, but does not save the resource clauses.

This way, the development environment uses the files (and modifications) to the [resource/3](program-resources.html#resource/3) declarations and/or files containing resource info, thus immediately affecting the running environment, while the runtime system quickly accesses the system resources.

### 14.4.3 Declaring resources

**resource**(`:Name, +FileSpec`)  
**resource**(`:Name, +FileSpec, +Options`)  
These predicates are defined as dynamic predicates in the module `user`. Clauses for them may be defined in any module, including the user module. `Name` is the name of the resource (an atom). A resource name may contain any character, except for \$ and :, which are reserved for internal usage by the resource library. `FileSpec` is a file specification that may exploit [file_search_path/2](consulting.html#file_search_path/2) (see [absolute_file_name/2](files.html#absolute_file_name/2)).

Often, resources are defined as unit clauses (facts), but the definition of this predicate also allows for rules. For proper generation of the saved state, it must be possible to enumerate the available resources by calling this predicate with all its arguments unbound.

If `FileSpec` points at a directory, the content of the directory is recursively added below `Name`. If `FileSpec` a term of the form `Alias(Name)`, all directories that match this specification are enumerated and their content is added to the resource database. If an file appears in multiple results of this search path only the first file is added. Note that this is consistent with the normal behaviour where [absolute_file_name/3](files.html#absolute_file_name/3) returns the first match. The `Options` can be used to control what is saved from a directory.

**include**(`+Patterns`)  
Only include a file from a directory if it matches at least one of the members of `Patterns`.

**exclude**(`+Patterns`)  
Excludes a file from a directory if it matches at least one of the members of `Patterns`.

### 14.4.4 Managing resource files

As of version 7.7.13, SWI-Prolog resource files are zip(1) files. Prolog creates and accesses its resource files using the [minizip](http://www.winimage.com/zLibDll/minizip.html) project. The resource files may be examined and modified using any tool that can process zip files.
