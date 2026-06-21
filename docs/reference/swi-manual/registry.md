
## A.51 library(registry): Manipulating the Windows registry

The `library(registry)` is only available on the MS-Windows version of SWI-Prolog. It loads the foreign extension `plregtry.dll`, providing the predicates described below. This library only makes the most common operations on the registry available through the Prolog user. The underlying DLL provides a more complete coverage of the Windows registry API. Please consult the sources in `pl/src/win32/foreign/plregtry.c` for further details.

In all these predicates, `Path` refers to a‘/’separated path into the registry. This is *not* an atom containing‘/’-characters as used for filenames, but a term using the functor `/``/2`. Windows defines the following roots for the registry: `classes_root`, `current_user`, `local_machine` and `users`.

**registry_get_key**(`+Path, -Value`)  
Get the principal (default) value associated to this key. Fails silently if the key does not exist.

**registry_get_key**(`+Path, +Name, -Value`)  
Get a named value associated to this key.

**registry_set_key**(`+Path, +Value`)  
Set the principal (default) value of this key. Creates (a path to) the key if it does not already exist.

**registry_set_key**(`+Path, +Name, +Value`)  
Associate a named value to this key. Creates (a path to) the key if it does not already exist.

**registry_delete_key**(`+Path`)  
Delete the indicated key.

**shell_register_file_type**(`+Ext, +Type, +Name, +OpenAction`)  
Register a file-type. `Ext` is the extension to associate. `Type` is the type name, often something like `prolog.type`. `Name` is the name visible in the Windows file-type browser. Finally, `OpenAction` defines the action to execute when a file with this extension is opened in the Windows explorer.

**shell_register_dde**(`+Type, +Action, +Service, +Topic, +Command, +IfNotRunning`)  
Associate DDE actions to a type. `Type` is the same type as used for the 2nd argument of [shell_register_file_type/4](registry.html#shell_register_file_type/4), `Action` is the action to perform, `Service` and `Topic` specify the DDE topic to address, and `Command` is the command to execute on this topic. Finally, `IfNotRunning` defines the command to execute if the required DDE server is not present.

**shell_register_prolog**(`+Ext`)  
Default registration of SWI-Prolog, which is invoked as part of the initialisation process on Windows systems. As the source also includes the above predicates, it is given as an example:

``` code
shell_register_prolog(Ext) :-
        current_prolog_flag(argv, [Me|_]),
        atomic_list_concat(['"', Me, '" "%1"'], OpenCommand),
        shell_register_file_type(
            Ext, 'prolog.type', 'Prolog Source', OpenCommand),
        shell_register_dde(
            'prolog.type', consult,
            prolog, control, 'consult(''%1'')', Me),
        shell_register_dde(
            'prolog.type', edit,
            prolog, control, 'edit(''%1'')', Me).
```
