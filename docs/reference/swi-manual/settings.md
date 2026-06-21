
## A.53 library(settings): Setting management

author  
Jan Wielemaker

See also  
`library(config)` distributed with XPCE provides an alternative aimed at graphical applications.

This library allows management of configuration settings for Prolog applications. Applications define settings in one or multiple files using the directive [setting/4](settings.html#setting/4) as illustrated below:

``` code
:- use_module(library(settings)).

:- setting(version, atom,   '1.0', 'Current version').
:- setting(timeout, number,    20, 'Timeout in seconds').
```

The directive is subject to [term_expansion/2](consulting.html#term_expansion/2), which guarantees proper synchronisation of the database if source-files are reloaded. This implies it is **not** possible to call [setting/4](settings.html#setting/4) as a predicate.

Settings are local to a module. This implies they are defined in a two-level namespace. Managing settings per module greatly simplifies assembling large applications from multiple modules that configuration through settings. This settings management library ensures proper access, loading and saving of settings.

\[det\]**setting**(`:Name, +Type, +Default, +Comment`)  
Define a setting. `Name` denotes the name of the setting, `Type` its type. `Default` is the value before it is modified. `Default` can refer to environment variables and can use arithmetic expressions as defined by eval_default/4.

If a second declaration for a setting is encountered, it is ignored if `Type` and `Default` are the same. Otherwise a permission_error is raised.

|  |  |
|----|----|
| `Name` | `Name` of the setting (an atom) |
| `Type` | `Type` for setting. One of `any` or a type defined by [must_be/2](error.html#must_be/2). |
| `Default` | `Default` value for the setting. |
| `Comment` | Atom containing a (short) descriptive note. |

\[nondet\]**setting**(`:Name, ?Value`)  
True when `Name` is a currently defined setting with `Value`. Note that `setting(Name, Value)` only enumerates the settings of the current module. All settings can be enumerated using `setting(Module:Name, Value)`. This predicate is `det` if `Name` is ground.

Errors  
`existence_error(setting, Name)`

\[det\]**env**(`+Name:atom, -Value:number`)  
\[det\]**env**(`+Name:atom, +Default:number, -Value:number`)  
Evaluate environment variables on behalf of arithmetic expressions.

\[det\]**set_setting**(`:Name, +Value`)  
Change a setting. Performs existence and type-checking for the setting. If the effective value of the setting is changed it broadcasts the event below.

``` code
settings(changed(Module:Name, Old, New))
```

Note that modified settings are **not** automatically persistent. The application should call [save_settings/0](settings.html#save_settings/0) to persist the changes.

Errors  
\- `existence_error(setting, Name)`  
- `type_error(Type, Value)`

\[det\]**restore_setting**(`:Name`)  
Restore the value of setting `Name` to its default. Broadcast a change like [set_setting/2](settings.html#set_setting/2) if the current value is not the default.

\[det\]**set_setting_default**(`:Name, +Default`)  
Change the default for a setting. The effect is the same as [set_setting/2](settings.html#set_setting/2), but the new value is considered the default when saving and restoring a setting. It is intended to change application defaults in a particular context.

\[det\]**load_settings**(`File`)  
\[det\]**load_settings**(`File, +Options`)  
Load local settings from `File`. Succeeds if `File` does not exist, setting the default save-file to `File`. `Options` are:

**undefined**(`+Action`)  
Define how to handle settings that are not defined. When `error`, an error is printed and the setting is ignored. when `load`, the setting is loaded anyway, waiting for a definition.

If possibly changed settings need to be persistent, the application must call [save_settings/0](settings.html#save_settings/0) as part of its shutdown. In simple cases calling `at_halt(save_settings)` is sufficient.

\[semidet\]**save_settings**  
\[semidet\]**save_settings**(`+File`)  
Save modified settings to `File`. Fails silently if the settings file cannot be written. The [save_settings/0](settings.html#save_settings/0) only attempts to save the settings file if some setting was modified using [set_setting/2](settings.html#set_setting/2).

Errors  
`context_error(settings, no_default_file)` for [save_settings/0](settings.html#save_settings/0) if no default location is known.

\[nondet\]**current_setting**(`?Setting`)  
True if `Setting` is a currently defined setting

\[det\]**setting_property**(`+Setting, +Property`)  
\[nondet\]**setting_property**(`?Setting, ?Property`)  
Query currently defined settings. `Property` is one of

**comment**(`-Atom`)  
**type**(`-Type`)  
`Type` of the setting.

**default**(`-Default`)  
`Default` value. If this is an expression, it is evaluated.

**source**(`-File:-Line`)  
Location where the setting is defined.

\[det\]**list_settings**  
\[det\]**list_settings**(`+Module`)  
List settings to `current_output`. The second form only lists settings on the matching module.

To be done  
Compute the required column widths

**convert_setting_text**(`+Type, +Text, -Value`)  
Converts from textual form to Prolog `Value`. Used to convert values obtained from the environment. Public to provide support in user-interfaces to this library.

Errors  
`type_error(Type, Value)`
