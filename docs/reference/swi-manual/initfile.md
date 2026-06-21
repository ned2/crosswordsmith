
## 2.2 The user's initialisation file

After the system initialisation, the system consults (see [consult/1](consulting.html#consult/1)) the user's *init* file. This file is searched using [absolute_file_name/3](files.html#absolute_file_name/3) using the path alias (see [file_search_path/2](consulting.html#file_search_path/2)) `app_config`. This is a directory named `swi-prolog` below the OS default name for placing application configuration data:

- On Windows, the CSIDL folder `CSIDL_APPDATA`, typically `C:\Documents and Settings\username\Application Data`.
- If the environment variable `XDG_DATA_HOME` is set, use this. This follows the [free desktop](https://standards.freedesktop.org) standard.
- The expansion of ` /.config`.

The directory can be found using this call:

``` code
?- absolute_file_name(app_config(.), Dir, [file_type(directory)]).
Dir = '/home/jan/.config/swi-prolog'.
```

After the first startup file is found it is loaded and Prolog stops looking for further startup files. The name of the startup file can be changed with the‘**-f** `file`’option. If `File` denotes an absolute path, this file is loaded, otherwise the file is searched for using the same conventions as for the default startup file. Finally, if `file` is `none`, no file is loaded.

The installation provides a file `customize/init.pl` with (commented) commands that are often used to customize the behaviour of Prolog, such as interfacing to the editor, color selection or history parameters. Many of the development tools provide menu entries for editing the startup file and starting a fresh startup file from the system skeleton.

See also the **-s** (script) and **-F** (system-wide initialisation) in [section 2.4](cmdline.html#sec:2.4) and [section 2.3](initgoal.html#sec:2.3).
