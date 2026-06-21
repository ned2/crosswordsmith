
## 2.5 UI Themes

UI (colour) themes play a role in two parts: when writing to the *console* and for the xpce-based development tools such as PceEmacs or the graphical debugger. Coloured console output is based on [ansi_format/3](ansiterm.html#ansi_format/3). The central message infra structure based on [print_message/2](printmsg.html#print_message/2) labels message (components) with a Prolog term that specifies the role. This is mapped to concrete colours by means of the hook [prolog:console_color/2](ansiterm.html#prolog:console_color/2). Theming the IDE uses xpce *class variables* that are initialised from Prolog when xpce is loaded.

Themes are implemented as a Prolog file in the file search path library/theme. A theme can be loaded using (for example) the directive below in the user's initialization file (see [section 2.2](initfile.html#sec:2.2)).

``` code
:- use_module(library(theme/dark)).
```

The theme file `library(theme/auto)` is provided to automatically choose a reasonable theme based on the environment. The current version detects the background color on *xterm* compatible terminal emulators (found on most Unix systems) and loads the `dark` theme if the background is‘darkish’.

The following notes apply to the different platforms on which SWI-Prolog is supported:

**Unix/Linux**  
If an xterm compatible terminal emulator is used to run Prolog you may wish to load either an explicit theme or `library(theme/auto)`.

**Epilog Prolog consoles**  
The **swipl-win** graphical application can be themed by loading a theme file. The theme file also sets the foreground and background colours for the Epilog console.

### 2.5.1 Status of theme support

Theme support was added in SWI-Prolog 8.1.11. Only part of the IDE tools are covered and the only additional theme (`dark`) is not net well balanced. The interfaces between the theme file and notably the IDE components is not very well established. Please contribute by improving the `dark` theme. Once that is complete and properly functioning we can start adding new themes.
