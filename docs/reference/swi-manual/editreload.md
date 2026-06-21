
## 3.3 The test-edit-reload cycle

SWI-Prolog does not enforce the use of a particular editor for writing Prolog source code. Editors are complicated programs that must be mastered in detail for real productive programming. If you are familiar with a specific editor you should not be forced to change. You may specify your favourite editor using the Prolog flag [editor](flags.html#flag:editor), the environment variable `EDITOR` or by defining rules for [prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1).

The use of a built-in editor, which is selected by setting the Prolog flag [editor](flags.html#flag:editor) to `pce_emacs`, has advantages. The XPCE *editor* object, around which the built-in PceEmacs is built, can be opened as a Prolog stream allowing analysis of your source by the real Prolog system.

### 3.3.1 Locating things to edit

The central predicate for editing something is [edit/1](edit.html#edit/1), an extensible front-end that searches for objects (files, predicates, modules, as well as XPCE classes and methods) in the Prolog database. If multiple matches are found it provides a choice. Together with the built-in completion on atoms bound to the `TAB` key this provides a quick way to edit objects:

``` code
?- edit(country).
Please select item to edit:

  1 chat:country/10  '/home/jan/.config/swi-prolog/lib/chat/countr.pl':16
  2 chat:country/1   '/home/jan/.config/swi-prolog/lib/chat/world0.pl':72

Your choice?
```

### 3.3.2 Editing and incremental compilation

One of the nice features of Prolog is that the code can be modified while the program is running. Using pure Prolog you can trace a program, find it is misbehaving, enter a *break environment*, modify the source code, reload it and finally do *retry* on the misbehaving predicate and try again. This sequence is not uncommon for long-running programs. For faster programs one will normally abort after understanding the misbehaviour, edit the source, reload it and try again.

One of the nice features of SWI-Prolog is the availability of [make/0](consulting.html#make/0), a simple predicate that checks all loaded source files to see which ones you have modified. It then reloads these files, considering the module from which the file was loaded originally. This greatly simplifies the trace-edit-verify development cycle. For example, after the tracer reveals there is something wrong with prove/3 , you do:

``` code
?- edit(prove).
```

Now edit the source, possibly switching to other files and making multiple changes. After finishing, invoke [make/0](consulting.html#make/0), either through the editor UI (**Compile/Make** (`Control-C Control-M`)) or on the top level, and watch the files being reloaded.^(46Watching these files is a good habit. If expected files are not reloaded you may have forgotten to save them from the editor or you may have been editing the wrong file (wrong directory).)

``` code
?- make.
% show compiled into photo_gallery 0.03 sec, 3,360 bytes
```
