
## 3.9 Summary of the IDE

The SWI-Prolog development environment consists of a number of interrelated but not (yet) integrated tools. Here is a list of the most important features and tips.

- *Atom completion*  
  The console completes a partial atom on the `TAB` key and shows alternatives on the command `Alt-?`.
- *Use [edit/1](edit.html#edit/1) for finding locations*  
  The command [edit/1](edit.html#edit/1) takes the name of a file, module, predicate or other entity registered through extensions and starts the user's preferred editor at the right location.
- *Select editor*  
  External editors are selected using the `EDITOR` environment variable, by setting the Prolog flag [editor](flags.html#flag:editor), or by defining the hook [prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1).
- *Update Prolog after editing*  
  Using [make/0](consulting.html#make/0), all files you have edited are re-loaded.
- *PceEmacs*  
  Offers syntax highlighting and checking based on real-time parsing of the editor's buffer, layout support and navigation support.
- *Using the graphical debugger*  
  The predicates [guitracer/0](guitracer.html#guitracer/0) and [noguitracer/0](guitracer.html#noguitracer/0) switch between traditional text-based and window-based debugging. The tracer is activated using the [trace/0](debugger.html#trace/0), [spy/1](debugger.html#spy/1) or menu items from PceEmacs or the Prolog Navigator.
- *The Prolog Navigator*  
  Shows the file structure and structure inside the file. It allows for loading files, editing, setting spy points, etc.
