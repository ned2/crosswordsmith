
## 2.13 An overview of hook predicates

SWI-Prolog provides a large number of hooks, mainly to control handling messages, debugging, startup, shut-down, macro-expansion, etc. Below is a summary of all defined hooks with an indication of their portability.

- *[portray/1](termrw.html#portray/1)*  
  Hook into [write_term/3](termrw.html#write_term/3) to alter the way terms are printed (ISO).
- *[message_hook/3](printmsg.html#message_hook/3)*  
  Hook into [print_message/2](printmsg.html#print_message/2) to alter the way system messages are printed (Quintus/SICStus).
- *[message_property/2](printmsg.html#message_property/2)*  
  Hook into [print_message/2](printmsg.html#print_message/2) that defines prefix, output stream, color, etc.
- *message_prefix_hook/2*  
  Hook into [print_message/2](printmsg.html#print_message/2) to add additional prefixes to the message such as the time and thread.
- *[library_directory/1](consulting.html#library_directory/1)*  
  Hook into [absolute_file_name/3](files.html#absolute_file_name/3) to define new library directories (most Prolog systems).
- *[file_search_path/2](consulting.html#file_search_path/2)*  
  Hook into [absolute_file_name/3](files.html#absolute_file_name/3) to define new search paths (Quintus/SICStus).
- *[term_expansion/2](consulting.html#term_expansion/2)*  
  Hook into [load_files/2](consulting.html#load_files/2) to modify read terms before they are compiled (macro-processing) (most Prolog systems).
- *[goal_expansion/2](consulting.html#goal_expansion/2)*  
  Same as [term_expansion/2](consulting.html#term_expansion/2) for individual goals (SICStus).
- *[prolog_load_file/2](loadfilehook.html#prolog_load_file/2)*  
  Hook into [load_files/2](consulting.html#load_files/2) to load other data formats for Prolog sources from‘non-file’resources. The [load_files/2](consulting.html#load_files/2) predicate is the ancestor of [consult/1](consulting.html#consult/1), [use_module/1](import.html#use_module/1), etc.
- *[prolog_edit:locate/3](edit.html#prolog_edit:locate/3)*  
  Hook into [edit/1](edit.html#edit/1) to locate objects (SWI).
- *[prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1)*  
  Hook into [edit/1](edit.html#edit/1) to call an internal editor (SWI).
- *[prolog_edit:edit_command/2](edit.html#prolog_edit:edit_command/2)*  
  Hook into [edit/1](edit.html#edit/1) to define the external editor to use (SWI).
- *[prolog_list_goal/1](intlibs.html#prolog_list_goal/1)*  
  Hook into the tracer to list the code associated to a particular goal (SWI).
- *[prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4)*  
  Hook into the tracer to handle trace events (SWI).
- *[prolog:debug_control_hook/1](prologdebug.html#prolog:debug_control_hook/1)*  
  Hook in [spy/1](debugger.html#spy/1), [nospy/1](debugger.html#nospy/1), [nospyall/0](debugger.html#nospyall/0) and [debugging/0](debugger.html#debugging/0) to extend these control predicates to higher-level libraries.
- *[prolog:help_hook/1](intlibs.html#prolog:help_hook/1)*  
  Hook in [help/0](online-help.html#help/0), [help/1](online-help.html#help/1) and [apropos/1](online-help.html#apropos/1) to extend the help system.
- *[resource/3](program-resources.html#resource/3)*  
  Define a new resource (not really a hook, but similar) (SWI).
- *[exception/3](exception3.html#exception/3)*  
  Old attempt to a generic hook mechanism. Handles undefined predicates (SWI).
- *[attr_unify_hook/2](attvar.html#attr_unify_hook/2)*  
  Unification hook for attributed variables. Can be defined in any module. See [section 8.1](attvar.html#sec:8.1) for details.
