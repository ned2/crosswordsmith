
## F.1 Predicates

The predicate summary is used by the Prolog predicate [apropos/1](online-help.html#apropos/1) to suggest predicates from a keyword.

[@/2](overrule.html#@/2)

Call using calling context

[!/0](control.html#!/0)

Cut (discard choicepoints)

[\$/0](debug-determinism.html#$/0)

Discard choicepoints and demand deterministic success

[\$/1](debug-determinism.html#$/1)

Verify goal succeeds deterministically

[,/2](control.html#,/2)

Conjunction of goals

[-\>/2](control.html#-%3E/2)

If-then-else

[\*-\>/2](control.html#*-%3E/2)

Soft-cut

./2

Consult. Also functional notation

[:\</2](bidicts.html#:%3C/2)

Select keys from a dict

[:=/2](wasm-js-call.html#:=/2)

WASM: Call JavaScript

[;/2](control.html#;/2)

Disjunction of two goals

[\</2](arith.html#%3C/2)

Arithmetic smaller

[=/2](compare.html#=/2)

True when arguments are unified

[=../2](manipterm.html#=../2)

“Univ.” Term to list conversion

[=:=/2](arith.html#=:=/2)

Arithmetic equality

[=\</2](arith.html#=%3C/2)

Arithmetic smaller or equal

[==/2](compare.html#==/2)

Test for strict equality

[=@=/2](compare.html#=@=/2)

Test for structural equality (variant)

[=\\/2](arith.html#=\=/2)

Arithmetic not equal

[\>/2](arith.html#%3E/2)

Arithmetic larger

[\>=/2](arith.html#%3E=/2)

Arithmetic larger or equal

[\>:\</2](bidicts.html#%3E:%3C/2)

Partial dict unification

[?=/2](compare.html#?=/2)

Test of terms can be compared now

[@\</2](compare.html#@%3C/2)

Standard order smaller

[@=\</2](compare.html#@=%3C/2)

Standard order smaller or equal

[@\>/2](compare.html#@%3E/2)

Standard order larger

[@\>=/2](compare.html#@%3E=/2)

Standard order larger or equal

[\\/1](control.html#\+/1)

Negation by failure. Same as [not/1](metacall.html#not/1)

[\\/2](compare.html#\=/2)

True if arguments cannot be unified

[\\=/2](compare.html#\==/2)

True if arguments are not strictly equal

[\\@=/2](compare.html#\=@=/2)

Not structural identical

^/2

Existential quantification ([bagof/3](allsolutions.html#bagof/3), [setof/3](allsolutions.html#setof/3))

[\|/2](control.html#%7C/2)

Disjunction in DCGs. Same as [;/2](control.html#;/2)

/1

DCG escape; constraints

[abolish/1](db.html#abolish/1)

Remove predicate definition from the database

[abolish/2](db.html#abolish/2)

Remove predicate definition from the database

[abolish_all_tables/0](tabling-preds.html#abolish_all_tables/0)

Abolish computed tables

[abolish_module_tables/1](tabling-preds.html#abolish_module_tables/1)

Abolish all tables in a module

[abolish_monotonic_tables/0](tabling-monotonic.html#abolish_monotonic_tables/0)

Abolish all monotonic tables

[abolish_nonincremental_tables/0](tabling-preds.html#abolish_nonincremental_tables/0)

Abolish non-automatic tables

[abolish_nonincremental_tables/1](tabling-preds.html#abolish_nonincremental_tables/1)

Abolish non-automatic tables

[abolish_private_tables/0](tabling-preds.html#abolish_private_tables/0)

Abolish tables of this thread

[abolish_shared_tables/0](tabling-preds.html#abolish_shared_tables/0)

Abolish tables shared between threads

[abolish_table_subgoals/1](tabling-preds.html#abolish_table_subgoals/1)

Abolish tables for a goal

[abort/0](toplevel.html#abort/0)

Abort execution, return to top level

[absolute_file_name/2](files.html#absolute_file_name/2)

Get absolute path name

[absolute_file_name/3](files.html#absolute_file_name/3)

Get absolute path name with options

[answer_count_restraint/0](tabling-restraints.html#answer_count_restraint/0)

Undefined answer due to `max_answers`

[access_file/2](files.html#access_file/2)

Check access permissions of a file

[acyclic_term/1](typetest.html#acyclic_term/1)

Test term for cycles

[add_import_module/3](importmodule.html#add_import_module/3)

Add module to the auto-import list

[add_nb_set/2](nb_set.html#add_nb_set/2)

Add term to a non-backtrackable set

[add_nb_set/3](nb_set.html#add_nb_set/3)

Add term to a non-backtrackable set

[append/1](IO.html#append/1)

Append to a file

[apple_current_locale_identifier/1](system.html#apple_current_locale_identifier/1)

Get Apple locale info

[apply/2](metacall.html#apply/2)

Call goal with additional arguments

[apropos/1](online-help.html#apropos/1)

`library(online_help)` Search manual

[arg/3](manipterm.html#arg/3)

Access argument of a term

[assoc_to_list/2](assoc.html#assoc_to_list/2)

Convert association tree to list

[assert/1](db.html#assert/1)

Add a clause to the database

[assert/2](db.html#assert/2)

Add a clause to the database, give reference

[asserta/1](db.html#asserta/1)

Add a clause to the database (first)

[asserta/2](db.html#asserta/2)

Add a clause to the database (first)

[assertion/1](debug.html#assertion/1)

Make assertions about your program

[assertz/1](db.html#assertz/1)

Add a clause to the database (last)

[assertz/2](db.html#assertz/2)

Add a clause to the database (last)

[attach_console/0](threadutil.html#attach_console/0)

Attach I/O console to thread

[attach_packs/0](pack-attach.html#attach_packs/0)

Attach add-ons

[attach_packs/1](pack-attach.html#attach_packs/1)

Attach add-ons from directory

[attach_packs/2](pack-attach.html#attach_packs/2)

Attach add-ons from directory

attribute_goals/3

Project attributes to goals

[attr_unify_hook/2](attvar.html#attr_unify_hook/2)

Attributed variable unification hook

[attr_portray_hook/2](attvar.html#attr_portray_hook/2)

Attributed variable print hook

[attvar/1](attvar.html#attvar/1)

Type test for attributed variable

[at_end_of_stream/0](chario.html#at_end_of_stream/0)

Test for end of file on input

[at_end_of_stream/1](chario.html#at_end_of_stream/1)

Test for end of file on stream

[at_halt/1](consulting.html#at_halt/1)

Register goal to run at [halt/1](toplevel.html#halt/1)

[atom/1](typetest.html#atom/1)

Type check for an atom

[atom_chars/2](manipatom.html#atom_chars/2)

Convert between atom and list of characters

[atom_codes/2](manipatom.html#atom_codes/2)

Convert between atom and list of characters codes

[atom_concat/3](manipatom.html#atom_concat/3)

Concatenate two atoms

[atom_length/2](manipatom.html#atom_length/2)

Determine length of an atom

[atom_number/2](manipatom.html#atom_number/2)

Convert between atom and number

[atom_prefix/2](manipatom.html#atom_prefix/2)

Test for start of atom

[atom_string/2](string.html#atom_string/2)

Conversion between atom and string

[atom_to_term/3](manipatom.html#atom_to_term/3)

Convert between atom and term

[atomic/1](typetest.html#atomic/1)

Type check for primitive

[atomic_concat/3](manipatom.html#atomic_concat/3)

Concatenate two atomic values to an atom

[atomic_list_concat/2](manipatom.html#atomic_list_concat/2)

Append a list of atomics

[atomic_list_concat/3](manipatom.html#atomic_list_concat/3)

Append a list of atomics with separator

[atomics_to_string/2](string.html#atomics_to_string/2)

Concatenate list of inputs to a string

[atomics_to_string/3](string.html#atomics_to_string/3)

Concatenate list of inputs to a string

[autoload/1](module-autoload.html#autoload/1)

Declare a file for autoloading

[autoload/2](module-autoload.html#autoload/2)

Declare a file for autoloading specific predicates

[autoload_all/0](saved-states.html#autoload_all/0)

Autoload all predicates now

[autoload_call/1](metacall.html#autoload_call/1)

Call after autoloading

[autoload_path/1](autoload.html#autoload_path/1)

Add directories for autoloading

[await/2](wasm-js-call.html#await/2)

WASM: Wait for a Promise

[b_getval/2](gvar.html#b_getval/2)

Fetch backtrackable global variable

[b_set_dict/3](bidicts.html#b_set_dict/3)

Destructive assignment on a dict

[b_setval/2](gvar.html#b_setval/2)

Assign backtrackable global variable

[bagof/3](allsolutions.html#bagof/3)

Find all solutions to a goal

[between/3](arith.html#between/3)

Integer range checking/generating

[blob/2](typetest.html#blob/2)

Type check for a blob

[bounded_number/3](arith.html#bounded_number/3)

Number between bounds

[break/0](toplevel.html#break/0)

Start interactive top level

break_hook/6

*(hook)* Debugger hook

[byte_count/2](streamstat.html#byte_count/2)

Byte-position in a stream

[call/1](metacall.html#call/1)

Call a goal

[call/\[2..\]](metacall.html#call/2)

Call with additional arguments

[call_cleanup/2](metacall.html#call_cleanup/2)

Guard a goal with a cleanup-handler

[call_dcg/3](DCG.html#call_dcg/3)

As [phrase/3](DCG.html#phrase/3) without type checking

[call_delays/2](WFS.html#call_delays/2)

Get the condition associated with an answer

[call_residue_vars/2](coroutining.html#call_residue_vars/2)

Find residual attributed variables

[call_residual_program/2](WFS.html#call_residual_program/2)

Get residual program associated with an answer

[call_shared_object_function/2](foreignlink.html#call_shared_object_function/2)

UNIX: Call C-function in shared (.so) file

[call_with_depth_limit/3](metacall.html#call_with_depth_limit/3)

Prove goal with bounded depth

[call_with_inference_limit/3](metacall.html#call_with_inference_limit/3)

Prove goal in limited inferences

[callable/1](typetest.html#callable/1)

Test for atom or compound term

[cancel_halt/1](consulting.html#cancel_halt/1)

Cancel [halt/0](toplevel.html#halt/0) from an [at_halt/1](consulting.html#at_halt/1) hook

[catch/3](exception.html#catch/3)

Call goal, watching for exceptions

[char_code/2](manipatom.html#char_code/2)

Convert between character and character code

[char_conversion/2](charconv.html#char_conversion/2)

Provide mapping of input characters

[char_type/2](chartype.html#char_type/2)

Classify characters

[character_count/2](streamstat.html#character_count/2)

Get character index on a stream

[chdir/1](files.html#chdir/1)

Compatibility: change working directory

[chr_constraint/1](practical.html#chr_constraint/1)

CHR Constraint declaration

[chr_show_store/1](chr-debugging.html#chr_show_store/1)

List suspended CHR constraints

[chr_trace/0](chr-debugging.html#chr_trace/0)

Start CHR tracer

[chr_type/1](practical.html#chr_type/1)

CHR Type declaration

[chr_notrace/0](chr-debugging.html#chr_notrace/0)

Stop CHR tracer

[chr_leash/1](chr-debugging.html#chr_leash/1)

Define CHR leashed ports

[chr_option/2](chr-syntaxandsemantics.html#chr_option/2)

Specify CHR compilation options

[clause/2](examineprog.html#clause/2)

Get clauses of a predicate

[clause/3](examineprog.html#clause/3)

Get clauses of a predicate

[clause_property/2](examineprog.html#clause_property/2)

Get properties of a clause

[close/1](IO.html#close/1)

Close stream

[close/2](IO.html#close/2)

Close stream (forced)

[close_dde_conversation/1](DDE.html#close_dde_conversation/1)

Win32: Close DDE channel

[close_shared_object/1](foreignlink.html#close_shared_object/1)

UNIX: Close shared library (.so file)

[collation_key/2](chartype.html#collation_key/2)

Sort key for locale dependent ordering

comment_hook/3

*(hook)* handle comments in sources

[compare/3](compare.html#compare/3)

Compare, using a predicate to determine the order

[compile_aux_clauses/1](consulting.html#compile_aux_clauses/1)

Compile predicates for [goal_expansion/2](consulting.html#goal_expansion/2)

[compile_predicates/1](dynamic.html#compile_predicates/1)

Compile dynamic code to static

[compiling/0](consulting.html#compiling/0)

Is this a compilation run?

[compound/1](typetest.html#compound/1)

Test for compound term

[compound_name_arity/3](manipterm.html#compound_name_arity/3)

Name and arity of a compound term

[compound_name_arguments/3](manipterm.html#compound_name_arguments/3)

Name and arguments of a compound term

[code_type/2](chartype.html#code_type/2)

Classify a character-code

[consult/1](consulting.html#consult/1)

Read (compile) a Prolog source file

[context_module/1](ctxmodule.html#context_module/1)

Get context module of current goal

convert_time/8

Break time stamp into fields

convert_time/2

Convert time stamp to string

[copy_stream_data/2](chario.html#copy_stream_data/2)

Copy all data from stream to stream

[copy_stream_data/3](chario.html#copy_stream_data/3)

Copy n bytes from stream to stream

[copy_predicate_clauses/2](db.html#copy_predicate_clauses/2)

Copy clauses between predicates

[copy_term/2](manipterm.html#copy_term/2)

Make a copy of a term

[copy_term/3](attvar.html#copy_term/3)

Copy a term and obtain attribute-goals

[copy_term/4](manipterm.html#copy_term/4)

Copy part of the variables in a term

[copy_term_nat/2](attvar.html#copy_term_nat/2)

Make a copy of a term without attributes

[copy_term_nat/4](manipterm.html#copy_term_nat/4)

Copy part of the variables in a term

[create_prolog_flag/3](flags.html#create_prolog_flag/3)

Create a new Prolog flag

[current_arithmetic_function/1](miscarith.html#current_arithmetic_function/1)

Examine evaluable functions

[current_atom/1](examineprog.html#current_atom/1)

Examine existing atoms

[current_blob/2](examineprog.html#current_blob/2)

Examine typed blobs

[current_char_conversion/2](charconv.html#current_char_conversion/2)

Query input character mapping

[current_engine/1](engine-predicates.html#current_engine/1)

Enumerate known engines

[current_flag/1](examineprog.html#current_flag/1)

Examine existing flags

[current_foreign_library/2](foreignlink.html#current_foreign_library/2)

`library(shlib)` Examine loaded shared libraries (.so files)

[current_format_predicate/2](format.html#current_format_predicate/2)

Enumerate user-defined format codes

[current_functor/2](examineprog.html#current_functor/2)

Examine existing name/arity pairs

[current_input/1](IO.html#current_input/1)

Get current input stream

[current_key/1](examineprog.html#current_key/1)

Examine existing database keys

[current_locale/1](locale.html#current_locale/1)

Get the current locale

[current_module/1](manipmodule.html#current_module/1)

Examine existing modules

[current_op/3](operators.html#current_op/3)

Examine current operator declarations

[current_output/1](IO.html#current_output/1)

Get the current output stream

[current_predicate/1](examineprog.html#current_predicate/1)

Examine existing predicates (ISO)

[current_predicate/2](examineprog.html#current_predicate/2)

Examine existing predicates

[current_signal/3](signal.html#current_signal/3)

Current software signal mapping

[current_stream/3](IO.html#current_stream/3)

Examine open streams

[current_table/2](tabling-preds.html#current_table/2)

Find answer table for a variant

[current_transaction/1](db.html#current_transaction/1)

Detect encapsulating transactions

[current_trie/1](db.html#current_trie/1)

Enumerate known tries

[cyclic_term/1](typetest.html#cyclic_term/1)

Test term for cycles

[day_of_the_week/2](system.html#day_of_the_week/2)

Determine ordinal-day from date

[date_time_stamp/2](system.html#date_time_stamp/2)

Convert date structure to time-stamp

[date_time_value/3](system.html#date_time_value/3)

Extract info from a date structure

[dcg_translate_rule/2](consulting.html#dcg_translate_rule/2)

Source translation of DCG rules

[dcg_translate_rule/4](consulting.html#dcg_translate_rule/4)

Source translation of DCG rules

[dde_current_connection/2](DDE.html#dde_current_connection/2)

Win32: Examine open DDE connections

[dde_current_service/2](DDE.html#dde_current_service/2)

Win32: Examine DDE services provided

[dde_execute/2](DDE.html#dde_execute/2)

Win32: Execute command on DDE server

[dde_register_service/2](DDE.html#dde_register_service/2)

Win32: Become a DDE server

[dde_request/3](DDE.html#dde_request/3)

Win32: Make a DDE request

dde_poke/3

Win32: POKE operation on DDE server

[dde_unregister_service/1](DDE.html#dde_unregister_service/1)

Win32: Terminate a DDE service

[debug/0](debugger.html#debug/0)

Test for debugging mode

[debug/1](debug.html#debug/1)

Select topic for debugging

[debug/3](debug.html#debug/3)

Print debugging message on topic

debug_control_hook/1

*(hook)* Extend [spy/1](debugger.html#spy/1), etc.

debug_reset_from_class/0

Restore debug mode

[debugging/0](debugger.html#debugging/0)

Show debugger status

[debugging/1](debug.html#debugging/1)

Test where we are debugging topic

[default_module/2](importmodule.html#default_module/2)

Query module inheritance

[del_attr/2](attvar.html#del_attr/2)

Delete attribute from variable

[del_attrs/1](attvar.html#del_attrs/1)

Delete all attributes from variable

[del_dict/4](bidicts.html#del_dict/4)

Delete Key-Value pair from a dict

[delays_residual_program/2](WFS.html#delays_residual_program/2)

Get the residual program for an answer

[delete_directory/1](files.html#delete_directory/1)

Remove a folder from the file system

[delete_file/1](files.html#delete_file/1)

Remove a file from the file system

[delete_import_module/2](importmodule.html#delete_import_module/2)

Remove module from import list

[det/1](debug-determinism.html#det/1)

Declare predicates as deterministic

[deterministic/1](manipstack.html#deterministic/1)

Test determinicity of current clause

[dict_create/3](bidicts.html#dict_create/3)

Create a dict from data

[dict_pairs/3](bidicts.html#dict_pairs/3)

Convert between dict and list of pairs

[dict_same_keys/2](bidicts.html#dict_same_keys/2)

True when dicts have the same keys

[dif/2](coroutining.html#dif/2)

Constrain two terms to be different

[directory_files/2](files.html#directory_files/2)

Get entries of a directory/folder

[discontiguous/1](dynamic.html#discontiguous/1)

Indicate distributed definition of a predicate

[divmod/4](arith.html#divmod/4)

Compute quotient and remainder of two integers

[downcase_atom/2](chartype.html#downcase_atom/2)

Convert atom to lower-case

[duplicate_term/2](manipterm.html#duplicate_term/2)

Create a copy of a term

[dwim_match/2](miscpreds.html#dwim_match/2)

Atoms match in “Do What I Mean” sense

[dwim_match/3](miscpreds.html#dwim_match/3)

Atoms match in “Do What I Mean” sense

[dwim_predicate/2](examineprog.html#dwim_predicate/2)

Find predicate in “Do What I Mean” sense

[dynamic/1](dynamic.html#dynamic/1)

Indicate predicate definition may change

[dynamic/2](dynamic.html#dynamic/2)

Indicate predicate definition may change

[edit/0](edit.html#edit/0)

Edit current script- or associated file

[edit/1](edit.html#edit/1)

Edit a file, predicate, module (extensible)

[elif/1](consulting.html#elif/1)

Part of conditional compilation (directive)

[else/0](consulting.html#else/0)

Part of conditional compilation (directive)

[empty_assoc/1](assoc.html#empty_assoc/1)

Create/test empty association tree

[empty_nb_set/1](nb_set.html#empty_nb_set/1)

Test/create an empty non-backtrackable set

[encoding/1](consulting.html#encoding/1)

Define encoding inside a source file

[endif/0](consulting.html#endif/0)

End of conditional compilation (directive)

[engine_create/3](engine-predicates.html#engine_create/3)

Create an interactor

[engine_create/4](engine-predicates.html#engine_create/4)

Create an interactor

[engine_destroy/1](engine-predicates.html#engine_destroy/1)

Destroy an interactor

[engine_fetch/1](engine-predicates.html#engine_fetch/1)

Get term from caller

[engine_next/2](engine-predicates.html#engine_next/2)

Ask interactor for next term

[engine_next_reified/2](engine-predicates.html#engine_next_reified/2)

Ask interactor for next term

[engine_post/2](engine-predicates.html#engine_post/2)

Send term to an interactor

[engine_post/3](engine-predicates.html#engine_post/3)

Send term to an interactor and wait for reply

[engine_self/1](engine-predicates.html#engine_self/1)

Get handle to running interactor

[engine_yield/1](engine-predicates.html#engine_yield/1)

Make term available to caller

[ensure_loaded/1](consulting.html#ensure_loaded/1)

Consult a file if that has not yet been done

[epilog/0](system.html#epilog/0)

Create a Prolog console

[epilog/1](system.html#epilog/1)

Create a Prolog console

[erase/1](db.html#erase/1)

Erase a database record or clause

[exception/3](exception3.html#exception/3)

*(hook)* Handle runtime exceptions

[exists_directory/1](files.html#exists_directory/1)

Check existence of directory

[exists_file/1](files.html#exists_file/1)

Check existence of file

[exists_source/1](consulting.html#exists_source/1)

Check existence of a Prolog source

[exists_source/2](consulting.html#exists_source/2)

Check existence of a Prolog source

expand_answer/2

Expand answer of query (deprecated)

expand_answer/3

Expand answer of query

[expand_file_name/2](files.html#expand_file_name/2)

Wildcard expansion of file names

[expand_file_search_path/2](consulting.html#expand_file_search_path/2)

Wildcard expansion of file paths

[expand_goal/2](consulting.html#expand_goal/2)

Compiler: expand goal in clause-body

[expand_goal/4](consulting.html#expand_goal/4)

Compiler: expand goal in clause-body

expand_query/4

Expanded entered query

[expand_term/2](consulting.html#expand_term/2)

Compiler: expand read term into clause(s)

[expand_term/4](consulting.html#expand_term/4)

Compiler: expand read term into clause(s)

[expects_dialect/1](dialect.html#expects_dialect/1)

For which Prolog dialect is this code written?

[explain/1](online-help.html#explain/1)

`library(explain)` Explain argument

[explain/2](online-help.html#explain/2)

`library(explain)` 2nd argument is explanation of first

[export/1](altmoduleapi.html#export/1)

Export a predicate from a module

[fail/0](control.html#fail/0)

Always false

[false/0](control.html#false/0)

Always false

[fast_term_serialized/2](IO.html#fast_term_serialized/2)

Fast term (de-)serialization

[fast_read/2](IO.html#fast_read/2)

Read binary term serialization

[fast_write/2](IO.html#fast_write/2)

Write binary term serialization

[current_prolog_flag/2](flags.html#current_prolog_flag/2)

Get system configuration parameters

[file_base_name/2](files.html#file_base_name/2)

Get file part of path

[file_directory_name/2](files.html#file_directory_name/2)

Get directory part of path

[file_name_extension/3](files.html#file_name_extension/3)

Add, remove or test file extensions

[file_search_path/2](consulting.html#file_search_path/2)

Define path-aliases for locating files

[find_chr_constraint/1](chr-debugging.html#find_chr_constraint/1)

Returns a constraint from the store

[findall/3](allsolutions.html#findall/3)

Find all solutions to a goal

[findall/4](allsolutions.html#findall/4)

Difference list version of [findall/3](allsolutions.html#findall/3)

[findnsols/4](allsolutions.html#findnsols/4)

Find first `N` solutions

[findnsols/5](allsolutions.html#findnsols/5)

Difference list version of [findnsols/4](allsolutions.html#findnsols/4)

[fill_buffer/1](chario.html#fill_buffer/1)

Fill the input buffer of a stream

[flag/3](db.html#flag/3)

Simple global variable system

[float/1](typetest.html#float/1)

Type check for a floating point number

[float_class/2](arith.html#float_class/2)

Classify (special) floats

[float_parts/4](arith.html#float_parts/4)

Get mantissa and exponent of a float

[flush_output/0](chario.html#flush_output/0)

Output pending characters on current stream

[flush_output/1](chario.html#flush_output/1)

Output pending characters on specified stream

[forall/2](forall2.html#forall/2)

Prove goal for all solutions of another goal

[format/1](format.html#format/1)

Formatted output

[format/2](format.html#format/2)

Formatted output with arguments

[format/3](format.html#format/3)

Formatted output on a stream

[format_time/3](system.html#format_time/3)

C **strftime()** like date/time formatter

[format_time/4](system.html#format_time/4)

date/time formatter with explicit locale

[format_predicate/2](format.html#format_predicate/2)

Program [format/\[1,2\]](format.html#format/1)

[term_attvars/2](attvar.html#term_attvars/2)

Find attributed variables in a term

[term_variables/2](manipterm.html#term_variables/2)

Find unbound variables in a term

[term_variables/3](manipterm.html#term_variables/3)

Find unbound variables in a term

[text_to_string/2](string.html#text_to_string/2)

Convert arbitrary text to a string

[freeze/2](coroutining.html#freeze/2)

Delay execution until variable is bound

[frozen/2](coroutining.html#frozen/2)

Query delayed goals on var

[functor/3](manipterm.html#functor/3)

Get name and arity of a term or construct a term

[functor/4](manipterm.html#functor/4)

Get name and arity of a term or construct a term

[garbage_collect/0](memory.html#garbage_collect/0)

Invoke the garbage collector

[garbage_collect_atoms/0](memory.html#garbage_collect_atoms/0)

Invoke the atom garbage collector

[garbage_collect_clauses/0](memory.html#garbage_collect_clauses/0)

Invoke clause garbage collector

[gen_assoc/3](assoc.html#gen_assoc/3)

Enumerate members of association tree

[gen_nb_set/2](nb_set.html#gen_nb_set/2)

Generate members of non-backtrackable set

[gensym/2](gensym.html#gensym/2)

Generate unique atoms from a base

[get/1](chario.html#get/1)

Read first non-blank character

[get/2](chario.html#get/2)

Read first non-blank character from a stream

[get_assoc/3](assoc.html#get_assoc/3)

Fetch key from association tree

[get_assoc/5](assoc.html#get_assoc/5)

Fetch key from association tree

[get_attr/3](attvar.html#get_attr/3)

Fetch named attribute from a variable

[get_attrs/2](attvar.html#get_attrs/2)

Fetch all attributes of a variable

[get_byte/1](chario.html#get_byte/1)

Read next byte (ISO)

[get_byte/2](chario.html#get_byte/2)

Read next byte from a stream (ISO)

[get_char/1](chario.html#get_char/1)

Read next character as an atom (ISO)

[get_char/2](chario.html#get_char/2)

Read next character from a stream (ISO)

[get_code/1](chario.html#get_code/1)

Read next character (ISO)

[get_code/2](chario.html#get_code/2)

Read next character from a stream (ISO)

[get_dict/3](bidicts.html#get_dict/3)

Get the value associated to a key from a dict

[get_dict/5](bidicts.html#get_dict/5)

Replace existing value in a dict

[get_flag/2](db.html#get_flag/2)

Get value of a flag

[get_single_char/1](chario.html#get_single_char/1)

Read next character from the terminal

[get_string_code/3](string.html#get_string_code/3)

Get character code at index in string

[get_time/1](system.html#get_time/1)

Get current time

[get0/1](chario.html#get0/1)

Read next character

[get0/2](chario.html#get0/2)

Read next character from a stream

[getenv/2](system.html#getenv/2)

Get shell environment variable

[goal_expansion/2](consulting.html#goal_expansion/2)

Hook for macro-expanding goals

[goal_expansion/4](consulting.html#goal_expansion/4)

Hook for macro-expanding goals

[ground/1](typetest.html#ground/1)

Verify term holds no unbound variables

[gdebug/0](guitracer.html#gdebug/0)

Debug using graphical tracer

[gspy/1](guitracer.html#gspy/1)

Spy using graphical tracer

[gtrace/0](guitracer.html#gtrace/0)

Trace using graphical tracer

[guitracer/0](guitracer.html#guitracer/0)

Install hooks for the graphical debugger

[gxref/0](xref.html#gxref/0)

Cross-reference loaded program

[halt/0](toplevel.html#halt/0)

Exit from Prolog

[halt/1](toplevel.html#halt/1)

Exit from Prolog with status

[term_hash/2](db.html#term_hash/2)

Hash-value of ground term

[term_hash/4](db.html#term_hash/4)

Hash-value of term with depth limit

[help/0](online-help.html#help/0)

Give help on help

[help/1](online-help.html#help/1)

Give help on predicates and show parts of manual

help_hook/1

*(hook)* User-hook in the help-system

[if/1](consulting.html#if/1)

Start conditional compilation (directive)

[ignore/1](metacall.html#ignore/1)

Call the argument, but always succeed

[import/1](altmoduleapi.html#import/1)

Import a predicate from a module

[import_module/2](importmodule.html#import_module/2)

Query import modules

[in_pce_thread/1](mt-xpce.html#in_pce_thread/1)

Run goal in XPCE thread

[in_pce_thread_sync/1](mt-xpce.html#in_pce_thread_sync/1)

Run goal in XPCE thread

[include/1](consulting.html#include/1)

Include a file with declarations

[initialization/1](consulting.html#initialization/1)

Initialization directive

[initialization/2](consulting.html#initialization/2)

Initialization directive

[initialize/0](consulting.html#initialize/0)

Run program initialization

[instance/2](db.html#instance/2)

Fetch clause or record from reference

[integer/1](typetest.html#integer/1)

Type check for integer

interactor/0

Start new thread with console and top level

[is/2](arith.html#is/2)

Evaluate arithmetic expression

[is_absolute_file_name/1](files.html#is_absolute_file_name/1)

True if arg defines an absolute path

[is_assoc/1](assoc.html#is_assoc/1)

Verify association list

[is_async/0](wasm-js-call.html#is_async/0)

WASM: Test Prolog can call [await/2](wasm-js-call.html#await/2)

[is_dict/1](bidicts.html#is_dict/1)

Type check for a dict

[is_dict/2](bidicts.html#is_dict/2)

Type check for a dict in a class

[is_engine/1](engine-predicates.html#is_engine/1)

Type check for an engine handle

[is_list/1](builtinlist.html#is_list/1)

Type check for a list

[is_message_queue/1](threadcom.html#is_message_queue/1)

Type check for a message queue

[is_most_general_term/1](manipterm.html#is_most_general_term/1)

Type check for general term

[is_object/1](wasm-js-call.html#is_object/1)

WASM: Test JavaScript object

[is_object/2](wasm-js-call.html#is_object/2)

WASM: Test JavaScript object and class

[is_stream/1](IO.html#is_stream/1)

Type check for a stream handle

[is_trie/1](db.html#is_trie/1)

Type check for a trie handle

[is_thread/1](thmonitor.html#is_thread/1)

Type check for an thread handle

[join_threads/0](threadutil.html#join_threads/0)

Join all terminated threads interactively

[keysort/2](builtinlist.html#keysort/2)

Sort, using a key

[known_licenses/0](softlicense.html#known_licenses/0)

Print known licenses

[last/2](lists.html#last/2)

Last element of a list

[leash/1](debugger.html#leash/1)

Change ports visited by the tracer

[length/2](builtinlist.html#length/2)

Length of a list

[library_directory/1](consulting.html#library_directory/1)

*(hook)* Directories holding Prolog libraries

[license/0](softlicense.html#license/0)

Evaluate licenses of loaded modules

[license/1](softlicense.html#license/1)

Define license for current file

[license/2](softlicense.html#license/2)

Define license for named module

[line_count/2](streamstat.html#line_count/2)

Line number on stream

[line_position/2](streamstat.html#line_position/2)

Character position in line on stream

[list_debug_topics/0](debug.html#list_debug_topics/0)

List registered topics for debugging

[list_to_assoc/2](assoc.html#list_to_assoc/2)

Create association tree from list

[list_to_set/2](lists.html#list_to_set/2)

Remove duplicates from a list

[list_strings/0](string.html#list_strings/0)

Help porting to version 7

[load_files/1](consulting.html#load_files/1)

Load source files

[load_files/2](consulting.html#load_files/2)

Load source files with options

[load_foreign_library/1](foreignlink.html#load_foreign_library/1)

`library(shlib)` Load shared library (.so file)

[load_foreign_library/2](foreignlink.html#load_foreign_library/2)

`library(shlib)` Load shared library (.so file)

[locale_create/3](locale.html#locale_create/3)

Create a new locale object

[locale_destroy/1](locale.html#locale_destroy/1)

Destroy a locale object

[locale_property/2](locale.html#locale_property/2)

Query properties of locale objects

[locale_sort/2](chartype.html#locale_sort/2)

Language dependent sort of atoms

[make/0](consulting.html#make/0)

Reconsult all changed source files

[make_directory/1](files.html#make_directory/1)

Create a folder on the file system

[make_library_index/1](autoload.html#make_library_index/1)

Create autoload file INDEX.pl

[malloc_property/1](memory.html#malloc_property/1)

Property of the allocator

[make_library_index/2](autoload.html#make_library_index/2)

Create selective autoload file INDEX.pl

[map_assoc/2](assoc.html#map_assoc/2)

Map association tree

[map_assoc/3](assoc.html#map_assoc/3)

Map association tree

[max_assoc/3](assoc.html#max_assoc/3)

Highest key in association tree

[memberchk/2](builtinlist.html#memberchk/2)

Deterministic [member/2](lists.html#member/2)

message_action/2

*(hook)* Associate actions with messages

[message_hook/3](printmsg.html#message_hook/3)

Intercept [print_message/2](printmsg.html#print_message/2)

message_line_element/2

*(hook)* Intercept [print_message_lines/3](printmsg.html#print_message_lines/3)

[message_property/2](printmsg.html#message_property/2)

*(hook)* Define display of a message

[message_queue_create/1](threadcom.html#message_queue_create/1)

Create queue for thread communication

[message_queue_create/2](threadcom.html#message_queue_create/2)

Create queue for thread communication

[message_queue_destroy/1](threadcom.html#message_queue_destroy/1)

Destroy queue for thread communication

[message_queue_property/2](threadcom.html#message_queue_property/2)

Query message queue properties

[message_queue_set/2](threadcom.html#message_queue_set/2)

Set a message queue property

[message_to_string/2](printmsg.html#message_to_string/2)

Translate message-term to string

[meta_predicate/1](metapred.html#meta_predicate/1)

Declare access to other predicates

[min_assoc/3](assoc.html#min_assoc/3)

Lowest key in association tree

[mode/1](dynamic.html#mode/1)

module

1Query/set current type-in module

[module/2](defmodule.html#module/2)

Declare a module

[module/3](defmodule.html#module/3)

Declare a module with language options

[module_property/2](manipmodule.html#module_property/2)

Find properties of a module

[module_transparent/1](ctxmodule.html#module_transparent/1)

Indicate module based meta-predicate

[msort/2](builtinlist.html#msort/2)

Sort, do not remove duplicates

[multifile/1](dynamic.html#multifile/1)

Indicate distributed definition of predicate

[mutex_create/1](threadsync.html#mutex_create/1)

Create a thread-synchronisation device

[mutex_create/2](threadsync.html#mutex_create/2)

Create a thread-synchronisation device

[mutex_destroy/1](threadsync.html#mutex_destroy/1)

Destroy a mutex

[mutex_lock/1](threadsync.html#mutex_lock/1)

Become owner of a mutex

[mutex_property/2](threadsync.html#mutex_property/2)

Query mutex properties

[mutex_statistics/0](thmonitor.html#mutex_statistics/0)

Print statistics on mutex usage

[mutex_trylock/1](threadsync.html#mutex_trylock/1)

Become owner of a mutex (non-blocking)

[mutex_unlock/1](threadsync.html#mutex_unlock/1)

Release ownership of mutex

[mutex_unlock_all/0](threadsync.html#mutex_unlock_all/0)

Release ownership of all mutexes

[name/2](manipatom.html#name/2)

Convert between atom and list of character codes

[nb_current/2](gvar.html#nb_current/2)

Enumerate non-backtrackable global variables

[nb_delete/1](gvar.html#nb_delete/1)

Delete a non-backtrackable global variable

[nb_getval/2](gvar.html#nb_getval/2)

Fetch non-backtrackable global variable

[nb_link_dict/3](bidicts.html#nb_link_dict/3)

Non-backtrackable assignment to dict

[nb_linkarg/3](manipterm.html#nb_linkarg/3)

Non-backtrackable assignment to term

[nb_linkval/2](gvar.html#nb_linkval/2)

Assign non-backtrackable global variable

[nb_set_to_list/2](nb_set.html#nb_set_to_list/2)

Convert non-backtrackable set to list

[nb_set_dict/3](bidicts.html#nb_set_dict/3)

Non-backtrackable assignment to dict

[nb_setarg/3](manipterm.html#nb_setarg/3)

Non-backtrackable assignment to term

[nb_setval/2](gvar.html#nb_setval/2)

Assign non-backtrackable global variable

[nl/0](chario.html#nl/0)

Generate a newline

[nl/1](chario.html#nl/1)

Generate a newline on a stream

[nodebug/0](debugger.html#nodebug/0)

Disable debugging

[nodebug/1](debug.html#nodebug/1)

Disable debug-topic

[noguitracer/0](guitracer.html#noguitracer/0)

Disable the graphical debugger

[nonground/2](manipterm.html#nonground/2)

Term is not ground due to witness

[nonvar/1](typetest.html#nonvar/1)

Type check for bound term

nonterminal/1

Set predicate property

noprofile/1

Hide (meta-) predicate for the profiler

[noprotocol/0](protocol.html#noprotocol/0)

Disable logging of user interaction

[normalize_space/2](chartype.html#normalize_space/2)

Normalize white space

[nospy/1](debugger.html#nospy/1)

Remove spy point

[nospyall/0](debugger.html#nospyall/0)

Remove all spy points

[not/1](metacall.html#not/1)

Negation by failure (argument not provable). Same as [\\/1](control.html#\+/1)

[not_exists/1](tabling-preds.html#not_exists/1)

Tabled negation for non-ground or non-tabled goals

[notrace/0](debugger.html#notrace/0)

Stop tracing

[notrace/1](debugger.html#notrace/1)

Do not debug argument goal

[nth_clause/3](examineprog.html#nth_clause/3)

N-th clause of a predicate

[nth_integer_root_and_remainder/4](arith.html#nth_integer_root_and_remainder/4)

Integer root and remainder

[number/1](typetest.html#number/1)

Type check for integer or float

[number_chars/2](manipatom.html#number_chars/2)

Convert between number and one-char atoms

[number_codes/2](manipatom.html#number_codes/2)

Convert between number and character codes

[number_string/2](string.html#number_string/2)

Convert between number and string

[numbervars/3](manipterm.html#numbervars/3)

Number unbound variables of a term

[numbervars/4](manipterm.html#numbervars/4)

Number unbound variables of a term

[on_signal/3](signal.html#on_signal/3)

Handle a software signal

[once/1](metacall.html#once/1)

Call a goal deterministically

[op/3](operators.html#op/3)

Declare an operator

[open/3](IO.html#open/3)

Open a file (creating a stream)

[open/4](IO.html#open/4)

Open a file (creating a stream)

[open_dde_conversation/3](DDE.html#open_dde_conversation/3)

Win32: Open DDE channel

[open_null_stream/1](IO.html#open_null_stream/1)

Open a stream to discard output

[open_resource/3](program-resources.html#open_resource/3)

Open a program resource as a stream

[open_shared_object/2](foreignlink.html#open_shared_object/2)

UNIX: Open shared library (.so file)

[open_shared_object/3](foreignlink.html#open_shared_object/3)

UNIX: Open shared library (.so file)

open_source_hook/3

*(hook)* Open a source file

[open_string/2](string.html#open_string/2)

Open a string as a stream

[ord_list_to_assoc/2](assoc.html#ord_list_to_assoc/2)

Convert ordered list to assoc

[parse_time/2](system.html#parse_time/2)

Parse text to a time-stamp

[parse_time/3](system.html#parse_time/3)

Parse text to a time-stamp

[pce_dispatch/1](mt-xpce.html#pce_dispatch/1)

Run XPCE GUI in separate thread

pce_call/1

Run goal in XPCE GUI thread

[peek_byte/1](chario.html#peek_byte/1)

Read byte without removing

[peek_byte/2](chario.html#peek_byte/2)

Read byte without removing

[peek_char/1](chario.html#peek_char/1)

Read character without removing

[peek_char/2](chario.html#peek_char/2)

Read character without removing

[peek_code/1](chario.html#peek_code/1)

Read character-code without removing

[peek_code/2](chario.html#peek_code/2)

Read character-code without removing

[peek_string/3](chario.html#peek_string/3)

Read a string without removing

[phrase/2](DCG.html#phrase/2)

Activate grammar-rule set

[phrase/3](DCG.html#phrase/3)

Activate grammar-rule set (returning rest)

[phrase_from_quasi_quotation/2](quasiquotations.html#phrase_from_quasi_quotation/2)

Parse quasi quotation with DCG

please/3

Query/change environment parameters

[plus/3](arith.html#plus/3)

Logical integer addition

[portray/1](termrw.html#portray/1)

*(hook)* Modify behaviour of [print/1](termrw.html#print/1)

[predicate_property/2](examineprog.html#predicate_property/2)

Query predicate attributes

[predsort/3](builtinlist.html#predsort/3)

Sort, using a predicate to determine the order

[print/1](termrw.html#print/1)

Print a term

[print/2](termrw.html#print/2)

Print a term on a stream

[print_message/2](printmsg.html#print_message/2)

Print message from (exception) term

[print_message_lines/3](printmsg.html#print_message_lines/3)

Print message to stream

[profile/1](profile.html#profile/1)

Obtain execution statistics

[profile/2](profile.html#profile/2)

Obtain execution statistics

profile_count/3

Obtain profile results on a predicate

profiler/2

Obtain/change status of the profiler

[prolog/0](toplevel.html#prolog/0)

Run interactive top level

[prolog_alert_signal/2](signal.html#prolog_alert_signal/2)

Query/set unblock signal

[prolog_choice_attribute/3](manipstack.html#prolog_choice_attribute/3)

Examine the choice point stack

[prolog_current_choice/1](manipstack.html#prolog_current_choice/1)

Reference to most recent choice point

[prolog_current_frame/1](manipstack.html#prolog_current_frame/1)

Reference to goal's environment stack

[prolog_cut_to/1](ancestral-cut.html#prolog_cut_to/1)

Realise global cuts

[prolog_edit:locate/2](edit.html#prolog_edit:locate/2)

Locate targets for [edit/1](edit.html#edit/1)

[prolog_edit:locate/3](edit.html#prolog_edit:locate/3)

Locate targets for [edit/1](edit.html#edit/1)

[prolog_edit:edit_source/1](edit.html#prolog_edit:edit_source/1)

Call editor for [edit/1](edit.html#edit/1)

[prolog_edit:edit_command/2](edit.html#prolog_edit:edit_command/2)

Specify editor activation

[prolog_edit:load/0](edit.html#prolog_edit:load/0)

Load [edit/1](edit.html#edit/1) extensions

prolog_exception_hook/5

Rewrite exceptions

[prolog_file_type/2](consulting.html#prolog_file_type/2)

Define meaning of file extension

[prolog_frame_attribute/3](manipstack.html#prolog_frame_attribute/3)

Obtain information on a goal environment

[prolog_ide/1](idepreds.html#prolog_ide/1)

Program access to the development environment

[prolog_interrupt/0](interrupt.html#prolog_interrupt/0)

Allow debugging a thread

[prolog_list_goal/1](intlibs.html#prolog_list_goal/1)

*(hook)* Intercept tracer’L’command

[prolog_listen/2](prolog-event.html#prolog_listen/2)

Listen to Prolog events

[prolog_listen/3](prolog-event.html#prolog_listen/3)

Listen to Prolog events

[prolog_load_context/2](consulting.html#prolog_load_context/2)

Context information for directives

[prolog_load_file/2](loadfilehook.html#prolog_load_file/2)

*(hook)* Program [load_files/2](consulting.html#load_files/2)

[prolog_skip_level/2](tracehook.html#prolog_skip_level/2)

Indicate deepest recursion to trace

[prolog_stack_property/2](memory.html#prolog_stack_property/2)

Query properties of the stacks

[prolog_to_os_filename/2](files.html#prolog_to_os_filename/2)

Convert between Prolog and OS filenames

[prolog_trace_interception/4](tracehook.html#prolog_trace_interception/4)

`library(user)` Intercept the Prolog tracer

[prolog_unlisten/2](prolog-event.html#prolog_unlisten/2)

Stop listening to Prolog events

[project_attributes/2](attvar.html#project_attributes/2)

Project constraints to query variables

[prompt1/1](termrw.html#prompt1/1)

Change prompt for 1 line

[prompt/2](termrw.html#prompt/2)

Change the prompt used by [read/1](termrw.html#read/1)

[protocol/1](protocol.html#protocol/1)

Make a log of the user interaction

[protocola/1](protocol.html#protocola/1)

Append log of the user interaction to file

[protocolling/1](protocol.html#protocolling/1)

On what file is user interaction logged

[public/1](dynamic.html#public/1)

Declaration that a predicate may be called

[put/1](chario.html#put/1)

Write a character

[put/2](chario.html#put/2)

Write a character on a stream

[put_assoc/4](assoc.html#put_assoc/4)

Add Key-Value to association tree

[put_attr/3](attvar.html#put_attr/3)

Put attribute on a variable

[put_attrs/2](attvar.html#put_attrs/2)

Set/replace all attributes on a variable

[put_byte/1](chario.html#put_byte/1)

Write a byte

[put_byte/2](chario.html#put_byte/2)

Write a byte on a stream

[put_char/1](chario.html#put_char/1)

Write a character

[put_char/2](chario.html#put_char/2)

Write a character on a stream

[put_code/1](chario.html#put_code/1)

Write a character-code

[put_code/2](chario.html#put_code/2)

Write a character-code on a stream

[put_dict/3](bidicts.html#put_dict/3)

Add/replace multiple keys in a dict

[put_dict/4](bidicts.html#put_dict/4)

Add/replace a single key in a dict

[qcompile/1](consulting.html#qcompile/1)

Compile source to Quick Load File

[qcompile/2](consulting.html#qcompile/2)

Compile source to Quick Load File

[qsave_program/1](saved-states.html#qsave_program/1)

Create runtime application

[qsave_program/2](saved-states.html#qsave_program/2)

Create runtime application

[quasi_quotation_syntax/1](quasiquotations.html#quasi_quotation_syntax/1)

Declare quasi quotation syntax

[quasi_quotation_syntax_error/1](quasiquotations.html#quasi_quotation_syntax_error/1)

Raise syntax error

[radial_restraint/0](tabling-restraints.html#radial_restraint/0)

Tabling radial restraint was violated

[random_property/1](miscarith.html#random_property/1)

Query properties of random generation

[rational/1](typetest.html#rational/1)

Type check for a rational number

[rational/3](typetest.html#rational/3)

Decompose a rational

[read/1](termrw.html#read/1)

Read Prolog term

[read/2](termrw.html#read/2)

Read Prolog term from stream

[read_clause/3](termrw.html#read_clause/3)

Read clause from stream

[read_link/3](files.html#read_link/3)

Read a symbolic link

[read_pending_codes/3](chario.html#read_pending_codes/3)

Fetch buffered input from a stream

[read_pending_chars/3](chario.html#read_pending_chars/3)

Fetch buffered input from a stream

[read_string/3](string.html#read_string/3)

Read a number of characters into a string

[read_string/5](string.html#read_string/5)

Read string up to a delimiter

[read_term/2](termrw.html#read_term/2)

Read term with options

[read_term/3](termrw.html#read_term/3)

Read term with options from stream

[read_term_from_atom/3](termrw.html#read_term_from_atom/3)

Read term with options from atom

[read_term_with_history/2](termrw.html#read_term_with_history/2)

Read term with command line history

[recorda/2](db.html#recorda/2)

Record term in the database (first)

[recorda/3](db.html#recorda/3)

Record term in the database (first)

[recorded/2](db.html#recorded/2)

Obtain term from the database

[recorded/3](db.html#recorded/3)

Obtain term from the database

[recordz/2](db.html#recordz/2)

Record term in the database (last)

[recordz/3](db.html#recordz/3)

Record term in the database (last)

[redefine_system_predicate/1](db.html#redefine_system_predicate/1)

Abolish system definition

[reexport/1](reexport.html#reexport/1)

Load files and re-export the imported predicates

[reexport/2](reexport.html#reexport/2)

Load predicates from a file and re-export it

[reload_foreign_libraries/0](foreignlink.html#reload_foreign_libraries/0)

Reload DLLs/shared objects

[reload_library_index/0](autoload.html#reload_library_index/0)

Force reloading the autoload index

[rename_file/2](files.html#rename_file/2)

Change name of file

[repeat/0](control.html#repeat/0)

Succeed, leaving infinite backtrack points

[require/1](consulting.html#require/1)

This file requires these predicates

[reset/3](delcont.html#reset/3)

Wrapper for delimited continuations

[reset_gensym/1](gensym.html#reset_gensym/1)

Reset a gensym key

[reset_gensym/0](gensym.html#reset_gensym/0)

Reset all gensym keys

reset_profiler/0

Clear statistics obtained by the profiler

[resource/2](program-resources.html#resource/2)

Declare a program resource

[resource/3](program-resources.html#resource/3)

Declare a program resource

[retract/1](db.html#retract/1)

Remove clause from the database

[retractall/1](db.html#retractall/1)

Remove unifying clauses from the database

[same_file/2](files.html#same_file/2)

Succeeds if arguments refer to same file

[same_term/2](manipterm.html#same_term/2)

Test terms to be at the same address

[see/1](IO.html#see/1)

Change the current input stream

[seeing/1](IO.html#seeing/1)

Query the current input stream

[seek/4](IO.html#seek/4)

Modify the current position in a stream

[seen/0](IO.html#seen/0)

Close the current input stream

select_dict/2

Select matching attributes from a dict

[select_dict/3](bidicts.html#select_dict/3)

Select matching attributes from a dict

[set_end_of_stream/1](chario.html#set_end_of_stream/1)

Set physical end of an open file

[set_epilog/1](system.html#set_epilog/1)

Control the SWI-Prolog console

[set_flag/2](db.html#set_flag/2)

Set value of a flag

[set_input/1](IO.html#set_input/1)

Set current input stream from a stream

[set_locale/1](locale.html#set_locale/1)

Set the default local

[set_malloc/1](memory.html#set_malloc/1)

Set memory allocator property

[set_module/1](manipmodule.html#set_module/1)

Set properties of a module

[set_output/1](IO.html#set_output/1)

Set current output stream from a stream

[set_prolog_IO/3](IO.html#set_prolog_IO/3)

Prepare streams for interactive session

[set_prolog_flag/2](flags.html#set_prolog_flag/2)

Define a system feature

[set_prolog_gc_thread/1](memory.html#set_prolog_gc_thread/1)

Control the gc thread

[set_prolog_stack/2](memory.html#set_prolog_stack/2)

Modify stack characteristics

[set_random/1](miscarith.html#set_random/1)

Control random number generation

[set_stream/2](IO.html#set_stream/2)

Set stream attribute

[set_stream_position/2](IO.html#set_stream_position/2)

Seek stream to position

[set_system_IO/3](IO.html#set_system_IO/3)

Rebind stdin/stderr/stdout

[setup_call_cleanup/3](metacall.html#setup_call_cleanup/3)

Undo side-effects safely

[setup_call_catcher_cleanup/4](metacall.html#setup_call_catcher_cleanup/4)

Undo side-effects safely

[setarg/3](manipterm.html#setarg/3)

Destructive assignment on term

[setenv/2](system.html#setenv/2)

Set shell environment variable

[setlocale/3](system.html#setlocale/3)

Set/query C-library regional information

[setof/3](allsolutions.html#setof/3)

Find all unique solutions to a goal

[shell/1](system.html#shell/1)

Execute OS command

[shell/2](system.html#shell/2)

Execute OS command

[shift/1](delcont.html#shift/1)

Shift control to the closest [reset/3](delcont.html#reset/3)

[shift_for_copy/1](delcont.html#shift_for_copy/1)

Shift control to the closest [reset/3](delcont.html#reset/3)

[show_profile/1](profile.html#show_profile/1)

Show results of the profiler

[sig_atomic/1](threadcom.html#sig_atomic/1)

Run goal without handling signals

[sig_block/1](threadcom.html#sig_block/1)

Block matching thread signals

[sig_pending/1](threadcom.html#sig_pending/1)

Query pending signals

[sig_remove/2](threadcom.html#sig_remove/2)

Remove pending signals

[sig_unblock/1](threadcom.html#sig_unblock/1)

Unblock matching thread signals

[size_abstract_term/3](tabling-restraints.html#size_abstract_term/3)

Abstract a term (tabling support)

[size_file/2](files.html#size_file/2)

Get size of a file in characters

[size_nb_set/2](nb_set.html#size_nb_set/2)

Determine size of non-backtrackable set

[skip/1](chario.html#skip/1)

Skip to character in current input

[skip/2](chario.html#skip/2)

Skip to character on stream

[sleep/1](miscpreds.html#sleep/1)

Suspend execution for specified time

[snapshot/1](db.html#snapshot/1)

Run goal in isolation

[sort/2](builtinlist.html#sort/2)

Sort elements in a list

[sort/4](builtinlist.html#sort/4)

Sort elements in a list

[source_exports/2](dialect.html#source_exports/2)

Check whether source exports a predicate

[source_file/1](consulting.html#source_file/1)

Examine currently loaded source files

[source_file/2](consulting.html#source_file/2)

Obtain source file of predicate

[source_file_property/2](consulting.html#source_file_property/2)

Information about loaded files

[source_location/2](consulting.html#source_location/2)

Location of last read term

[split_string/4](string.html#split_string/4)

Break a string into substrings

[spy/1](debugger.html#spy/1)

Force tracer on specified predicate

[stamp_date_time/3](system.html#stamp_date_time/3)

Convert time-stamp to date structure

[statistics/2](builtin-statistics.html#statistics/2)

Obtain collected statistics

[stream_pair/3](IO.html#stream_pair/3)

Create/examine a bi-directional stream

[stream_position_data/3](IO.html#stream_position_data/3)

Access fields from stream position

[stream_property/2](IO.html#stream_property/2)

Get stream properties

[string/1](typetest.html#string/1)

Type check for string

[string_bytes/3](string.html#string_bytes/3)

Translates between text and bytes in encoding

[string_concat/3](string.html#string_concat/3)

[atom_concat/3](manipatom.html#atom_concat/3) for strings

[string_length/2](string.html#string_length/2)

Determine length of a string

[string_chars/2](string.html#string_chars/2)

Conversion between string and list of characters

[string_codes/2](string.html#string_codes/2)

Conversion between string and list of character codes

[string_code/3](string.html#string_code/3)

Get or find a character code in a string

[string_lower/2](string.html#string_lower/2)

Case conversion to lower case

[string_upper/2](string.html#string_upper/2)

Case conversion to upper case

[string_predicate/1](check.html#string_predicate/1)

*(hook)* Predicate contains strings

[strip_module/3](ctxmodule.html#strip_module/3)

Extract context module and term

[style_check/1](debugger.html#style_check/1)

Change level of warnings

[sub_atom/5](manipatom.html#sub_atom/5)

Take a substring from an atom

[sub_atom_icasechk/3](manipatom.html#sub_atom_icasechk/3)

Case insensitive substring match

[sub_string/5](string.html#sub_string/5)

Take a substring from a string

[subsumes_term/2](compare.html#subsumes_term/2)

One-sided unification test

[succ/2](arith.html#succ/2)

Logical integer successor relation

[tab/1](chario.html#tab/1)

Output number of spaces

[tab/2](chario.html#tab/2)

Output number of spaces on a stream

[table/1](tabling-preds.html#table/1)

Declare predicate to be tabled

[tabled_call/1](tabling-preds.html#tabled_call/1)

Helper for [not_exists/1](tabling-preds.html#not_exists/1)

[tdebug/0](threadutil.html#tdebug/0)

Switch all threads into debug mode

[tdebug/1](threadutil.html#tdebug/1)

Switch a thread into debug mode

[tell/1](IO.html#tell/1)

Change current output stream

[telling/1](IO.html#telling/1)

Query current output stream

[term_expansion/2](consulting.html#term_expansion/2)

*(hook)* Convert term before compilation

[term_expansion/4](consulting.html#term_expansion/4)

*(hook)* Convert term before compilation

[term_singletons/2](manipterm.html#term_singletons/2)

Find singleton variables in a term

[term_string/2](string.html#term_string/2)

Read/write a term from/to a string

[term_string/3](string.html#term_string/3)

Read/write a term from/to a string

[term_subsumer/3](compare.html#term_subsumer/3)

Most specific generalization of two terms

[term_to_atom/2](manipatom.html#term_to_atom/2)

Convert between term and atom

[thread_affinity/3](threadcreate.html#thread_affinity/3)

Query and control the *affinity* mask

[thread_alias/1](threadutil.html#thread_alias/1)

Set the alias name of a thread

[thread_at_exit/1](threadcreate.html#thread_at_exit/1)

Register goal to be called at exit

[thread_create/2](threadcreate.html#thread_create/2)

Create a new Prolog task

[thread_create/3](threadcreate.html#thread_create/3)

Create a new Prolog task

[thread_detach/1](threadcreate.html#thread_detach/1)

Make thread cleanup after completion

[thread_exit/1](threadcreate.html#thread_exit/1)

Terminate Prolog task with value

[thread_get_message/1](threadcom.html#thread_get_message/1)

Wait for message

[thread_get_message/2](threadcom.html#thread_get_message/2)

Wait for message in a queue

[thread_get_message/3](threadcom.html#thread_get_message/3)

Wait for message in a queue

[thread_idle/2](memory.html#thread_idle/2)

Reduce footprint while waiting

[thread_initialization/1](threadcreate.html#thread_initialization/1)

Run action at start of thread

[thread_join/1](threadcreate.html#thread_join/1)

Wait for Prolog task-completion

[thread_join/2](threadcreate.html#thread_join/2)

Wait for Prolog task-completion

[thread_local/1](threadcom.html#thread_local/1)

Declare thread-specific clauses for a predicate

[thread_message_hook/3](printmsg.html#thread_message_hook/3)

Thread local [message_hook/3](printmsg.html#message_hook/3)

[thread_peek_message/1](threadcom.html#thread_peek_message/1)

Test for message

[thread_peek_message/2](threadcom.html#thread_peek_message/2)

Test for message in a queue

[thread_property/2](thmonitor.html#thread_property/2)

Examine Prolog threads

[thread_self/1](threadcreate.html#thread_self/1)

Get identifier of current thread

[thread_send_message/2](threadcom.html#thread_send_message/2)

Send message to another thread

[thread_send_message/3](threadcom.html#thread_send_message/3)

Send message to another thread

[thread_setconcurrency/2](threadcreate.html#thread_setconcurrency/2)

Number of active threads

[thread_signal/2](threadcom.html#thread_signal/2)

Execute goal in another thread

[thread_statistics/3](thmonitor.html#thread_statistics/3)

Get statistics of another thread

[thread_update/2](threadcom.html#thread_update/2)

Update a module and signal waiters

[thread_wait/2](threadcom.html#thread_wait/2)

Wait for a goal to become true

[threads/0](threadutil.html#threads/0)

List running threads

[throw/1](exception.html#throw/1)

Raise an exception (see [catch/3](exception.html#catch/3))

[time/1](statistics.html#time/1)

Determine time needed to execute goal

[time_file/2](files.html#time_file/2)

Get last modification time of file

[tmp_file/2](files.html#tmp_file/2)

Create a temporary filename

[tmp_file_stream/3](files.html#tmp_file_stream/3)

Create a temporary file and open it

[tnodebug/0](threadutil.html#tnodebug/0)

Switch off debug mode in all threads

[tnodebug/1](threadutil.html#tnodebug/1)

Switch off debug mode in a thread

[tnot/1](tabling-preds.html#tnot/1)

Tabled negation

[told/0](IO.html#told/0)

Close current output

[tprofile/1](threadutil.html#tprofile/1)

Profile a thread for some period

[trace/0](debugger.html#trace/0)

Start the tracer

[tracing/0](debugger.html#tracing/0)

Query status of the tracer

[transaction/1](db.html#transaction/1)

Run goal in a transaction

[transaction/2](db.html#transaction/2)

Run goal in a transaction

[transaction/3](db.html#transaction/3)

Run goal in a transaction

[transaction_updates/1](db.html#transaction_updates/1)

Updates to be committed in a transaction

[trie_delete/3](db.html#trie_delete/3)

Remove term from trie

[trie_destroy/1](db.html#trie_destroy/1)

Destroy a trie

[trie_gen/3](db.html#trie_gen/3)

Get all terms from a trie

[trie_gen_compiled/2](db.html#trie_gen_compiled/2)

Get all terms from a trie

[trie_gen_compiled/3](db.html#trie_gen_compiled/3)

Get all terms from a trie

[trie_insert/2](db.html#trie_insert/2)

Insert term into a trie

[trie_insert/3](db.html#trie_insert/3)

Insert term into a trie

[trie_insert/4](db.html#trie_insert/4)

Insert term into a trie

[trie_lookup/3](db.html#trie_lookup/3)

Lookup a term in a trie

[trie_new/1](db.html#trie_new/1)

Create a trie

[trie_property/2](db.html#trie_property/2)

Examine a trie's properties

[trie_update/3](db.html#trie_update/3)

Update associated value in trie

[trie_term/2](db.html#trie_term/2)

Get term from a trie by handle

[trim_heap/0](memory.html#trim_heap/0)

Release unused **malloc()** resources

[trim_stacks/0](memory.html#trim_stacks/0)

Release unused stack resources

tripwire/2

*(hook)* Handle a tabling tripwire event

[true/0](control.html#true/0)

Succeed

[tspy/1](threadutil.html#tspy/1)

Set spy point and enable debugging in all threads

[tspy/2](threadutil.html#tspy/2)

Set spy point and enable debugging in a thread

[tty_get_capability/3](tty.html#tty_get_capability/3)

Get terminal parameter

[tty_goto/2](tty.html#tty_goto/2)

Goto position on screen

[tty_put/2](tty.html#tty_put/2)

Write control string to terminal

[tty_size/2](tty.html#tty_size/2)

Get row/column size of the terminal

[ttyflush/0](chario.html#ttyflush/0)

Flush output on terminal

[undefined/0](WFS.html#undefined/0)

Well Founded Semantics: true nor false

[undo/1](metacall.html#undo/1)

Schedule goal for backtracking

[unify_with_occurs_check/2](compare.html#unify_with_occurs_check/2)

Logically sound unification

[unifiable/3](compare.html#unifiable/3)

Determining binding required for unification

[unknown/2](debugger.html#unknown/2)

Trap undefined predicates

[unload_file/1](consulting.html#unload_file/1)

Unload a source file

[unload_foreign_library/1](foreignlink.html#unload_foreign_library/1)

`library(shlib)` Detach shared library (.so file)

[unload_foreign_library/2](foreignlink.html#unload_foreign_library/2)

`library(shlib)` Detach shared library (.so file)

[unsetenv/1](system.html#unsetenv/1)

Delete shell environment variable

[untable/1](tabling-preds.html#untable/1)

Remove tabling instrumentation

[upcase_atom/2](chartype.html#upcase_atom/2)

Convert atom to upper-case

[use_foreign_library/1](foreignlink.html#use_foreign_library/1)

Load DLL/shared object (directive)

[use_foreign_library/2](foreignlink.html#use_foreign_library/2)

Load DLL/shared object (directive)

[use_module/1](import.html#use_module/1)

Import a module

[use_module/2](import.html#use_module/2)

Import predicates from a module

[valid_string_goal/1](check.html#valid_string_goal/1)

*(hook)* Goal handles strings

[var/1](typetest.html#var/1)

Type check for unbound variable

[var_number/2](manipterm.html#var_number/2)

Check that var is numbered by numbervars

[var_property/2](consulting.html#var_property/2)

Variable properties during macro expansion

[variant_sha1/2](db.html#variant_sha1/2)

Term-hash for term-variants

[variant_hash/2](db.html#variant_hash/2)

Term-hash for term-variants

[version/0](printmsg.html#version/0)

Print system banner message

[version/1](printmsg.html#version/1)

Add messages to the system banner

[visible/1](debugger.html#visible/1)

Ports that are visible in the tracer

[volatile/1](saved-states.html#volatile/1)

Predicates that are not saved

[wait_for_input/3](streamstat.html#wait_for_input/3)

Wait for input with optional timeout

[when/2](coroutining.html#when/2)

Execute goal when condition becomes true

[wildcard_match/2](miscpreds.html#wildcard_match/2)

POSIX style glob pattern matching

[wildcard_match/3](miscpreds.html#wildcard_match/3)

POSIX style glob pattern matching

[win_add_dll_directory/1](system.html#win_add_dll_directory/1)

Add directory to DLL search path

[win_add_dll_directory/2](system.html#win_add_dll_directory/2)

Add directory to DLL search path

[win_remove_dll_directory/1](system.html#win_remove_dll_directory/1)

Remove directory from DLL search path

[win_exec/2](system.html#win_exec/2)

Win32: spawn Windows task

win_has_menu/0

Win32: true if console menu is available

[win_folder/2](system.html#win_folder/2)

Win32: get special folder by CSIDL

[win_process_modules/1](system.html#win_process_modules/1)

Win32 get .exe and .dll files of the process

[win_shell/2](system.html#win_shell/2)

Win32: open document through Shell

[win_shell/3](system.html#win_shell/3)

Win32: open document through Shell

[win_registry_get_value/3](system.html#win_registry_get_value/3)

Win32: get registry value

[win_get_user_preferred_ui_languages/2](system.html#win_get_user_preferred_ui_languages/2)

Win32: get language preferences

[with_mutex/2](threadsync.html#with_mutex/2)

Run goal while holding mutex

[with_output_to/2](IO.html#with_output_to/2)

Write to strings and more

[with_quasi_quotation_input/3](quasiquotations.html#with_quasi_quotation_input/3)

Parse quasi quotation from stream

[with_tty_raw/1](chario.html#with_tty_raw/1)

Run goal with terminal in raw mode

[working_directory/2](files.html#working_directory/2)

Query/change CWD

[write/1](termrw.html#write/1)

Write term

[write/2](termrw.html#write/2)

Write term to stream

[writeln/1](termrw.html#writeln/1)

Write term, followed by a newline

[writeln/2](termrw.html#writeln/2)

Write term, followed by a newline to a stream

[write_canonical/1](termrw.html#write_canonical/1)

Write a term with quotes, ignore operators

[write_canonical/2](termrw.html#write_canonical/2)

Write a term with quotes, ignore operators on a stream

[write_length/3](termrw.html#write_length/3)

Determine \#characters to output a term

[write_term/2](termrw.html#write_term/2)

Write term with options

[write_term/3](termrw.html#write_term/3)

Write term with options to stream

[writeq/1](termrw.html#writeq/1)

Write term, insert quotes

[writeq/2](termrw.html#writeq/2)

Write term, insert quotes on stream
