
## F.2 Library predicates

### F.2.1 library(aggregate)

|  |  |
|----|----|
| [aggregate/3](aggregate.html#aggregate/3) | Aggregate bindings in Goal according to Template. |
| [aggregate/4](aggregate.html#aggregate/4) | Aggregate bindings in Goal according to Template. |
| [aggregate_all/3](aggregate.html#aggregate_all/3) | Aggregate bindings in Goal according to Template. |
| [aggregate_all/4](aggregate.html#aggregate_all/4) | Aggregate bindings in Goal according to Template. |
| [foldall/4](aggregate.html#foldall/4) | Use Folder to fold V0 to V using all answers of Goal. |
| [foreach/2](aggregate.html#foreach/2) | True when the conjunction of \_instances\_ of Goal created from solutions for Generator is true. |
| [free_variables/4](aggregate.html#free_variables/4) | Find free variables in bagof/setof template. |

### F.2.2 library(ansi_term)

|  |  |
|----|----|
| [ansi_format/3](ansiterm.html#ansi_format/3) | Format text with ANSI attributes. |
| [ansi_format/4](ansiterm.html#ansi_format/4) | Format text with ANSI attributes. |
| [ansi_get_color/2](ansiterm.html#ansi_get_color/2) | Obtain the RGB color for an ANSI color parameter. |
| [ansi_hyperlink/2](ansiterm.html#ansi_hyperlink/2) | Create a hyperlink for a terminal emulator. |
| [ansi_hyperlink/3](ansiterm.html#ansi_hyperlink/3) | Create a hyperlink for a terminal emulator. |
| console_color/2 | Hook that allows for mapping abstract terms to concrete ANSI attributes. |
| [tty_url_hook/2](ansiterm.html#tty_url_hook/2) | Hook for location_url/2. |

### F.2.3 library(apply)

|  |  |
|----|----|
| [convlist/3](apply.html#convlist/3) | Similar to maplist/3, but elements for which call(Goal, ElemIn, \_) fails are omitted from ListOut. |
| [exclude/3](apply.html#exclude/3) | Filter elements for which Goal fails. |
| [foldl/4](apply.html#foldl/4) | Fold an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("fold-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [foldl/5](apply.html#foldl/5) | Fold an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("fold-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [foldl/6](apply.html#foldl/6) | Fold an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("fold-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [foldl/7](apply.html#foldl/7) | Fold an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("fold-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [include/3](apply.html#include/3) | Filter elements for which Goal succeeds. |
| [maplist/2](apply.html#maplist/2) | True if Goal is successfully applied on all matching elements of the list. |
| [maplist/3](apply.html#maplist/3) | True if Goal is successfully applied on all matching elements of the list. |
| [maplist/4](apply.html#maplist/4) | True if Goal is successfully applied on all matching elements of the list. |
| [maplist/5](apply.html#maplist/5) | True if Goal is successfully applied on all matching elements of the list. |
| [partition/4](apply.html#partition/4) | Filter elements of List according to Pred. |
| [partition/5](apply.html#partition/5) | Filter List according to Pred in three sets. |
| [scanl/4](apply.html#scanl/4) | Scan an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("scan-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [scanl/5](apply.html#scanl/5) | Scan an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("scan-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [scanl/6](apply.html#scanl/6) | Scan an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("scan-left"), using columns of \_m\_ list elements as arguments for Goal. |
| [scanl/7](apply.html#scanl/7) | Scan an ensemble of \_m\_ (0 `<=` \_m\_ `<=` 4) lists of length \_n\_ head-to-tail ("scan-left"), using columns of \_m\_ list elements as arguments for Goal. |

### F.2.4 library(assoc)

|  |  |
|----|----|
| [assoc_to_list/2](assoc.html#assoc_to_list/2) | Translate assoc into a pairs list |
| [assoc_to_keys/2](assoc.html#assoc_to_keys/2) | Translate assoc into a key list |
| [assoc_to_values/2](assoc.html#assoc_to_values/2) | Translate assoc into a value list |
| [empty_assoc/1](assoc.html#empty_assoc/1) | Test/create an empty assoc |
| [gen_assoc/3](assoc.html#gen_assoc/3) | Non-deterministic enumeration of assoc |
| [get_assoc/3](assoc.html#get_assoc/3) | Get associated value |
| [get_assoc/5](assoc.html#get_assoc/5) | Get and replace associated value |
| [list_to_assoc/2](assoc.html#list_to_assoc/2) | Translate pair list to assoc |
| [map_assoc/2](assoc.html#map_assoc/2) | Test assoc values |
| [map_assoc/3](assoc.html#map_assoc/3) | Map assoc values |
| [max_assoc/3](assoc.html#max_assoc/3) | Max key-value of an assoc |
| [min_assoc/3](assoc.html#min_assoc/3) | Min key-value of an assoc |
| [ord_list_to_assoc/2](assoc.html#ord_list_to_assoc/2) | Translate ordered list into an assoc |
| [put_assoc/4](assoc.html#put_assoc/4) | Add association to an assoc |

### F.2.5 library(broadcast)

|  |  |
|----|----|
| [broadcast/1](broadcast.html#broadcast/1) | Send event notification |
| [broadcast_request/1](broadcast.html#broadcast_request/1) | Request all agents |
| [listen/2](broadcast.html#listen/2) | Listen to event notifications |
| [listen/3](broadcast.html#listen/3) | Listen to event notifications |
| [unlisten/1](broadcast.html#unlisten/1) | Stop listening to event notifications |
| [unlisten/2](broadcast.html#unlisten/2) | Stop listening to event notifications |
| [unlisten/3](broadcast.html#unlisten/3) | Stop listening to event notifications |
| [listening/3](broadcast.html#listening/3) | Who is listening to event notifications? |

### F.2.6 library(charsio)

|  |  |
|----|----|
| [atom_to_chars/2](charsio.html#atom_to_chars/2) | Convert Atom into a list of character codes. |
| [atom_to_chars/3](charsio.html#atom_to_chars/3) | Convert Atom into a difference list of character codes. |
| [format_to_chars/3](charsio.html#format_to_chars/3) | Use format/2 to write to a list of character codes. |
| [format_to_chars/4](charsio.html#format_to_chars/4) | Use format/2 to write to a difference list of character codes. |
| [number_to_chars/2](charsio.html#number_to_chars/2) | Convert Atom into a list of character codes. |
| [number_to_chars/3](charsio.html#number_to_chars/3) | Convert Number into a difference list of character codes. |
| [open_chars_stream/2](charsio.html#open_chars_stream/2) | Open Codes as an input stream. |
| [read_from_chars/2](charsio.html#read_from_chars/2) | Read Codes into Term. |
| [read_term_from_chars/3](charsio.html#read_term_from_chars/3) | Read Codes into Term. |
| [with_output_to_chars/2](charsio.html#with_output_to_chars/2) | Run Goal as with once/1. |
| [with_output_to_chars/3](charsio.html#with_output_to_chars/3) | Run Goal as with once/1. |
| [with_output_to_chars/4](charsio.html#with_output_to_chars/4) | Same as with_output_to_chars/3 using an explicit stream. |
| [write_to_chars/2](charsio.html#write_to_chars/2) | Write a term to a code list. |
| [write_to_chars/3](charsio.html#write_to_chars/3) | Write a term to a code list. |

### F.2.7 library(check)

|  |  |
|----|----|
| [check/0](check.html#check/0) | Run all consistency checks defined by checker/2. |
| [checker/2](check.html#checker/2) | Register code validation routines. |
| [list_autoload/0](check.html#list_autoload/0) | Report predicates that may be auto-loaded. |
| [list_cross_module_calls/0](check.html#list_cross_module_calls/0) | List calls from one module to another using Module:Goal where the callee is not defined exported, public or multifile, i.e., where the callee should be considered \_private\_. |
| [list_format_errors/0](check.html#list_format_errors/0) | List argument errors for format/2,3. |
| [list_format_errors/1](check.html#list_format_errors/1) | List argument errors for format/2,3. |
| [list_rationals/0](check.html#list_rationals/0) | List rational numbers that appear in clauses. |
| [list_rationals/1](check.html#list_rationals/1) | List rational numbers that appear in clauses. |
| [list_redefined/0](check.html#list_redefined/0) | Lists predicates that are defined in the global module =user= as well as in a normal module; that is, predicates for which the local definition overrules the global default definition. |
| [list_strings/0](string.html#list_strings/0) | List strings that appear in clauses. |
| [list_strings/1](check.html#list_strings/1) | List strings that appear in clauses. |
| [list_trivial_fails/0](check.html#list_trivial_fails/0) | List goals that trivially fail because there is no matching clause. |
| [list_trivial_fails/1](check.html#list_trivial_fails/1) | List goals that trivially fail because there is no matching clause. |
| [list_undefined/0](check.html#list_undefined/0) | Report undefined predicates. |
| [list_undefined/1](check.html#list_undefined/1) | Report undefined predicates. |
| [list_void_declarations/0](check.html#list_void_declarations/0) | List predicates that have declared attributes, but no clauses. |
| [string_predicate/1](check.html#string_predicate/1) | Multifile hook to disable list_strings/0 on the given predicate. |
| [trivial_fail_goal/1](check.html#trivial_fail_goal/1) | Multifile hook that tells list_trivial_fails/0 to accept Goal as valid. |
| [valid_string_goal/1](check.html#valid_string_goal/1) | Multifile hook that qualifies Goal as valid for list_strings/0. |

### F.2.8 library(clpb)

|  |  |
|----|----|
| [labeling/1](clpb.html#labeling/1) | Enumerate concrete solutions. |
| [random_labeling/2](clpb.html#random_labeling/2) | Select a single random solution. |
| [sat/1](clpb.html#sat/1) | True iff Expr is a satisfiable Boolean expression. |
| [sat_count/2](clpb.html#sat_count/2) | Count the number of admissible assignments. |
| [taut/2](clpb.html#taut/2) | Tautology check. |
| [weighted_maximum/3](clpb.html#weighted_maximum/3) | Enumerate weighted optima over admissible assignments. |

### F.2.9 library(clpfd)

|  |  |
|----|----|
| [\#/\\2](clpfd.html##/\/2) | P and Q hold. |
| [\#\</2](clpfd.html##%3C/2) | The arithmetic expression X is less than Y. |
| [\#\<==/2](clpfd.html##%3C==/2) | Q implies P. |
| [\#\<==\>/2](clpfd.html##%3C==%3E/2) | P and Q are equivalent. |
| [\#=/2](clpfd.html##=/2) | The arithmetic expression X equals Y. |
| [\#=\</2](clpfd.html##=%3C/2) | The arithmetic expression X is less than or equal to Y. |
| [\#==\>/2](clpfd.html##==%3E/2) | P implies Q. |
| [\#\>/2](clpfd.html##%3E/2) | Same as Y `#<` X. |
| [\#\>=/2](clpfd.html##%3E=/2) | Same as Y `#=<` X. |
| [\#\\1](clpfd.html##\/1) | Q does \_not\_ hold. |
| [\#\\2](clpfd.html##\/2) | Either P holds or Q holds, but not both. |
| [\#\\/2](clpfd.html##\//2) | P or Q holds. |
| [\#\\/2](clpfd.html##\=/2) | The arithmetic expressions X and Y evaluate to distinct integers. |
| [all_different/1](clpfd.html#all_different/1) | Like all_distinct/1, but with weaker propagation. |
| [all_distinct/1](clpfd.html#all_distinct/1) | True iff Vars are pairwise distinct. |
| [automaton/3](clpfd.html#automaton/3) | Describes a list of finite domain variables with a finite automaton. |
| [automaton/8](clpfd.html#automaton/8) | Describes a list of finite domain variables with a finite automaton. |
| [chain/2](clpfd.html#chain/2) | Zs form a chain with respect to Relation. |
| [circuit/1](clpfd.html#circuit/1) | True iff the list Vs of finite domain variables induces a Hamiltonian circuit. |
| [cumulative/1](clpfd.html#cumulative/1) | Equivalent to cumulative(Tasks, \[limit(1)\]). |
| [cumulative/2](clpfd.html#cumulative/2) | Schedule with a limited resource. |
| [disjoint2/1](clpfd.html#disjoint2/1) | True iff Rectangles are not overlapping. |
| [element/3](clpfd.html#element/3) | The N-th element of the list of finite domain variables Vs is V. |
| [empty_fdset/1](clpfd.html#empty_fdset/1) | Set is the empty FD set. |
| [empty_interval/2](clpfd.html#empty_interval/2) | Min..Max is an empty interval. |
| [fd_degree/2](clpfd.html#fd_degree/2) | Degree is the number of constraints currently attached to Var. |
| [fd_dom/2](clpfd.html#fd_dom/2) | Dom is the current domain (see in/2) of Var. |
| [fd_inf/2](clpfd.html#fd_inf/2) | Inf is the infimum of the current domain of Var. |
| [fd_set/2](clpfd.html#fd_set/2) | Set is the FD set representation of the current domain of Var. |
| [fd_size/2](clpfd.html#fd_size/2) | Reflect the current size of a domain. |
| [fd_sup/2](clpfd.html#fd_sup/2) | Sup is the supremum of the current domain of Var. |
| [fd_var/1](clpfd.html#fd_var/1) | True iff Var is a CLP(FD) variable. |
| [fdset_add_element/3](clpfd.html#fdset_add_element/3) | Set2 is the same FD set as Set1, but with the integer Elt added. |
| [fdset_complement/2](clpfd.html#fdset_complement/2) | The FD set Complement is the complement of the FD set Set. |
| [fdset_del_element/3](clpfd.html#fdset_del_element/3) | Set2 is the same FD set as Set1, but with the integer Elt removed. |
| [fdset_disjoint/2](clpfd.html#fdset_disjoint/2) | The FD sets Set1 and Set2 have no elements in common. |
| [fdset_eq/2](clpfd.html#fdset_eq/2) | True if the FD sets Set1 and Set2 are equal, i. |
| [fdset_intersect/2](clpfd.html#fdset_intersect/2) | The FD sets Set1 and Set2 have at least one element in common. |
| [fdset_intersection/3](clpfd.html#fdset_intersection/3) | Intersection is an FD set (possibly empty) of all elements that the FD sets Set1 and Set2 have in common. |
| [fdset_interval/3](clpfd.html#fdset_interval/3) | Interval is a non-empty FD set consisting of the single interval Min..Max. |
| [fdset_max/2](clpfd.html#fdset_max/2) | Max is the upper bound (supremum) of the non-empty FD set Set. |
| [fdset_member/2](clpfd.html#fdset_member/2) | The integer Elt is a member of the FD set Set. |
| [fdset_min/2](clpfd.html#fdset_min/2) | Min is the lower bound (infimum) of the non-empty FD set Set. |
| [fdset_parts/4](clpfd.html#fdset_parts/4) | Set is a non-empty FD set representing the domain Min..Max `\/` Rest, where Min..Max is a non-empty interval (see fdset_interval/3) and Rest is another FD set (possibly empty). |
| [fdset_singleton/2](clpfd.html#fdset_singleton/2) | Set is the FD set containing the single integer Elt. |
| [fdset_size/2](clpfd.html#fdset_size/2) | Size is the number of elements of the FD set Set, or the atom \*sup\* if Set is infinite. |
| [fdset_subset/2](clpfd.html#fdset_subset/2) | The FD set Set1 is a (non-strict) subset of Set2, i. |
| [fdset_subtract/3](clpfd.html#fdset_subtract/3) | The FD set Difference is Set1 with all elements of Set2 removed, i. |
| [fdset_to_list/2](clpfd.html#fdset_to_list/2) | List is a list containing all elements of the finite FD set Set, in ascending order. |
| [fdset_to_range/2](clpfd.html#fdset_to_range/2) | Domain is a domain equivalent to the FD set Set. |
| [fdset_union/2](clpfd.html#fdset_union/2) | The FD set Union is the n-ary union of all FD sets in the list Sets. |
| [fdset_union/3](clpfd.html#fdset_union/3) | The FD set Union is the union of FD sets Set1 and Set2. |
| [global_cardinality/2](clpfd.html#global_cardinality/2) | Global Cardinality constraint. |
| [global_cardinality/3](clpfd.html#global_cardinality/3) | Global Cardinality constraint. |
| [in/2](clpfd.html#in/2) | Var is an element of Domain. |
| [in_set/2](clpfd.html#in_set/2) | Var is an element of the FD set Set. |
| [indomain/1](clpfd.html#indomain/1) | Bind Var to all feasible values of its domain on backtracking. |
| [ins/2](clpfd.html#ins/2) | The variables in the list Vars are elements of Domain. |
| [is_fdset/1](clpfd.html#is_fdset/1) | Set is currently bound to a valid FD set. |
| [label/1](clpfd.html#label/1) | Equivalent to labeling(\[\], Vars). |
| [labeling/2](clpfd.html#labeling/2) | Assign a value to each variable in Vars. |
| [lex_chain/1](clpfd.html#lex_chain/1) | Lists are lexicographically non-decreasing. |
| [list_to_fdset/2](clpfd.html#list_to_fdset/2) | Set is an FD set containing all elements of List, which must be a list of integers. |
| [range_to_fdset/2](clpfd.html#range_to_fdset/2) | Set is an FD set equivalent to the domain Domain. |
| [scalar_product/4](clpfd.html#scalar_product/4) | True iff the scalar product of Cs and Vs is in relation Rel to Expr. |
| [serialized/2](clpfd.html#serialized/2) | Describes a set of non-overlapping tasks. |
| [sum/3](clpfd.html#sum/3) | The sum of elements of the list Vars is in relation Rel to Expr. |
| [transpose/2](clpfd.html#transpose/2) | Transpose a list of lists of the same length. |
| [tuples_in/2](clpfd.html#tuples_in/2) | True iff all Tuples are elements of Relation. |
| [zcompare/3](clpfd.html#zcompare/3) | Analogous to compare/3, with finite domain variables A and B. |

### F.2.10 library(clpqr)

|  |  |
|----|----|
| [entailed/1](clpqr.html#entailed/1) | Check if constraint is entailed |
| [inf/2](clpqr.html#inf/2) | Find the infimum of an expression |
| [sup/2](clpqr.html#sup/2) | Find the supremum of an expression |
| [minimize/1](clpqr.html#minimize/1) | Minimizes an expression |
| [maximize/1](clpqr.html#maximize/1) | Maximizes an expression |
| [bb_inf/3](clpqr.html#bb_inf/3) | Infimum of expression for mixed-integer problems |
| [bb_inf/4](clpqr.html#bb_inf/4) | Infimum of expression for mixed-integer problems |
| [bb_inf/5](clpqr.html#bb_inf/5) | Infimum of expression for mixed-integer problems |
| [dump/3](clpqr.html#dump/3) | Dump constraints on variables |

### F.2.11 library(csv)

|  |  |
|----|----|
| [csv_options/2](csv.html#csv_options/2) | Compiled is the compiled representation of the CSV processing options as they may be passed into csv`//`2, etc. |
| [csv_read_file/2](csv.html#csv_read_file/2) | Read a CSV file into a list of rows. |
| [csv_read_file/3](csv.html#csv_read_file/3) | Read a CSV file into a list of rows. |
| [csv_read_file_row/3](csv.html#csv_read_file_row/3) | True when Row is a row in File. |
| [csv_read_row/3](csv.html#csv_read_row/3) | Read the next CSV record from Stream and unify the result with Row. |
| [csv_read_stream/3](csv.html#csv_read_stream/3) | Read CSV data from Stream. |
| [csv_write_file/2](csv.html#csv_write_file/2) | Write a list of Prolog terms to a CSV file. |
| [csv_write_file/3](csv.html#csv_write_file/3) | Write a list of Prolog terms to a CSV file. |
| [csv_write_stream/3](csv.html#csv_write_stream/3) | Write the rows in Data to Stream. |
| [csv//1](csv.html#csv//1) | Prolog DCG to‘read/write’CSV data. |
| [csv//2](csv.html#csv//2) | Prolog DCG to‘read/write’CSV data. |

### F.2.12 library(dcgbasics)

|  |  |
|----|----|
| [alpha_to_lower//1](basics.html#alpha_to_lower//1) | Read a letter (class =alpha=) and return it as a lowercase letter. |
| [atom//1](basics.html#atom//1) | Generate codes of Atom. |
| [blank//0](basics.html#blank//0) | Take next =space= character from input. |
| [blanks//0](basics.html#blanks//0) | Skip zero or more white-space characters. |
| [blanks_to_nl//0](basics.html#blanks_to_nl//0) | Take a sequence of blank`//`0 codes if blanks are followed by a newline or end of the input. |
| [csym//1](basics.html#csym//1) | Recognise a C symbol according to the‘csymf\` and‘csym\` code type classification provided by the C library. |
| [digit//1](basics.html#digit//1) | Number processing. |
| [digits//1](basics.html#digits//1) | Number processing. |
| [eol//0](basics.html#eol//0) | Matches end-of-line. |
| [eos//0](basics.html#eos//0) | Matches end-of-input. |
| [float//1](basics.html#float//1) | Process a floating point number. |
| [integer//1](basics.html#integer//1) | Number processing. |
| [nonblank//1](basics.html#nonblank//1) | Code is the next non-blank (=graph=) character. |
| [nonblanks//1](basics.html#nonblanks//1) | Take all =graph= characters. |
| [number//1](basics.html#number//1) | Generate extract a number. |
| [prolog_var_name//1](basics.html#prolog_var_name//1) | Matches a Prolog variable name. |
| [remainder//1](basics.html#remainder//1) | Unify List with the remainder of the input. |
| [string//1](basics.html#string//1) | Take as few as possible tokens from the input, taking one more each time on backtracking. |
| [string_without//2](basics.html#string_without//2) | Take as many codes from the input until the next character code appears in the list EndCodes. |
| [white//0](basics.html#white//0) | Take next =white= character from input. |
| [whites//0](basics.html#whites//0) | Skip white space \_inside\_ a line. |
| [xdigit//1](basics.html#xdigit//1) | True if the next code is a hexdecimal digit with Weight. |
| [xdigits//1](basics.html#xdigits//1) | List of weights of a sequence of hexadecimal codes. |
| [xinteger//1](basics.html#xinteger//1) | Generate or extract an integer from a sequence of hexadecimal digits. |

### F.2.13 library(dcghighorder)

|  |  |
|----|----|
| [foreach//2](highorder.html#foreach//2) | Generate a list from the solutions of Generator. |
| [foreach//3](highorder.html#foreach//3) | Generate a list from the solutions of Generator. |
| [optional//2](highorder.html#optional//2) | Perform an optional match, executing Default if Match is not matched. |
| [sequence//2](highorder.html#sequence//2) | Match or generate a sequence of Element. |
| [sequence//3](highorder.html#sequence//3) | Match or generate a sequence of Element where each pair of elements is separated by Sep. |
| [sequence//5](highorder.html#sequence//5) | Match or generate a sequence of Element enclosed by Start end End, where each pair of elements is separated by Sep. |

### F.2.14 library(debug)

|  |  |
|----|----|
| [assertion/1](debug.html#assertion/1) | Acts similar to C assert() macro. |
| assertion_failed/2 | This hook is called if the Goal of assertion/1 fails. |
| [debug/1](debug.html#debug/1) | Add/remove a topic from being printed. |
| [debug/3](debug.html#debug/3) | Format a message if debug topic is enabled. |
| [debug_message_context/1](debug.html#debug_message_context/1) | Specify additional context for debug messages. |
| debug_print_hook/3 | Hook called by debug/3. |
| [debugging/1](debug.html#debugging/1) | Examine debug topics. |
| [debugging/2](debug.html#debugging/2) | Examine debug topics. |
| [list_debug_topics/0](debug.html#list_debug_topics/0) | List currently known topics for debug/3 and their setting. |
| [list_debug_topics/1](debug.html#list_debug_topics/1) | List currently known topics for debug/3 and their setting. |
| [nodebug/1](debug.html#nodebug/1) | Add/remove a topic from being printed. |

### F.2.15 library(dicts)

|  |  |
|----|----|
| [dict_fill/4](dicts.html#dict_fill/4) | Implementation for the dicts_to_same_keys/3‘OnEmpty\` closure that fills new cells with a copy of ValueIn. |
| [dict_keys/2](dicts.html#dict_keys/2) | True when Keys is an ordered set of the keys appearing in Dict. |
| [dict_size/2](dicts.html#dict_size/2) | True when KeyCount is the number of keys in Dict. |
| [dicts_join/3](dicts.html#dicts_join/3) | Join dicts in Dicts that have the same value for Key, provided they do not have conflicting values on other keys. |
| [dicts_join/4](dicts.html#dicts_join/4) | Join two lists of dicts (Dicts1 and Dicts2) on Key. |
| [dicts_same_keys/2](dicts.html#dicts_same_keys/2) | True if List is a list of dicts that all have the same keys and Keys is an ordered set of these keys. |
| [dicts_same_tag/2](dicts.html#dicts_same_tag/2) | True when List is a list of dicts that all have the tag Tag. |
| [dicts_slice/3](dicts.html#dicts_slice/3) | DictsOut is a list of Dicts only containing values for Keys. |
| [dicts_to_compounds/4](dicts.html#dicts_to_compounds/4) | True when Dicts and Compounds are lists of the same length and each element of Compounds is a compound term whose arguments represent the values associated with the corresponding keys in Keys. |
| [dicts_to_same_keys/3](dicts.html#dicts_to_same_keys/3) | DictsOut is a copy of DictsIn, where each dict contains all keys appearing in all dicts of DictsIn. |
| [mapdict/2](dicts.html#mapdict/2) | True when all dicts have the same set of keys and call(Goal, Key, V1, ...) is true for all keys in the dicts. |
| [mapdict/3](dicts.html#mapdict/3) | True when all dicts have the same set of keys and call(Goal, Key, V1, ...) is true for all keys in the dicts. |
| [mapdict/4](dicts.html#mapdict/4) | True when all dicts have the same set of keys and call(Goal, Key, V1, ...) is true for all keys in the dicts. |

### F.2.16 library(error)

|  |  |
|----|----|
| [current_encoding/1](error.html#current_encoding/1) | True if Name is the name of a supported encoding. |
| [current_type/3](error.html#current_type/3) | True when Type is a currently defined type and Var satisfies Type of the body term Body succeeds. |
| [domain_error/2](error.html#domain_error/2) | The argument is of the proper type, but has a value that is outside the supported values. |
| [existence_error/2](error.html#existence_error/2) | Culprit is of the correct type and correct domain, but there is no existing (external) resource of type ObjectType that is represented by it. |
| [existence_error/3](error.html#existence_error/3) | Culprit is of the correct type and correct domain, but there is no existing (external) resource of type ObjectType that is represented by it in the provided set. |
| [has_type/2](error.html#has_type/2) | True if Term satisfies Type. |
| [instantiation_error/1](error.html#instantiation_error/1) | An argument is under-instantiated. |
| [is_of_type/2](error.html#is_of_type/2) | True if Term satisfies Type. |
| [must_be/2](error.html#must_be/2) | True if Term satisfies the type constraints for Type. |
| [permission_error/3](error.html#permission_error/3) | It is not allowed to perform Operation on (whatever is represented by) Culprit that is of the given PermissionType (in fact, the ISO Standard is confusing and vague about these terms’meaning). |
| [representation_error/1](error.html#representation_error/1) | A representation error indicates a limitation of the implementation. |
| [resource_error/1](error.html#resource_error/1) | A goal cannot be completed due to lack of resources. |
| [syntax_error/1](error.html#syntax_error/1) | A text has invalid syntax. |
| [type_error/2](error.html#type_error/2) | Tell the user that Culprit is not of the expected ValidType. |
| [uninstantiation_error/1](error.html#uninstantiation_error/1) | An argument is over-instantiated. |

### F.2.17 library(exceptions)

|  |  |
|----|----|
| [catch/4](exceptions.html#catch/4) | As catch/3, only catching exceptions for which exception(ErrorType,Ball) is true. |
| [error_term/2](exceptions.html#error_term/2) | Describe the formal part of error(Formal,ImplDefined) exceptions. |
| [exception/2](exceptions.html#exception/2) | If Ball is unbound, adds a delayed goal that tests the error belongs to Type when Ball is instantiated (by catch/3). |
| [exception_term/2](exceptions.html#exception_term/2) | Describe exceptions that are not error(Formal, \_) terms. |
| [exception_type/2](exceptions.html#exception_type/2) | Declare all exceptions subsumed by Term to be an exception of Type. |

### F.2.18 library(fastrw)

|  |  |
|----|----|
| [fast_read/1](fastrw.html#fast_read/1) | The next term is read from current standard input and is unified with Term. |
| [fast_write/1](fastrw.html#fast_write/1) | Output Term in a way that fast_read/1 and fast_read/2 will be able to read it back. |
| [fast_write_to_string/3](fastrw.html#fast_write_to_string/3) | Perform a fast-write to the difference-slist String`\`Tail. |

### F.2.19 library(explain)

|  |  |
|----|----|
| [explain/1](online-help.html#explain/1) | Give an explanation on Term. |
| [explain/2](online-help.html#explain/2) | True when Explanation is an explanation of Term. |

### F.2.20 library(help)

|  |  |
|----|----|
| [apropos/1](online-help.html#apropos/1) | Print objects from the manual whose name or summary match with Query. |
| [help/0](online-help.html#help/0) | Show help for What. |
| [help/1](online-help.html#help/1) | Show help for What. |
| [help_text/2](online-help.html#help_text/2) | When =Predicate= is a term of the form =Name/Arity= for which documentation exists, =HelpText= is the documentation in textual format (parsed from the HTML help). |
| [show_html_hook/1](online-help.html#show_html_hook/1) | Hook called to display the extracted HTML document. |

### F.2.21 library(gensym)

|  |  |
|----|----|
| [gensym/2](gensym.html#gensym/2) | Generate `<`Base`>`1, `<`Base`>`2, etc atoms on each subsequent call. |
| [reset_gensym/0](gensym.html#reset_gensym/0) | Reset gensym for all registered keys. |
| [reset_gensym/1](gensym.html#reset_gensym/1) | Restart generation of identifiers from Base at `<`Base`>`1. |

### F.2.22 library(heaps)

|  |  |
|----|----|
| [add_to_heap/4](heaps.html#add_to_heap/4) | Adds Value with priority Key to Heap0, constructing a new heap in Heap. |
| [delete_from_heap/4](heaps.html#delete_from_heap/4) | Deletes Value from Heap0, leaving its priority in Key and the resulting data structure in Heap. |
| [empty_heap/1](heaps.html#empty_heap/1) | True if Heap is an empty heap. |
| [get_from_heap/4](heaps.html#get_from_heap/4) | Retrieves the minimum-priority pair Key-Value from Heap0. |
| [heap_size/2](heaps.html#heap_size/2) | Determines the number of elements in Heap. |
| [heap_to_list/2](heaps.html#heap_to_list/2) | Constructs a list List of Key-Value terms, ordered by (ascending) priority. |
| [is_heap/1](heaps.html#is_heap/1) | Returns true if X is a heap. |
| [list_to_heap/2](heaps.html#list_to_heap/2) | If List is a list of Key-Value terms, constructs a heap out of List. |
| [merge_heaps/3](heaps.html#merge_heaps/3) | Merge the two heaps Heap0 and Heap1 in Heap. |
| [min_of_heap/3](heaps.html#min_of_heap/3) | Unifies Value with the minimum-priority element of Heap and Key with its priority value. |
| [min_of_heap/5](heaps.html#min_of_heap/5) | Gets the two minimum-priority elements from Heap. |
| [singleton_heap/3](heaps.html#singleton_heap/3) | True if Heap is a heap with the single element Key-Value. |

### F.2.23 library(increval)

|  |  |
|----|----|
| [incr_directly_depends/2](increval.html#incr_directly_depends/2) | True if Goal1 depends on Goal2 in the IDG. |
| [incr_invalid_subgoals/1](increval.html#incr_invalid_subgoals/1) | List is a sorted list (set) of the incremental subgoals that are currently invalid. |
| [incr_invalidate_call/1](increval.html#incr_invalidate_call/1) | This is the XSB name, but the manual says incr_invalidate_calls/1 and the comment with the code suggests this is misnamed. |
| [incr_invalidate_calls/1](increval.html#incr_invalidate_calls/1) | Invalidate all tables for subgoals of Goal as well as tables that are affected by these. |
| [incr_is_invalid/1](increval.html#incr_is_invalid/1) | True when Subgoal's table is marked as invalid. |
| [incr_propagate_calls/1](tabling-monotonic.html#incr_propagate_calls/1) | Activate the monotonic answer propagation similarly to when a new fact is asserted for a monotonic dynamic predicate. |
| [incr_table_update/0](increval.html#incr_table_update/0) | Updated all invalid tables. |
| [incr_trans_depends/2](increval.html#incr_trans_depends/2) | True for each pair in the transitive closure of incr_directly_depends(G1, G2). |
| [is_incremental_subgoal/1](increval.html#is_incremental_subgoal/1) | This predicate non-deterministically unifies Subgoal with incrementally tabled subgoals that are currently table entries. |

### F.2.24 library(intercept)

|  |  |
|----|----|
| [intercept/3](intercept.html#intercept/3) | Run Goal as call/1. |
| [intercept/4](intercept.html#intercept/4) | Similar to intercept/3, but the copy of Handler is called as call(Copy,Arg), which allows passing large context arguments or arguments subject to unification or \_destructive assignment\_. |
| [intercept_all/4](intercept.html#intercept_all/4) | True when List contains all instances of Template that have been sent using send_signal/1 where the argument unifies with Ball. |
| [nb_intercept_all/4](intercept.html#nb_intercept_all/4) | As intercept_all/4, but backtracing inside Goal does not reset List. |
| [send_signal/1](intercept.html#send_signal/1) | If this predicate is called from a sub-goal of intercept/3, execute the associated \_Handler\_ of the intercept/3 environment. |
| [send_silent_signal/1](intercept.html#send_silent_signal/1) | As send_signal/1, but succeed silently if there is no matching intercept environment. |

### F.2.25 library(iostream)

|  |  |
|----|----|
| [close_any/1](iostream.html#close_any/1) | Execute the‘Close\` closure returned by open_any/5. |
| [open_any/5](iostream.html#open_any/5) | Establish a stream from Specification that should be closed using Close, which can either be called or passed to close_any/1. |
| [open_hook/6](iostream.html#open_hook/6) | Open Spec in Mode, producing Stream. |

### F.2.26 library(listing)

|  |  |
|----|----|
| [listing/0](listing.html#listing/0) | Lists all predicates defined in the calling module. |
| [listing/1](listing.html#listing/1) | List matching clauses. |
| [listing/2](listing.html#listing/2) | List matching clauses. |
| [portray_clause/1](listing.html#portray_clause/1) | Portray‘Clause’on the current output stream. |
| [portray_clause/2](listing.html#portray_clause/2) | Portray‘Clause’on the current output stream. |
| [portray_clause/3](listing.html#portray_clause/3) | Portray‘Clause’on the current output stream. |

### F.2.27 library(lists)

|  |  |
|----|----|
| [append/2](lists.html#append/2) | Concatenate a list of lists. |
| [append/3](lists.html#append/3) | List1AndList2 is the concatenation of List1 and List2. |
| [clumped/2](lists.html#clumped/2) | Pairs is a list of‘Item-Count\` pairs that represents the \_run length encoding\_ of Items. |
| [delete/3](lists.html#delete/3) | Delete matching elements from a list. |
| [flatten/2](lists.html#flatten/2) | Is true if FlatList is a non-nested version of NestedList. |
| [intersection/3](lists.html#intersection/3) | True if Set3 unifies with the intersection of Set1 and Set2. |
| [is_set/1](lists.html#is_set/1) | True if Set is a proper list without duplicates. |
| [last/2](lists.html#last/2) | Succeeds when Last is the last element of List. |
| [list_to_set/2](lists.html#list_to_set/2) | True when Set has the same elements as List in the same order. |
| [max_list/2](lists.html#max_list/2) | True if Max is the largest number in List. |
| [max_member/2](lists.html#max_member/2) | True when Max is the largest member in the standard order of terms. |
| [max_member/3](lists.html#max_member/3) | True when Max is the largest member according to Pred, which must be a 2-argument callable that behaves like (`@=<`)/2. |
| [member/2](lists.html#member/2) | True if Elem is a member of List. |
| [min_list/2](lists.html#min_list/2) | True if Min is the smallest number in List. |
| [min_member/2](lists.html#min_member/2) | True when Min is the smallest member in the standard order of terms. |
| [min_member/3](lists.html#min_member/3) | True when Min is the smallest member according to Pred, which must be a 2-argument callable that behaves like (`@=<`)/2. |
| [nextto/3](lists.html#nextto/3) | True if Y directly follows X in List. |
| [nth0/3](lists.html#nth0/3) | True when Elem is the Index'th element of List. |
| [nth0/4](lists.html#nth0/4) | Select/insert element at index. |
| [nth1/3](lists.html#nth1/3) | Is true when Elem is the Index'th element of List. |
| [nth1/4](lists.html#nth1/4) | As nth0/4, but counting starts at 1. |
| [numlist/3](lists.html#numlist/3) | List is a list \[Low, Low+1, ... High\]. |
| [permutation/2](lists.html#permutation/2) | True when Xs is a permutation of Ys. |
| [prefix/2](lists.html#prefix/2) | True iff Part is a leading substring of Whole. |
| [proper_length/2](lists.html#proper_length/2) | True when Length is the number of elements in the proper list List. |
| [reverse/2](lists.html#reverse/2) | Is true when the elements of List2 are in reverse order compared to List1. |
| [same_length/2](lists.html#same_length/2) | Is true when List1 and List2 are lists with the same number of elements. |
| [select/3](lists.html#select/3) | Is true when List1, with Elem removed, results in List2. |
| [select/4](lists.html#select/4) | Select from two lists at the same position. |
| [selectchk/3](lists.html#selectchk/3) | Semi-deterministic removal of first element in List that unifies with Elem. |
| [selectchk/4](lists.html#selectchk/4) | Semi-deterministic version of select/4. |
| [subseq/3](lists.html#subseq/3) | Is true when SubList contains a subset of the elements of List in the same order and Complement contains all elements of List not in SubList, also in the order they appear in List. |
| [subset/2](lists.html#subset/2) | True if all elements of SubSet belong to Set as well. |
| [subtract/3](lists.html#subtract/3) | Delete all elements in Delete from Set. |
| [sum_list/2](lists.html#sum_list/2) | Sum is the result of adding all numbers in List. |
| [union/3](lists.html#union/3) | True if Set3 unifies with the union of the lists Set1 and Set2. |

### F.2.28 library(macros)

|  |  |
|----|----|
| [expand_macros/5](macros.html#expand_macros/5) | Perform macro expansion on TermIn with layout PosIn to produce TermOut with layout PosOut. |
| [include_macros/3](macros.html#include_macros/3) | Include macros from another module. |
| [macro_position/1](macros.html#macro_position/1) | True when Position is the position of the macro. |

### F.2.29 library(main)

|  |  |
|----|----|
| [argv_options/3](main.html#argv_options/3) | Parse command line arguments. |
| [argv_options/4](main.html#argv_options/4) | As argv_options/3 in \_\_guided\_\_ mode, Currently this version allows parsing argument options throwing an exception rather than calling halt/1 by passing an empty list to ParseOptions. |
| [argv_usage/1](main.html#argv_usage/1) | Use print_message/2 to print a usage message at Level. |
| [cli_debug_opt_help/2](main.html#cli_debug_opt_help/2) | Implements opt_type/3, opt_help/2 and opt_meta/2 for debug arguments. |
| [cli_debug_opt_meta/2](main.html#cli_debug_opt_meta/2) | Implements opt_type/3, opt_help/2 and opt_meta/2 for debug arguments. |
| [cli_debug_opt_type/3](main.html#cli_debug_opt_type/3) | Implements opt_type/3, opt_help/2 and opt_meta/2 for debug arguments. |
| [cli_enable_development_system/0](main.html#cli_enable_development_system/0) | Re-enable the development environment. |
| [cli_parse_debug_options/2](main.html#cli_parse_debug_options/2) | Parse certain commandline options for debugging and development purposes. |
| [main/0](main.html#main/0) | Call main/1 using the passed command-line arguments. |

### F.2.30 library(occurs)

|  |  |
|----|----|
| [contains_term/2](occurs.html#contains_term/2) | Succeeds if Sub is contained in Term (=, deterministically). |
| [contains_var/2](occurs.html#contains_var/2) | Succeeds if Sub is contained in Term (==, deterministically). |
| [free_of_term/2](occurs.html#free_of_term/2) | Succeeds of Sub does not unify to any subterm of Term. |
| [free_of_var/2](occurs.html#free_of_var/2) | Succeeds of Sub is not equal (`==`) to any subterm of Term. |
| [occurrences_of_term/3](occurs.html#occurrences_of_term/3) | Count the number of SubTerms in Term that \_unify\_ with SubTerm. |
| [occurrences_of_var/3](occurs.html#occurrences_of_var/3) | Count the number of SubTerms in Term that are \_equal\_ to SubTerm. |
| [sub_term/2](occurs.html#sub_term/2) | Generates (on backtracking) all subterms of Term. |
| [sub_term_shared_variables/3](occurs.html#sub_term_shared_variables/3) | If Sub is a sub term of Term, Vars is bound to the list of variables in Sub that also appear outside Sub in Term. |
| [sub_var/2](occurs.html#sub_var/2) | Generates (on backtracking) all subterms (`==`) of Term. |

### F.2.31 library(option)

|  |  |
|----|----|
| [dict_options/2](option.html#dict_options/2) | Convert between an option list and a dictionary. |
| [merge_options/3](option.html#merge_options/3) | Merge two option sets. |
| [meta_options/3](option.html#meta_options/3) | Perform meta-expansion on options that are module-sensitive. |
| [option/2](option.html#option/2) | Get an Option from Options. |
| [option/3](option.html#option/3) | Get an Option from Options. |
| [select_option/3](option.html#select_option/3) | Get and remove Option from Options. |
| [select_option/4](option.html#select_option/4) | Get and remove Option with default value. |

### F.2.32 library(optparse)

|  |  |
|----|----|
| [opt_arguments/3](optparse.html#opt_arguments/3) | Extract commandline options according to a specification. |
| [opt_help/2](optparse.html#opt_help/2) | True when Help is a help string synthesized from OptsSpec. |
| [opt_parse/4](optparse.html#opt_parse/4) | Equivalent to opt_parse(OptsSpec, ApplArgs, Opts, PositionalArgs, `[]`). |
| [opt_parse/5](optparse.html#opt_parse/5) | Parse the arguments Args (as list of atoms) according to OptsSpec. |
| [parse_type/3](optparse.html#parse_type/3) | Hook to parse option text Codes to an object of type Type. |

### F.2.33 library(ordsets)

|  |  |
|----|----|
| [is_ordset/1](ordsets.html#is_ordset/1) | True if Term is an ordered set. |
| [list_to_ord_set/2](ordsets.html#list_to_ord_set/2) | Transform a list into an ordered set. |
| [ord_add_element/3](ordsets.html#ord_add_element/3) | Insert an element into the set. |
| [ord_del_element/3](ordsets.html#ord_del_element/3) | Delete an element from an ordered set. |
| [ord_disjoint/2](ordsets.html#ord_disjoint/2) | True if Set1 and Set2 have no common elements. |
| [ord_empty/1](ordsets.html#ord_empty/1) | True when List is the empty ordered set. |
| [ord_intersect/2](ordsets.html#ord_intersect/2) | True if both ordered sets have a non-empty intersection. |
| [ord_intersect/3](ordsets.html#ord_intersect/3) | Intersection holds the common elements of Set1 and Set2. |
| [ord_intersection/2](ordsets.html#ord_intersection/2) | Intersection of a powerset. |
| [ord_intersection/3](ordsets.html#ord_intersection/3) | Intersection holds the common elements of Set1 and Set2. |
| [ord_intersection/4](ordsets.html#ord_intersection/4) | Intersection and difference between two ordered sets. |
| [ord_memberchk/2](ordsets.html#ord_memberchk/2) | True if Element is a member of OrdSet, compared using ==. |
| [ord_selectchk/3](ordsets.html#ord_selectchk/3) | Selectchk/3, specialised for ordered sets. |
| [ord_seteq/2](ordsets.html#ord_seteq/2) | True if Set1 and Set2 have the same elements. |
| [ord_subset/2](ordsets.html#ord_subset/2) | Is true if all elements of Sub are in Super. |
| [ord_subtract/3](ordsets.html#ord_subtract/3) | Diff is the set holding all elements of InOSet that are not in NotInOSet. |
| [ord_symdiff/3](ordsets.html#ord_symdiff/3) | Is true when Difference is the symmetric difference of Set1 and Set2. |
| [ord_union/2](ordsets.html#ord_union/2) | True if Union is the union of all elements in the superset SetOfSets. |
| [ord_union/3](ordsets.html#ord_union/3) | Union is the union of Set1 and Set2. |
| [ord_union/4](ordsets.html#ord_union/4) | True iff ord_union(Set1, Set2, Union) and ord_subtract(Set2, Set1, New). |

### F.2.34 library(persistency)

|  |  |
|----|----|
| [current_persistent_predicate/1](persistency.html#current_persistent_predicate/1) | True if PI is a predicate that provides access to the persistent database DB. |
| [db_assert/1](persistency.html#db_assert/1) | Assert Term into the database and record it for persistency. |
| [db_attach/2](persistency.html#db_attach/2) | Use File as persistent database for the calling module. |
| [db_attached/1](persistency.html#db_attached/1) | True if the context module attached to the persistent database File. |
| [db_detach/0](persistency.html#db_detach/0) | Detach persistency from the calling module and delete all persistent clauses from the Prolog database. |
| [db_retract/1](persistency.html#db_retract/1) | Retract terms from the database one-by-one. |
| [db_retractall/1](persistency.html#db_retractall/1) | Retract all matching facts and do the same in the database. |
| [db_sync/1](persistency.html#db_sync/1) | Synchronise database with the associated file. |
| [db_sync_all/1](persistency.html#db_sync_all/1) | Sync all registered databases. |
| [persistent/1](persistency.html#persistent/1) | Declare dynamic database terms. |

### F.2.35 library(portraytext)

|  |  |
|----|----|
| [is_text_code/1](portraytext.html#is_text_code/1) | Multifile hook that can be used to extend the set of character codes that is recognised as likely text. |
| [portray_text/1](portraytext.html#portray_text/1) | Switch portraying on or off. |
| [set_portray_text/2](portraytext.html#set_portray_text/2) | Set options for portraying. |
| [set_portray_text/3](portraytext.html#set_portray_text/3) | Set options for portraying. |

### F.2.36 library(predicate_options)

|  |  |
|----|----|
| [assert_predicate_options/4](predicate_options.html#assert_predicate_options/4) | As predicate_options(:PI, +Arg, +Options). |
| [check_predicate_option/3](predicate_options.html#check_predicate_option/3) | Verify predicate options at runtime. |
| [check_predicate_options/0](predicate_options.html#check_predicate_options/0) | Analyse loaded program for erroneous options. |
| [current_option_arg/2](predicate_options.html#current_option_arg/2) | True when Arg of PI processes predicate options. |
| [current_predicate_option/3](predicate_options.html#current_predicate_option/3) | True when Arg of PI processes Option. |
| [current_predicate_options/3](predicate_options.html#current_predicate_options/3) | True when Options is the current active option declaration for PI on Arg. |
| [derive_predicate_options/0](predicate_options.html#derive_predicate_options/0) | Derive new predicate option declarations. |
| [derived_predicate_options/1](predicate_options.html#derived_predicate_options/1) | Derive predicate option declarations for a module. |
| [derived_predicate_options/3](predicate_options.html#derived_predicate_options/3) | Derive option arguments using static analysis. |
| [predicate_options/3](predicate_options.html#predicate_options/3) | Declare that the predicate PI processes options on Arg. |
| [retractall_predicate_options/0](predicate_options.html#retractall_predicate_options/0) | Remove all dynamically (derived) predicate options. |

### F.2.37 library(prologcoverage)

|  |  |
|----|----|
| [cov_load_data/2](prologcoverage.html#cov_load_data/2) | Reload coverage data from File. |
| [cov_property/1](prologcoverage.html#cov_property/1) | True when coverage analysis satisfies Property. |
| [cov_reset/0](prologcoverage.html#cov_reset/0) | Discard all collected coverage data. |
| [cov_save_data/2](prologcoverage.html#cov_save_data/2) | Save the coverage information to File. |
| [coverage/1](prologcoverage.html#coverage/1) | As call(Goal), collecting coverage information while Goal is running. |
| [coverage/2](prologcoverage.html#coverage/2) | Collect and optionally report coverage by Goal. |
| [report_hook/2](prologcoverage.html#report_hook/2) | This hook is called after the data collection. |
| [show_coverage/1](prologcoverage.html#show_coverage/1) | Show collected coverage data. |

### F.2.38 library(prologdebug)

|  |  |
|----|----|
| [debugging/0](debugger.html#debugging/0) | Report current status of the debugger. |
| [debugging_hook/1](prologdebug.html#debugging_hook/1) | Multifile hook that is called as forall(debugging_hook(DebugMode), true) and that may be used to extend the information printed from other debugging libraries. |
| [exception_hook/5](prologdebug.html#exception_hook/5) | Trap exceptions and consider whether or not to start the tracer. |
| [nospy/1](debugger.html#nospy/1) | Set/clear spy-points. |
| [nospyall/0](debugger.html#nospyall/0) | Set/clear spy-points. |
| [notrap/1](prologdebug.html#notrap/1) | Install a trap on error(Formal, Context) exceptions that unify. |
| [spy/1](debugger.html#spy/1) | Set/clear spy-points. |
| [trap/1](prologdebug.html#trap/1) | Install a trap on error(Formal, Context) exceptions that unify. |
| [trap_alias/2](prologdebug.html#trap_alias/2) | Define short hands for commonly used exceptions. |

### F.2.39 library(prologjiti)

|  |  |
|----|----|
| [jiti_list/0](prologjiti.html#jiti_list/0) | List the JITI (Just In Time Indexes) of selected predicates. |
| [jiti_list/1](prologjiti.html#jiti_list/1) | List the JITI (Just In Time Indexes) of selected predicates. |
| [jiti_suggest_modes/0](prologjiti.html#jiti_suggest_modes/0) | Propose modes for the predicates referenced by Spec. |
| [jiti_suggest_modes/1](prologjiti.html#jiti_suggest_modes/1) | Propose modes for the predicates referenced by Spec. |

### F.2.40 library(prologpack)

|  |  |
|----|----|
| [pack_info/1](prologpack.html#pack_info/1) | Print more detailed information about Pack. |
| [pack_install/1](prologpack.html#pack_install/1) | Install one or more packs from SpecOrList. |
| [pack_install/2](prologpack.html#pack_install/2) | Install one or more packs from SpecOrList. |
| [pack_install_local/3](prologpack.html#pack_install_local/3) | Install a number of packages in a local directory. |
| [pack_list/1](prologpack.html#pack_list/1) | Query package server and installed packages and display results. |
| [pack_list/2](prologpack.html#pack_list/2) | Query package server and installed packages and display results. |
| [pack_list_installed/0](prologpack.html#pack_list_installed/0) | List currently installed packages and report possible dependency issues. |
| [pack_property/2](prologpack.html#pack_property/2) | True when Property is a property of an installed Pack. |
| [pack_publish/2](prologpack.html#pack_publish/2) | Publish a package. |
| [pack_rebuild/0](prologpack.html#pack_rebuild/0) | Rebuild possible foreign components of Pack. |
| [pack_rebuild/1](prologpack.html#pack_rebuild/1) | Rebuild possible foreign components of Pack. |
| [pack_remove/1](prologpack.html#pack_remove/1) | Remove the indicated package. |
| [pack_remove/2](prologpack.html#pack_remove/2) | Remove the indicated package. |
| [pack_search/1](prologpack.html#pack_search/1) | Query package server and installed packages and display results. |
| [pack_upgrade/1](prologpack.html#pack_upgrade/1) | Upgrade Pack. |
| [pack_url_file/2](prologpack.html#pack_url_file/2) | True if File is a unique id for the referenced pack and version. |

### F.2.41 library(prologversions)

|  |  |
|----|----|
| [cmp_versions/3](prologversions.html#cmp_versions/3) | Compare to versions. |
| [require_prolog_version/2](prologversions.html#require_prolog_version/2) | Claim that the running Prolog version is at least version Required and provides the requested Features. |
| [require_version/3](prologversions.html#require_version/3) | Require Component to have version CmpRequired, while Component is know to have version Available. |

### F.2.42 library(prologtrace)

|  |  |
|----|----|
| [list_tracing/0](prologtrace.html#list_tracing/0) | List predicates we are currently tracing. |
| [notraceall/0](prologtrace.html#notraceall/0) | Remove all trace points. |
| [trace/1](prologtrace.html#trace/1) | Print passes through \_ports\_ of specified predicates. |
| [trace/2](prologtrace.html#trace/2) | Print passes through \_ports\_ of specified predicates. |
| [tracing/2](prologtrace.html#tracing/2) | True if Spec is traced using Ports. |

### F.2.43 library(prologxref)

|  |  |
|----|----|
| [prolog:called_by/2](prologxref.html#prolog:called_by/2) | *(hook)* Extend cross-referencer |
| xref_built_in/1 | Examine defined built-ins |
| [xref_called/3](prologxref.html#xref_called/3) | Examine called predicates |
| [xref_clean/1](prologxref.html#xref_clean/1) | Remove analysis of source |
| [xref_current_source/1](prologxref.html#xref_current_source/1) | Examine cross-referenced sources |
| [xref_defined/3](prologxref.html#xref_defined/3) | Examine defined predicates |
| [xref_exported/2](prologxref.html#xref_exported/2) | Examine exported predicates |
| [xref_module/2](prologxref.html#xref_module/2) | Module defined by source |
| [xref_source/1](prologxref.html#xref_source/1) | Cross-reference analysis of source |

### F.2.44 library(pairs)

|  |  |
|----|----|
| [group_pairs_by_key/2](pairs.html#group_pairs_by_key/2) | Group values with equivalent (==/2) consecutive keys. |
| [map_list_to_pairs/3](pairs.html#map_list_to_pairs/3) | Create a Key-Value list by mapping each element of List. |
| [pairs_keys/2](pairs.html#pairs_keys/2) | Remove the values from a list of Key-Value pairs. |
| [pairs_keys_values/3](pairs.html#pairs_keys_values/3) | True if Keys holds the keys of Pairs and Values the values. |
| [pairs_values/2](pairs.html#pairs_values/2) | Remove the keys from a list of Key-Value pairs. |
| [transpose_pairs/2](pairs.html#transpose_pairs/2) | Swap Key-Value to Value-Key. |

### F.2.45 library(pio)

#### F.2.45.1 library(pure_input)

|  |  |
|----|----|
| [phrase_from_file/2](pio.html#phrase_from_file/2) | Process the content of File using the DCG rule Grammar. |
| [phrase_from_file/3](pio.html#phrase_from_file/3) | As phrase_from_file/2, providing additional Options. |
| [phrase_from_stream/2](pio.html#phrase_from_stream/2) | Run Grammer against the character codes on Stream. |
| [stream_to_lazy_list/2](pio.html#stream_to_lazy_list/2) | Create a lazy list representing the character codes in Stream. |
| [lazy_list_character_count//1](pio.html#lazy_list_character_count//1) | True when CharCount is the current character count in the Lazy list. |
| [lazy_list_location//1](pio.html#lazy_list_location//1) | Determine current (error) location in a lazy list. |
| [syntax_error//1](pio.html#syntax_error//1) | Throw the syntax error Error at the current location of the input. |

### F.2.46 library(random)

|  |  |
|----|----|
| [getrand/1](random.html#getrand/1) | Query/set the state of the random generator. |
| [maybe/0](random.html#maybe/0) | Succeed/fail with equal probability (variant of maybe/1). |
| [maybe/1](random.html#maybe/1) | Succeed with probability P, fail with probability 1-P. |
| [maybe/2](random.html#maybe/2) | Succeed with probability K/N (variant of maybe/1). |
| [random/1](random.html#random/1) | Binds R to a new random float in the \_open\_ interval (0.0,1.0). |
| [random/3](random.html#random/3) | Generate a random integer or float in a range. |
| [random_between/3](random.html#random_between/3) | Binds R to a random integer in \[L,U\] (i.e., including both L and U). |
| [random_member/2](random.html#random_member/2) | X is a random member of List. |
| [random_numlist/4](random.html#random_numlist/4) | Unify List with an ascending list of integers between L and U (inclusive). |
| [random_perm2/4](random.html#random_perm2/4) | Does X=A,Y=B or X=B,Y=A with equal probability. |
| [random_permutation/2](random.html#random_permutation/2) | Permutation is a random permutation of List. |
| [random_select/3](random.html#random_select/3) | Randomly select or insert an element. |
| [random_subseq/3](random.html#random_subseq/3) | Selects a random subsequence Subseq of List, with Complement containing all elements of List that were not selected. |
| [randseq/3](random.html#randseq/3) | S is a list of K unique random integers in the range 1..N. |
| [randset/3](random.html#randset/3) | S is a sorted list of K unique random integers in the range 1..N. |
| [setrand/1](random.html#setrand/1) | Query/set the state of the random generator. |

### F.2.47 library(rbtrees)

|  |  |
|----|----|
| [is_rbtree/1](rbtrees.html#is_rbtree/1) | True if Term is a valid Red-Black tree. |
| [list_to_rbtree/2](rbtrees.html#list_to_rbtree/2) | Tree is the red-black tree corresponding to the mapping in List, which should be a list of Key-Value pairs. |
| [ord_list_to_rbtree/2](rbtrees.html#ord_list_to_rbtree/2) | Tree is the red-black tree corresponding to the mapping in list List, which should be a list of Key-Value pairs. |
| [rb_apply/4](rbtrees.html#rb_apply/4) | If the value associated with key Key is Val0 in Tree, and if call(G,Val0,ValF) holds, then NewTree differs from Tree only in that Key is associated with value ValF in tree NewTree. |
| [rb_clone/3](rbtrees.html#rb_clone/3) | ‘Clone’the red-back tree TreeIn into a new tree TreeOut with the same keys as the original but with all values set to unbound values. |
| [rb_del_max/4](rbtrees.html#rb_del_max/4) | Delete the largest element from the tree Tree, returning the key Key, the value Val associated with the key and a new tree NewTree. |
| [rb_del_min/4](rbtrees.html#rb_del_min/4) | Delete the least element from the tree Tree, returning the key Key, the value Val associated with the key and a new tree NewTree. |
| [rb_delete/3](rbtrees.html#rb_delete/3) | Delete element with key Key from the tree Tree, returning the value Val associated with the key and a new tree NewTree. |
| [rb_delete/4](rbtrees.html#rb_delete/4) | Same as rb_delete(Tree, Key, NewTree), but also unifies Val with the value associated with Key in Tree. |
| [rb_empty/1](rbtrees.html#rb_empty/1) | Succeeds if Tree is an empty Red-Black tree. |
| [rb_fold/4](rbtrees.html#rb_fold/4) | Fold the given predicate over all the key-value pairs in Tree, starting with initial state State0 and returning the final state State. |
| [rb_in/3](rbtrees.html#rb_in/3) | True when Key-Value is a key-value pair in red-black tree Tree. |
| [rb_insert/4](rbtrees.html#rb_insert/4) | Add an element with key Key and Value to the tree Tree creating a new red-black tree NewTree. |
| [rb_insert_new/4](rbtrees.html#rb_insert_new/4) | Add a new element with key Key and Value to the tree Tree creating a new red-black tree NewTree. |
| [rb_keys/2](rbtrees.html#rb_keys/2) | Keys is unified with an ordered list of all keys in the Red-Black tree Tree. |
| [rb_lookup/3](rbtrees.html#rb_lookup/3) | True when Value is associated with Key in the Red-Black tree Tree. |
| [rb_map/2](rbtrees.html#rb_map/2) | True if call(Goal, Value) is true for all nodes in T. |
| [rb_map/3](rbtrees.html#rb_map/3) | For all nodes Key in the tree Tree, if the value associated with key Key is Val0 in tree Tree, and if call(G,Val0,ValF) holds, then the value associated with Key in NewTree is ValF. |
| [rb_max/3](rbtrees.html#rb_max/3) | Key is the maximal key in Tree, and is associated with Val. |
| [rb_min/3](rbtrees.html#rb_min/3) | Key is the minimum key in Tree, and is associated with Val. |
| [rb_new/1](rbtrees.html#rb_new/1) | Create a new Red-Black tree Tree. |
| [rb_next/4](rbtrees.html#rb_next/4) | Next is the next element after Key in Tree, and is associated with Val. |
| [rb_partial_map/4](rbtrees.html#rb_partial_map/4) | For all nodes Key in Keys, if the value associated with key Key is Val0 in tree Tree, and if call(G,Val0,ValF) holds, then the value associated with Key in NewTree is ValF, otherwise it is the value associated with the key in Tree. |
| [rb_previous/4](rbtrees.html#rb_previous/4) | Previous is the previous element after Key in Tree, and is associated with Val. |
| [rb_size/2](rbtrees.html#rb_size/2) | Size is the number of elements in Tree. |
| [rb_update/4](rbtrees.html#rb_update/4) | Tree NewTree is tree Tree, but with value for Key associated with NewVal. |
| [rb_update/5](rbtrees.html#rb_update/5) | Same as =`|`rb_update(Tree, Key, NewVal, NewTree)`|`= but also unifies OldVal with the value associated with Key in Tree. |
| [rb_visit/2](rbtrees.html#rb_visit/2) | Pairs is an infix visit of tree Tree, where each element of Pairs is of the form Key-Value. |

### F.2.48 library(readutil)

|  |  |
|----|----|
| [read_file_to_codes/3](readutil.html#read_file_to_codes/3) | Read the file Spec into a list of Codes. |
| [read_file_to_string/3](readutil.html#read_file_to_string/3) | Read the file Spec into a the string String. |
| [read_file_to_terms/3](readutil.html#read_file_to_terms/3) | Read the file Spec into a list of terms. |
| [read_line_to_codes/2](readutil.html#read_line_to_codes/2) | Read the next line of input from Stream. |
| [read_line_to_codes/3](readutil.html#read_line_to_codes/3) | Difference-list version to read an input line to a list of character codes. |
| [read_line_to_string/2](readutil.html#read_line_to_string/2) | Read the next line from Stream into String. |
| [read_stream_to_codes/2](readutil.html#read_stream_to_codes/2) | Read input from Stream to a list of character codes. |
| [read_stream_to_codes/3](readutil.html#read_stream_to_codes/3) | Read input from Stream to a list of character codes. |

### F.2.49 library(record)

|                                  |                               |
|----------------------------------|-------------------------------|
| [record/1](record.html#record/1) | Define named fields in a term |

### F.2.50 library(registry)

This library is only available on Windows systems.

|  |  |
|----|----|
| [registry_get_key/2](registry.html#registry_get_key/2) | Get principal value of key |
| [registry_get_key/3](registry.html#registry_get_key/3) | Get associated value of key |
| [registry_set_key/2](registry.html#registry_set_key/2) | Set principal value of key |
| [registry_set_key/3](registry.html#registry_set_key/3) | Set associated value of key |
| [registry_delete_key/1](registry.html#registry_delete_key/1) | Remove a key |
| [shell_register_file_type/4](registry.html#shell_register_file_type/4) | Register a file-type |
| [shell_register_dde/6](registry.html#shell_register_dde/6) | Register DDE action |
| [shell_register_prolog/1](registry.html#shell_register_prolog/1) | Register Prolog |

### F.2.51 library(rwlocks)

|  |  |
|----|----|
| [with_rwlock/3](rwlocks.html#with_rwlock/3) | Run Goal, synchronized with LockId in ModeSpec. |
| [with_rwlock/4](rwlocks.html#with_rwlock/4) | Run Goal, synchronized with LockId in ModeSpec. |

### F.2.52 library(settings)

|  |  |
|----|----|
| [convert_setting_text/3](settings.html#convert_setting_text/3) | Converts from textual form to Prolog Value. |
| [current_setting/1](settings.html#current_setting/1) | True if Setting is a currently defined setting. |
| [env/2](settings.html#env/2) | Evaluate environment variables on behalf of arithmetic expressions. |
| [env/3](settings.html#env/3) | Evaluate environment variables on behalf of arithmetic expressions. |
| [list_settings/0](settings.html#list_settings/0) | List settings to =current_output=. |
| [list_settings/1](settings.html#list_settings/1) | List settings to =current_output=. |
| [load_settings/1](settings.html#load_settings/1) | Load local settings from File. |
| [load_settings/2](settings.html#load_settings/2) | Load local settings from File. |
| [restore_setting/1](settings.html#restore_setting/1) | Restore the value of setting Name to its default. |
| [save_settings/0](settings.html#save_settings/0) | Save modified settings to File. |
| [save_settings/1](settings.html#save_settings/1) | Save modified settings to File. |
| [set_setting/2](settings.html#set_setting/2) | Change a setting. |
| [set_setting_default/2](settings.html#set_setting_default/2) | Change the default for a setting. |
| [setting/2](settings.html#setting/2) | True when Name is a currently defined setting with Value. |
| [setting/4](settings.html#setting/4) | Define a setting. |
| [setting_property/2](settings.html#setting_property/2) | Query currently defined settings. |

### F.2.53 library(simplex)

|  |  |
|----|----|
| [assignment/2](simplex.html#assignment/2) | Solve assignment problem |
| [constraint/3](simplex.html#constraint/3) | Add linear constraint to state |
| [constraint/4](simplex.html#constraint/4) | Add named linear constraint to state |
| [constraint_add/4](simplex.html#constraint_add/4) | Extend a named constraint |
| [gen_state/1](simplex.html#gen_state/1) | Create empty linear program |
| [maximize/3](simplex.html#maximize/3) | Maximize objective function in to linear constraints |
| [minimize/3](simplex.html#minimize/3) | Minimize objective function in to linear constraints |
| [objective/2](simplex.html#objective/2) | Fetch value of objective function |
| [shadow_price/3](simplex.html#shadow_price/3) | Fetch shadow price in solved state |
| [transportation/4](simplex.html#transportation/4) | Solve transportation problem |
| [variable_value/3](simplex.html#variable_value/3) | Fetch value of variable in solved state |

### F.2.54 library(statistics)

|  |  |
|----|----|
| [call_time/2](statistics.html#call_time/2) | Call Goal as call/1, unifying Time with a dict that provides information on the resource usage. |
| [call_time/3](statistics.html#call_time/3) | Call Goal as call/1, unifying Time with a dict that provides information on the resource usage. |
| [statistics/0](statistics.html#statistics/0) | Print information about resource usage using print_message/2. |
| [statistics/1](statistics.html#statistics/1) | Stats is a dict representing the same information as statistics/0. |
| [thread_statistics/2](statistics.html#thread_statistics/2) | Obtain statistical information about a single thread. |
| [time/1](statistics.html#time/1) | Execute Goal, reporting statistics to the user. |

### F.2.55 library(tableutils)

|  |  |
|----|----|
| [summarize_idg/0](tableutil.html#summarize_idg/0) | Implements XSB's statistics(summarize_idg). |
| [summarize_idg/1](tableutil.html#summarize_idg/1) | Implements XSB's statistics(summarize_idg). |
| [table_statistics/0](tableutil.html#table_statistics/0) | Print a summary of statistics relevant to tabling. |
| [table_statistics/1](tableutil.html#table_statistics/1) | Print a summary for the statistics of all tables for subgoals of Variant. |
| [table_statistics/2](tableutil.html#table_statistics/2) | Give summary statistics for the tables associated with all subgoals of Variant. |
| [table_statistics/3](tableutil.html#table_statistics/3) | Give summary statistics for the tables associated with all subgoals of Variant. |
| [table_statistics_by_predicate/0](tableutil.html#table_statistics_by_predicate/0) | Print statistics on memory usage and lookups per predicate. |
| [table_statistics_by_predicate/1](tableutil.html#table_statistics_by_predicate/1) | Print statistics on memory usage and lookups per predicate. |
| [tdump/0](tableutil.html#tdump/0) | Dump all tables and their status that \_unify\_ with Goal. |
| [tdump/1](tableutil.html#tdump/1) | Dump all tables and their status that \_unify\_ with Goal. |
| [tdump/2](tableutil.html#tdump/2) | Dump all tables and their status that \_unify\_ with Goal. |
| [tidg/0](tableutil.html#tidg/0) | Dump the incremental dependency graph. |
| [tidg/1](tableutil.html#tidg/1) | Dump the incremental dependency graph. |
| [tstat/2](tableutil.html#tstat/2) | Print the top-N (for positive Top) or bottom-N (for negative Top) for‘Stat\` for all tabled subgoals of Variant (or all tabled subgoals for tstat/2). |
| [tstat/3](tableutil.html#tstat/3) | Print the top-N (for positive Top) or bottom-N (for negative Top) for‘Stat\` for all tabled subgoals of Variant (or all tabled subgoals for tstat/2). |

### F.2.56 library(terms)

|  |  |
|----|----|
| [foldsubterms/4](terms.html#foldsubterms/4) | The predicate foldsubterms/5 calls call(Goal4, SubTerm1, SubTerm2, StateIn, StateOut) for each subterm, including variables, in Term1. |
| [foldsubterms/5](terms.html#foldsubterms/5) | The predicate foldsubterms/5 calls call(Goal4, SubTerm1, SubTerm2, StateIn, StateOut) for each subterm, including variables, in Term1. |
| [mapargs/3](terms.html#mapargs/3) | Term1 and Term2 have the same functor (name/arity) and for each matching pair of arguments call(Goal, A1, A2) is true. |
| [mapsubterms/3](terms.html#mapsubterms/3) | Recursively map sub terms of Term1 into subterms of Term2 for every pair for which call(Goal, ST1, ST2) succeeds. |
| [mapsubterms_var/3](terms.html#mapsubterms_var/3) | Recursively map sub terms of Term1 into subterms of Term2 for every pair for which call(Goal, ST1, ST2) succeeds. |
| [same_functor/2](terms.html#same_functor/2) | True when Term1 and Term2 are terms that have the same functor (Name/Arity). |
| [same_functor/3](terms.html#same_functor/3) | True when Term1 and Term2 are terms that have the same functor (Name/Arity). |
| [same_functor/4](terms.html#same_functor/4) | True when Term1 and Term2 are terms that have the same functor (Name/Arity). |
| [subsumes/2](terms.html#subsumes/2) | True if Generic is unified to Specific without changing Specific. |
| [subsumes_chk/2](terms.html#subsumes_chk/2) | True if Generic can be made equivalent to Specific without changing Specific. |
| [term_factorized/3](terms.html#term_factorized/3) | Is true when Skeleton is Term where all subterms that appear multiple times are replaced by a variable and Substitution is a list of Var=Value that provides the subterm at the location Var. |
| [term_size/2](terms.html#term_size/2) | True if Size is the size in \_cells\_ occupied by Term on the global (term) stack. |
| [term_subsumer/3](compare.html#term_subsumer/3) | General is the most specific term that is a generalisation of Special1 and Special2. |
| [variant/2](terms.html#variant/2) | Same as SWI-Prolog =`|`Term1 `=@=` Term2`|`=. |

### F.2.57 library(ugraphs)

|  |  |
|----|----|
| [add_edges/3](ugraphs.html#add_edges/3) | Unify NewGraph with a new graph obtained by adding the list of Edges to Graph. |
| [add_vertices/3](ugraphs.html#add_vertices/3) | Unify NewGraph with a new graph obtained by adding the list of Vertices to Graph. |
| [complement/2](ugraphs.html#complement/2) | UGraphOut is a ugraph with an edge between all vertices that are \_not\_ connected in UGraphIn and all edges from UGraphIn removed. |
| [compose/3](ugraphs.html#compose/3) | Compose NewGraph by connecting the \_drains\_ of LeftGraph to the \_sources\_ of RightGraph. |
| [connect_ugraph/3](ugraphs.html#connect_ugraph/3) | Adds Start as an additional vertex that is connected to all vertices in UGraphIn. |
| [del_edges/3](ugraphs.html#del_edges/3) | Unify NewGraph with a new graph obtained by removing the list of Edges from Graph. |
| [del_vertices/3](ugraphs.html#del_vertices/3) | Unify NewGraph with a new graph obtained by deleting the list of Vertices and all the edges that start from or go to a vertex in Vertices to the Graph. |
| [edges/2](ugraphs.html#edges/2) | Unify Edges with all edges appearing in Graph. |
| [neighbors/3](ugraphs.html#neighbors/3) | Neigbours is a sorted list of the neighbours of Vertex in Graph. |
| [neighbours/3](ugraphs.html#neighbours/3) | Neigbours is a sorted list of the neighbours of Vertex in Graph. |
| [reachable/3](ugraphs.html#reachable/3) | True when Vertices is an ordered set of vertices reachable in UGraph, including Vertex. |
| [top_sort/2](ugraphs.html#top_sort/2) | Sort vertices topologically. |
| [transitive_closure/2](ugraphs.html#transitive_closure/2) | Generate the graph Closure as the transitive closure of Graph. |
| [transpose_ugraph/2](ugraphs.html#transpose_ugraph/2) | Unify NewGraph with a new graph obtained from Graph by replacing all edges of the form V1-V2 by edges of the form V2-V1. |
| [ugraph_layers/2](ugraphs.html#ugraph_layers/2) | Sort vertices topologically. |
| [ugraph_union/3](ugraphs.html#ugraph_union/3) | NewGraph is the union of Graph1 and Graph2. |
| [vertices/2](ugraphs.html#vertices/2) | Unify Vertices with all vertices appearing in Graph. |
| [vertices_edges_to_ugraph/3](ugraphs.html#vertices_edges_to_ugraph/3) | Create a UGraph from Vertices and Edges. |

|  |  |
|----|----|
| [vertices_edges_to_ugraph/3](ugraphs.html#vertices_edges_to_ugraph/3) | Create unweighted graph |
| [vertices/2](ugraphs.html#vertices/2) | Find vertices in graph |
| [edges/2](ugraphs.html#edges/2) | Find edges in graph |
| [add_vertices/3](ugraphs.html#add_vertices/3) | Add vertices to graph |
| [del_vertices/3](ugraphs.html#del_vertices/3) | Delete vertices from graph |
| [add_edges/3](ugraphs.html#add_edges/3) | Add edges to graph |
| [del_edges/3](ugraphs.html#del_edges/3) | Delete edges from graph |
| [transpose_ugraph/2](ugraphs.html#transpose_ugraph/2) | Invert the direction of all edges |
| [neighbors/3](ugraphs.html#neighbors/3) | Find neighbors of vertice |
| [neighbours/3](ugraphs.html#neighbours/3) | Find neighbors of vertice |
| [complement/2](ugraphs.html#complement/2) | Inverse presense of edges |
| [compose/3](ugraphs.html#compose/3) |  |
| [top_sort/2](ugraphs.html#top_sort/2) | Sort graph topologically |
| top_sort/3 | Sort graph topologically |
| [transitive_closure/2](ugraphs.html#transitive_closure/2) | Create transitive closure of graph |
| [reachable/3](ugraphs.html#reachable/3) | Find all reachable vertices |
| [ugraph_union/3](ugraphs.html#ugraph_union/3) | Union of two graphs |

### F.2.58 library(url)

|  |  |
|----|----|
| [file_name_to_url/2](url.html#file_name_to_url/2) | Translate between a filename and a file:`//` URL. |
| [global_url/3](url.html#global_url/3) | Translate a possibly relative URL into an absolute one. |
| [http_location/2](url.html#http_location/2) | Construct or analyze an HTTP location. |
| [is_absolute_url/1](url.html#is_absolute_url/1) | True if URL is an absolute URL. |
| [parse_url/2](url.html#parse_url/2) | Construct or analyse a URL. |
| [parse_url/3](url.html#parse_url/3) | Similar to parse_url/2 for relative URLs. |
| [parse_url_search/2](url.html#parse_url_search/2) | Construct or analyze an HTTP search specification. |
| [set_url_encoding/2](url.html#set_url_encoding/2) | Query and set the encoding for URLs. |
| [url_iri/2](url.html#url_iri/2) | Convert between a URL, encoding in US-ASCII and an IRI. |
| [www_form_encode/2](url.html#www_form_encode/2) | En/decode to/from application/x-www-form-encoded. |

### F.2.59 library(writef)

|  |  |
|----|----|
| [swritef/2](writef.html#swritef/2) | Use writef/1 or writef/2 and write the result to a \_string\_. |
| [swritef/3](writef.html#swritef/3) | Use writef/1 or writef/2 and write the result to a \_string\_. |
| [writef/1](writef.html#writef/1) | Formatted write to the =current_output=. |
| [writef/2](writef.html#writef/2) | Formatted write to the =current_output=. |

### F.2.60 library(www_browser)

|  |  |
|----|----|
| [expand_url_path/2](wwwbrowser.html#expand_url_path/2) | Expand URL specifications similar to absolute_file_name/3. |
| [known_browser/2](wwwbrowser.html#known_browser/2) | True if browser FileBaseName has a remote protocol compatible to Compatible. |
| [www_open_url/1](wwwbrowser.html#www_open_url/1) | Open URL in running version of the users’browser or start a new browser. |

### F.2.61 library(solution_sequences)

|  |  |
|----|----|
| [call_nth/2](solutionsequences.html#call_nth/2) | True when Goal succeeded for the Nth time. |
| [distinct/1](solutionsequences.html#distinct/1) | True if Goal is true and no previous solution of Goal bound Witness to the same value. |
| [distinct/2](solutionsequences.html#distinct/2) | True if Goal is true and no previous solution of Goal bound Witness to the same value. |
| [group_by/4](solutionsequences.html#group_by/4) | Group bindings of Template that have the same value for By. |
| [limit/2](solutionsequences.html#limit/2) | Limit the number of solutions. |
| [offset/2](solutionsequences.html#offset/2) | Ignore the first Count solutions. |
| [order_by/2](solutionsequences.html#order_by/2) | Order solutions according to Spec. |
| [reduced/1](solutionsequences.html#reduced/1) | Similar to distinct/1, but does not guarantee unique results in return for using a limited amount of memory. |
| [reduced/3](solutionsequences.html#reduced/3) | Similar to distinct/1, but does not guarantee unique results in return for using a limited amount of memory. |

### F.2.62 library(thread)

|  |  |
|----|----|
| [call_in_thread/2](thread.html#call_in_thread/2) | Run Goal as an interrupt in the context of Thread. |
| [call_in_thread/3](thread.html#call_in_thread/3) | Run Goal as an interrupt in the context of Thread. |
| [concurrent/3](thread.html#concurrent/3) | Run Goals in parallel using N threads. |
| [concurrent_and/2](thread.html#concurrent_and/2) | Concurrent version of‘(Generator,Test)\`. |
| [concurrent_and/3](thread.html#concurrent_and/3) | Concurrent version of‘(Generator,Test)\`. |
| [concurrent_forall/2](thread.html#concurrent_forall/2) | True when Action is true for all solutions of Generate. |
| [concurrent_forall/3](thread.html#concurrent_forall/3) | True when Action is true for all solutions of Generate. |
| [concurrent_maplist/2](thread.html#concurrent_maplist/2) | Concurrent version of maplist/2. |
| [concurrent_maplist/3](thread.html#concurrent_maplist/3) | Concurrent version of maplist/2. |
| [concurrent_maplist/4](thread.html#concurrent_maplist/4) | Concurrent version of maplist/2. |
| [first_solution/3](thread.html#first_solution/3) | Try alternative solvers concurrently, returning the first answer. |

### F.2.63 library(thread_pool)

|  |  |
|----|----|
| [create_pool/1](threadpool.html#create_pool/1) | Hook to create a thread pool lazily. |
| [current_thread_pool/1](threadpool.html#current_thread_pool/1) | True if Name refers to a defined thread pool. |
| [thread_create_in_pool/4](threadpool.html#thread_create_in_pool/4) | Create a thread in Pool. |
| [thread_pool_create/3](threadpool.html#thread_pool_create/3) | Create a pool of threads. |
| [thread_pool_destroy/1](threadpool.html#thread_pool_destroy/1) | Destroy the thread pool named Name. |
| [thread_pool_property/2](threadpool.html#thread_pool_property/2) | True if Property is a property of thread pool Name. |
| [worker_exitted/3](threadpool.html#worker_exitted/3) | It is possible that’\_\_thread_pool_manager’no longer exists while closing down the process because the manager was killed before the worker. |

### F.2.64 library(varnumbers)

|  |  |
|----|----|
| [max_var_number/3](varnumbers.html#max_var_number/3) | True when Max is the max of Start and the highest numbered \$VAR(N) term. |
| [numbervars/1](varnumbers.html#numbervars/1) | Number variables in Term using \$VAR(N). |
| [varnumbers/2](varnumbers.html#varnumbers/2) | Inverse of numbervars/1. |
| [varnumbers/3](varnumbers.html#varnumbers/3) | Inverse of numbervars/3. |
| [varnumbers_names/3](varnumbers.html#varnumbers_names/3) | If Term is a term with numbered and named variables using the reserved term’\$VAR’(X), Copy is a copy of Term where each’\$VAR’(X) is consistently replaced by a fresh variable and Bindings is a list‘X = Var\`, relating the‘X\` terms with the variable it is mapped to. |

### F.2.65 library(yall)

|  |  |
|----|----|
| [//2](yall.html#//2) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//3](yall.html#//3) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//4](yall.html#//4) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//5](yall.html#//5) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//6](yall.html#//6) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//7](yall.html#//7) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//8](yall.html#//8) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [//9](yall.html#//9) | Shorthand for‘Free/\[\]`>>`Lambda\`. |
| [\>\>/2](yall.html#%3E%3E/2) | Calls a copy of Lambda. |
| [\>\>/3](yall.html#%3E%3E/3) | Calls a copy of Lambda. |
| [\>\>/4](yall.html#%3E%3E/4) | Calls a copy of Lambda. |
| [\>\>/5](yall.html#%3E%3E/5) | Calls a copy of Lambda. |
| [\>\>/6](yall.html#%3E%3E/6) | Calls a copy of Lambda. |
| [\>\>/7](yall.html#%3E%3E/7) | Calls a copy of Lambda. |
| [\>\>/8](yall.html#%3E%3E/8) | Calls a copy of Lambda. |
| [\>\>/9](yall.html#%3E%3E/9) | Calls a copy of Lambda. |
| [is_lambda/1](yall.html#is_lambda/1) | True if Term is a valid Lambda expression. |
| [lambda_calls/2](yall.html#lambda_calls/2) | Goal is the goal called if call/N is applied to LambdaExpression, where ExtraArgs are the additional arguments to call/N. |
| [lambda_calls/3](yall.html#lambda_calls/3) | Goal is the goal called if call/N is applied to LambdaExpression, where ExtraArgs are the additional arguments to call/N. |
