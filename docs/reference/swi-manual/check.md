
## A.7 library(check): Consistency checking

See also  
\- [gxref/0](xref.html#gxref/0) provides a graphical cross referencer  
- PceEmacs performs real time consistency checks while you edit  
- `library(prolog_xref)` implements‘offline’cross-referencing  
- `library(prolog_codewalk)` implements‘online’analysis

This library provides some consistency checks for the loaded Prolog program. The predicate [make/0](consulting.html#make/0) runs [list_undefined/0](check.html#list_undefined/0) to find undefined predicates in‘user’modules.

\[det\]**check**  
Run all consistency checks defined by [checker/2](check.html#checker/2). Checks enabled by default are:

- [list_undefined/0](check.html#list_undefined/0) reports undefined predicates
- [list_trivial_fails/0](check.html#list_trivial_fails/0) reports calls for which there is no matching clause.
- [list_format_errors/0](check.html#list_format_errors/0) reports mismatches in [format/2](format.html#format/2),3 templates and the list of arguments.
- [list_redefined/0](check.html#list_redefined/0) reports predicates that have a local definition and a global definition. Note that these are **not** errors.
- [list_void_declarations/0](check.html#list_void_declarations/0) reports on predicates with defined properties, but no clauses.
- [list_autoload/0](check.html#list_autoload/0) lists predicates that will be defined at runtime using the autoloader.
- [check_predicate_options/0](predicate_options.html#check_predicate_options/0) tests for options passed to predicates such as [open/4](IO.html#open/4) that are unknown or are used with an invalid argument.

The checker can be expanded or restricted by modifying the dynamic multifile hook [checker/2](check.html#checker/2).

The checker may be used in batch, e.g., for CI workflows by calling SWI-Prolog as below. Note that by using `-l` to load the program, the program is not started if it used [initialization/2](consulting.html#initialization/2) of type `main` to start the program.

``` code
swipl -q --on-warning=status --on-error=status \
      -g check -t halt -l myprogram.pl
```

\[det\]**list_undefined**  
\[det\]**list_undefined**(`+Options`)  
Report undefined predicates. This predicate finds undefined predicates by decompiling and analyzing the body of all clauses. `Options`:

**module_class**(`+Classes`)  
Process modules of the given `Classes`. The default for classes is `[user]`. For example, to include the libraries into the examination, use `[user,library]`.

See also  
\- [gxref/0](xref.html#gxref/0) provides a graphical cross-referencer.  
- [make/0](consulting.html#make/0) calls [list_undefined/0](check.html#list_undefined/0)

\[det\]**list_autoload**  
Report predicates that may be auto-loaded. These are predicates that are not defined, but will be loaded on demand if referenced.

See also  
autoload/0

To be done  
This predicate uses an older mechanism for finding undefined predicates. Should be synchronized with list undefined.

**list_redefined**  
Lists predicates that are defined in the global module `user` as well as in a normal module; that is, predicates for which the local definition overrules the global default definition.

\[det\]**list_cross_module_calls**  
List calls from one module to another using Module:Goal where the callee is not defined exported, public or multifile, i.e., where the callee should be considered *private*.

\[det\]**list_void_declarations**  
List predicates that have declared attributes, but no clauses.

\[det\]**list_trivial_fails**  
\[det\]**list_trivial_fails**(`+Options`)  
List goals that trivially fail because there is no matching clause. `Options`:

**module_class**(`+Classes`)  
Process modules of the given `Classes`. The default for classes is `[user]`. For example, to include the libraries into the examination, use `[user,library]`.

\[multifile\]**trivial_fail_goal**(`:Goal`)  
Multifile hook that tells [list_trivial_fails/0](check.html#list_trivial_fails/0) to accept `Goal` as valid.

\[det\]**list_strings**  
\[det\]**list_strings**(`+Options`)  
List strings that appear in clauses. This predicate is used to find portability issues for changing the Prolog flag `double_quotes` from `codes` to `string`, creating packed string objects. Warnings may be suppressed using the following multifile hooks:

- [string_predicate/1](check.html#string_predicate/1) to stop checking certain predicates
- [valid_string_goal/1](check.html#valid_string_goal/1) to tell the checker that a goal is safe.

See also  
Prolog flag `double_quotes`.

\[det\]**list_rationals**  
\[det\]**list_rationals**(`+Options`)  
List rational numbers that appear in clauses. This predicate is used to find portability issues for changing the Prolog flag `rational_syntax` to `natural`, creating rational numbers from \<integer\>/\<nonneg\>. `Options`:

**module_class**(`+Classes`)  
Determines the modules classes processed. By default only user code is processed. See prolog_program_clause/2.

**arithmetic**(`+Bool`)  
If `true` (default `false`) also warn on rationals appearing in arithmetic expressions.

See also  
Prolog flag [rational_syntax](flags.html#flag:rational_syntax) and `prefer_rationals`.

\[det\]**list_format_errors**  
\[det\]**list_format_errors**(`+Options`)  
List argument errors for [format/2](format.html#format/2),3.

\[multifile\]**string_predicate**(`:PredicateIndicator`)  
Multifile hook to disable [list_strings/0](string.html#list_strings/0) on the given predicate. This is typically used for facts that store strings.

\[semidet,multifile\]**valid_string_goal**(`+Goal`)  
Multifile hook that qualifies `Goal` as valid for [list_strings/0](string.html#list_strings/0). For example, `format("Hello world~n")` is considered proper use of string constants.

\[nondet,multifile\]**checker**(`:Goal, +Message:text`)  
Register code validation routines. Each clause defines a `Goal` which performs a consistency check executed by [check/0](check.html#check/0). `Message` is a short description of the check. For example, assuming the `my_checks` module defines a predicate list_format_mistakes/0:

``` code
:- multifile check:checker/2.
check:checker(my_checks:list_format_mistakes,
              "errors with format/2 arguments").
```

The predicate is dynamic, so you can disable checks with [retract/1](db.html#retract/1). For example, to stop reporting redefined predicates:

``` code
retract(check:checker(list_redefined,_)).
```

\[det\]**list_confusable_identifiers**  
\[det\]**list_confusable_identifiers**(`+Options`)  
Walk every clause in the modules selected by `module_class` (default `[user]`) and warn on atoms whose written form is a possible UTS \#39 spoof. Two issues are reported:

- an atom whose unicode_restriction_level/2 is not safe (worse than `single_script`), reported at the first source location where it occurs;

- groups of distinct atoms with equal unicode_skeleton/2 — two identifiers in the same program that look the same.

  Only atoms whose written form does not require quotes are considered (see’\$needs_quotes’/1). Loaded via the multifile check:checker/2 hook, so a plain [check/0](check.html#check/0) runs this check too.
