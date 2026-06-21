
## 9.2 CHR Syntax and Semantics

### 9.2.1 Syntax of CHR rules

``` code
rules --> rule, rules ; [].

rule --> name, actual_rule, pragma, [atom('.')].

name --> atom, [atom('@')] ; [].

actual_rule --> simplification_rule.
actual_rule --> propagation_rule.
actual_rule --> simpagation_rule.

simplification_rule --> head, [atom('<=>')], guard, body.
propagation_rule --> head, [atom('==>')], guard, body.
simpagation_rule --> head, [atom('\')], head, [atom('<=>')],
                     guard, body.

head --> constraints.

constraints --> constraint, constraint_id.
constraints --> constraint, constraint_id,
                [atom(',')], constraints.

constraint --> compound_term.

constraint_id --> [].
constraint_id --> [atom('#')], variable.
constraint_id --> [atom('#')], [atom('passive')] .

guard --> [] ; goal, [atom('|')].

body --> goal.

pragma --> [].
pragma --> [atom('pragma')], actual_pragmas.

actual_pragmas --> actual_pragma.
actual_pragmas --> actual_pragma, [atom(',')], actual_pragmas.

actual_pragma --> [atom('passive(')], variable, [atom(')')].
```

Note that the guard of a rule may not contain any goal that binds a variable in the head of the rule with a non-variable or with another variable in the head of the rule. It may, however, bind variables that do not appear in the head of the rule, e.g. an auxiliary variable introduced in the guard.

### 9.2.2 Semantics of CHR

In this subsection the operational semantics of CHR in Prolog are presented informally. They do not differ essentially from other CHR systems.

When a constraint is called, it is considered an active constraint and the system will try to apply the rules to it. Rules are tried and executed sequentially in the order they are written.

A rule is conceptually tried for an active constraint in the following way. The active constraint is matched with a constraint in the head of the rule. If more constraints appear in the head, they are looked for among the suspended constraints, which are called passive constraints in this context. If the necessary passive constraints can be found and all match with the head of the rule and the guard of the rule succeeds, then the rule is committed and the body of the rule executed. If not all the necessary passive constraints can be found, or the matching or the guard fails, then the body is not executed and the process of trying and executing simply continues with the following rules. If for a rule there are multiple constraints in the head, the active constraint will try the rule sequentially multiple times, each time trying to match with another constraint.

This process ends either when the active constraint disappears, i.e. it is removed by some rule, or after the last rule has been processed. In the latter case the active constraint becomes suspended.

A suspended constraint is eligible as a passive constraint for an active constraint. The other way it may interact again with the rules is when a variable appearing in the constraint becomes bound to either a non-variable or another variable involved in one or more constraints. In that case the constraint is triggered, i.e. it becomes an active constraint and all the rules are tried.

**Rule Types** There are three different kinds of rules, each with its specific semantics:

- *simplification*  
  The simplification rule removes the constraints in its head and calls its body.
- *propagation*  
  The propagation rule calls its body exactly once for the constraints in its head.
- *simpagation*  
  The simpagation rule removes the constraints in its head after the `\` and then calls its body. It is an optimization of simplification rules of the form: \[constraints_1, constraints_2 \<=\> constraints_1, body \] Namely, in the simpagation form: \[ constraints_1 `\`constraints_2 \<=\> body \] The *`constraints`*`_1` constraints are not called in the body.

**Rule Names** Naming a rule is optional and has no semantic meaning. It only functions as documentation for the programmer.

**Pragmas** The semantics of the pragmas are:

**passive**(`Identifier`)  
The constraint in the head of a rule `Identifier` can only match a passive constraint in that rule. There is an abbreviated syntax for this pragma. Instead of:

``` code
                ..., c # Id, ... <=> ... pragma passive(Id)
```

you can also write

``` code
                ..., c # passive, ... <=> ...
```

Additional pragmas may be released in the future.

:- **chr_option**(`+Option, +Value`)  
It is possible to specify options that apply to all the CHR rules in the module. Options are specified with the [`chr_option/2`](chr-syntaxandsemantics.html#chr_option/2) declaration:

``` code
:- chr_option(Option,Value).
```

and may appear in the file anywhere after the first constraints declaration.

Available options are:

**check_guard_bindings**  
This option controls whether guards should be checked for (illegal) variable bindings or not. Possible values for this option are `on` to enable the checks, and `off` to disable the checks. If this option is on, any guard fails when it binds a variable that appears in the head of the rule. When the option is off (default), the behaviour of a binding in the guard is undefined.

**optimize**  
This option controls the degree of optimization. Possible values are `full` to enable all available optimizations, and `off` (default) to disable all optimizations. The default is derived from the SWI-Prolog flag [optimise](flags.html#flag:optimise), where `true` is mapped to `full`. Therefore the command line option **-O** provides full CHR optimization. If optimization is enabled, debugging must be disabled.

**debug**  
This option enables or disables the possibility to debug the CHR code. Possible values are `on` (default) and `off`. See [section 9.4](chr-debugging.html#sec:9.4) for more details on debugging. The default is derived from the Prolog flag [generate_debug_info](flags.html#flag:generate_debug_info), which is `true` by default. See **--no-debug**. If debugging is enabled, optimization must be disabled.
