
## 9.8 CHR Compiler Errors and Warnings

In this section we summarize the most important error and warning messages of the CHR compiler.

### 9.8.1 CHR Compiler Errors

**Type clash**  
for variable ... in rule ...

This error indicates an inconsistency between declared types; a variable can not belong to two types. See static type checking.

**Invalid functor**  
in head ... of rule ...

This error indicates an inconsistency between a declared type and the use of a functor in a rule. See static type checking.

**Cyclic alias**  
definition: ... == ...

You have defined a type alias in terms of itself, either directly or indirectly.

**Ambiguous type aliases**  
You have defined two overlapping type aliases.

**Multiple definitions**  
for type

You have defined the same type multiple times.

**Non-ground type**  
in constraint definition: ...

You have declared a non-ground type for a constraint argument.

**Could not find type definition**  
for ...

You have used an undefined type in a type declaration.

**Illegal mode/type declaration**  
You have used invalid syntax in a constraint declaration.

**Constraint multiply defined**  
There is more than one declaration for the same constraint.

**Undeclared constraint**  
... in head of ...

You have used an undeclared constraint in the head of a rule. This often indicates a misspelled constraint name or wrong number of arguments.

**Invalid pragma**  
... in ... Pragma should not be a variable.

You have used a variable as a pragma in a rule. This is not allowed.

**Invalid identifier**  
... in pragma passive in ...

You have used an identifier in a passive pragma that does not correspond to an identifier in the head of the rule. Likely the identifier name is misspelled.

**Unknown pragma**  
... in ...

You have used an unknown pragma in a rule. Likely the pragma is misspelled or not supported.

**Something unexpected**  
happened in the CHR compiler

You have most likely bumped into a bug in the CHR compiler. Please contact Tom Schrijvers to notify him of this error.
