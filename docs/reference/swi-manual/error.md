
## A.16 library(error): Error generating support

author  
\- Jan Wielemaker  
- Richard O'Keefe  
- Ulrich Neumerkel

See also  
\- `library(debug)` and `library(prolog_stack)`.  
- [print_message/2](printmsg.html#print_message/2) is used to print (uncaught) error terms.

This module provides predicates to simplify error generation and checking. It's implementation is based on a discussion on the SWI-Prolog mailinglist on best practices in error handling. The utility predicate [must_be/2](error.html#must_be/2) provides simple run-time type validation. The \*\_error predicates are simple wrappers around [throw/1](exception.html#throw/1) to simplify throwing the most common ISO error terms.

**type_error**(`+ValidType, +Culprit`)  
Tell the user that `Culprit` is not of the expected `ValidType`. This error is closely related to [domain_error/2](error.html#domain_error/2) because the notion of types is not really set in stone in Prolog. We introduce the difference using a simple example.

Suppose an argument must be a non-negative integer. If the actual argument is not an integer, this is a *type_error*. If it is a negative integer, it is a *domain_error*.

Typical borderline cases are predicates accepting a compound term, e.g., `point(X,Y)`. One could argue that the basic type is a compound-term and any other compound term is a domain error. Most Prolog programmers consider each compound as a type and would consider a compound that is not `point(_,_)` a *type_error*.

**domain_error**(`+ValidDomain, +Culprit`)  
The argument is of the proper type, but has a value that is outside the supported values. See [type_error/2](error.html#type_error/2) for a more elaborate discussion of the distinction between type- and domain-errors.

**existence_error**(`+ObjectType, +Culprit`)  
`Culprit` is of the correct type and correct domain, but there is no existing (external) resource of type `ObjectType` that is represented by it.

**existence_error**(`+ObjectType, +Culprit, +Set`)  
`Culprit` is of the correct type and correct domain, but there is no existing (external) resource of type `ObjectType` that is represented by it in the provided set. The thrown exception term carries a formal term structured as follows: `existence_error(ObjectType, Culprit, Set)`

Compatibility  
This error is outside the ISO Standard.

**permission_error**(`+Operation, +PermissionType, +Culprit`)  
It is not allowed to perform `Operation` on (whatever is represented by) `Culprit` that is of the given `PermissionType` (in fact, the ISO Standard is confusing and vague about these terms’meaning).

**instantiation_error**(`+FormalSubTerm`)  
An argument is under-instantiated. I.e. it is not acceptable as it is, but if some variables are bound to appropriate values it would be acceptable.

|  |  |
|----|----|
| `FormalSubTerm` | is the term that needs (further) instantiation. Unfortunately, the ISO error does not allow for passing this term along with the error, but we pass it to this predicate for documentation purposes and to allow for future enhancement. |

**uninstantiation_error**(`+Culprit`)  
An argument is over-instantiated. This error is used for output arguments whose value cannot be known upfront. For example, the goal `open(File, read, input)` cannot succeed because the system will allocate a new unique stream handle that will never unify with `input`.

**representation_error**(`+Flag`)  
A representation error indicates a limitation of the implementation. SWI-Prolog has no such limits that are not covered by other errors, but an example of a representation error in another Prolog implementation could be an attempt to create a term with an arity higher than supported by the system.

**syntax_error**(`+Culprit`)  
A text has invalid syntax. The error is described by `Culprit`. According to the ISO Standard, `Culprit` should be an implementation-dependent atom.

To be done  
Deal with proper description of the location of the error. For short texts, we allow for Type(Text), meaning Text is not a valid Type. E.g. `syntax_error(number('1a'))` means that `1a` is not a valid number.

**resource_error**(`+Resource`)  
A goal cannot be completed due to lack of resources. According to the ISO Standard, `Resource` should be an implementation-dependent atom.

\[det\]**must_be**(`+Type, @Term`)  
True if `Term` satisfies the type constraints for `Type`. Defined types are `atom`, `atomic`, `between`, `boolean`, `callable`, `chars`, `codes`, `text`, `compound`, `constant`, `float`, `integer`, `nonneg`, `positive_integer`, `negative_integer`, `nonvar`, `number`, `oneof`, `list`, `list_or_partial_list`, `symbol`, `var`, `rational`, `encoding`, `dict` and `string`.

Most of these types are defined by an arity-1 built-in predicate of the same name. Below is a brief definition of the other types.

> |  |  |
> |----|----|
> | acyclic | Acyclic term (tree); see [acyclic_term/1](typetest.html#acyclic_term/1) |
> | any | any term |
> | `between(FloatL,FloatU)` | Number \[FloatL..FloatU\] |
> | `between(IntL,IntU)` | Integer \[IntL..IntU\] |
> | boolean | One of `true` or `false` |
> | callable | Atom or compound term |
> | char | Atom of length 1 |
> | chars | Proper list of 1-character atoms |
> | code | Representation Unicode code point (0..0x10ffff) |
> | codes | Proper list of Unicode character codes |
> | compound | compound term |
> | `compound(Term)` | Compound with same name/arity as term; checks arguments |
> | constant | Same as `atomic` |
> | cyclic | Cyclic term (rational tree); see [cyclic_term/1](typetest.html#cyclic_term/1) |
> | dict | A dictionary term; see [is_dict/1](bidicts.html#is_dict/1) |
> | encoding | Valid name for a character encoding; see [current_encoding/1](error.html#current_encoding/1) |
> | list | A (non-open) list; see [is_list/1](builtinlist.html#is_list/1) |
> | `list(Type)` | Proper list with elements of `Type` |
> | list_or_partial_list | A list or an open list (ending in a variable); see is_list_or_partial_list/1 |
> | negative_integer | Integer `<` 0 |
> | nonneg | Integer `>=` 0 |
> | `oneof(L)` | Ground term that is member of L |
> | pair | Key-Value pair. Same as `compound(any-any)` |
> | positive_integer | Integer `>` 0 |
> | proper_list | Same as list |
> | stream | A stream name or valid stream handle; see [is_stream/1](IO.html#is_stream/1) |
> | symbol | Same as `atom` |
> | text | One of `atom`, `string`, `chars` or `codes` |
> | type | `Term` is a valid type specification |

In addition, types may be composed using `TypeA,TypeB`, `TypeA;TypeB` and negated using `\Type`.

Errors  
instantiation_error if `Term` is insufficiently instantiated and `type_error(Type, Term)` if `Term` is not of `Type`.

\[semidet\]**is_of_type**(`+Type, @Term`)  
True if `Term` satisfies `Type`.

\[semidet,multifile\]**has_type**(`+Type, @Term`)  
True if `Term` satisfies `Type`.

\[nondet\]**current_encoding**(`?Name`)  
True if `Name` is the name of a supported encoding. See encoding option of e.g., [open/4](IO.html#open/4).

\[nondet\]**current_type**(`?Type, @Var, -Body`)  
True when `Type` is a currently defined type and `Var` satisfies `Type` of the body term `Body` succeeds.
