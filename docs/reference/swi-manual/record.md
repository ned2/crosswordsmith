
## A.50 library(record): Access named fields in a term

The library `library(record)` provides named access to fields in a record represented as a compound term such as `point(X, Y)`. The Prolog world knows various approaches to solve this problem, unfortunately with no consensus. The approach taken by this library is proposed by Richard O'Keefe on the SWI-Prolog mailinglist.

The approach automates a technique commonly described in Prolog text-books, where access and modification predicates are defined for the record type. Such predicates are subject to normal import/export as well as analysis by cross-referencers. Given the simple nature of the access predicates, an optimizing compiler can easily inline them for optimal performance.

A record is defined using the directive [record/1](record.html#record/1). We introduce the library with a short example:

``` code
:- record point(x:integer=0, y:integer=0).

        ...,
        default_point(Point),
        point_x(Point, X),
        set_x_of_point(10, Point, Point1),

        make_point([y(20)], YPoint),
```

The principal functor and arity of the term used defines the name and arity of the compound used as records. Each argument is described using a term of the format below.

> \<`name`\>\[:\<`type`\>\]\[=\<`default`\>\]

In this definition, \<`name`\> is an atom defining the name of the argument, \<`type`\> is an optional type specification as defined by [must_be/2](error.html#must_be/2) from library `library(error)`, and \<`default`\> is the default initial value. The \<`type`\> defaults to `any`. If no default value is specified the default is an unbound variable.

A record declaration creates a set of predicates through *term-expansion*. We describe these predicates below. In this description, \<`constructor`\> refers to the name of the record (\`point’in the example above) and \<`name`\> to the name of an argument (field).

- *default\_\<`constructor`\>(-Record)*  
  Create a new record where all fields have their default values. This is the same as make\_\<`constructor`\>(\[\], Record) .
- *make\_\<`constructor`\>(+Fields, -Record)*  
  Create a new record where specified fields have the specified values and remaining fields have their default value. Each field is specified as a term \<`name`\>(\<`value`\>). See example in the introduction.
- *make\_\<`constructor`\>(+Fields, -Record, -RestFields)*  
  Same as make\_\<`constructor`\>/2, but named fields that do not appear in `Record` are returned in `RestFields`. This predicate is motivated by option-list processing. See library `library(option)`.
- *\<`constructor`\>\_\<`name`\>(Record, Value)*  
  Unify `Value` with argument in `Record` named \<`name`\>.^(253Note this is not called‘get\_’as it performs unification and can perfectly well instantiate the argument.)
- *\<`constructor`\>\_data(?Name, +Record, ?Value)*  
  True when `Value` is the value for the field named `Name` in `Record`. This predicate does not perform type-checking.
- *set\_\<`name`\>\_of\_\<`constructor`\>(+Value, +OldRecord, -NewRecord)*  
  Replace the value for \<`name`\> in `OldRecord` by `Value` and unify the result with `NewRecord`.
- *set\_\<`name`\>\_of\_\<`constructor`\>(+Value, !Record)*  
  Destructively replace the argument \<`name`\> in `Record` by `Value` based on [setarg/3](manipterm.html#setarg/3). Use with care.
- *nb_set\_\<`name`\>\_of\_\<`constructor`\>(+Value, !Record)*  
  As above, but using non-backtrackable assignment based on [nb_setarg/3](manipterm.html#nb_setarg/3). Use with *extreme* care.
- *set\_\<`constructor`\>\_fields(+Fields, +Record0, -Record)*  
  Set multiple fields using the same syntax as make\_\<`constructor`\>/2, but starting with `Record0` rather than the default record.
- *set\_\<`constructor`\>\_fields(+Fields, +Record0, -Record, -RestFields)*  
  Similar to set\_\<`constructor`\>\_fields/4, but fields not defined by \<`constructor`\> are returned in `RestFields`.
- *set\_\<`constructor`\>\_field(+Field, +Record0, -Record)*  
  Set a single field specified as a term \<`name`\>(\<`value`\>).

**record**(`+Spec`)  
The construct `:- record Spec, ...` is used to define access to named fields in a compound. It is subject to term-expansion (see [expand_term/2](consulting.html#expand_term/2)) and cannot be called as a predicate. See [section A.50](record.html#sec:A.50) for details.
