
## 5.4 Dicts: structures with named arguments

SWI-Prolog version 7 introduces dicts as an abstract object with a concrete modern syntax and functional notation for accessing members and as well as access functions defined by the user. The syntax for a dict is illustrated below. `Tag` is an atom.^(178That `Tag` can be a variable. This used to be common for jargonanonymous dicts, but is now deprecated in favour of the reserved tag `#`.) As with compound terms, there is **no** space between the tag and the opening brace. The keys are either atoms or small integers (up to [max_tagged_integer](flags.html#flag:max_tagged_integer)). The values are arbitrary Prolog terms which are parsed using the same rules as used for arguments in compound terms.

> Tag{Key1:Value1, Key2:Value2, ...}

A dict can *not* hold duplicate keys. The dict is transformed into an opaque internal representation that does *not* respect the order in which the key-value pairs appear in the input text. If a dict is written, the keys are written according to the standard order of terms (see [section 4.6.1](compare.html#sec:4.6.1)). Here are some examples, where the second example illustrates that the order is not maintained and the third illustrates an anonymous dict.

``` code
?- A = point{x:1, y:2}.
A = point{x:1, y:2}.

?- A = point{y:2, x:1}.
A = point{x:1, y:2}.

?- A = #{first_name:"Mel", last_name:"Smith"}.
A = #{first_name:"Mel", last_name:"Smith"}.
```

Dicts can be unified following the standard symmetric Prolog unification rules. As dicts use an internal canonical form, the order in which the named keys are represented is not relevant. This behaviour is illustrated by the following example.

``` code
?- point{x:1, y:2} = point{y:2, x:X}.
X = 1.
```

Two dicts unify only if they have the same tag and the same set of keys and values associated with the keys unify. As normal in Prolog, after successful unification the two dicts are indistingsable according to [==/2](compare.html#==/2) and [compare/3](compare.html#compare/3).

### 5.4.1 Dict types and tags

The [PIP](https://prolog-lang.org/ImplementersForum/PIPs.html) implementors forum is working on three PIPs related to dicts and their syntax: 0102, 0104 and 0109. These PIPs are converging to the following consensus regarding the `Tag{...}` syntax:

- The ... in `Tag{...}` is a set of `Key`:`Value`, where `Key` is an atom (and possibly (small) integer) and `Value` is any term. This is compatible to SWI-Prolog.
- The Tag is either
  - The atom `#`. This denotes a *dynamic dict*, which is what our dicts are.
  - Any other atom. This denotes a *named arguments* or *struct*. This syntax is mapped to a compound term whose arity and allowed keys is defined by a declaration.
  - A variable. This is parsed into an *attributed variable* (see [section 8.1](attvar.html#sec:8.1)).

Our aim is to support these PIPs. Next to mapping “any other atom” to a named argument structure, we will support mapping these to a dynamic dict, optionally with user-defined functions (see [section 5.4.2.1](bidicts.html#sec:5.4.2.1)). The transition is guided by the Prolog flag [var_tag](flags.html#flag:var_tag). Its current default (`dict`) is compatible with the original dict implementation of SWI-Prolog.

The use of unbound tags is deprecated. It is recommended to use `#` for anonymous dynamic dicts. The main reason for using unbound tags is in using them with the predicates [:\</2](bidicts.html#:%3C/2) and [\>:\</2](bidicts.html#%3E:%3C/2) to access a set of keys in a dict with unknown tag. In the current version the `#` tag matches any tag in these predicates. For example, to extract the `x` and `y` from a dict, we can now use

``` code
    ...,
    #{x:X, y:Y} :< Dict,
```

### 5.4.2 Functions on dicts

The infix operator dot (`op(100, yfx, .)` is used to extract values and evaluate functions on dicts. Functions are recognised if they appear in the argument of a *goal* in the source text, possibly nested in a term. The keys act as field selector, which is illustrated in this example.

``` code
?- X = point{x:1,y:2}.x.
X = 1.

?- Pt = point{x:1,y:2}, write(Pt.y).
2
Pt = point{x:1,y:2}.

?- X = point{x:1,y:2}.C.
X = 1,
C = x ;
X = 2,
C = y.
```

The compiler translates a goal that contains `.``/2` terms in its arguments into a conjunction of calls to [./3](bidicts.html#./3) defined in the `system` module. Terms functor`.`2 that appears in the head are replaced with a variable and calls to [./3](bidicts.html#./3) are inserted at the start of the body. Below are two examples, where the first extracts the `x` key from a dict and the second extends a dict containing an address with the postal code, given a find_postal_code/4 predicate.

``` code
dict_x(X, X.x).

add_postal_code(Dict, Dict.put(postal_code, Code)) :-
        find_postal_code(Dict.city,
                         Dict.street,
                         Dict.house_number,
                         Code).
```

Note that expansion of `.``/2` terms implies that such terms cannot be created by writing them explicitly in your source code. Such terms can still be created with [functor/3](manipterm.html#functor/3), [=../2](manipterm.html#=../2), [compound_name_arity/3](manipterm.html#compound_name_arity/3) and [compound_name_arguments/3](manipterm.html#compound_name_arguments/3).^(179Traditional code is unlikely to use `.``/2` terms because they were practically reserved for usage in lists. We do not provide a quoting mechanism as found in functional languages because it would only be needed to quote `.``/2` terms, such terms are rare and term manipulation provides an escape route.)

**.**(`+Dict, +Function, -Result`)  
This predicate is called to evaluate `.``/2` terms found in the arguments of a goal. This predicate evaluates the field extraction described above, raising an exception if `Function` is an atom (*key*) and `Dict` does not contain the requested key. If `Function` is a compound term, it checks for the predefined functions on dicts described in [section 5.4.2.2](bidicts.html#sec:5.4.2.2) or executes a user defined function as described in [section 5.4.2.1](bidicts.html#sec:5.4.2.1).

#### 5.4.2.1 User defined functions on dicts

The tag of a dict associates the dict to a module. If the dot notation uses a compound term, this calls the goal below.

> \<`module`\>:\<`name`\>(Arg1, ..., +Dict, -Value)

Functions are normal Prolog predicates. The dict infrastructure provides a more convenient syntax for representing the head of such predicates without worrying about the argument calling conventions. The code below defines a function `multiply(Times)` on a point that creates a new point by multiplying both coordinates. and `len`^(180as `length` would result in a predicate [length/2](builtinlist.html#length/2), this name cannot be used. This might change in future versions.) to compute the length from the origin. The . and `:=` operators are used to abstract the location of the predicate arguments. It is allowed to define multiple a function with multiple clauses, providing overloading and non-determinism.

``` code
:- module(point, []).

M.multiply(F) := point{x:X, y:Y} :-
        X is M.x*F,
        Y is M.y*F.

M.len() := Len :-
        Len is sqrt(M.x**2 + M.y**2).
```

After these definitions, we can evaluate the following functions:

``` code
?- X = point{x:1, y:2}.multiply(2).
X = point{x:2, y:4}.

?- X = point{x:1, y:2}.multiply(2).len().
X = 4.47213595499958.
```

#### 5.4.2.2 Predefined functions on dicts

Dicts currently define the following reserved functions:

**get**(`?KeyPath`)  
Return the value associates with `KeyPath`. `KeyPath` is either a single key or a term `Key1/Key2/...`. Each key is either an atom, small integer or a variable. While `Dict.Key` throws an existence error, this function *fails* silently if a key does not exist in the target dict. See also [:\</2](bidicts.html#:%3C/2), which can be used to test for existence and unify multiple key values from a dict. For example:

``` code
?- write(t{a:x}.get(a)).
x
?- write(t{a:x}.get(b)).
false.
?- write(t{a:t{b:x}}.get(a/b)).
x
```

**put**(`+New`)  
Evaluates to a new dict where the key-values in `New` replace or extend the key-values in the original dict. See [put_dict/3](bidicts.html#put_dict/3).

**get**(`?KeyPath, +Default`)  
Same as get/1 , but if no match is found the function evaluates to `Default`. If `KeyPath` contains variables possible choice points are respected and the function only evaluates to `Default` if the pattern has no matches.

**put**(`+KeyPath, +Value`)  
Evaluates to a new dict where the `KeyPath`-`Value` replaces or extends the key-values in the original dict. `KeyPath` is either a key or a term `KeyPath`/`Key`,^(181Note that we do not use the’.’functor here, because the `.``/2` would *evaluate*.) replacing the value associated with `Key` in a sub-dict of the dict on which the function operates. See [put_dict/4](bidicts.html#put_dict/4). Below are some examples:

``` code
?- A = #{}.put(a, 1).
A = #{a:1}.

?- A = #{a:1}.put(a, 2).
A = #{a:2}.

?- A = #{a:1}.put(b/c, 2).
A = #{a:1, b: #{c:2}}.

?- A = #{a: #{b:1}}.put(a/b, 2).
A = #{a: #{b:2}}.

?- A = #{a:1}.put(a/b, 2).
A = #{a: #{b:2}}.
```

### 5.4.3 Predicates for managing dicts

This section documents the predicates that are defined on dicts. We use the naming and argument conventions of the traditional `library(assoc)`.

**is_dict**(`@Term`)  
True if `Term` is a dict. This is the same as `is_dict(Term,_)`.

**is_dict**(`@Term, -Tag`)  
True if `Term` is a dict of `Tag`.

**get_dict**(`?Key, +Dict, -Value`)  
Unify the value associated with `Key` in dict with `Value`. If `Key` is unbound, all associations in `Dict` are returned on backtracking. The order in which the associations are returned is undefined. This predicate is normally accessed using the functional notation `Dict.Key`. See [section 5.4.2](bidicts.html#sec:5.4.2).

Fails silently if Key does not appear in Dict. This is different from the behavior of the functional‘.\`-notation, which throws an existence error in that case.

\[semidet\]**get_dict**(`+Key, +Dict, -Value, -NewDict, +NewValue`)  
Create a new dict after updating the value for `Key`. Fails if `Value` does not unify with the current value associated with `Key`. `Dict` is either a dict or a list the can be converted into a dict.

Has the behavior as if defined in the following way:

``` code
get_dict(Key, Dict, Value, NewDict, NewValue) :-
        get_dict(Key, Dict, Value),
        put_dict(Key, Dict, NewValue, NewDict).
```

**dict_create**(`-Dict, +Tag, +Data`)  
Create a dict in `Tag` from `Data`. `Data` is a list of attribute-value pairs using the syntax `Key:Value`, `Key=Value`, `Key-Value` or `Key(Value)`. An exception is raised if `Data` is not a proper list, one of the elements is not of the shape above, a key is neither an atom nor a small integer or there is a duplicate key.

**dict_pairs**(`?Dict, ?Tag, ?Pairs`)  
Bi-directional mapping between a dict and an ordered list of pairs (see [section A.35](pairs.html#sec:A.35)).

**dict_same_keys**(`?Dict1, ?Dict2`)  
True when `Dict1` and `Dict2` have the same set of keys. The *tag* is not considered. This predicate is semidet if both arguments are instantiated to a dict. If one is instantiated to a dict and the other is unbound, it generates a dict with the same keys and unbound tag and values. The predicate [is_dict/2](bidicts.html#is_dict/2) may be used test tag equivalence or unify the tags. This predicate raises an `instantiation_error` if both argument are unbound or a `type_error` if one of the arguments is neither a dict nor a variable.

**put_dict**(`+New, +DictIn, -DictOut`)  
`DictOut` is a new dict created by replacing or adding key-value pairs from `New` to `Dict`. `New` is either a dict or a valid input for [dict_create/3](bidicts.html#dict_create/3). This predicate is normally accessed using the functional notation. Below are some examples:

``` code
?- A = point{x:1, y:2}.put(#{x:3}).
A = point{x:3, y:2}.

?- A = point{x:1, y:2}.put([x=3]).
A = point{x:3, y:2}.

?- A = point{x:1, y:2}.put([x=3,z=0]).
A = point{x:3, y:2, z:0}.
```

**put_dict**(`+Key, +DictIn, +Value, -DictOut`)  
`DictOut` is a new dict created by replacing or adding `Key`-`Value` to `DictIn`. For example:

``` code
?- A = point{x:1, y:2}.put(x, 3).
A = point{x:3, y:2}.
```

This predicate can also be accessed by using the functional notation, in which case Key can also be a \*path\* of keys. For example:

``` code
?- Dict = #{}.put(a/b, c).
Dict = #{a: #{b:c}}.
```

**del_dict**(`+Key, +DictIn, ?Value, -DictOut`)  
True when `Key`-`Value` is in `DictIn` and `DictOut` contains all associations of `DictIn` except for `Key`.

\[semidet\]`+Select` **:\<** `+From`  
True when `Select` is a‘sub dict’of `From`: the tags must match and all keys in `Select` must appear with unifying values in `From`. Two tags match if the unify or one of the tags is the reserved `#` tag. `From` may contain keys that are not in `Select`. This operation is frequently used to *match* a dict and at the same time extract relevant values from it. For example:

``` code
plot(Dict, On),
  #{x:X, y:Y, z:Z} :< Dict =>
    plot_xyz(X, Y, Z, On).
plot(Dict, On) :-
  #{x:X, y:Y} :< Dict =>
    plot_xy(X, Y, On).
```

The goal `Select :< From` is equivalent to `select_dict(Select, From, _)`.

\[semidet\]**select_dict**(`+Select, +From, -Rest`)  
True when the tags of `Select` and `From` match, all keys in `Select` appear in `From` and the corresponding values have been unified. The key-value pairs of `From` that do not appear in `Select` are used to form an anonymous dict, which is unified with `Rest`. For example:

``` code
?- select_dict(#{x:0, y:Y}, point{x:0, y:1, z:2}, R).
Y = 1,
R = #{z:2}.
```

See also [:\</2](bidicts.html#:%3C/2) to ignore `Rest` and [\>:\</2](bidicts.html#%3E:%3C/2) for a symmetric partial unification of two dicts.

`+Dict1` **\>:\<** `+Dict2`  
This operator specifies a *partial unification* between `Dict1` and `Dict2`. It is true when the tags match and the values associated with all *common* keys have been unified. The values associated to keys that do not appear in the other dict are ignored. Partial unification is symmetric. For example, given a list of dicts, find dicts that represent a point with X equal to zero:

``` code
    member(Dict, List),
    Dict >:< point{x:0, y:Y}.
```

See also [:\</2](bidicts.html#:%3C/2) and [select_dict/3](bidicts.html#select_dict/3).

#### 5.4.3.1 Destructive assignment in dicts

This section describes the destructive update operations defined on dicts. These actions can only *update* keys and not add or remove keys. If the requested key does not exist the predicate raises `existence_error(key, Key, Dict)`. Note the additional argument.

Destructive assignment is a non-logical operation and should be used with care because the system may copy or share identical Prolog terms at any time. Some of this behaviour can be avoided by adding an additional unbound value to the dict. This prevents unwanted sharing and ensures that [copy_term/2](manipterm.html#copy_term/2) actually copies the dict. This pitfall is demonstrated in the example below:

``` code
?- A = a{a:1}, copy_term(A,B), b_set_dict(a, A, 2).
A = B, B = a{a:2}.

?- A = a{a:1,dummy:_}, copy_term(A,B), b_set_dict(a, A, 2).
A = a{a:2, dummy:_G3195},
B = a{a:1, dummy:_G3391}.
```

\[det\]**b_set_dict**(`+Key, !Dict, +Value`)  
Destructively update the value associated with `Key` in `Dict` to `Value`. The update is trailed and undone on backtracking. This predicate raises an existence error if `Key` does not appear in `Dict`. The update semantics are equivalent to [setarg/3](manipterm.html#setarg/3) and [b_setval/2](gvar.html#b_setval/2).

\[det\]**nb_set_dict**(`+Key, !Dict, +Value`)  
Destructively update the value associated with `Key` in `Dict` to a copy of `Value`. The update is *not* undone on backtracking. This predicate raises an existence error if `Key` does not appear in `Dict`. The update semantics are equivalent to [nb_setarg/3](manipterm.html#nb_setarg/3) and [nb_setval/2](gvar.html#nb_setval/2).

\[det\]**nb_link_dict**(`+Key, !Dict, +Value`)  
Destructively update the value associated with `Key` in `Dict` to `Value`. The update is *not* undone on backtracking. This predicate raises an existence error if `Key` does not appear in `Dict`. The update semantics are equivalent to [nb_linkarg/3](manipterm.html#nb_linkarg/3) and [nb_linkval/2](gvar.html#nb_linkval/2). Use with extreme care and consult the documentation of [nb_linkval/2](gvar.html#nb_linkval/2) before use.

### 5.4.4 When to use dicts?

Dicts are a new type in the Prolog world. They compete with several other types and libraries. In the list below we have a closer look at these relations. We will see that dicts are first of all a good replacement for compound terms with a high or not clearly fixed arity, library `library(record)` and option processing.

**Compound terms**  
Compound terms with positional arguments form the traditional way to package data in Prolog. This representation is well understood, fast and compound terms are stored efficiently. Compound terms are still the representation of choice, provided that the number of arguments is low and fixed or compactness or performance are of utmost importance.

A good example of a compound term is the representation of RDF triples using the term `rdf(Subject, Predicate, Object)` because RDF triples are defined to have precisely these three arguments and they are always referred to in this order. An application processing information about persons should probably use dicts because the information that is related to a person is not so fixed. Typically we see first and last name. But there may also be title, middle name, gender, date of birth, etc. The number of arguments becomes unmanageable when using a compound term, while adding or removing an argument leads to many changes in the program.

**Library `library(record)`**  
Using library `library(record)` relieves the maintenance issues associated with using compound terms significantly. The library generates access and modification predicates for each field in a compound term from a declaration. The library provides sound access to compound terms with many arguments. One of its problems is the verbose syntax needed to access or modify fields which results from long names for the generated predicates and the restriction that each field needs to be extracted with a separate goal. Consider the example below, where the first uses library `library(record)` and the second uses dicts.

``` code
    ...,
    person_first_name(P, FirstName),
    person_last_name(P, LastName),
    format('Dear ~w ~w,~n~n', [FirstName, LastName]).

    ...,
    format('Dear ~w ~w,~n~n', [Dict.first_name, Dict.last_name]).
```

Records have a fixed number of arguments and (non-)existence of an argument must be represented using a value that is outside the normal domain. This lead to unnatural code. For example, suppose our person also has a title. If we know the first name we use this and else we use the title. The code samples below illustrate this.

``` code
salutation(P) :-
    person_first_name(P, FirstName), nonvar(FirstName), !,
    person_last_name(P, LastName),
    format('Dear ~w ~w,~n~n', [FirstName, LastName]).
salutation(P) :-
    person_title(P, Title), nonvar(Title), !,
    person_last_name(P, LastName),
    format('Dear ~w ~w,~n~n', [Title, LastName]).

salutation(P) :-
    #{first_name:FirstName, last_name:LastName} :< P, !,
    format('Dear ~w ~w,~n~n', [FirstName, LastName]).
salutation(P) :-
    #{title:Title, last_name:LastName} :< P, !,
    format('Dear ~w ~w,~n~n', [Title, LastName]).
```

**Library `library(assoc)`**  
This library implements a balanced binary tree. Dicts can replace the use of this library if the association is fairly static (i.e., there are few update operations), all keys are atoms or (small) integers and the code does not rely on ordered operations.

**Library `library(option)`**  
Option lists are introduced by ISO Prolog, for example for [read_term/3](termrw.html#read_term/3), [open/4](IO.html#open/4), etc. The `library(option)` library provides operations to extract options, merge options lists, etc. Dicts are well suited to replace option lists because they are cheaper, can be processed faster and have a more natural syntax.

**Library `library(pairs)`**  
This library is commonly used to process large name-value associations. In many cases this concerns short-lived data structures that result from [findall/3](allsolutions.html#findall/3), [maplist/3](apply.html#maplist/3) and similar list processing predicates. Dicts may play a role if frequent random key lookups are needed on the resulting association. For example, the skeleton‘create a pairs list’,‘use [list_to_assoc/2](assoc.html#list_to_assoc/2) to create an assoc’, followed by frequent usage of [get_assoc/3](assoc.html#get_assoc/3) to extract key values can be replaced using [dict_pairs/3](bidicts.html#dict_pairs/3) and the dict access functions. Using dicts in this scenario is more efficient and provides a more pleasant access syntax.

### 5.4.5 A motivation for dicts as primary citizens

Dicts, or key-value associations, are a common data structure. A good old example are *property lists* as found in Lisp, while a good recent example is formed by JavaScript *objects*. Traditional Prolog does not offer native property lists. As a result, people are using a wide range of data structures for key-value associations:

- Using compound terms and positional arguments, e.g., `point(1,2)`.
- Using compound terms with library `library(record)`, which generates access predicates for a term using positional arguments from a description.
- Using lists of terms `Name=Value`, `Name-Value`, `Name:Value` or `Name(Value)`.
- Using library `library(assoc)` which represents the associations as a balanced binary tree.

This situation is unfortunate. Each of these have their advantages and disadvantages. E.g., compound terms are compact and fast, but inflexible and using positional arguments quickly breaks down. Library `library(record)` fixes this, but the syntax is considered hard to use. Lists are flexible, but expensive and the alternative key-value representations that are used complicate the matter even more. Library `library(assoc)` allows for efficient manipulation of changing associations, but the syntactical representation of an assoc is complex, which makes them unsuitable for e.g., *options lists* as seen in predicates such as [open/4](IO.html#open/4).

### 5.4.6 Implementation notes about dicts

Although dicts are designed as an abstract data type and we deliberately reserve the possibility to change the representation and even use multiple representations, this section describes the current implementation.

Dicts are currently represented as a compound term using the functor `` `dict` ``. The first argument is the tag. The remaining arguments create an array of sorted key-value pairs. This representation is compact and guarantees good locality. Lookup is order `log( N )`, while adding values, deleting values and merging with other dicts has order `N`. The main disadvantage is that changing values in large dicts is costly, both in terms of memory and time.

Future versions may share keys in a separate structure or use a binary trees to allow for cheaper updates. One of the issues is that the representation must either be kept canonical or unification must be extended to compensate for alternate representations.
