JPL: A bidirectional Prolog/Java interface

Paul Singleton & Fred Dushin

Abstract

JPL is a library using the SWI-Prolog foreign interface and the Java jni interface providing a bidirectional interface between Java and Prolog that can be used to embed Prolog in Java as well as for embedding Java in Prolog. In both setups it provides a reentrant bidirectional interface.

This document is a reference for the Prolog API. The overall documentation is maintained on the [GitHub Wiki](https://github.com/ssardina-research/packages-jpl/wiki)

### 1 library(jpl): A Java interface for SWI Prolog 7.x

See also  
[http://jpl7.org/](http://jpl7.org/)

The `library(jpl)` provides a bidirectional interface to a Java Virtual Machine.

\[det\]**jpl_new**(`+X, +Params, -V`)  
`X` can be:

- an atomic classname, e.g. `'java.lang.String'`
- or an atomic descriptor, e.g. `'[I'` or `'Ljava.lang.String;'`
- or a suitable type, i.e. any `class(_,_)` or `array(_)`, e.g. `class([java,util],['Date'])`

If `X` is an object (non-array) type or descriptor and `Params` is a list of values or references, then `V` is the result of an invocation of that type's most specifically-typed constructor to whose respective formal parameters the actual `Params` are assignable (and assigned).

If `X` is an array type or descriptor and `Params` is a list of values or references, each of which is (independently) assignable to the array element type, then `V` is a new array of as many elements as `Params` has members, initialised with the respective members of `Params`.

If `X` is an array type or descriptor and `Params` is a non-negative integer N, then `V` is a new array of that type, with N elements, each initialised to Java's appropriate default value for the type.

If `V` is literally `{Term}` then we attempt to convert a `new org.jpl7.Term` instance to a corresponding term; this is of little obvious use here, but is consistent with [jpl_call/4](#jpl_call/4) and [jpl_get/3](#jpl_get/3).

\[det\]**jpl_call**(`+X, +MethodName:atom, +Params:list(datum), -Result:datum`)  
`X` should be either

- an object reference, e.g. `<jref>(1552320)` (for static or instance methods)
- or a classname, e.g. `'java.util.Date'` (for static methods only)
- or a descriptor, e.g. `'Ljava.util.Date;'` (for static methods only)
- or type, e.g. `class([java,util],['Date'])` (for static methods only)

`MethodName` should be a method name (as an atom) (may involve dynamic overload resolution based on inferred types of params)

`Params` should be a proper list (perhaps empty) of suitable actual parameters for the named method.

The class or object may have several methods with the given name; JPL will resolve (per call) to the most appropriate method based on the quantity and inferred types of `Params`. This resolution mimics the corresponding static resolution performed by Java compilers.

Finally, an attempt will be made to unify `Result` with the method's returned value, or with `@(void)` (the compound term with name `@` and argument `void`) if it has none.

\[det\]**jpl_get**(`+X, +Fspec, -V:datum`)  
`X` can be

- a classname
- or a descriptor
- or an (object or array) type (for static fields)
- or a non-array object (for static and non-static fields)
- or an array (for’length’pseudo field, or indexed element retrieval)

`Fspec` can be

- an atomic field name
- or an integral array index (to get an element from an array)
- or a pair I-J of integers (to get a subrange of an array).

Finally, an attempt will be made to unify `V` with the retrieved value or object reference.

Examples

``` code
jpl_get('java.awt.Cursor', 'NE_RESIZE_CURSOR', Q).
Q = 7.

jpl_new(array(class([java,lang],['String'])), [for,while,do,if,then,else,try,catch,finally], A),
jpl_get(A, 3-5, B).
B = [if, then, else].
```

\[det\]**jpl_set**(`+X, +Fspec, +V`)  
sets the `Fspec`-th field of (class or object) `X` to value `V` iff it is assignable

`X` can be

- a class instance (for static or non-static fields)
- or an array (for indexed element or subrange assignment)
- or a classname, or a `class(_,_)` or `array(_)` type (for static fields)
- but not a String (no fields to retrieve)

`Fspec` can be

- an atomic field name (overloading through shadowing has yet to be handled properly)
- or an array index I (`X` must be an array object: `V` is assigned to `X`\[I\])
- or a pair I-J of integers (`X` must be an array object, `V` must be a list of values: successive members of `V` are assigned to `X`\[I..J\])

`V` must be a suitable value or object.

\[det\]**jpl_get_default_jvm_opts**(`-Opts:list(atom)`)  
Returns (as a list of atoms) the options which will be passed to the JVM when it is initialised, e.g. `['-Xrs']`

\[det\]**jpl_set_default_jvm_opts**(`+Opts:list(atom)`)  
Replaces the default JVM initialisation options with those supplied.

\[semidet\]**jpl_get_actual_jvm_opts**(`-Opts:list(atom)`)  
Returns (as a list of atoms) the options with which the JVM was initialised.

Fails silently if a JVM has not yet been started, and can thus be used to test for this.

**jpl_pl_lib_version**(`-Version`)  
`Version` is the fully qualified version identifier of the in-use Prolog component (`jpl.pl`) of JPL.

It should exactly match the version identifiers of JPL's C (jpl.c) and Java (jpl.jar) components.

Example

``` code
?- jpl_pl_lib_version(V).
V = '7.6.1'.
```

**jpl_c_lib_version**(`-Version`)  
`Version` is the fully qualified version identifier of the in-use C component (jpl.c) of JPL.

It should exactly match the version identifiers of JPL's Prolog (`jpl.pl`) and Java (jpl.jar) components.

Example

``` code
?- jpl_c_lib_version(V).
V = '7.4.0-alpha'.
```

**jpl_class_to_classname**(`+Class:jref, -ClassName:entityName`)  
`Class` is a reference to a class object.

`ClassName` is its canonical (?) source-syntax (dotted) name, e.g. `'java.util.Date'`

NB not used outside jni_junk and jpl_test (is this (still) true?)

NB oughta use the available caches (but their indexing doesn't suit)

`TODO` This shouldn't exist as we have jpl_class_to_entityname/2 ???

The implementation actually just calls `Class.getName()` to get the entity name (dotted name)

**jpl_class_to_type**(`+Class:jref, -Type:jpl_type`)  
The `Class` is a reference to a (Java Universe) instance of `java.lang.Class`. The `Type` is the (Prolog Universe) JPL type term denoting the same type as does the instance of `Class`.

NB should ensure that, if not found in cache, then cache is updated.

Intriguingly, getParameterTypes returns class objects (undocumented AFAIK) with names’boolean’,’byte’etc. and even’void’(?!)

**jpl_classname_to_class**(`+EntityName:atom, -Class:jref`)  
`EntityName` is the entity name to be mapped to a class reference.

`Class` is a (canonical) reference to the corresponding class object.

NB uses caches where the class has already been mapped once before.

**jpl_entityname_to_type**(`+EntityName:atom, -Type:jpl_type`)  
`EntityName` is the entity name (an atom) denoting a Java type, to be mapped to a JPL type. This is the string returned by `java.lang.Class.getName()`.

`Type` is the JPL type (a ground term) denoting the same Java type as `EntityName` does.

The Java type in question may be a reference type (class, abstract class, interface), and array type or a primitive, including "void".

Examples:

``` code
int                       int
integer                   class([],[integer])
void                      void
char                      char
double                    double
[D                        array(double)
[[I                       array(array(int))
java.lang.String          class([java,lang],['String'])
[Ljava.lang.String;       array(class([java,lang],['String']))
[[Ljava.lang.String;      array(array(class([java, lang], ['String'])))
[[[Ljava.util.Calendar;   array(array(array(class([java,util],['Calendar']))))
foo.bar.Bling$Blong       class([foo,bar],['Bling','Blong'])
```

NB uses caches where the class has already been mapped once before.

See also  
[https://docs.oracle.com/en/java/javase/14/docs/api/java.base/java/lang/Class.html\\getName()](https://docs.oracle.com/en/java/javase/14/docs/api/java.base/java/lang/Class.html\#getName())

**jpl_type_to_entityname**(`+Type:jpl_type, -EntityName:atom`)  
This is the converse of [jpl_entityname_to_type/2](#jpl_entityname_to_type/2)

**jpl_classname_to_type**(`+EntityName:atom, -Type:jpl_type`)  
This is a wrapper around [jpl_entityname_to_type/2](#jpl_entityname_to_type/2) to keep the old exported predicate alive. The name of this predicate does not fully reflect that it actually deals in entity names instead of just class names.

Use [jpl_entityname_to_type/2](#jpl_entityname_to_type/2) in preference.

**jpl_type_to_classname**(`+Type:jpl_type, -EntityName:atom`)  
This is a wrapper around [jpl_type_to_entityname/2](#jpl_type_to_entityname/2) to keep the old exported predicate alive. The name of this predicate does not fully reflect that it actually deals in entity names instead of just class names.

Use [jpl_type_to_entityname/2](#jpl_type_to_entityname/2) in preference.

**jpl_datum_to_type**(`+Datum:datum, -Type:type`)  
`Datum` must be a JPL representation of an instance of one (or more) Java types;

`Type` is the unique most specialised type of which `Datum` denotes an instance;

NB 3 is an instance of byte, char, short, int and long, of which byte and char are the joint, overlapping most specialised types, so this relates 3 to the pseudo subtype’char_byte’;

See also  
jpl_type_to_preferred_concrete_type/2 for converting inferred types to instantiable types

**jpl_object_to_class**(`+Object:jref, -Class:jref`)  
fails silently if `Object` is not a valid reference to a Java object

`Class` is a (canonical) reference to the (canonical) class object which represents the class of `Object`

NB what's the point of caching the type if we don't look there first?

**jpl_object_to_type**(`+Object:jref, -Type:type`)  
`Object` must be a proper JPL reference to a Java object (i.e. a class or array instance, but not null, void or String).

`Type` is the JPL type of that object.

\[nondet\]**jpl_primitive_type**(`-Type:atom`)  
`Type` is an atomic JPL representation of one of Java's primitive types. N.B: `void` is not included.

``` code
?- setof(Type, jpl_primitive_type(Type), Types).
Types = [boolean, byte, char, double, float, int, long, short].
```

**jpl_ref_to_type**(`+Ref:jref, -Type:type`)  
`Ref` must be a proper JPL reference (to an object, null or void).

`Type` is its type.

**jpl_type_to_class**(`+Type:jpl_type, -Class:jref`)  
`Type` is the JPL type, a ground term designating a class or an array type.

Incomplete types are now never cached (or otherwise passed around).

jFindClass throws an exception if FCN can't be found.

**jpl_is_class**(`@Term`)  
True if `Term` is a JPL reference to an instance of `java.lang.Class`.

**jpl_is_false**(`@Term`)  
True if `Term` is `@(false)`, the JPL representation of the Java boolean value’false’.

**jpl_is_null**(`@Term`)  
True if `Term` is `@(null)`, the JPL representation of Java's’null’reference.

**jpl_is_object**(`@Term`)  
True if `Term` is a well-formed JPL object reference.

NB this checks only syntax, not whether the object exists.

**jpl_is_object_type**(`@Term`)  
True if `Term` is an object (class or array) type, not e.g. a primitive, null or void.

**jpl_is_ref**(`@Term`)  
True if `Term` is a well-formed JPL reference, either to a Java object or to Java's notional but important’null’non-object.

**jpl_is_true**(`@Term`)  
True if `Term` is `@(true)`, the JPL representation of the Java boolean value’true’.

**jpl_is_type**(`@Term`)  
True if `Term` is a well-formed JPL type structure.

**jpl_is_void**(`@Term`)  
True if `Term` is `@(void)`, the JPL representation of the pseudo Java value’void’(which is returned by [jpl_call/4](#jpl_call/4) when invoked on void methods).

NB you can try passing’void’back to Java, but it won't ever be interested.

\[semidet\]**jpl_false**(`-X:datum`)  
`X` is `@(false)`, the JPL representation of the Java boolean value’false’.

See also  
[jpl_is_false/1](#jpl_is_false/1)

\[semidet\]**jpl_null**(`-X:datum`)  
`X` is `@(null)`, the JPL representation of Java's’null’reference.

See also  
[jpl_is_null/1](#jpl_is_null/1)

\[semidet\]**jpl_true**(`-X:datum`)  
`X` is `@(true)`, the JPL representation of the Java boolean value’true’.

See also  
[jpl_is_true/1](#jpl_is_true/1)

\[semidet\]**jpl_void**(`-X:datum`)  
`X` is `@(void)`, the JPL representation of the pseudo Java value’void’.

See also  
[jpl_is_void/1](#jpl_is_void/1)

**jpl_array_to_length**(`+Array:jref, -Length:integer`)  
`Array` should be a JPL reference to a Java array of any type.

`Length` is the length of that array. This is a utility predicate, defined thus:

``` code
jpl_array_to_length(A, N) :-
    (   jpl_ref_to_type(A, array(_))
    ->  jGetArrayLength(A, N)
    ).
```

**jpl_array_to_list**(`+Array:jref, -Elements:list(datum)`)  
`Array` should be a JPL reference to a Java array of any type.

`Elements` is a Prolog list of JPL representations of the array's elements (values or references, as appropriate). This is a utility predicate, defined thus:

``` code
jpl_array_to_list(A, Es) :-
    jpl_array_to_length(A, Len),
    (   Len > 0
    ->  LoBound is 0,
        HiBound is Len-1,
        jpl_get(A, LoBound-HiBound, Es)
    ;   Es = []
    ).
```

**jpl_datums_to_array**(`+Datums:list(datum), -A:jref`)  
`A` will be a JPL reference to a new Java array, whose base type is the most specific Java type of which each member of `Datums` is (directly or indirectly) an instance.

NB this fails silently if

- `Datums` is an empty list (no base type can be inferred)
- `Datums` contains both a primitive value and an object (including array) reference (no common supertype)

**jpl_enumeration_element**(`+Enumeration:jref, -Element:datum`)  
Generates each `Element` from `Enumeration`.

- if the element is a java.lang.String then `Element` will be an atom
- if the element is null then `Element` will (oughta) be null
- otherwise I reckon it has to be an object ref

**jpl_enumeration_to_list**(`+Enumeration:jref, -Elements:list(datum)`)  
`Enumeration` should be a JPL reference to an object which implements the `Enumeration` interface.

`Elements` is a Prolog list of JPL references to the enumerated objects. This is a utility predicate, defined thus:

``` code
jpl_enumeration_to_list(Enumeration, Es) :-
    (   jpl_call(Enumeration, hasMoreElements, [], @(true))
    ->  jpl_call(Enumeration, nextElement, [], E),
        Es = [E|Es1],
        jpl_enumeration_to_list(Enumeration, Es1)
    ;   Es = []
    ).
```

\[nondet\]**jpl_hashtable_pair**(`+HashTable:jref, -KeyValuePair:pair(datum,datum)`)  
Generates Key-Value pairs from the given `HashTable`.

NB String is converted to atom but Integer is presumably returned as an object ref (i.e. as elsewhere, no auto unboxing);

NB this is anachronistic: the Map interface is preferred.

**jpl_iterator_element**(`+Iterator:jref, -Element:datum`)  
`Iterator` should be a JPL reference to an object which implements the `java.util.Iterator` interface.

`Element` is the JPL representation of the next element in the iteration. This is a utility predicate, defined thus:

``` code
jpl_iterator_element(I, E) :-
    (   jpl_call(I, hasNext, [], @(true))
    ->  (   jpl_call(I, next, [], E)
        ;   jpl_iterator_element(I, E)
        )
    ).
```

**jpl_list_to_array**(`+Datums:list(datum), -Array:jref`)  
`Datums` should be a proper Prolog list of JPL datums (values or references).

If `Datums` have a most specific common supertype, then `Array` is a JPL reference to a new Java array, whose base type is that common supertype, and whose respective elements are the Java values or objects represented by `Datums`.

\[semidet\]**jpl_terms_to_array**(`+Terms:list(term), -Array:jref`)  
`Terms` should be a proper Prolog list of arbitrary terms.

`Array` is a JPL reference to a new Java array of `org.jpl7.Term`, whose elements represent the respective members of the list.

**jpl_array_to_terms**(`+JRef:jref, -Terms:list(term)`)  
`JRef` should be a JPL reference to a Java array of org.jpl7.Term instances (or ots subtypes); `Terms` will be a list of the terms which the respective array elements represent.

\[nondet\]**jpl_map_element**(`+Map:jref, -KeyValue:pair(datum,datum)`)  
`Map` must be a JPL Reference to an object which implements the `java.util.Map` interface

This generates each Key-Value pair from the `Map`, e.g.

``` code
?- jpl_call('java.lang.System', getProperties, [], Map), jpl_map_element(Map, E).
Map = @<jref>(0x20b5c38),
E = 'java.runtime.name'-'Java(TM) SE Runtime Environment' ;
Map = @<jref>(0x20b5c38),
E = 'sun.boot.library.path'-'C:\\Program Files\\Java\\jre7\\bin'
etc.
```

This is a utility predicate, defined thus:

``` code
jpl_map_element(Map, K-V) :-
    jpl_call(Map, entrySet, [], ES),
    jpl_set_element(ES, E),
    jpl_call(E, getKey, [], K),
    jpl_call(E, getValue, [], V).
```

\[nondet\]**jpl_set_element**(`+Set:jref, -Element:datum`)  
`Set` must be a JPL reference to an object which implements the `java.util.Set` interface.

On backtracking, `Element` is bound to a JPL representation of each element of `Set`. This is a utility predicate, defined thus:

``` code
jpl_set_element(S, E) :-
    jpl_call(S, iterator, [], I),
    jpl_iterator_element(I, E).
```

**jpl_servlet_byref**(`+Config, +Request, +Response`)  
This serves the *byref* servlet demo, exemplifying one tactic for implementing a servlet in Prolog by accepting the `Request` and `Response` objects as JPL references and accessing their members via JPL as required;

See also  
[jpl_servlet_byval/3](#jpl_servlet_byval/3)

**jpl_servlet_byval**(`+MultiMap, -ContentType:atom, -Body:atom`)  
This exemplifies an alternative (to jpl_servlet_byref) tactic for implementing a servlet in Prolog; most Request fields are extracted in Java before this is called, and passed in as a multimap (a map, some of whose values are maps).

**jpl_pl_syntax**(`-Syntax:atom`)  
Unifies `Syntax` with’traditional’or’modern’according to the mode in which SWI Prolog 7.x was started

# Index

?  
[jpl_array_to_length/2](#jpl_array_to_length/2)  
[jpl_array_to_list/2](#jpl_array_to_list/2)  
[jpl_array_to_terms/2](#jpl_array_to_terms/2)  
[jpl_c_lib_version/1](#jpl_c_lib_version/1)  
[jpl_call/4](#jpl_call/4)  
[jpl_class_to_classname/2](#jpl_class_to_classname/2)  
[jpl_class_to_type/2](#jpl_class_to_type/2)  
[jpl_classname_to_class/2](#jpl_classname_to_class/2)  
[jpl_classname_to_type/2](#jpl_classname_to_type/2)  
[jpl_datum_to_type/2](#jpl_datum_to_type/2)  
[jpl_datums_to_array/2](#jpl_datums_to_array/2)  
[jpl_entityname_to_type/2](#jpl_entityname_to_type/2)  
[jpl_enumeration_element/2](#jpl_enumeration_element/2)  
[jpl_enumeration_to_list/2](#jpl_enumeration_to_list/2)  
[jpl_false/1](#jpl_false/1)  
[jpl_get/3](#jpl_get/3)  
[jpl_get_actual_jvm_opts/1](#jpl_get_actual_jvm_opts/1)  
[jpl_get_default_jvm_opts/1](#jpl_get_default_jvm_opts/1)  
[jpl_hashtable_pair/2](#jpl_hashtable_pair/2)  
[jpl_is_class/1](#jpl_is_class/1)  
[jpl_is_false/1](#jpl_is_false/1)  
[jpl_is_null/1](#jpl_is_null/1)  
[jpl_is_object/1](#jpl_is_object/1)  
[jpl_is_object_type/1](#jpl_is_object_type/1)  
[jpl_is_ref/1](#jpl_is_ref/1)  
[jpl_is_true/1](#jpl_is_true/1)  
[jpl_is_type/1](#jpl_is_type/1)  
[jpl_is_void/1](#jpl_is_void/1)  
[jpl_iterator_element/2](#jpl_iterator_element/2)  
[jpl_list_to_array/2](#jpl_list_to_array/2)  
[jpl_map_element/2](#jpl_map_element/2)  
[jpl_new/3](#jpl_new/3)  
[jpl_null/1](#jpl_null/1)  
[jpl_object_to_class/2](#jpl_object_to_class/2)  
[jpl_object_to_type/2](#jpl_object_to_type/2)  
[jpl_pl_lib_version/1](#jpl_pl_lib_version/1)  
[jpl_pl_syntax/1](#jpl_pl_syntax/1)  
[jpl_primitive_type/1](#jpl_primitive_type/1)  
[jpl_ref_to_type/2](#jpl_ref_to_type/2)  
[jpl_servlet_byref/3](#jpl_servlet_byref/3)  
[jpl_servlet_byval/3](#jpl_servlet_byval/3)  
[jpl_set/3](#jpl_set/3)  
[jpl_set_default_jvm_opts/1](#jpl_set_default_jvm_opts/1)  
[jpl_set_element/2](#jpl_set_element/2)  
[jpl_terms_to_array/2](#jpl_terms_to_array/2)  
[jpl_true/1](#jpl_true/1)  
[jpl_type_to_class/2](#jpl_type_to_class/2)  
[jpl_type_to_classname/2](#jpl_type_to_classname/2)  
[jpl_type_to_entityname/2](#jpl_type_to_entityname/2)  
[jpl_void/1](#jpl_void/1)  
