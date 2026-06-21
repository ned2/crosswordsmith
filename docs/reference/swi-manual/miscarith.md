
## 4.28 Misc arithmetic support predicates

**set_random**(`+Option`)  
Controls the random number generator accessible through the *functions* [random/1](arith.html#f-random/1) and [random_float/0](arith.html#f-random_float/0). Note that the library `library(random)` provides an alternative API to the same random primitives.

**seed**(`+Seed`)  
Set the seed of the random generator for this thread. `Seed` is an integer or the atom `random`. If `random`, repeat the initialization procedure described with the function [random/1](arith.html#f-random/1). Here is an example:

``` code
?- set_random(seed(111)), A is random(6).
A = 5.
?- set_random(seed(111)), A is random(6).
A = 5.
```

**state**(`+State`)  
Set the generator to a state fetched using the state property of [random_property/1](miscarith.html#random_property/1). Using other values may lead to undefined behaviour.^(135The limitations of the underlying (GMP) library are unknown, which makes it impossible to validate the `State`.)

**random_property**(`?Option`)  
True when `Option` is a current property of the random generator. Currently, this predicate provides access to the state. This predicate is not present on systems where the state is inaccessible.

**state**(`-State`)  
Describes the current state of the random generator. State is a normal Prolog term that can be asserted or written to a file. Applications should make no other assumptions about its representation. The only meaningful operation is to use as argument to [set_random/1](miscarith.html#set_random/1) using the `state(State)` option.^(bugGMP provides no portable mechanism to fetch and restore the state. The current implementation works, but the state depends on the platform. I.e., it is generally not possible to reuse the state with another version of GMP or on a CPU with different datasizes or endian-ness.)

**current_arithmetic_function**(`?Head`)  
True when `Head` is an evaluable function. For example:

``` code
?- current_arithmetic_function(sin(_)).
true.
```
