
## A.18 library(fastrw): Fast reading and writing of terms

Compatibility  
The format is not compatible to SICStus/Ciao (which are not compatible either). Future versions of this library might implement a different encoding.

bug  
The current implementation of [fast_read/1](fastrw.html#fast_read/1) **is not safe**. It is guaranteed to safely read terms written using [fast_write/1](fastrw.html#fast_write/1), but may crash on arbitrary input. The implementation does perform some basic sanity checks, including validation of the magic start byte.

To be done  
Establish a portable binary format.

This library provides the SICStus and Ciao `library(fastrw)` interface. The idea behind this library is to design a fast serialization for Prolog terms. Ideally, this should be portable between Prolog implementation. Unfortunately there is no portably binary term format defined.

The current implementation is based on PL_record_external(), which provides a binary representation of terms that is processed efficiently and can handle subterm sharing, cycles and attributed variables. In other words, this library can handle any Prolog term except *blobs* such as stream handles, database references, etc. We try to keep the format compatible between versions, but this is not guaranteed. Conversion is always possible by reading a database using the old version, dump it using [write_canonical/1](termrw.html#write_canonical/1) and read it into the new version.

This library is built upon the following built in predicates:

- [fast_term_serialized/2](IO.html#fast_term_serialized/2) translates between a term and its serialization as a byte string.
- [fast_read/2](IO.html#fast_read/2) and [fast_write/2](IO.html#fast_write/2) read/write binary serializations.

**fast_read**(`-Term`)  
The next term is read from current standard input and is unified with `Term`. The syntax of the term must agree with fast_read / fast_write format. If the end of the input has been reached, `Term` is unified with the term `end_of_file`.

**fast_write**(`+Term`)  
Output `Term` in a way that [fast_read/1](fastrw.html#fast_read/1) and [fast_read/2](IO.html#fast_read/2) will be able to read it back.

**fast_write_to_string**(`+Term, -String, ?Tail`)  
Perform a fast-write to the difference-slist `String``\``Tail`.
