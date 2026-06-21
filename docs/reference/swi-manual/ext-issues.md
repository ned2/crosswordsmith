
## 5.7 Remaining issues

The changes and extensions described in this chapter resolve many limitations of the Prolog language we have encountered. Still, there are remaining issues for which we seek solutions in the future.

**Text representation**

Although strings resolve this issue for many applications, we are still faced with the representation of text as lists of characters which we need for parsing using DCGs. The ISO standard provides two representations, a list of *character codes* (\`codes’for short) and a list of *one-character atoms* (\`chars’for short). There are two sets of predicates, named \*\_code(s) and \*\_char(s) that provide the same functionality (e.g., [atom_codes/2](manipatom.html#atom_codes/2) and [atom_chars/2](manipatom.html#atom_chars/2)) using their own representation of characters. Codes can be used in arithmetic expressions, while chars are more readable. Neither can unambiguously be interpreted as a representation for text because codes can be interpreted as a list of integers and chars as a list of atoms.

We have not found a convincing way out. One of the options could be the introduction of a‘char’type. This type can be allowed in arithmetic and with the 0’\<`char`\> syntax we have a concrete syntax for it.

**Arrays**

Although lists are generally a much cleaner alternative for Prolog, real arrays with direct access to elements can be useful for particular tasks. The problem of integrating arrays is twofold. First of all, there is no good one-size-fits-all data representation for arrays. Many tasks that involve arrays require *mutable* arrays, while Prolog data is immutable by design. Second, standard Prolog has no good syntax support for arrays. SWI-Prolog version 7 has‘block operators’(see [section 5.3.3](ext-syntax.html#sec:5.3.3)) which can resolve the syntactic issues. Block operators have been adopted by YAP.

**Lambda expressions**

Although many alternatives^(184See e.g., [http://www.complang.tuwien.ac.at/ulrich/Prolog-inedit/ISO-Hiord](http://www.complang.tuwien.ac.at/ulrich/Prolog-inedit/ISO-Hiord)) have been proposed, we still feel uneasy with them.

**Loops**

Many people have explored routes to avoid the need for recursion in Prolog for simple iterations over data. ECLiPSe have proposed *logical loops* [Schimpf, 2002](Bibliography.html#logicalloops:2002), while B-Prolog introduced *declarative loops* and *list comprehension* [Zhou, 2010](Bibliography.html#declarativeloops:2010). The above mentioned lambda expressions, combined with [maplist/2](apply.html#maplist/2) can achieve similar results.
