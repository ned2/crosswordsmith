
## 14.6 Protecting your code

Prolog in general, but SWI-Prolog in particular is an transparent environment. Prolog's “code is data” point of view makes this natural as it simplifies development and debugging. Some users though want or need to protect their code against copying or reverse engineering.

There are three ways to distribute code: as source, as `.qlf` file and in a saved state. Both QLF files and saved states contain the code as *virtual machine code*. QLF files capture the predicates and directives, while saved state capture the current state of the program. From the viewpoint of protecting code there is no significant difference.

There are two aspects to protection. One is to make sure the attacker has no access to the code in any format and the other is to provide access to a non-human-readable version of the code. The second approach is known as code obfuscation. Code obfuscation typically remove layout and comments and rename all internal identifiers. If an attacker gets access to the SWI-Prolog virtual machine code this can be *decompiled*. The decompiled code does not include layout information variable names and comments. Other identifiers, notably predicate and module names are maintained. This provides some protection against understanding the source as Prolog code without meaningful variable names and comments is generally hard to follow.

For further protecting the code, there are several scenarios.

- If the user has unrestricted access to the file system on which the application is installed the user can always access the state or QLF file. This data can be loaded into a compatible emulator and be *decompiled*.
- If the user can run arbitrary Prolog code or shell commands the state can be protected by embedding it as a string in the executable deny read access to the executable. This requires a small C program that includes the string and uses [PL_set_resource_db_mem()](foreigninclude.html#PL_set_resource_db_mem()) to register the string as the resource database. See [PL_set_resource_db_mem()](foreigninclude.html#PL_set_resource_db_mem()) for details. This protection should be combined with the [protect_static_code](flags.html#flag:protect_static_code) described below.
- Some extra protection can be provided using the Prolog flag [protect_static_code](flags.html#flag:protect_static_code), which disables decompilation of *static* predicates. Note that most Prolog implementations cannot decompile static code. Various SWI-Prolog tools depend on this ability though. Examples are [list_undefined/0](check.html#list_undefined/0), autoload/0, [show_coverage/1](prologcoverage.html#show_coverage/1), etc.

### 14.6.1 Obfuscating code in saved states

If the option `obfuscate(true)` is used with [qsave_program/2](saved-states.html#qsave_program/2), certain atoms in the saved state are renamed. The renaming is performed by library `library(obfuscate)`. The current implementation is rather conservative, renaming atoms that are used only to define the functor that names a predicate. This is a safe operation, provided the application does not create new references to renamed predicates by reading additional source code or constructing the atom that names the predicate dynamically in some other way such as using [atom_concat/3](manipatom.html#atom_concat/3). Predicates that are called this way must be declared using [public/1](dynamic.html#public/1).

Note that more aggressive renaming is possible, but this requires more detailed analysis of the various roles played by some atom. Helpful and descriptive predicate names tend to be unique and are thus subject to this transformation. More general names tend to collide with other roles of the same atom and thus prevent renaming.
