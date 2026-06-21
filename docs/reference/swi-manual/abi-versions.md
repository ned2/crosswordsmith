
## 2.21 Binary compatibility

SWI-Prolog first of all attempts to maintain *source code* compatibility between versions. Data and programs can often be represented in binary form. This touches a number of interfaces with varying degrees of compatibility. The relevant version numbers and signatures are made available by [PL_version_info()](foreigninclude.html#PL_version_info()), the **--abi-version** and the Prolog flag [abi_version](flags.html#flag:abi_version).

**Foreign extensions**  
Dynamically loadable foreign extensions have the usual dependencies on the architecture, ABI model of the (C) compiler, dynamic link library format, etc. They also depend on the backward compatibility of the PL\_\* API functions provided lib `libswipl`.

A compatible API allows distribution of foreign extensions in binary form, notably for platforms on which compilation is complicated (e.g., Windows). This compatibility is therefore high on the priority list, but must infrequently be compromised.

[PL_version_info()](foreigninclude.html#PL_version_info()): `PL_VERSION_FLI`, [abi_version](flags.html#flag:abi_version) key: `foreign_interface`

**Binary terms**  
Terms may be represented in binary format using [PL_record_external()](foreigninclude.html#PL_record_external()) and [fast_write/2](IO.html#fast_write/2). As these formats are used for storing binary terms in databases or communicate terms between Prolog processes in binary form, great care is taken to maintain compatibility.

[PL_version_info()](foreigninclude.html#PL_version_info()): `PL_VERSION_REC`, [abi_version](flags.html#flag:abi_version) key: `record`

**QLF files**  
QLF files (see [qcompile/1](consulting.html#qcompile/1)) are binary representation of Prolog file or module. They represent clauses as sequences of *virtual machine* (VM) instructions. Their compatibility relies on the QLF file format and the ABI of the VM. Some care is taken to maintain compatibility.

[PL_version_info()](foreigninclude.html#PL_version_info()): `PL_VERSION_QLF`, `PL_VERSION_QLF_LOAD` and `PL_VERSION_VM`, [abi_version](flags.html#flag:abi_version) key: `qlf`, `qlf_min_load`, `vmi`

**Saved states**  
Saved states (see **-c** and [qsave_program/2](saved-states.html#qsave_program/2)) is a zip file that contains the entire Prolog database using the same representation as QLF files. A saved state may contain additional resources, such as foreign extensions, data files, etc. In addition to the dependency concerns of QLF files, built-in and core library predicates may call *internal* foreign predicates. The interface between the public built-ins and internal foreign predicates changes frequently. Patch level releases in the *stable branch* will as much as possible maintain compatibility.

The relevant ABI version keys are the same as for QLF files with one addition: [PL_version_info()](foreigninclude.html#PL_version_info()): `PL_VERSION_BUILT_IN`, [abi_version](flags.html#flag:abi_version) key: `built_in`
