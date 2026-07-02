
## 12.5 Linking embedded applications using swipl-ld

The utility program **swipl-ld** (Win32: swipl-ld.exe) may be used to link a combination of C files and Prolog files into a stand-alone executable. **swipl-ld** automates most of what is described in the previous sections.

In normal usage, a copy is made of the default embedding template `.../swipl/include/stub.c`. The **main()** routine is modified to suit your application. [PL_initialise()](foreigninclude.html#PL_initialise()) **must** be passed the program name (`argv[0]`) (Win32: the executing program can be obtained using **GetModuleFileName()**). The other elements of the command line may be modified. Next, **swipl-ld** is typically invoked as:

``` code
swipl-ld -o output stubfile.c [other-c-or-o-files] [plfiles]
```

**swipl-ld** will first split the options into various groups for both the C compiler and the Prolog compiler. Next, it will add various default options to the C compiler and call it to create an executable holding the user's C code and the Prolog kernel. Then, it will call the SWI-Prolog compiler to create a saved state from the provided Prolog files and finally, it will attach this saved state to the created emulator to create the requested executable.

Below, it is described how the options are split and which additional options are passed.

**-help**  
Print brief synopsis.

**-pl** `prolog`  
Select the Prolog to use. This Prolog is used for two purposes: get the home directory as well as the compiler/linker options and create a saved state of the Prolog code.

**-ld** `linker`  
Linker used to link the raw executable. Default is to use the C compiler (Win32: **link.exe**).

**-cc** `C compiler`  
Compiler for `.c` files found on the command line. Default is the compiler used to build SWI-Prolog accessible through the Prolog flag [c_cc](flags.html#flag:c_cc) (Win32: **cl.exe**).

**-c++** `C++-compiler`  
Compiler for C++ source file (extensions `.cpp`, `.cxx`, `.cc` or `.C`) found on the command line. Default is **c++** or **g++** if the C compiler is **gcc** (Win32: cl.exe).

**-nostate**  
Just relink the kernel, do not add any Prolog code to the new kernel. This is used to create a new kernel holding additional foreign predicates on machines that do not support the shared-library (DLL) interface, or if building the state cannot be handled by the default procedure used by **swipl-ld**. In the latter case the state is created separately and appended to the kernel using `cat <``kernel``> <``state``> > <``out``>` (Win32: `copy /b <``kernel``>+<``state``> <``out``>`).

**-shared**  
Link C, C++ or object files into a shared object (DLL) that can be loaded by the [load_foreign_library/1](foreignlink.html#load_foreign_library/1) predicate. If used with **-c** it sets the proper options to compile a C or C++ file ready for linking into a shared object.

**-dll**  
*Windows only*. Embed SWI-Prolog into a DLL rather than an executable.

**-c**  
Compile C or C++ source files into object files. This turns **swipl-ld** into a replacement for the C or C++ compiler, where proper options such as the location of the include directory are passed automatically to the compiler.

**-E**  
Invoke the C preprocessor. Used to make **swipl-ld** a replacement for the C or C++ compiler.

**-pl-options** `, ...`  
Additional options passed to Prolog when creating the saved state. The first character immediately following `pl-options` is used as separator and translated to spaces when the argument is built. Example: `-pl-options,-F,xpce` passes `-F xpce` as additional flags to Prolog.

**-ld-options** `, ...`  
Passes options to the linker, similar to **-pl-options**.

**-cc-options** `, ...`  
Passes options to the C/C++ compiler, similar to **-pl-options**.

**-v**  
Select verbose operation, showing the various programs and their options.

**-o** `outfile`  
Reserved to specify the final output file.

**-l**`library`  
Specifies a library for the C compiler. By default, `-lswipl` (Win32: libpl.lib) and the libraries needed by the Prolog kernel are given.

**-L**`library-directory`  
Specifies a library directory for the C compiler. By default the directory containing the Prolog C library for the current architecture is passed.

****-g** \| **-I`include-directory`** \| **-D`definition`****  
These options are passed to the C compiler. By default, the include directory containing `SWI-Prolog.h` is passed. **swipl-ld** adds three additional **\*** `-D`def flags:

**-D**`__SWI_PROLOG__`  
Indicates the code is to be connected to SWI-Prolog.

**-D**`__SWI_EMBEDDED__`  
Indicates the creation of an embedded program.

**-D**`_SWIPL_HOME=...`  
Provides the current SWI-Prolog home. This is used by `SWI-Prolog.h` to define the `SWIPL_HOME` macro to a string holding the home directory. This may be used to construct the commandline argument `--home=<dir>`. For example

``` code
int
main(int argc, char **argv)
{ char *av[] = {argv[0], "--home=" SWIPL_HOME, "-q", NULL};
  int ac = sizeof(av)/sizeof(*av)-1;
  if ( !PL_initialise(ac, argv) )
    PL_halt(1);
}
```

 `*.o | *.c | *.C | *.cxx | *.cpp`  
Passed as input files to the C compiler.

 `*.pl | *.qlf`  
Passed as input files to the Prolog compiler to create the saved state.

 `*`  
All other options. These are passed as linker options to the C compiler.

### 12.5.1 A simple example

The following is a very simple example going through all the steps outlined above. It provides an arithmetic expression evaluator. We will call the application **calc** and define it in the files `calc.c`^(241A similar C++ program is in [C++ interface to SWI-Prolog (Version 2)](https://www.swi-prolog.org/pldoc/man?section=cpp2).) and `calc.pl`. The Prolog file is simple:

``` code
calc(String) :-
    term_string(Expr, String),
    A is Expr,
    writeln(A).
```

The C part of the application parses the command line options, initialises the Prolog engine, locates the `calc/1` predicate and calls it. The code is in [figure 8](plld.html#fig:calc).

``` code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SWI-Prolog.h>

int
main(int argc, char **argv)
{ char *expression;
  char *program = argv[0];

  if ( argc < 2 )
  { fprintf(stderr, "Usage: %s expression\n", program);
    exit(1);
  }

  /* Determine length of joined arguments */
  size_t len = 1 + argc - 2;
  for(int n=1; n<argc; n++)
    len += strlen(argv[n]);
  if ( !(expression = malloc(len)) ) // freed on exit
  { perror("allocate");
    exit(1);
  }

  /* combine all the arguments as a single string */
  char *e = expression;
  for(int n=1; n<argc; n++)
  { if ( n != 1 )
      *e++ = ' ';
    strcpy(e, argv[n]);
    e += strlen(e);
  }

  /* make the argument vector for Prolog */

  int   plac=0;
  char *plav[2];
  plav[plac++] = program;
  plav[plac]   = NULL;

  /* initialise Prolog */

  if ( !PL_initialise(plac, plav) )
    PL_halt(1);

  /* Lookup calc/1 and make the arguments and call */

  { predicate_t pred = PL_predicate("calc", 1, "user");
    bool rc;
    term_t h0;

    rc = ( (h0=PL_new_term_refs(1)) &&
           PL_unify_chars(h0, PL_STRING|REP_MB, (size_t)-1, expression) &&
           PL_call_predicate(NULL, PL_Q_NORMAL, pred, h0) );

    PL_halt(!rc);
  }

  return 0;
}
```

**Figure 8 :** C source for the calc application

The application is now created using the command line below. The option `-goal true` sets the Prolog initialization goal to suppress the banner. Note that the `-o calc` does not specify an extension. If the platform uses a file extension for executables, **swipl-ld** will add this (e.g., `.exe` on Windows). For more details on the `swipl-ld` command, see [section 12.5](plld.html#sec:12.5).

``` code
% swipl-ld -goal true -o calc calc.c calc.pl
```

The created program **calc** is a native executable with the Prolog code attached to it. Note that the program typically depends on the shared object `libswipl` and, depending on the platform and configuration, on several external shared objects.

``` code
% ./calc pi/2
1.5708
```
