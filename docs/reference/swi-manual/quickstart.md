
## 2.1 Getting started quickly

### 2.1.1 Starting SWI-Prolog

#### 2.1.1.1 Starting SWI-Prolog on Unix

By default, SWI-Prolog is installed as **swipl**. The command line arguments of SWI-Prolog itself and its utility programs are documented using standard Unix **man** pages. SWI-Prolog is normally operated as an interactive application simply by starting the program:

``` code
$ swipl
Welcome to SWI-Prolog ...
...

1 ?-
```

After starting Prolog, one normally loads a program into it using [consult/1](consulting.html#consult/1), which may be abbreviated by putting the name of the program file between square brackets. The following goal loads the file [likes.pl](https://raw.githubusercontent.com/SWI-Prolog/swipl-devel/master/demo/likes.pl) containing clauses for the predicates likes/2 :

``` code
?- [likes].
true.

?-
```

Alternatively, the source file may be given as command line arguments:

``` code
$ swipl likes.pl
Welcome to SWI-Prolog ...
...

1 ?-
```

> Both the above assume `likes.pl` is in your *working directory*. If you use the command line version **swipl** the working directory is the same as the shell from which you started SWI-Prolog. If you started the GUI version (**swipl-win**) this depends largely on the OS. You can use pwd/0 and cd/0 to find and change the working directory. The utility ls/0 lists the contents of the working directory.
>
> ``` code
> ?- pwd.
> % /home/janw/src/swipl-devel/linux/
> true.
> ?- cd('~/tmp').
> true.
>
> ?- pwd.
> % /home/janw/tmp/
> true.
> ```
>
> The file `likes.pl` is also installed in a subdirectory `demo` insides SWI-Prolog's installation directory and may be loaded regardless of the working directory using the command below. See [absolute_file_name/3](files.html#absolute_file_name/3) and [file_search_path/2](consulting.html#file_search_path/2) for details on how SWI-Prolog specifies file locations.
>
> ``` code
> ?- [swi(demo/likes)].
> true.
> ```

After this point, Unix and Windows users unite, so if you are using Unix please continue at [section 2.1.2](quickstart.html#sec:2.1.2).

#### 2.1.1.2 Starting SWI-Prolog on Windows

After SWI-Prolog has been installed on a Windows system, the following important new things are available to the user:

- A folder (called *directory* in the remainder of this document) called `swipl` containing the executables, libraries, etc., of the system. No files are installed outside this directory.
- A program **swipl-win.exe**, providing a window for interaction with Prolog. The program **swipl.exe** is a version of SWI-Prolog that runs in a console window.
- The file extension `.pl` is associated with the program **swipl-win.exe**. Opening a `.pl` file will cause **swipl-win.exe** to start, change directory to the directory in which the file to open resides, and load this file.

The normal way to start the `likes.pl` file mentioned in [section 2.1.1.1](quickstart.html#sec:2.1.1.1) is by simply double-clicking this file in the Windows explorer.

### 2.1.2 Adding rules from the console

Although we strongly advice to put your program in a file, optionally edit it and use [make/0](consulting.html#make/0) to reload it (see [section 2.1.4](quickstart.html#sec:2.1.4)), it is possible to manage facts and rules from the terminal. The most convenient way to add a few clauses is by consulting the pseudo file `user`. The input is ended using the system end-of-file character.

``` code
?- [user].
|: hello :- format('Hello world~n').
|: ^D
true.

?- hello.
Hello world
true.
```

The predicates [assertz/1](db.html#assertz/1) and [retract/1](db.html#retract/1) are alternatives to add and remove rules and facts.

### 2.1.3 Executing a query

After loading a program, one can ask Prolog queries about the program. The query below asks Prolog what food‘sam’likes. The system responds with `X = <``value``>` if it can prove the goal for a certain `X`. The user can type the semi-colon (;) or spacebar^(7On most installations, single-character commands are executed without waiting for the RETURN key.) if (s)he wants another solution. Use the return key if you do not want to see more answers. Prolog completes the output with a full stop (.) if the user uses the return key or Prolog *knows* there are no more answers. If Prolog cannot find (more) answers, it writes **false.** Finally, Prolog answers using an error message to indicate the query or program contains an error.

``` code
?- likes(sam, X).
X = dahl ;
X = tandoori ;
...
X = chips.

?-
```

Note that the answer written by Prolog is a valid Prolog program that, when executed, produces the same set of answers as the original program.^(8The SWI-Prolog top level differs in several ways from traditional Prolog top level. The current top level was designed in cooperation with Ulrich Neumerkel.)

### 2.1.4 Examining and modifying your program

If properly configured, the predicate [edit/1](edit.html#edit/1) starts the built-in or user configured editor on the argument. The argument can be anything that can be linked to a location: a file name, predicate name, module name, etc. If the argument resolves to only one location the editor is started on this location, otherwise the user is presented a choice.

If a graphical user interface is available, the editor normally creates a new window and the system prompts for the next command. The user may edit the source file, save it and run [make/0](consulting.html#make/0) to update any modified source file. If the editor cannot be opened in a window, it opens in the same console and leaving the editor runs [make/0](consulting.html#make/0) to reload any source files that have been modified.

``` code
?- edit(likes).

true.
?- make.
% /home/jan/src/pl-devel/linux/likes compiled 0.00 sec, 0 clauses

?- likes(sam, X).
...
```

The program can also be *decompiled* using [listing/1](listing.html#listing/1) as below. The argument of [listing/1](listing.html#listing/1) is just a predicate name, a predicate *indicator* of the form `Name/Arity`, e.g., `?- listing(``mild/1).` or a *head*, e.g., `?- listing(likes(sam, _)).`, listing all *matching* clauses. The predicate [listing/0](listing.html#listing/0), i.e., without arguments lists the entire program.^(9This lists several *hook* predicates that are defined by default and is typically not very informative.)

``` code
?- listing(mild).
mild(dahl).
mild(tandoori).
mild(kurma).

true.
```

### 2.1.5 Stopping Prolog

The interactive toplevel can be stopped in two ways: enter the system end-of-file character (typically *Control-D*) or by executing the [halt/0](toplevel.html#halt/0) predicate:

``` code
?- halt.
$
```
