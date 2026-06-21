
## 2.3 Initialisation files and goals

Using command line arguments (see [section 2.4](cmdline.html#sec:2.4)), SWI-Prolog can be forced to load files and execute queries for initialisation purposes or non-interactive operation. The most commonly used options are **-f** `file` or **-s** `file` to make Prolog load a file, **-g** `goal` to define initialisation goals and **-t** `goal` to define the *top-level goal*. The following is a typical example for starting an application directly from the command line.

``` code
machine% swipl -s load.pl -g go -t halt
```

It tells SWI-Prolog to load `load.pl`, start the application using the *entry point* go/0 and ---instead of entering the interactive top level--- exit after completing go/0 .

The command line may have multiple **-g** `goal` occurrences. The goals are executed in order. Possible choice points of individual goals are pruned. If a `goal` fails execution stops with exit status `1`. If a `goal` raises an exception, the exception is printed and the process stops with exit code `2`.

The **-q** may be used to suppress all informational messages as well as the error message that is normally printed if an initialisation goal *fails*.

In MS-Windows, the same can be achieved using a short-cut with appropriately defined command line arguments. A typically seen alternative is to write a file `run.pl` with content as illustrated below. Double-clicking `run.pl` will start the application.

``` code
:- [load].                      % load program
:- go.                          % run it
:- halt.                        % and exit
```

[Section 2.11.1.1](compilation.html#sec:2.11.1.1) discusses further scripting options, and [chapter 14](runtime.html#sec:14) discusses the generation of runtime executables. Runtime executables are a means to deliver executables that do not require the Prolog system.
