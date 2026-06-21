Prolog Unit Tests

Jan Wielemaker  
University of Amsterdam  
VU University Amsterdam  
The Netherlands  
E-mail: [jan@swi-prolog.org](mailto:jan@swi-prolog.org)

Abstract

This document describes a Prolog unit-test framework. This framework was initially developed for [SWI-Prolog](http://www.swi-prolog.org). The current version also runs on [SICStus Prolog](http://www.sics.se/sicstus/), providing a portable testing framework. See [section 8.1](#sec:8.1).

# Table of Contents

[1 Introduction](#sec:1)

[2 A Unit Test box](#sec:2)

[2.1 Test Unit options](#sec:2.1)

[2.2 Writing the test body](#sec:2.2)

[2.2.1 Testing deterministic predicates](#sec:2.2.1)

[2.2.2 Testing semi-deterministic predicates](#sec:2.2.2)

[2.2.3 Testing non-deterministic predicates](#sec:2.2.3)

[2.2.4 Testing error conditions](#sec:2.2.4)

[2.2.5 One body with multiple tests using assertions](#sec:2.2.5)

[3 Using separate test files](#sec:3)

[4 Running the test-suite](#sec:4)

[4.1 Running the test suite from Prolog](#sec:4.1)

[4.2 Running the test suite from the command line](#sec:4.2)

[5 Tests and production systems](#sec:5)

[6 Controlling the test suite](#sec:6)

[7 Auto-generating tests](#sec:7)

[8 Portability of the test-suite](#sec:8)

[8.1 PlUnit on SICStus](#sec:8.1)

[9 Motivation of choices](#sec:9)

## 1 Introduction

There is really no excuse not to write tests!

Automatic testing of software during development is probably the most important Quality Assurance measure. Tests can validate the final system, which is nice for your users. However, most (Prolog) developers forget that it is not just a burden during development.

- Tests document how the code is supposed to be used.
- Tests can validate claims you make on the Prolog implementation. Writing a test makes the claim explicit.
- Tests avoid big applications saying‘No’after modifications. This saves time during development, and it saves *a lot* of time if you must return to the application a few years later or you must modify and debug someone else's application.

## 2 A Unit Test box

Tests are written in pure Prolog and enclosed within the directives [begin_tests/1](#begin_tests/1),2 and end_tests/1. They can be embedded inside a normal source module, or be placed in a separate test-file that loads the files to be tested. Code inside a test box is normal Prolog code. The entry points are defined by rules using the head `test(Name)` or `test(Name, Options)`, where `Name` is a ground term and `Options` is a list describing additional properties of the test. Here is a very simple example:

``` code
:- begin_tests(lists).
:- use_module(library(lists)).

test(reverse) :-
        reverse([a,b], [b,a]).

:- end_tests(lists).
```

The optional second argument of the test-head defines additional processing options. Defined options are:

**blocked**(`+Reason:atom`)  
The test is currently disabled. Tests are flagged as blocked if they cannot be run for some reason. E.g. they crash Prolog, they rely on some service that is not available, they take too much resources, etc. Tests that fail but do not crash, etc. should be flagged using `fixme(Fixme)`.

**fixme**(`+Reason:atom`)  
Similar to `blocked(Reason)`, but the test it executed anyway. If it fails, a `-` is printed instead of the `.` character. If it passes a `+` and if it passes with a choicepoint, `!`. A summary is printed at the end of the test run and the goal `test_report(fixme)` can be used to get details.

**condition**(`:Goal`)  
Pre-condition for running the test. If the condition fails the test is skipped. The condition can be used as an alternative to the `setup` option. The only difference is that failure of a condition skips the test and is considered an error when using the `setup` option.

**cleanup**(`:Goal`)  
`Goal` is always called after completion of the test-body, regardless of whether it fails, succeeds or throws an exception. This option or call_cleanup/2 must be used by tests that require side-effects that must be reverted after the test completes. `Goal` may share variables with the test body.

``` code
create_file(Tmp) :-
        tmp_file(plunit, Tmp),
        open(Tmp, write, Out),
        write(Out, 'hello(World).\n'),
        close(Out).

test(read, [ setup(create_file(Tmp)),
             cleanup(delete_file(Tmp))
           ]) :-
        read_file_to_terms(Tmp, Terms, []),
        Term = hello(_).
```

**setup**(`:Goal`)  
`Goal` is run before the test-body. Typically used together with the `cleanup` option to create and destroy the required execution environment.

**forall**(`:Generator`)  
Run the same test for each solution of `Generator`. Each run invokes the setup and cleanup handlers. This can be used to run the same test with different inputs. If an error occurs, the test is reported as `name (forall bindings =` \<`vars`\> `)`, where \<`vars`\> indicates the bindings of variables in `Generator`.

**true**(`AnswerTerm Cmp Value`)  
Body should succeed deterministically. If a choicepoint is left open, a warning is printed to STDERR ("Test succeeded with choicepoint"). That warning can be suppressed by adding the `nondet` keyword. `AnswerTerm` is compared to `Value` using the comparison operator `Cmp`. `Cmp` is typically one of =/2, ==/2, =:=/2 or =@=/2,^(1The =@= predicate (denoted *structural equivalence*) is the same as variant/2 in SICStus.) but any test can be used. This is the same as inserting the test at the end of the conjunction, but it allows the test engine to distinguish between failure of copy_term/2 and producing the wrong value. Multiple variables must be combined in an arbitrary compound term. E.g. `A1-A2 == v1-v2`

``` code
test(copy, [ true(Copy =@= hello(X,X))
           ]) :-
        copy_term(hello(Y,Y), Copy).
```

**AnswerTerm Cmp Value**(`E`)  
quivalent to `true(AnswerTerm Cmp Value)` if `Cmp` is one of the comparison operators given above.

**fail**  
Body must fail.

**throws**(`Error`)  
Body must throw `Error`. The thrown error term is matched against term `Error` using `subsumes_term(Error, ThrownError)`. I.e., the thrown error must be more specific than the specified `Error`. See subsumes_term/2.

**error**(`Error`)  
Body must throw `error(Error, _Context)`. See keyword `throws` (as well as predicate throw/1 and library(error)) for details.

**all**(`AnswerTerm Cmp Instances`)  
Similar to `true(AnswerTerm Cmp Values)`, but used for non-deterministic predicates. Each element is compared using `Cmp`. Order matters. For example:

``` code
test(or, all(X == [1,2])) :-
        ( X = 1 ; X = 2 ).
```

**set**(`AnswerTerm Cmp Instances`)  
Similar to `all(AnswerTerm Cmp Instances)`, but before testing both the bindings of `AnswerTerm` and `Instances` are sorted using sort/2. This removes duplicates and places both sets in the same order.^(2The result is only well-defined of `Cmp` is `==`.)

**nondet**  
If this keyword appears in the option list, non-deterministic success of the body is not considered an error.

**occurs_check**(`Mode`)  
Run the test in a particular *occurs check mode*. Mode is one of `false` (default), `true` or `error`. See the Prolog flag **occurs_check** for details.

### 2.1 Test Unit options

**begin_tests**(`+Name`)  
Start named test-unit. Same as `begin_tests(Name, [])`.

**begin_tests**(`+Name, +Options`)  
Start named test-unit with options. Options provide conditional processing, setup and cleanup similar to individual tests (second argument of test/2 rules).

Defined options are:

**blocked**(`+Reason`)  
Test-unit has been blocked for the given `Reason`.

**condition**(`:Goal`)  
Executed before executing any of the tests. If `Goal` fails, the test of this unit is skipped.

**setup**(`:Goal`)  
Executed before executing any of the tests.

**cleanup**(`:Goal`)  
Executed after completion of all tests in the unit.

**occurs_check**(`+Mode`)  
Specify default for subject-to-occurs-check mode. See [section 2](#sec:2) for details on the `occurs_check` option.

### 2.2 Writing the test body

The test-body is ordinary Prolog code. Without any options, the body must be designed to succeed *deterministically*. Any other result is considered a failure. One of the options `fail`, `true`, `throws`, `all` or `set` can be used to specify a different expected result. See [section 2](#sec:2) for details. In this section we illustrate typical test-scenarios by testing SWI-Prolog built-in and library predicates.

#### 2.2.1 Testing deterministic predicates

Deterministic predicates are predicates that must succeed exactly once and, for well behaved predicates, leave no choicepoints. Typically they have zero or more input- and zero or more output arguments. The test goal supplies proper values for the input arguments and verifies the output arguments. Verification can use test-options or be explicit in the body. The tests in the example below are equivalent.

``` code
test(add) :-
        A is 1 + 2,
        A =:= 3.

test(add, [true(A =:= 3)]) :-
        A is 1 + 2.
```

The test engine verifies that the test-body does not leave a choicepoint. We illustrate that using the test below:

``` code
test(member) :-
        member(b, [a,b,c]).
```

Although this test succeeds, member/2 leaves a choicepoint which is reported by the test subsystem. To make the test silent, use one of the alternatives below.

``` code
test(member) :-
        member(b, [a,b,c]), !.

test(member, [nondet]) :-
        member(b, [a,b,c]).
```

#### 2.2.2 Testing semi-deterministic predicates

Semi-deterministic predicates are predicates that either fail or succeed exactly once and, for well behaved predicates, leave no choicepoints. Testing such predicates is the same as testing deterministic predicates. Negative tests must be specified using the option `fail` or by negating the body using `\+/1`.

``` code
test(is_set) :-
        \+ is_set([a,a]).

test(is_set, [fail]) :-
        is_set([a,a]).
```

#### 2.2.3 Testing non-deterministic predicates

Non-deterministic predicates succeed zero or more times. Their results are tested either using findall/3 or setof/3 followed by a value-check or using the `all` or `set` options. The following are equivalent tests:

``` code
test(member) :-
        findall(X, member(X, [a,b,c]), Xs),
        Xs == [a,b,c].

test(member, all(X == [a,b,c])) :-
        member(X, [a,b,c]).
```

#### 2.2.4 Testing error conditions

Error-conditions are tested using the option `throws(Error)` or by wrapping the test in a catch/3. The following tests are equivalent:

``` code
test(div0) :-
     catch(A is 1/0, error(E, _), true),
     E =@= evaluation_error(zero_divisor).

test(div0, [error(evaluation_error(zero_divisor))]) :-
     A is 1/0.
```

#### 2.2.5 One body with multiple tests using assertions

PlUnit is designed to cooperate with the assertion/1 test provided by library(debug).^(3This integration was suggested by Günter Kniesel.) If an assertion fails in the context of a test, the test framework reports this and considers the test failed, but does not trap the debugger. Using assertion/1 in the test-body is attractive for two scenarios:

- Confirm that multiple claims hold. Where multiple claims about variable bindings can be tested using the == option in the test header, arbitrary boolean tests, notably about the state of the database, are harder to combine. Simply adding them in the body of the test has two disadvantages: it is less obvious to distinguish the tested code from the test and if one of the tests fails there is no easy way to find out which one.
- Testing‘scenarios’or sequences of actions. If one step in such a sequence fails there is again no easy way to find out which one. By inserting assertions into the sequence this becomes obvious.

Below is a simple example, showing two failing assertions. The first line of the failure message gives the test. The second reports the location of the assertion.^(4If known. The location is determined by analysing the stack. The second failure shows a case where this does not work because last-call optimization has already removed the context of the test-body.) If the assertion call originates from a different file this is reported appropriately. The last line gives the actually failed goal.

``` code
:- begin_tests(test).

test(a) :-
        A is 2^3,
        assertion(float(A)),
        assertion(A == 9).

:- end_tests(test).
```

``` code
?- run_tests.
% PL-Unit: test
ERROR: /home/jan/src/pl-devel/linux/t.pl:5:
        test a: assertion at line 7 failed
        Assertion: float(8)
ERROR: /home/jan/src/pl-devel/linux/t.pl:5:
        test a: assertion failed
        Assertion: 8==9
. done
% 2 assertions failed
```

## 3 Using separate test files

Test-units can be embedded in normal Prolog source-files. Alternatively, tests for a source-file can be placed in another file alongside the file to be tested. Test files use the extension `.plt`. The predicate [load_test_files/1](#load_test_files/1) can load all files that are related to source-files loaded into the current project.

## 4 Running the test-suite

### 4.1 Running the test suite from Prolog

To run tests from the Prolog prompt, first load the program and then run [run_tests/0](#run_tests/0) or `run_tests(+Unit)`.

**run_tests**  
Run all test-units.

**run_tests**(`+Spec`)  
Run only the specified tests. `Spec` can be a list to run multiple tests. A single specification is either the name of a test unit or a term \<`Unit`\>:\<`Tests`\>, running only the specified test. \<`Tests`\> is either the name of a test or a list of names. Running particular tests is particularly useful for tracing a test:^(5Unfortunately the body of the test is called through meta-calling, so it cannot be traced. The called user-code can be traced normally though.)

``` code
?- gtrace, run_tests(lists:member).
```

To identify nonterminating tests, interrupt the looping process with *Control-C*. The test name and location will be displayed.

### 4.2 Running the test suite from the command line

To run a file's tests from the command line, run the following command, replacing `your/file.pl` with the path to your file.

``` code
swipl -g run_tests -t halt your/file.pl
```

Prolog will (1) load the file you specify, as well as any modules it depends on; (2) run all tests in those files, and (3) exit with status 0 or 1 depending on whether the test suite succeeds or fails.

If you want to test multiple files, you can pass multiple `..pl` files.

## 5 Tests and production systems

Most applications do not want the test-suite to end up in the final application. There are several ways to achieve this. One is to place all tests in separate files and not to load the tests when creating the production environment. Alternatively, use the directive below before loading the application.

``` code
:- set_test_options([load(never)]).
```

## 6 Controlling the test suite

**set_test_options**(`+Options`)  
Defined options are:

**load**(`+Load`)  
Determines whether or not tests are loaded. When `never`, everything between [begin_tests/1](#begin_tests/1) and end_tests/1 is simply ignored. When `always`, tests are always loaded. Finally, when using the default value `normal`, tests are loaded if the code is not compiled with optimisation turned on.

**run**(`+Run`)  
Specifies when tests are run. Using `manual`, tests can only be run using [run_tests/0](#run_tests/0) or [run_tests/1](#run_tests/1). Using `make`, tests will be run for reloaded files, but not for files loaded the first time. Using `make(all)` make/0 will run all test-suites, not only those that belong to files that are reloaded.

**format**(`+Format`)  
Currently one of `tty` (default if there is a console) or `log`. `tty` uses terminal control to overwrite successful tests, allowing the user to see the currently running tests and output from failed tests. This is the default of the output is a tty. `log` prints a full log of the executed tests and their result and is intended for non-interactive usage.

**output**(`+When`)  
If `always`, emit all output as it is produced, if `never`, suppress all output and if `on_failure`, emit the output if the test fails.

**show_blocked**(`+Bool`)  
Show individual blocked tests during the report.

**occurs_check**(`+Mode`)  
Defines the default for the `occurs_check` flag during testing.

**cleanup**(`+Bool`)  
If `true` (default `false`), cleanup report at the end of [run_tests/1](#run_tests/1). Used to improve cooperation with memory debuggers such as dmalloc.

**jobs**(`Num`)  
Number of jobs to use for concurrent testing. Default is one, implying sequential testing.

**timeout**(`+Seconds`)  
Set timeout for each individual test. This acts as a default that may be overuled at the level of units or individual tests. A timeout of 0 or negative is handled as *inifinite*.

**load_test_files**(`+Options`)  
Load `.plt` test-files that belong to the currently loaded sources.

**running_tests**  
Print all currently running tests to the terminal. It can be used to find running thread in multi-threaded test operation or find the currently running test if a test appears to be blocking.

**test_report**(`+What`)  
Print report on the executed tests. `What` defines the type of report. Currently this only supports `fixme`, providing details on how the fixme-flagged tests proceeded.

## 7 Auto-generating tests

Prolog is an interactive environment. Where users of non-interactive systems tend to write tests as code, Prolog developers tend to run queries interactively during development. This interactive testing is generally faster, but the disadvantage is that the tests are lost at the end of the session. The test-wizard tries to combine the advantages. It collects toplevel queries and saves them to a specified file. Later, it extracts these queries from the file and locates the predicates that are tested by the queries. It runs the query and creates a test clause from the query.

Auto-generating test cases is experimentally supported through the library `library(test_wizard)`. We briefly introduce the functionality using examples. First step is to log the queries into a file. This is accomplished with the commands below. `Queries.pl` is the name in which to store all queries. The user can choose any filename for this purpose. Multiple Prolog instances can share the same name, as data is appended to this file and write is properly locked to avoid file corruption.

``` code
:- use_module(library(test_wizard)).
:- set_prolog_flag(log_query_file, 'Queries.pl').
```

Next, we will illustrate using the library by testing the predicates from library `library(lists)`. To generate test cases we just make calls on the terminal. Note that all queries are recorded and the system will select the appropriate ones when generating the test unit for a particular module.

``` code
?- member(b, [a,b]).
Yes
?- reverse([a,b], [b|A]).
A = [a] ;
No
```

Now we can generate the test-cases for the module list using make_tests/3:

``` code
?- make_tests(lists, 'Queries.pl', current_output).
:- begin_tests(lists).

test(member, [nondet]) :-
        member(b, [a, b]).
test(reverse, [true(A==[a])]) :-
        reverse([a, b], [b|A]).

:- end_tests(lists).
```

## 8 Portability of the test-suite

One of the reasons to have tests is to simplify migrating code between Prolog implementations. Unfortunately creating a portable test-suite implies a poor integration into the development environment. Luckily, the specification of the test-system proposed here can be ported quite easily to most Prolog systems sufficiently compatible to SWI-Prolog to consider porting your application. Most important is to have support for term_expansion/2.

In the current system, test units are compiled into sub-modules of the module in which they appear. Few Prolog systems allow for sub-modules and therefore ports may have to fall-back to inject the code in the surrounding module. This implies that support predicates used inside the test unit should not conflict with predicates of the module being tested.

### 8.1 PlUnit on SICStus

The directory of `plunit.pl` and `swi.pl` must be in the `library` search-path. With PLUNITDIR replaced accordingly, add the following into your `.sicstusrc` or `sicstus.ini`.

``` code
:- set_prolog_flag(language, iso). % for maximal compatibility
library_directory('PLUNITDIR').
```

The current version runs under SICStus 3. Open issues:

- Some messages are unformatted because SICStus 3 reports all ISO errors as instantiation errors.
- Only `plunit.pl`. Both coverage analysis and the test generation wizard currently require SWI-Prolog.
- The `load` option `normal` is the same as `always`. Use `set_test_options(load, never)` to avoid loading the test suites.
- The `run` option is not supported.
- Tests are loaded into the enclosing module instead of a separate test module. This means that predicates in the test module must not conflict with the enclosing module, nor with other test modules loaded into the same module.

## 9 Motivation of choices

### Easy to understand and flexible

There are two approaches for testing. In one extreme the tests are written using declarations dealing with setup, cleanup, running and testing the result. In the other extreme a test is simply a Prolog goal that is supposed to succeed. We have chosen to allow for any mixture of these approaches. Written down as test/1 we opt for the simple succeeding goal approach. Using options to the test the user can choose for a more declarative specification. The user can mix both approaches.

The body of the test appears at the position of a clause-body. This simplifies identification of the test body and ensures proper layout and colouring support from the editor without the need for explicit support of the unit test module. Only clauses of test/1 and test/2 may be marked as non-called in environments that perform cross-referencing.

# Index

?  
assertion/1  
[2.2.5](#idx:assertion1:14) [2.2.5](#idx:assertion1:15)

[begin_tests/1](#begin_tests/1)  
[2](#idx:begintests1:1) [6](#idx:begintests1:18)

[begin_tests/2](#begin_tests/2)  
call_cleanup/2  
[2](#idx:callcleanup2:3)

catch/3  
[2.2.4](#idx:catch3:13)

copy_term/2  
[2](#idx:copyterm2:5)

end_tests/1  
[2](#idx:endtests1:2) [6](#idx:endtests1:19)

findall/3  
[2.2.3](#idx:findall3:11)

[load_test_files/1](#load_test_files/1)  
[3](#idx:loadtestfiles1:16)

make/0  
[6](#idx:make0:22)

make_tests/3  
[7](#idx:maketests3:24)

member/2  
[2.2.1](#idx:member2:10)

[run_tests/0](#run_tests/0)  
[4.1](#idx:runtests0:17) [6](#idx:runtests0:20)

[run_tests/1](#run_tests/1)  
[6](#idx:runtests1:21) [6](#idx:runtests1:23)

[running_tests/0](#running_tests/0)  
[set_test_options/1](#set_test_options/1)  
setof/3  
[2.2.3](#idx:setof3:12)

sort/2  
[2](#idx:sort2:8)

subsumes_term/2  
[2](#idx:subsumesterm2:6)

term_expansion/2  
[8](#idx:termexpansion2:25)

test/1  
[9](#idx:test1:26) [9](#idx:test1:27)

test/2  
[2.1](#idx:test2:9) [9](#idx:test2:28)

[test_report/1](#test_report/1)  
throw/1  
[2](#idx:throw1:7)

variant/2  
[2](#idx:variant2:4)
