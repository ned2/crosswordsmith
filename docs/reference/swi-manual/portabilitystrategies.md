
## C.1 Some considerations for writing portable code

The traditional way to write portable code is to define custom predicates for all potentially non-portable code and define these separately for all Prolog dialects one wishes to support. Here are some considerations.

- Probably the best reason for this is that it allows to define minimal semantics required by the application for the portability predicates. Such functionality can often be mapped efficiently to the target dialect. Contrary, if code was written for dialect `X`, the defined semantics are those of dialect `X`. Emulating all extreme cases and full error handling compatibility may be tedious and result in a much slower implementation than needed. Take for example [call_cleanup/2](metacall.html#call_cleanup/2). The SICStus definition is fundamentally different from the SWI definition, but 99% of the applications just want to make calls like below to guarantee `StreamIn` is closed, even if process/1 misbehaves.

  ``` code
          call_cleanup(process(StreamIn), close(In))
  ```

- As a drawback, the code becomes full of *my_call_cleanup*, etc. and every potential portability conflict needs to be abstracted. It is hard for people who have to maintain such code later to grasp the exact semantics of the *my\_\** predicates and applications that combine multiple libraries using this compatibility approach are likely to encounter conflicts between the portability layers. A good start is not to use *my\_\**, but a prefix derived from the library or application name or names that explain the intended semantics more precisely.

- Another problem is that most code is initially not written with portability in mind. Instead, ports are requested by users or arise from the desire to switch Prolog dialect. Typically, we want to achieve compatibility with the new Prolog dialect with minimal changes, often keeping compatibility with the original dialect(s). This problem is well known from the C/Unix world and we advise anyone to study the philosophy of [GNU autoconf](http://www.gnu.org/software/autoconf/), from which we will illustrate some highlights below.

The GNU autoconf suite, known to most people as **configure**, was an answer to the frustrating life of Unix/C programmers when Unix dialects were about as abundant and poorly standardised as Prolog dialects today. Writing a portable C program can only be achieved using cpp, the C preprocessor. The C preprocessor performs two tasks: macro expansion and conditional compilation. Prolog realises macro expansion through [term_expansion/2](consulting.html#term_expansion/2) and [goal_expansion/2](consulting.html#goal_expansion/2). Conditional compilation is achieved using `:- if(Condition)` as explained in [section 4.3.1.2](consulting.html#sec:4.3.1.2). The situation appears similar.

The important lesson learned from GNU autoconf is that the *last* resort for conditional compilation to achieve portability is to switch on the platform or dialect. Instead, GNU autoconf allows you to write tests for specific properties of the platform. Most of these are whether or not some function or file is available. Then there are some standard tests for difficult-to-write-portable situations and finally there is a framework that allows you to write arbitrary C programs and check whether they can be compiled and/or whether they show the intended behaviour. Using a separate **configure** program is needed in C, as you cannot perform C compilation step or run C programs from the C preprocessor. In most Prolog environments we do not need this distinction as the compiler is integrated into the runtime environment and Prolog has excellent reflexion capabilities.

We must learn from the distinction to test for features instead of platform (dialect), as this makes the platform-specific code robust for future changes of the dialect. Suppose we need [compare/3](compare.html#compare/3) as defined in this manual. The [compare/3](compare.html#compare/3) predicate is not part of the ISO standard, but many systems support it and it is not unlikely it will become ISO standard or the intended dialect will start supporting it. GNU autoconf strongly advises to test for the availability:

``` code
:- if(\+current_predicate(_, compare(_,_,_))).
compare(<, Term1, Term2) :-
        Term1 @< Term2, !.
compare(>, Term1, Term2) :-
        Term1 @> Term2, !.
compare(=, Term1, Term2) :-
        Term1 == Term2.
:- endif.
```

This code is **much** more robust against changes to the intended dialect and, possibly at least as important, will provide compatibility with dialects you didn't even consider porting to right now.

In a more challenging case, the target Prolog has [compare/3](compare.html#compare/3), but the semantics are different. What to do? One option is to write a my_compare/3 and change all occurrences in the code. Alternatively you can rename calls using [goal_expansion/2](consulting.html#goal_expansion/2) like below. This construct will not only deal with Prolog dialects lacking [compare/3](compare.html#compare/3) as well as those that only implement it for numeric comparison or have changed the argument order. Of course, writing rock-solid code would require a complete test-suite, but this example will probably cover all Prolog dialects that allow for conditional compilation, have core ISO facilities and provide [goal_expansion/2](consulting.html#goal_expansion/2), the things we claim a Prolog dialect should have to start writing portable code for it.

``` code
:- if(\+catch(compare(<,a,b), _, fail)).
compare_standard_order(<, Term1, Term2) :-
        Term1 @< Term2, !.
compare_standard_order(>, Term1, Term2) :-
        Term1 @> Term2, !.
compare_standard_order(=, Term1, Term2) :-
        Term1 == Term2.

goal_expansion(compare(Order, Term1, Term2),
               compare_standard_order(Order, Term1, Term2)).
:- endif.
```
