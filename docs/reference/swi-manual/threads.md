
# 10 Multithreaded applications

SWI-Prolog multithreading is based on standard C language multithreading support. It is not like *ParLog* or other parallel implementations of the Prolog language. Prolog threads have their own stacks and only share the Prolog *heap*: predicates, records, flags and other global non-backtrackable data. SWI-Prolog thread support is designed with the following goals in mind.

- *Multithreaded server applications*  
  Today's computing services often focus on (internet) server applications. Such applications often have need for communication between services and/or fast non-blocking service to multiple concurrent clients. The shared heap provides fast communication, and thread creation is relatively cheap.^(196On an Intel i7-2600K, running Ubuntu Linux 12.04, SWI-Prolog 6.2 creates and joins 32,000 threads per second elapsed time.)
- *Interactive applications*  
  Interactive applications often need to perform extensive computation. If such computations are executed in a new thread, the main thread can process events and allow the user to cancel the ongoing computation. User interfaces can also use multiple threads, each thread dealing with input from a distinct group of windows. See also [section 10.7](mt-xpce.html#sec:10.7).
- *Natural integration with foreign code*  
  Each Prolog thread runs in a native thread of the operating system, automatically making them cooperate with *MT-safe* foreign code. In addition, any foreign thread can create its own Prolog engine for dealing with calling Prolog from C code.

SWI-Prolog multithreading is based on the POSIX thread standard [Butenhof, 1997](Bibliography.html#Butenhof:1997:PPT) used on most popular systems except for MS-Windows. On Windows it uses the [pthread-win32](http://sources.redhat.com/pthreads-win32/) emulation of POSIX threads mixed with the Windows native API for smoother and faster operation. The SWI-Prolog thread implementation has been discussed in the ISO WG17 working group and is largely adopted by YAP and XSB Prolog.^(197The latest version of the ISO draft can be found at [http://logtalk.org/plstd/threads.pdf](http://logtalk.org/plstd/threads.pdf). It appears to have dropped from the ISO WG17 agenda.)

------------------------------------------------------------------------

## Section Index

------------------------------------------------------------------------

[10.1 Creating and destroying Prolog threads](threadcreate.html)

[10.2 Monitoring threads](thmonitor.html)

[10.3 Thread communication](threadcom.html)

[10.3.1 Message queues](threadcom.html#sec:10.3.1)

[10.3.2 Waiting for events](threadcom.html#sec:10.3.2)

[10.3.3 Signalling threads](threadcom.html#sec:10.3.3)

[10.3.4 Threads and dynamic predicates](threadcom.html#sec:10.3.4)

[10.4 Thread synchronisation](threadsync.html)

[10.5 library(threadutil): Interactive thread utilities](threadutil.html)

[10.6 Multithreaded mixed C and Prolog applications](foreignthread.html)

[10.6.1 A Prolog thread for each native thread (one-to-one)](foreignthread.html#sec:10.6.1)

[10.6.2 Using Prolog engines from C](foreignthread.html#sec:10.6.2)

[10.7 Multithreading and the XPCE graphics system](mt-xpce.html)
