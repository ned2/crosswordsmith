
## 9.4 Debugging CHR programs

The CHR debugging facilities are currently rather limited. Only tracing is currently available. To use the CHR debugging facilities for a CHR file it must be compiled for debugging. Generating debug info is controlled by the CHR option [debug](flags.html#flag:debug), whose default is derived from the SWI-Prolog flag [generate_debug_info](flags.html#flag:generate_debug_info). Therefore debug info is provided unless the **--no-debug** is used.

### 9.4.1 CHR debug ports

For CHR constraints the four standard ports are defined:

**call**  
A new constraint is called and becomes active.

**exit**  
An active constraint exits: it has either been inserted in the store after trying all rules or has been removed from the constraint store.

**fail**  
An active constraint fails.

**redo**  
An active constraint starts looking for an alternative solution.

In addition to the above ports, CHR constraints have five additional ports:

**wake**  
A suspended constraint is woken and becomes active.

**insert**  
An active constraint has tried all rules and is suspended in the constraint store.

**remove**  
An active or passive constraint is removed from the constraint store.

**try**  
An active constraint tries a rule with possibly some passive constraints. The try port is entered just before committing to the rule.

**apply**  
An active constraint commits to a rule with possibly some passive constraints. The apply port is entered just after committing to the rule.

### 9.4.2 Tracing CHR programs

Tracing is enabled with the [chr_trace/0](chr-debugging.html#chr_trace/0) predicate and disabled with the [chr_notrace/0](chr-debugging.html#chr_notrace/0) predicate.

When enabled the tracer will step through the `call`, `exit`, `fail`, `wake` and `apply` ports, accepting debug commands, and simply write out the other ports.

The following debug commands are currently supported:

``` code
        CHR debug options:

                <cr>    creep           c       creep
                s       skip
                g       ancestors
                n       nodebug
                b       break
                a       abort
                f       fail
                ?       help            h       help
```

Their meaning is:

**creep**  
Step to the next port.

**skip**  
Skip to exit port of this call or wake port.

**ancestors**  
Print list of ancestor call and wake ports.

**nodebug**  
Disable the tracer.

**break**  
Enter a recursive Prolog top level. See [break/0](toplevel.html#break/0).

**abort**  
Exit to the top level. See [abort/0](toplevel.html#abort/0).

**fail**  
Insert failure in execution.

**help**  
Print the above available debug options.

### 9.4.3 CHR Debugging Predicates

The `library(chr)` module contains several predicates that allow inspecting and printing the content of the constraint store.

**chr_trace**  
Activate the CHR tracer. By default the CHR tracer is activated and deactivated automatically by the Prolog predicates [trace/0](debugger.html#trace/0) and [notrace/0](debugger.html#notrace/0).

**chr_notrace**  
Deactivate the CHR tracer. By default the CHR tracer is activated and deactivated automatically by the Prolog predicates [trace/0](debugger.html#trace/0) and [notrace/0](debugger.html#notrace/0).

**chr_leash**(`+Spec`)  
Define the set of CHR ports on which the CHR tracer asks for user intervention (i.e. stops). `Spec` is either a list of ports as defined in [section 9.4.1](chr-debugging.html#sec:9.4.1) or a predefined‘alias’. Defined aliases are: `full` to stop at all ports, `none` or `off` to never stop, and `default` to stop at the `call`, `exit`, `fail`, `wake` and `apply` ports. See also [leash/1](debugger.html#leash/1).

**chr_show_store**(`+Mod`)  
Prints all suspended constraints of module `Mod` to the standard output. This predicate is automatically called by the SWI-Prolog top level at the end of each query for every CHR module currently loaded. The Prolog flag `chr_toplevel_show_store` controls whether the top level shows the constraint stores. The value `true` enables it. Any other value disables it.

**find_chr_constraint**(`-Constraint`)  
Returns a constraint in the constraint store. Via backtracking, all constraints in the store can be enumerated.
