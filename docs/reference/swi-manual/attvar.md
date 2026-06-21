
## 8.1 Attributed variables

*Attributed variables* provide a technique for extending the Prolog unification algorithm [Holzbaur, 1992](Bibliography.html#holzbaur:1992) by hooking the binding of attributed variables. There is no consensus in the Prolog community on the exact definition and interface to attributed variables. The SWI-Prolog interface is identical to the one realised by Bart Demoen for hProlog [Demoen, 2002](Bibliography.html#Demoen:CW350). This interface is simple and available on all Prolog systems that can run the Leuven CHR system (see [chapter 9](chr.html#sec:9) and the Leuven [CHR page](https://dtai.cs.kuleuven.be/CHR/)).

Binding an attributed variable schedules a goal to be executed at the first possible opportunity. In the current implementation the hooks are executed immediately after a successful unification of the clause-head or successful completion of a foreign language (built-in) predicate. Each attribute is associated to a module, and the hook ([attr_unify_hook/2](attvar.html#attr_unify_hook/2)) is executed in this module. The example below realises a very simple and incomplete finite domain reasoner:

``` code
:- module(domain,
          [ domain/2                    % Var, ?Domain
          ]).
:- use_module(library(ordsets)).

domain(X, Dom) :-
        var(Dom), !,
        get_attr(X, domain, Dom).
domain(X, List) :-
        list_to_ord_set(List, Domain),
        put_attr(Y, domain, Domain),
        X = Y.

%       An attributed variable with attribute value Domain has been
%       assigned the value Y

attr_unify_hook(Domain, Y) :-
        (   get_attr(Y, domain, Dom2)
        ->  ord_intersection(Domain, Dom2, NewDomain),
            (   NewDomain == []
            ->  fail
            ;   NewDomain = [Value]
            ->  Y = Value
            ;   put_attr(Y, domain, NewDomain)
            )
        ;   var(Y)
        ->  put_attr( Y, domain, Domain )
        ;   ord_memberchk(Y, Domain)
        ).

%       Translate attributes from this module to residual goals

attribute_goals(X) -->
        { get_attr(X, domain, List) },
        [domain(X, List)].
```

Before explaining the code we give some example queries:

> |                                            |                     |
> |--------------------------------------------|---------------------|
> | `?- domain(X, [a,b]), X = c`               | fail                |
> | `?- domain(X, [a,b]), domain(X, [a,c]).`   | X = a               |
> | `?- domain(X, [a,b,c]), domain(X, [a,c]).` | domain(X, \[a, c\]) |

The predicate domain/2 fetches (first clause) or assigns (second clause) the variable a *domain*, a set of values the variable can be unified with. In the second clause, domain/2 first associates the domain with a fresh variable (Y) and then unifies X to this variable to deal with the possibility that X already has a domain. The predicate [attr_unify_hook/2](attvar.html#attr_unify_hook/2) (see below) is a hook called after a variable with a domain is assigned a value. In the simple case where the variable is bound to a concrete value, we simply check whether this value is in the domain. Otherwise we take the intersection of the domains and either fail if the intersection is empty (first example), assign the value if there is only one value in the intersection (second example), or assign the intersection as the new domain of the variable (third example). The nonterminal attribute_goals//1 is used to translate remaining attributes to user-readable goals that, when called, reinstate these attributes or attributes that correspond to equivalent constraints.

Implementing constraint solvers ([chapter 8](clp.html#sec:8)) is the most common, but not the only use case for attributed variables: If you implement algorithms that require efficient destructive modifications, then using attributed variables is often a more convenient and somewhat more declarative alternative for [setarg/3](manipterm.html#setarg/3) and related predicates whose sharing semantics are harder to understand. In particular, attributed variables make it easy to express graph networks and graph-oriented algorithms, since each variable can store pointers to further variables in its attributes. In such cases, the use of attributed variables should be confined within a module that exposes its functionality via more declarative interface predicates.

### 8.1.1 Attribute manipulation predicates

**attvar**(`@Term`)  
Succeeds if `Term` is an attributed variable. Note that [var/1](typetest.html#var/1) also succeeds on attributed variables. Attributed variables are created with [put_attr/3](attvar.html#put_attr/3).

**put_attr**(`+Var, +Module, +Value`)  
If `Var` is a variable or attributed variable, set the value for the attribute named `Module` to `Value`. If an attribute with this name is already associated with `Var`, the old value is replaced. Backtracking will restore the old value (i.e., an attribute is a mutable term; see also [setarg/3](manipterm.html#setarg/3)). This predicate raises an uninstantiation error if `Var` is not a variable, and a type error if `Module` is not an atom.

**get_attr**(`+Var, +Module, -Value`)  
Request the current `value` for the attribute named `Module`. If `Var` is not an attributed variable or the named attribute is not associated to `Var` this predicate fails silently. If `Module` is not an atom, a type error is raised.

**del_attr**(`+Var, +Module`)  
Delete the named attribute. If `Var` loses its last attribute it is transformed back into a traditional Prolog variable. If `Module` is not an atom, a type error is raised. In all other cases this predicate succeeds regardless of whether or not the named attribute is present.

### 8.1.2 Attributed variable hooks

Attribute names are linked to modules. This means that certain operations on attributed variables cause *hooks* to be called in the module whose name matches the attribute name.

**attr_unify_hook**(`+AttValue, +VarValue`)  
A hook that must be defined in the module to which an attributed variable refers. It is called *after* the attributed variable has been unified with a non-var term, possibly another attributed variable. `AttValue` is the attribute that was associated to the variable in this module and `VarValue` is the new value of the variable. If this predicate fails, the unification fails. If `VarValue` is another attributed variable the hook often combines the two attributes and associates the combined attribute with `VarValue` using [put_attr/3](attvar.html#put_attr/3).

To be done  
The way in which this hook currently works makes the implementation of important classes of constraint solvers impossible or at least extremely impractical. For increased generality and convenience, simultaneous unifications as in `[X,Y]=[0,1]` should be processed sequentially by the Prolog engine, or a more general hook should be provided in the future. See [Triska, 2016](Bibliography.html#clpb:Triska2016) for more information.

**attribute_goals**(`+Var`)`//`  
This nonterminal is the main mechanism in which residual constraints are obtained. It is called in every module where it is defined, and `Var` has an attribute. Its argument is that variable. In each module, attribute_goals//1 must describe a list of Prolog goals that are declaratively equivalent to the goals that caused the attributes of that module to be present and in their current state. It is always possible to do this (since these attributes stem from such goals), and it is the responsibility of constraint library authors to provide this mapping without exposing any library internals. Ideally and typically, remaining relevant attributes are mapped to *pure* and potentially simplified Prolog goals that satisfy both of the following:

- They are declaratively equivalent to the constraints that were originally posted.
- They use only predicates that are themselves exported and documented in the modules they stem from.

The latter property ensures that users can reason about residual goals, and see for themselves whether a constraint library behaves correctly. It is this property that makes it possible to thoroughly test constraint solvers by contrasting obtained residual goals with expected answers.

This nonterminal is used by [copy_term/3](attvar.html#copy_term/3), on which the Prolog top level relies to ensure the basic invariant of pure Prolog programs: The answer is *declaratively equivalent* to the query.

The [copy_term/3](attvar.html#copy_term/3) primitive uses attribute_goals//1 inside a [findall/3](allsolutions.html#findall/3) call. This implies that attribute_goals//1 can unify variables and modify attributes, for example, to tell other hooks that some attribute has already been taken care of. This nonterminal is also used by [frozen/2](coroutining.html#frozen/2) which does *not* create a copy. Ideally attribute_goals//1 should not modify anything to allow direct application in [frozen/2](coroutining.html#frozen/2). In the current implementation [frozen/2](coroutining.html#frozen/2) backtracks over attribute_goals//1 to tolerate the current behavior. This work-around harms the performance of [frozen/2](coroutining.html#frozen/2). New implementations of attribute_goals//1 should avoid relying on backtracking when feasible. Future versions of [frozen/2](coroutining.html#frozen/2) and [copy_term/3](attvar.html#copy_term/3) may require attribute_goals//1 not to modify any variables or attributes.

Note that instead of *defaulty* representations, a Prolog *list* is used to represent residual goals. This simplifies processing and reasoning about residual goals throughout all programs that need this functionality.

**project_attributes**(`+QueryVars, +ResidualVars`)  
A hook that can be defined in each module to project constraints on newly introduced variables back to the query variables. `QueryVars` is the list of variables occurring in the query and `ResidualVars` is a list of variables that have attributes attached. There may be variables that occur in both lists. If possible, [project_attributes/2](attvar.html#project_attributes/2) should change the attributes so that all constraints are expressed as residual goals that refer only to `QueryVars`, while other variables are existentially quantified.

\[deprecated\]**attr_portray_hook**(`+AttValue, +Var`)  
Called by [write_term/2](termrw.html#write_term/2) and friends for each attribute if the option `attributes(portray)` is in effect. If the hook succeeds the attribute is considered printed. Otherwise `Module = ...` is printed to indicate the existence of a variable. This predicate is deprecated because it cannot work with pure interface predicates like [copy_term/3](attvar.html#copy_term/3). Use attribute_goals//1 instead to map attributes to residual goals.

### 8.1.3 Operations on terms with attributed variables

**copy_term**(`+Term, -Copy, -Gs`)  
Create a regular term `Copy` as a copy of `Term` (without any attributes), and a list `Gs` of goals that represents the attributes. The goal `maplist(call, Gs)` recreates the attributes for `Copy`. The nonterminal attribute_goals//1, as defined in the modules the attributes stem from, is used to convert attributes to lists of goals.

This building block is used by the top level to report pending attributes in a portable and understandable fashion. This predicate is the preferred way to reason about and communicate terms with constraints.

The form `copy_term(Term, Term, Gs)` can be used to reason about the constraints in which `Term` is involved.

**copy_term_nat**(`+Term, -Copy`)  
As [copy_term/2](manipterm.html#copy_term/2). Attributes, however, are *not* copied but replaced by fresh variables.

**term_attvars**(`+Term, -AttVars`)  
`AttVars` is a list of all attributed variables in `Term` and its attributes. That is, [term_attvars/2](attvar.html#term_attvars/2) works recursively through attributes. This predicate is cycle-safe. The goal `term_attvars(Term,[])` is an efficient test that `Term` has *no* attributes; scanning the term is aborted after the first attributed variable is found.

### 8.1.4 Special purpose predicates for attributes

Normal user code should deal with [put_attr/3](attvar.html#put_attr/3), [get_attr/3](attvar.html#get_attr/3) and [del_attr/2](attvar.html#del_attr/2). The routines in this section fetch or set the entire attribute list of a variable. Use of these predicates is anticipated to be restricted to printing and other special purpose operations.

**get_attrs**(`+Var, -Attributes`)  
Get all attributes of `Var`. `Attributes` is a term of the form `att(Module, Value, MoreAttributes)`, where `MoreAttributes` is `[]` for the last attribute.

**put_attrs**(`+Var, -Attributes`)  
Set all attributes of `Var`. See [get_attrs/2](attvar.html#get_attrs/2) for a description of `Attributes`.

**del_attrs**(`+Var`)  
If `Var` is an attributed variable, delete *all* its attributes. In all other cases, this predicate succeeds without side-effects.
