
## A.56 library(simplex): Solve linear programming problems

author  
[Markus Triska](https://www.metalevel.at)

### A.56.1 Introduction

A **linear programming problem** or simply **linear program** (LP) consists of:

- a set of *linear* **constraints**
- a set of **variables**
- a *linear* **objective function**.

The goal is to assign values to the variables so as to *maximize* (or minimize) the value of the objective function while satisfying all constraints.

Many optimization problems can be modeled in this way. As one basic example, consider a knapsack with fixed capacity C, and a number of items with sizes `s(i)` and values `v(i)`. The goal is to put as many items as possible in the knapsack (not exceeding its capacity) while maximizing the sum of their values.

As another example, suppose you are given a set of *coins* with certain values, and you are to find the minimum number of coins such that their values sum up to a fixed amount. Instances of these problems are solved below.

Solving an LP or integer linear program (ILP) with this library typically comprises 4 stages:

1.  an initial state is generated with [gen_state/1](simplex.html#gen_state/1)
2.  all relevant constraints are added with [constraint/3](simplex.html#constraint/3)
3.  [maximize/3](simplex.html#maximize/3) or [minimize/3](simplex.html#minimize/3) are used to obtain a *solved state* that represents an optimum solution
4.  [variable_value/3](simplex.html#variable_value/3) and [objective/2](simplex.html#objective/2) are used on the solved state to obtain variable values and the objective function at the optimum.

The most frequently used predicates are thus:

**gen_state**(`-State`)  
Generates an initial state corresponding to an empty linear program.

**constraint**(`+Constraint, +S0, -S`)  
Adds a linear or integrality constraint to the linear program corresponding to state `S0`. A linear constraint is of the form `Left Op C`, where `Left` is a list of `Coefficient*Variable` terms (variables in the context of linear programs can be atoms or compound terms) and `C` is a non-negative numeric constant. The list represents the sum of its elements. `Op` can be `=`, `=<` or `>=`. The coefficient `1` can be omitted. An integrality constraint is of the form `integral(Variable)` and constrains Variable to an integral value.

**maximize**(`+Objective, +S0, -S`)  
Maximizes the objective function, stated as a list of `Coefficient*Variable` terms that represents the sum of its elements, with respect to the linear program corresponding to state `S0`. `\`arg{`S`} is unified with an internal representation of the solved instance.

**minimize**(`+Objective, +S0, -S`)  
Analogous to [maximize/3](simplex.html#maximize/3).

**variable_value**(`+State, +Variable, -Value`)  
`Value` is unified with the value obtained for `Variable`. `State` must correspond to a solved instance.

**objective**(`+State, -Objective`)  
Unifies `Objective` with the result of the objective function at the obtained extremum. `State` must correspond to a solved instance.

All numeric quantities are converted to rationals via rationalize/1, and rational arithmetic is used throughout solving linear programs. In the current implementation, all variables are implicitly constrained to be *non-negative*. This may change in future versions, and non-negativity constraints should therefore be stated explicitly.

### A.56.2 Delayed column generation

*Delayed column generation* means that more constraint columns are added to an existing LP. The following predicates are frequently used when this method is applied:

**constraint**(`+Name, +Constraint, +S0, -S`)  
Like [constraint/3](simplex.html#constraint/3), and attaches the name `Name` (an atom or compound term) to the new constraint.

**shadow_price**(`+State, +Name, -Value`)  
Unifies `Value` with the shadow price corresponding to the linear constraint whose name is `Name`. `State` must correspond to a solved instance.

**constraint_add**(`+Name, +Left, +S0, -S`)  
`Left` is a list of `Coefficient*Variable` terms. The terms are added to the left-hand side of the constraint named `Name`. `S` is unified with the resulting state.

An example application of *delayed column generation* to solve a *bin packing* task is available from: [**metalevel.at/various/colgen/**](https://www.metalevel.at/various/colgen/)

### A.56.3 Solving LPs with special structure

The following predicates allow you to solve specific kinds of LPs more efficiently:

**transportation**(`+Supplies, +Demands, +Costs, -Transport`)  
Solves a transportation problem. `Supplies` and `Demands` must be lists of non-negative integers. Their respective sums must be equal. `Costs` is a list of lists representing the cost matrix, where an entry (*i*,*j*) denotes the integer cost of transporting one unit from *i* to *j*. A transportation plan having minimum cost is computed and unified with `Transport` in the form of a list of lists that represents the transportation matrix, where element (*i*,*j*) denotes how many units to ship from *i* to *j*.

**assignment**(`+Cost, -Assignment`)  
Solves a linear assignment problem. `Cost` is a list of lists representing the quadratic cost matrix, where element (i,j) denotes the integer cost of assigning entity \$i\$ to entity \$j\$. An assignment with minimal cost is computed and unified with `Assignment` as a list of lists, representing an adjacency matrix.

### A.56.4 Examples

We include a few examples for solving LPs with this library.

#### A.56.4.1 Example 1

This is the "radiation therapy" example, taken from *Introduction to Operations Research* by Hillier and Lieberman.

[**Prolog DCG notation**](https://www.metalevel.at/prolog/dcg) is used to *implicitly* thread the state through posting the constraints:

``` code
:- use_module(library(simplex)).

radiation(S) :-
        gen_state(S0),
        post_constraints(S0, S1),
        minimize([0.4*x1, 0.5*x2], S1, S).

post_constraints -->
        constraint([0.3*x1, 0.1*x2] =< 2.7),
        constraint([0.5*x1, 0.5*x2] = 6),
        constraint([0.6*x1, 0.4*x2] >= 6),
        constraint([x1] >= 0),
        constraint([x2] >= 0).
```

An example query:

``` code
?- radiation(S), variable_value(S, x1, Val1),
                 variable_value(S, x2, Val2).
Val1 = 15 rdiv 2,
Val2 = 9 rdiv 2.
```

#### A.56.4.2 Example 2

Here is an instance of the knapsack problem described above, where `C = 8`, and we have two types of items: One item with value 7 and size 6, and 2 items each having size 4 and value 4. We introduce two variables, `x(1)` and `x(2)` that denote how many items to take of each type.

``` code
:- use_module(library(simplex)).

knapsack(S) :-
        knapsack_constraints(S0),
        maximize([7*x(1), 4*x(2)], S0, S).

knapsack_constraints(S) :-
        gen_state(S0),
        constraint([6*x(1), 4*x(2)] =< 8, S0, S1),
        constraint([x(1)] =< 1, S1, S2),
        constraint([x(2)] =< 2, S2, S).
```

An example query yields:

``` code
?- knapsack(S), variable_value(S, x(1), X1),
                variable_value(S, x(2), X2).
X1 = 1
X2 = 1 rdiv 2.
```

That is, we are to take the one item of the first type, and half of one of the items of the other type to maximize the total value of items in the knapsack.

If items can not be split, integrality constraints have to be imposed:

``` code
knapsack_integral(S) :-
        knapsack_constraints(S0),
        constraint(integral(x(1)), S0, S1),
        constraint(integral(x(2)), S1, S2),
        maximize([7*x(1), 4*x(2)], S2, S).
```

Now the result is different:

``` code
?- knapsack_integral(S), variable_value(S, x(1), X1),
                         variable_value(S, x(2), X2).

X1 = 0
X2 = 2
```

That is, we are to take only the *two* items of the second type. Notice in particular that always choosing the remaining item with best performance (ratio of value to size) that still fits in the knapsack does not necessarily yield an optimal solution in the presence of integrality constraints.

#### A.56.4.3 Example 3

We are given:

- 3 coins each worth 1 unit
- 20 coins each worth 5 units and
- 10 coins each worth 20 units.

The task is to find a *minimal* number of these coins that amount to 111 units in total. We introduce variables `c(1)`, `c(5)` and `c(20)` denoting how many coins to take of the respective type:

``` code
:- use_module(library(simplex)).

coins(S) :-
        gen_state(S0),
        coins(S0, S).

coins -->
        constraint([c(1), 5*c(5), 20*c(20)] = 111),
        constraint([c(1)] =< 3),
        constraint([c(5)] =< 20),
        constraint([c(20)] =< 10),
        constraint([c(1)] >= 0),
        constraint([c(5)] >= 0),
        constraint([c(20)] >= 0),
        constraint(integral(c(1))),
        constraint(integral(c(5))),
        constraint(integral(c(20))),
        minimize([c(1), c(5), c(20)]).
```

An example query:

``` code
?- coins(S), variable_value(S, c(1), C1),
             variable_value(S, c(5), C5),
             variable_value(S, c(20), C20).

C1 = 1,
C5 = 2,
C20 = 5.
```
