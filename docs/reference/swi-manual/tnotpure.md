
## 7.4 Tabling for impure programs

Tabling guarantees logically correct results and termination provided the computation only involves terms of bounded size on *pure* Prolog programs, i.e., Prolog programs without side effects or pruning of choice points (cut, [-\>/2](control.html#-%3E/2), etc.). Notably pruning choice points of an incomplete tabled goal may cause an incomplete table and thus cause subsequent queries for the same goal to return an incomplete set of answers. The current SWI-Prolog implementation provides several mechanisms to improve on this situation.

- *Dynamic Strongly Connected Components (SCC)*  
  Tabled goals are *completed* as soon as possible. Each fresh tabled goal creates a scheduling component which the system attempts to solve immediately. If a subgoal of the fresh goal refers to an incomplete tabled goal the scheduling components for both goals are merged such that the related goals are completed together. Dynamic rather than static determination of strongly connected components guarantees that the components are minimal because only actually reached code needs to be considered rather than maximally reachable code.

  Minimal SCCs imply that goals are completed as early as possible. This implies that tabled goals may be embedded in e.g., [findall/3](allsolutions.html#findall/3) or be used as a condition as long as there is no dependency (*loop*) with goals outside the [findall/3](allsolutions.html#findall/3) or condition. For example, the code below misbehaves when called as `p(X)` because the argument of [findall/3](allsolutions.html#findall/3) calls a *variant* of the goal that initiated the findall goal. A call `p(1)` however is ok as `p(1)` is not a variant of `p(X)`.

  ``` code
  p(X) :-
      findall(Y, p(Y), Ys),
      ...
  ```

- *Early completion*  
  Ground goals, i.e., goals without variables, are subject to early completion. This implies they are considered completed after the first solution.
