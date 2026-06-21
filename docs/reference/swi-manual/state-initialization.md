
## 14.3 State initialization

The [initialization/1](consulting.html#initialization/1) and [initialization/2](consulting.html#initialization/2) directive may be used to register goals to be executed at various points in the life cycle of an executable. Alternatively, one may consider *lazy initialization* which typically follows the pattern below. Single threaded code can avoid using [with_mutex/2](threadsync.html#with_mutex/2).

``` code
:- dynamic x_done/0.
:- volatile x_done/0.

x(X) :-
    x_done,
    !,
    use_x(X).
x(X) :-
    with_mutex(x, create_x),
    use_x(X).

create_x :-
    x_done,
    !.
create_x :-
    <create x>
    asserta(x_done).
```
