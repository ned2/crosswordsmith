
## A.36 library(persistency): Provide persistent dynamic predicates

To be done  
\- Provide type safety while loading  
- Thread safety must now be provided at the user-level. Can we provide generic thread safety? Basically, this means that we must wrap all exported predicates. That might better be done outside this library.  
- Transaction management?  
- Should assert\_`<`name`>` only assert if the database does not contain a variant?  
- Since we have [prolog_listen/2](prolog-event.html#prolog_listen/2), we could use direct [assert/1](db.html#assert/1) and [retract/1](db.html#retract/1) and use the system hooks to deal with the updates.

This module provides simple persistent storage for one or more dynamic predicates. A database is always associated with a module. A module that wishes to maintain a database must declare the terms that can be placed in the database using the directive [persistent/1](persistency.html#persistent/1).

The [persistent/1](persistency.html#persistent/1) expands each declaration into five predicates:

- `name(Arg, ...)`
- `assert_name(Arg, ...)`
- `asserta_name(Arg, ...)`
- `retract_name(Arg, ...)`
- `retractall_name(Arg, ...)`

As mentioned, a database can only be accessed from within a single module. This limitation is on purpose, forcing the user to provide a proper API for accessing the shared persistent data.

This module requires the same thread-synchronization as the normal Prolog database. This implies that if each individual assert or retract takes the database from one consistent state to the next, no additional locking is required. If more than one elementary database operation is required to get from one consistent state to the next, both updating and querying the database must be locked using [with_mutex/2](threadsync.html#with_mutex/2).

Below is a simple example, where adding a user does not need locking as it is a single *assert*, while modifying a user requires both a retract and assert and thus needs to be locked.

``` code
:- module(user_db,
          [ attach_user_db/1,           % +File
            current_user_role/2,        % ?User, ?Role
            add_user/2,                 % +User, +Role
            set_user_role/2             % +User, +Role
          ]).
:- use_module(library(persistency)).

:- persistent
        user_role(name:atom, role:oneof([user,administrator])).

attach_user_db(File) :-
        db_attach(File, []).

%%      current_user_role(+Name, -Role) is semidet.

current_user_role(Name, Role) :-
        with_mutex(user_db, user_role(Name, Role)).

add_user(Name, Role) :-
        assert_user_role(Name, Role).

set_user_role(Name, Role) :-
        user_role(Name, Role), !.
set_user_role(Name, Role) :-
        with_mutex(user_db,
                   (  retractall_user_role(Name, _),
                      assert_user_role(Name, Role))).
```

**persistent** `+Spec`  
Declare dynamic database terms. Declarations appear in a directive and have the following format:

``` code
:- persistent
        <callable>,
        <callable>,
        ...
```

Each specification is a callable term, following the conventions of `library(record)`, where each argument is of the form

``` code
name:type
```

Types are defined by `library(error)`.

\[nondet\]**current_persistent_predicate**(`:PI`)  
True if `PI` is a predicate that provides access to the persistent database DB.

**db_attach**(`:File, +Options`)  
Use `File` as persistent database for the calling module. The calling module must defined [persistent/1](persistency.html#persistent/1) to declare the database terms. Defined options:

**sync**(`+Sync`)  
One of `close` (close journal after write), `flush` (default, flush journal after write) or `none` (handle as fully buffered stream).

If `File` is already attached this operation may change the `sync` behaviour.

\[semidet\]**db_attached**(`:File`)  
True if the context module attached to the persistent database `File`.

\[det\]**db_assert**(`:Term`)  
Assert `Term` into the database and record it for persistency. Note that if the on-disk file has been modified it is first reloaded.

\[det\]**db_detach**  
Detach persistency from the calling module and delete all persistent clauses from the Prolog database. Note that the file is not affected. After this operation another file may be attached, providing it satisfies the same persistency declaration.

\[det\]**db_retractall**(`:Term`)  
Retract all matching facts and do the same in the database. If `Term` is unbound, [persistent/1](persistency.html#persistent/1) from the calling module is used as generator.

\[nondet\]**db_retract**(`:Term`)  
Retract terms from the database one-by-one.

**db_sync**(`:What`)  
Synchronise database with the associated file. `What` is one of:

**reload**  
Database is reloaded from file if the file was modified since loaded.

**update**  
As `reload`, but use incremental loading if possible. This allows for two processes to examine the same database file, where one writes the database and the other periodycally calls `db_sync(update)` to follow the modified data.

**gc**  
Database was re-written, deleting all retractall statements. This is the same as `gc(50)`.

**gc**(`Percentage`)  
GC DB if the number of deleted terms is greater than the given percentage of the total number of terms.

**gc**(`always`)  
GC DB without checking the percentage.

**close**  
Database stream was closed

**detach**  
Remove all registered persistency for the calling module

**nop**  
No-operation performed

With unbound `What`, [db_sync/1](persistency.html#db_sync/1) reloads the database if it was modified on disk, gc it if it is dirty and close it if it is opened.

**db_sync_all**(`+What`)  
Sync all registered databases.
