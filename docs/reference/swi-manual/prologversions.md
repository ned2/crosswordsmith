
## A.44 library(prolog_versions): Demand specific (Prolog) versions

To be done  
\- Not only provide a minimal version but a more version ranges, exclude certain versions, etc.  
- More features and better messages to help the user resolving problems.

This library is provided to make it easier to reason about software versions, in particular require that that hosting Prolog system is of the right version and provides the features required by the library or application.

\[det\]**require_prolog_version**(`+Required, +Features:list`)  
Claim that the running Prolog version is at least version `Required` and provides the requested `Features`. `Required` is an expression of versions. At the lowest level, a version is an atom or string that provides the version as

``` code
Major.Minor[[.Patch][[-GitRev],-GitHash]]]
```

Example strings are’8.5’,’8.5.0’,’8.5.0-50’,’8.5.0-69-gad38e8ad8\`. The last two require fetching the sources from git or using the Windows daily builds.

Versions may be embedded in a comparison operator (`<`, `=<`, `=`, `>=` or `>`), e.g., `=<('9.1')`. Versions are considered to compare equal only on the components of the `Required` version. I.e., `'9.1'` compares equal to `'9.1.2'`.

Version expressions can be constructed from the Prolog operators’,’/2,’;’/2 and’`\+`’/1. An example of a complicated expression is below, which demands major version 9, but considers 9.1.2 not suitable.

``` code
(>=('9'), \+(=('9.1.2')))
```

`Features` is a list of required or preferred features. Individual features are:

**warning**(`Feature`)  
Only print a warning instead of throwing an error.

**library**(`Lib`)  
Demand `library(Lib)` to be present. Thde library not being there may indicate an incomplete installation. For example `library(pce)` to demand xpce graphics support.

**`Flag`**  
Demand `current_prolog_flag(Flag, true)` to be true.

**`FlagValue`**  
If `FlagValue` is Flag(Value), demand `current_prolog_flag(Flag, Value)` to be true.

Errors  
\- `version_error('SWI-Prolog', PrologVersion, Cmp, Required)`  
- `existence_error(prolog_feature, Feature)`

\[det\]**require_version**(`+Component, +Available, +CmpRequired`)  
Require `Component` to have version `CmpRequired`, while `Component` is know to have version `Available`.

Errors  
`version_error(Component, Required, Cmp, Available)`

\[semidet\]**cmp_versions**(`?Cmp, +V1, +V2`)  
Compare to versions. `Cmp` is one of `<`, `=<`, `=`, `>=` or `>`. If `Cmp` is unbound we check whether `<` or `>` hold or else bind `Cmp` to `=`.

When comparing for equality (`=`), the versions are considered equal if they compare equal up to the detail level of the least specified. E.g,’9.1.2’is considered equal to’9.1’.
