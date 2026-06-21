
## E.2 Software support to keep track of license conditions

Given the above, it is possible that SWI-Prolog packages and extensions rely on the GPL, LGPL or other licenses. The predicates below allow for registering license requirements for Prolog files and foreign modules. The predicate [license/0](softlicense.html#license/0) reports which components from the currently configured system are distributed under non-permissive open source licenses and therefore may need to be replaced to suit your requirements.

**license**  
Evaluate the license conditions of all loaded components. If the system contains one or more components that are licensed under GPL-like restrictions the system indicates this program may only be distributed under the `GPL` license as well as which components prohibit the use of other license conditions. Likewise for for LGPL components.

**license**(`+LicenseId, +Component`)  
Register the fact that `Component` is distributed under a license identified by `LicenseId`. Known license identifiers can be listed using [known_licenses/0](softlicense.html#known_licenses/0). A new license can be registered as a known language using a declaration like below. The second argument defines the *category* if the license, which is one of `gpl`, `lgpl`, `permissive` or `proprietary`.

``` code
:- multifile license:license/3.

license:license(mylicense, permissive,
                [ comment('My personal license'),
                  url('http://www.mine.org/license.html')
                ]).

:- license(mylicense).
```

**license**(`+LicenseId`)  
Intended as a directive in Prolog source files. It takes the current filename and calls [license/2](softlicense.html#license/2).

`void` **PL_license**(`const char *LicenseId, const char *Component`)  
Intended for the **install()** procedure of foreign libraries. This call can be made *before* [PL_initialise()](foreigninclude.html#PL_initialise()).

**known_licenses**  
List all licenses *known* to the system. This does not imply the system contains code covered by the listed licenses. See [license/2](softlicense.html#license/2).
