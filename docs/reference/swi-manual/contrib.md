
## E.1 Contributing to the SWI-Prolog project

To reach maximal coherence we will, as a rule of thumb, only accept new code that has the Simplified BSD license and existing code with a *permissive* license such as MIT, Apache, BSD-3, etc. In exceptional cases we may accept code with GPL or LGPL conditions. Such code must be tagged using a [license/1](softlicense.html#license/1) directive (Prolog) or a call to [PL_license()](softlicense.html#PL_license()) for foreign code and, if they are part of the core, the code must be excluded using the `--without-gpl` or `--without-lgpl` option.
