
## 14.5 Debugging and updating deployed systems

SWI-Prolog provides several facilities to debug and update running (server) applications. The core to these facilities are:

- Hot-swap recompilation ([section 4.3.2](consulting.html#sec:4.3.2) and the library `library(hotswap)`) allow, with some limitation, making modifications to running services. This includes adding debugging and logging statements.
- To make this useful some form of interaction is required. This can be implemented using signal handlers (Unix), specific HTTP services, generic HTTP services (e.g., [SWISH](https://swish.swi-prolog.org)) or networked interaction using the library `library(prolog_server)` that allow interaction using netcat (**nc**) or **telnet**.
