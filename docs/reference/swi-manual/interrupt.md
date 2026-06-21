
## B.4 Simulating a debugger interrupt

**prolog_interrupt**  
Calls the functionality that allows for debugging after receiving (normally) `SIGINT`. This may be used in IDE environments to start debugging a toplevel thread by injecting this into the target thread using [thread_signal/2](threadcom.html#thread_signal/2).
