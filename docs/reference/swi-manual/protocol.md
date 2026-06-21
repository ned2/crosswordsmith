
## 4.38 Creating a Protocol of the User Interaction

SWI-Prolog offers the possibility to log the interaction with the user on a file.^(162A similar facility was added to Edinburgh C-Prolog by Wouter Jansweijer.) All Prolog interaction, including warnings and tracer output, are written to the protocol file.

**protocol**(`+File`)  
Start protocolling on file `File`. If there is already a protocol file open, then close it first. If `File` exists it is truncated.

**protocola**(`+File`)  
Equivalent to [protocol/1](protocol.html#protocol/1), but does not truncate the `File` if it exists.

**noprotocol**  
Stop making a protocol of the user interaction. Pending output is flushed on the file.

**protocolling**(`-File`)  
True if a protocol was started with [protocol/1](protocol.html#protocol/1) or [protocola/1](protocol.html#protocola/1) and unifies `File` with the current protocol output file.
