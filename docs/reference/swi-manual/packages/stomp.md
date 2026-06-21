STOMP -- a SWI-Prolog STOMP client

Hongxin Liang and Jan Wielemaker  
SWI-Prolog Solutions b.v.  
The Netherlands  
E-mail: [jan@swi-prolog.org](mailto:jan@swi-prolog.org)

Abstract

This package provide a client to STOMP, The *Simple Text Oriented Messaging Protocol*. See the [STOMP specification](http://stomp.github.io/) for details. We used STOMP to connect to the [RabbitMQ](https://www.rabbitmq.com/) message broker.

# Table of Contents

[1 library(stomp): STOMP client.](#sec:1)

[1.1 Threading](#sec:1.1)

[1.2 Reconnecting](#sec:1.2)

## 1 library(stomp): STOMP client.

author  
Hongxin Liang and Jan Wielemaker

See also  
\- [http://stomp.github.io/index.html](http://stomp.github.io/index.html)  
- [https://github.com/jasonrbriggs/stomp.py](https://github.com/jasonrbriggs/stomp.py)

license  
BSD-2

To be done  
TSL support

This module provides a STOMP (Simple (or Streaming) Text Orientated Messaging Protocol) client. This client is based on work by Hongxin Liang. The current version is a major rewrite, both changing the API and the low-level STOMP frame (de)serialization.

The predicate [stomp_connection/5](#stomp_connection/5) is used to register a connection. The connection is established by [stomp_connect/1](#stomp_connect/1), which is lazily called from any of the predicates that send a STOMP frame. After establishing the connection two threads are created. One receives STOMP frames and the other manages and watches the *heart beat*.

### 1.1 Threading

Upon receiving a frame the callback registered with [stomp_connection/5](#stomp_connection/5) is called in the context of the receiving thread. More demanding applications may decide to send incomming frames to a SWI-Prolog message queue and have one or more *worker threads* processing the input. Alternatively, frames may be inspected by the receiving thread and either processed immediately or be dispatched to either new or running threads. The best scenario depends on the message rate, processing time and concurrency of the Prolog application.

All message sending predicates of this library are *thread safe*. If two threads send a frame to the same connection the library ensures that both frames are properly serialized.

### 1.2 Reconnecting

By default this library tries to establish the connection and the user gets notified by means of a `disconnected` pseudo frame that the connection is lost. Using the Options argument of [stomp_connection/6](#stomp_connection/6) the system can be configured to try and keep connecting if the server is not available and reconnect if the connection is lost. See the `pong.pl` example.

\[det\]**stomp_connection**(`+Address, +Host, +Headers, :Callback, -Connection`)  
\[det\]**stomp_connection**(`+Address, +Host, +Headers, :Callback, -Connection, +Options`)  
Create a connection reference. The connection is not set up yet by this predicate. `Callback` is called on any received frame except for *heart beat* frames as below.

``` code
call(Callback, Command, Connection, Header, Body)
```

Where command is one of the commands below. `Header` is a dict holding the STOMP frame header, where all values are strings except for the `'content-length'` key value which is passed as an integer.

Body is a string or, if the data is of the type `application/json`, a dict.

**connected**  
A connection was established. `Connection` and Header are valid.

**disconnected**  
The connection was lost. Only `Connection` is valid.

**message**  
A message arrived. All three arguments are valid. Body is a dict if the `content-type` of the message is `application/json` and a string otherwise.

**heartbeat**  
A heartbeat was received. Only `Connection` is valid. This callback is also called for each newline that follows a frame. These newlines can be a heartbeat, but can also be additional newlines that follow a frame.

**heartbeat_timeout**  
No heartbeat was received. Only `Connection` is valid.

**error**  
An error happened. All three arguments are valid and handled as `message`.

Note that [stomp_teardown/1](#stomp_teardown/1) causes the receiving and heartbeat thread to be signalled with abort/0.

`Options` processed:

**reconnect**(`+Bool`)  
Try to reestablish the connection to the server if the connection is lost. Default is `false`

**connect_timeout**(`+Seconds`)  
Maximum time to try and reestablish a connection. The default is `600` (10 minutes).

**json_options**(`+Options`)  
`Options` passed to json_read_dict/3 to translate `application/json` content to Prolog. Default is `[]`.

|  |  |
|----|----|
| `Address` | is a valid address for tcp_connect/3. Normally a term `Host`:Port, e.g., `localhost:32772`. |
| `Host` | is a path denoting the STOMP host. Often just `/`. |
| `Headers` | is a dict with STOMP headers used for the `CONNECT` request. |
| `Connection` | is an opaque ground term that identifies the connection. |

\[nondet\]**stomp_connection_property**(`?Connection, ?Property`)  
True when `Property`, is a property of `Connection`. Defined properties are:

**address**(`Address`)  
**callback**(`Callback`)  
**host**(`Host`)  
**headers**(`Headers`)  
**reconnect**(`Bool`)  
**connect_timeout**(`Seconds`)  
All the above properties result from the [stomp_connection/6](#stomp_connection/6) registration.

**receiver_thread_id**(`Thread`)  
**stream**(`Stream`)  
**heartbeat_thread_id**(`Thread`)  
**received_heartbeat**(`TimeStamp`)  
These describe an active STOMP connection.

**stomp_destroy_connection**(`+Connection`)  
Destroy a connection. If it is active, first use [stomp_teardown/1](#stomp_teardown/1) to disconnect.

\[det\]**stomp_setup**(`+Connection, +Options`)  
Set up the actual socket connection and start receiving thread. This is a no-op if the connection has already been created. The `Options` processed are below. Other options are passed to tcp_connect/3.

**timeout**(`+Seconds`)  
If tcp_connect/3 fails, retry until the timeout is reached. If `Seconds` is `inf` or `infinite`, keep retrying forever.

\[semidet\]**stomp_teardown**(`+Connection`)  
Tear down the socket connection, stop receiving thread and heartbeat thread (if applicable). The registration of the connection as created by [stomp_connection/5](#stomp_connection/5) is preserved and the connection may be reconnected using [stomp_connect/1](#stomp_connect/1).

\[det\]**stomp_reconnect**(`+Connection`)  
Teardown the connection and try to reconnect.

\[det\]**stomp_connect**(`+Connection`)  
\[det\]**stomp_connect**(`+Connection, +Options`)  
Setup the connection. First ensures a socket connection and if this is new send a `CONNECT` frame. Protocol version and heartbeat negotiation will be handled. `STOMP` frame is not used for backward compatibility.

This predicate waits for the connection handshake to have completed. Currently it waits for a maximum of 10 seconds after establishing the socket for the server reply.

Calling this on an established connection has no effect.

Errors  
`stomp_error(connect, Connection, Detail)` if no connection could be established.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\CONNECT_or_STOMP_Frame).](http://stomp.github.io/stomp-specification-1.2.html\#CONNECT_or_STOMP_Frame).)

\[det\]**stomp_send**(`+Connection, +Destination, +Headers, +Body`)  
Send a `SEND` frame. If `content-type` is not provided, `text/plain` will be used. `content-length` will be filled in automatically.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\SEND](http://stomp.github.io/stomp-specification-1.2.html\#SEND)

\[det\]**stomp_send_json**(`+Connection, +Destination, +Headers, +JSON`)  
Send a `SEND` frame. `JSON` can be either a `JSON` term or a dict. `content-type` is filled in automatically as `application/json` and `content-length` will be filled in automatically as well.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\SEND](http://stomp.github.io/stomp-specification-1.2.html\#SEND)

\[det\]**stomp_subscribe**(`+Connection, +Destination, +Id, +Headers`)  
Send a `SUBSCRIBE` frame.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\SUBSCRIBE](http://stomp.github.io/stomp-specification-1.2.html\#SUBSCRIBE)

\[det\]**stomp_unsubscribe**(`+Connection, +Id`)  
Send an `UNSUBSCRIBE` frame.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\UNSUBSCRIBE](http://stomp.github.io/stomp-specification-1.2.html\#UNSUBSCRIBE)

\[det\]**stomp_ack**(`+Connection, +MessageId, +Headers`)  
Send an `ACK` frame. See [stomp_ack/2](#stomp_ack/2) for simply passing the header received with the message we acknowledge.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\ACK](http://stomp.github.io/stomp-specification-1.2.html\#ACK)

\[det\]**stomp_nack**(`+Connection, +MessageId, +Headers`)  
Send a `NACK` frame. See [stomp_nack/2](#stomp_nack/2) for simply passing the header received with the message we acknowledge.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\NACK](http://stomp.github.io/stomp-specification-1.2.html\#NACK)

\[det\]**stomp_ack**(`+Connection, +MsgHeader`)  
\[det\]**stomp_nack**(`+Connection, +MsgHeader`)  
Reply with an ACK or NACK based on the received message header. On a STOMP 1.1 request we get an `ack` field in the header and reply with an `id`. For STOMP 1.2 we reply with the `message-id` and `subscription` that we received with the message.

\[det\]**stomp_begin**(`+Connection, +Transaction`)  
Send a `BEGIN` frame. @see [http://stomp.github.io/stomp-specification-1.2.html\\BEGIN](http://stomp.github.io/stomp-specification-1.2.html\#BEGIN)

\[det\]**stomp_commit**(`+Connection, +Transaction`)  
Send a `COMMIT` frame.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\COMMIT](http://stomp.github.io/stomp-specification-1.2.html\#COMMIT)

\[det\]**stomp_abort**(`+Connection, +Transaction`)  
Send a `ABORT` frame.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\ABORT](http://stomp.github.io/stomp-specification-1.2.html\#ABORT)

\[semidet\]**stomp_transaction**(`+Connection, :Goal`)  
Run `Goal` as once/1, tagging all `SEND` messages inside the transaction with the transaction id. If `Goal` fails or raises an exception the transaction is aborted. Failure or exceptions cause the transaction to be aborted using [stomp_abort/2](#stomp_abort/2), after which the result is forwarded.

\[det\]**stomp_disconnect**(`+Connection, +Headers`)  
Send a `DISCONNECT` frame. If the connection has the `reconnect` property set to `true`, this will be reset to `disconnected` to avoid reconnecting. A subsequent [stomp_connect/2](#stomp_connect/2) resets the reconnect status.

See also  
[http://stomp.github.io/stomp-specification-1.2.html\\DISCONNECT](http://stomp.github.io/stomp-specification-1.2.html\#DISCONNECT)

# Index

?  
[stomp_abort/2](#stomp_abort/2)  
[stomp_ack/2](#stomp_ack/2)  
[stomp_ack/3](#stomp_ack/3)  
[stomp_begin/2](#stomp_begin/2)  
[stomp_commit/2](#stomp_commit/2)  
[stomp_connect/1](#stomp_connect/1)  
[stomp_connect/2](#stomp_connect/2)  
[stomp_connection/5](#stomp_connection/5)  
[stomp_connection/6](#stomp_connection/6)  
[stomp_connection_property/2](#stomp_connection_property/2)  
[stomp_destroy_connection/1](#stomp_destroy_connection/1)  
[stomp_disconnect/2](#stomp_disconnect/2)  
[stomp_nack/2](#stomp_nack/2)  
[stomp_nack/3](#stomp_nack/3)  
[stomp_reconnect/1](#stomp_reconnect/1)  
[stomp_send/4](#stomp_send/4)  
[stomp_send_json/4](#stomp_send_json/4)  
[stomp_setup/2](#stomp_setup/2)  
[stomp_subscribe/4](#stomp_subscribe/4)  
[stomp_teardown/1](#stomp_teardown/1)  
[stomp_transaction/2](#stomp_transaction/2)  
[stomp_unsubscribe/2](#stomp_unsubscribe/2)  
