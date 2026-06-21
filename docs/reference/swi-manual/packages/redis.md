Redis -- a SWI-Prolog client for redis

Jan Wielemaker  
SWI-Prolog Solutions b.v.  
The Netherlands  
E-mail: [jan@swi-prolog.org](mailto:jan@swi-prolog.org)

Abstract

This package provides the library `redis.pl` that implements a client for the \[redis\](https://redis.io/), *” An open source (BSD licensed), in-memory data structure store, used as a database, cache and message broker’*

# Table of Contents

[1 About the SWI-Prolog Redis client](#sec:1)

[1.1 Redis and threads](#sec:1.1)

[1.2 Redis TLS support](#sec:1.2)

[1.3 Using Redis sentinels](#sec:1.3)

[1.4 About versions](#sec:1.4)

[1.5 Redis as a message brokering system](#sec:1.5)

[1.6 History](#sec:1.6)

[2 library(redis): Redis client](#sec:2)

[3 library(redis_streams): Using Redis streams](#sec:3)

## 1 About the SWI-Prolog Redis client

[Redis](https://redis.io) is an in-memory key-value store. Redis can be operated as a simple store for managing (notably) volatile persistent data. Redis can operate in serveral modes, ranging from a single server to clusters organised in several different ways to provide high availability, resilience and replicating the data close to the involved servers. In addition to being a key-value store, Redis enables additional communication between clients such as *publish/subscribe* to message channels, *streams*, etc.

These features can be used to connect *micro services*, both for sharing state, notifications and distributing tasks.

### 1.1 Redis and threads

The connection between the redis client and server uses a *stream pair*. Although SWI-Prolog stream I/O is thread-safe, having multiple threads using this same connection will mixup writes and their replies.

At the moment, the following locking is in place.

- Connections created using [redis_connect/3](#redis_connect/3) are *not* locked. This implies the connection handle may be used from a single thread only, or [redis/3](#redis/3) requests must be protected using with_mutex/2.
- Redis/3 request using a *server name* established using [redis_server/3](#redis_server/3) are locked using a mutex with the same name as the server name.

### 1.2 Redis TLS support

If SWI-Prolog includes the `ssl` library, the Redis client can connect to the server using TLS (SSL). Connecting requires the same three files as `redis-cli` requires: the root certificate file, a client certificate and the private key of the client certificate. Below is an example call to [redis_server/3](#redis_server/3):

``` code
:- redis_server(swish, localhost:6379,
                [ user(bob),
                  password("topsecret"),
                  version(3),
                  tls(true),
                  cacert('ca.crt'),
                  key('client.key'),
                  cert('client.cert')
                ]).
```

### 1.3 Using Redis sentinels

Redis sentinels is one of the two options to create a high availability service. It consists of minimally three Redis servers and mininally three sentinel servers. The sentinel servers monitor the Redis servers and will initiate a fail-over when the master becomes disfunctional and certain safety constraints are satisfied. A client needs to be aware of this setup. It is given an initial list with (a subset of) the known sentinels. The client attempts to connect to one of the sentinels and ask it for the current Redis master server. Details are described in [Sentinel client spec](https://redis.io/docs/reference/sentinel-clients/). The SWI-Prolog client maintains the actual list of sentinels dynamically after successful discovery of the first sentinel. Below is an example [redis_server/3](#redis_server/3) to connect to a sentinel network. The *Address* specification `sentinel(swish)` tells the library we want to connect to a sentinel network that is monitored under the name `swish`.

``` code
:- redis_server(swish_sentinel, sentinel(swish),
                [ user(janbob),
                  password("topsecret"),
                  version(3),
                  sentinels([ host1:26379,
                                  host2:26379,
                                  host3:26379
                                ])
                ]).
```

### 1.4 About versions

The current stable version of Redis is 6. Many Linux distros still ship with version 5. Both talk protocol version 2. Version 6 also supports protocol version 3. The main differences are:

- The version 3 protocol has several improvements that notably improvement passing large objects using a *streaming* protocol.
- Hashes (maps) in the version 3 protocol are exchanged as lists of *pairs* (`Name-Value`), while version 2 exchanges hashes as a list of alternating names and values. This is visible to the user. High level predicates such as [redis_get_hash/3](#redis_get_hash/3) deal with both representations.
- The version 3 protocol supports *push messages* to deal with *monitor* and *subscribe* events on the same connection as used for handling normal requests.

New projects are encouraged to use Redis version 6 with the version 3 protocol. See [redis_server/3](#redis_server/3).

### 1.5 Redis as a message brokering system

Starting with Redis 5, redis supports *streams*. A stream is a list of messages. Streams can be used as a reliable alternative to the older Redis PUB/SUB (Publish Subscribe) mechanism that has no memory, i.e., if a node is down when a message arrives the message is missed. In addition, they can be used to have each message processed by a *consumer* that belongs to a *consumer group*. Both facilities are supported by `library(redis_streams)` ([section 3](#sec:3))

Redis streams provide all the low-level primitives to realise message brokering. Putting it all together is non-trivial though. Notably:

- We must take care of messages that have been sent to some consumer but the consumer fails to process the message and (thus) ACK it is processed. This is handled by [xlisten_group/5](#xlisten_group/5) using several options. Good defaults for these options are hard to give as it depends on the required processing time for a message, how common failures are and an acceptable delay time in case of a failure, what to do in case of a persistent failure, etc.
- Streams are independent from consumer groups and acknowledged messages remain in the stream. [xstream_set/3](#xstream_set/3) can be used to limit the length of the stream, discarding the oldest messages. However, it is hard to give a sensible default. The required queue length depends on the the size of the messages, whether messages come in more or less randomly or in bursts (that cause the stream to grow for a while), available memory, how bad it is if some messages get lost, etc.

The directory `doc/packages/examples/redis` in the installation provides an example using streams and consumer groups to realise one or more clients connected to one or more compute nodes.

### 1.6 History

This module is based on the `gpredis.pl` by Sean Charles for GNU-Prolog. This file greatly helped me understanding what had to be done, although, eventually, not much of the original interface is left. The main difference to the original client are:

- Replies are not wrapped by type in a compound term.
- String replies use the SWI-Prolog string type.
- Values can be specified as `Value as prolog`, after which they are returns as a (copy of) Value. This prefixes the value using "`\`u0000T`\`u0000".
- Strings are in UTF-8 encoding to support full Unicode.
- Using [redis_server/3](#redis_server/3), actual connections are established lazily and when a connection is lost it is automatically restarted.
- This library allows for using the Redis publish/subscribe interface. Messages are propagated using broadcast/1.

## 2 library(redis): Redis client

This library is a client to [Redis](https://redis.io), a popular key value store to deal with caching and communication between micro services.

In the typical use case we register the details of one or more Redis servers using [redis_server/3](#redis_server/3). Subsequenly, [redis/2](#redis/2)-3 is used to issue commands on the server. For example:

``` code
?- redis_server(default, redis:6379, [password("secret")]).
?- redis(default, set(user, "Bob")).
?- redis(default, get(user), User).
User = "Bob"
```

\[det\]**redis_server**(`+ServerName, +Address, +Options`)  
Register a redis server without connecting to it. The `ServerName` acts as a lazy connection alias. Initially the `ServerName` `default` points at `localhost:6379` with no connect options. The `default` server is used for [redis/1](#redis/1) and [redis/2](#redis/2) and may be changed using this predicate. `Options` are described with [redis_connect/3](#redis_connect/3).

Connections established this way are by default automatically reconnected if the connection is lost for some reason unless a `reconnect(false)` option is specified.

\[det\]**redis_connect**(`-Connection`)  
\[det\]**redis_connect**(`+Address, -Connection, +Options`)  
\[det\]**redis_connect**(`-Connection, +Host, +Port`)  
Connect to a redis server. The main mode is `redis_connect(+Address, -Connection, +Options)`. [redis_connect/1](#redis_connect/1) is equivalent to `redis_connect(localhost:6379, Connection, [])`. `Options`:

**reconnect**(`+Boolean`)  
If `true`, try to reconnect to the service when the connection seems lost. Default is `true` for connections specified using [redis_server/3](#redis_server/3) and `false` for explictly opened connections.

**user**(`+User`)  
If `version(3)` and `password(Password)` are specified, these are used to authenticate using the `HELLO` command.

**password**(`+Password`)  
Authenticate using `Password`

**version**(`+Version`)  
Specify the connection protocol version. Initially this is version 2. Redis 6 also supports version 3. When specified as `3`, the `HELLO` command is used to upgrade the protocol.

**tls**(`true`)  
When specified, initiate a TLS connection. If this option is specified we must also specify the `cacert`, `key` and `cert` options.

**cacert**(`+File`)  
CA Certificate file to verify with.

**cert**(`+File`)  
Client certificate to authenticate with.

**key**(`+File`)  
Private key file to authenticate with.

**sentinels**(`+ListOfAddresses`)  
Used together with an `Address` of the form `sentinel(MasterName)` to enable contacting a network of Redis servers guarded by a sentinel network.

**sentinel_user**(`+User`)  
**sentinel_password**(`+Password`)  
Authentication information for the senitels. When omitted we try to connect withour authentication.

Instead of using these predicates, [redis/2](#redis/2) and [redis/3](#redis/3) are normally used with a *server name* argument registered using [redis_server/3](#redis_server/3). These predicates are meant for creating a temporary paralel connection or using a connection with a *blocking* call.

|  |  |
|----|----|
| `Address` | is a term `Host`:`Port`, `unix(File)` or the name of a server registered using [redis_server/3](#redis_server/3). The latter realises a *new* connection that is typically used for blocking redis commands such as listening for published messages, waiting on a list or stream. |

Compatibility  
`redis_connect(-Connection, +Host, +Port)` provides compatibility to the original GNU-Prolog interface and is equivalent to `redis_connect(Host:Port, Connection, [])`.

\[semidet\]**tls_verify**(`+SSL, +ProblemCert, +AllCerts, +FirstCert, +Status`)  
Accept or reject the certificate verification. Similar to the Redis command line client (`redis-cli`), we accept the certificate as long as it is signed, not verifying the hostname.

\[nondet\]**sentinel_slave**(`+ServerId, +Pool, -Slave, +Options`)  
True when `Slave` is a slave server in the sentinel cluster. `Slave` is a dict holding the keys and values as described by the Redis command

``` code
SENTINEL SLAVES mastername
```

\[det\]**redis_disconnect**(`+Connection`)  
\[det\]**redis_disconnect**(`+Connection, +Options`)  
Disconnect from a redis server. The second form takes one option, similar to close/2:

**force**(`Force`)  
When `true` (default `false`), do not raise any errors if `Connection` does not exist or closing the connection raises a network or I/O related exception. This version is used internally if a connection is in a broken state, either due to a protocol error or a network issue.

\[semidet\]**redis**(`+Connection, +Request`)  
This predicate is overloaded to handle two types of requests. First, it is a shorthand for `redis(Connection, Command, _)` and second, it can be used to exploit Redis *pipelines* and *transactions*. The second form is acticated if `Request` is a *list*. In that case, each element of the list is either a term `Command -> Reply` or a simple `Command`. Semantically this represents a sequence of [redis/3](#redis/3) and [redis/2](#redis/2) calls. It differs in the following aspects:

- All commands are sent in one batch, after which all replies are read. This reduces the number of *round trips* and typically greatly improves performance.
- If the first command is `multi` and the last `exec`, the commands are executed as a Redis *transaction*, i.e., they are executed *atomically*.
- If one of the commands returns an error, the subsequent commands **are still executed**.
- You can not use variables from commands earlier in the list for commands later in the list as a result of the above execution order.

Procedurally, the process takes the following steps:

1.  Send all commands
2.  Read all replies and push messages
3.  Handle all callbacks from push messages
4.  Check whether one of the replies is an error. If so, raise this error (subsequent errors are lost)
5.  Bind all replies for the `Command -> Reply` terms.

Examples

``` code
?- redis(default,
         [ lpush(li,1),
           lpush(li,2),
           lrange(li,0,-1) -> List
         ]).
List = ["2", "1"].
```

\[semidet\]**redis**(`+Connection, +Command, -Reply`)  
Execute a redis `Command` on Connnection. Next, bind `Reply` to the returned result. `Command` is a callable term whose functor is the name of the Redis command and whose arguments are translated to Redis arguments according to the rules below. Note that all text is always represented using UTF-8 encoding.

- Atomic values are emitted verbatim
- A term A:B:... where all arguments are either atoms, strings or integers (**no floats**) is translated into a string `"A:B:..."`. This is a common shorthand for representing Redis keys.
- A term Term as prolog is emitted as "`\`u0000T`\`u0000" followed by Term in canonical form.
- Any other term is emitted as write/1.

`Reply` is either a plain term (often a variable) or a term `Value as Type`. In the latter form, `Type` dictates how the Redis *bulk* reply is translated to Prolog. The default equals to `auto`, i.e., as a number of the content satisfies the Prolog number syntax and as an atom otherwise.

- `status(Atom)` Returned if the server replies with `+ Status`. Atom is the textual value of `Status` converted to lower case, e.g., `status(ok)` or `status(pong)`.
- `nil` This atom is returned for a NIL/NULL value. Note that if the reply is only `nil`, [redis/3](#redis/3) *fails*. The `nil` value may be embedded inside lists or maps.
- A number Returned if the server replies an integer (":Int"), double (",Num") or big integer ("(Num")
- A string Returned on a *bulk* reply. Bulk replies are supposed to be in UTF-8 encoding. The the bulk reply starts with "`\`u0000T`\`u0000" it is supposed to be a Prolog term. Note that this intepretation means it is **not** possible to read arbitrary binary blobs.
- A list of replies. A list may also contain `nil`. If `Reply` as a whole would be `nil` the call fails.
- A list of *pairs*. This is returned for the redis version 3 protocol "%Map". Both the key and value respect the same rules as above.

Redis *bulk* replies are translated depending on the `as` `Type` as explained above.

**string**  
**string**(`Encoding`)  
Create a SWI-Prolog string object interpreting the blob as following `Encoding`. `Encoding` is a restricted set of SWI-Prolog's encodings: `bytes` (`iso_latin_1`), `utf8` and `text` (the current locale translation).

**atom**  
**atom**(`Encoding`)  
As above, producing an atom.

**codes**  
**codes**(`Encoding`)  
As above, producing a list of integers (Unicode code points)

**chars**  
**chars**(`Encoding`)  
As above, producing a list of one-character atoms.

**integer**  
**float**  
**rational**  
**number**  
Interpret the bytes as a string representing a number. If the string does not represent a number of the requested type a `type_error(Type, String)` is raised.

**tagged_integer**  
Same as integer, but demands the value to be between the Prolog flags `min_tagged_integer` and `max_tagged_integer`, allowing the value to be used as a dict key.

**auto**  
Same as `auto(atom, number)`

**auto**(`AsText, AsNumber`)  
If the bulk string confirms the syntax of `AsNumber`, convert the value to the requested numberical type. Else convert the value to text according to `AsText`. This is similar to the Prolog predicate name/2.

**dict_key**  
Alias for `auto(atom,tagged_integer)`. This allows the value to be used as a key for a SWI-Prolog dict.

**pairs**(`AsKey, AsValue`)  
Convert a map or array of even length into pairs for which the key satisfies `AsKey` and the value `AsValue`. The `pairs` type can also be applied to a Redis array. In this case the array length must be even. This notably allows fetching a Redis *hash* as pairs using `HGETALL` using version 2 of the Redis protocol.

**dict**(`AsKey, AsValue`)  
Similar to `pairs(AsKey, AsValue)`, but convert the resulting pair list into a SWI-Prolog dict. `AsKey` must convert to a valid dict key, i.e., an atom or tagged integer. See `dict_key`.

**dict**(`AsValue`)  
Shorthand for `dict(dict_key, AsValue)`.

Here are some simple examples

``` code
?- redis(default, set(a, 42), X).
X = status("OK").
?- redis(default, get(a), X).
X = "42".
?- redis(default, get(a), X as integer).
X = 42.
?- redis(default, get(a), X as float).
X = 42.0.
?- redis(default, set(swipl:version, 8)).
true.
?- redis(default, incr(swipl:version), X).
X = 9.
```

Errors  
`redis_error(Code, String)`

**redis**(`+Request`)  
Connect to the default redis server, call redist/3 using `Request`, disconnect and print the result. This predicate is intended for interactive usage.

\[det\]**redis_write**(`+Redis, +Command`)  
\[det\]**redis_read**(`+Redis, -Reply`)  
Write command and read replies from a `Redis` server. These are building blocks for subscribing to event streams.

\[det\]**redis_get_list**(`+Redis, +Key, -List`)  
\[det\]**redis_get_list**(`+Redis, +Key, +ChunkSize, -List`)  
Get the content of a `Redis` list in `List`. If `ChunkSize` is given and smaller than the list length, `List` is returned as a *lazy list*. The actual values are requested using redis `LRANGE` requests. Note that this results in O(N`^`2) complexity. Using a lazy list is most useful for relatively short lists holding possibly large items.

Note that values retrieved are *strings*, unless the value was added using `Term as prolog`.

It seems possible for `LLEN` to return `OK`. I don't know why. As a work-around we return the empty list rather than an error.

See also  
lazy_list/2 for a discussion on the difference between lazy lists and normal lists.

\[det\]**redis_set_list**(`+Redis, +Key, +List`)  
Associate a `Redis` key with a list. As `Redis` has no concept of an empty list, if `List` is `[]`, `Key` is *deleted*. Note that key values are always strings in `Redis`. The same conversion rules as for [redis/1](#redis/1)-3 apply.

\[det\]**redis_get_hash**(`+Redis, +Key, -Data:dict`)  
\[det\]**redis_set_hash**(`+Redis, +Key, +Data:dict`)  
Put/get a `Redis` hash as a Prolog dict. Putting a dict first deletes `Key`. Note that in many cases applications will manage `Redis` hashes by key. [redis_get_hash/3](#redis_get_hash/3) is notably a user friendly alternative to the `Redis` `HGETALL` command. If the `Redis` hash is not used by other (non-Prolog) applications one may also consider using the `Term as prolog` syntax to store the Prolog dict as-is.

\[det\]**redis_array_dict**(`?Array, ?Tag, ?Dict`)  
Translate a Redis reply representing hash data into a SWI-Prolog dict. `Array` is either a list of alternating keys and values or a list of *pairs*. When translating to an array, this is always a list of alternating keys and values.

|       |                             |
|-------|-----------------------------|
| `Tag` | is the SWI-Prolog dict tag. |

\[det\]**redis_scan**(`+Redis, -LazyList, +Options`)  
\[det\]**redis_sscan**(`+Redis, +Set, -LazyList, +Options`)  
\[det\]**redis_hscan**(`+Redis, +Hash, -LazyList, +Options`)  
\[det\]**redis_zscan**(`+Redis, +Set, -LazyList, +Options`)  
Map the `Redis` `SCAN`, `SSCAN`, `HSCAN` and `ZSCAN`‘commands into a *lazy list*. For [redis_scan/3](#redis_scan/3) and [redis_sscan/4](#redis_sscan/4) the result is a list of strings. For [redis_hscan/4](#redis_hscan/4) and [redis_zscan/4](#redis_zscan/4), the result is a list of *pairs*. `Options` processed:

**match**(`Pattern`)  
Adds the `MATCH` subcommand, only returning matches for `Pattern`.

**count**(`Count`)  
Adds the `COUNT` subcommand, giving a hint to the size of the chunks fetched.

**type**(`Type`)  
Adds the `TYPE` subcommand, only returning answers of the indicated type.

See also  
lazy_list/2.

\[nondet\]**redis_current_command**(`+Redis, ?Command`)  
\[nondet\]**redis_current_command**(`+Redis, ?Command, -Properties`)  
True when `Command` has `Properties`. Fails if `Command` is not defined. The [redis_current_command/3](#redis_current_command/3) version returns the command argument specification. See `Redis` documentation for an explanation.

\[nondet\]**redis_property**(`+Redis, ?Property`)  
True if `Property` is a property of the `Redis` server. Currently uses `redis(info, String)` and parses the result. As this is for machine usage, properties names \*\_human are skipped.

\[det\]**redis_subscribe**(`+Redis, +Channels, -Id, +Options`)  
Subscribe to one or more `Redis` PUB/SUB channels. This predicate creates a thread using thread_create/3 with the given `Options`. Once running, the thread listens for messages. The message content is a string or Prolog term as described in [redis/3](#redis/3). On receiving a message, the following message is broadcasted:

``` code
redis(Id, Channel, Data)
```

If [redis_unsubscribe/2](#redis_unsubscribe/2) removes the last subscription, the thread terminates.

To simply print the incomming messages use e.g.

``` code
?- listen(redis(_, Channel, Data),
          format('Channel ~p got ~p~n', [Channel,Data])).
true.
?- redis_subscribe(default, test, Id, []).
Id = redis_pubsub_3,
?- redis(publish(test, "Hello world")).
Channel test got "Hello world"
1
true.
```

|  |  |
|----|----|
| `Id` | is the thread identifier of the listening thread. Note that the `Options` `alias(Name)` can be used to get a system wide name. |

\[det\]**redis_subscribe**(`+Id, +Channels`)  
\[det\]**redis_unsubscribe**(`+Id, +Channels`)  
Add/remove channels from for the subscription. If no subscriptions remain, the listening thread terminates.

|  |  |
|----|----|
| `Channels` | is either a single channel or a list thereof. Each channel specification is either an atom or a term‘A:B:...\`, where all parts are atoms. |

**redis_current_subscription**(`?Id, ?Channels`)  
True when a PUB/SUB subscription with `Id` is listening on `Channels`.

## 3 library(redis_streams): Using Redis streams

See also  
[https://redis.io/topics/streams-intro](https://redis.io/topics/streams-intro)

A Redis stream is a set of messages consisting of key-value pairs that are identified by a time and sequence number. Streams are powerful objects that can roughly be used for three purposes:

- Maintain and query a log of events, i.e., a *timeline*.
- Provide an alternative to Redis’publish/subscribe API that ensures messages get delivered by all clients even if they are offline at the moment an event is published.
- Distribute messages over a group of clients. This mode assigns messages to clients in a round-robin fashion. Clients confirm a specific message is handled. Living clients can inspect the stream for possibly dead clients and migrate the pending messages to other clients.

This library abstracts the latter two scenarios. The main predicates are

- [xadd/4](#xadd/4) to add to a stream
- [xlisten/3](#xlisten/3) to read and broadcast messages from a stream
- [xlisten_group/5](#xlisten_group/5) to act as a *consumer* in a consumer group.

**xstream_set**(`+Redis, +Key, +Option`)  
Set an option on for `Key` on `Redis`. Currently supports:

**maxlen**(`+Count`)  
Make [xadd/4](#xadd/4) add a `MAXLEN ~` `Count` option to the `XADD` command, capping the length of the stream. See also `Redis` as a message brokering system ([section 1.5](#sec:1.5))

\[det\]**xadd**(`+Redis, +Key, ?Id, +Data:dict`)  
Add a message to a the stream `Key` on `Redis`. The length of the stream can be capped using the [xstream_set/3](#xstream_set/3) option `maxlen(Count)`. If `Id` is unbound, generating the id is left to the server and `Id` is unified with the returned id. The returned id is a string consisting of the time stamp in milliseconds and a sequence number. See `Redis` docs for details.

**xlisten**(`+Redis, +Streams, +Options`)  
Listen using `XREAD` on one or more `Streams` on the server `Redis`. For each message that arrives, call broadcast/1, where Data is a dict representing the message.

``` code
broadcast(redis(Redis, Stream, Id, Data))
```

`Options`:

**count**(`+Count`)  
Process at most `Count` messages per stream for each request.

**start**(`+Start`)  
Normally either `0` to start get all messages from the epoch or `$` to get messages starting with the last. Default is `$`.

**starts**(`+List`)  
May be used as an alternative to the start/1 option to specify the start for each stream. This may be used to restart listening if the application remembers the last processed id.

Note that this predicate does **not terminate**. It is normally executed in a thread. The following call listens to the streams `key1` and `key2` on the default `Redis` server. Using `reconnect(true)`, the client will try to re-establish a connection if the collection got lost.

``` code
?- redis_connect(default, C, [reconnect(true)]),
   thread_create(xlisten(C, [key1, key2], [start($)]),
                 _, [detached(true)]).
```

|  |  |
|----|----|
| `Redis` | is either a `Redis` server name (see [redis_server/3](#redis_server/3)) or an open connection. If it is a server name, a new connection is opened that is closed if [xlisten/3](#xlisten/3) completes. |

See also  
[redis_subscribe/2](#redis_subscribe/2) implements the classical pub/sub system of `Redis` that does not have any memory.

**xlisten_group**(`+Redis, +Group, +Consumer, +Streams, +Options`)  
Listen as `Consumer` to `Group`. This is similar to [xlisten/3](#xlisten/3), with the following differences:

- Instead of using broadcast/1, broadcast_request/1 is used and the message is only considered processed if broadcast_request/1 succeeds. If the message is handled with success, an `XACK` is sent to the server.

`Options` processed:

**block**(`+Seconds`)  
Causes `XREADGROUP` to return with timeout when no messages arrive within `Seconds`. On a timeout, xidle_group/5 is called which will try to handle messages to other consumers pending longer than `Seconds`. Choosing the time depends on the application. Notably:

- Using a time shorter than the required processing time will make the job migrate from consumer to consumer until `max_deliveries(Count)` is exceeded. Note that the original receiver does not notice that the job is claimed and thus multiple consumers may ultimately answer the message.
- Using a too long time causes an unnecessarily long delay if a node fails.

**max_deliveries**(`+Count`)  
Re-deliver (using `XCLAIM`) a message max `Count` times. Exceeding this calls [xhook/2](#xhook/2). Default `Count` is `3`.

**max_claim**(`+Count`)  
Do not claim more than `Count` messages during a single idle action. Default is `10`.

**xconsumer_stop**(`+Leave`)  
May be called from a consumer listener to stop the consumer. This predicate throws the exception `redis(stop(Leave))`, which is caught by [xlisten_group/5](#xlisten_group/5).

\[multifile\]**xhook**(`+Stream, +Event`)  
This multifile predicate is called on certain stream events. Defined events are:

**delivery_failed**(`Id, Group, Delivered`)  
A message was delivered more than specified by max_deliveries/1 of [xlisten_group/5](#xlisten_group/5). `Id` is the message id, `Group` the group and `Delivered` the current delivery count. If the hooks fails, the message is acknowledged using `XACK`. From [introduction to streams](https://redis.io/topics/streams-intro):

> "So once the deliveries counter reaches a given large number that you chose, it is probably wiser to put such messages in another stream and send a notification to the system administrator. This is basically the way that Redis streams implement the concept of the dead letter."

# Index

?  
[redis/1](#redis/1)  
[redis/2](#redis/2)  
[redis/3](#redis/3)  
[redis_array_dict/3](#redis_array_dict/3)  
[redis_connect/1](#redis_connect/1)  
[redis_connect/3](#redis_connect/3)  
[redis_current_command/2](#redis_current_command/2)  
[redis_current_command/3](#redis_current_command/3)  
[redis_current_subscription/2](#redis_current_subscription/2)  
[redis_disconnect/1](#redis_disconnect/1)  
[redis_disconnect/2](#redis_disconnect/2)  
[redis_get_hash/3](#redis_get_hash/3)  
[redis_get_list/3](#redis_get_list/3)  
[redis_get_list/4](#redis_get_list/4)  
[redis_hscan/4](#redis_hscan/4)  
[redis_property/2](#redis_property/2)  
[redis_read/2](#redis_read/2)  
[redis_scan/3](#redis_scan/3)  
[redis_server/3](#redis_server/3)  
[redis_set_hash/3](#redis_set_hash/3)  
[redis_set_list/3](#redis_set_list/3)  
[redis_sscan/4](#redis_sscan/4)  
[redis_subscribe/2](#redis_subscribe/2)  
[redis_subscribe/4](#redis_subscribe/4)  
[redis_unsubscribe/2](#redis_unsubscribe/2)  
[redis_write/2](#redis_write/2)  
[redis_zscan/4](#redis_zscan/4)  
[sentinel_slave/4](#sentinel_slave/4)  
[tls_verify/5](#tls_verify/5)  
[xadd/4](#xadd/4)  
[xconsumer_stop/1](#xconsumer_stop/1)  
[xhook/2](#xhook/2)  
[xlisten/3](#xlisten/3)  
[xlisten_group/5](#xlisten_group/5)  
[xstream_set/3](#xstream_set/3)  
