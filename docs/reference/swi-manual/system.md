
## 4.35 Operating System Interaction

The predicates in this section provide basic access to the operating system that has been part of the Prolog legacy tradition. Note that more advanced access to low-level OS features is provided by several libraries from the `clib` package, notably library `library(process)`, `library(socket)`, `library(unix)` and `library(filesex)`.

**shell**(`+Command`)  
Equivalent to‘`shell(Command, 0)`’. See [shell/2](system.html#shell/2) for details.

**shell**(`+Command, -Status`)  
Execute `Command` on the operating system. `Command` is given to the Bourne shell (/bin/sh). `Status` is unified with the exit status of the command.

On Windows, [shell/\[1,2\]](system.html#shell/1) executes the command using the **CreateProcess()** API and waits for the command to terminate. If the command ends with a `&` sign, the command is handed to the **WinExec()** API, which does not wait for the new task to terminate. See also [win_exec/2](system.html#win_exec/2) and [win_shell/2](system.html#win_shell/2). Please note that the **CreateProcess()** API does **not** imply the Windows command interpreter (**cmd.exe** and therefore commands that are built in the command interpreter can only be activated using the command interpreter. For example, a file can be copied using the command below.

``` code
?- shell('cmd.exe /C copy file1.txt file2.txt').
```

Note that many of the operations that can be achieved using the shell built-in commands can easily be achieved using Prolog primitives. See [make_directory/1](files.html#make_directory/1), [delete_file/1](files.html#delete_file/1), [rename_file/2](files.html#rename_file/2), etc. The clib package provides `library(filesex)`, implementing various high level file operations such as copy_file/2. Using Prolog primitives instead of shell commands improves the portability of your program.

The library `library(process)` provides process_create/3 and several related primitives that support more fine-grained interaction with processes, including I/O redirection and management of asynchronous processes.

**getenv**(`+Name, -Value`)  
Get environment variable. Fails silently if the variable does not exist. Please note that environment variable names are case-sensitive on Unix systems and case-insensitive on Windows.

**setenv**(`+Name, +Value`)  
Set an environment variable. `Name` and `Value` must be instantiated to atoms or integers. The environment variable will be passed to shell/\[0-2\] and can be requested using [getenv/2](system.html#getenv/2). They also influence [expand_file_name/2](files.html#expand_file_name/2). Environment variables are shared between threads. Depending on the underlying C library, [setenv/2](system.html#setenv/2) and [unsetenv/1](system.html#unsetenv/1) may not be thread-safe and may cause memory leaks. Only changing the environment once and before starting threads is safe in all versions of SWI-Prolog.

**unsetenv**(`+Name`)  
Remove an environment variable from the environment. Some systems lack the underlying **unsetenv()** library function. On these systems [unsetenv/1](system.html#unsetenv/1) sets the variable to the empty string.

**setlocale**(`+Category, -Old, +New`)  
Set/Query the *locale* setting which tells the C library how to interpret text files, write numbers, dates, etc. Category is one of `all`, `collate`, `ctype`, `messages`, `monetary`, `numeric` or `time`. For details, please consult the C library locale documentation. See also [section 2.18.1](widechars.html#sec:2.18.1). Please note that the locale is shared between all threads and thread-safe usage of [setlocale/3](system.html#setlocale/3) is in general not possible. Do locale operations before starting threads or thoroughly study threading aspects of locale support in your environment before using in multithreaded environments. Locale settings are used by [format_time/3](system.html#format_time/3), [collation_key/2](chartype.html#collation_key/2) and [locale_sort/2](chartype.html#locale_sort/2).

The `messages` locale defines the language used by [print_message/2](printmsg.html#print_message/2). Note that this locale is not available on all operating system. Notably Windows does not support this category. See [win_get_user_preferred_ui_languages/2](system.html#win_get_user_preferred_ui_languages/2).

### 4.35.1 Windows-specific Operating System Interaction

The predicates in this section are only available on the Windows version of SWI-Prolog. Their use is discouraged if there are portable alternatives. For example, [win_exec/2](system.html#win_exec/2) and [win_shell/2](system.html#win_shell/2) can often be replaced by the more portable [shell/2](system.html#shell/2) or the more powerful process_create/3.

**win_exec**(`+Command, +Show`)  
Windows only. Spawns a Windows task without waiting for its completion. `Show` is one of the Win32 `SW_*` constants written in lowercase without the `SW_*`: `hide` `maximize` `minimize` `restore` `show` `showdefault` `showmaximized` `showminimized` `showminnoactive` `showna` `shownoactive` `shownormal`. In addition, `iconic` is a synonym for `minimize` and `normal` for `shownormal`.

**win_shell**(`+Operation, +File, +Show`)  
Windows only. Opens the document `File` using the Windows shell rules for doing so. `Operation` is one of `open`, `print` or `explore` or another operation registered with the shell for the given document type. On modern systems it is also possible to pass a URL as `File`, opening the URL in Windows default browser. This call interfaces to the Win32 API **ShellExecute()**. The `Show` argument determines the initial state of the opened window (if any). See [win_exec/2](system.html#win_exec/2) for defined values.

**win_shell**(`+Operation, +File`)  
Same as `win_shell(Operation, File, normal)`.

**win_registry_get_value**(`+Key, +Name, -Value`)  
Windows only. Fetches the value of a Windows registry key. `Key` is an atom formed as a path name describing the desired registry key. `Name` is the desired attribute name of the key. `Value` is unified with the value. If the value is of type `DWORD`, the value is returned as an integer. If the value is a string, it is returned as a Prolog atom. Other types are currently not supported. The default‘root’is `HKEY_CURRENT_USER`. Other roots can be specified explicitly as `HKEY_CLASSES_ROOT`, `HKEY_CURRENT_USER`, `HKEY_LOCAL_MACHINE` or `HKEY_USERS`. The example below fetches the extension to use for Prolog files (see `README.TXT` on the Windows version):

``` code
?- win_registry_get_value(
       'HKEY_LOCAL_MACHINE/Software/SWI/Prolog',
       fileExtension,
       Ext).

Ext = pl
```

**win_folder**(`?Name, -Directory`)  
True if `Name` is the Windows‘CSIDL’of `Directory`. If `Name` is unbound, all known Windows special paths are generated. `Name` is the CSIDL after deleting the leading `CSIDL_` and mapping the constant to lowercase. Check the Windows documentation for the function **SHGetSpecialFolderPath()** for a description of the defined constants. This example extracts the‘My Documents’folder:

``` code
?- win_folder(personal, MyDocuments).

MyDocuments = 'C:/Documents and Settings/jan/My Documents'
```

**win_add_dll_directory**(`+AbsDir`)  
This predicate adds a directory to the search path for dependent DLL files. If possible, this is achieved with [win_add_dll_directory/2](system.html#win_add_dll_directory/2). Otherwise, `%PATH%` is extended with the provided directory. `AbsDir` may be specified in the Prolog canonical syntax. See [prolog_to_os_filename/2](files.html#prolog_to_os_filename/2). Note that [use_foreign_library/1](foreignlink.html#use_foreign_library/1) passes an absolute path to the DLL if the destination DLL can be located from the specification using [absolute_file_name/3](files.html#absolute_file_name/3). This predicate is available from library `library(shlib)` and can be autoloaded.

**win_add_dll_directory**(`+AbsDir, -Cookie`)  
This predicate adds a directory to the search path for dependent DLL files. If the call is successful it unifies `Cookie` with a handle that must be passed to [win_remove_dll_directory/1](system.html#win_remove_dll_directory/1) to remove the directory from the search path. Error conditions:

- This predicate *fails* if Windows does not yet support the underlying primitives. These are available in recently patched Windows 7 systems and later.
- This predicate throws an exception if the provided path is invalid or the underlying Windows API returns an error.

If [open_shared_object/2](foreignlink.html#open_shared_object/2) is passed an *absolute* path to a DLL on a Windows installation that supports **AddDllDirectory()** and friends,^(148Windows 7 with up-to-date patches or Windows 8.) SWI-Prolog uses **LoadLibraryEx()** with the flags `LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR` and `LOAD_LIBRARY_SEARCH_DEFAULT_DIRS`. In this scenario, directories from `%PATH%` are *not* searched. Additional directories can be added using [win_add_dll_directory/2](system.html#win_add_dll_directory/2).

**win_remove_dll_directory**(`-Cookie`)  
Remove a DLL search directory installed using [win_add_dll_directory/2](system.html#win_add_dll_directory/2).

**win_process_modules**(`-FileNames`)  
This predicate is a wrapper around **EnumProcessModules()**. `FileNames` is unified with a list of absolute paths for all *modules* of the Windows process. Modules are the main executable file and all DLLs loaded into the process, *except* *data DLLs*. The returned file names are in canonical Prolog representation. This predicate may be used to debug loading a DLL from an unexpected location and as a helper for packaging all dependencies when creating a distribution. According to the Windows documentation this API may return incorrect results if DLLs are loaded or unloaded while **EnumProcessModules()** is in progress. See also [qsave_program/2](saved-states.html#qsave_program/2).

**win_get_user_preferred_ui_languages**(`+Format, -Languages`)  
Unifies `Languages` with a list of the user preferred languages (*Windows Display Languages*) in order of preference. If `Format` is `name`, the list elements are atoms. See [Language Names](https://docs.microsoft.com/en-us/windows/win32/intl/language-names) for details. If `Format` is `id`, `Languages` is a list of numeric language ids represented as Prolog integers. This predicate provides Windows alternative to [setlocale/3](system.html#setlocale/3) using the category `messages`.

### 4.35.2 Apple specific Operating System Interaction

Non-portable Apple MacOS specific predicates are prefixed with `apple_`.

**apple_current_locale_identifier**(`-Identifier`)  
Unify `Identifier` with the value for **CFLocaleGetIdentifier()** of the Apple current locale. The `Identifier` is an atom that consists of the primary language identifier, e.g., `en` for English followed by an underscore and an identifier for the *Region* in the MacOS *Language & Region* preferences. For example, with the primary language set to “English (UK)” and the *Region* to “United Kingdom” we get `en_GB`. This relates to the locale identifier `en_GB.UTF-8`. Unfortunately it is not that simple. For example, we can combine the primary language “English (UK)” with the *Region* “Netherlands” to end up with `en_NL` which is not a valid MacOS locale.

### 4.35.3 Dealing with time and date

Representing time in a computer system is surprisingly complicated. There are a large number of time representations in use, and the correct choice depends on factors such as compactness, resolution and desired operations. Humans tend to think about time in hours, days, months, years or centuries. Physicists think about time in seconds. But, a month does not have a defined number of seconds. Even a day does not have a defined number of seconds as sometimes a leap-second is introduced to synchronise properly with our earth's rotation. At the same time, resolution demands a range from better than pico-seconds to millions of years. Finally, civilizations have a wide range of calendars. Although there exist libraries dealing with most of this complexity, our desire to keep Prolog clean and lean stops us from fully supporting these.

For human-oriented tasks, time can be broken into years, months, days, hours, minutes, seconds and a timezone. Physicists prefer to have time in an arithmetic type representing seconds or fraction thereof, so basic arithmetic deals with comparison and durations. An additional advantage of the physicist's approach is that it requires much less space. For these reasons, SWI-Prolog uses an arithmetic type as its prime time representation.

Many C libraries deal with time using fixed-point arithmetic, dealing with a large but finite time interval at constant resolution. In our opinion, using a floating point number is a more natural choice as we can use a natural unit and the interface does not need to be changed if a higher resolution is required in the future. Our unit of choice is the second as it is the scientific unit.^(149Using Julian days is a choice made by the Eclipse team. As conversion to dates is needed for a human readable notation of time and Julian days cannot deal naturally with leap seconds, we decided for the second as our unit.) We have placed our origin at 1970-01-01T0:0:0Z for compatibility with the POSIX notion of time as well as with older time support provided by SWI-Prolog.

Where older versions of SWI-Prolog relied on the POSIX conversion functions, the current implementation uses [libtai](http://cr.yp.to/libtai.html) to realise conversion between time-stamps and calendar dates for a period of 10 million years.

#### 4.35.3.1 Time and date data structures

We use the following time representations

**TimeStamp**  
A TimeStamp is a floating point number expressing the time in seconds since the Epoch at 1970-01-01.

**date**(`Y,M,D,H,Mn,S,Off,TZ,DST`)  
We call this term a *date-time* structure. The first 5 fields are integers expressing the year, month (1..12), day (1..31), hour (0..23) and minute (0..59). The `S` field holds the seconds as a floating point number between 0.0 and 60.0. `Off` is an integer representing the offset relative to UTC in seconds, where positive values are west of Greenwich. If converted from local time (see [stamp_date_time/3](system.html#stamp_date_time/3)), `TZ` holds the name of the local timezone. If the timezone is not known, `TZ` is the atom `-`. `DST` is `true` if daylight saving time applies to the current time, `false` if daylight saving time is relevant but not effective, and `-` if unknown or the timezone has no daylight saving time.

**date**(`Y,M,D`)  
Date using the same values as described above. Extracted using [date_time_value/3](system.html#date_time_value/3).

**time**(`H,Mn,S`)  
Time using the same values as described above. Extracted using [date_time_value/3](system.html#date_time_value/3).

#### 4.35.3.2 Time and date predicates

**get_time**(`-TimeStamp`)  
Return the current time as a `TimeStamp`. The granularity is system-dependent. See [section 4.35.3.1](system.html#sec:4.35.3.1).

**stamp_date_time**(`+TimeStamp, -DateTime, +TimeZone`)  
Convert a `TimeStamp` to a `DateTime` in the given timezone. See [section 4.35.3.1](system.html#sec:4.35.3.1) for details on the data types. `TimeZone` describes the timezone for the conversion. It is one of `local` to extract the local time, `’UTC’` to extract a UTC time or an integer describing the seconds west of Greenwich.

**date_time_stamp**(`+DateTime, -TimeStamp`)  
Compute the timestamp from a **date/9** term. Values for month, day, hour, minute or second may be left unbound from the right, i.e., seconds may be unbound, or seconds and minutes, etc. It is not allowed to specify the seconds and leave one of the more significant fields unbound. Unbound values unified to their lowest possible value.

Values for month, day, hour, minute or second need not be *normalized*. This flexibility allows for easy computation of the time at any given number of these units from a given timestamp. Normalization can be achieved following this call with [stamp_date_time/3](system.html#stamp_date_time/3). This example computes the date 200 days after 2006-07-14:

``` code
?- date_time_stamp(date(2006,7,214,0,0,0,0,-,-), Stamp),
   stamp_date_time(Stamp, D, 0),
   date_time_value(date, D, Date).
Date = date(2007, 1, 30)
```

When computing a time stamp from a local time specification, the UTC offset (arg 7), TZ (arg 8) and DST (arg 9) argument may be left unbound and are unified with the proper information. The example below, executed in Amsterdam, illustrates this behaviour. On the 25th of March at 01:00, DST does not apply. At 02.00, the clock is advanced by one hour and thus both 02:00 and 03:00 represent the same time stamp.

``` code
1 ?- date_time_stamp(date(2012,3,25,1,0,0,UTCOff,TZ,DST),
                     Stamp).
UTCOff = -3600,
TZ = 'CET',
DST = false,
Stamp = 1332633600.0.

2 ?- date_time_stamp(date(2012,3,25,2,0,0,UTCOff,TZ,DST),
                     Stamp).
UTCOff = -7200,
TZ = 'CEST',
DST = true,
Stamp = 1332637200.0.

3 ?- date_time_stamp(date(2012,3,25,3,0,0,UTCOff,TZ,DST),
                     Stamp).
UTCOff = -7200,
TZ = 'CEST',
DST = true,
Stamp = 1332637200.0.
```

Note that DST and offset calculation are based on the POSIX function **mktime()**. If **mktime()** returns an error, a representation_error `dst` is generated.

**date_time_value**(`?Key, +DateTime, ?Value`)  
Extract values from a date/9 term. Provided keys are:

|                   |                                             |
|-------------------|---------------------------------------------|
| **key**           | **value**                                   |
| `year`            | Calendar year as an integer                 |
| `month`           | Calendar month as an integer 1..12          |
| `day`             | Calendar day as an integer 1..31            |
| `hour`            | Clock hour as an integer 0..23              |
| `minute`          | Clock minute as an integer 0..59            |
| `second`          | Clock second as a float 0.0..60.0           |
| `utc_offset`      | Offset to UTC in seconds (positive is west) |
| `time_zone`       | Name of timezone; fails if unknown          |
| `daylight_saving` | Bool `(`true) if dst is in effect           |
| `date`            | Term `date(Y,M,D)`                          |
| `time`            | Term `time(H,M,S)`                          |

**format_time**(`+Out, +Format, +StampOrDateTime`)  
Modelled after POSIX **strftime()**, using GNU extensions. `Out` is a destination as specified with [with_output_to/2](IO.html#with_output_to/2). `Format` is an atom or string with the following conversions. Conversions start with a percent (%) character.^(150Descriptions taken from Linux Programmer's Manual) `StampOrDateTime` is either a numeric time-stamp, a term `date(Y,M,D,H,M,S,O,TZ,DST)` or a term `date(Y,M,D)`.

- `a`  
  The abbreviated weekday name according to the current locale. Use [format_time/4](system.html#format_time/4) for POSIX locale.
- `A`  
  The full weekday name according to the current locale. Use [format_time/4](system.html#format_time/4) for POSIX locale.
- `b`  
  The abbreviated month name according to the current locale. Use [format_time/4](system.html#format_time/4) for POSIX locale.
- `B`  
  The full month name according to the current locale. Use [format_time/4](system.html#format_time/4) for POSIX locale.
- `c`  
  The preferred date and time representation for the current locale.
- `C`  
  The century number (year/100) as a 2-digit integer.
- `d`  
  The day of the month as a decimal number (range 01 to 31).
- `D`  
  Equivalent to %m/%d/%y. (For Americans only. Americans should note that in other countries %d/%m/%y is rather common. This means that in an international context this format is ambiguous and should not be used.)
- `e`  
  Like %d, the day of the month as a decimal number, but a leading zero is replaced by a space.
- `E`  
  Modifier. Not implemented.
- `f`  
  Number of microseconds. The `f` can be prefixed by an integer to print the desired number of digits. E.g., `%3f` prints milliseconds. This format is not covered by any standard, but available with different format specifiers in various incarnations of the **strftime()** function.
- `F`  
  Equivalent to %Y-%m-%d (the ISO 8601 date format).
- `g`  
  Like %G, but without century, i.e., with a 2-digit year (00-99).
- `G`  
  The ISO 8601 year with century as a decimal number. The 4-digit year corresponding to the ISO week number (see %V). This has the same format and value as %y, except that if the ISO week number belongs to the previous or next year, that year is used instead.
- `V`  
  The ISO 8601:1988 week number of the current year as a decimal number, range 01 to 53, where week 1 is the first week that has at least 4 days in the current year, and with Monday as the first day of the week. See also %U and %W.
- `h`  
  Equivalent to %b.
- `H`  
  The hour as a decimal number using a 24-hour clock (range 00 to 23).
- `I`  
  The hour as a decimal number using a 12-hour clock (range 01 to 12).
- `j`  
  The day of the year as a decimal number (range 001 to 366).
- `k`  
  The hour (24-hour clock) as a decimal number (range 0 to 23); single digits are preceded by a blank. (See also %H.)
- `l`  
  The hour (12-hour clock) as a decimal number (range 1 to 12); single digits are preceded by a blank. (See also %I.)
- `m`  
  The month as a decimal number (range 01 to 12).
- `M`  
  The minute as a decimal number (range 00 to 59).
- `n`  
  A newline character.
- `O`  
  Modifier to select locale-specific output. Not implemented.
- `p`  
  Either‘AM’or‘PM’according to the given time value, or the corresponding strings for the current locale. Noon is treated as‘pm’and midnight as‘am’.^(151Despite the above claim, some locales yield `am` or `pm` in lower case.)
- `P`  
  Like %p but in lowercase:‘am’or‘pm’or a corresponding string for the current locale.
- `r`  
  The time in a.m. or p.m. notation. In the POSIX locale this is equivalent to‘%I:%M:%S %p’.
- `R`  
  The time in 24-hour notation (%H:%M). For a version including the seconds, see %T below.
- `s`  
  The number of seconds since the Epoch, i.e., since 1970-01-01 00:00:00 UTC.
- `S`  
  The second as a decimal number (range 00 to 60). (The range is up to 60 to allow for occasional leap seconds.)
- `t`  
  A tab character.
- `T`  
  The time in 24-hour notation (%H:%M:%S).
- `u`  
  The day of the week as a decimal, range 1 to 7, Monday being 1. See also %w.
- `U`  
  The week number of the current year as a decimal number, range 00 to 53, starting with the first Sunday as the first day of week 01. See also %V and %W.
- `w`  
  The day of the week as a decimal, range 0 to 6, Sunday being 0. See also %u.
- `W`  
  The week number of the current year as a decimal number, range 00 to 53, starting with the first Monday as the first day of week 01.
- `x`  
  The preferred date representation for the current locale without the time.
- `X`  
  The preferred time representation for the current locale without the date.
- `y`  
  The year as a decimal number without a century (range 00 to 99).
- `Y`  
  The year as a decimal number including the century.
- `z`  
  The timezone as hour offset from GMT using the format HHmm. Required to emit RFC822-conforming dates (using `’%a, %d %b %Y %T %z’`). Our implementation supports `%:z`, which modifies the output to HH:mm as required by XML-Schema. Note that both notations are valid in ISO 8601. The sequence `%:z` is compatible to the GNU date(1) command.
- `Z`  
  The timezone or name or abbreviation.
- `+`  
  The date and time in date(1) format.
- `%`  
  A literal‘%’character.

The table below gives some format strings for popular time representations. RFC1123 is used by HTTP. The full implementation of http_timestamp/2 as available from `library(http/http_header)` is here.

``` code
http_timestamp(Time, Atom) :-
        stamp_date_time(Time, Date, 'UTC'),
        format_time(atom(Atom),
                    '%a, %d %b %Y %T GMT',
                    Date, posix).
```

|              |                         |
|--------------|-------------------------|
| **Standard** | **Format string**       |
| **xsd**      | `’%FT%T%:z’`            |
| **ISO8601**  | `’%FT%T%z’`             |
| **RFC822**   | `’%a, %d %b %Y %T %z’`  |
| **RFC1123**  | `’%a, %d %b %Y %T GMT’` |

**format_time**(`+Out, +Format, +StampOrDateTime, +Locale`)  
Format time given a specified `Locale`. This predicate is a work-around for lacking proper portable and thread-safe time and locale handling in current C libraries. In its current implementation the only value allowed for `Locale` is `posix`, which currently only modifies the behaviour of the `a`, `A`, `b` and `B` format specifiers. The predicate is used to be able to emit POSIX locale week and month names for emitting standardised time-stamps such as RFC1123.

**parse_time**(`+Text, -Stamp`)  
Same as `parse_time(Text, _Format, Stamp)`. See [parse_time/3](system.html#parse_time/3).

**parse_time**(`+Text, ?Format, -Stamp`)  
Parse a textual time representation, producing a time-stamp. Supported formats for `Text` are in the table below. If the format is known, it may be given to reduce parse time and avoid ambiguities. Otherwise, `Format` is unified with the format encountered.

|          |                                    |
|----------|------------------------------------|
| **Name** | **Example**                        |
| rfc_1123 | `Fri, 08 Dec 2006 15:29:44 GMT `   |
|          | `Fri, 08 Dec 2006 15:29:44 +0000 ` |
| iso_8601 | `2006-12-08T17:29:44+02:00 `       |
|          | `20061208T172944+0200 `            |
|          | `2006-12-08T15:29Z `               |
|          | `2006-12-08 `                      |
|          | `20061208 `                        |
|          | `2006-12 `                         |
|          | `2006-W49-5 `                      |
|          | `2006-342 `                        |

**day_of_the_week**(`+Date,-DayOfTheWeek`)  
Computes the day of the week for a given date. `Date`` = date(``Year``,``Month``,``Day``)`. Days of the week are numbered from one to seven: Monday = 1, Tuesday = 2, ... , Sunday = 7.

### 4.35.4 Controlling the **swipl-win** (Epilog) console window

The SWI-Prolog package `xpce` provides the SWI-Prolog GUI capabilities. Part of it is a console window that provides a menu, called *Epilog*. If Prolog is started using **swipl-win**, it starts as a GUI application running the Prolog REPL loop in an Epilog window. SWI-Prolog allows opening multiple such console windows, each running the REPL loop in a new thread. New Epilog windows can be started from several menus in the GUI as well as using the predicate [epilog/0](system.html#epilog/0).

**epilog**  
pen a new Epilog console. SWI-Prolog creates a new thread that runs the REPL loop and has its input and output connected to the new console.

**epilog**(`+Options`)  
Open an Epilog console with additional `Options`.

**title**(`+String`)  
Set the window title.

**rows**(`+Integer`)  
Height of the initial terminal in lines (default 25).

**cols**(`+Integer`)  
Width of the initial terminal in characters (default 80).

**init**(`:Goal`)  
Run Goal as initialization goal. Default calls [version/0](printmsg.html#version/0) for the first and [true/0](control.html#true/0) for subsequent terminals.

**goal**(`:Goal`)  
Run Goal as REPL loop. Default runs [prolog/0](toplevel.html#prolog/0).

**main**(`+Bool`)  
If `true`, act as application main window. In this case [epilog/1](system.html#epilog/1) runs the main thread and returns after all windows have been closed.

**set_epilog**(`+Option`)  
Modify the Epilog console attached to the calling thread. Option is one of:

**title**(`+Title`)  
**foreground**(`+Color`)  
**background**(`+Color`)  
**selection_foreground**(`+Color`)  
**selection_background**(`+Color`)  
**menu**(`+Label, +Before`)  
Add a new popup to the Epilog menu. The popus is added before Before. If Before is‘-\`, the new popup is added to the right.

**menu_item**(`+PopupName, +Item, +Before, :Goal`)  
Insert an item in the Epilog console menu. PopupName is the popup in which to insert the item. Item is the name for the new item. If Item is `--`, a *separator* is inserted. Before is the name of the item before which to insert the new item. If this is `-`, the item is appended.

Goal is *injected* into the current terminal of the Epilog window. This implies that we assume that the console is waiting for the user. Eventually, we probably want a more flexible solution.

This predicate raises `existence_error(epilog, Thread)` if the (calling) thread has no attached Epilog window. Threads are started in an Epilog window have the Prolog flag [console_menu](flags.html#flag:console_menu) set to `true`.
