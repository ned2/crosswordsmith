
## A.2 library(ansi_term): Print decorated text to ANSI consoles

See also  
[http://en.wikipedia.org/wiki/ANSI_escape_code](http://en.wikipedia.org/wiki/ANSI_escape_code)

This library allows for exploiting the color and attribute facilities of most modern terminals using ANSI escape sequences. This library provides the following:

- [ansi_format/3](ansiterm.html#ansi_format/3) allows writing messages to the terminal with ansi attributes.
- It defines the hook [prolog:message_line_element/2](printmsg.html#prolog:message_line_element/2), which provides ansi attributes and hyperlinks for [print_message/2](printmsg.html#print_message/2).

The behavior of this library is controlled by two Prolog flags:

**`[|]`**(`99, [111,108,111,114,95,116,101,114,109]`)  
When `true`, activate the color output for this library. Otherwise simply call [format/3](format.html#format/3).

**`[|]`**(`104, [121,112,101,114,108,105,110,107,95,116,101,114,109]`)  
Emit terminal hyperlinks for `url(Location)` and `url(URL, Label)` elements of Prolog messages.

&nbsp;

\[det\]**ansi_format**(`+ClassOrAttributes, +Format, +Args`)  
\[det\]**ansi_format**(`+Stream, +ClassOrAttributes, +Format, +Args`)  
`Format` text with ANSI attributes. This predicate behaves as [format/2](format.html#format/2) using `Format` and `Args`, but if the `current_output` is a terminal, it adds ANSI escape sequences according to Attributes. For example, to print a text in bold cyan, do

``` code
?- ansi_format([bold,fg(cyan)], 'Hello ~w', [world]).
```

Attributes is either a single attribute, a list thereof or a term that is mapped to concrete attributes based on the current theme (see [prolog:console_color/2](ansiterm.html#prolog:console_color/2)). The attribute names are derived from the ANSI specification. See the source for sgr_code/2 for details. Some commonly used attributes are:

**bold**  
**underline**  
`fg(Color)` **`,`** `bg(Color)``,``hfg(Color)``,``hbg(Color)`  
For `fg(Color)` and `bg(Color)`, the colour name can be’#RGB’or’#RRGGBB’

`fg8(Spec)` **`,`** `bg8(Spec)`  
8-bit color specification. `Spec` is a colour name, `h(Color)` or an integer 0..255.

`fg(R, G, B)` **`,`** `bg(R, G, B)`  
24-bit (direct color) specification. The components are integers in the range 0..255.

Defined color constants are below. `default` can be used to access the default color of the terminal.

- black, red, green, yellow, blue, magenta, cyan, white

ANSI sequences are sent if and only if

- The `current_output` has the property `tty(true)` (see [stream_property/2](IO.html#stream_property/2)).
- The Prolog flag [color_term](flags.html#flag:color_term) is `true`.

\[semidet,multifile\]prolog:**console_color**(`+Term, -AnsiAttributes`)  
Hook that allows for mapping abstract terms to concrete ANSI attributes. This hook is used by *theme* files to adjust the rendering based on user preferences and context. Defaults are defined in the file `boot/messages.pl`, default_theme/2.

See also  
`library(theme/dark)` for an example implementation and the `Term` values used by the system messages.

\[semidet,multifile\]prolog:**message_line_element**(`+Stream, +Term`)  
Hook implementation that deals with `ansi(+Attr, +Fmt, +Args)` in message specifications.

\[det\]**ansi_hyperlink**(`+Stream, +Location`)  
\[det\]**ansi_hyperlink**(`+Stream, +URL, +Label`)  
Create a hyperlink for a terminal emulator. The file is fairly easy, but getting the line and column across is not as there seems to be no established standard. The current implementation emits, i.e., inserting a capital `L` before the line.

``` code
``file://AbsFileName[#LLine[:Column]]``
```

See also  
[https://gist.github.com/egmontkob/eb114294efbcd5adb1944c9f3cb5feda](https://gist.github.com/egmontkob/eb114294efbcd5adb1944c9f3cb5feda)

\[multifile\]**tty_url_hook**(`+Location, -URL`)  
Hook for location_url/2.

\[semidet\]**ansi_get_color**(`+Which, -RGB`)  
Obtain the `RGB` color for an ANSI color parameter. `Which` is either a color alias or an integer ANSI color id. Defined aliases are `foreground` and `background`. This predicate sends a request to the console (`user_output`) and reads the reply. This assumes an `xterm` compatible terminal.

|  |  |
|----|----|
| `RGB` | is a term `rgb(Red,Green,Blue)`. The color components are integers in the range 0..65535. |
