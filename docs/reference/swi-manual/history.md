
## 2.8 Command line history

SWI-Prolog offers a query substitution mechanism similar to what is seen in Unix shells. The availability of this feature is controlled by [set_prolog_flag/2](flags.html#set_prolog_flag/2), using the **history** Prolog flag. By default, history is available if no interactive command line editor is available. To enable history, remembering the last 50 commands, put the following into your startup file (see [section 2.2](initfile.html#sec:2.2)):

``` code
:- set_prolog_flag(history, 50).
```

The history system allows the user to compose new queries from those typed before and remembered by the system. The available history commands are shown in [table 1](history.html#tab:history). History expansion is not done if these sequences appear in quoted atoms or strings.

|         |                                           |
|---------|-------------------------------------------|
| `!!.`   | Repeat last query                         |
| `!nr.`  | Repeat query numbered \<`nr`\>            |
| `!str.` | Repeat last query starting with \<`str`\> |
| `h.`    | Show history of commands                  |
| `!h.`   | Show this list                            |

**Table 1 :** History commands
