
## A.38 library(portray_text): Portray text

SWI-Prolog has the special string data type. However, in Prolog, text may be represented more traditionally as a list of character-codes, i.e. (small) integers (in SWI-Prolog specifically, those are Unicode code points). This results in output like the following (here using the backquote notation which maps text to a list of codes):

``` code
?- writeln(`hello`).
[104, 101, 108, 108, 111]

?- atom_codes("hello",X).
X = [104,101,108,108,111].
```

Unless you know the Unicode tables by heart, this is pretty unpleasant for debugging. Loading `library(portray_text)` makes the toplevel and debugger consider certain lists of integers as text and print them as text. This is called "portraying". Of course, interpretation is imperfect as there is no way to tell in general whether `[65,66]` should written as `` `AB` `` or as `[65,66]`. Therefore it is important that the user be aware of the fact that this conversion is enabled. This is why this library must be loaded explicitly.

To be able to copy the printed representation and paste it back, printed text is enclosed in *back quotes* if [current_prolog_flag/2](flags.html#current_prolog_flag/2) for the flag `back_quotes` is `codes` (the default), and enclosed in *double quotes* otherwise. Certain control characters are printed out in backslash-escaped form.

The default heuristic only considers list of codes as text if the codes are all from the set of 7-bit ASCII without most of the control characters. A code is classified as text by text_code/1, which in turn calls [is_text_code/1](portraytext.html#is_text_code/1). Define portray_text:is_text_code/1 to succeed on additional codes for more flexibility (by default, that predicate succeeds nowhere). For example:

``` code
?- maplist([C,R]>>(portray_text:text_code(C)->R=y;R=n),
           `G\u00e9n\u00e9rateur`,Results).
Results = [y,n,y,n,y,y,y,y,y,y].
```

Now make [is_text_code/1](portraytext.html#is_text_code/1) accept anything:

``` code
?- [user].
|: portray_text:is_text_code(_).
|: ^D
% user://3 compiled 0.00 sec, 1 clauses
true.
```

Then:

``` code
?- maplist([C,R]>>(portray_text:text_code(C)->R=y;R=n),
           `G\u00e9n\u00e9rateur`,Results).
Results = [y,y,y,y,y,y,y,y,y,y].
```

\[det\]**portray_text**(`+OnOff:boolean`)  
Switch portraying on or off. If `true`, consider lists of integers as list of Unicode code points and print them as corresponding text inside quotes: `` `text` `` or `"text"`. Quoting depends on the value of [current_prolog_flag/2](flags.html#current_prolog_flag/2) `back_quotes`. Same as

``` code
?- set_portray_text(enabled, true).
```

\[det\]**set_portray_text**(`+Key, +Value`)  
\[det\]**set_portray_text**(`+Key, ?Old, +New`)  
Set options for portraying. Defined Keys are:

**enabled**  
Enable/disable portray text

**min_length**  
Only consider for conversion lists of integers that have a length of at least `Value`. Default is 3.

**ellipsis**  
When converting a list that is longer than `Value`, display the output as `start...end`.

\[semidet,multifile\]**is_text_code**(`+Code:nonneg`)  
Multifile hook that can be used to extend the set of character codes that is recognised as likely text. By default, [is_text_code/1](portraytext.html#is_text_code/1) fails everywhere and internally, only non-control ASCII characters (32-126) and the the control codes (9,10,13) are accepted.

To be done  
we might be able to use the current locale to include the appropriate code page. (Does that really make sense?)
