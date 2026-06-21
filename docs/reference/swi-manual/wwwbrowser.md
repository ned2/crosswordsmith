
## A.30 library(www_browser): Open a URL in the users browser

This library deals with the highly platform specific task of opening a web page. In addition, is provides a mechanism similar to [absolute_file_name/3](files.html#absolute_file_name/3) that expands compound terms to concrete URLs. For example, the SWI-Prolog home page can be opened using:

``` code
?- www_open_url(swipl(.)).
```

**www_open_url**(`+Url`)  
Open URL in running version of the users’browser or start a new browser. This predicate tries the following steps:

1.  If a prolog flag (see [set_prolog_flag/2](flags.html#set_prolog_flag/2)) `browser` is set and this is the name of a known executable, use this. The flag may be set to `Command-Mode`, where mode is one of `fg` or `bg`, requesting Command to run in foreground or background mode. Default is `bg`.
2.  On Windows, use `win_shell(open, URL)`
3.  Find a generic‘open’comment. Candidates are `xdg-open`, `open` or `gnome-open`.
4.  If a environment variable `BROWSER` is set and this is the name of a known executable, use this.
5.  Try to find a known browser. @tbd Figure out the right tool in step 3 as it is not uncommon that multiple are installed.

\[multifile\]**known_browser**(`+FileBaseName, -Compatible`)  
True if browser `FileBaseName` has a remote protocol compatible to `Compatible`.

**expand_url_path**(`+Spec, -URL`)  
Expand `URL` specifications similar to [absolute_file_name/3](files.html#absolute_file_name/3). The predicate url_path/2 plays the role of [file_search_path/2](consulting.html#file_search_path/2).

Errors  
`existence_error(url_path, Spec)` if the location is not defined.
