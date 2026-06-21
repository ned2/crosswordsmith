
## 3.4 Using the PceEmacs built-in editor

### 3.4.1 Activating PceEmacs

Initially [edit/1](edit.html#edit/1) uses the editor specified in the `EDITOR` environment variable. There are two ways to force it to use the built-in editor. One is to set the Prolog flag [editor](flags.html#flag:editor) to `pce_emacs` and the other is by starting the editor explicitly using the emacs/\[0,1\] predicates.

### 3.4.2 Bluffing through PceEmacs

PceEmacs closely mimics Richard Stallman's GNU-Emacs commands, adding features from modern window-based editors to make it more acceptable for beginners.^(47Decent merging with MS-Windows control-key conventions is difficult as many conflict with GNU-Emacs. Especially the cut/copy/paste commands conflict with important GNU-Emacs commands.)

At the basis, PceEmacs maps keyboard sequences to methods defined on the extended *editor* object. Some frequently used commands are, with their key-binding, presented in the menu bar above each editor window. A complete overview of the bindings for the current *mode* is provided through **Help/Show key bindings** (`Control-h Control-b`).

#### 3.4.2.1 Edit modes

Modes are the heart of (Pce)Emacs. Modes define dedicated editing support for a particular kind of (source) text. For our purpose we want *Prolog mode*. There are various ways to make PceEmacs use Prolog mode for a file.

- *Using the proper extension*  
  If the file ends in `.pl` or the selected alternative (e.g. `.pro`) extension, Prolog mode is selected.
- *Using `#!/path/to/.../swipl`*  
  If the file is a *Prolog Script* file, starting with the line `#!/path/to/swipl` `options`, Prolog mode is selected regardless of the extension.
- *Using `-*- Prolog -*-`*  
  If the above sequence appears in the first line of the file (inside a Prolog comment) Prolog mode is selected.
- *Explicit selection*  
  Finally, using **File/Mode/Prolog** you can switch to Prolog mode explicitly.

#### 3.4.2.2 Frequently used editor commands

Below we list a few important commands and how to activate them.

- *Cut/Copy/Paste*  
  These commands follow Unix/X11 traditions. You're best suited with a three-button mouse. After selecting using the left-mouse (double-click uses word-mode and triple line-mode), the selected text is *automatically* copied to the clipboard (X11 primary selection on Unix). *Cut* is achieved using the `DEL` key or by typing something else at the location. *Paste* is achieved using the middle-mouse (or wheel) button. If you don't have a middle-mouse button, pressing the left- and right-button at the same time is interpreted as a middle-button click. If nothing helps, there is the **Edit/Paste** menu entry. Text is pasted at the caret location.

- *Undo*  
  Undo is bound to the GNU-Emacs `Control-_` as well as the MS-Windows `Control-Z` sequence.

- *Abort*  
  Multi-key sequences can be aborted at any stage using `Control-G`.

- *Find*  
  Find (Search) is started using `Control-S` (forward) or `Control-R` (backward). PceEmacs implements *incremental search*. This is difficult to use for novices, but very powerful once you get the clue. After one of the above start keys, the system indicates search mode in the status line. As you are typing the search string, the system searches for it, extending the search with every character you type. It illustrates the current match using a green background.

  If the target cannot be found, PceEmacs warns you and no longer extends the search string.^(48GNU-Emacs keeps extending the string, but why? Adding more text will not make it match.) During search, some characters have special meaning. Typing anything but these characters commits the search, re-starting normal edit mode. Special commands are:

  **`Control-S`**  
  Search forwards for next.

  **`Control-R`**  
  Search backwards for next.

  **`Control-W`**  
  Extend search to next word boundary.

  **`Control-G`**  
  Cancel search, go back to where it started.

  **`ESC`**  
  Commit search, leaving caret at found location.

  **`Backspace`**  
  Remove a character from the search string.

- *Dynamic Abbreviation*  
  Also called *dabbrev*, dynamic abbreviation is an important feature of Emacs clones to support programming. After typing the first few letters of an identifier, you may press `Alt-/`, causing PceEmacs to search backwards for identifiers that start the same and use it to complete the text you typed. A second `Alt-/` searches further backwards. If there are no hits before the caret, it starts searching forwards. With some practice, this system allows for entering code very fast with nice and readable identifiers (or other difficult long words).

- *Open (a file)*  
  Is called **File/Find file** (`Control-x Control-f`). By default the file is loaded into the current window. If you want to keep this window, press `Alt-s` or click the little icon at the bottom left to make the window *sticky*.

- *Split view*  
  Sometimes you want to look at two places in the same file. To do this, use `Control-x 2` to create a new window pointing to the same file. Do not worry, you can edit as well as move around in both. `Control-x 1` kills all other windows running on the same file.

These are the most commonly used commands. In [section 3.4.3](pceemacs.html#sec:3.4.3) we discuss specific support for dealing with Prolog source code.

### 3.4.3 Prolog Mode

In the previous section ([section 3.4.2](pceemacs.html#sec:3.4.2)) we explained the basics of PceEmacs. Here we continue with Prolog-specific functionality. Possibly the most interesting is *Syntax highlighting*. Unlike most editors where this is based on simple patterns, PceEmacs syntax highlighting is achieved by Prolog itself actually reading and interpreting the source as you type it. There are three moments at which PceEmacs checks (part of) the syntax.

- *After typing a `.`*  
  After typing a `.` that is not preceded by a *symbol* character, the system assumes you completed a clause, tries to find the start of this clause and verifies the syntax. If this process succeeds it colours the elements of the clause according to the rules given below. Colouring is done using information from the last full check on this file. If it fails, the syntax error is displayed in the status line and the clause is not coloured.
- *After the command `Control-c Control-s`*  
  Acronym for **C**heck **S**yntax, it performs the same checks as above for the clause surrounding the caret. On a syntax error, however, the caret is moved to the expected location of the error.^(49In most cases the location where the parser cannot proceed is further down the file than the actual error location.)
- *After pausing for two seconds*  
  After a short pause (2 seconds), PceEmacs opens the edit buffer and reads it as a whole, creating an index of defined, called, dynamic, imported and exported predicates. After completing this, it re-reads the file and colours all clauses and calls with valid syntax.
- *After typing `Control-l Control-l`*  
  The `Control-l` command re-centers the window (scrolls the window to make the caret the center of the window). Typing this command twice starts the same process as above.

**The colour schema** itself is defined in `library(emacs/prolog_colour)`. The colouring can be extended and modified using multifile predicates. Please check this source file for details. In general, underlined objects have a popup (right-mouse button) associated with common commands such as viewing the documentation or source. **Bold** text is used to indicate the definition of objects (typically predicates when using plain Prolog). Other colours follow intuitive conventions. See [table 3](pceemacs.html#tab:plcolour).

Clauses

Blue bold

Head of an exported predicate

Red bold

Head of a predicate that is not called

Black bold

Head of remaining predicates

Calls in the clause body

Blue

Call to built-in or imported predicate

Red

Call to undefined predicate

Purple

Call to dynamic predicate

Other entities

Dark green

Comment

Dark blue

Quoted atom or string

Brown

Variable

**Table 3 :** Colour conventions

**Layout support**

Layout is not‘just nice’, it is *essential* for writing readable code. There is much debate on the proper layout of Prolog. PceEmacs, being a rather small project, supports only one particular style for layout.^(50Defined in Prolog in the file `library(emacs/prolog_mode)`, you may wish to extend this. Please contribute your extensions!) Below are examples of typical constructs.

``` code
head(arg1, arg2).

head(arg1, arg2) :- !.

head(Arg1, arg2) :- !,
        call1(Arg1).

head(Arg1, arg2) :-
        (   if(Arg1)
        ->  then
        ;   else
        ).

head(Arg1) :-
        (   a
        ;   b
        ).

head :-
        a(many,
          long,
          arguments(with,
                    many,
                    more),
          and([ a,
                long,
                list,
                with,
                a,
              | tail
              ])).
```

PceEmacs uses the same conventions as GNU-Emacs. The `TAB` key indents the current line according to the syntax rules. `Alt-q` indents all lines of the current clause. It provides support for head, calls (indented 1 tab), if-then-else, disjunction and argument lists broken across multiple lines as illustrated above.

#### 3.4.3.1 Finding your way around

The command `Alt-.` extracts name and arity from the caret location and jumps (after conformation or edit) to the definition of the predicate. It does so based on the source-location database of loaded predicates also used by [edit/1](edit.html#edit/1). This makes locating predicates reliable if all sources are loaded and up-to-date (see [make/0](consulting.html#make/0)).

In addition, references to files in [use_module/\[1,2\]](import.html#use_module/1), [consult/1](consulting.html#consult/1), etc. are red if the file cannot be found and underlined blue if the file can be loaded. A popup allows for opening the referenced file.
