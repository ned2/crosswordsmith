
## 2.15 The SWI-Prolog syntax

SWI-Prolog syntax is close to ISO-Prolog standard syntax, which is based on the Edinburgh Prolog syntax. A formal description can be found in the ISO standard document. For an informal introduction we refer to Prolog text books (see [section 1](intro.html#sec:1)) and [online tutorials](http://www.swi-prolog.org/Links.html). In addition to the differences from the ISO standard documented here, SWI-Prolog offers several extensions, some of which also extend the syntax. See [section 5](extensions.html#sec:5) for more information.

### 2.15.1 ISO Syntax Support

This section lists various extensions w.r.t. the ISO Prolog syntax.

#### 2.15.1.1 Processor Character Set

The processor character set specifies the class of each character used for parsing Prolog source text. Character classification is fixed to [Unicode](http://www.unicode.org/). See also [section 2.18](widechars.html#sec:2.18).

#### 2.15.1.2 Nested comments

SWI-Prolog allows for nesting `/* ... */` comments. Where the ISO standard accepts `/* ... /* ... */` as a comment, SWI-Prolog will search for a terminating `*/`. This is useful if some code with `/* ... */` comment statements in it should be commented out. This modification also avoids unintended commenting in the example below, where the closing `*/` of the first comment has been forgotten.^(28Recent copies of GCC give a style warning if `/*` is encountered in a comment, which suggests that this problem has been recognised more widely.)

``` code
/* comment

code

/* second comment */

code
```

#### 2.15.1.3 Character Escape Syntax

Within quoted atoms (using single quotes: `’<``atom``>’`) special characters are represented using escape sequences. An escape sequence is led in by the backslash (`\`) character. The list of escape sequences is compatible with the ISO standard but contains some extensions, and the interpretation of numerically specified characters is slightly more flexible to improve compatibility. Undefined escape characters raise a `syntax_error` exception.^(29Up to SWI-Prolog 6.1.9, undefined escape characters were copied verbatim, i.e., removing the backslash.)

`\a`  
Alert character. Normally the ASCII character 7 (beep).

`\b`  
Backspace character.

`\c`  
No output. All input characters up to but not including the first non-layout character are skipped. This allows for the specification of pretty-looking long lines. Not supported by ISO. Example:

``` code
format('This is a long line that looks better if it was \c
       split across multiple physical lines in the input')
```

`\<``NEWLINE``>`  
When in ISO mode (see the Prolog flag [iso](flags.html#flag:iso)), only skip this sequence. In native mode, white space that follows the newline is skipped as well and a warning is printed, indicating that this construct is deprecated and advising to use `\c`. We advise using `\c` or putting the layout *before* the `\`, as shown below. Using `\c` is supported by various other Prolog implementations and will remain supported by SWI-Prolog. The style shown below is the most compatible solution.^(30Future versions will interpret `\`\<`return`\> according to ISO.)

``` code
format('This is a long line that looks better if it was \
split across multiple physical lines in the input')
```

instead of

``` code
format('This is a long line that looks better if it was\
 split across multiple physical lines in the input')
```

Note that SWI-Prolog also allows unescaped newlines to appear in quoted material. This is not allowed by the ISO standard, but used to be common practice before.

`\e`  
Escape character (ASCII 27). Not ISO, but widely supported.

`\f`  
Form-feed character.

`\n`  
Next-line character.

`\r`  
Carriage-return only (i.e., go back to the start of the line).

`\s`  
Space character. Intended to allow writing `0'\s` to get the character code of the space character. Not ISO.

`\t`  
Horizontal tab character.

`\v`  
Vertical tab character (ASCII 11).

`\``xXX..\`  
Hexadecimal specification of a character. The closing `\` is obligatory according to the ISO standard, but optional in SWI-Prolog to enhance compatibility with the older Edinburgh standard. The code `\xa\3` emits the character 10 (hexadecimal‘a’) followed by‘3’. Characters specified this way are interpreted as Unicode characters. See also `\u`.

`\uXXXX`  
Unicode character specification where the character is specified using *exactly* 4 hexadecimal digits. This is an extension to the ISO standard, fixing two problems. First, where `\x` defines a numeric character code, it doesn't specify the character set in which the character should be interpreted. Second, it is not needed to use the idiosyncratic closing `\` ISO Prolog syntax.

`\UXXXXXXXX`  
Same as `\uXXXX`, but using 8 digits to cover the whole Unicode set.

`\40`  
Octal character specification. The rules and remarks for hexadecimal specifications apply to octal specifications as well.

`\``\`  
Escapes the backslash itself. Thus, `'\\'` is an atom consisting of a single `\`.

`\’`  
Single quote. Note that `'\''` and `''''` both describe the atom with a single `’`, i.e., `'\'' == ''''` is true.

`\"`  
Double quote.

`\‘`  
Back quote.

Character escaping is only available if `current_prolog_flag(character_escapes, true)` is active (default). See [current_prolog_flag/2](flags.html#current_prolog_flag/2). Character escapes conflict with [writef/2](writef.html#writef/2) in two ways: `\40` is interpreted as decimal 40 by [writef/2](writef.html#writef/2), but as octal 40 (decimal 32) by `read`. Also, the [writef/2](writef.html#writef/2) sequence `\l` is illegal. It is advised to use the more widely supported [format/\[2,3\]](format.html#format/2) predicate instead. If you insist upon using [writef/2](writef.html#writef/2), either switch [character_escapes](flags.html#flag:character_escapes) to `false`, or use double `\\`, as in `writef('\\l')`.

#### 2.15.1.4 Syntax for non-decimal numbers

SWI-Prolog implements both Edinburgh and ISO representations for non-decimal numbers. According to Edinburgh syntax, such numbers are written as `<``radix``>’<``number``>`, where \<`radix`\> is a number between 2 and 36. ISO defines binary, octal and hexadecimal numbers using `0`*`[bxo]`*`<``number``>`. For example: `A is 0b100 \/ 0xf00` is a valid expression. Such numbers are always unsigned.

#### 2.15.1.5 Using digit groups in large integers

SWI-Prolog supports splitting long integers into *digit groups*. Digit groups can be separated with the sequence \<`underscore`\>, \<`optional white space`\>. If the \<`radix`\> is 10 or lower, they may also be separated with exactly one space. The following all express the integer 1 million:

``` code
1_000_000
1 000 000
1_000_/*more*/000
```

Integers can be printed using this notation with [format/2](format.html#format/2), using the `~I` format specifier. For example:

``` code
?- format('~I', [1000000]).
1_000_000
```

The current syntax has been proposed by Ulrich Neumerkel on the SWI-Prolog mailinglist.

#### 2.15.1.6 Rational number syntax

As of version 8.1.22, SWI-Prolog supports rational numbers as a primary citizen atomic data type if SWI-Prolog is compiled with the GMP library. This can be tested using the [bounded](flags.html#flag:bounded) Prolog flag. An atomic type also requires a syntax. Unfortunately there are few options for adding rational numbers without breaking the ISO standard.^(31ECLiPSe uses `numerator`\_`denominator`. This syntax conflicts with SWI-Prolog digit groups (see [section 2.15.1.5](syntax.html#sec:2.15.1.5)) and does not have a recognised link to rational numbers. The notation `1/3r` and `1/3R` have also been proposed. The `1/3r` is compatible to Ruby, but is hard to parse due to the required look-ahead and not very natural. See also [https://en.wikipedia.org/wiki/Rational_data_type](https://en.wikipedia.org/wiki/Rational_data_type).)

ECLiPSe and SWI-Prolog have agreed to define the canonical syntax for rational numbers to be e.g., `1r3`. In addition, ECLiPSe accepts `1_3` and SWI-Prolog can be asked to accept `1/3` using the module sensitive Prolog flag [rational_syntax](flags.html#flag:rational_syntax), which has the values below. Note that [write_canonical/1](termrw.html#write_canonical/1) always uses the compatible `1r3` syntax.

**natural**  
This is the default mode where we ignore the ambiguity issue and follow the most natural \<`integer`\>/\<`nonneg`\> alternative. Here, \<`integer`\> follows the normal rules for Prolog decimal integers and \<`nonneg`\> does the same, but does not allows for a sign. Note that the parser translates a rational number to its canonical form which implies there are no common divisors in the resulting numerator and denominator. Examples of ration numbers are:

|                  |         |
|------------------|---------|
| 1/2              | 1/2     |
| 2/4              | 1/2     |
| 1 000 000/33 000 | 1000/33 |
| -3/5             | -3/5    |

We expect very few programs to have text parsed into a rational number while a term was expected. Note that for rationals appearing in an arithmetic expression the only difference is that evaluation moves from runtime to compiletime. The utility [list_rationals/0](check.html#list_rationals/0) may be used on a loaded program to check whether the program contains rational numbers inside clauses and thus may be subject to compatibility issues. If a term is intended this can be written as `/(1,2)`, `(1)/2`, `1 / 2` or some variation thereof.

**compatibility**  
Read and write rational numbers as e.g., `1r3`. In other words, this adheres to the same rules as `natural` above, but using the‘`r`’instead of‘`/`’. Note that this may conflict with traditional Prolog as‘`r`’can be defined as an infix operator. The same argument holds for `0x23` and similar syntax for numbers that are part of the ISO standard.

While the syntax is controlled by the flag [rational_syntax](flags.html#flag:rational_syntax), behavior on integer division and exponentiation is controlled by the flag [prefer_rationals](flags.html#flag:prefer_rationals). See section [section 4.27.2.2](arith.html#sec:4.27.2.2) for arithmetic on rational numbers.

#### 2.15.1.7 NaN and Infinity floats and their syntax

SWI-Prolog supports reading and printing‘special’floating point values according to [Proposal for Prolog Standard core update wrt floating point arithmetic](http://eclipseclp.org/Specs/core_update_float.html) by Joachim Schimpf and available in ECLiPSe Prolog. In particular,

- Infinity is printed as `1.0Inf` or `-1.0Inf`. Any sequence matching the regular expression `[+-]?\sd+[.]\sd+Inf` is mapped to plus or minus infinity.

- `NaN` (Not a Number) is printed as `1.xxxNaN`, where *1.xxx* is the float after replacing the exponent by‘1’. Such numbers are read, resulting in the same `NaN`. The `NaN` constant can also be produced using the function [nan/0](arith.html#f-nan/0), e.g.,

  ``` code
  ?- A is nan.
  A = 1.5NaN.
  ```

By default SWI-Prolog arithmetic (see [section 4.27](arith.html#sec:4.27)) follows the ISO standard with describes that floating point operations either produce a *normal* floating point number or raise an exception. [section 4.27.2.4](arith.html#sec:4.27.2.4) describes the Prolog flags that can be used to support the IEEE special float values. The ability to create, read and write such values facilitates the exchange of data with languages that can represent the full range of IEEE doubles.

#### 2.15.1.8 Force only underscore to introduce a variable

According to the ISO standard and most Prolog systems, identifiers that start with an uppercase letter or an underscore are variables. In the past, *Prolog by BIM* provided an alternative syntax, where only the underscore (`_`) introduces a variable. As of SWI-Prolog 7.3.27 SWI-Prolog supports this alternative syntax, controlled by the Prolog flag [var_prefix](flags.html#flag:var_prefix). As the [character_escapes](flags.html#flag:character_escapes) flag, this flag is maintained per module, where the default is `false`, supporting standard syntax.

Having only the underscore introduce a variable is particularly useful if code contains identifiers for case sensitive external languages. Examples are the RDF library where code frequently specifies property and class names^(32Samer Abdallah suggested this feature based on experience with non-Prolog users using the RDF library.) and the R interface for specifying functions or variables that start with an uppercase character. Lexical databases where part of the terms start with an uppercase letter is another category were the readability of the code improves using this option.

#### 2.15.1.9 Unicode Prolog source

The ISO standard specifies the Prolog syntax in ASCII characters. As SWI-Prolog supports Unicode in source files we must extend the syntax. This section describes the implication for the source files, while writing international source files is described in [section 3.1.3](projectfiles.html#sec:3.1.3).

The SWI-Prolog Unicode character classification follows the Unicode release reported by the read-only Prolog flag [unicode_syntax_version](flags.html#flag:unicode_syntax_version). Note that [char_type/2](chartype.html#char_type/2) and friends, intended for processing arbitrary text rather than Prolog source code, are based on the C library locale-based classification routines, and that the predicates in `library(unicode)` report the version of the bundled *utf8proc* data, which may differ from the syntax classifier's version (see unicode_version/1).

- *Quoted atoms and strings*  
  Any character of any script can be used in quoted atoms and strings. The escape sequences `\uXXXX` and `\UXXXXXXXX` (see [section 2.15.1.3](syntax.html#sec:2.15.1.3)) were introduced to specify Unicode code points in ASCII files.

- *Atoms and Variables*  
  We handle them in one item as they are closely related. The Unicode standard defines a syntax for identifiers in computer languages.^(33[http://www.unicode.org/reports/tr31/](http://www.unicode.org/reports/tr31/)) SWI-Prolog uses the `XID_Start` and `XID_Continue` sets: identifiers start with an `XID_Start` code point followed by a sequence of `XID_Continue` code points. As a profile addition, the superscript digits (², ³, ¹, and ⁰--⁹, i.e. U+00B2, U+00B3, U+00B9, U+2070, U+2074..U+2079) and the subscript digits (₀--₉, U+2080..U+2089) are also accepted as `XID_Continue`, allowing variables such as X² and X₁. Such sequences are handled as a single token. The token is a *variable* iff it starts with an underscore (`_`) or with a code point in general category `Lu` (uppercase letter). Otherwise it is an atom. Note that titlecase letters (general category `Lt`, e.g. ǅ) start an atom, not a variable; this differs from earlier releases that used the broader derived `Uppercase` property. Many languages do not have the notion of character case; in such languages variables *must* be written as `_name`.

- *Numbers*  

  In source code ([read_term/2](termrw.html#read_term/2)), numeric literals use ASCII digits `0`--`9` only. Conversion via [atom_number/2](manipatom.html#atom_number/2), [number_codes/2](manipatom.html#number_codes/2), and [number_chars/2](manipatom.html#number_chars/2) additionally accepts any Unicode `Nd` block for integers, rational numbers (see [section 2.15.1.6](syntax.html#sec:2.15.1.6)) and floating point numbers; in a single number all digits must come from the same block, i.e., if the numerator of a rational uses Indian script the denominator must too. The sign, rational separator, floating point `.`, and floating point exponent are always ASCII.

- *White space*  
  The layout characters are exactly the Unicode `Pattern_White_Space` set defined by UAX #31: U+0009..U+000D, U+0020, U+0085, U+200E, U+200F, U+2028, and U+2029. NBSP (U+00A0) is deliberately excluded from `Pattern_White_Space`; appearing outside quoted material it raises a stray-character syntax error. Programs that paste from word processors will occasionally encounter NBSP in the wrong place, and reporting it explicitly is preferable to silently treating it as a separator.

- *Line termination*  
  Seven of the Pattern_White_Space code points end a line: U+000A (LF), U+000B (VT), U+000C (FF), U+000D (CR), U+0085 (NEL), U+2028 (LINE SEPARATOR), and U+2029 (PARAGRAPH SEPARATOR). They terminate `%` line comments, drive the source-position line counter, and act as the *newline* for backslash-newline continuation inside quoted strings (`\<EOL>` followed by zero or more blanks is consumed). The same set is exposed to user code through `prolog_end_of_line` in [code_type/2](chartype.html#code_type/2) and [char_type/2](chartype.html#char_type/2); the unprefixed `end_of_line` stays restricted to the four ISO/POSIX control codes (LF, VT, FF, CR), and the eleven-member Pattern_White_Space set itself is `prolog_layout`.

- *Stray characters in source text*  
  At token-start position (where layout is allowed), a code point that is not in any of the recognised syntax classes --- layout (see above), decimal digit, identifier-start, identifier-continue, solo, bracket open or quote open --- raises `syntax_error(illegal_character)`. This includes the C0 and C1 control range, unassigned and noncharacter code points, surrogates, the `Zs` / `Zl` / `Zp` separator classes that are not in `Pattern_White_Space` (NBSP, OGHAM SPACE MARK, NARROW NO-BREAK SPACE, IDEOGRAPHIC SPACE, ...), `Cf` format characters that are not in `Pattern_White_Space` and not in `Other_ID_Continue` (SOFT HYPHEN, ZERO WIDTH SPACE, ...), enclosing combining marks (`Me`), and other-number characters (`No`, e.g. the vulgar fractions and Roman-numeral form U+00BC..U+00BE) that are not in the explicit super- or subscript-digit profile.

  Non-spacing combining marks (`Mn`, `Mc`) are likewise rejected at token-start position --- they do not start an identifier --- but *are* in `XID_Continue` so they absorb into a preceding identifier (the sequence U+0061 followed by U+0300 COMBINING GRAVE reads as a single-token identifier two code points long).

- *Inside quoted material*  
  Inside single-quoted atoms (`'...'`), double-quoted strings (`"..."`), back-quoted text (`` `...` ``), Unicode quote pairs (see above), `%` comments, and `/* ... */` comments, *any* Unicode scalar value (U+0000 to U+10FFFF, excluding surrogates which UTF-8 cannot encode) is accepted verbatim. The escape sequences `\uXXXX` and `\UXXXXXXXX` ([section 2.15.1.3](syntax.html#sec:2.15.1.3)) are available for portability and explicit clarity, not as a gate. The single exception is the bidirectional override / isolate range (U+202A..U+202E and U+2066..U+2069), which is rejected as a Trojan-source defense (see [unicode_atoms](flags.html#flag:unicode_atoms)). Quoted writing (e.g. [writeq/1](termrw.html#writeq/1)) of an atom or string that contains a control or zero-width code point causes the atom or string to be quoted and the offending code points to be written using an escape sequence; see [section 2.15.1.3](syntax.html#sec:2.15.1.3).

- *Other characters*  
  The first 128 characters follow the ISO Prolog standard. In particular, the ASCII symbol characters glue into compound atoms, giving the familiar operator tokens `==`, `=..`, `:-` and so on. Beyond ASCII, all Unicode symbol characters (general categories `Sm`, `Sc`, `Sk`, `So`) and the connector, dash, and other-punctuation classes (`Pc`, `Pd`, `Po`) are treated as *solo*: each forms an atom on its own and does not glue with adjacent symbols. This is a deliberate change from earlier releases, in which Unicode symbols glued into compound atoms in the same way as ASCII symbols; the change ensures that characters such as ≤, €, and · keep their per-character meaning. Operators built from Unicode symbols must be declared explicitly with [op/3](operators.html#op/3). Numeric characters of other type (general category `No`, e.g. fractions and circled digits) are not part of the identifier set; only the explicitly listed super- and subscript digits extend identifiers.

- *Brackets (paired delimiters)*  
  The opening / closing punctuation classes `Ps` and `Pe` form *bracket pairs*: an opening character followed by a *Prolog term* and the matching closing character reads as a unary compound whose functor is the two delimiter characters joined. This is the same shape as `{Term}` becoming `'{}'(Term)`, generalised to the full Unicode `Ps`/`Pe` set (64 pairs, sourced from Unicode `BidiMirroring.txt` filtered by general category). Operators inside brackets are honoured; nesting works as expected. Mismatched or stray closes raise `syntax_error`. The analogy with `{}` is complete: an empty pair, optionally containing only layout, reads as the two-character atom rather than a compound; that atom followed by `(` is a functor; and on output `'<open><close>'(X)` is written as `<open>X<close>` and the bare atom unquoted, both subject to the `brace_terms(true)` write option.

- *Quotes (paired literal-text delimiters)*  
  The initial / final quotation classes `Pi` and `Pf` form *quote pairs*: an opening character followed by literal text and the matching closing character reads as a unary compound whose functor is the two delimiter characters joined and whose argument is the contained text in the form selected by [double_quotes](flags.html#flag:double_quotes) (string by default; also atom, codes, or chars). The contained text is *not* parsed as a Prolog term; escape sequences (`\n`, `\uXXXX`, ...) are processed as in ASCII quoted strings. For example, with [double_quotes](flags.html#flag:double_quotes) set to `string`, the source text «hello, world» reads as a compound whose functor is the two-character atom «» and whose single argument is the string `"hello, world"`. The quote pairs come from `Pi`/`Pf` entries of `BidiMirroring.txt` (8 pairs) plus the standard left/right curly quote pairs U+2018/U+2019 and U+201C/U+201D, which are absent from `BidiMirroring.txt`. Mismatched closes and unmatched opens raise `syntax_error`.

The features above let source text contain Unicode without escapes:

``` code
p(X⁰, X) :-                   % Superscript variable profile for
    q(X⁰, X¹),                % threaded variables.
    r(X¹, X).

?- atom_number('१२३', N).     % Devanagari Nd via atom_number/2
N = 123.

?- atom_codes(≤, Cs).         % Unicode symbol stays solo
Cs = [8804].

?- term_string(T, "⟨a, b⟩"),  % bracket pair (Ps/Pe)
   display(T).
⟨⟩(','(a,b))
T = ⟨a, b⟩.

?- term_string(T, "⟨ ⟩").     % empty pair is the atom
T = ⟨⟩.

?- term_string(T, "«hello").  % quote pair (Pi/Pf)
T = '«»'("hello").
```

#### 2.15.1.10 Singleton variable checking

A *singleton variable* is a variable that appears only one time in a clause. It can always be replaced by `_`, the *anonymous* variable. In some cases, however, people prefer to give the variable a name. As mistyping a variable is a common mistake, Prolog systems generally give a warning (controlled by [style_check/1](debugger.html#style_check/1)) if a variable is used only once. The system can be informed that a variable is meant to appear once by *starting* it with an underscore, e.g., `_Name`. Please note that any variable, except plain `_`, shares with variables of the same name. The term `t(_X, _X)` is equivalent to `t(X, X)`, which is *different* from `t(_, _)`.

As Unicode requires variables to start with an underscore in many languages, this schema has been extended.^(34After a proposal by Richard O'Keefe.) First we define the two classes of named variables.

- *Named singleton variables*  
  Named singletons start with a double underscore (`__`) or a single underscore followed by an uppercase letter, e.g., `__var` or `_Var`.
- *Normal variables*  
  All other variables are‘normal’variables. Note this makes `_var` a normal variable.^(35Some Prolog dialects write variables this way.)

Any normal variable appearing exactly once in the clause *and* any named singleton variables appearing more than once are reported. Below are some examples with warnings in the right column. Singleton messages can be suppressed using the [style_check/1](debugger.html#style_check/1) directive.

Finally, variables named `_<digit>` are never subject to style checking. These variable names are emitted by [write/1](termrw.html#write/1) and friends. This exception can also be used to pass singletons to [debug/3](debug.html#debug/3). Passing singletons to [debug/3](debug.html#debug/3) is otherwise problematic as using a normal variable results in a singleton warning when optimization removes the [debug/3](debug.html#debug/3) statement while using a named anonymous variable results in a *multiton* warning. For example:

``` code
p(X) :-
    q(X,_0Y),
    debug(demo, 'q/2 says ~p', [_0Y]).
```

|  |  |
|----|----|
| test(\_). |  |
| test(\_a). | Singleton variables: \[\_a\] |
| test(A). | Singleton variables: \[A\] |
| test(\_12). |  |
| test(\_A). |  |
| test(\_\_a). |  |
| test(\_, \_). |  |
| test(\_a, \_a). |  |
| test(\_\_a, \_\_a). | Singleton-marked variables appearing more than once: \[\_\_a\] |
| test(\_A, \_A). | Singleton-marked variables appearing more than once: \[\_A\] |
| test(A, A). |  |

**Semantic singletons**

Starting with version 6.5.1, SWI-Prolog has *syntactic singletons* and *semantic singletons*. The first are checked by [read_clause/3](termrw.html#read_clause/3) (and [read_term/3](termrw.html#read_term/3) using the option `singletons(warning)`). The latter are generated by the compiler for variables that appear alone in a *branch*. For example, in the code below the variable `X` is not a *syntactic* singleton, but the variable `X` does not communicate any bindings and replacing `X` with `_` does not change the semantics.

``` code
test :-
        (   test_1(X)
        ;   test_2(X)
        ).
```
