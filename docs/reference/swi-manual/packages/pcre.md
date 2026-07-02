SWI-Prolog Regular Expression library

Jan Wielemaker and Peter Ludemann  
VU University Amsterdam  
The Netherlands  
E-mail: [J.Wielemaker@vu.nl](mailto:J.Wielemaker@vu.nl)

Abstract

The library `library(pcre)` provides access to Perl Compatible Regular Expressions.

# Table of Contents

[1 Motivation](#sec:1)

[2 library(pcre): Perl compatible regular expression matching for SWI-Prolog](#sec:2)

## 1 Motivation

The core facility for string matching in Prolog is provided by DCG (*Definite Clause Grammars*). Using DCGs is typically more verbose but gives reuse, modularity, readability and mixing with arbitrary Prolog code in return. Supporting regular expressions has some advantages: (1) in simple cases, the terse specification of a regular expression is more comfortable; (2) many programmers are familar with them; and (3) regular expressions are part of domain specific languages one may wish to implement in Prolog, e.g., SPARQL.

There are roughly three options for adding regular expressions to Prolog. One is to simply interpret them in Prolog. Given Prolog's unification and backtracking facilities this is remarkable simple and performs quite reasonably. Still, implementing all facilities of modern regular expression engines requires significant effort. Alternatively, we can *compile* them into DCGs. This brings terse expressions to DCGs while staying in the same framework. The disadvantage is that regular expressions become programs that are hard to work with, making this approach less attractive for applications that potentially execute many different regular expressions. The final option is to wrap an existing regular expression engine. This provides access to a robust implementation for which we only have to document the Prolog binding. That is the option taken by library `library(pcre)`.

## 2 library(pcre): Perl compatible regular expression matching for SWI-Prolog

See also  
‘man pcre2api\` or [https://www.pcre.org/current/doc/html/pcre2api.html](https://www.pcre.org/current/doc/html/pcre2api.html) for details of the PCRE2 syntax and options.

This module provides an interface to the [PCRE2](http://www.pcre.org/) (Perl Compatible Regular Expression) library. This Prolog interface provides an almost complete wrapper around PCRE2 (the successor to PCRE) with as much backward compatibility to PCRE as possible, because the original implementation was for PCRE (also known as PCRE1).

Regular expressions are created from a pattern and options and represented as a SWI-Prolog *blob*. This implies they are subject to (atom) garbage collection. Compiled regular expressions can safely be used in multiple threads. Most predicates accept both an explicitly compiled regular expression, a pattern, or a term Pattern/Flags. The semantics of the pattern can be additionally modified by options. In the latter two cases a regular expression *blob* is created and stored in a cache. The cache can be cleared using [re_flush/0](#re_flush/0).

Most of the predicates in this library take both a regular expression represented as a string with optional flags, e.g., `'aap'/i` or a *compiled regular* expression. If a string (+flags) alternative is used, the library maintains a cache of compiled regular expressions. See also [re_flush/0](#re_flush/0). The library can be asked to rewrite the [re_match/2](#re_match/2) and [re_match/3](#re_match/3) goals to use inlined compiled regular expression objects using

``` code
:- set_prolog_flag(re_compile, true).
```

This has some consequences:

- Performance is considerable better.
- Compiled regular expressions are currently incompatible with *Quick Load Files* (\`.qlf\`, see qcompile/1) and *Saved States* (see qsave_program/2 and the `-c` command line option.
- Debugging may be harder.

\[semidet\]**re_match**(`+Regex, +String`)  
\[semidet\]**re_match**(`+Regex, +String, +Options`)  
Succeeds if `String` matches `Regex`. For example:

``` code
?- re_match("^needle"/i, "Needle in a haystack").
true.
```

Defined `Options` are given below. For details, see the PCRE documentation. If an option is repeated, the first value is used and subsequent values are ignored. Unrecognized options are ignored. Unless otherwise specified, boolean options default to `false`.

If `Regex` is a text pattern (optionally with flags), then any of the `Options` for [re_compile/3](#re_compile/3) can be used, in addition to the `Options` listed below. If `Regex` is the result of [re_compile/3](#re_compile/3), then only the following execution-time `Options` are recognized and any others are ignored. Some options may not exist on your system, depending on the PCRE2 version and how it was built - these unsupported options are silently ignored.

- `start(From)` Start at the given character index
- `anchored(Bool)` If `true`, match only at the first position
- `bol(Bool)` `String` is the beginning of a line (default `true`) - affects behavior of circumflex metacharacter (`^`).
- `empty(Bool)` An empty string is a valid match (default `true`)
- `empty_atstart(Bool)` An empty string at the start of the subject is a valid match (default `true`)
- `eol(Bool)` `String` is the end of a line - affects behavior of dollar metacharacter (`$`) (default `true`).
- `newline(Mode)` If `any`, recognize any Unicode newline sequence, if `anycrlf`, recognize CR, LF, and CRLF as newline sequences, if `cr`, recognize CR, if `lf`, recognize LF, if `crlf` recognize CRLF as newline. The default is determined by how PCRE was built, and can be found by `re_config(newline2(NewlineDefault))`.
- `newline2(Mode)` - synonym for `newline(Mode)`.
- `utf_check(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) You should not need this because SWI-Prolog ensures that the UTF8 strings are valid, so the default is `false`.
- `endanchored(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `partial_soft(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `partial_hard(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `dfa_restart(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `dfa_shortest(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)

[TABLE]

\[semidet\]**re_matchsub**(`+Regex, +String, -Sub:dict`)  
\[semidet\]**re_matchsub**(`+Regex, +String, -Sub:dict, +Options`)  
Match `String` against `Regex`. On success, `Sub` is a dict containing integer keys for the numbered capture group and atom keys for the named capture groups. The entire match string has the key `0`. The associated value is determined by the `capture_type(Type)` option passed to [re_compile/3](#re_compile/3), or by flags if `Regex` is of the form Pattern/Flags; and may be specified at the level of individual captures using a naming convention for the caption name. See [re_compile/3](#re_compile/3) for details.

The example below exploits the typed groups to parse a date specification:

``` code
?- re_matchsub("(?<date> (?<year_I>(?:\\d\\d)?\\d\\d) -
                (?<month_I>\\d\\d) - (?<day_I>\\d\\d) )"/x,
               "2017-04-20", Sub, []).
Sub = re_match{0:"2017-04-20", date:"2017-04-20",
               day:20, month:4, year:2017}.
```

|  |  |
|----|----|
| `Both` | compilation and execution options are processed. See [re_compile/3](#re_compile/3) and [re_match/3](#re_match/3) for the set of options. In addition, some compilation options may passed as `/Flags` to `Regex` - see [re_match/3](#re_match/3) for the list of flags. |
| `Regex` | See [re_match/2](#re_match/2) for a description of this argument. |

\[semidet\]**re_foldl**(`:Goal, +Regex, +String, ?V0, ?V, +Options`)  
Fold all matches of `Regex` on `String`. Each match is represented by a dict as specified for [re_matchsub/4](#re_matchsub/4). `V0` and `V` are related using a sequence of invocations of `Goal` as illustrated below.

``` code
call(Goal, Dict1, V0, V1),
call(Goal, Dict2, V1, V2),
...
call(Goal, Dictn, Vn, V).
```

This predicate is used to implement [re_split/4](#re_split/4) and [re_replace/4](#re_replace/4). For example, we can count all matches of a `Regex` on `String` using this code:

``` code
re_match_count(Regex, String, Count) :-
    re_foldl(increment, Regex, String, 0, Count, []).

increment(_Match, V0, V1) :-
    V1 is V0+1.
```

After which we can query

``` code
?- re_match_count("a", "aap", X).
X = 2.
```

Here is an example `Goal` for extracting all the matches with their offsets within the string:

``` code
range_match(Dict, StringIndex-[MatchStart-Substring|List], StringIndex-List) :-
    Dict.(StringIndex.index) = MatchStart-MatchLen,
    sub_string(StringIndex.string, MatchStart, MatchLen, _, Substring).
```

And can be used with this query (note the `capture_type(range)` option, which is needed by range_match/3, and `greedy(false)` to invert the meaning of `*?`):

``` code
?- String = "{START} Mary {END} had a {START} little lamb {END}",
   re_foldl(range_match,
            "{START} *?(?<piece>.*) *?{END}",
            String, _{string:String,index:piece}-Matches, _-[],
            [capture_type(range),greedy(false)]).
Matches = [8-"Mary", 33-"little lamb"].
```

\[det\]**re_split**(`+Pattern, +String, -Splits:list`)  
\[det\]**re_split**(`+Pattern, +String, -Splits:list, +Options`)  
Split `String` using the regular expression `Pattern`. `Splits` is a list of strings holding alternating matches of `Pattern` and skipped parts of the `String`, starting with a skipped part. The `Splits` lists ends with a string of the content of `String` after the last match. If `Pattern` does not appear in `String`, `Splits` is a list holding a copy of `String`. This implies the number of elements in `Splits` is *always* odd. For example:

``` code
?- re_split("a+", "abaac", Splits, []).
Splits = ["","a","b","aa","c"].
?- re_split(":\\s*"/n, "Age: 33", Splits, []).
Splits = ['Age', ': ', 33].
```

|  |  |
|----|----|
| `Pattern` | is the pattern text, optionally follows by /Flags. Similar to [re_matchsub/4](#re_matchsub/4), the final output type can be controlled by a flag `a` (atom), `s` (string, default) or `n` (number if possible, atom otherwise). |

\[det\]**re_replace**(`+Pattern, +With, +String, -NewString`)  
\[det\]**re_replace**(`+Pattern, +With, +String, -NewString, +Options`)  
Replace matches of the regular expression `Pattern` in `String` with `With` (possibly containing references to captured substrings).

Throws an error if `With` uses a name that doesn't exist in the `Pattern`.

[TABLE]

\[det\]**re_compile**(`+Pattern, -Regex, +Options`)  
Compiles `Pattern` to a `Regex` *blob* of type `regex` (see blob/2). Defined `Options` are given below. Please consult the [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) for details. If an option is repeated, the first value is used and subsequent values are ignored. Unrecognized options are ignored. Unless otherwise specified, boolean options default to `false`. Some options may not exist on your system, depending on the PCRE2 version and how it was built - these unsupported options are silently ignored.

The various matching predicates can take either a `Regex` *blob* or a string pattern; if they are given a string pattern, they call [re_compile/3](#re_compile/3) and cache the result; so, there is little reason to use [re_compile/3](#re_compile/3) directly.

- `anchored(Bool)` If `true`, match only at the first position
- `auto_capture(Bool)` Enable use of numbered capturing parentheses. (default `true`)
- `bsr(Mode)` If `anycrlf`, `\`R only matches CR, LF or CRLF; if `unicode`, `\`R matches all Unicode line endings.
- `bsr2(Mode)` - synonym for `bsr(Mode)`.
- `caseless(Bool)` If `true`, do caseless matching.
- `compat(With)` Error - PCRE1 had `compat(javascript)` for JavaScript compatibility, but PCRE2 has removed that.
- `dollar_endonly(Bool)` If `true`, \$ not to match newline at end
- `dotall(Bool)` If `true`, . matches anything including NL
- `dupnames(Bool)` If `true`, allow duplicate names for subpatterns
- `extended(Bool)` If `true`, ignore white space and \# comments
- `firstline(Bool)` If `true`, force matching to be before newline
- `greedy(Bool)` If `true`, operators such as `+` and `*` are greedy unless followed by `?`; if `false`, the operators are not greedy and `?` has the opposite meaning. It can also beset by a‘(?U)\` within the pattern - see the [PCRE2 pattern internal option setting documentation](https://www.pcre.org/current/doc/html/pcre2pattern.html\#SEC13) for details and note that the PCRE2 option is `UNGREEDY`, which is the inverse of this packages `greedy` options. (default `true`)
- `compat(With)` Raises an errr - PCRE1 had `compat(javascript)` for JavaScript compatibility, but PCRE2 has removed that option . Consider using the `alt_bsux` and `extra_alt_bsux` options.
- `multiline(Bool)` If `true`, `^` and \$ match newlines within data
- `newline(Mode)` If `any`, recognize any Unicode newline sequence; if `anycrlf` (default), recognize CR, LF, and CRLF as newline sequences; if `cr`, recognize CR; if `lf`, recognize LF; `crlf` recognize CRLF as newline; if `nul`, recognize the NULL character (0x00) as newline.
- `newline2(Mode)` - synonym for `newline(Mode)`.
- `ucp(Bool)` If `true`, use Unicode properties for `\`d, `\`w, etc.
- `utf_check(Bool)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) You should not need this because SWI-Prolog ensures that the UTF8 strings are valid,
- `endanchored(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `allow_empty_class(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `alt_bsux(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `auto_callout(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `match_unset_backref(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `never_ucp(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `never_utf(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `auto_possess(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) (default `true`)
- `dotstar_anchor(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) (default `true`)
- `start_optimize(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) (default `true`)
- `utf(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `never_backslash_c(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `alt_circumflex(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `alt_verbnames(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `use_offset_limit(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `extended_more(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `literal(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `match_invalid_utf(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `jit_complete(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `jit_partial_soft(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `jit_partial_hard(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `jit_invalid_utf(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)
- `jit(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html) (default `true`)
- `copy_matched_subject(boolean)` - see [PCRE2 API documentation](https://www.pcre.org/current/doc/html/pcre2api.html)

In addition to the options above that directly map to PCRE flags the following options are processed:

- `optimise(Bool)` or `optimize(Bool)` Turns on the JIT compiler for additional optimization that greatly that speeds up the matching performance of many patterns. (Note that he meaning has changed slightly from the PCRE1 implementation
- PCRE2 always optimises where possible; this is an additional optimisation.)
- `capture_type(+Type)` How to return the matched part of the input and possibly captured groups in there. Possible values are:
  **string**  
  Return the captured string as a string (default).

  **atom**  
  Return the captured string as an atom.

  **range**  
  Return the captured string as a pair `Start-Length`. Note that we use `Start-Length` rather than the more conventional `Start-End` to allow for immediate use with sub_atom/5 and sub_string/5.

  **term**  
  Parse the captured string as a Prolog term. This is notably practical if you capture a number.

The `capture_type` specifies the default for this pattern. The interface supports a different type for each *named* group using the syntax‘(?\<name_T\>...)\`, where `T` is one of `S` (string), `A` (atom), `I` (integer), `F` (float), `N` (number), `T` (term) and `R` (range). In the current implementation `I`, `F` and `N` are synonyms for `T`. Future versions may act different if the parsed value is not of the requested numeric type.

Note that [re_compile/3](#re_compile/3) does not support the `Pattern`/Flags form that is supported by [re_match/3](#re_match/3), [re_replace/4](#re_replace/4), etc.; the `Pattern` must be text and all compile options specified in `Options`.

**re_flush**  
Clean pattern and replacement caches.

To be done  
Flush automatically if the cache becomes too large.

**re_config**(`?Term`)  
Extract configuration information from the pcre library. `Term` is of the form `Name(Value)`. Name is derived from the `PCRE_CONFIG_*` constant after removing `PCRE_CONFIG_` and mapping the name to lower case, e.g. `utf8`, `unicode_properties`, etc. Value is a Prolog boolean, integer, or atom. For boolean (1 or 0) values, `true` or `false` is returned.

[re_config/1](#re_config/1) will backtrack through all the possible configuration values if its argument is a variable. If an unknown option is specified, [re_config/1](#re_config/1) fails.

Non-compatible changes between PCRE1 and PCRE2 because numeric values changed: `bsr` and `newline` have been replaced by `bsr2` and `newline2`:

- `bsr2` - previously `bsr` returned 0 or 1; now returns `unicode` or `anycrlf`
- `newline2` - previously `newline` returned an integer, now returns `cr`, `lf`, `crlf`, `any`, `anycrlf`, `nul`

`Term` values are as follows. Some values might not exist, depending on the version of PCRE2 and the options it was built with.

- bsr2 The character sequences that the `\R` escape sequence matches by default. Replaces `bsr` option from PCRE1, which is not compatible.

- compiled_widths An integer whose lower bits indicate which code unit widths were selected when PCRE2 was built. The 1-bit indicates 8-bit support, and the 2-bit and 4-bit indicate 16-bit and 32-bit support, respectively. The 1 bit should always be set because the wrapper code requires 8 bit support.

- depthlimit

- heaplimit

- jit `true` if just-in-time compiling is available.

- jittarget A string containing the name of the architecture for which the JIT compiler is configured. e.g.,’x86 64bit (little endian + unaligned)’.

- linksize

- matchlimit

- never_backslash_c

- newline2 An atom whose value specifies the default character sequence that is recognized as meaning "newline" (`cr`, `lf`, `crlf`, `any`, `anycrlf`, `nul`). Replaces `newline` option from PCRE1, which is not compatible.

- parenslimit

- stackrecurse

- unicode Always `true`

- unicode_version The unicode version as an atom, e.g.’12.1.0’.

- utf8 - synonym for `unicode`

- parens_limit

- version The version information as an atom, containing the PCRE version number and release date, e.g.’10.34 2019-11-21’.

  For backwards compatibility with PCRE1, the following are accepted, but are deprecated:

  - `utf8` - synonym for `unicode`
  - `link_size` - synonym for `linksize`
  - `match_limit` - synonym for `matchlimit`
  - `parens_limit` - synonym for `parenslimit`
  - `unicode_properties` - always true

  The following have been removed because they don't exist in PCRE2 and don't seem to have any meaningful use in PCRE1:

  - `posix_malloc_threshold`
  - `match_limit_recursion`

# Index

?  
[re_compile/3](#re_compile/3)  
[re_config/1](#re_config/1)  
[re_flush/0](#re_flush/0)  
[re_foldl/6](#re_foldl/6)  
[re_match/2](#re_match/2)  
[re_match/3](#re_match/3)  
[re_matchsub/3](#re_matchsub/3)  
[re_matchsub/4](#re_matchsub/4)  
[re_replace/4](#re_replace/4)  
[re_replace/5](#re_replace/5)  
[re_split/3](#re_split/3)  
[re_split/4](#re_split/4)  
