SWI-Prolog Unicode library

Jan Wielemaker  
VU University, Amsterdam  
The Netherlands  
E-mail: [J.Wielemaker@cs.vu.nl](mailto:J.Wielemaker@cs.vu.nl)

Abstract

This package wraps [utf8proc](https://github.com/JuliaStrings/utf8proc) unicode routines. This library provides the four unicode normalization forms (NFC, NFD, NFKC, NFKD) as well as access to the Unicode character properties.

# Table of Contents

[1 library(unicode): Unicode string handling](#sec:1)

[2 library(unicode_security): Unicode security helpers (UTS \#39, UAX \#24)](#sec:2)

[3 License](#sec:3)

## 1 library(unicode): Unicode string handling

This library wraps the [utf8proc](https://github.com/JuliaStrings/utf8proc) library, giving Prolog code access to Unicode character properties, string normalization, case folding, and grapheme-cluster iteration.

Three levels of API are provided:

1.  **Normalization**: [unicode_nfd/2](#unicode_nfd/2), [unicode_nfc/2](#unicode_nfc/2), [unicode_nfkd/2](#unicode_nfkd/2), [unicode_nfkc/2](#unicode_nfkc/2) implement the four standard Unicode normalization forms (NFD, NFC, NFKD, NFKC; see UAX#15). [unicode_nfkc_casefold/2](#unicode_nfkc_casefold/2) combines NFKC with case folding for caseless identifier matching (see UAX#31).

2.  **Per-codepoint properties**: [unicode_property/2](#unicode_property/2) queries the Unicode property database (general category, bidi class, decomposition type, display width, case mappings, grapheme boundary class, ...).

3.  **Mixed string-level transformations**: [unicode_map/3](#unicode_map/3) is the workhorse; it accepts a list of flags chosen from a fixed set and performs the corresponding composition of decompose / compose / strip / lump / case-fold / grapheme-boundary-mark operations in a single pass. [unicode_casefold/2](#unicode_casefold/2) is a convenience wrapper.

    Grapheme clusters (user-perceived characters) can be iterated with [atom_graphemes/2](#atom_graphemes/2) and [string_graphemes/2](#string_graphemes/2).

    Loading this library also installs a Unicode NFC normalisation hook into the SWI-Prolog kernel. The kernel's `unicode_atoms` policy (Prolog flag, stream property and `read_term/2,3` option) uses this hook for its `nfc` and `error` modes; without the library loaded those modes raise `existence_error(hook, unicode_normalize)`. The kernel's quoted-write rule that force-quotes atoms containing combining marks is independent and works even without this library loaded.

    Lump handling:

    ``` code
    U+0020      <-- all space characters (general category Zs)
    U+0027  '   <-- left/right single quotation mark U+2018..2019,
                    modifier letter apostrophe U+02BC,
                    modifier letter vertical line U+02C8
    U+002D  -   <-- all dash characters (general category Pd),
                    minus U+2212
    U+002F  /   <-- fraction slash U+2044,
                    division slash U+2215
    U+003A  :   <-- ratio U+2236
    U+003C  <   <-- single left-pointing angle quotation mark U+2039,
                    left-pointing angle bracket U+2329,
                    left angle bracket U+3008
    U+003E  >   <-- single right-pointing angle quotation mark U+203A,
                    right-pointing angle bracket U+232A,
                    right angle bracket U+3009
    U+005C  \   <-- set minus U+2216
    U+005E  ^   <-- modifier letter up arrowhead U+02C4,
                    modifier letter circumflex accent U+02C6,
                    caret U+2038,
                    up arrowhead U+2303
    U+005F  _   <-- all connector characters (general category Pc),
                    modifier letter low macron U+02CD
    U+0060  `   <-- modifier letter grave accent U+02CB
    U+007C  |   <-- divides U+2223
    U+007E  ~   <-- tilde operator U+223C
    ```

    @see [http://www.unicode.org/reports/tr15/](http://www.unicode.org/reports/tr15/) (UAX#15 Normalization) @see [http://www.unicode.org/reports/tr29/](http://www.unicode.org/reports/tr29/) (UAX#29 Grapheme Clusters) @see [http://www.unicode.org/reports/tr31/](http://www.unicode.org/reports/tr31/) (UAX#31 Identifiers) @see [https://github.com/JuliaStrings/utf8proc](https://github.com/JuliaStrings/utf8proc)

\[det\]**unicode_map**(`+In, -Out, +Options`)  
Perform a Unicode mapping on `In`, returning `Out`. `Options` is a list that may contain any combination of the flags below; a call is roughly equivalent to `utf8proc_map(In, Options)` in the C API.

**stable**  
Respect Unicode versioning stability --- the result does not depend on which (recent) version of Unicode is in use.

**compat**  
Use compatibility decomposition (i.e. formatting information is lost).

**compose**  
Produce a composed result (e.g. NFC or NFKC, depending on the presence of `compat`).

**decompose**  
Produce a decomposed result (NFD/NFKD).

**ignore**  
Strip "default ignorable" characters (e.g. soft hyphen, zero-width space).

**rejectna**  
Raise an error instead of returning output when the input contains unassigned code points.

**nlf2ls**  
Convert all NLF-sequences (LF, CRLF, CR, NEL) to U+2028 LINE SEPARATOR.

**nlf2ps**  
Convert all NLF-sequences to U+2029 PARAGRAPH SEPARATOR.

**nlf2lf**  
Convert all NLF-sequences to U+000A LINE FEED.

**stripcc**  
Strip or convert control characters. NLF-sequences become a space, except if one of the NLF-conversion flags is set; HT and FF are treated as NLF in this case. All other control characters are removed.

**casefold**  
Apply Unicode case folding (for caseless comparison).

**charbound**  
Insert a U+00FF byte at the beginning of every grapheme cluster (UAX#29). The result can be split on 0xFF to recover individual graphemes; [atom_graphemes/2](#atom_graphemes/2) wraps this pattern.

**lump**  
Normalise typographic variants to their ASCII equivalents (see module header for the full list). Combined with `nlf2lf`, paragraph and line separators become U+000A as well.

**stripmark**  
Strip all combining marks (non-spacing, spacing, enclosing). Must be combined with `compose` or `decompose`.

\[det\]**unicode_nfd**(`+In, -Out`)  
Characters in `In` are decomposed by canonical equivalence (NFD). Precomposed characters expand into base + combining marks. For example U+00C5 (LATIN CAPITAL LETTER A WITH RING ABOVE) becomes the two-code sequence `A` + U+030A.

See also  
[http://www.unicode.org/reports/tr15/](http://www.unicode.org/reports/tr15/)

\[det\]**unicode_nfc**(`+In, -Out`)  
Characters in `In` are decomposed and then recomposed by canonical equivalence (NFC). Precomposed code points are preferred; for example `A` + U+030A becomes the single code point U+00C5.

See also  
[http://en.wikipedia.org/wiki/Unicode_equivalence\\Normal_forms](http://en.wikipedia.org/wiki/Unicode_equivalence\#Normal_forms)

\[det\]**unicode_nfkd**(`+In, -Out`)  
Characters in `In` are decomposed by compatibility equivalence (NFKD). Compatibility decomposition expands presentation forms (ligatures, subscripts, fullwidth letters) into their base forms; for example U+FB03 (LATIN SMALL LIGATURE FFI) becomes‘f f i\`.

\[det\]**unicode_nfkc**(`+In, -Out`)  
Characters in `In` are decomposed by compatibility equivalence, then recomposed by canonical equivalence (NFKC).

\[det\]**unicode_nfkc_casefold**(`+In, -Out`)  
Equivalent to [unicode_nfkc/2](#unicode_nfkc/2) followed by [unicode_casefold/2](#unicode_casefold/2) done in a single pass. This is the normalisation form recommended by UAX#31 for caseless identifier matching. For example German `'Strasse'` written with U+00DF (LATIN SMALL LETTER SHARP S) maps to `'strasse'`, and U+FB03 (LATIN SMALL LIGATURE FFI) maps to `'ffi'`.

\[det\]**unicode_casefold**(`+In, -Out`)  
`Out` is the case-folded form of `In`. Use this for caseless comparison that does not require NFKC; otherwise prefer [unicode_nfkc_casefold/2](#unicode_nfkc_casefold/2).

\[nondet\]**unicode_property**(`?Code, ?Property`)  
Query the Unicode character database for `Code`. `Code` is an integer code point (0 .. 0x10FFFF) or a single-character atom; `Property` is a term of the form Name(Value) drawn from the list below.

This predicate is a thin wrapper over utf8proc's property struct, so its vocabulary matches the utf8proc documentation. In the modes (+,?) and (-,?) the predicate enumerates properties for the given code (or the code for the given property); in (+,+) it is a deterministic test.

Supported properties:

**category**(`Atom`)  
Unicode general category. `Atom` is one of `Cc`, `Cf`, `Cn`, `Co`, `Cs`, `Ll`, `Lm`, `Lo`, `Lt`, `Lu`, `Mc`, `Me`, `Mn`, `Nd`, `Nl`, `No`, `Pc`, `Pd`, `Pe`, `Pf`, `Pi`, `Po`, `Ps`, `Sc`, `Sk`, `Sm`, `So`, `Zl`, `Zp`, `Zs`. When querying, the single capital letter of a subcategory stands for all its subcategories; e.g.

``` code
?- unicode_property(0'A, category('L')).
true.
```

**combining_class**(`Integer`)  
Canonical combining class (0 for base characters, 230 for accents above, etc.).

**bidi_class**(`Atom`)  
Bidirectional class. One of `l`, `lre`, `lro`, `r`, `al`, `rle`, `rlo`, `pdf`, `en`, `es`, `et`, `an`, `cs`, `nsm`, `bn`, `b`, `s`, `ws`, `on`.

**bidi_mirrored**(`Bool`)  
`true` if the character is mirrored for bidi (parentheses, brackets, math operators, ...).

**decomp_type**(`Atom`)  
Compatibility decomposition type. One of `font`, `nobreak`, `initial`, `medial`, `final`, `isolated`, `circle`, `super`, `sub`, `vertical`, `wide`, `narrow`, `small`, `square`, `fraction`, `compat`. Fails when there is no decomposition.

**ignorable**(`Bool`)  
`true` if the character is a "default ignorable" code point.

**boundclass**(`Atom`)  
UAX#29 grapheme-cluster break class. One of `start`, `other`, `cr`, `lf`, `control`, `extend`, `l`, `v`, `t`, `lv`, `lvt`, `regional_indicator`, `spacingmark`, `prepend`, `zwj`, `extended_pictographic`, `e_zwg`.

**width**(`Integer`)  
Display width in fixed-width cells, 0..3. Zero for combining marks and control characters, 1 for most, 2 for "wide" characters (CJK, emoji).

**ambiguous_width**(`Bool`)  
`true` if the character has East-Asian Ambiguous width --- normally one column, but two in a legacy CJK context.

**uppercase**(`Code`)  
**lowercase**(`Code`)  
**titlecase**(`Code`)  
Single-code-point case mapping. Fails when the code point has no mapping of that kind (e.g. `unicode_property(0'A, uppercase(_))` fails because `'A'` is already upper-case). For characters whose case mapping produces more than one code point (e.g. U+00DF LATIN SMALL LETTER SHARP S maps to "SS"), use [unicode_map/3](#unicode_map/3) with the `[casefold]` option or [unicode_casefold/2](#unicode_casefold/2) for a full string-level transformation.

**indic_conjunct_break**(`Atom`)  
Indic_Conjunct_Break property (Unicode 15+; UAX#44). One of `none`, `linker`, `consonant`, `extend`. Used by the grapheme-cluster-break algorithm for Devanagari, Bengali, etc.

See also  
[http://www.unicode.org/reports/tr44/](http://www.unicode.org/reports/tr44/) (Unicode property database)

\[det\]**atom_graphemes**(`?Atom, ?Graphemes`)  
Relate `Atom` to a list of its grapheme clusters. Grapheme clusters are "user-perceived characters" as defined by UAX#29 --- e.g. the precomposed U+00E9 (LATIN SMALL LETTER E WITH ACUTE) and the decomposed sequence `e` + U+0301 are both one grapheme, an emoji ZWJ sequence such as `MAN + ZWJ + WOMAN + ZWJ + GIRL` is one grapheme, and a regional-indicator pair (e.g. U+1F1F3 U+1F1F1, rendered as the Dutch flag) is one grapheme.

In the forward mode (+`Atom`, ?`Graphemes`), `Atom` is decomposed into a list of atoms, each covering one cluster. In the reverse mode (?`Atom`, +`Graphemes`), the elements of `Graphemes` are concatenated into `Atom`. Both arguments instantiated means both modes run and the result must agree.

``` code
?- atom_codes(A, [0'c, 0'a, 0'f, 0'e, 0x0301]),
   atom_graphemes(A, Gs).
Gs = [c, a, f, G],
atom_codes(G, [0'e, 0x0301]).

?- atom_graphemes(A, [a, b, c]).
A = abc.
```

See also  
[string_graphemes/2](#string_graphemes/2) for the string analogue.

\[det\]**string_graphemes**(`?String, ?Graphemes`)  
As [atom_graphemes/2](#atom_graphemes/2), but the elements of `Graphemes` are strings.

\[det\]**unicode_version**(`-Version`)  
`Version` is an atom describing the Unicode version implemented by the linked utf8proc library, e.g. `'15.1.0'`. This drives the normalisation, case-folding and grapheme-cluster predicates in this module, and may differ from the Unicode version of the SWI-Prolog source syntax classifier reported by the read-only Prolog flag **unicode_syntax_version**.

\[semidet\]**unicode_codepoint_valid**(`+Code`)  
True when `Code` is a non-negative integer that is a valid and *assigned* Unicode code point. Unassigned code points (general category `Cn`), surrogate halves (`Cs`), and integers outside‘0..0x10FFFF\` all fail.

## 2 library(unicode_security): Unicode security helpers (UTS \#39, UAX \#24)

This library implements helpers from [UTS \#39 (Unicode Security Mechanisms)](https://www.unicode.org/reports/tr39/) and the script properties of [UAX \#24](https://www.unicode.org/reports/tr24/). It is intended for linters, identifier validators and any code that needs to reason about confusable look-alike text or mixed-script identifiers. It does **not** alter the Prolog reader; UTS \#39 is deliberately a library-level facility.

The library ships its own UCD-derived tables and is independent of `library(unicode)` (which wraps libutf8proc for normalisation and per-code-point properties). See `etc/gen_uts39.pl` in the package directory to regenerate the tables on a Unicode-version bump.

Predicates fall into three groups:

- Per-code-point lookups: [unicode_script/2](#unicode_script/2), [unicode_script_extensions/2](#unicode_script_extensions/2), [unicode_identifier_status/2](#unicode_identifier_status/2), [unicode_identifier_type/2](#unicode_identifier_type/2).
- Skeleton and confusable test (UTS \#39 §4): [unicode_skeleton/2](#unicode_skeleton/2), [unicode_confusable/2](#unicode_confusable/2), [unicode_confusable/3](#unicode_confusable/3).
- String-level identifier checks (UTS \#39 §5): [unicode_resolved_scripts/2](#unicode_resolved_scripts/2), [unicode_restriction_level/2](#unicode_restriction_level/2).

\[semidet\]**unicode_script**(`+Code:integer, -Script:atom`)  
True when `Script` is the UAX \#24 Script_Property of `Code`. `Script` is a lower-case atom of the long property value (`latin`, `cyrillic`, `han`, `common`, `inherited`, ...). Fails for code points outside the Unicode range or with no entry in `Scripts.txt`.

\[semidet\]**unicode_script_extensions**(`+Code:integer, -Scripts:list(atom)`)  
`Scripts` is the sorted list of UAX \#24 Script_Extensions of `Code`. For most code points this is a singleton `[Script]`. Fails for code points outside the Unicode range and for code points with no entry in either `ScriptExtensions.txt` or `Scripts`.txt.

\[semidet\]**unicode_identifier_status**(`+Code:integer, -Status:atom`)  
Succeeds, unifying `Status` with `allowed`, when `Code` is listed in UTS \#39 `IdentifierStatus.txt` with status `Allowed`. Fails otherwise — per UTS \#39 every code point not listed there is Restricted by default; rather than return `restricted` for everything else, this predicate simply fails.

\[semidet\]**unicode_identifier_type**(`+Code:integer, -Types:list(atom)`)  
`Types` is the sorted list of UTS \#39 Identifier_Type atoms for `Code` (`recommended`, `inclusion`, `technical`, `obsolete`, `limited_use`, `exclusion`, `not_nfkc`, `not_xid`, `default_ignorable`, `deprecated`, `uncommon_use`). Fails for code points outside the Unicode range or with no entry in `IdentifierType.txt`.

\[det\]**unicode_skeleton**(`+Text, -Skeleton:atom`)  
Compute the UTS \#39 §4 skeleton of `Text`: apply NFD, substitute each code point with its `confusables.txt` prototype string, then apply NFD again. Two strings are confusable iff their skeletons compare equal.

\[semidet\]**unicode_confusable**(`+T1, +T2`)  
True when [unicode_skeleton/2](#unicode_skeleton/2) of `T1` and `T2` are equal.

\[semidet\]**unicode_confusable**(`+T1, +T2, +Options`)  
As [unicode_confusable/2](#unicode_confusable/2). `Options`:

**ignore_intentional**(`+Bool`)  
If `true`, skip the per-character substitution when the source and target form a pair listed in UTS \#39 `intentional.txt` (e.g. Latin A versus Greek capital Alpha). Default `false`.

\[det\]**unicode_resolved_scripts**(`+Text, -Scripts:list(atom)`)  
`Scripts` is the UTS \#39 §5.1 resolved augmented Script_Extensions set of `Text`: the intersection of `augscx(c)` over all non-Common/non-Inherited characters, with the augmentation rules for Han, Hiragana, Katakana, Hangul and Bopomofo applied. The empty list signals a mixed-script string.

\[det\]**unicode_restriction_level**(`+Text, -Level:atom`)  
Classify `Text` under UTS \#39 §5.2 at the most restrictive level for which it qualifies. `Level` is one of:

- `ascii_only` — every code point in U+0020..U+007E and Allowed.
- `single_script` — augmented resolved-script-set non-empty and every code point Allowed.
- `highly_restrictive` — covered by Latin plus one of `Hanb`, `Jpan` or `Kore` (UTS \#39 §5.1 augmented profiles).
- `moderately_restrictive` — covered by Latin plus a single non-Latin Recommended script (`Cyrl` or `Grek`).
- `minimally_restrictive` — every code point has Identifier_Type in `{recommended, inclusion}`.
- `unrestricted` — otherwise. A linter that walks source clauses and reports atoms with the confusability issues above is registered in `library(check)` itself (predicate list_confusable_identifiers/0); see the `library(check)` documentation for details.

## 3 License

Copyright (c) 2009 Public Software Group e. V., Berlin, Germany

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

This software distribution contains derived data from a modified version of the Unicode data files. The following license applies to that data:

COPYRIGHT AND PERMISSION NOTICE

Copyright (c) 1991-2007 Unicode, Inc. All rights reserved. Distributed under the Terms of Use in http://www.unicode.org/copyright.html.

Permission is hereby granted, free of charge, to any person obtaining a copy of the Unicode data files and any associated documentation (the "Data Files") or Unicode software and any associated documentation (the "Software") to deal in the Data Files or Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, and/or sell copies of the Data Files or Software, and to permit persons to whom the Data Files or Software are furnished to do so, provided that (a) the above copyright notice(s) and this permission notice appear with all copies of the Data Files or Software, (b) both the above copyright notice(s) and this permission notice appear in associated documentation, and (c) there is clear notice in each modified Data File or in the Software as well as in the documentation associated with the Data File(s) or Software that the data or software has been modified.

THE DATA FILES AND SOFTWARE ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THE DATA FILES OR SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall not be used in advertising or otherwise to promote the sale, use or other dealings in these Data Files or Software without prior written authorization of the copyright holder.

Unicode and the Unicode logo are trademarks of Unicode, Inc., and may be registered in some jurisdictions. All other trademarks and registered trademarks mentioned herein are the property of their respective owners.

# Index

?  
[atom_graphemes/2](#atom_graphemes/2)  
[string_graphemes/2](#string_graphemes/2)  
[unicode_casefold/2](#unicode_casefold/2)  
[unicode_codepoint_valid/1](#unicode_codepoint_valid/1)  
[unicode_confusable/2](#unicode_confusable/2)  
[unicode_confusable/3](#unicode_confusable/3)  
[unicode_identifier_status/2](#unicode_identifier_status/2)  
[unicode_identifier_type/2](#unicode_identifier_type/2)  
[unicode_map/3](#unicode_map/3)  
[unicode_nfc/2](#unicode_nfc/2)  
[unicode_nfd/2](#unicode_nfd/2)  
[unicode_nfkc/2](#unicode_nfkc/2)  
[unicode_nfkc_casefold/2](#unicode_nfkc_casefold/2)  
[unicode_nfkd/2](#unicode_nfkd/2)  
[unicode_property/2](#unicode_property/2)  
[unicode_resolved_scripts/2](#unicode_resolved_scripts/2)  
[unicode_restriction_level/2](#unicode_restriction_level/2)  
[unicode_script/2](#unicode_script/2)  
[unicode_script_extensions/2](#unicode_script_extensions/2)  
[unicode_skeleton/2](#unicode_skeleton/2)  
[unicode_version/1](#unicode_version/1)  
