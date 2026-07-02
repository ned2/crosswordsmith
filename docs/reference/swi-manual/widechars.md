
## 2.18 Wide character support

SWI-Prolog represents characters using [Unicode](https://home.unicode.org). Unicode defines *code points* in the range `0 ... 0x10FFFF`. These code points represent virtually any character in any language. In addition, the Unicode standard defines character classes (letter, digit, punctuation, etc.), case conversion and much more. Unicode is a super set of ISO 8859-1 (ISO Latin-1), which is a superset of US-ASCII.

SWI-Prolog has two representations for atoms and string objects (see [section 5.2](string.html#sec:5.2)). If the text fits in ISO Latin-1, it is represented as an array of 8-bit characters. Otherwise the text is represented as an array of `wchar_t` characters. On virtually all systems except for MS-Windows, `wchar_t` is a 32-bit unsigned integer and thus capable of representing all Unicode code points. On MS-Windows `wchar_t` is a 16-bit unsigned integer and thus only capable of representing the code points `0 ... 0xFFFF`. As of SWI-Prolog version 8.5.14, the `wchar_t` is (on Windows) interpreted as a UTF-16 string. The UTF-16 encoding uses *surrogate pairs* to represent the code points `0x10000 ... 0x10FFFF` as two *code units* in the `0xD800 ... 0xDFFF`. As Unicode *code points*, this range is *unassigned*. For consistency, SWI-Prolog does not accept integers in the surrogate pair range as valid code points, e.g.

``` code
?- char_code(X, 0xD800).
ERROR: Type error: `character_code' expected, found `55296' (an integer)
```

The internal character representation is completely transparent to the Prolog user. Users of the foreign language interface as described in [chapter 12](foreign.html#sec:12) sometimes need to be aware of these issues though.

Character coding comes into view when characters of strings need to be read from or written to file or when they have to be communicated to other software components using the foreign language interface. In this section we only deal with I/O through streams, which includes file I/O as well as I/O through network sockets.

### 2.18.1 Wide character encodings on streams

Although characters are uniquely coded using the Unicode standard internally, streams and files are byte (8-bit) oriented and there are a variety of ways to represent the larger Unicode codes in an 8-bit octet stream. The most popular one, especially in the context of the web, is UTF-8. Bytes 0 ... 127 represent simply the corresponding US-ASCII character, while bytes 128 ... 255 are used for multi-byte encoding of characters placed higher in the Unicode space. Especially on MS-Windows the 16-bit UTF-16 standard, represented by pairs of bytes, is also popular.

Prolog I/O streams have a property called *encoding* which specifies the used encoding that influences [get_code/2](chario.html#get_code/2) and [put_code/2](chario.html#put_code/2) as well as all the other text I/O predicates.

The default encoding for files is derived from the Prolog flag [encoding](flags.html#flag:encoding), which is initialised from `setlocale(LC_CTYPE, NULL)` to one of `text`, `utf8` or `iso_latin_1`. One of the latter two is used if the encoding name is recognized, while `text` is used as default. Using `text`, the translation is left to the wide-character functions of the C library.^(42The Prolog native UTF-8 mode is considerably faster than the generic **mbrtowc()** one.) On MS-Windows the default is unconditionally `utf8`, irrespective of the system code page, as UTF-8 is the de-facto encoding for source files and the Windows C runtime's locale-based wide-character functions provide weaker Unicode coverage than Prolog's own tables. The encoding can be specified explicitly in [load_files/2](consulting.html#load_files/2) for loading Prolog source with an alternative encoding, [open/4](IO.html#open/4) when opening files or using [set_stream/2](IO.html#set_stream/2) on any open stream. For Prolog source files we also provide the [encoding/1](consulting.html#encoding/1) directive that can be used to switch between encodings that are compatible with US-ASCII (`ascii`, `iso_latin_1`, `utf8` and many locales). See also [section 3.1.3](projectfiles.html#sec:3.1.3) for writing Prolog files with non-US-ASCII characters and [section 2.15.1.9](syntax.html#sec:2.15.1.9) for syntax issues. For additional information and Unicode resources, please visit [http://www.unicode.org/](http://www.unicode.org/).

SWI-Prolog currently defines and supports the following encodings:

**octet**  
Default encoding for `binary` streams. This causes the stream to be read and written fully untranslated.

**ascii**  
7-bit encoding in 8-bit bytes. Equivalent to `iso_latin_1`, but generates errors and warnings on encountering values above 127.

**iso_latin_1**  
8-bit encoding supporting many Western languages. This causes the stream to be read and written fully untranslated. The above is the SWI-Prolog native name. This encoding may be specified using the official [IANA](https://www.iana.org) name `ISO-8859-1`.

**text**  
C library default locale encoding for text files. Files are read and written using the C library functions **mbrtowc()** and **wcrtomb()**. This may be the same as one of the other encodings, notably it may be the same as `iso_latin_1` for Western languages and `utf8` in a UTF-8 context.

**utf8**  
Multi-byte encoding of full Unicode, compatible with `ascii`. See above. The above is the SWI-Prolog native name. This encoding may be specified using the official [IANA](https://www.iana.org) name `UTF-8`.

**utf16be**  
**utf16le**  
UTF-16 encoding. Reads input in pairs of bytes. `utf16be` is *Big Endian*, putting the most significant byte first and `utf16le` is *Little Endian*, putting the most significant byte second. UTF-16 can represent full Unicode using *surrogate pairs*. The above are the SWI-Prolog native names. These encodings may be specified using the official [IANA](https://www.iana.org) names `UTF-16BE` and `UTF-16LE`. For backward compatibility we also support `unicode_be` and `unicode_le`.

Note that not all encodings can represent all characters. This implies that writing text to a stream may cause errors because the stream cannot represent these characters. The behaviour of a stream on these errors can be controlled using [set_stream/2](IO.html#set_stream/2). Initially the terminal stream writes the characters using Prolog escape sequences while other streams generate an I/O exception.

#### 2.18.1.1 BOM: Byte Order Mark

From [section 2.18.1](widechars.html#sec:2.18.1), you may have got the impression that text files are complicated. This section deals with a related topic, making life often easier for the user, but providing another worry to the programmer. **BOM** or *Byte Order Marker* is a technique for identifying Unicode text files as well as the encoding they use. Such files start with the Unicode character 0xFEFF, a non-breaking, zero-width space character. This is a pretty unique sequence that is not likely to be the start of a non-Unicode file and uniquely distinguishes the various Unicode file formats. As it is a zero-width blank, it even doesn't produce any output. This solves all problems, or ... Some formats start off as US-ASCII and may contain some encoding mark to switch to UTF-8, such as the `encoding="UTF-8"` in an XML header. Such formats often explicitly forbid the use of a UTF-8 BOM. In other cases there is additional information revealing the encoding, making the use of a BOM redundant or even illegal.

The BOM is handled by SWI-Prolog [open/4](IO.html#open/4) predicate. By default, text files are probed for the BOM when opened for reading. If a BOM is found, the encoding is set accordingly and the property `bom(true)` is available through [stream_property/2](IO.html#stream_property/2). When opening a file for writing, writing a BOM can be requested using the option `bom(true)` with [open/4](IO.html#open/4).
