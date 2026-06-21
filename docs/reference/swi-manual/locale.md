
## 4.23 Localization (locale) support

SWI-Prolog provides (currently limited) support for localized applications.

- The predicates [char_type/2](chartype.html#char_type/2) and [code_type/2](chartype.html#code_type/2) query character classes depending on the locale.
- The predicates [collation_key/2](chartype.html#collation_key/2) and [locale_sort/2](chartype.html#locale_sort/2) can be used for locale dependent sorting of atoms.
- The predicate [format_time/3](system.html#format_time/3) can be used to format time and date representations, where some of the specifiers are locale dependent.
- The predicate [format/2](format.html#format/2) provides locale-specific formatting of numbers. This functionality is based on a more fine-grained localization model that is the subject of this section.

A locale is a (optionally named) read-only object that provides information to locale specific functions.^(121The locale interface described in this section and its effect on [format/2](format.html#format/2) and reading integers from digit groups was discussed on the SWI-Prolog mailinglist. Most input in this discussion is from Ulrich Neumerkel and Richard O'Keefe. The predicates in this section were designed by Jan Wielemaker.) The system creates a default locale object named `default` from the system locale. This locale is used as the initial locale for the three standard streams as well as the `main` thread. Locale sensitive output predicates such as [format/3](format.html#format/3) get their locale from the stream to which they deliver their output. New streams get their locale from the thread that created the stream. Threads get their locale from the thread that created them.

**locale_create**(`-Locale, +Default, +Options`)  
Create a new locale object. `Default` is either an existing locale or a string that denotes the name of a locale provided by the system, such as `"en_EN.UTF-8"`. The values read from the default locale can be modified using `Options`. `Options` provided are:

**alias**(`+Atom`)  
Give the locale a name.

**decimal_point**(`+Atom`)  
Specify the decimal point to use.

**thousands_sep**(`+Atom`)  
Specify the string that delimits digit groups. Only effective is `grouping` is also specified.

**grouping**(`+List`)  
Specify the grouping of digits. Groups are created from the right (least significant) digits, left of the decimal point. `List` is a list of integers, specifying the number of digits in each group, counting from the right. If the last element is `repeat(Count)`, the remaining digits are grouped in groups of size `Count`. If the last element is a normal integer, digits further to the left are not grouped.

For example, the English locale uses

``` code
[ decimal_point('.'), thousands_sep(','), grouping([repeat(3)]) ]
```

Named locales exists until they are destroyed using [locale_destroy/1](locale.html#locale_destroy/1) and they are no longer referenced. Unnamed locales are subject to (atom) garbage collection.

**locale_destroy**(`+Locale`)  
Destroy a locale. If the locale is named, this removes the name association from the locale, after which the locale is left to be reclaimed by garbage collection.

**locale_property**(`?Locale, ?Property`)  
True when `Locale` has `Property`. Properties are the same as the `Options` described with [locale_create/3](locale.html#locale_create/3).

**set_locale**(`+Locale`)  
Set the default locale for the current thread, as well as the locale for the standard streams (`user_input`, `user_output`, `user_error`, `current_output` and `current_input`. This locale is used for new streams, unless overruled using the `locale(Locale)` option of [open/4](IO.html#open/4) or [set_stream/2](IO.html#set_stream/2).

**current_locale**(`-Locale`)  
True when `Locale` is the locale of the calling thread.
