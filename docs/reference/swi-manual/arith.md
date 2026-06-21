
## 4.27 Arithmetic

Arithmetic can be divided into some special purpose integer predicates and a series of general predicates for integer, floating point and rational arithmetic as appropriate. The general arithmetic predicates all handle `expressions`. An expression is either a simple number or a `function`. The arguments of a function are expressions. The functions are described in [section 4.27.2.6](arith.html#sec:4.27.2.6).

### 4.27.1 Special purpose integer arithmetic

The predicates in this section provide more logical operations between integers. They are not covered by the ISO standard, although they are‘part of the community’and found as either library or built-in in many other Prolog systems.

**between**(`+Low, +High, ?Value`)  
`Low` and `High` are integers, `High`` ≥``Low`. If `Value` is an integer, `Low`` ≤``Value`` ≤``High`. When `Value` is a variable it is successively bound to all integers between `Low` and `High`. If `High` is `inf` or `infinite`^(122We prefer `infinite`, but some other Prolog systems already use `inf` for infinity; we accept both for the time being.) [between/3](arith.html#between/3) is true iff `Value`` ≥``Low`, a feature that is particularly interesting for generating integers from a certain value.

**succ**(`?Int1, ?Int2`)  
True if `Int2`` = ``Int1`` + 1` and `Int1`` ≥`. At least one of the arguments must be instantiated to a natural number. This predicate raises the domain error `not_less_than_zero` if called with a negative integer. E.g. `succ(X, 0)` fails silently and `succ(X, -1)` raises a domain error.^(123The behaviour to deal with natural numbers only was defined by Richard O'Keefe to support the common count-down-to-zero in a natural way. Up to 5.1.8, [succ/2](arith.html#succ/2) also accepted negative integers.)

**plus**(`?Int1, ?Int2, ?Int3`)  
True if `Int3`` = ``Int1`` + ``Int2`. At least two of the three arguments must be instantiated to integers.

**divmod**(`+Dividend, +Divisor, -Quotient, -Remainder`)  
This predicate is a shorthand for computing both the `Quotient` and `Remainder` of two integers in a single operation. This allows for exploiting the fact that the low level implementation for computing the quotient also produces the remainder. Timing confirms that this predicate is almost twice as fast as performing the steps independently. Semantically, [divmod/4](arith.html#divmod/4) is defined as below.

``` code
divmod(Dividend, Divisor, Quotient, Remainder) :-
        Quotient  is Dividend div Divisor,
        Remainder is Dividend mod Divisor.
```

Note that this predicate is only available if SWI-Prolog is compiled with unbounded integer support. This is the case for all packaged versions.

**nth_integer_root_and_remainder**(`+N, +I, -Root, -Remainder`)  
True when `Root ** N + Remainder = I`. `N` and `I` must be integers.^(124This predicate was suggested by Markus Triska. The final name and argument order is by Richard O'Keefe. The decision to include the remainder is by Jan Wielemaker. Including the remainder makes this predicate about twice as slow if `Root` is not exact.) `N` must be one or more. If `I` is negative and `N` is *odd*, `Root` and `Remainder` are negative, i.e., the following holds for `I`` < 0`:

``` code
%   I < 0,
%   N mod 2 =\= 0,
    nth_integer_root_and_remainder(
        N, I, Root, Remainder),
    IPos is -I,
    nth_integer_root_and_remainder(
        N, IPos, RootPos, RemainderPos),
    Root =:= -RootPos,
    Remainder =:= -RemainderPos.
```

### 4.27.2 General purpose arithmetic

The general arithmetic predicates are optionally compiled (see [set_prolog_flag/2](flags.html#set_prolog_flag/2) and the **-O** command line option). Compiled arithmetic reduces global stack requirements and improves performance. Unfortunately compiled arithmetic cannot be traced, which is why it is optional.

\[ISO\]`+Expr1` **\>** `+Expr2`  
True if expression `Expr1` evaluates to a larger number than `Expr2`.

\[ISO\]`+Expr1` **\<** `+Expr2`  
True if expression `Expr1` evaluates to a smaller number than `Expr2`.

\[ISO\]`+Expr1` **=\<** `+Expr2`  
True if expression `Expr1` evaluates to a smaller or equal number to `Expr2`.

\[ISO\]`+Expr1` **\>=** `+Expr2`  
True if expression `Expr1` evaluates to a larger or equal number to `Expr2`.

\[ISO\]`+Expr1` **=\\** `+Expr2`  
True if expression `Expr1` evaluates to a number non-equal to `Expr2`.

\[ISO\]`+Expr1` **=:=** `+Expr2`  
True if expression `Expr1` evaluates to a number equal to ` Expr2`.

\[ISO\]`-Number` **is** `+Expr`  
True when `Number` is the value to which `Expr` evaluates. Typically, [is/2](arith.html#is/2) should be used with unbound left operand. If equality is to be tested, [=:=/2](arith.html#=:=/2) should be used. For example:

|  |  |
|----|----|
| `?- 1 is sin(pi/2).` | Fails! sin(pi/2) evaluates to the float 1.0, which does not unify with the integer 1. |
| `?- 1 =:= sin(pi/2).` | Succeeds as expected. |

#### 4.27.2.1 Arithmetic types

SWI-Prolog defines the following numeric types:

- *integer*  
  If SWI-Prolog is built using the *GNU multiple precision arithmetic library* (GMP), integer arithmetic is *unbounded*, which means that the size of integers is limited by available memory only. Without GMP, SWI-Prolog integers are 64-bits, regardless of the native integer size of the platform. The type of integer support can be detected using the Prolog flags [bounded](flags.html#flag:bounded), [min_integer](flags.html#flag:min_integer) and [max_integer](flags.html#flag:max_integer). As the use of GMP is default, most of the following descriptions assume unbounded integer arithmetic.

  Internally, SWI-Prolog has three integer representations. Small integers (defined by the Prolog flag [max_tagged_integer](flags.html#flag:max_tagged_integer)) are encoded directly. Larger integers are represented as 64-bit values on the global stack. Integers that do not fit in 64 bits are represented as serialised GNU MPZ structures on the global stack.

- *rational number*  
  Rational numbers (`Q`) are quotients of two integers (`N/M`). Rational arithmetic is only provided if GMP is used (see above). Rational numbers satisfy the type tests [rational/1](typetest.html#rational/1), [number/1](typetest.html#number/1) and [atomic/1](typetest.html#atomic/1) and may satisfy the type test [integer/1](typetest.html#integer/1), i.e., integers are considered rational numbers. Rational numbers are always kept in *canonical representation*, which means `M` is positive and `N` and `M` have no common divisors. Rational numbers are introduced into the computation using the functions [rational/1](arith.html#f-rational/1), [rationalize/1](arith.html#f-rationalize/1) or the [rdiv/2](arith.html#f-rdiv/2) (rational division) function. If the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals) is `true` (default), division ([//2](arith.html#f-//2)) and integer power ([^/2](arith.html#f-%5E/2)) also produce a rational number.

- *float*  
  Floating point numbers are represented using the C type `double`. On most of today's platforms these are 64-bit IEEE floating point numbers.

Arithmetic functions that require integer arguments accept, in addition to integers, rational numbers with (canonical) denominator‘1’. If the required argument is a float the argument is converted to float. Note that conversion of integers to floating point numbers may raise an overflow exception. In all other cases, arguments are converted to the same type using the order below.

> integer `→` rational number `→` floating point number

#### 4.27.2.2 Rational number examples

The use of rational numbers with unbounded integers allows for exact integer or *fixed point* arithmetic under addition, subtraction, multiplication, division and exponentiation ([^/2](arith.html#f-%5E/2)). Support for rational numbers depends on the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals). If this is `true`, the number division function ([//2](arith.html#f-//2)) and exponentiation function ([^/2](arith.html#f-%5E/2)) generate a rational number on integer and rational arguments and [read/1](termrw.html#read/1) and friends read `[-+][0-9_ ]+/[0-9_ ]+` into a rational number. See also [section 2.15.1.6](syntax.html#sec:2.15.1.6). Here are some examples.

|                             |             |
|-----------------------------|-------------|
| A is 2/6                    | A = 1/3     |
| A is 4/3 + 1                | A = 7/3     |
| A is 4/3 + 1.5              | A = 2.83333 |
| A is 4/3 + rationalize(1.5) | A = 17/6    |

Note that floats cannot represent all decimal numbers exactly. The function [rational/1](arith.html#f-rational/1) creates an *exact* equivalent of the float, while [rationalize/1](arith.html#f-rationalize/1) creates a rational number that is within the float rounding error from the original float. Please check the documentation of these functions for details and examples.

Rational numbers can be printed as decimal numbers with arbitrary precision using the [format/3](format.html#format/3) floating point conversion:

``` code
?- A is 4/3 + rational(1.5),
   format('~50f~n', [A]).
2.83333333333333333333333333333333333333333333333333

A = 17/6
```

#### 4.27.2.3 Rational numbers or floats

SWI-Prolog uses rational number arithmetic if the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals) is `true` and if this is defined for a function on the given operands. This results in perfectly precise answers. Unfortunately rational numbers can get really large and, if a precise answer is not needed, a big waste of memory and CPU time. In such cases one should use floating point arithmetic. The Prolog flag [max_rational_size](flags.html#flag:max_rational_size) provides a *tripwire* to detect cases where rational numbers get big and react on these events.

Floating point arithmetic can be forced by forcing a float into an argument at any point, i.e., the result of a function with at least one float is always float except for the float-to-integer rounding and truncating functions such as [round/1](arith.html#f-round/1), [rational/1](arith.html#f-rational/1) or [float_integer_part/1](arith.html#f-float_integer_part/1).

Float arithmetic is typically forced by using a floating point constant as initial value or operand. Alternatively, the [float/1](arith.html#f-float/1) function forces conversion of the argument.

#### 4.27.2.4 IEEE 754 floating point arithmetic

The Prolog ISO standard defines that floating point arithmetic returns a valid floating point number or raises an exception. IEEE floating point arithmetic defines two modes: raising exceptions and propagating the special float values `NaN`, `Inf`, `-Inf` and `-0.0`. SWI-Prolog implements a part of the [ECLiPSe proposal](http://eclipseclp.org/Specs/core_update_float.html) to support non-exception based processing of floating point numbers. There are four flags that define handling the four exceptional events in floating point arithmetic, providing the choice between `error` and returning the IEEE special value. Note that these flags *only* apply for floating point arithmetic. For example rational division by zero always raises an exception.

|  |  |  |
|----|:--:|:--:|
| **Flag** | **Default** | **Alternative** |
| [float_overflow](flags.html#flag:float_overflow) | error | infinity |
| [float_zero_div](flags.html#flag:float_zero_div) | error | infinity |
| [float_undefined](flags.html#flag:float_undefined) | error | nan |
| [float_underflow](flags.html#flag:float_underflow) | ignore | error |

The Prolog flag [float_rounding](flags.html#flag:float_rounding) and the function [roundtoward/2](arith.html#f-roundtoward/2) control the rounding mode for floating point arithmetic. The default rounding is `to_nearest` and the following alternatives are provided: `to_positive`, `to_negative` and `to_zero`.

\[det\]**float_class**(`+Float, -Class`)  
Wraps C99 **fpclassify()** to access the class of a floating point number. Raises a type error if `Float` is not a float. Defined classes are below.

**nan**  
`Float` is “Not a number” . See [nan/0](arith.html#f-nan/0). May be produced if the Prolog flag [float_undefined](flags.html#flag:float_undefined) is set to `nan`. Although IEEE 754 allows NaN to carry a *payload* and have a sign, SWI-Prolog has only a single NaN values. Note that two NaN *terms* compare equal in the standard order of terms ([==/2](compare.html#==/2), etc.), they compare non-equal for arithmetic ([=:=/2](arith.html#=:=/2), etc.).

**infinite**  
`Float` is positive or negative infinity. See [inf/0](arith.html#f-inf/0). May be produced if the Prolog flag [float_overflow](flags.html#flag:float_overflow) or the flag [float_zero_div](flags.html#flag:float_zero_div) is set to `infinity`.

**zero**  
`Float` is zero (0.0 or -0.0)

**subnormal**  
`Float` is too small to be represented in normalized format. May **not** be produced if the Prolog flag [float_underflow](flags.html#flag:float_underflow) is set to `error`.

**normal**  
`Float` is a normal floating point number.

\[det\]**float_parts**(`+Float, -Mantissa, -Base, -Exponent`)  
True when `Mantissa` is the normalized fraction of `Float`, `Base` is the *radix* and `Exponent` is the exponent. This uses the C function **frexp()**. If `Float` is NaN or `±`Inf `Mantissa` has the same value and `Exponent` is 0 (zero). In the current implementation `Base` is always 2. The following relation is always true:

> `Float =:= Mantissa × Base^Exponent`

\[det\]**bounded_number**(`?Low, ?High, +Num`)  
True if `Low` \< `Num` \< `High`. Raises a type error if `Num` is not a number. This predicate can be used both to check and generate bounds across the various numeric types. Note that a number cannot be bounded by itself and `NaN`, `Inf`, and `-Inf` are not bounded numbers.

If `Low` and/or `High` are variables they will be unified with *tightest* values that still meet the bounds criteria. The generated bounds will be integers if `Num` is an integer; otherwise they will be floats (also see [nexttoward/2](arith.html#f-nexttoward/2) for generating float bounds). Some examples:

``` code
?- bounded_number(0,10,1).
true.

?- bounded_number(0.0,1.0,1r2).
true.

?- bounded_number(L,H,1.0).
L = 0.9999999999999999,
H = 1.0000000000000002.

?- bounded_number(L,H,-1).
L = -2,
H = 0.

?- bounded_number(0,1r2,1).
false.

?- bounded_number(L,H,1.0Inf).
false.
```

#### 4.27.2.5 Floating point arithmetic precision

SWI-Prolog represents floats using the C `double` type. On virtually all modern hardware this implies it uses 64-bit IEEE 754 floating point numbers. See also [section 4.27.2.4](arith.html#sec:4.27.2.4). All floating point arithmetic is performed using C. Different C compilers, different C math libraries and different hardware floating point support may yield different results for the same expression on different instances of SWI-Prolog.

#### 4.27.2.6 Arithmetic Functions

Arithmetic functions are terms which are evaluated by the arithmetic predicates described in [section 4.27.2](arith.html#sec:4.27.2). There are four types of arguments to functions:

|  |  |
|----|----|
| `Expr` | Arbitrary expression, returning either a floating point value or an integer. |
| `IntExpr` | Arbitrary expression that must evaluate to an integer. |
| `RatExpr` | Arbitrary expression that must evaluate to a rational number. |
| `FloatExpr` | Arbitrary expression that must evaluate to a floating point. |

For systems using bounded integer arithmetic (default is unbounded, see [section 4.27.2.1](arith.html#sec:4.27.2.1) for details), integer operations that would cause overflow automatically convert to floating point arithmetic.

SWI-Prolog provides many extensions to the set of floating point functions defined by the ISO standard. The current policy is to provide such functions on‘as-needed’basis if the function is widely supported elsewhere and notably if it is part of the [C99](http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1124.pdf) mathematical library. In addition, we try to maintain compatibility with other Prolog implementations.

\[ISO\]**-** `+Expr`  
`Result`` = -``Expr`

\[ISO\]**+** `+Expr`  
`Result`` = ``Expr`. Note that if `+` is followed by a number, the parser discards the `+`. I.e. `?- integer(+1)` succeeds.

\[ISO\]`+Expr1` **+** `+Expr2`  
`Result`` = ``Expr1`` + ``Expr2`

\[ISO\]`+Expr1` **-** `+Expr2`  
`Result`` = ``Expr1`` - ``Expr2`

\[ISO\]`+Expr1` **\*** `+Expr2`  
`Result`` = ``Expr1`` × ``Expr2`

\[ISO\]`+Expr1` **/** `+Expr2`  
`Result`` = ``Expr1``/``Expr2`. If the flag [iso](flags.html#flag:iso) is `true` or one of the arguments is a float, both arguments are converted to float and the return value is a float. Otherwise the result type depends on the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals). If `true`, the result is always a rational number. If `false` the result is rational if at least one of the arguments is rational. Otherwise (both arguments are integer) the result is integer if the division is exact and float otherwise. See also [section 4.27.2.2](arith.html#sec:4.27.2.2), [///2](arith.html#f-///2), and [rdiv/2](arith.html#f-rdiv/2).

The current default for the Prolog flag [prefer_rationals](flags.html#flag:prefer_rationals) is `false`. Future version may switch this to `true`, providing precise results when possible. The pitfall is that in general rational arithmetic is slower and can become very slow and produce huge numbers that require a lot of (global stack) memory. Code for which the exact results provided by rational numbers is not needed should force float results by making one of the operands float, for example by dividing by `10.0` rather than `10` or by using [float/1](arith.html#f-float/1). Note that when one of the arguments is forced to a float the division is a float operation while if the result is forced to the float the division is done using rational arithmetic.

\[ISO\]`+IntExpr1` **mod** `+IntExpr2`  
Modulo, defined as `Result` = `IntExpr1` - (`IntExpr1` div `IntExpr2`) ` × ` `IntExpr2`, where `div` is *floored* division.

\[ISO\]`+IntExpr1` **rem** `+IntExpr2`  
Remainder of integer division. Behaves as if defined by `Result` is `IntExpr1` - (`IntExpr1` // `IntExpr2`) ` × ` `IntExpr2`

\[ISO\]`+IntExpr1` **//** `+IntExpr2`  
Integer division, defined as `Result` is `rnd_I`(`Expr1`/`Expr2`) . The function `rnd_I` is the default rounding used by the C compiler and available through the Prolog flag [integer_rounding_function](flags.html#flag:integer_rounding_function). In the C99 standard, C-rounding is defined as `towards_zero`.^(125Future versions might guarantee rounding towards zero.)

\[ISO\]**div**(`+IntExpr1, +IntExpr2`)  
Integer division, defined as `Result` is (`IntExpr1` - `IntExpr1` `mod` `IntExpr2`) // `IntExpr2`. In other words, this is integer division that rounds towards -infinity. This function guarantees behaviour that is consistent with [mod/2](arith.html#f-mod/2), i.e., the following holds for every pair of integers `X,Y` where `Y =\= 0`.

``` code
        Q is div(X, Y),
        M is mod(X, Y),
        X =:= Y*Q+M.
```

`+RatExpr` **rdiv** `+RatExpr`  
Rational number division. This function is only available if SWI-Prolog has been compiled with rational number support. See [section 4.27.2.2](arith.html#sec:4.27.2.2) for details.

**gcd**(`+IntExpr1, +IntExpr2`)  
Result is the greatest common divisor of `IntExpr1` and `IntExpr2`. The GCD is always a positive integer. If either expression evaluates to zero the GCD is the result of the other expression.

**lcm**(`+IntExpr1, +IntExpr2`)  
Result is the least common multiple of `IntExpr1`, `IntExpr2`.^(bugIf the system is compiled for bounded integers only [lcm/2](arith.html#f-lcm/2) produces an integer overflow if the product of the two expressions does not fit in a 64 bit signed integer. The default build with unbounded integer support has no such limit.) If either expression evaluates to zero the LCM is zero.

\[ISO\]**abs**(`+Expr`)  
Evaluate `Expr` and return the absolute value of it.

\[ISO\]**sign**(`+Expr`)  
Evaluate to -1 if `Expr`` < 0`, 1 if `Expr`` > 0` and 0 if `Expr`` = 0`. If `Expr` evaluates to a float, the return value is a float (e.g., -1.0, 0.0 or 1.0). In particular, note that sign(-0.0) evaluates to 0.0. See also [copysign/2](arith.html#f-copysign/2).

**cmpr**(`+Expr1, +Expr2`)  
Exactly compares the values `Expr1` and `Expr2` and returns -1 if `Expr1` \< `Expr2`, 0 if they are equal, and 1 if `Expr1` \> `Expr2`. Evaluates to NaN if either or both `Expr1` and `Expr2` are NaN and the Prolog flag [float_undefined](flags.html#flag:float_undefined) is set to `nan`. See also [minr/2](arith.html#f-minr/2) and [maxr/2](arith.html#f-maxr/2).

This function relates to the Prolog numerical comparison predicates [\>/2](arith.html#%3E/2), [=:=/2](arith.html#=:=/2), etc. The Prolog numerical comparison converts the rational in a mixed rational/float comparison to a float, possibly rounding the value. This function converts the float to a rational, comparing the exact values.

\[ISO\]**copysign**(`+Expr1, +Expr2`)  
Evaluate to `X`, where the absolute value of `X` equals the absolute value of `Expr1` and the sign of `X` matches the sign of `Expr2`. This function is based on **copysign()** from C99, which works on double precision floats and deals with handling the sign of special floating point values such as -0.0. Our implementation follows C99 if both arguments are floats. Otherwise, [copysign/2](arith.html#f-copysign/2) evaluates to `Expr1` if the sign of both expressions matches or -`Expr1` if the signs do not match. Here, we use the extended notion of signs for floating point numbers, where the sign of -0.0 and other special floats is negative.

**nexttoward**(`+Expr1, +Expr2`)  
Evaluates to floating point number following `Expr1` in the direction of `Expr2`. This relates to [epsilon/0](arith.html#f-epsilon/0) in the following way:

``` code
?- epsilon =:= nexttoward(1,2)-1.
true.
```

**roundtoward**(`+Expr1, +RoundMode`)  
Evaluate `Expr1` using the floating point rounding mode `RoundMode`. This provides a local alternative to the Prolog flag [float_rounding](flags.html#flag:float_rounding). This function can be nested. The supported values for `RoundMode` are the same as the flag values: `to_nearest`, `to_positive`, `to_negative` or `to_zero`.

Note that floating point arithmetic is provided by the C compiler and C runtime library. Unfortunately most C libraries do not correctly implement the rounding modes for notably the trigonometry and exponential functions. There exist correct libraries such as [crlibm](https://github.com/taschini/crlibm), but these libraries are large, most of them are poorly maintained or have an incompatible license. C runtime libraries do a better job using the default *to nearest* rounding mode. SWI-Prolog now assumes this mode is correct and translates upward rounding to be the [nexttoward/2](arith.html#f-nexttoward/2) infinity and downward rounding [nexttoward/2](arith.html#f-nexttoward/2) -infinity. If the “to nearest” rounding mode is correct, this ensures that the true value is between the downward and upward rounded values, although the generated interval is larger than needed. Unfortunately this is not the case as shown in [Accuracy of Mathematical Functions in Single, Double, Extended Double and Quadruple Precision](https://hal.inria.fr/hal-03141101) by *Vincenzo Innocente and Paul Zimmermann*.

\[ISO\]**max**(`+Expr1, +Expr2`)  
Evaluate to the larger of `Expr1` and `Expr2`. Both arguments are compared after converting to the same type, but the return value is in the original type. For example, max(2.5, 3) compares the two values after converting to float, but returns the integer 3. If both values are numerical equal the returned max is of the type used for the comparison. For example, the max of 1 and 1.0 is 1.0 because both numbers are converted to float for the comparison. However, the special float -0.0 is smaller than 0.0 as well as the integer 0. If the Prolog flag [float_undefined](flags.html#flag:float_undefined) is set to `nan` and one of the arguments evaluates to NaN, the result is NaN.

The function [maxr/2](arith.html#f-maxr/2) is similar, but uses exact (rational) comparison if `Expr1` and `Expr2` have a different type, propagate the rational (integer) rather and the float if the two compare equal and propagate the non-NaN value in case one is NaN.

**maxr**(`+Expr1, +Expr2`)  
Evaluate to the larger of `Expr1` and `Expr2` using exact comparison (see [cmpr/2](arith.html#f-cmpr/2)). If the two values are exactly equal, and one of the values is rational, the result will be that value; the objective being to avoid "pollution" of any precise calculation with a potentially imprecise float. So `max(1,1.0)` evaluates to 1.0 while `maxr(1,1.0)` evaluates to 1. This also means that 0 is preferred over 0.0 or -0.0; -0.0 is still considered smaller than 0.0.

[maxr/2](arith.html#f-maxr/2) also treats NaN's as missing values so `maxr(1,nan)` evaluates to 1.

\[ISO\]**min**(`+Expr1, +Expr2`)  
Evaluate to the smaller of `Expr1` and `Expr2`. See [max/2](arith.html#f-max/2) for a description of type handling.

**minr**(`+Expr1, +Expr2`)  
Evaluate to the smaller of `Expr1` and `Expr2` using exact comparison (see [cmpr/2](arith.html#f-cmpr/2)). See [maxr/2](arith.html#f-maxr/2) for a description of type handling.

\[deprecated\]**.**(`+Char,[]`)  
A list of one element evaluates to the character code of this element.^(126The function is documented as `.``/2`. Using SWI-Prolog v7 and later the actual functor is `[|]``/2`.) This implies `"a"` evaluates to the character code of the letter‘a’(97) using the traditional mapping of double quoted string to a list of character codes. `Char` is either a valid code point (non-negative integer up to the Prolog flag [max_char_code](flags.html#flag:max_char_code)) or a one-character atom. Arithmetic evaluation also translates a string object (see [section 5.2](string.html#sec:5.2)) of one character length into the character code for that character. This implies that expression `"a"` works if the Prolog flag [double_quotes](flags.html#flag:double_quotes) is set to one of `codes`, `chars` or `string`.

Getting access to character codes this way originates from DEC10 Prolog. ISO has the `0'` syntax and the predicate [char_code/2](manipatom.html#char_code/2). Future versions may drop support for `X is "a"`.

**random**(`+IntExpr`)  
Evaluate to a random integer `i` for which `0 ≤i < ``IntExpr`. The system has two implementations. If it is compiled with support for unbounded arithmetic (default) it uses the GMP library random functions. In this case, each thread keeps its own random state. The default algorithm is the *Mersenne Twister* algorithm. The seed is set when the first random number in a thread is generated. If available, it is set from `/dev/random`.^(127On Windows the state is initialised from **CryptGenRandom()**.) Otherwise it is set from the system clock. If unbounded arithmetic is not supported, random numbers are shared between threads and the seed is initialised from the clock when SWI-Prolog was started. The predicate [set_random/1](miscarith.html#set_random/1) can be used to control the random number generator.

**Warning!** Although properly seeded (if supported on the OS), the Mersenne Twister algorithm does *not* produce cryptographically secure random numbers. To generate cryptographically secure random numbers, use crypto_n_random_bytes/2 from library `library(crypto)` provided by the `ssl` package.

**random_float**  
Evaluate to a random float `I` for which `0.0 < i < 1.0`. This function shares the random state with [random/1](arith.html#f-random/1). All remarks with the function [random/1](arith.html#f-random/1) also apply for [random_float/0](arith.html#f-random_float/0). Note that both sides of the domain are *open*. This avoids evaluation errors on, e.g., [log/1](arith.html#f-log/1) or [//2](arith.html#f-//2) while no practical application can expect 0.0.^(128Richard O'Keefe said: “If you *are* generating IEEE doubles with the claimed uniformity, then 0 has a 1 in `2^53 = 1 in 9,007,199,254,740,992` chance of turning up. No program that expects \[0.0,1.0) is going to be surprised when 0.0 fails to turn up in a few millions of millions of trials, now is it? But a program that expects (0.0,1.0) could be devastated if 0.0 did turn up.’)

\[ISO\]**round**(`+Expr`)  
Evaluate `Expr` and round the result to the nearest integer. According to ISO, [round/1](arith.html#f-round/1) is defined as `floor(Expr+1/2)`, i.e., rounding *down*. This is an unconventional choice under which the relation `round(Expr) == -round(-Expr)` does not hold. SWI-Prolog rounds *outward*, e.g., `round(1.5) =:= 2` and `round(-1.5) =:= -2`.

**integer**(`+Expr`)  
Same as [round/1](arith.html#f-round/1) (backward compatibility).

\[ISO\]**float**(`+Expr`)  
Translate the result to a floating point number. Normally, Prolog will use integers whenever possible. When used around the 2nd argument of [is/2](arith.html#is/2), the result will be returned as a floating point number. In other contexts, the operation has no effect.

**rational**(`+Expr`)  
Convert the `Expr` to a rational number or integer. The function returns the input on integers and rational numbers. For floating point numbers, the returned rational number *exactly* represents the float. As floats cannot exactly represent all decimal numbers the results may be surprising. In the examples below, doubles can represent 0.25 and the result is as expected, in contrast to the result of `rational(0.1)`. The function [rationalize/1](arith.html#f-rationalize/1) remedies this. See [section 4.27.2.2](arith.html#sec:4.27.2.2) for more information on rational number support.

``` code
?- A is rational(0.25).

A is 1r4
?- A is rational(0.1).
A = 3602879701896397r36028797018963968
```

For every *normal* float `X` the relation `X` `=:=` rational(`X`) holds.

This function raises an `evaluation_error(undefined)` if `Expr` is NaN and `evaluation_error(rational_overflow)` if `Expr` is Inf.

**rationalize**(`+Expr`)  
Convert the `Expr` to a rational number or integer. The function is similar to [rational/1](arith.html#f-rational/1), but the result is only accurate within the rounding error of floating point numbers, generally producing a much smaller denominator.^(129The names [rational/1](arith.html#f-rational/1) and [rationalize/1](arith.html#f-rationalize/1) as well as their semantics are inspired by Common Lisp.130The implementation of rationalize as well as converting a rational number into a float is copied from ECLiPSe and covered by the *Cisco-style Mozilla Public License Version 1.1*.)

``` code
?- A is rationalize(0.25).

A = 1r4
?- A is rationalize(0.1).

A = 1r10
```

For every *normal* float `X` the relation `X` `=:=` rationalize(`X`) holds.

This function raises the same exceptions as [rational/1](arith.html#f-rational/1) on non-normal floating point numbers.

**numerator**(`+RationalExpr`)  
If `RationalExpr` evaluates to a rational number or integer, evaluate to the top/left value. Evaluates to itself if `RationalExpr` evaluates to an integer. See also [denominator/1](arith.html#f-denominator/1). The following is true for any rational `X`.

``` code
X =:= numerator(X)/denominator(X).
```

**denominator**(`+RationalExpr`)  
If `RationalExpr` evaluates to a rational number or integer, evaluate to the bottom/right value. Evaluates to 1 (one) if `RationalExpr` evaluates to an integer. See also [numerator/1](arith.html#f-numerator/1). The following is true for any rational `X`.

``` code
X =:= numerator(X)/denominator(X).
```

\[ISO\]**float_fractional_part**(`+Expr`)  
Fractional part of a floating point number. Negative if `Expr` is negative, rational if `Expr` is rational and 0 if `Expr` is integer. The following relation is always true: `X is float_fractional_part(X) + float_integer_part(X)`.

\[ISO\]**float_integer_part**(`+Expr`)  
Integer part of floating point number. Negative if `Expr` is negative, `Expr` if `Expr` is integer.

\[ISO\]**truncate**(`+Expr`)  
Truncate `Expr` to an integer. If `Expr`` ≥` this is the same as `floor(Expr)`. For `Expr`` < 0` this is the same as `ceil(Expr)`. That is, [truncate/1](arith.html#f-truncate/1) rounds towards zero.

\[ISO\]**floor**(`+Expr`)  
Evaluate `Expr` and return the largest integer smaller or equal to the result of the evaluation.

\[ISO\]**ceiling**(`+Expr`)  
Evaluate `Expr` and return the smallest integer larger or equal to the result of the evaluation.

**ceil**(`+Expr`)  
Same as [ceiling/1](arith.html#f-ceiling/1) (backward compatibility).

\[ISO\]`+IntExpr1` **\>\>** `+IntExpr2`  
Bitwise shift `IntExpr1` by `IntExpr2` bits to the right. The ISO standard dictates shifting a negative value is *implementation defined*. SWI-Prolog defines shifting negative integers to be defined as `-(-Int>>Shift)`. Shifting positive integers by more than their size results in 0 (zero). Shifting negative integers by more then their size results in -1. I.e., `A is -3464 >> 100` binds `A` to -1. If `IntExpr2` is negative, a right shift (see [\>\>/2](arith.html#f-%3E%3E/2)) is performed with the negated value of `IntExpr2`.

\[ISO\]`+IntExpr1` **\<\<** `+IntExpr2`  
Bitwise shift `IntExpr1` by `IntExpr2` bits to the left. The ISO standard dictates shifting a negative value is *implementation defined*. SWI-Prolog defines shifting negative integers to be defined as `-(-Int<<Shift)`. If `IntExpr2` is negative, a left shift (see [\<\</2](arith.html#f-%3C%3C/2)) is performed with the negated value of `IntExpr2`.

\[ISO\]`+IntExpr1` **\\** `+IntExpr2`  
Bitwise‘or’ `IntExpr1` and `IntExpr2`.

\[ISO\]`+IntExpr1` **/\\** `+IntExpr2`  
Bitwise‘and’ `IntExpr1` and `IntExpr2`.

\[ISO\]`+IntExpr1` **xor** `+IntExpr2`  
Bitwise‘exclusive or’ `IntExpr1` and `IntExpr2`.

\[ISO\]**\\** `+IntExpr`  
Bitwise negation. The returned value is the one's complement of `IntExpr`.

\[ISO\]**sqrt**(`+Expr`)  
`Result`` = &#8730;(``Expr``)`.

\[ISO\]**sin**(`+Expr`)  
`Result`` = sin(``Expr``)`. `Expr` is the angle in radians.

\[ISO\]**cos**(`+Expr`)  
`Result`` = cos(``Expr``)`. `Expr` is the angle in radians.

\[ISO\]**tan**(`+Expr`)  
`Result`` = tan(``Expr``)`. `Expr` is the angle in radians.

\[ISO\]**asin**(`+Expr`)  
`Result`` = arcsin(``Expr``)`. `Result` is the angle in radians.

\[ISO\]**acos**(`+Expr`)  
`Result`` = arccos(``Expr``)`. `Result` is the angle in radians.

\[ISO\]**atan**(`+Expr`)  
`Result`` = arctan(``Expr``)`. `Result` is the angle in radians.

\[ISO\]**atan2**(`+YExpr, +XExpr`)  
`Result`` = arctan(``YExpr``/``XExpr``)`. `Result` is the angle in radians. The return value is in the range `[-π...π`. Used to convert between rectangular and polar coordinate system.

Note that the ISO Prolog standard demands `atan2(0.0,0.0)` to raise an evaluation error, whereas the C99 and POSIX standards demand this to evaluate to 0.0. SWI-Prolog follows C99 and POSIX.

**atan**(`+YExpr, +XExpr`)  
Same as [atan2/2](arith.html#f-atan2/2) (backward compatibility).

**sinh**(`+Expr`)  
`Result`` = sinh(``Expr``)`. The hyperbolic sine of `X` is defined as `e ** X - e ** -X / 2`.

**cosh**(`+Expr`)  
`Result`` = cosh(``Expr``)`. The hyperbolic cosine of `X` is defined as `e ** X + e ** -X / 2`.

**tanh**(`+Expr`)  
`Result`` = tanh(``Expr``)`. The hyperbolic tangent of `X` is defined as `sinh( X ) / cosh( X )`.

**asinh**(`+Expr`)  
`Result`` = arcsinh(``Expr``)` (inverse hyperbolic sine).

**acosh**(`+Expr`)  
`Result`` = arccosh(``Expr``)` (inverse hyperbolic cosine).

**atanh**(`+Expr`)  
`Result`` = arctanh(``Expr``)`. (inverse hyperbolic tangent).

\[ISO\]**log**(`+Expr`)  
Natural logarithm. `Result`` = ln(``Expr``)`

**log10**(`+Expr`)  
Base-10 logarithm. `Result`` = log10(``Expr``)`

\[ISO\]**exp**(`+Expr`)  
`Result`` = e **``Expr`

\[ISO\]`+Expr1` **\*\*** `+Expr2`  
`Result`` = ``Expr1``**``Expr2`. The result is a float, unless SWI-Prolog is compiled with unbounded integer support and the inputs are integers and produce an integer result. The integer expressions `0 ** I`, `1 ** I` and `-1 ** I` are guaranteed to work for any integer `I`. Other integer base values generate a `resource` error if the result does not fit in memory.

The ISO standard demands a float result for all inputs and introduces [^/2](arith.html#f-%5E/2) for integer exponentiation. The function [float/1](arith.html#f-float/1) can be used on one or both arguments to force a floating point result. Note that casting the *input* result in a floating point computation, while casting the *output* performs integer exponentiation followed by a conversion to float.

\[ISO\]`+Expr1` **^** `+Expr2`  
In SWI-Prolog, [^/2](arith.html#f-%5E/2) is equivalent to [\*\*/2](arith.html#f-**/2). The ISO version is similar, except that it produces a evaluation error if both `Expr1` and `Expr2` are integers and the result is not an integer. The table below illustrates the behaviour of the exponentiation functions in ISO and SWI. Note that if the exponent is negative the behavior of `Int``^``Int` depends on the flag [prefer_rationals](flags.html#flag:prefer_rationals), producing either a rational number or a floating point number.

|          |         |                             |                 |              |
|----------|---------|-----------------------------|-----------------|--------------|
| `Expr1`  | `Expr2` | Function                    | SWI             | ISO          |
| Int      | Int     | [\*\*/2](arith.html#f-**/2) | Int or Rational | Float        |
| Int      | Float   | [\*\*/2](arith.html#f-**/2) | Float           | Float        |
| Rational | Int     | [\*\*/2](arith.html#f-**/2) | Rational        | \-           |
| Float    | Int     | [\*\*/2](arith.html#f-**/2) | Float           | Float        |
| Float    | Float   | [\*\*/2](arith.html#f-**/2) | Float           | Float        |
| Int      | Int     | [^/2](arith.html#f-%5E/2)   | Int or Rational | Int or error |
| Int      | Float   | [^/2](arith.html#f-%5E/2)   | Float           | Float        |
| Rational | Int     | [^/2](arith.html#f-%5E/2)   | Rational        | \-           |
| Float    | Int     | [^/2](arith.html#f-%5E/2)   | Float           | Float        |
| Float    | Float   | [^/2](arith.html#f-%5E/2)   | Float           | Float        |

**powm**(`+IntExprBase, +IntExprExp, +IntExprMod`)  
`Result`` = (``IntExprBase``**``IntExprExp``) modulo ``IntExprMod`. Only available when compiled with unbounded integer support. This formula is required for Diffie-Hellman key-exchange, a technique where two parties can establish a secret key over a public network. `IntExprBase` and `IntExprExp` must be non-negative (`>=0`), `IntExprMod` must be positive (`>0`).^(131The underlying GMP **mpz_powm()** function allows negative values under some conditions. As the conditions are expensive to pre-compute, error handling from GMP is non-trivial and negative values are not needed for Diffie-Hellman key-exchange we do not support these.)

**lgamma**(`+Expr`)  
Return the natural logarithm of the absolute value of the Gamma function.^(132Some interfaces also provide the sign of the Gamma function. We cannot do that in an arithmetic function. Future versions may provide a *predicate* lgamma/3 that returns both the value and the sign.)

**erf**(`+Expr`)  
[Wikipedia](https://en.wikipedia.org/wiki/Error_function): “In mathematics, the error function (also called the Gauss error function) is a special function (non-elementary) of sigmoid shape which occurs in probability, statistics and partial differential equations.”

**erfc**(`+Expr`)  
[Wikipedia](https://en.wikipedia.org/wiki/Error_function): “The complementary error function.”

\[ISO\]**pi**  
Evaluate to the mathematical constant `π` (3.14159 ... ).

**e**  
Evaluate to the mathematical constant `e` (2.71828 ... ).

**epsilon**  
Evaluate to the difference between the float 1.0 and the first larger floating point number. Deprecated. The function [nexttoward/2](arith.html#f-nexttoward/2) provides a better alternative.

**inf**  
Evaluate to positive infinity. See [section 2.15.1.7](syntax.html#sec:2.15.1.7) and [section 4.27.2.4](arith.html#sec:4.27.2.4). This value can be negated using [-/1](arith.html#f--/1).

**nan**  
Evaluate to *Not a Number*. See [section 2.15.1.7](syntax.html#sec:2.15.1.7) and [section 4.27.2.4](arith.html#sec:4.27.2.4).

**cputime**  
Evaluate to a floating point number expressing the CPU time (in seconds) used by Prolog up till now. See also [statistics/2](builtin-statistics.html#statistics/2) and [time/1](statistics.html#time/1).

**eval**(`+Expr`)  
Evaluate `Expr`. Although ISO standard dictates that‘`A`=1+2, `B` is `A`’works and unifies `B` to 3, it is widely felt that source level variables in arithmetic expressions should have been limited to numbers. In this view the eval function can be used to evaluate arbitrary expressions.^(133The [eval/1](arith.html#f-eval/1) function was first introduced by ECLiPSe and is under consideration for YAP.)

**Bitvector functions**

The functions below are not covered by the standard. The [msb/1](arith.html#f-msb/1) function also appears in hProlog and SICStus Prolog. The [getbit/2](arith.html#f-getbit/2) function also appears in ECLiPSe, which also provides `setbit(Vector,Index)` and `clrbit(Vector,Index)`. The others are SWI-Prolog extensions that improve handling of ---unbounded--- integers as bit-vectors.

**msb**(`+IntExpr`)  
Return the largest integer `N` such that `(IntExpr >> N) /\ 1 =:= 1`. This is the (zero-origin) index of the most significant 1 bit in the value of `IntExpr`, which must evaluate to a positive integer. Errors for 0, negative integers, and non-integers.

**lsb**(`+IntExpr`)  
Return the smallest integer `N` such that `(IntExpr >> N) /\ 1 =:= 1`. This is the (zero-origin) index of the least significant 1 bit in the value of `IntExpr`, which must evaluate to a positive integer. Errors for 0, negative integers, and non-integers.

**popcount**(`+IntExpr`)  
Return the number of 1s in the binary representation of the non-negative integer `IntExpr`.

**getbit**(`+IntExprV, +IntExprI`)  
Evaluates to the bit value (0 or 1) of the `IntExprI`-th bit of `IntExprV`. Both arguments must evaluate to non-negative integers. The result is equivalent to `(IntExprV >> IntExprI)/\1`, but more efficient because materialization of the shifted value is avoided. Future versions will optimise `(IntExprV >> IntExprI)/\1` to a call to [getbit/2](arith.html#f-getbit/2), providing both portability and performance.^(134This issue was fiercely debated at the ISO standard mailinglist. The name *getbit* was selected for compatibility with ECLiPSe, the only system providing this support. Richard O'Keefe disliked the name and argued that efficient handling of the above implementation is the best choice for this functionality.)
