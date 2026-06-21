sweep: SWI-Prolog Embedded in Emacs

Eshel Yaron[me@eshelyaron.com](mailto:me@eshelyaron.com)

Abstract

Sweep is an embedding of SWI-Prolog in GNU Emacs. It provides an interface for executing Prolog queries and consuming their results from Emacs Lisp. Sweep further builds on top of this interface and on top of the standard Emacs facilities to provide advanced features for developing SWI-Prolog programs in Emacs.

# Table of Contents

[1 Installation](#sec:1)

[2 Getting started](#sec:2)

## 1 Installation

Installing Sweep requires:

- Emacs 27 or later, and
- SWI-Prolog 8.5.18 or later.

Sweep is available from NonGNU ELPA, to install it simply type in Emacs `M-x package-install RET sweeprolog RET`.

Note that in Emacs prior to version 28, you need to explicitly enable NonGNU ELPA by adding something like the following to your Emacs configuration:

``` code
  (with-eval-after-load 'package
    (add-to-list 'package-archives '("nongnu" . "https://elpa.nongnu.org/nongnu/")))
```

## 2 Getting started

After installing the `sweeprolog` Elisp library, load it into Emacs:

``` code
  (require 'sweeprolog)
```

All set! You can now use Sweep for Prolog development and for integrating Prolog into your Emacs Lisp code. For a full description of the different features of Sweep, see [the Sweep manual](https://eshelyaron.com/sweep.html).
