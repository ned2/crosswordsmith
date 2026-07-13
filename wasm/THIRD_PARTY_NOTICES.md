# Third-party notices — crosswordsmith browser bundle

This file is the redistribution notice for the built browser bundle in
`wasm/client/` — `swipl-web.js`, `swipl-web.wasm`, `swipl-web.data`,
`crosswordsmith.qlf`, `worker.js`, and their content-hashed copies. Ship this
file alongside those artifacts (the build stamps it into `wasm/client/`, and
`build-manifest.json`'s `licenses` field names it); make it reachable from any
site that serves them.

The bundle is compiled from the following pinned sources (also recorded
per-build in `build-manifest.json`):

| Pinned input | Version / commit |
|---|---|
| SWI-Prolog (`swipl-devel` + its package submodules) | commit `aa6289399` (V10.1.10-17) |
| zlib | 1.3.2 (sha256 `bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16`) |
| Emscripten (emsdk) | 6.0.1 |
| crosswordsmith | the repo this file lives in (`buildId` in `build-manifest.json`) |

Every licence text below was copied verbatim from those pinned sources on
2026-07-13. The inventory was derived from the build's actual link line and
object list (not from documentation): the engine built with the
crosswordsmith-web package profile (payload plan Phase 3 — packages
clib+http+json+sgml configured; foreign extensions uri, readutil and json are
the ONLY package C code linked, `build-manifest.json`'s
`profile.staticExtensions`), the C code the shipped objects vendor, the one
cross-compiled C library (zlib), and the Emscripten runtime components
embedded in the emitted `.js`/`.wasm`.

Components inspected and ruled OUT (present in the source tree but not in the
shipped artifacts): every plugin archive outside the static-extension
keep-list (sgml2pl, http_stream, sha4pl, hashstream, md54pl, crypt, memfile,
files, prolog_stream — dropped from the swipl-web link by build-wasm.sh step
2.6, so the MD5/SHA/Gladman notices those would carry do not apply) and the
Info-ZIP decryption code inside minizip (`NOUNCRYPT` is force-defined, so it
is compiled out). Earlier full-image builds also linked PCRE2 and the nlp/
semweb packages' vendored code (Double Metaphone, Porter, Snowball, i-Sub);
those packages are no longer configured, so their notices were removed —
restore them from git history if the package profile ever grows back.

## Summary

| # | Component | What of it ships | Licence |
|---|---|---|---|
| 1 | SWI-Prolog | the engine (`swipl-web.*`) built with the crosswordsmith-web profile: foreign extensions uri/readutil/json, boot + the 22-file preload library in `.data` | BSD-2-Clause |
| 2.1 | LibBF (Fabrice Bellard) | arbitrary-precision arithmetic, compiled into the engine (`USE_GMP=OFF`) | MIT |
| 2.2 | MT19937 (Christian Stigen Larsen) | PRNG used by LibBF | modified BSD |
| 2.3 | minizip (Gilles Vollant et al.) | zip archive reading, compiled into the engine | zlib |
| 2.4 | libtai (D. J. Bernstein) | TAI/UTC time handling, compiled into the engine | public domain |
| 2.5 | dtoa (David M. Gay / Lucent) | float↔string conversion, compiled into the engine | permissive notice |
| 2.6 | MurmurHash (Austin Appleby) | hashing in the engine | public domain |
| 3 | zlib | statically linked into the engine | zlib |
| 4.1 | Emscripten | the JS glue in `swipl-web.js`, runtime support in `.wasm` | MIT / NCSA (dual) |
| 4.2 | musl libc | libc compiled into `.wasm` | MIT |
| 4.3 | LLVM compiler-rt builtins | compiler runtime compiled into `.wasm` | NCSA / MIT (dual, MIT elected) |
| 4.4 | dlmalloc (Doug Lea) | the allocator compiled into `.wasm` | public domain (CC0) |
| 4.5 | MiniLZ4 (Pierre Curto) | `.data` decompression code embedded in `swipl-web.js` (`-sLZ4=1`) | MIT |
| 5 | crosswordsmith | `crosswordsmith.qlf` (the qcompiled app), `worker.js`, the SDK facade | MIT |

---

## 1. SWI-Prolog — BSD-2-Clause

Copyright (c) 1985–2026, University of Amsterdam, VU University Amsterdam,
SWI-Prolog Solutions b.v., and contributors. Verbatim from `LICENSE` at the
pinned commit:

````text
# SWI-Prolog copyright status

SWI-Prolog is covered by the Simplified BSD license:

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in
       the documentation and/or other materials provided with the
       distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
    ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

## Possible additional license requirements

Note that SWI-Prolog may be linked with libraries covered by more
restrictive license and therefore the above conditions may not suffice
for a particular version of SWI-Prolog or a program loaded into
SWI-Prolog. Run the following predicate to find which components with
additional requirements are loaded into a particular version.

  ```
  ?- license.
  ```

````

The additional-requirement components that predicate refers to are enumerated
in sections 2.1–2.6 below for exactly the set compiled into this bundle.

## 2. Components vendored inside SWI-Prolog

### 2.1 LibBF — MIT

Arbitrary-precision arithmetic (the wasm build uses `USE_GMP=OFF`, so LibBF —
not GMP — ships). Verbatim from `src/libbf/libbf.h`:

````text
/*
 * Tiny arbitrary precision floating point library
 *
 * Copyright (c) 2017-2020 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
````

### 2.2 Mersenne Twister (MT19937) — modified BSD

PRNG compiled with LibBF (`src/libbf/mersenne-twister.c`). The vendored source
carries this statement (verbatim; the file embeds no longer licence text):

````text
/*
 * The Mersenne Twister pseudo-random number generator (PRNG)
 *
 * This is an implementation of fast PRNG called MT19937, meaning it has a
 * period of 2^19937-1, which is a Mersenne prime.
 *
 * This PRNG is fast and suitable for non-cryptographic code.  For instance, it
 * would be perfect for Monte Carlo simulations, etc.
 *
 * For all the details on this algorithm, see the original paper:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf
 *
 * Written by Christian Stigen Larsen
 * Distributed under the modified BSD license.
 * 2015-02-17, 2017-12-06
 */
````

Upstream (https://github.com/cslarsen/mersenne-twister) distributes this code
under the modified (3-clause) BSD licence; the BSD-3-Clause text is reproduced
in §2.11.

### 2.3 minizip — zlib licence

Zip archive I/O (`src/minizip/{ioapi,unzip,zip}.c`), used by the engine to
read `.data`/saved states. Copyright (C) 1998–2010 Gilles Vollant, (C)
2007–2008 Even Rouault, (C) 2009–2010 Mathias Svensson. Upstream minizip is
distributed under the same conditions as zlib; the zlib licence text is
reproduced in §3. The Info-ZIP decryption code referenced in its headers is
compiled out (`NOUNCRYPT`).

### 2.4 libtai — public domain

TAI/UTC date-time arithmetic (`src/libtai/`), by D. J. Bernstein, who placed
libtai in the public domain. The vendored copy carries no licence text.

### 2.5 dtoa — permissive notice

Float/string conversion (`src/os/dtoa.c`, included by `pl-dtoa.c`). Verbatim:

````text
/****************************************************************
 *
 * The author of this software is David M. Gay.
 *
 * Copyright (c) 1991, 2000, 2001 by Lucent Technologies.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 *
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR LUCENT MAKES ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 *
 ***************************************************************/
````

### 2.6 MurmurHash — public domain

Hash functions in the engine. Verbatim from `src/pl-hash.c`:

````text
/*  Part of SWI-Prolog

    Author:	Austin Appleby
    License:	Public domain
    See:	http://murmurhash.googlepages.com/
*/
````

## 3. zlib — zlib licence

Statically linked into the engine (and the licence §2.3 minizip is
distributed under). Verbatim from the pinned 1.3.2 tarball's `LICENSE`:

````text
Copyright notice:

 (C) 1995-2026 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu
````

## 4. Emscripten toolchain runtime

The emitted `swipl-web.js` and `swipl-web.wasm` embed Emscripten's JS glue and
the system libraries emsdk 6.0.1 links statically.

### 4.1 Emscripten — MIT / University of Illinois-NCSA (dual)

Verbatim `LICENSE` from the pinned emsdk's Emscripten:

````text
Emscripten is available under 2 licenses, the MIT license and the
University of Illinois/NCSA Open Source License.

Both are permissive open source licenses, with little if any
practical difference between them.

The reason for offering both is that (1) the MIT license is
well-known, while (2) the University of Illinois/NCSA Open Source
License allows Emscripten's code to be integrated upstream into
LLVM, which uses that license, should the opportunity arise.

The full text of both licenses follows.

==============================================================================

Copyright (c) 2010-2014 Emscripten authors, see AUTHORS file.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

==============================================================================

Copyright (c) 2010-2014 Emscripten authors, see AUTHORS file.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal with the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimers.

    Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following disclaimers
    in the documentation and/or other materials provided with the
    distribution.

    Neither the names of Mozilla,
    nor the names of its contributors may be used to endorse
    or promote products derived from this Software without specific prior
    written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

==============================================================================

This program uses portions of Node.js source code located in src/library_path.js,
in accordance with the terms of the MIT license. Node's license follows:

    """
        Copyright Joyent, Inc. and other Node contributors. All rights reserved.
        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to
        deal in the Software without restriction, including without limitation the
        rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
        sell copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in
        all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.
    """

The musl libc project is bundled in this repo, and it has the MIT license, see
system/lib/libc/musl/COPYRIGHT

The third_party/ subdirectory contains code with other licenses. None of it is
used by default, but certain options use it (e.g., the optional closure compiler
flag will run closure compiler from third_party/).

````

### 4.2 musl libc — MIT

Verbatim from the bundled musl's `COPYRIGHT` (the file's full contributor
enumeration is in `system/lib/libc/musl/COPYRIGHT` of the pinned Emscripten):

````text
musl as a whole is licensed under the following standard MIT license:

----------------------------------------------------------------------
Copyright © 2005-2020 Rich Felker, et al.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
----------------------------------------------------------------------
````

### 4.3 LLVM compiler-rt builtins — NCSA / MIT (dual, MIT elected)

Compiler runtime support compiled into `.wasm`. Verbatim legacy dual-licence
section from the bundled `compiler-rt/LICENSE.TXT` (the library is dual
licensed; this bundle elects the MIT terms, and both texts follow):

````text
==============================================================================
Legacy LLVM License (https://llvm.org/docs/DeveloperPolicy.html#legacy):
==============================================================================

The compiler_rt library is dual licensed under both the University of Illinois
"BSD-Like" license and the MIT license.  As a user of this code you may choose
to use it under either license.  As a contributor, you agree to allow your code
to be used under both.

Full text of the relevant licenses is included below.

==============================================================================

University of Illinois/NCSA
Open Source License

Copyright (c) 2009-2019 by the contributors listed in CREDITS.TXT

All rights reserved.

Developed by:

    LLVM Team

    University of Illinois at Urbana-Champaign

    http://llvm.org

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal with
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimers.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimers in the
      documentation and/or other materials provided with the distribution.

    * Neither the names of the LLVM Team, University of Illinois at
      Urbana-Champaign, nor the names of its contributors may be used to
      endorse or promote products derived from this Software without specific
      prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
SOFTWARE.

==============================================================================

Copyright (c) 2009-2015 by the contributors listed in CREDITS.TXT

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
````

### 4.4 dlmalloc — public domain (CC0)

The allocator compiled into `.wasm`. The bundled `system/lib/dlmalloc.c`
states (verbatim):

```text
This is a version (aka dlmalloc) of malloc/free/realloc written by
Doug Lea and released to the public domain, as explained at
http://creativecommons.org/publicdomain/zero/1.0/
```

### 4.5 MiniLZ4 — MIT

`-sLZ4=1` embeds LZ4 block decoding (minified) in `swipl-web.js` to unpack
`swipl-web.data`. Verbatim header from the bundled `third_party/mini-lz4.js`:

````text
/*
MiniLZ4: Minimal LZ4 block decoding and encoding.

based off of node-lz4, https://github.com/pierrec/node-lz4

====
Copyright (c) 2012 Pierre Curto

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
====

````

## 5. crosswordsmith — MIT

`crosswordsmith.qlf` (the qcompiled application), `worker.js`, and the SDK
facade. Verbatim from the repository `LICENSE`:

````text
The MIT License (MIT)

Copyright (c) 2011 Ned Letcher

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
````

## 6. Recorded future obligation — UKACD18

The browser bundle deliberately does **not** include the `fill` verb or any
dictionary today. The day `fill` browserifies with the UKACD18 word list, that
dictionary's **verbatim freeware licence must ship with the bundle** (the repo
README already carries this rule for native use: "redistributable freeware —
ship its license verbatim"). Add it to the summary table and to the packaged
assets in the same change that bundles the dictionary.
