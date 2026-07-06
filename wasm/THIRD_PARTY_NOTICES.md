# Third-party notices — wasm bundle redistribution manifest

> **Status: draft until publish** (wasm-sdk-strategy OQ-6). The SDK is consumed
> in-repo today (OQ-1); nothing here ships to npm/CDN yet. This file records
> the redistribution obligations the built `wasm/client/` bundle carries so
> the eventual publish inherits a manifest instead of an archaeology task.
> **Before any npm/CDN publish:** replace each summary below with the verbatim
> upstream licence text (or bundle the upstream files) and re-verify the
> licence of each pinned version. `build-manifest.json` points here.

The build outputs (`swipl-web.{js,wasm,data}`, `crosswordsmith.qlf`) are
compiled from, and therefore redistribute, the following (exact versions/
commits pinned in [`wasm/build/pins.sh`](build/pins.sh) and recorded per-build
in `wasm/client/build-manifest.json`):

| Component | What of it ships | Licence | Obligation |
|---|---|---|---|
| SWI-Prolog (swipl-devel + its bundled packages) | the whole engine compiled into `swipl-web.*`; boot + library code in `.data` | **BSD-2-Clause** | reproduce copyright notice + licence text in distributions |
| zlib | statically linked into the engine | **zlib licence** | notice must not misrepresent origin; no attribution required in binaries, but the notice text customarily ships |
| PCRE2 | statically linked into the engine (`library(pcre)`) | **BSD-3-Clause** | reproduce copyright notice + licence text; no endorsement use of contributor names |
| Emscripten runtime | the JS glue emitted into `swipl-web.js` and support code in `.wasm` | **MIT / University of Illinois-NCSA** (dual) | reproduce copyright notice + licence text |
| crosswordsmith (this repo) | `crosswordsmith.qlf` (the qcompiled app) | **MIT** (repo `LICENSE`) | reproduce copyright notice |

## Recorded future obligation — UKACD18

The browser bundle deliberately does **not** include the `fill` verb or any
dictionary today (strategy §6 phase 4). The day `fill` browserifies with the
UKACD18 word list, that dictionary's **verbatim freeware licence must ship
with the bundle** (the repo README already carries this rule for native use:
"redistributable freeware — ship its license verbatim"). Add it to the table
above and to the packaged assets in the same change that bundles the
dictionary.
