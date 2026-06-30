# Reference docs

## SWI-Prolog manual (local Markdown copy)

`swi-manual/` is a complete, text-searchable copy of the **SWI-Prolog reference
manual**, converted to Markdown from the HTML docs shipped with the SWI-Prolog
installed on this machine — so it matches the exact version the project runs
against (see the version header in [`swi-manual/INDEX.md`](swi-manual/INDEX.md)).

It exists so agents (and humans) can look up predicate signatures, modes, and
semantics from the repo, without guessing or going online.

> **License (third-party).** The `swi-manual/` tree is **not** covered by this
> repo's top-level MIT `LICENSE`. The SWI-Prolog reference manual is © the
> SWI-Prolog project and contributors and licensed under the **Creative Commons
> Attribution-ShareAlike 3.0 Unported License** (CC BY-SA 3.0,
> <https://creativecommons.org/licenses/by-sa/3.0/>); that notice is preserved
> in [`swi-manual/index.md`](swi-manual/index.md). These files are a Markdown
> *adaptation* (pandoc-converted; see the build script) and so remain CC BY-SA
> 3.0 — share-alike applies to this manual copy and adaptations of it, not to
> the rest of the repository (crosswordsmith's own code and data stay MIT).

- **Find a predicate** — grep the tree, e.g. `grep -rn "between(" docs/reference/swi-manual`,
  or grep the bold definition: `grep -rn '\*\*between\*\*' docs/reference/swi-manual`.
- **Read a topic** — open the relevant file; names are topic-based (`arith.md`,
  `lists.md`, `clpfd.md`, `assoc.md`, `apply.md`, …). The full topic→file map is
  in [`swi-manual/INDEX.md`](swi-manual/INDEX.md).
- **Packages** (clpfd, http, etc.) live under `swi-manual/packages/`.

These files are **generated** — don't hand-edit them; change the build script
instead.

### Regenerating

Run [`build-swi-manual.sh`](build-swi-manual.sh) to rebuild the tree from the
installed docs (requires `pandoc`). Re-run after an SWI-Prolog upgrade so the
manual tracks the version in use. The script reads the version and source path
straight from `swipl`, so there is nothing to configure.
