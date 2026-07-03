# `xword` Phase 0 — technical spikes

Companion to [`xword-spec.md`](xword-spec.md). Phase 0 de-risks Phase 1 by
running small **throwaway** probes and folding **durable findings** back into
concrete stories and spec corrections. This file is the backlog + the findings
log; it is the **system of record** for Phase 0.

## Operating rules

- **Code is ephemeral; findings are durable.** Run each spike in `scratch/` or
  an isolated worktree and **discard the code**. The only persistent output is a
  Findings entry here (and any correction folded into `xword-spec.md`). No spike
  code is committed to the package.
- **Spikes may overturn the spec.** A spike can conclude a tool/approach is
  wrong (e.g. "use resvg, not cairosvg"); record it and edit the spec.
- **Phase 0 gates Phase 1.** Phase 1 does not start until every spike marked
  **gates: P1** is `done`.
- **Agents may add spikes.** On completing a spike, if you surface an unknown
  that must be resolved first, **append a new spike** to the backlog with
  `gates: P1`. The seed set is a starting point, not a ceiling.
- **Deliverable shape.** Each finding must give the next agent enough to build
  without re-discovering: a working **snippet**, the **gotchas** hit, the
  **pro-tips**, and the **story-ready** recommendation.

## Findings entry template

```
### Sn — <title>   [status: done]
- **Question(s):** …
- **Tried:** … (fixtures/commands used)
- **Findings:** the answer(s).
- **Gotchas:** traps to warn agents about.
- **Pro-tips:** what made it work.
- **Snippet:** minimal copy-pasteable illustration.
- **Spec changes:** sections edited in xword-spec.md (or "none").
- **New spikes filed:** Sx, … (or "none").
```

## Backlog

| # | Spike | Gates | Status |
|---|---|---|---|
| S1 | CLI & app skeleton (cyclopts) | P1 | **done** |
| S2 | Terminal grid rendering (Rich) | P1 | **done** |
| S3 | Determinism & golden capture | P1 | **done** |
| S4 | ipuz archaeology + round-trip | P1 | **done** |
| S5 | Exolve archaeology + round-trip | P1 | **done** |
| S6 | SVG→raster backend (recon) | P3 | **done** |

`native` gets no spike: it is our own well-specified format (`json-output-spec.md`)
and is exercised as a fixture inside S2–S5.

## Spike definitions

### S1 — CLI & app skeleton (cyclopts)  ·  gates: P1
- **Questions:** How do `view`/`convert`/`render` map to `cyclopts` subcommands
  with **stdin as default input** (positional file overrides), `--from`/`--to`
  as enums, boolean flags (`--blank`, `-q/--quiet`), per-verb `--help`, non-zero
  exit on error, and clean error rendering? What is the package layout (`xword/`
  src, `pyproject.toml`, `xword` console entry point), the `[raster]` optional
  extra stub, and the dev/run story via `uvx`/`pipx`/editable + `pytest`?
- **Try:** a minimal cyclopts app whose three verbs echo parsed args; wire
  TTY/stdin detection; a `pyproject.toml` with entry point + extra; run through
  `uvx` and an editable install.
- **Fold back:** the skeleton, cyclopts idioms + gotchas, dev-run commands;
  resolves the repo-layout open item (§14) and §12.

### S2 — Terminal grid rendering (Rich)  ·  gates: P1
- **Questions:** Draw a grid that **reads as square** despite the terminal's
  ~2:1 char aspect ratio, with corner numbers, block cells, and solved letters.
  Rich `Table` vs manual `Segment`/`Text`? Cell geometry (chars w×h), box-join
  artifacts, block styling (reverse-video vs `█`), number placement, colour +
  `NO_COLOR`, `wcwidth`/Unicode traps, wide-grid overflow.
- **Try:** render `bundled_17` solved and `--blank` both ways; eyeball
  squareness; test non-TTY and `NO_COLOR`.
- **Fold back:** chosen primitive + why, cell-geometry decision (resolves §14),
  gotchas, snippet. (Also validate the Rich↔Textual renderer seam, S7, here.)

### S3 — Determinism & golden capture  ·  gates: P1
- **Questions:** Produce byte-stable output for goldens across sinks — Rich→text
  with **pinned `Console(width=…)` + colour disabled** in a non-TTY CI, JSON key
  ordering for `convert`, and the SVG no-timestamp / no-random-id rule. How to
  capture Rich deterministically (record console / `export_text`)?
- **Try:** capture S2's render to text under fixed width + `NO_COLOR`; assert
  stable across runs; dump JSON with stable ordering.
- **Fold back:** the golden-capture recipe + pitfalls (width, ANSI, locale);
  confirms D6/§11.

### S4 — ipuz archaeology + round-trip  ·  gates: P1
- **Questions:** ipuz v2's real encodings — blocks, numbering, `puzzle` vs
  `solution`, enumeration, title/author, rebus, styling, custom extensions — and
  detection sniff (`version`/`kind`). Does the §10 capability table hold?
- **Try:** parse real ipuz samples **and** `crosswordsmith export --to ipuz`
  output into a prototype `Board`; round-trip `native→ipuz→native`; confirm
  `link` drops with a warning and structure maps; diff structurally vs the
  engine.
- **Fold back:** field-level native↔ipuz mapping, encoding gotchas, any §10
  correction, detection snippet, a sample fixture.

### S5 — Exolve archaeology + round-trip  ·  gates: P1
- **Questions:** the `exolve-*` directive grammar (begin/width/height/grid/
  across/down), block/bar/circle encoding, title/setter, and the
  **enumeration-embedded-in-clue-text** parse-back ambiguity. Detection; §10
  validation.
- **Try:** parse `crosswordsmith export --to exolve` output + real Exet-exported
  samples; round-trip; probe the enumeration-extraction ambiguity for a rule.
- **Fold back:** field mapping, the enumeration-parse rule + gotcha, detection
  snippet, a sample fixture.

### S6 — SVG→raster backend (recon)  ·  gates: P3 (front-loaded)
- **Questions:** which `[raster]` backend — `cairosvg` (needs cairo **system**
  libs) vs `resvg` (static wrapper) vs headless browser — for install friction,
  determinism, and PDF support? Must the SVG be shaped raster-friendly?
- **Try:** a tiny SVG → PNG/PDF through each candidate; note system-dep pain;
  check deterministic bytes.
- **Fold back:** backend recommendation + the `[raster]` install story; any SVG
  design constraints for Phase 3.

## Findings

### S1 — CLI & app skeleton (cyclopts)   [status: done]
- **Question(s):** How do view/convert/render map to cyclopts subcommands with
  stdin-as-default (positional file overrides), `--out` on success, `--from`/
  `--to` as enums, boolean flags (`--blank`, `-q/--quiet`), per-verb `--help`,
  non-zero exit + clean error rendering (diagnostics→stderr, data→stdout), and a
  bare `xword` that prints usage and exits non-zero? Package layout, `[raster]`
  extra, dev/run story via uvx/pipx/editable + pytest?
- **Tried:** Built a real src-layout package (hatchling) in `scratch/phase0/s1/`.
  Versions: **cyclopts 4.20.0, rich 15.0.0, uv 0.11.23, Python 3.14.6**. Ran
  every path (bare `xword`; `echo … | xword view`; positional override; enum +
  stderr warn + `-q`; per-verb `--help`; missing `--to`; bad enum; missing file)
  checking `echo $?` and the stdout/stderr split; forced the png→TTY and
  no-input→TTY guards under a real pty; `uv run pytest`; `uvx --from . xword`;
  `uv run --extra raster` (imported cairosvg 2.9.0); `uv tool install --force .`.
- **Findings:**
  - **Subcommands:** `app = App(name="xword", version=…)`; `@app.command` per
    verb; docstring → verb help; `--help`/`-h` per verb + root `--version` free.
  - **Enums:** plain `enum.Enum` for `--from`/`--to`; cyclopts validates, lists
    choices in help, errors with exit 1 on a bad value. Required `--to` = a
    keyword-only param with no default.
  - **`--from` keyword clash:** name the param **`from_`** — cyclopts strips one
    trailing underscore, so it maps to `--from` automatically.
  - **Short aliases / negation:** `Parameter(name=("--quiet","-q"))`; booleans
    auto-get a `--no-<flag>` negation.
  - **stdin/positional/--out:** ordinary code, not a cyclopts feature. Type the
    positional as `cyclopts.types.ExistingFile` (existence validated → exit 1 on
    missing). Build the **whole output string, write once at the end** → the
    "written only on success, no partial file" guarantee is trivial.
  - **Errors/streams:** cyclopts renders all parse/validation errors as a Rich
    panel **to stderr** and exits **1** automatically. Domain errors need a
    manual hook: wrap `app()` in `main()`, catch your exception → print to
    stderr, `raise SystemExit(1)`.
  - **Bare `xword`:** default prints help and exits **0** — NOT what we want.
    Register `@app.default` that prints help to a stderr `Console` and raises
    `SystemExit(2)`; does not break subcommands/`--help`/`--version`.
  - **Layout/entry point:** src-layout — `src/xword/{__init__,cli}.py`,
    `[project.scripts] xword = "xword.cli:main"`, `[tool.hatch.build.targets.
    wheel] packages=["src/xword"]`. Entry point is `:main` (a wrapper), not
    `:app`, so domain errors map to exit 1.
  - **`[raster]` extra:** `[project.optional-dependencies] raster = ["cairosvg"]`
    (backend = S6); `uv run --extra raster` resolves/imports it, base omits it.
  - **Dev/run:** `uv sync` (editable), `uv run xword …`, `uv run pytest`,
    `uvx --from . xword …` (ephemeral from source), `uv tool install .` /
    `pipx install .` (PATH command) — all confirmed.
- **Gotchas:** bare-invocation exit is 0 by default (non-zero is NOT free — add
  the `@app.default` handler). Calling `app(tokens)` in-process raises
  `SystemExit` even on success → **test the CLI via `subprocess`**, not by
  calling `app()`. `--from` is a Python keyword → use `from_`. The auto
  `--no-<flag>` negation always exists. The png/pdf→TTY guard only fires on a
  real terminal stdout — test under a pty, not a normal pipe.
- **Pro-tips:** let cyclopts own parse/validation errors (already stderr + exit 1
  + formatted); hand-handle only your domain exception in a thin `main()`. Type
  the positional as `ExistingFile` for free existence validation.
  `if __name__ == "__main__": main()` enables `python -m xword.cli` for
  pty/subprocess tests.
- **Snippet:**
  ```python
  # src/xword/cli.py
  import sys
  from enum import Enum
  from pathlib import Path
  from typing import Annotated, Optional
  from cyclopts import App, Parameter
  from cyclopts.types import ExistingFile
  from rich.console import Console

  app = App(name="xword", help="Crossword layout multitool.", version="0.0.0")
  _stderr = Console(stderr=True)

  class DataFormat(Enum):
      native = "native"; ipuz = "ipuz"; exolve = "exolve"

  class XwordError(Exception): ...

  @app.default
  def _bare():                       # bare `xword` -> usage on stderr, non-zero
      app.help_print(console=_stderr); raise SystemExit(2)

  def _read(file):                   # positional file overrides stdin
      if file is not None: return file.read_text()
      if sys.stdin.isatty(): raise XwordError("no input: pass a FILE or pipe stdin")
      return sys.stdin.read()

  def _write(text, out):             # --out on success, else stdout
      out.write_text(text) if out is not None else sys.stdout.write(text)

  @app.command
  def convert(
      file: Annotated[Optional[ExistingFile], Parameter(help="Layout; omit for stdin.")] = None,
      *,
      to: Annotated[DataFormat, Parameter(name="--to", help="Target (required).")],
      from_: Annotated[Optional[DataFormat], Parameter(help="Override detection.")] = None,  # -> --from
      out: Annotated[Optional[Path], Parameter(name=("--out", "-o"))] = None,
      quiet: Annotated[bool, Parameter(name=("--quiet", "-q"))] = False,
  ):
      """Convert between data formats (any -> any)."""
      data = _read(file)
      if not quiet: print("[warn] dropping native link/meta", file=sys.stderr)
      _write(f"[convert] from={from_} to={to.value} bytes={len(data)}\n", out)

  def main():                        # console entry: xword = "xword.cli:main"
      try: app()
      except XwordError as e:
          print(f"xword: error: {e}", file=sys.stderr); raise SystemExit(1)
  ```
  ```toml
  [project]
  name = "xword"; requires-python = ">=3.11"
  dependencies = ["cyclopts", "rich"]
  [project.optional-dependencies]
  raster = ["cairosvg"]
  [project.scripts]
  xword = "xword.cli:main"
  [build-system]
  requires = ["hatchling"]; build-backend = "hatchling.build"
  [tool.hatch.build.targets.wheel]
  packages = ["src/xword"]
  ```
- **Spec changes:** §8 (bare `xword` non-zero needs `@app.default`; cyclopts
  auto-renders parse errors to stderr/exit 1 — done); §12 (src-layout, entry
  point `:main` wrapper, `[raster]` extra + dev/run commands — done); §6.2
  (`--to` required keyword enum, `--from` = param `from_` — captured in §12/§8).
- **New spikes filed:** none. (Determinism of Rich error/help panels under fixed
  width + `NO_COLOR` for goldens folds into S3.)

### S2 — Terminal grid rendering (Rich)   [status: done]
- **Question(s):** How to draw a terminal crossword that reads square at ~2:1
  char aspect, with corner numbers, blocks, and letters? `Table` vs manual
  `Segment`/`Text`? Cell geometry, block styling, number placement, colour/
  `NO_COLOR`, Unicode-width traps, wide-grid overflow (resolving §14). Does it
  port to a Textual widget (S7 seam)?
- **Tried:** Rendered the real `bundled_17.native.json` (17×17) **solved and
  blank** via **both** a Rich `Table` renderer and a manual box-drawing
  `Text`/`Segment` renderer over a shared renderer-agnostic `BoardGeom`.
  Measured squareness (chars × ~2.05:1 pixel aspect), captured colour ANSI +
  plain via `Console(record=True).export_text()`, probed `wcwidth` on every
  glyph, tested overflow into `width=40` with/without `soft_wrap`.
- **Findings:**
  - **Primitive: manual `Segment`/`Text`, NOT `Table`.** Both can render a
    correct square grid (69w×36h → pixel-aspect ~0.93, reads square), but
    `Table` (1) omits inter-row lines unless `show_lines=True`, and decisively
    (2) **always reflows to console width** — at `width=40` it shrinks columns
    to 1 char and drops letters. You cannot pin a Table to fixed geometry. Manual
    `Text` renders at true geometry regardless of width and, with
    `soft_wrap=True`, overflows cleanly (the spec's "grid never wraps" rule),
    gives full junction/background control, and is byte-stable for goldens (S3).
  - **Cell geometry (resolves §14): 3w×1h interior, single-line box borders.**
    Grid pitch 4w×2h = 2:1 cancels the terminal's 2:1 aspect → square. 17-grid =
    69 chars wide (fits 80 cols), 36 rows tall. (`w5×h2` "number-above-letter" is
    also square but overflows 80 cols and is very tall — reject as default, keep
    as a future `--big`.)
  - **Number placement:** superscript digits, top-left, in-line before the
    letter (`str(n).translate(superscript-map)`); `¹O `, blank `¹  `. Two digits
    fit exactly (`¹⁷O` = 3 cells). Non-start cells centre the letter (` O `).
  - **Blocks: solid `█`×3, NOT reverse-video (default).** Reverse-video is an SGR
    attribute that **vanishes when styles are stripped** (piped/`NO_COLOR`/golden)
    → a block then looks identical to empty. `█` is a real width-1 glyph that
    survives plain capture. (Offer reverse-video only as a colour-TTY flourish.)
  - **Colour + `NO_COLOR`/non-TTY:** numbers `bold cyan`, letters default/white,
    borders `grey42`. `enabled = stdout.isatty() and not NO_COLOR` →
    `Console(no_color=not enabled, force_terminal=enabled)`. Confirmed colour
    path has escapes, no-colour path none, plain output byte-identical across
    runs.
  - **Unicode widths:** no traps — `wcwidth` = 1 for `█`, box-drawing, and every
    superscript digit; `¹⁷O` = 3 cells. (Superscripts are EA-ambiguous width, so
    a rare CJK-ambiguous-wide terminal could widen them — note it.)
  - **S7 Rich↔Textual seam:** clean/low-risk. The renderer is a pure
    `BoardGeom → Rich renderable` (only `Text`/`Segment`); a Textual widget's
    `render()` returns Rich renderables → drops in over the identical
    `BoardGeom`. Only interactivity (cursor/highlight) is new = restyling cells.
- **Gotchas:** `Table` needs `show_header=False, padding=0, pad_edge=False`
  **and `show_lines=True`** (miss the last → no row separators). `Console.print`
  word-wraps by default — a 69-char line in a 40-col console becomes 71 lines;
  pass **`soft_wrap=True`** (or `no_wrap`/`overflow="crop"`). Reverse-video blocks
  vanish in stripped output — use `█`. For deterministic capture use
  `Console(width=…, record=True, file=io.StringIO(), legacy_windows=False)` +
  `export_text(styles=False)`.
- **Pro-tips:** keep the presentation split real — one `BoardGeom`
  (renderer-agnostic), one pure `geom → Text`. Render at true geometry, compare
  grid width to `console.width`; if wider, scroll (`soft_wrap`) or stderr-note —
  never let Rich reflow it.
- **Snippet:** (condensed; full in `scratch/phase0/s2/render_segment.py`)
  ```python
  from rich.text import Text
  SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
  JUNC = {(0,1,0,1):"┌",(0,1,1,0):"┐",(1,0,0,1):"└",(1,0,1,0):"┘",
          (0,1,1,1):"┬",(1,0,1,1):"┴",(1,1,0,1):"├",(1,1,1,0):"┤",(1,1,1,1):"┼"}
  def render(geom, *, w=3, grid="grey42", num="bold cyan", let="white"):
      R, C = geom.rows, geom.cols
      out = Text()
      def hline(ra, rb):
          for c in range(C):
              out.append(JUNC[(int(0<=ra<R), int(0<=rb<R), int(c>0), 1)], style=grid)
              out.append("─" * w, style=grid)
          out.append(JUNC[(int(0<=ra<R), int(0<=rb<R), 1, 0)], style=grid); out.append("\n")
      def cell(cl):
          if cl.block: return Text("█" * w, style=grid)
          letter = cl.letter or " "
          if cl.number is not None:
              s = (str(cl.number).translate(SUP) + letter)[:w].ljust(w)
              t = Text(); k = len(str(cl.number))
              t.append(s[:k], style=num); t.append(s[k:], style=let if cl.letter else ""); return t
          return Text(letter.center(w), style=let if cl.letter else "")
      hline(-1, 0)
      for r in range(R):
          out.append("│", style=grid)
          for c in range(C): out.append(cell(geom.cells[r][c])); out.append("│", style=grid)
          out.append("\n"); hline(r, r + 1)
      return out
  # console.print(render(geom), soft_wrap=True)   # <- never wraps
  ```
  Rendered sample (solved, plain default; top-left):
  ```
  ┌───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┐
  │¹O │ M │ E │²G │ A │ P │ O │ I │ N │ T │███│███│███│███│███│███│███│
  ├───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┼───┤
  │███│███│⁴E │ T │ E │ R │ N │ A │ L │ R │ E │ T │ U │ R │ N │███│███│
  └── (17×17, 69w×36h chars, reads square) ──┘
  ```
- **Spec changes:** §14 (remove colour/geometry open Q, record decision — done);
  §6.1 (concrete never-wrap rule + geometry/colour decision — done); §5.2
  (affirm the pure `BoardGeom → Rich renderable` seam — done).
- **New spikes filed:** none. S7 (Rich↔Textual seam) validated low-risk and
  pre-closed at spec level; the full Textual widget stays a Phase-4 spike.

### S3 — Determinism & golden capture   [status: done]
- **Question(s):** Byte-stable goldens across all three sinks — (1) Rich→text for
  `view`, (2) JSON for `convert`, (3) SVG for `render`? Is S2's Rich capture
  stable across runs *and* decoupled from the ambient terminal
  (`COLUMNS`/`TERM`/`NO_COLOR`/`LANG`/TTY)? What JSON key-ordering rule matches D6
  and the engine? SVG no-timestamp/no-random-id discipline?
- **Tried:** `uv init` (Python 3.14, rich 15.0.0). Minimal manual `Text` grid
  renderer over `bundled_17.native.json`; `capture.py` (recipe), `emit.py`, and
  `harness.sh` running each sink **twice + `sha256sum` compare**, then re-running
  the plain sink under `COLUMNS=40/200`, `NO_COLOR=1`, `TERM=dumb/xterm-256color`,
  `LANG=C`/`en_US.UTF-8`, `FORCE_COLOR=1`, and piped-through-`cat`; same sweep on
  the ANSI sink. JSON: walked every object for sort-order + compared the engine's
  on-disk key-token sequence to a `sort_keys=True` dump + re-dump. Every claim is
  an observed `sha256sum`, not asserted.
- **Findings:** recipe per sink, all byte-verified:
  - **Rich→text (`view`)** — S2's recipe holds and is **fully decoupled** from
    the terminal. *Plain golden:* `Console(width=<pinned>, no_color=True,
    force_terminal=False, record=True, file=io.StringIO(), legacy_windows=False)`
    → `print(r)` → `export_text(styles=False)`: same sha in all 8 perturbation
    cases (invariant under `COLUMNS`/`TERM`/`NO_COLOR`/`LANG`/`FORCE_COLOR`/
    piped). *ANSI golden:* same but `force_terminal=True, color_system=
    "truecolor"`, `export_text(styles=True)` — stable, contains `ESC[`, and also
    invariant under `NO_COLOR=1` (because `export_text` reconstructs ANSI from the
    recorded **segment buffer**, not live colour gating). The two decoupling
    levers: **pinned `width=`** and **`file=io.StringIO()`** (capture never
    touches real stdout → piping can't perturb it).
  - **JSON (`convert`)** — rule: **`json.dumps(obj, sort_keys=True, indent=2,
    ensure_ascii=False) + "\n"`**. This **matches the engine at the key level**:
    both `bundled_17.native.json` and `.ipuz.json` already have every object's
    keys sorted at every depth (0 unsorted objects), and the engine's key-token
    sequence is identical to the `sort_keys=True` dump. Deterministic across runs
    and **idempotent** (dump==redump fixpoint). So `sort_keys=True`, not insertion
    order. (Key-order parity, not byte parity — byte-parity still gated on the
    engine's whitespace style + invented ipuz title, both deferred.)
  - **SVG (`render`, Phase 3)** — element ids are **pure functions of grid
    coordinates**: `cell_id(r,c)="cell-r{r}-c{c}"`, `word_id(n,dir)="word-{n}-
    {dir}"`. No `uuid`, no counters, no `datetime` — those are simply banned from
    the SVG serializer. (S6 confirmed the cairosvg raster stage adds none either.)
- **Gotchas:** **unpinned width leaks the environment (the #1 trap)** — omitting
  `width=` makes `Console.width` follow `COLUMNS`, so two goldens differ; the
  golden's width is a *deliberate pinned constant* (S2/S3 used 200), not the
  ambient terminal. `NO_COLOR` does **not** strip an ANSI golden captured via
  `styles=True` (styles live in the record buffer) — use `styles=False`/
  `no_color=True` for a plain golden. Set `legacy_windows=False` so captures
  don't depend on host OS. `LANG` had **zero** effect (Rich operates on `str`) —
  but write goldens to disk with explicit `encoding="utf-8"`. JSON needs
  `ensure_ascii=False` to match the engine's UTF-8 literals.
- **Pro-tips:** capture to `io.StringIO()`, never real stdout — severs all
  TTY/pipe coupling in one move. Define **one project-wide pinned width constant**
  and feed every `view` golden so grids never wrap. Verify goldens with
  `sha256sum` on double-runs in CI. For JSON assert `dump==redump` (idempotence).
  For SVG a one-line lint — grep output for `uuid`/`urn:`/4-digit years — guards
  the no-random-id rule.
- **Snippet:**
  ```python
  import io, json, os, sys
  from rich.console import Console
  VIEW_WIDTH = 200  # pinned project-wide; grids never wrap at this width

  def capture_plain(renderable, *, width=VIEW_WIDTH):     # byte-stable, env-invariant
      con = Console(width=width, no_color=True, force_terminal=False,
                    record=True, file=io.StringIO(), legacy_windows=False)
      con.print(renderable); return con.export_text(styles=False)

  def capture_ansi(renderable, *, width=VIEW_WIDTH):      # keeps ANSI, env-invariant
      con = Console(width=width, force_terminal=True, color_system="truecolor",
                    record=True, file=io.StringIO(), legacy_windows=False)
      con.print(renderable); return con.export_text(styles=True)

  def make_view_console():                                 # LIVE view, not goldens
      enabled = sys.stdout.isatty() and not os.environ.get("NO_COLOR")
      return Console(no_color=not enabled, force_terminal=enabled or None)

  def dump_json(obj):                                      # matches engine key order
      return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"

  def cell_id(r, c): return f"cell-r{r}-c{c}"              # deterministic SVG ids
  def word_id(n, d): return f"word-{n}-{d.lower()}"
  ```
- **Spec changes:** §11 (name the exact `view` capture recipe + pinned-width
  requirement + ANSI variant — done); §11/D6 (JSON rule `sort_keys=True,indent=2,
  ensure_ascii=False`, reproduces engine key order — done); D6/§6.3 (SVG ids =
  pure functions of coords, `uuid`/`datetime` banned — done); §8 (confirm runtime
  colour gating `enabled = isatty() and not NO_COLOR` — already recorded via S2).
- **New spikes filed:** none. (D6 fully de-risked across all three sinks; the
  only open items — engine JSON *byte* parity via whitespace/title, and SVG
  `<text>` vs `<path>` for cross-machine raster goldens — are tracked in §10/§14
  and deferred, not Phase-0 gates.)

### S4 — ipuz archaeology + round-trip   [status: done]
- **Question(s):** ipuz v2 encodings for blocks, numbering, `puzzle` vs
  `solution`, enumeration, title/author, rebus, styling, extensions? Does the §10
  table hold? Reliable detection sniff? Does native→Board→ipuz agree structurally
  with the engine's `export --to ipuz`, and is the invented-title divergence real?
- **Tried:** Prototype ipuz↔Board parser/serializer (`scratch/phase0/s4/`),
  detection sniff, structural (not byte) comparison. Parsed the engine's
  `bundled_17.ipuz.json` and the native fixture; hand-built a real-world-shaped
  `sample.ipuz.json` exercising a block, numbering, title/author, a **circle**
  cell, a **rebus** (`"PH"`), and both clue-item encodings. Read the engine
  emitter `prolog/crosswordsmith/export.pl:93-150`.
- **Findings:** native↔ipuz mapping —
  | native (§6) | ipuz v2 |
  |---|---|
  | `gridLength` N | `dimensions.{width,height}` (independent ⇒ rectangular) |
  | block (`null`) | `puzzle` and `solution` both `block` (default `"#"`; engine emits `"#"`) |
  | numbered start cell | `puzzle[r][c] = <int>` |
  | non-start white | `puzzle[r][c] = empty` (default `"0"`; **engine emits int `0`**) |
  | `cell.letter` | `solution[r][c] = "<letter>"` (multi-char = **rebus** `"PH"`) |
  | across/down/number | *(re-derived by row-major scan; ipuz ints advisory)* |
  | `words[].answer` | *(derived from `solution` letters along cells)* |
  | `words[].meta.clue` | `clues.Across/Down[i]` — item forms `{number,clue,enumeration}` **or** `[number,"text"]` **or** bare string |
  | enumeration | `clues[].enumeration` `"(5,5)"` |
  | title/author | top-level `title`/`author` (+ copyright/publisher/notes/date…) |
  | `link`/`meta` | **no home** — dropped + warned |
  Encodings: styling = a cell **object** `{"cell":<n|empty>,"style":{...}}` with
  `shapebg:"circle"`, `color`/`highlight`, `barred:"TL…"`; extensions =
  namespaced top-level keys / `volatile` (droppable). **§10 table HOLDS — no
  correction.** native→Board→ipuz gave **full structural equality** with the
  engine export (puzzle, solution, dimensions, clues incl. enumeration all `==`).
  **Invented-title CONFIRMED:** engine ipuz has `title:"Untitled"`, **no
  `author`** (`export.pl:104` hardcodes it — the only structural-view diff).
  **Detection:** check ipuz **first** — `version` is an `ipuz.org` URL or `kind`
  is a list containing one — else native (`grid`+`words`+`gridLength`).
- **Gotchas:** `empty`/`block` are *declared* fields with defaults `"0"`/`"#"` —
  engine emits **int** `0`, real-world uses **str** `"0"`; compare against the
  declared value, never hardcode. ipuz `null` cell (void) ≠ block cell. Clue
  items come in 3 shapes; only the object form carries `enumeration` (else
  re-derive from the answer). Rebus inflates answers beyond cell count — don't
  map answer chars to cells positionally. Styled cells are dicts — a block test
  must special-case them (`cell.cell` = label/empty, `cell.style` = styling).
- **Pro-tips:** block-ness from `puzzle`, letters from `solution` (tolerates a
  missing `solution` = blank template, matching §5.1). Compare `convert` to the
  engine on a **structural view** `{dimensions,puzzle,solution,clues}` (strips
  the title divergence + key-order noise). Port the engine's `enum_clean`
  zero-run handling (`export.pl:80-90`) for pathological separators.
- **Snippet:**
  ```python
  def sniff(obj):  # obj = json.loads(text)
      ipuz_url = lambda v: isinstance(v, str) and "ipuz.org" in v
      k = obj.get("kind")
      if ipuz_url(obj.get("version")) or (isinstance(k, list) and any(ipuz_url(x) for x in k)):
          return "ipuz"
      if {"grid","words","gridLength"} <= obj.keys(): return "native"
      return "unknown-json"
  BLOCK = lambda o: o.get("block", "#")     # engine: "#"
  EMPTY = lambda o: o.get("empty", "0")     # engine: 0 (int!)  real-world: "0"
  def puzzle_cell(pc, block, empty):        # -> (kind, number, style)
      if pc is None: return ("void", None, None)
      if pc == block: return ("block", None, None)
      if isinstance(pc, dict):
          n = pc.get("cell")
          return ("block" if n == block else "cell", n if isinstance(n, int) else None, pc.get("style"))
      return ("cell", pc if isinstance(pc, int) else None, None)
  ```
  Keep `scratch/phase0/s4/sample.ipuz.json` (3×3: block, circle, `"PH"` rebus,
  title/author, both clue-item forms) as the Phase-1/2 ipuz fixture.
- **Spec changes:** none required — §10 validated as-is; invented-title
  (§6.2/§11/§14) confirmed and existing wording correct. Applied optional
  clarifications: §7 (sniff checks ipuz first; signal = a value containing the
  `ipuz.org` host) and §5.1 (ipuz `empty`/`block` are declared fields with
  defaults, read not hardcode — engine emits int `0`; styled cell = object).
- **New spikes filed:** none. (A bars/shaded end-to-end ipuz round-trip is a
  small Phase-2 test task, not a gating spike.)

### S5 — Exolve archaeology + round-trip   [status: done]
- **Question(s):** `exolve-*` grammar; block/bar/circle/unfilled encoding in
  `exolve-grid`; title/setter; and the load-bearing enumeration-parse-back rule
  (separate real clue text from a trailing `(5,5)` without misfiring on
  clue-internal parens). Detection; does §10 hold?
- **Tried:** Read the engine emitter (`export.pl` `layout_to_exolve/2`,
  `exolve_cell_char/2`, `answer_enumeration/2`) and the real-engine ingest test
  `tests/exolve_ingest_check.sh` (which drives `exolve-m.js` and already contains
  a parse-back rule); fetched the authoritative grid grammar from the Exolve
  README. Prototype `Exolve↔Board` (`scratch/phase0/s5/proto.py`), parsed the
  engine's `bundled_17.exolve.txt`, ran `native→Board→exolve→Board` structurally,
  stress-tested the enum extractor on ~20 adversarial clues. (No real Exet sample
  exists in-repo — only the engine golden `tests/golden/export_bundled_17.exolve`;
  Exet round-trip is manual per `docs/exet-verification.md`.)
- **Findings:** field mapping — width/height ↔ `exolve-width/height` (may differ
  ⇒ rectangular); block `null` ↔ `.`; `letter` ↔ letter char, **unfilled light =
  `0`**; numbering **re-derived by Exolve from geometry** (not stored per-cell);
  clue+enum ↔ one line `N clue (enum)` under `exolve-across:`/`-down:`;
  title/author ↔ `exolve-title`/`exolve-setter`; `link`/`meta` = **no home**,
  dropped. **Decorators are cell-char *suffixes*:** `.`=block, `0`=unfilled,
  `A–Z`=solution, `?`=undecided; suffix `|`=east-bar, `_`=south-bar, `@`=circle,
  `!`=prefilled, `*`=diagramless, `~`=skip-number, `{k}`=SVG-deco; rebus via
  `exolve-option: rebus-cells`. So `0@` = circled empty, `R|` = filled + east
  bar. **§10 HOLDS — no correction.**
  - **Enumeration-parse rule (normative):** the enum is the **last**
    parenthesised group anchored at end-of-string whose interior is only the enum
    alphabet — digits + separators `, - . ' ` and space:
    `\s*\(([0-9]+(?:[,\-. ']*[0-9]+)*)\)\s*$`. Text before it (rstripped) = clue.
    Keeps clue-internal parens (`Send (a letter) to HQ (4)` → clue `Send (a
    letter) to HQ`, enum `(4)`); `Ring (2 wds)` → no enum. Verified on hyphen
    `(4-5)`, apostrophe `(1'4)`, multi `(3,4,5)`, nested `Nested (a (b) c) (3)`.
  - **Where it still fails + the fix:** irreducibly ambiguous when the clue text
    itself ends in a bare `(n)` (`Beethoven's Fifth (5)`, `See 12 (4)`). The
    **disambiguator is the grid** — accept the enum only if its integer runs sum
    to the word's cell count, else treat the parens as clue text (and derive enum
    from cell count). Reduces but doesn't eliminate (`(10)` on a real 10-cell
    entry still collides).
- **Gotchas:** Exolve grids store **no spaces** — word boundaries live only in
  the enumeration; so `native→exolve→native` is not byte-lossless on the
  `answer` display string unless the enum is re-applied (puzzle-lossless per §10,
  expected). Grid rows are **not one-char-per-column** once decorators/rebus
  appear (`0@BR|S` = 6 chars for 4 cols) — tokenise main-char + trailing
  decorator run, never column-index. The engine **invents `exolve-title:
  Untitled`** (`export.pl` ~line 172); `xword` (invent-nothing) emits none —
  `exolve-title` is optional so this round-trips fine, **but**
  `docs/exet-verification.md` records Exet's *Save* crashing on a null title, so
  a downstream Exet re-save may hit it. Enum can be `(?)` (undecided) — strict
  digit rule treats it as clue text (extend alphabet with `?` if needed).
- **Pro-tips:** `tests/exolve_ingest_check.sh` is gold — it encodes both
  enum-derivation and a parse-back `strip_enum` and drives the real
  `exolve-m.js`; mirror its regex. `.strip()` each line before matching (sub
  directives are indented). Build numbering/words from grid geometry, not
  per-cell numbers (Exolve stores none).
- **Snippet:**
  ```python
  def is_exolve(text):  # §7 step 3, robust to indentation
      return any(s.startswith(("exolve-begin", "exolve-width"))
                 for s in (ln.strip() for ln in text.splitlines()))
  import re
  ENUM_RE = re.compile(r"\s*\(([0-9]+(?:[,\-. ']*[0-9]+)*)\)\s*$")
  def split_clue_enum(s):
      s = s.rstrip(); m = ENUM_RE.search(s)
      return (s.strip(), None) if not m else (s[:m.start()].rstrip(), "(" + m.group(1) + ")")
  def enum_letter_count(enum): return sum(int(n) for n in re.findall(r"[0-9]+", enum))
  # accept enum iff enum_letter_count(enum) == len(word.cells); else clue text
  ```
  Keepsakes for the Phase-2 corpus: `scratch/phase0/s5/sample_decorated.exolve`
  (4×4 with `.` block, `0` unfilled, `0@` circle, `R|`/`E_` bars, and a
  bracketed-ref collision clue) + the engine golden
  `tests/golden/export_bundled_17.exolve`.
- **Spec changes:** §5.1 (Exolve enum parse-back rule = end-anchored regex +
  cell-count cross-check; unfilled `0`; decorators-as-suffixes ⇒ tokenise —
  done); §6.2/§11/§14 (generalise the invented-title to **both** ipuz and Exolve
  + the Exet-Save-crash note — done). §10 table itself correct — no change.
- **New spikes filed:** none. (The Exet UI save-back remains a manual procedure
  per `docs/exet-verification.md`, not a code spike.)

### S6 — SVG→raster backend (recon)   [status: done]
- **Question(s):** Which `[raster]` backend converts the master SVG → PNG/PDF —
  `cairosvg` (system cairo libs) vs `resvg` (static Rust wrapper) vs headless
  browser — judged on install friction, byte-determinism (D6), and PDF support?
  Must the master SVG be shaped raster-friendly?
- **Tried:** Fresh `uv init` (Python 3.14, Linux). Authored a self-contained
  150×150 crossword-ish SVG (rects + corner numbers + centred letters, explicit
  `width`/`height` + `viewBox`, no external refs). Installs: `uv add cairosvg`
  (→ cairosvg 2.9.0 + cairocffi 1.7.1 + cffi/pillow/tinycss2/cssselect2/
  defusedxml), ran `svg2png`/`svg2pdf` twice each; `uv add resvg-py`
  (→ resvg-py 0.3.3, single static Rust wheel, zero transitive deps), ran
  `svg_to_bytes` twice. Determinism via `sha256sum`+`cmp`; magic/dims via `file`;
  PDF/PNG scanned for dates/ids. Headless browser assessed only.
- **Findings:** **`cairosvg` is the sole `xword[raster]` backend.** It is the
  only candidate emitting **both PNG and PDF** from one library, its output is
  **byte-identical across runs**, and its PDF is true **vector** (embedded glyph
  outlines / `/FontFile`, print-quality). `resvg-py` has the best install story
  (zero system deps) but emits **PNG only** — and PDF is required (D4/§6.3), so
  it cannot stand alone. Headless browser is heaviest (~130–170 MB chromium) and
  its PDFs carry non-deterministic metadata.

  | Backend | Install friction | Determinism (same machine) | PDF | PNG |
  |---|---|---|---|---|
  | **cairosvg** | needs **system libcairo** (cairocffi dlopens `libcairo.so.2`); clean box = `apt install libcairo2`; pure-Python otherwise | ✅ byte-identical (no `CreationDate`/`ModDate`/`/ID`, no PNG `tIME`) | ✅ vector | ✅ RGB |
  | resvg-py | ✅ zero system deps — static Rust wheel | ✅ byte-identical PNG | ❌ none | ✅ RGBA |
  | headless browser | ❌ heaviest (chromium download) | ⚠️ weak (dates/version drift) | ✅ | ✅ |

  **Install story:** the extra declares just `cairosvg`; the one non-pip
  prerequisite is the system cairo library (`apt install libcairo2` / distro
  equivalent). `uv add cairosvg` + convert works out of the box wherever cairo
  is already present (most dev machines / desktop distros); minimal CI/containers
  add the apt line.

  **SVG design constraints for Phase 3:** text rendering is **host-font
  dependent** (resvg with `skip_system_fonts=True` and no supplied font dropped
  all text). Same-machine output is deterministic (satisfies D6, §11's
  dimension/magic goldens); cross-machine pixel identity is **not** guaranteed
  with `<text>`. If cross-machine byte-goldens are ever wanted, emit glyphs as
  `<path>` or bundle+pin a font. Keep the master SVG self-contained (no external
  href/font/image refs), with explicit `width`/`height` + `viewBox`, and honour
  D6 at the SVG layer (stable ids, no timestamps) — the raster stage adds none.
- **Gotchas:** pypi package is **`resvg-py`** (import `resvg_py`), not `resvg`.
  `cairocffi` does **not** bundle cairo — it `dlopen`s `libcairo.so.2` at import;
  a box without it fails with `OSError`, not a pip error (so "it just worked"
  here only because this box has cairo). `resvg_py.svg_to_bytes` returns a
  Rust-side bytes-like — wrap in `bytes(...)` before writing. `skip_system_fonts
  =True` silently drops text unless you pass `font_files`/`font_dirs`. resvg PNG
  is RGBA, cairosvg PNG is RGB — the two backends' PNGs are never byte-comparable
  to each other (only within a backend).
- **Pro-tips:** cairosvg PDF is genuinely deterministic — nothing to scrub for
  D6. `svg2png(..., output_width=/scale=)` gives hi-DPI PNGs from one SVG. Verify
  raster goldens with `file` (dims + magic), not sha, for cross-machine CI.
- **Snippet:**
  ```python
  import cairosvg
  cairosvg.svg2png(url="grid.svg", write_to="out.png")   # + output_width/scale for hi-DPI
  cairosvg.svg2pdf(url="grid.svg", write_to="out.pdf")   # vector, print-ready, deterministic
  # in-memory: pass bytestring=..., omit write_to to get bytes back
  ```
  ```toml
  [project.optional-dependencies]
  raster = ["cairosvg>=2.9"]   # runtime prereq: system cairo (apt install libcairo2)
  ```
- **Spec changes:** §12 (record the cairosvg decision — done); §8/§6.3 (the
  "missing raster extra" error should also cover missing-system-cairo — done);
  §11 (note *why* raster goldens are dims/magic not byte: host-font-dependent
  `<text>` — done); §14 (add the `<text>`-vs-`<path>` question for the SVG→PDF
  path — done).
- **New spikes filed:** none. (One optional Phase-3 sub-task, not a Phase-0
  gate: decide `<text>` vs `<path>` glyphs for the SVG→PDF path — noted in §14.)
