# Native schema uplift — title/author anchor, styling + enumeration riders; lossless-only conversion policy

Status: DRAFT (pre-review). Ratified in session 2026-07-07: extend the native
schema where the engine can compute-or-carry, keep every conversion lossless,
and REJECT the lossy best-effort mode (xword-spec §10 hatch (a)) as
unnecessary. This plan supersedes the §10 "deferred escape hatches" paragraph:
its decisions land in the spec as numbered Dn entries as each phase ships.

Decisions (fixed; the phases implement them):
- **Uplift**: optional `title` + `author` (anchor, P1+P2); closed-subset cell
  styling (rider, P3); enumeration separator charset `'` `.` (rider, P4).
- **Permanent fail-strict, no uplift, no lossy fallback**: rectangular grids
  (the engine's whole model is square-N×N arithmetic — a schema wider than
  every engine verb is a lie; ipuz is the rectangular carrier format) and
  rebus (breaks the single-letter crossing model; flattening produces
  *invalid* puzzles, not diminished ones). Better error text is allowed.
- **Bars: gated, not designed here** — the engine has barred-Ximenean lint
  RULES (lint.pl:64–266) but no bar representation; the representation belongs
  to docs/cryptic-layout-spec.md (status: proposed, unimplemented). Bars enter
  the schema when that spec lands, on its terms.
- **Compatibility rule**: every schema addition is additive-optional — absent
  fields mean today's semantics; old layouts stay valid; nothing new is ever
  required. (json-output-spec has no version field; additive-optional IS the
  compatibility story — state it in the spec.)

References verified 2026-07-07 @ 3df11bc (grep to re-confirm):
json-output-spec §6 output format (:123–268); export.pl json_read_dict :40,
ipuz invented title :157–161, exolve-id :230, exolve-title :235;
exet-verification.md:144 (Exet Save crashes on null title — the reason the
exolve default exists); xword board.py:60 (Board.meta models title?/author?),
native.py:8/:106 (native's title/author rejection), ipuz.py :126/:223
(title/author both directions); xword-spec §10 table + closing paragraph, §11
~:400-408, §14; the TWO affected engine goldens:
tests/golden/export_bundled_17.{ipuz,exolve}.

---

## P1 — engine: schema fields + export policy convergence (M)
Branch: `feat/layout-title-author`. Plan doc = this file, committed first.

1. **json-output-spec**: new §6.x — optional top-level `title` and `author`
   (strings) on the layout object. `arrange` never emits them (it has no
   source for them — non-goal below); they exist so downstream tools (xword,
   hand edits) can carry them and so `export` can pass them through. Record
   the additive-optional rule.
2. **export.pl — converge on the Q5 policy (invent nothing, except exolve's
   ecosystem default)**:
   - ipuz: emit `title`/`author` only when present in the layout dict; STOP
     inventing `"Untitled"` (:161). Keep the comment history — the Exet
     rationale moves to the exolve side where it belongs.
   - exolve: `exolve-title:` from the layout title, defaulting to `Untitled`
     ONLY when absent (the Exet Save-crash defense stays, cite
     exet-verification.md:144); `exolve-setter:` from author when present.
   - **Drop `exolve-id: crosswordsmith-export` (:230)** — exolve documents id
     as optional/auto-derived from a content signature; xword-spec ~:242-245
     already records that a constant id doesn't preserve solving state across
     edits anyway. RECOMMENDED; if review overturns, record why in the spec.
3. Engine plunit (tests/export.plt): pass-through present/absent × both
   targets; absent-field layouts (all existing fixtures) emit no title/author
   in ipuz and the default title in exolve.
4. **Goldens: exactly two change, intentionally** —
   export_bundled_17.ipuz loses `"title": "Untitled"`;
   export_bundled_17.exolve loses `exolve-id` (keeps `exolve-title: Untitled`
   — the bundled layout is titleless). Update via the sanctioned
   `make update-golden` path; document + date in the commit message per
   discipline. Add ONE new golden with a titled+authored layout fixture so
   pass-through is byte-pinned.
5. Verify: `make unit` / `make test` / `make golden` (diffs confined to the
   two goldens + the new one); arrange + fill ratchets as +0.00% sentinels
   (export is outside the gated paths — prove it, don't assume); grep the
   wasm surface for export usage (browser.pl is arrange-only — confirm; the
   shipped qlf embeds old export.pl until the next wasm build — state it).
6. Tracker: STATUS.md export row + a de-accretion-style note that the
   invented-title hack is retired; exet-verification.md gets a dated addendum
   line (behavior it recorded has changed: default now exolve-only-when-absent).

## P2 — xword: native carries title/author (S–M)
Branch: `xword/native-title-author`, AFTER P1 merges (the engine cross-check
tests consume engine export output/goldens).

1. native.py: parse + serialize optional `title`/`author` (Board.meta already
   models them, board.py:60); delete the title→native structural-failure path
   (:106 region) and re-pin its tests.
2. xword-spec: §10 table cells title/author ✗→✓ for native (D7 rule text
   unchanged — the table drives it); the "Deferred escape hatches" closing
   paragraph → resolved decisions block (uplift shipped; lossy mode rejected;
   rectangular/rebus permanent-strict; bars gated on cryptic-layout-spec);
   §11 cross-check upgraded: title/author lines/fields assert byte-parity;
   §14 byte-parity entry updated (see Stretch). xword-status.md rows same
   commit.
3. Round-trip contracts (test_convert.py): titled boards are now IDENTITY
   through both ipuz and exolve hops (the Q5 fail-strict re-pin is replaced —
   the way back now succeeds carrying the title). Titleless→exolve still
   gains the injected default with the Q5 warning; the return trip now
   succeeds carrying `title: "Untitled"` — document as "gains the ecosystem
   default, loses nothing" (inherent to Exet's requirement, NOT lossiness;
   assert it explicitly in a test).
4. Verify: `cd xword && uv run pytest` before/after counts (baseline 111
   passed / 3 skipped); CLI demo end-to-end on a titled and a titleless
   engine export — paste the emitted lines; `make unit` once (engine
   untouched in P2).

## P3 — rider: closed-subset cell styling (M; own go/no-go, after P1+P2)
Branch: `feat/layout-cell-styles` (engine) + xword changes — sequence like
P1→P2 if the engine side lands first, or one branch if small enough; the
phase OPENS with its own mini plan-commit deciding the schema shape.

- Subset rule: carry ONLY what all three formats express. Background shade
  colour qualifies today (ipuz StyleSpec `color` ↔ exolve-colour — xword
  already serializes that direction, exolve.py ~:228-234). Circled cells join
  ONLY if exolve expressibility is verified against xword's own exolve
  parser/emitter knowledge (no guessing from memory); otherwise shade-only.
- Schema shape proposal (to be confirmed in the phase's plan-commit): optional
  top-level `styles: [{row, col, color}, …]` — matches xword's internal
  (colour,row,col) model; additive-optional.
- Engine: carry-and-ignore everywhere EXCEPT export, which passes the subset
  through to both targets (full pass-through, not carry-but-drop — a
  carry-but-drop export would recreate the engine-vs-xword divergence this
  campaign exists to close).
- xword: §10 table styling row ✗→✓ for native (subset-scoped footnote);
  exolve→ipuz StyleSpec coarseness note stays (unchanged by this).

## P4 — rider: enumeration separator charset (S; own go/no-go)
- First check whether `'`/`.` enumerations appear in any real corpus the repo
  can see (xword test fixtures, docs/reference examples); if none, decide on
  principle and say so in the commit.
- Widen the native answer display-form charset to include `'` and `.`
  (json-output-spec §6.2 display rules; xword's enumeration reconstruction —
  the §10 bullet: `,`→space, `-`→hyphen, now `'`→`'`, `.`→`.`); anything else
  stays fail-strict. Tests both sides.

## Stretch — byte-parity endgame (decision, likely tiny)
After P1+P2 (+exolve-id drop), the ONLY remaining engine↔xword divergence is
JSON whitespace style. Probe whether xword's ipuz/native emit can match SWI's
json_write_dict style cheaply (indent/separator tuning); if yes, do it and
upgrade §11's cross-check to full byte-parity; if not, record whitespace as
the one permanent, documented gap and pin cross-checks at
byte-parity-modulo-whitespace. Either way §14's byte-parity entry closes.

## Non-goals
- Lossy conversion mode of any kind (hatch (a) is REJECTED, not deferred).
- Rectangular/rebus support in native, in any phase.
- Bars (await cryptic-layout-spec).
- An `arrange --title` input flag / json-input-spec changes — title/author
  enter via the layout JSON (hand edit or xword native emit). Revisit only on
  user demand.
- Engine arrange/fill/lint behavior — export + schema docs only.

## Verification matrix (every phase reports real results)
| Check | P1 | P2 | P3 | P4 |
|---|---|---|---|---|
| make unit + make test (engine) | ✓ | once | ✓ | ✓ |
| make golden (intentional diffs enumerated) | ✓ | — | ✓ | ✓ |
| uv run pytest (xword, before/after counts) | — | ✓ | ✓ | ✓ |
| arrange+fill ratchets +0.00% sentinel | ✓ | — | ✓ | — |
| CLI end-to-end demo pasted | ✓ | ✓ | ✓ | ✓ |
| tracker/spec rows same commit | ✓ | ✓ | ✓ | ✓ |
