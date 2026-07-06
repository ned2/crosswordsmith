# Fill dictionary Unicode/locale hardening

Status: REVIEWED 2026-07-07 (v2) — adversarially reviewed; four load-bearing
corrections folded (⊳). Resolves the deferral in
`docs/plans/fill-perf-campaign.md:176` ("Unicode/locale hardening of
dictionary normalization (noted, deferred)").

Branch: `harden/fill-dict-unicode` from current main. This plan lands at
`docs/plans/fill-dict-unicode-normalization.md` as the first commit.

## 1. The defect (probes 2026-07-07; rows as seen under a UTF-8 locale)

`load_dict/3` → `normalize_word/2` (fill.pl:236–238) is `string_upper` +
keep-`char_type(C,alpha)` (alpha_chars/2, :240–246):

| Input line | Today | Wart |
|---|---|---|
| `café` (NFC) | `C,A,F,É` | no fold; É matches only É |
| `café` (NFD) | `C,A,F,E` | same text, different letters by byte encoding |
| `Straße` | `S,T,R,A,ß,E` | ß survives string_upper in ALL locales (probed) |
| `наöquery` | `Н,А,Ö,Q,U,E,R,Y` | any Unicode alphabet passes |
| `don't` | `D,O,N,T` | correct — ASCII punctuation squeeze stays |

⊳ Deeper than the draft knew (review findings 5, 7): BOTH retained primitives
are locale-dependent for non-ASCII — `string_upper(é)`=É under UTF-8 locales
but unchanged under LC_ALL=C; `char_type(é, alpha)` flips alpha↔not-alpha by
locale — and `read_file_to_string(File, Str, [])` (fill.pl:219–222) decodes
by locale (UTF-8 dict under LC_ALL=C → U+FFFD mojibake). The hardening must
remove every locale-sensitive decision from the pipeline, including the
decode.

## 2. Policy: fold-then-strict, fully locale-independent

Per line (replacing normalize_word/2's body; all decisions table/range-based;
`string_upper` retained only for its locale-stable ASCII a–z→A–Z mapping —
non-ASCII case is handled by double-casing the fold table):

1. ⊳ Decode pinned: `read_file_to_string(File, Str, [encoding(utf8)])`
   (finding 7). Dictionary files are DEFINED as UTF-8; document in PlDoc +
   fixtures/dict/README.md.
2. `string_upper`, then per char, first match wins:
   - `A`–`Z` → keep. ⊳ Fast path = a 26-fact first-arg-indexed `az/1` lookup
     — measured 2.10 inf/char on full ENABLE, inference-count-IDENTICAL to
     today's char_type scan; an @>=/@=< range chain costs 3.00 inf/char, ~+15%
     load_inf — use az/1 (finding 13).
   - Explicit DOUBLE-CASED fold table (Latin-1 + Latin Extended-A precomposed
     letters, both cases → A–Z sequences): `é/É→E`-class base folds; multi-char
     `ß→SS` (Unicode casefold-anchored: nfkc_casefold("Straße")="strasse"),
     `Æ/æ→AE`, `Œ/œ→OE` (established English orthographic equivalence),
     `Ĳ/ĳ→IJ` (NFKC-anchored), `Ø/ø→O`. ⊳ `Þ/Ð` DEMOTED to drop-word — no
     normative crossword-convention source found in docs/; pure-invention
     transliterations don't ship (finding 10).
   - ⊳ Combining marks: squeeze ONLY marks in an explicit allowed set — the
     decomposition marks of the table's letters (grave/acute/circumflex/
     tilde/diaeresis/ring/cedilla/macron/breve/... — enumerate from the
     table's NFD forms) and only directly after a kept A–Z char. Any OTHER
     mark, or a mark in any other position → drop the WORD. This kills the
     draft's NFC≠NFD inconsistency (finding 9b): Romanian `ș` (Ext-B, not in
     table) now drops the word in BOTH NFC (unfoldable letter) and NFD
     (disallowed mark U+0326) forms.
   - ASCII non-letters (digits, punctuation, space) → squeeze the char
     (unchanged: don't→DONT, 3d→D). ⊳ Plus a small explicit non-ASCII
     squeeze-list for typographic punctuation (U+2019 ’, U+2018, U+2013/14,
     U+00AD soft hyphen) so "don’t"→DONT matches "don't" (finding 9a).
   - Any other char (unfoldable letters — Cyrillic/Greek/CJK/unmapped Latin —
     or unknown ≥0x80) → **drop the whole word**, count it.
3. Empty-after-normalization → skip line (unchanged).
4. ⊳ Dropped-word count: **unconditional stderr warning when count > 0**
   (single line, e.g. "fill: N dictionary words dropped (unrepresentable in
   A-Z)"). The draft's verbose-gating was WRONG per the repo's own INV-3
   contract, which names "dropped words" verbatim as unconditional
   (core.pl:200–205, design-spec.md:109/:143–146; finding 11). ASCII dicts →
   0 drops → no output → the quiet-by-default stderr checks stay green.

Post-fold: café(NFC)=café(NFD)→CAFE; Straße→STRASSE; наöquery→dropped+counted;
don't = don’t → DONT; ș-words → dropped (both encodings).

⊳ Rejected-alternative rationale corrected (finding 8): this SWI install DOES
bundle `library(unicode)` (utf8proc ext; unicode_nfc/2, unicode_nfkc_casefold/2;
Unicode 16.0.0; probed locale-independent) — but the WASM build does NOT ship
utf8proc and fill.pl is qcompiled into crosswordsmith.qlf, so importing it in
fill.pl would break the browser build. Hand-rolled table it is — and
library(unicode) becomes the TEST ORACLE instead: a plunit compares every fold
table entry against unicode_nfd + mark-stripping / nfkc_casefold so the table
can't silently disagree with Unicode (native tests only; never imported by
fill.pl).

## 3. Where + what changes

- `prolog/crosswordsmith/fill.pl`: normalize_word/alpha_chars → the pipeline
  above (internal, prose comments; keep the F-L1 first-order rationale);
  line_word/2 drop-word case + counting; the stderr warning at the load_dict
  level (INV-3 style, mirroring core.pl:200–205's examples); encoding(utf8)
  pin in read_file_lines; `load_dict/3`'s PlDoc (:195–204) gains the
  UTF-8-defined-input + normalization-policy + unconditional-warning
  sentences. Determinism tags re-verified by probe.
- ⊳ F-L2 index artifact (finding 12): `fill_save_index` (fill.pl:744, :804)
  bakes normalization into saved indexes; its dict_sha256 check cannot see a
  normalization change. If the artifact header has a version field, bump it;
  if not, add one or extend the integrity input — implementer picks the
  minimal mechanism; for shipped ASCII dicts old artifacts remain semantically
  identical either way (document this).
- Scope note: answers/clues/seeds input paths are NOT touched (own validation
  paths); browser.pl does not import fill (verified :49–63); xword consumes
  JSON output only. Dictionary-line hardening only.
- Docs, same commit: README `--dict` (policy + warning line);
  ⊳ `fixtures/dict/README.md` (describes normalize_word with stale line refs;
  add UTF-8 definition + policy; finding 14c); design-spec fill/dict text if
  any live sentence describes normalization (grep first);
  `docs/plans/fill-perf-campaign.md:176` annotated "(resolved 2026-07-07 →
  fill-dict-unicode-normalization.md)".

## 4. Tests (tests/fill.plt; all additive — no existing test pins old
normalize behavior, verified finding 14)

- One plunit per probe-table row incl. the NFC=NFD equality pair and
  don't=don’t; ș both-encodings drop; drop-count value; empty-after-norm skip.
- ⊳ Oracle test: fold table ≡ library(unicode) derivation (see §2).
- ⊳ ASCII-freeze: word lists from `fixtures/dict/enable1.txt` AND
  `fixtures/wordlist_sample.txt` (the CLI's DEFAULT dict, crosswordsmith:487–489
  — missed by the draft, finding 2) byte-identical pre/post change.
- ⊳ Encoding-pin test: a UTF-8 dict loads identically under LC_ALL=C (would
  mojibake without the pin; finding 7).
- ⊳ Locale matrix (finding 14a): run the new plunit group under LC_ALL=C,
  C.UTF-8, en_US.UTF-8 — identical results. The draft's tr_TR probe is
  REPLACED: tr_TR is not installed on this host and SWI silently falls back
  (flipping default encoding to latin-1 — a false pass; finding 5).
- Stderr contract: existing quiet-by-default checks stay green (0 drops on
  ASCII); a plunit asserts the warning fires on a drops>0 dict (run_tests.sh
  :198–199 only checks /filled/ absence — don't rely on it for this).

## 5. Benchmarks

- ⊳ With the az/1 fast path the ASCII load_inf delta is expected ~0 (measured
  identical at probe scale) — the GOAL is landing inside the ratchet's 0.5%
  tolerance with NO re-record (finding 13; the draft's "expect FAIL then
  re-record" is retired). Report the real deltas from
  `swipl -q benchmarks/check_fill_baseline.pl` and `--heavy`.
- If load_inf still moves past tolerance: re-record ⊳ WITH `--heavy` (g09_full
  is tier-heavy; a core-only record leaves heavy rungs stale — finding 4),
  knowing `--record` rewrites all five metrics per rung + host/swi and appends
  fill_history.jsonl (finding 4): diff the file and assert search_inf/grid_inf
  numerically identical; label the commit.
- ⊳ Either way, fix `fill_baseline.json`'s stale `generated_note` (still says
  load_inf is "INFORMATIONAL/reporting-only"; it has been gated since F-L1 —
  finding 3) — hand-edit if no re-record happens.
- `bash benchmarks/check_fill_identity.sh` — all digests must match (ASCII
  dicts ⇒ identical word lists ⇒ identical fills).
- `swipl -q benchmarks/check_baseline.pl` — arrange sentinel, +0.00%.

## 6. Verification recap (report real results, incl. anything that failed)
make unit + make test (before/after counts; delta = +new tests only) · fill
ratchet core + --heavy · fill identity · arrange ratchet · locale-matrix plunit
run · encoding-pin test · list_undefined · state explicitly that no WASM run
is required (inference_parity is arrange-only; the shipped qlf embeds old
fill.pl until the next wasm build — say so, don't skip silently).

## 7. Non-goals
Pure-ASCII dictionaries' word lists change ZERO (hard acceptance). No ICU/
utf8proc import in fill.pl (WASM). No answer/clue-path normalization. No
arrange changes. No golden changes. Þ/Ð transliteration does not ship.
