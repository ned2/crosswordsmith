# Plan — Purity/Cut Audit Doc Consistency Fix

**Status:** planned (not started)
**Effort:** S — ~15 min for the summary/severity/header/count core + C25 relocation; ~half a day for full internal consistency (the ~63 per-finding inline relabels).
**Type:** documentation only. **No code changes.** The one lurking code change — exporting `crosswordsmith_arrange:arrange_best_layout/6` — is **explicitly deferred and OUT OF SCOPE** here (see §5).
**Files touched:** `docs/prolog-purity-audit-findings.md`, `docs/STATUS.md`. No source or test files.

## Problem

`docs/prolog-purity-audit-findings.md` doubles as its own remediation tracker. The tracker + log at the bottom record the true end state — `:1859` "**The audit backlog is now fully dispositioned: 63/63 findings closed.**" — but the top of the document, and the per-finding bodies, still describe an all-open, nothing-applied backlog. Concretely stale:

- Executive summary bullet + severity table say "47 open items … none applied."
- Section-4 header/preamble say "NOT applied" / "Everything here is `status: open`."
- **61 literal `**open**` inline status markers** remain in the per-finding bodies (only C1 and C48 were flipped inline). Any grep/skim reads these as live.
- The `[ ] post-commit: make bench-fill-record` checklist box is unticked, but the record WAS run.
- The backlog also grew 47 → 63 (browser.pl addendum added C48–C63), so any "47" total is doubly stale.
- Two system-of-record drifts in `docs/STATUS.md`: no section exists for this (fourth) audit, and the one genuinely-deferred sub-item (C25's `arrange_best_layout/6` export) is buried inside a closed finding row instead of the de-accretion roadmap.

## Evidence base (verified 2026-07-06)

- True end state: `docs/prolog-purity-audit-findings.md:1859` (63/63 closed); tracker checklist `:1675–1742` (all `[x]`); tracker header `:1665` "Applied (Phase-2 APPLY steps — **committed** 2026-07-06)".
- `make bench-fill-record` was run: git `86ac895 bench(fill): record post-audit baseline (6 WINs)`, then re-recorded in `efa4289` (post-C7–C16) and `1701c78` (post-lows); `benchmarks/fill_baseline.json` + `benchmarks/fill_history.jsonl` mtimes 2026-07-06 16:59, last history row commit `085a44c` @ 16:23:24; doc's own log `:1769` "recorded (`make bench-fill-record`, baseline + history ledger updated)".
- C25 export genuinely still missing: `prolog/crosswordsmith/arrange.pl:20–33` export list omits `arrange_best_layout/6` (defined unexported at `arrange.pl:265`); `benchmarks/subjects.pl:57–65` reaches it white-box under a deliberate C25 annotation; `tests/arrange.plt` reaches `arrange_best_layout/5,6` module-qualified throughout.
- No other doc uses the "47 open" framing (grep of `docs/` + `docs/STATUS.md` clean). `docs/STATUS.md` has sections for the revamp audit (`:92–109`) and idiom audit (`:111–120`) but none for the purity/cut audit.

---

## 1. Executive summary + counts (the ~15-min core)

### 1a. Exec summary backlog bullet — `docs/prolog-purity-audit-findings.md:102–104`
OLD:
> - **Findings backlog: 47 open items** (C1–C47: **16 med · 24 low · 7 nit**), none applied in Phase 2. The top three (C1–C3) are the determinism agent's MEDs — all about *measurement/state integrity in persistent processes*, not solver correctness.

NEW:
> - **Findings backlog: 47 raised in the Phase-1 sweep** (C1–C47: **16 med · 24 low · 7 nit**), none applied *during Phase 2*; the browser.pl addendum later added 16 more (C48–C63). **All 63 are now dispositioned and closed** — see the [Remediation tracker](#remediation-tracker). The top three (C1–C3) are the determinism agent's MEDs — all about *measurement/state integrity in persistent processes*, not solver correctness.

### 1b. Severity table — `docs/prolog-purity-audit-findings.md:110–115`
Relabel the header cell and append a closure note. OLD header row:
> `| open severity | count | ids |`

NEW header row:
> `| severity (as raised) | count | ids |`

Leave the four data rows (`:112–115`) unchanged (they are an accurate as-raised snapshot of C1–C47), and insert immediately after `:115`:
> `> All 47 above, plus the 16 addendum findings (C48–C63), are now closed — **63/63** (see [Remediation tracker](#remediation-tracker).)`

### 1c. Related staleness — "staged, uncommitted" (same drift class; recommended in the same pass)
The tracker header `:1665` and git `fb4cfd9` show the APPLY steps were committed; two spots still say uncommitted.

- `:90–91` OLD: `- **18 experiment steps earned APPLY and are applied to this worktree** (staged, uncommitted; the overlapping X4.3/X5.1 pair merged via git 3-way into 17 distinct patches — see the [count note](#applied-step-count-note)).`
  NEW: replace `applied to this worktree` → `applied and committed 2026-07-06 (`fb4cfd9`)` and drop `(staged, uncommitted; ` → `(the overlapping X4.3/X5.1 pair merged via git 3-way into 17 distinct patches — see the [count note](#applied-step-count-note)).`
- `:460–461` OLD: `All steps below are **applied to this worktree (staged, uncommitted)**.`
  NEW: `All steps below are **applied and committed 2026-07-06 (`fb4cfd9`)**.`

## 2. bench-fill-record — mark done (record WAS run)

- `:97` OLD: `... inside the 0.5% tolerance). `make bench-fill-record` is owed after commit.`
  NEW: `... inside the 0.5% tolerance). `make bench-fill-record` was run post-commit (baseline + history recorded, `86ac895`).`
- `:463–464` OLD: `**Owed after commit: `make bench-fill-record`** to ratchet the fill baseline down to the new numbers.`
  NEW: `**Done post-commit:** `make bench-fill-record` ratcheted the fill baseline to the new numbers (`86ac895`; re-recorded after later batches in `efa4289` / `1701c78`).`
- Checklist box `:1672` OLD: `- [ ] **post-commit**: `make bench-fill-record` (ratchet the 6-WIN fill baseline down)`
  NEW: `- [x] **post-commit**: `make bench-fill-record` — ratcheted the 6-WIN fill baseline down, 2026-07-06 (`86ac895`)`

**Leave as-is** (append-only historical log entries; `:1769` already records the record being done): `:1758`, `:1769`. **Leave as-is** (per-finding fix-recipes, not status claims): the `make bench-record` / `make bench-fill-record` references at `:699`, `:730`, `:808`.

## 3. Section-4 header + preamble — `docs/prolog-purity-audit-findings.md:510` and `:512–514`

- `:510` OLD: `## 4. Findings backlog (C1–C47) — found, verified, NOT applied`
  NEW: `## 4. Findings backlog (C1–C47) — found & verified in Phase 1 (all now closed — see [Remediation tracker](#remediation-tracker))`
- `:512–514` OLD (first sentence): `Deduplicated across all eight Phase-1 reports; severity-ordered. Everything here is `status: open`. Sources are noted as (census / style-core / style-arrange / style-fill / style-small / stdlib / determinism).`
  NEW: `Deduplicated across all eight Phase-1 reports; severity-ordered. Each finding below is recorded **as raised** in Phase 1; its per-finding status line has been updated to the live disposition, and the authoritative record (with commit) is the [Remediation tracker](#remediation-tracker) and [log](#remediation-log) — **all findings are now closed**. Sources are noted as (census / style-core / style-arrange / style-fill / style-small / stdlib / determinism).`

> NOTE: The preamble deliberately does NOT claim "everything is open" nor install a banner in lieu of relabeling. It commits us to §4's inline relabel below.

## 4. Disposition-aware inline relabel of the ~63 per-finding status markers (the half-day step)

There are **61** remaining literal `**open**` markers (C1 and C48 already carry inline dispositions). This is a **per-finding, disposition-aware** edit, NOT a blind global replace: each finding's `**open**` token in its metadata line is replaced with the disposition recorded in that finding's Remediation-tracker row. The tracker (`docs/prolog-purity-audit-findings.md:1675–1742`) is the source of truth for every disposition.

Procedure: for each finding `Cn` (n = 2..47, 49..63), locate its `#### Cn —` header and the single `**open**` in its metadata line; replace `**open**` with the mapped token below (mirroring C1's existing inline form `**fixed 2026-07-06**`). Do not otherwise alter finding bodies — forward-looking "Fix:" prose and "→ re-baseline" notes stay as historical record.

Default mapping: `**open**` → `**fixed 2026-07-06**`. Explicit non-default cases (must be respected — verify against the cited tracker row):

| Finding | Replace `**open**` with | Tracker source |
|---|---|---|
| C2  | `**relocated → C48 2026-07-06**` (closed as relocated, not fixed) | `:1677` |
| C9  | `**closed via C51 2026-07-06**` (original site deleted) | `:1684` |
| C42 | `**wont-fix 2026-07-06**` (`clumped/2` churn) | `:1732` |
| C44 | `**fixed 2026-07-06 (partial — (c) reverted on measurement)**` | `:1734` |
| C46 | `**wont-fix 2026-07-06**` (recorded-only notes) | `:1736` |
| all other C2–C47, C49–C63 | `**fixed 2026-07-06**` | respective tracker row `:1679–1742` |

Notes:
- C1 (`:521`) and C48 (`:1185`) already carry inline dispositions — leave them.
- C6 (`:609`), C25 (`:824`), C61 all resolve to `**fixed 2026-07-06**` per their tracker rows (`:1681`, `:1707`, `:1740`); C25 additionally gets the STATUS.md cross-ref in §5.
- The addendum findings C49–C63 have their `**open**` markers in §5 (approx lines `:1317`–`:1464`); same mapping applies.

Acceptance for this step: `grep -c '\*\*open\*\*' docs/prolog-purity-audit-findings.md` returns **0**.

## 5. C25 deferred export — relocate to STATUS.md de-accretion roadmap (keep the trigger + cross-ref)

The genuinely-deferred sub-item (export `arrange_best_layout/6`) is currently buried in a `[x]`-closed row. Surface it where export-surface cleanup is tracked, and keep the round-trip reference. **The export itself is NOT performed here — it is a future `arrange.pl` change, out of scope.**

5a. Add one bullet to `docs/STATUS.md`, in the **De-accretion / retirement roadmap** section (after `:145`):
> - **Deferred (trigger: next time `arrange.pl`'s export surface changes):** export `crosswordsmith_arrange:arrange_best_layout/6`. The benchmark search sampler (`benchmarks/subjects.pl:57–65`) reaches it white-box because it is not on the export list (`arrange.pl:20–33`); every recorded `search_inf` baseline is defined against `/6`, so re-pointing at a public seam (e.g. `arrange_outcome/5`) would silently move gated counts. A sanctioned white-box annotation is in place (C25, commit `085a44c`); promote to a plain `use_module` import when the export list is next touched. Ref: [`prolog-purity-audit-findings.md`](../prolog-purity-audit-findings.md) C25.

5b. Keep the C25 tracker row closed but add the cross-ref. `docs/prolog-purity-audit-findings.md:1707`, append to the row's tail (currently ends `... arrange_best_layout/6 export deferred to a future arrange change*`):
> ` · tracked in [STATUS.md de-accretion roadmap](./STATUS.md#de-accretion--retirement-roadmap)`

Do NOT reopen C25 and do NOT open a new C-finding (the audit is closed).

## 6. Add the missing purity/cut-audit section to STATUS.md

`docs/STATUS.md` records the revamp audit (`:92–109`) and idiom audit (`:111–120`) but not this fourth audit — a system-of-record with no row for a whole audit is the same drift class we are fixing. Add a section mirroring the existing two (place it after the idiom-audit section, keep chronological order, so after `:120`):

> ### SWI-Prolog purity & cut audit (2026-07-06) — remediation done
>
> A whole-core cut/impurity census plus six controlled, benchmark-gated experiments (the first audit in the series to *measure* its recommendations), grounded in the version-matched SWI 10.1.10 manual. Verdict: **zero high-severity findings** — the cleanest surface in the series. Applied 18 measured cut/purity steps (14 cuts + 9 `once/1` deleted; −6% RSS on 21×21 heavy rungs, −17…−18% greedy wall, −0.78…−2.9% fill `search_inf` where counted work shrank), with a published negative control proving where cuts must stay. Backlog: **63 findings** (C1–C47 sweep + C48–C63 browser.pl addendum), now **63/63 dispositioned** (fixed, with C42/C46 wont-fix and C44(c) reverted on measurement). One export cleanup is intentionally deferred — see the de-accretion roadmap. Full per-finding record + remediation log: [`prolog-purity-audit-findings.md`](../prolog-purity-audit-findings.md). **Status: done as of 2026-07-06.**

(Verify the med/low/nit severity split across both waves against the addendum log `:1785` before writing, and adjust the parenthetical counts to match.)

## 7. Verification (run after edits; all must hold)

```
# No stale live-status framing remains in the audit doc:
grep -n '47 open\|none applied\|NOT applied\|is owed after commit\|Owed after commit' docs/prolog-purity-audit-findings.md   # -> no matches
grep -n 'staged, uncommitted' docs/prolog-purity-audit-findings.md                                                          # -> no matches
grep -n '\[ \] \*\*post-commit\*\*' docs/prolog-purity-audit-findings.md                                                     # -> no matches (box ticked)

# Every per-finding marker relabeled:
grep -c '\*\*open\*\*' docs/prolog-purity-audit-findings.md            # -> 0
grep -c 'wont-fix 2026-07-06' docs/prolog-purity-audit-findings.md     # -> >= 2 (C42, C46)
grep -n 'relocated → C48' docs/prolog-purity-audit-findings.md         # -> C2 row present

# No other doc carries the stale framing:
grep -rn '47 open' docs/ docs/STATUS.md                                # -> nothing

# STATUS.md drifts closed:
grep -n 'purity & cut audit' docs/STATUS.md                            # -> new section present
grep -n 'arrange_best_layout/6' docs/STATUS.md                         # -> de-accretion "Deferred" bullet present

# Deferral intact by design (NOT a regression):
grep -n 'arrange_best_layout' prolog/crosswordsmith/arrange.pl         # -> still no export-list hit
```

## 8. Scope boundary

- **In scope (doc only):** all edits in §§1–6.
- **Out of scope:** actually exporting `arrange_best_layout/6` and re-pointing `benchmarks/subjects.pl` / `tests/arrange.plt` off the white-box reach. That is a source change deferred to the next `arrange.pl` export-surface change; §5 only *records* it so it is not lost.
- **Do not** edit any `.pl` / `.plt` / benchmark file as part of this task.
