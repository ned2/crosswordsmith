# F-L1 — normalize parse loop: yall lambda → first-order filter (2026-07-05)

First product-code experiment of the fill campaign's Phase 3 (load layer).
Branch `experiment/f-l1-normalize` (base ed5d556). Subject:
`prolog/crosswordsmith/fill.pl` `normalize_word/2` (single caller:
`load_dict/3`, fill.pl:116; repo-wide grep confirms no other callers).

## Change

`include([C]>>char_type(C, alpha), Cs, Letters)` → `alpha_chars(Cs, Letters)`,
a first-order recursion making the IDENTICAL `char_type(C, alpha)` decision
per character; `string_upper`/`string_chars` untouched (frozen semantics:
same case-mapping path, same per-char decision, order preserved).

## Mechanism (measured, not taken on faith)

fill.pl never imports `library(yall)`, so its lambda is not compile-expanded:
`include/3` meta-calls the `>>` term per element — a lambda copy + yall
interpretation per character. Micro-decomposition on the real ENABLE char
lists (172,823 lines / 1,570,540 post-upper chars; prep outside the timed
goal), in a module that (like fill.pl) does NOT import yall:

| filter kernel                          | inf        | inf/char | wall ms |
|----------------------------------------|-----------:|---------:|--------:|
| include + yall lambda (status quo)     | 19,364,948 | 12.33    | 2,088   |
| include + named helper (V1 kernel)     |  5,230,087 |  3.33    |   202   |
| first-order alpha_chars (V2 kernel)    |  3,486,724 |  2.22    |   124   |

Control: with `library(yall)` imported, the lambda compile-expands and costs
the same as the named helper (5,230,088) — the premium is the RUNTIME
meta-call, ~9.0 inf/char × 1.57M chars ≈ 14.13M inf. The per-char model
predicts product `load_dict/3` totals to 1 inference: predicted 12,468,069
(V1) / 10,724,706 (V2) vs measured 12,468,070 / 10,724,707 (warm; the fixed
+1 first-warm-call JIT event as documented in Phase 0).

## Variants (warm load_dict/3 inferences via call_time, product path)

| variant | 10k load_inf (Δ vs 1,603,141) | 172k load_inf (Δ vs 26,602,930) | equivalence |
|---------|------------------------------:|--------------------------------:|-------------|
| V1 include(alpha_char) | 787,255 (−50.89%) | 12,468,070 (−53.13%) | same goal per element as the lambda — trivially identical |
| **V2 alpha_chars (SHIPPED)** | **686,601 (−57.17%)** | **10,724,707 (−59.69%)** | proven below |

V2 also removes include/3's per-element call/3 overhead — strictly better at
both scales, so it ships. PROVENANCE: V2's draft came from a terminated
duplicate session; it was treated as a hypothesis and verified independently
here (equivalence + determinism + measurement), not trusted as reviewed code.

V1/V2 equivalence vs the frozen reference `include([C]>>char_type(C,alpha))`:
elementwise `==` on 42 edge strings through the full pipeline (empty,
all-non-alpha, digits, punctuation, Greek/Cyrillic/CJK, ß→SS, combining
marks, soft hyphen/NBSP/ZWSP, Arabic-Indic + circled digits, emoji,
ligatures, Turkish dotless i, titlecase digraphs), 20,000 random char lists ×
(raw + upper-cased) over a wide code-point pool, and every line of all four
frozen dicts (172,823 + 50k + 25k + 10k). V2 additionally verified
deterministic (exactly one solution per input, like include/3).

## Ratchet (all 11 rungs, `--heavy`; three independent runs byte-identical)

| rung | search_inf Δ | load_inf before → after (Δ) | cmd_wall_med_ms before → after |
|------|-------------:|------------------------------:|-------------------------------:|
| sq04_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3910 → 2580 (−34.0%) |
| g11_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3810 → 2630 (−31.0%) |
| g11_full_seed | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3840 → 2630 (−31.5%) |
| sq05_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3980 → 2710 (−31.9%) |
| g17_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3980 → 2770 (−30.4%) |
| g21_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 3970 → 2850 (−28.2%) |
| g13_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 4130 → 2930 (−29.1%) |
| sq04_50k | +0.00% | 7,757,017 → 3,164,117 (−59.21%) | 1340 → 1000 (−25.4%) |
| g15_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 4350 → 2935 (−32.5%) |
| g17_50k | +0.00% | 7,757,017 → 3,164,117 (−59.21%) | 2310 → 2035 (−11.9%) |
| g09_full | +0.00% | 26,602,930 → 10,724,707 (−59.69%) | 5180 → 4030 (−22.2%) |

(Walls from the first `--heavy` pass on the final branch; run-to-run spread
across the three passes ≈ ±150 ms. Full-ENABLE rungs gain ~1.1–1.5 s of
command wall — the in-process warm load_dict wall dropped 3.47 s → 1.96 s.
cmd_rss ±0.4% max — unchanged.)

## Validation

- `./run_tests.sh` fully green (plunit + all goldens incl. `fill 15 bench`
  byte-identical + exit/stderr contracts); no golden regenerated.
- `check_fill_identity.sh`: 11/11 rung digests OK (run before AND after the
  final branch rebase-relocation; both clean).
- `check_fill_baseline.pl` core and `--heavy`: PASS, search_inf +0.00% on
  every rung, exit 0 — three independent pairs, count tables md5-identical.
- Baseline files untouched: no `--record`, fill_baseline.json /
  fill_history.jsonl / fill_identity.sha256 / fixtures/dict/* unmodified.

Beats the pre-registered expectation (−53% at 172k was V1's number; shipped
V2 lands −59.7%). Exceeds it because V2 also deletes include/3's per-element
call/3 overhead, exactly as pre-registered for V2.
