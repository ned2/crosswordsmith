# A-G2 greedy transpose partner synthesis

Date: 2026-07-17
Base: `d7b5dfb5c22be70f6c4fc24ca1c3ee3b5a5f31a2`
Branch: `experiment/arrange-g2`
Verdict: **KEEP**

## Change and soundness

The greedy best-effort/candidate pool now directly searches only
`topleft_across` and `topright`. Each completed source construction is rebuilt
as its literal square-grid transpose, with a fresh clue-number variable in every
synthetic `pw/8` and a fresh copy of the original dropped terms. Blocks retain
the historical order: all TLA seeds, all TLD partners, all TR seeds, all BL
partners. A failed source contributes neither source nor partner. The existing
single top-level memo reset remains shared by both direct blocks.

The independent A-G2 observer still directly searches all four corners. Before
the change it proved the premise over 108 manifest slots and 28 generated slots.
After the change it compared those direct records with the synthesized product
raw pools before sorting/diversity selection: all seven manifest rungs and all
three generated cases passed with zero setup, eligibility, placed/dropped order,
reward, geometry, clue-variable, normalized-assoc, raw-pool, or involution
mismatches.

## Inference result

Construction remains the manifest's one pinned direct construction and is not
the attempt-count metric. Postprocess starts from the same raw pool and is also
unchanged. Sweep is the complete product pool construction.

| rung | construction old -> new | delta | sweep old -> new | delta | postprocess old -> new | delta |
|---|---:|---:|---:|---:|---:|---:|
| bundled_17_candidates | 12,811 -> 12,811 | 0.00% | 227,945 -> 127,545 | -44.05% | 6,910 -> 6,910 | 0.00% |
| bundled_11_best_effort | 2,098 -> 2,098 | 0.00% | 15,550 -> 10,037 | -35.45% | 29 -> 29 | 0.00% |
| benchmark_08_candidates | 43,676 -> 43,676 | 0.00% | 293,482 -> 164,166 | -44.06% | 5,974 -> 5,974 | 0.00% |
| real_13x13_12w_best_effort | 36,087 -> 36,087 | 0.00% | 642,508 -> 358,826 | -44.15% | 90 -> 90 | 0.00% |
| real_15x15_18w_best_effort | 98,467 -> 98,467 | 0.00% | 1,136,472 -> 627,096 | -44.82% | 130 -> 130 | 0.00% |
| ladder_15x15_32w_best_effort | 368,867 -> 368,867 | 0.00% | 2,335,326 -> 1,241,744 | -46.83% | 233 -> 233 | 0.00% |
| ladder_21x21_80w_best_effort | 239,344 -> 239,344 | 0.00% | 9,362,822 -> 4,910,962 | -47.55% | 593 -> 593 | 0.00% |

The greedy core+heavy ratchet reports seven wins and zero regressions.

## Direct attempts and differential

| rung | old direct attempts | new direct attempts | delta | direct-vs-derived |
|---|---:|---:|---:|---:|
| bundled_17_candidates | 20 | 10 | -50% | 10 paired slots, 0 mismatches |
| bundled_11_best_effort | 20 | 10 | -50% | 10 paired slots, 0 mismatches |
| benchmark_08_candidates | 20 | 10 | -50% | 10 paired slots, 0 mismatches |
| real_13x13_12w_best_effort | 20 | 10 | -50% | 10 paired slots, 0 mismatches |
| real_15x15_18w_best_effort | 16 | 8 | -50% | 8 paired slots, 0 mismatches |
| ladder_15x15_32w_best_effort | 8 | 4 | -50% | 4 paired slots, 0 mismatches |
| ladder_21x21_80w_best_effort | 4 | 2 | -50% | 2 paired slots, 0 mismatches |

The generated setup-failure, dropped-word, and all-fit sets contributed another
28 old direct slots, all matched. The setup-failure set retained four symmetric
old failures as two failed source attempts and no synthetic entries. Semantic
replay reports completed direct constructions `10,4,10,10,8,4,2`; the second
value is lower only because six of that rung's ten direct setups fail.

## Dense alternating wall

`benchmarks/run_arrange_g2_wall.pl` serialized accepted A-G1 baseline
`/tmp/opencode/crosswordsmith-p0-integration` at `c584962` against this candidate.
Each subprocess warmed the in-process sweep once; the outer protocol discarded
two alternating warmup pairs and measured 21 alternating pairs. No samples ran
concurrently.

| rung | baseline median | candidate median | paired median ratio | paired p10-p90 | speedup |
|---|---:|---:|---:|---:|---:|
| ladder_15x15_32w_best_effort | 0.102292 s | 0.055150 s | 0.543596 | 0.530757-0.558346 | 1.84x |
| ladder_21x21_80w_best_effort | 0.432068 s | 0.236204 s | 0.550034 | 0.535498-0.566627 | 1.82x |

## Verification

- Baseline `swipl -q benchmarks/run_arrange_g2_probe.pl`: PASS.
- Candidate A-G2 probe: all 108 manifest and 28 generated old direct slots,
  zero mismatches; synthesized raw pools exact.
- Focused greedy tests: 13/13 PASS.
- `make test`: 418 assertions PASS, every golden and CLI contract PASS.
- Greedy identity: exact
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.
- Strict identity: all 15 core+heavy rows exact; latency row remains
  `90289af7db529b0132bae8bb910a18e90daa6f7b1c19de7cb97508883a56c81b`.
- Strict core+heavy ratchet: all 14 gated counts exact at `+0.00%`:
  `25,006; 662,138; 38,503,164; 95,381; 341,666; 929,466;
  13,601,499; 38,215,934; 10,338,363; 280,528; 4,246,246;
  6,555,246; 3,545,338; 248,068` (listed in campaign-request order).
- Latency-only inference reference `500,000,191` remains informational; the
  measured row was `500,008,138` and its identity digest was unchanged.
- No baseline, history, identity manifest, golden, or authoritative ledger was
  modified or re-recorded.

## Draft ledger entry

### A-G2 - derive greedy transpose partners - KEPT

- **Change/soundness:** directly search the TLA and TR greedy seed blocks under
  the existing shared memo reset, rebuild literal square-grid transpose partners
  with fresh clue variables and copied original dropped terms, and concatenate
  the four blocks in historical TLA/TLD/TR/BL order. Failed source setup omits
  both symmetric entries; strict eligibility remains after construction.
- **Identity:** the independent four-corner observer matched all 108 manifest
  slots and 28 generated slots against synthesized product pools with zero
  mismatches. Ordered raw pools, ties, dropped order, rewards, selected layouts,
  normalized assocs, candidate counts/distances, and CLI bytes are unchanged;
  greedy identity remains `ffc32b71...f1c29` and all strict identities remain
  exact.
- **Result:** direct attempts halve exactly by rung:
  `20->10,20->10,20->10,20->10,16->8,8->4,4->2`. Sweep inferences improve
  35.45%-47.55% with zero construction/postprocess or strict regressions. Dense
  serialized wall ratios are 0.5436 and 0.5500 (paired p10-p90 0.5308-0.5583
  and 0.5355-0.5666), approximately 1.84x and 1.82x.
- **Verification/verdict:** focused tests, full native suite, exact replay,
  greedy/strict identities, and both core+heavy ratchets pass. **KEEP.** The
  residual gap from an exact 2x inference/wall result is expected helper,
  transpose, scoring, raw-pool, and shared pure-memo work, not surviving direct
  partner searches.
