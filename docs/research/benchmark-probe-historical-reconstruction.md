# Historical benchmark probe reconstruction

Status: ACTIVE INDEX (2026-07-18). This file records how to recover source for
retired benchmark probes whose results remain in `benchmarks/results/`. It is a
recipe index, not a supported-command list or a second source archive. Use a
disposable checkout; do not run these recipes in the main worktree.

## P1 arrange backtracking

- Result: `benchmarks/results/2026-07-04-p1-retreat-wipeout-tables.txt`.
- Driver: commit `22e8c21`, path `benchmarks/probe_backtrack.pl`.
- Product instrumentation:
  `docs/experiments/p1-backtrack-instrumentation.patch`.

The previously cited commits `a74ffcf` and `0b32369` are absent from this clone.
The source-preservation commit contains the complete 183-line driver:

```sh
git show 22e8c21:benchmarks/probe_backtrack.pl
```

## Arrange performance campaign

The accepted/rejected mechanism reports are under `benchmarks/results/` and the
durable outcome summary is
`docs/research/arrange-performance-campaign-2026-07b-summary.md`. The complete
integrated probe tree, its 12 generated cliff fixtures, and its 64-seed manifest
remain available at the last source-bearing snapshot:

```sh
git worktree add /tmp/crosswordsmith-arrange-probes --detach fbea282
```

The campaign-era reconstruction points are:

| group | measurement base | source / accepted snapshot |
|---|---|---|
| Phase-0 authority and twin | `aa4d3c8` | integrated probe snapshot `3ce0b8d` |
| P-D0 support delta | `1bccf47` | probe `d99ea4f`; integrated `9ee66ed` |
| P-C0 duplicate work | `1bccf47` | probe `4e4ee15`; integrated `55a0fde` |
| P-R0 fixed-instance pilot | `1bccf47` | runner/analysis patch `b0d53e7` |
| A-G1 premise / product | `3ce0b8d` | `ea65c1a` / `802a1ea`; integrated `88074b2` / `ef2e0b0` |
| A-G2 premise / product | `c584962` / `d7b5dfb` | `7d82081` / `c03624a`; accepted `7a14082` |
| A-D1 direct buckets | `362b399` | probe `b3e965e`; accepted `6bb50c8` |
| A-D2 newest-source delta | `83b7d89` | probe `c8ff053`; accepted `9d7055a` |
| Closeout | `5dcf975` | measurement code `c65df0a`; report `76e956b` |

The A-G1 semantic-counter twin survived beyond the campaign in
`benchmarks/greedy_subjects.pl`; its final source-bearing snapshot is `2a8108a`.

### P-R0 base-locked runner

`run_pr0.py` deliberately requires branch `probe/a-r0` with HEAD exactly at
`1bccf47917fea074a404bf467e51865ce988b8b8`. Its measured source was a working
patch over that base, later captured by `b0d53e7`. Reconstruct it in a disposable
clone without advancing HEAD:

```sh
git clone . /tmp/crosswordsmith-p-r0
git -C /tmp/crosswordsmith-p-r0 switch -C probe/a-r0 1bccf47
git diff --binary 1bccf47 b0d53e7 -- \
  benchmarks/probe_arrange/run_pr0.py \
  benchmarks/probe_arrange/analyze_pr0.py \
  benchmarks/probe_arrange/schema.py \
  benchmarks/probe_arrange/test_schema.py |
  git -C /tmp/crosswordsmith-p-r0 apply
git -C /tmp/crosswordsmith-p-r0 status --short
python3 /tmp/crosswordsmith-p-r0/benchmarks/probe_arrange/run_pr0.py
python3 /tmp/crosswordsmith-p-r0/benchmarks/probe_arrange/analyze_pr0.py \
  /tmp/crosswordsmith-p-r0/benchmarks/results/2026-07-17-p-r0-pilot.jsonl
```

The 131-row raw result
`benchmarks/results/2026-07-17-p-r0-pilot.jsonl` remains tracked for independent
cutoff-curve, paired-corner, reward, layout, and timing reanalysis. The generated
cliff fixtures and seed manifest do not: their parameters, hashes, and master
seed are recorded in the P-R0 report, and their exact bytes remain in `fbea282`.

### A-G2 transpose probes

- Premise report:
  `benchmarks/results/2026-07-17-a-g2-transpose-premise.md`.
- Premise launcher: commit `5500ebd`, path
  `benchmarks/run_arrange_g2_probe.pl`.
- Product report:
  `benchmarks/results/2026-07-17-a-g2-transpose-product.md`.
- Paired-wall launcher: commit `7a14082`, path
  `benchmarks/run_arrange_g2_wall.pl`.
- Four-corner oracle: commit `7d82081`, path
  `benchmarks/probe_arrange/g2_transpose.pl`.

The oracle's permanent geometry, fresh-variable, block-order, setup-failure,
and dropped-term invariants now live in `tests/arrange.plt`.

## P-F1 attribution

- Measurement base: `a178923`.
- Verbatim six-file probe snapshot: `22e8c21`.
- Result: `benchmarks/results/2026-07-05-p-f1-attribution.md`, landed at
  `93972d1`.

The original report-era object `0329b2f` is absent. Reconstruct the measured
base plus later-preserved probe source in a disposable checkout:

```sh
git worktree add /tmp/crosswordsmith-p-f1 --detach a178923
git diff --binary 22e8c21^ 22e8c21 -- benchmarks/probe_f1/ |
  git -C /tmp/crosswordsmith-p-f1 apply
```

## F-H2 gate probe

- Measurement base: `397615d`.
- Verbatim `index_stats.pl` and `kernel_bench.pl` snapshot: `22e8c21`.
- Result: `benchmarks/results/2026-07-05-f-h2-gate-probe.md`, landed at
  `ed5d556`.

The original object `62dc2b4` is absent. Reconstruct with:

```sh
git worktree add /tmp/crosswordsmith-f-h2-gate --detach 397615d
git diff --binary 22e8c21^ 22e8c21 -- \
  benchmarks/probe_fh2/index_stats.pl \
  benchmarks/probe_fh2/kernel_bench.pl |
  git -C /tmp/crosswordsmith-f-h2-gate apply
```

Use the SWI native/WASM versions pinned in the result report; current runtimes
are not equivalent evidence.

## F-H2 Phase A/B

- Measurement base: `7181337`.
- Initial implementation and Phase A/B probes: `a2248bd`.
- Optional-mask follow-up and final `build_probe.pl`: `0aa13ca`.
- Accepted merge: `2cb7c0f`.
- Result: `benchmarks/results/2026-07-05-f-h2-bitset-count.md`.

Phase A needs the full implementation patch, not only `phase_a.pl`, because its
archived source calls the post-F-H2 fill seam:

```sh
git worktree add /tmp/crosswordsmith-f-h2 --detach 7181337
git diff --binary 7181337 a2248bd |
  git -C /tmp/crosswordsmith-f-h2 apply
```

For the final optional-mask follow-up, apply through `0aa13ca` instead:

```sh
git worktree add /tmp/crosswordsmith-f-h2-final --detach 7181337
git diff --binary 7181337 0aa13ca |
  git -C /tmp/crosswordsmith-f-h2-final apply
```

At `a2248bd`, build fresh schema-v2 artifacts at the probe's historical
`*_v1.idx`/`*_v2.idx` names. For the `0aa13ca` follow-up, pass `--masks` when
building the mask-bearing Phase B artifacts; default artifacts intentionally
omit masks.

## Fill-quality campaigns

The complete pre-cleanup narrative and tables remain in the last source-bearing
README snapshot:

```sh
git show 0e7d124:benchmarks/fill_quality/README.md
```

### FS-3(b) / FS-4 frontier

- Driver: commit `23e3772`, path `benchmarks/fill_quality/matrix.sh`.
- Durable decision summary: `docs/design-spec.md` DP-6.

```sh
git show 23e3772:benchmarks/fill_quality/matrix.sh
```

### DP-7 MAC / dom-wdeg adoption

- Final probe: commit `4653996`, path
  `benchmarks/fill_quality/probe_mac.pl`.
- Durable adoption record: `docs/plans/fill-mac-dwd-implementation.md`.

```sh
git show 4653996:benchmarks/fill_quality/probe_mac.pl
```

### B0 shipped-MAC instrumentation

- Initial rig: commit `42c7eac`.
- Corrected budget-boundary rig: commit `9a88f36`, paths
  `benchmarks/fill_quality/probe_mac_b0.pl` and
  `benchmarks/fill_quality/run_b0_instrument.sh`.
- Result: `benchmarks/results/2026-07-16-fill-b0-mac-instrumentation.md`.

```sh
git show 9a88f36:benchmarks/fill_quality/probe_mac_b0.pl
git show 9a88f36:benchmarks/fill_quality/run_b0_instrument.sh
```

Later B1/B2/B3/C1 variants were never merged into the main source tree. Their
probe commits and verdicts remain in `docs/experiments.md` and
`docs/research/fill-perf-program-closeout-2026-07.md`.
