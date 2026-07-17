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

## A-G2 transpose probes

- Premise report:
  `benchmarks/results/2026-07-17-a-g2-transpose-premise.md`.
- Premise launcher: commit `5500ebd`, path
  `benchmarks/run_arrange_g2_probe.pl`.
- Product report:
  `benchmarks/results/2026-07-17-a-g2-transpose-product.md`.
- Paired-wall launcher: commit `7a14082`, path
  `benchmarks/run_arrange_g2_wall.pl`.

The separate oracle `benchmarks/probe_arrange/g2_transpose.pl` remains live
until its focused product tests migrate in Phase 3.

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
