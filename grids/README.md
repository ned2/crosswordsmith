# Stock-grid library

A small, curated set of **pre-validated legal blocked grid templates**
(design-spec §8.3). These are bundled assets — not generated — that feed the
`lint` profiles and, later, the `fill` engine. There is no CLI verb.

## Schema (OD-5: a black-square mask is the single source of truth)

```json
{
  "name": "blocked-13a",
  "size": 13,
  "symmetry": "rot180",
  "mask": ["...#.....#...", "...#.....#...", "..."]
}
```

- `size` — N (square, odd).
- `mask` — N strings of length N. `#` is a block; any other character (we use
  `.`) is a light (white) cell.
- Slots (across/down lights) are **derived** on load by run-scanning the mask
  (`stockgrid.pl`); they are not stored — no redundant slot list.

## The shipped set (OD-6)

| file | size | notes |
|---|---|---|
| `blocked_13a.json` | 13×13 | 180°-symmetric, lights 3/5/13 |
| `blocked_13b.json` | 13×13 | 180°-symmetric, lights 3/4/5/13 |
| `blocked_15a.json` | 15×15 | 180°-symmetric, lights 3/4/5/7/15 |
| `amer11.json` | 11×11 | American-style, fully checked, 180°-symmetric, spanning 11-letter slot — promoted from the fill-quality benchmark (DP-9) as the demo pairing for the bundled `dicts/cwl50.dict` (fills in ~6s, mean 78.1) |

All are original, license-clean, and **PASS `lint --profile blocked-uk`**. This
is a starter set (real publications curate only dozens); it grows by adding more
validated masks here.

## Adding / validating a grid

Each grid must be legal under the blocked-uk rules (≥3-cell lights, ≥half
checking, no triple-unch, 180° symmetry, connectivity, no isolated cells). The
validator reuses `lint`:

```prolog
?- stockgrid_report('grids/your_grid.json').   % prints PASS/WARN/FAIL + reasons
```

`tests/stockgrid.plt` asserts every bundled grid stays legal, so a regression
fails the suite.
