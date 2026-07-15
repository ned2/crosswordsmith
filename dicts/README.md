# Bundled scored lexicon — Collaborative Word List (clean-floor derivative)

`cwl50.dict` is the **opt-in scored dictionary** for `fill` (design-spec
§8.4a; bundled by decision pass **DP-9**, §10). It is a derivative of the
[Collaborative Word List](https://github.com/Crossword-Nexus/collaborative-word-list)
(CWL, Crossword Nexus) — **MIT-licensed**; the upstream license ships
verbatim as [`CWL-LICENSE`](CWL-LICENSE). It is *not* the `--dict` default:
pass it explicitly.

    ./crosswordsmith fill --grid grids/amer11.json --dict dicts/cwl50.dict

## Provenance (pinned — the artifact is byte-reproducible)

| | |
|---|---|
| Upstream commit | `2efe76e11ef311315e76d59700752733d69733d7` (`xwordlist.dict` last changed **2023-02-12**) |
| Raw `xwordlist.dict` sha256 | `a945a839a5f1e6f48caf9c8de446e5cd85f3567d7f62afcf54c6b738e8906ff4` (567,657 lines) |
| `cwl50.dict` sha256 | `a677198f46ce888eb20fbd5c165b37a5eaa882600fe337753d2ca8e390751667` (252,200 lines) |
| Regenerate | `scripts/fetch-cwl.sh` (default flags reproduce this file byte-identically) |

Two filters produce the derivative (rationale measured in
[`../benchmarks/fill_quality/README.md`](../benchmarks/fill_quality/README.md),
CWL sections):

1. **Digit-bearing words dropped** (992 lines, e.g. `MP3J;100`): the
   engine's A–Z fold squeezes digits, turning them into wrong short entries
   that score-descending order then plays eagerly.
2. **`score ≥ 50` floor baked in**: 50 is the documented clean-floor
   convention for 0–100 lists, and — decisive here — the largest CWL floor
   inside the engine's measured capacity envelope (the full list crashes
   dict load at SWI's default 1GB stack; the ≥30 band loads but can blow
   the global stack in §8.4c search on grids with full-length slots).
   `scripts/fetch-cwl.sh --min-score N` emits other floors from the same
   pin, with those caveats printed.

## Using it

- **No `--min-score` needed** — the floor is pre-applied, and the default
  `score ≥ 1` prune is a no-op on this file. `--min-score 75` or `90`
  select the stricter in-file bands.
- Scores are CWL's own 0–100 units. Band *meanings* are observed, not
  documented upstream (research: `docs/research/wordlist-scoring-2026.md`).
- **Register/staleness:** CWL is American-style and its data has been
  frozen upstream since **2023-02-12** — fine for the measured grids
  (means 78–83, beating ingrid_core on the same list), but expect no
  post-2023 vocabulary. UK/cryptics fills may prefer a user-supplied
  UKACD18 (unscored) or Spread the Wordlist (CC BY-NC-SA — never bundled).
- **Honest limit:** no loadable CWL floor fills the UK stock grids
  (`blocked_13a` at this floor defeats ingrid_core too). The bundled demo
  pairing is `grids/amer11.json` (fills in ~6s, mean 78.1, 0 below-clean).
- **`--seed`/`--shuffle` are slow on this file** (~2m45s at load regardless
  of grid, vs ~6s for the default path): the seeded equal-score-band
  permutation is O(n²) and this list's score-50 band has 85,800 entries.
  Measured 2026-07-16; tracked as a design-spec §8.5 backlog row.
