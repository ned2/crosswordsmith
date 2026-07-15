# Fill benchmark dictionaries

Frozen wordlists for the `fill` performance ladder
(`benchmarks/fill_workloads.pl`). These are **benchmark fixtures**, not a
product asset — the product's bundled scored lexicon lives in
[`dicts/`](../../dicts/README.md) (DP-9), and real fills otherwise pass
their own `--dict` (e.g. UKACD18). The bench freezes them so the
dictionary-load inference count (a measured, tracked metric) is
reproducible byte-for-byte.

## Primary dictionary — ENABLE

`enable1.txt` is the **ENABLE** word list (Enhanced North American Benchmark
Lexicon), a widely used **public-domain** Scrabble/word-game lexicon.

- Source: <https://raw.githubusercontent.com/dolph/dictionary/master/enable1.txt>
- Fetched: 2026-07-05.
- `172823` words, one per line, lowercase ASCII `a`–`z`, LF line endings.
- sha256: `3f16130220645692ed49c7134e24a18504c2ca55b3c012f7290e3e77c63b1a89`
- **License:** public domain. ENABLE was released to the public domain by its
  authors (M. Kantrowitz / A. Truscott lineage); it carries no copyright and is
  redistributed freely (it is the basis of many open word-game dictionaries).

Dictionary files are **defined as UTF-8** (fill pins `encoding(utf8)` on the
read — the process locale never affects decoding). fill's `normalize_word/2`
folds each line to A–Z letters locale-independently — ASCII is upcased and
non-letters squeezed; accented Latin letters fold through an explicit,
Unicode-derived table; words with unfoldable letters are dropped with an
unconditional stderr count — then `sort/2` dedupes (see `load_dict/3` in
fill.pl and docs/plans/fill-dict-unicode-normalization.md for the policy).
ENABLE is already lowercase ASCII, de-duplicated, and sorted, so the load is a
clean 1:1 (172823 distinct words, zero drops, no folding path taken).

## Frozen subsets — a size-scaling knob

`enable_10k.txt`, `enable_25k.txt`, `enable_50k.txt` are frozen subsets used as
ladder rungs at smaller dictionary scales (the dict-load and search cost scale
with dictionary size; see the campaign plan's Phase 3 trigger).

They are **NOT** the alphabetical first-N words — that would skew the letter
distribution (the front of a sorted ENABLE is all `a…` words) and destroy
crossing-completion realism. Instead each is a **seeded shuffle** of the full
list truncated to size, then re-sorted for a diff-friendly on-disk order:

    python3 benchmarks/gen_fill_dict.py          # regenerate
    python3 benchmarks/gen_fill_dict.py --check   # verify committed == regenerated

- Deterministic: `random.Random(20260705)` (seed recorded in the script).
- Exact counts: 10000 / 25000 / 50000 words.
- sha256:
  - `enable_10k.txt`: `5d7404b571fca4e0db57acff2dd4f33aea9a984890342e28aaba49d570e51b58`
  - `enable_25k.txt`: `0b02cdaeb470f7d70a8eaaadf8d3adcb7f572b87b4276d620de68fbea7c556b8`
  - `enable_50k.txt`: `50b1934d6f704bc8bcdd86475689b9e80dd1f23b58c25fb987de03ad042b50de`

The subset FILES are the source of truth (load-layer inference counts depend on
the exact bytes); `gen_fill_dict.py` documents how they were made and re-checks
them, but the bench reads the committed files, never regenerates.
