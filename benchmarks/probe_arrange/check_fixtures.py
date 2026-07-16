#!/usr/bin/env python3
"""Regenerate all campaign cliff fixtures to temp files and byte-compare."""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import subprocess
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
GENERATOR = ROOT / "benchmarks" / "gen_mesh_fixture.py"
SPECS = [
    (9, 18, 5, seed) for seed in (11, 12, 13)
] + [
    (15, 44, 6, seed) for seed in (11, 12, 13)
] + [
    (21, 84, 6, seed) for seed in (11, 12, 13)
] + [
    (21, 88, 6, seed) for seed in (11, 12, 13)
]


def name(grid: int, words: int, seed: int) -> str:
    return f"cliff_{grid:02d}x{grid:02d}_{words}w_seed{seed}.pl"


def count_words(path: pathlib.Path) -> int:
    text = path.read_text(encoding="ascii")
    return len(re.findall(r"^\s*\['[A-Z]+'\](?:,)?$", text, re.MULTILINE))


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    rows = []
    with tempfile.TemporaryDirectory(prefix="crosswordsmith-p0-corpus-") as temp:
        tempdir = pathlib.Path(temp)
        for grid, words, alphabet, seed in SPECS:
            committed = HERE / "fixtures" / name(grid, words, seed)
            generated = tempdir / committed.name
            subprocess.run(
                ["python3", str(GENERATOR), str(grid), str(words), str(alphabet),
                 "3", "4", str(seed), str(generated)],
                check=True, stdout=subprocess.DEVNULL,
            )
            actual = count_words(generated)
            committed_count = count_words(committed)
            if actual != words or committed_count != words:
                raise SystemExit(
                    f"count failure {committed.name}: generated={actual}, "
                    f"committed={committed_count}, requested={words}"
                )
            if generated.read_bytes() != committed.read_bytes():
                raise SystemExit(f"byte mismatch: {committed.name}")
            rows.append({"fixture": committed.name, "words": words,
                         "sha256": digest(committed)})
    print(json.dumps(rows, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
