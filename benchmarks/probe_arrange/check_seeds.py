#!/usr/bin/env python3
"""Check frozen SplitMix64 seeds with Python and product Prolog implementations."""
from __future__ import annotations

import json
import pathlib
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MANIFEST = HERE / "search_seeds.json"
MASK = (1 << 64) - 1


def splitmix64(state: int) -> tuple[int, int]:
    state = (state + 0x9E3779B97F4A7C15) & MASK
    z = ((state ^ (state >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return z ^ (z >> 31), state


def values(constant: int) -> list[int]:
    out = []
    state = constant
    for _ in range(64):
        value, state = splitmix64(state)
        out.append(value)
    return out


def product_values(constant: int) -> list[int]:
    goal = (
        f"S0={constant},findall(V,(between(1,64,I),"
        "(I=:=1->S=S0;nb_getval(sm_s,S)),"
        "crosswordsmith_core:splitmix64(S,V,S1),nb_setval(sm_s,S1)),Vs),"
        "json_write_dict(current_output,_{values:Vs}),nl,halt"
    )
    result = subprocess.run(
        ["swipl", "-q", "-l", "load.pl", "-g", goal], cwd=ROOT,
        check=True, text=True, capture_output=True,
    )
    return json.loads(result.stdout)["values"]


def main() -> int:
    doc = json.loads(MANIFEST.read_text(encoding="ascii"))
    constant = int(doc["constant_hex"], 16)
    expected = doc["values"]
    if values(constant) != expected:
        raise SystemExit("manifest differs from independent Python SplitMix64")
    if product_values(constant) != expected:
        raise SystemExit("manifest differs from product splitmix64/3")
    partitions = doc["partitions"]
    if partitions != {"pilot": [0, 16], "tuning": [16, 32], "held_out": [32, 64]}:
        raise SystemExit(f"invalid partitions: {partitions}")
    print("seeds: OK (64 values; pilot=16 tuning=16 held_out=32; Python=product)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
