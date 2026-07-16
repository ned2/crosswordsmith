#!/usr/bin/env python3
"""Print the canonical frozen seed manifest. Used only for regeneration."""
import json
from check_seeds import values

CONSTANT_HEX = "0x43524f5353574f52"  # ASCII "CROSSWOR"

print(json.dumps({
    "algorithm": "splitmix64-product-compatible",
    "constant_name": "CROSSWOR_ASCII_U64",
    "constant_hex": CONSTANT_HEX,
    "partitions": {"pilot": [0, 16], "tuning": [16, 32], "held_out": [32, 64]},
    "values": values(int(CONSTANT_HEX, 16)),
}, indent=2, sort_keys=True))
