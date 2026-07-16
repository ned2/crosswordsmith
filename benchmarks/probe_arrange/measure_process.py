#!/usr/bin/env python3
"""Run one probe row under GNU time and append process wall/RSS metadata."""
from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent


def main() -> int:
    if not sys.argv[1:]:
        raise SystemExit("usage: measure_process.py COMMAND PROBE-ARGS...")
    with tempfile.NamedTemporaryFile(mode="r+", encoding="ascii") as timing:
        command = ["/usr/bin/time", "-f", "%e %M", "-o", timing.name,
                   "swipl", "-q", str(HERE / "run.pl"), "--", *sys.argv[1:]]
        result = subprocess.run(command, text=True, capture_output=True)
        sys.stderr.write(result.stderr)
        if result.returncode:
            raise SystemExit(result.returncode)
        row = json.loads(result.stdout)
        timing.seek(0)
        wall, rss = timing.read().split()
        row["process_wall_seconds"] = float(wall)
        row["process_max_rss_kib"] = int(rss)
        print(json.dumps(row, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
