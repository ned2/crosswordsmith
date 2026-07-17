#!/usr/bin/env python3
"""Run the bounded P-R0 fixed-instance pilot through counter-free authority."""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import resource
import subprocess
import sys
from dataclasses import dataclass

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RUN = HERE / "run.pl"
MANIFEST = HERE / "search_seeds.json"
EXPECTED_COMMIT = "1bccf47917fea074a404bf467e51865ce988b8b8"
EXPECTED_BRANCH = "probe/a-r0"
BUDGET = 500_000_000
OUTER_TIMEOUT = 120
CPU_CAP_SECONDS = 5_400.0

sys.path.insert(0, str(ROOT))
from benchmarks.probe_arrange.schema import validate  # noqa: E402


@dataclass(frozen=True)
class Subject:
    path: str
    grid: int
    count: int
    fixture_seed: int
    kind: str

    @property
    def stem(self) -> str:
        return pathlib.Path(self.path).stem


CLIFFS = (
    Subject("benchmarks/probe_arrange/fixtures/cliff_09x09_18w_seed12.pl", 9, 18, 12, "cliff"),
    Subject("benchmarks/probe_arrange/fixtures/cliff_15x15_44w_seed11.pl", 15, 44, 11, "cliff"),
    Subject("benchmarks/probe_arrange/fixtures/cliff_21x21_84w_seed12.pl", 21, 84, 12, "cliff"),
    Subject("benchmarks/probe_arrange/fixtures/cliff_21x21_88w_seed11.pl", 21, 88, 11, "cliff"),
)
CONTROLS = (
    Subject("fixtures/ladder_09x09_08w.pl", 9, 8, 11, "control"),
    Subject("fixtures/ladder_15x15_12w.pl", 15, 12, 11, "control"),
    Subject("fixtures/ladder_21x21_25w.pl", 21, 25, 11, "control"),
)
CORNERS = ("topleft_across", "topright")


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip()


def check_base() -> None:
    commit = git("rev-parse", "HEAD")
    branch = git("branch", "--show-current")
    if (commit, branch) != (EXPECTED_COMMIT, EXPECTED_BRANCH):
        raise SystemExit(
            f"refusing pilot from {branch}@{commit}; expected "
            f"{EXPECTED_BRANCH}@{EXPECTED_COMMIT}"
        )


def load_rows(path: pathlib.Path) -> tuple[list[dict], dict[str, dict]]:
    if not path.exists():
        return [], {}
    rows = []
    by_id = {}
    with path.open(encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            validate(row)
            operation_id = row["operation_id"]
            if operation_id in by_id:
                raise SystemExit(f"duplicate operation_id on line {line_no}: {operation_id}")
            rows.append(row)
            by_id[operation_id] = row
    return rows, by_id


def consumed_cpu(rows: list[dict]) -> float:
    return sum(
        float(row.get("runner_child_cpu_seconds", row.get("cpu_seconds") or 0.0))
        for row in rows
    )


def child_limit(remaining: float):
    seconds = max(1, math.floor(remaining))

    def install() -> None:
        resource.setrlimit(resource.RLIMIT_CPU, (seconds, seconds))

    return install


def operation_id(subject: Subject, corner: str, seed_index: int) -> str:
    return f"pr0-{subject.kind}-{subject.stem}-{corner}-seed-{seed_index:02d}"


def command(subject: Subject, corner: str, seed_index: int, search_seed: int) -> list[str]:
    arm = "full_perturbation" if subject.kind == "cliff" else "easy_control_full_perturbation"
    return [
        "swipl", "-q", str(RUN), "--", "authority-corner", subject.path,
        str(subject.grid), str(subject.count), str(subject.fixture_seed), str(search_seed),
        str(BUDGET), corner, operation_id(subject, corner, seed_index), "0", arm,
        str(OUTER_TIMEOUT),
    ]


def append_row(path: pathlib.Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def run_one(
    subject: Subject,
    corner: str,
    seed_index: int,
    search_seed: int,
    path: pathlib.Path,
    rows: list[dict],
    by_id: dict[str, dict],
) -> None:
    op_id = operation_id(subject, corner, seed_index)
    if op_id in by_id:
        print(f"pr0 heartbeat operation={op_id} phase=resume-skip", file=sys.stderr, flush=True)
        return
    used = consumed_cpu(rows)
    remaining = CPU_CAP_SECONDS - used
    if remaining < 1.0:
        raise SystemExit(f"P-R0 CPU hard stop reached ({used:.3f}s consumed)")
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    process = subprocess.run(
        command(subject, corner, seed_index, search_seed),
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        check=False,
        timeout=OUTER_TIMEOUT + 15,
        preexec_fn=child_limit(remaining),
    )
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    child_cpu = (after.ru_utime + after.ru_stime) - (before.ru_utime + before.ru_stime)
    if process.returncode != 0:
        raise SystemExit(
            f"operation {op_id} exited {process.returncode}; child CPU {child_cpu:.3f}s"
        )
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise SystemExit(f"operation {op_id} emitted {len(lines)} JSON lines")
    row = json.loads(lines[0])
    validate(row)
    if row["operation_id"] != op_id or row["commit"] != EXPECTED_COMMIT:
        raise SystemExit(f"operation provenance mismatch: {op_id}")
    row.update(
        seed_index=seed_index,
        subject_kind=subject.kind,
        runner_child_cpu_seconds=child_cpu,
        pilot_partition="[0,16)",
        authority="unchanged_product_construct_one/7",
    )
    append_row(path, row)
    rows.append(row)
    by_id[op_id] = row
    print(
        f"pr0 heartbeat operation={op_id} phase=recorded outcome={row['outcome']} "
        f"pilot_cpu={consumed_cpu(rows):.3f}s",
        file=sys.stderr,
        flush=True,
    )


def first_eight_complete(by_id: dict[str, dict]) -> bool:
    return all(
        operation_id(subject, corner, seed_index) in by_id
        for seed_index in range(8)
        for subject in CLIFFS
        for corner in CORNERS
    )


def early_stop(by_id: dict[str, dict]) -> bool:
    if not first_eight_complete(by_id):
        return False
    rows = [
        by_id[operation_id(subject, corner, seed_index)]
        for seed_index in range(8)
        for subject in CLIFFS
        for corner in CORNERS
    ]
    return all(row["censored"] for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results",
        type=pathlib.Path,
        default=pathlib.Path("benchmarks/results/2026-07-17-p-r0-pilot.jsonl"),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    check_base()
    seeds = json.loads(MANIFEST.read_text(encoding="ascii"))["values"][:16]
    rows, by_id = load_rows(ROOT / args.results)
    schedule = [
        (subject, "topleft_across", 0, seeds[0]) for subject in CONTROLS
    ]
    for seed_index, search_seed in enumerate(seeds):
        for subject in CLIFFS:
            for corner in CORNERS:
                schedule.append((subject, corner, seed_index, search_seed))
    if args.dry_run:
        for subject, corner, seed_index, search_seed in schedule:
            print(" ".join(command(subject, corner, seed_index, search_seed)))
        return 0
    for subject, corner, seed_index, search_seed in schedule:
        if subject.kind == "cliff" and seed_index >= 8 and early_stop(by_id):
            print(
                "pr0 heartbeat phase=early-stop reason=all-first-eight-cliff-rows-censored",
                file=sys.stderr,
                flush=True,
            )
            break
        run_one(subject, corner, seed_index, search_seed, ROOT / args.results, rows, by_id)
    print(
        f"pr0 complete rows={len(rows)} pilot_cpu={consumed_cpu(rows):.3f}s "
        f"early_stop={str(early_stop(by_id)).lower()}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
