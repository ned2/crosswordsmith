"""Determinism helpers (S3): pinned-width capture and canonical JSON."""

from __future__ import annotations

import io
import json
import os
import sys
from typing import Any

from rich.console import Console

# Project-wide pinned width for every captured/golden `view` (S3): pinning is
# load-bearing — an unpinned Console follows COLUMNS and de-stabilises goldens.
VIEW_WIDTH = 200


def dump_json(obj: Any) -> str:
    """Canonical JSON: deterministic, idempotent, engine-key-order (S3)."""
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"


def capture_console(*, width: int = VIEW_WIDTH) -> Console:
    """A recording console decoupled from the ambient terminal (S3 recipe)."""
    return Console(
        width=width,
        no_color=True,
        force_terminal=False,
        record=True,
        file=io.StringIO(),
        legacy_windows=False,
    )


def make_view_console() -> Console:
    """The live `view` console: colour only on a TTY without NO_COLOR (§6.1)."""
    enabled = sys.stdout.isatty() and not os.environ.get("NO_COLOR")
    return Console(no_color=not enabled, force_terminal=enabled or None)
