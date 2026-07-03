"""CLI conventions (spec §8), driven via subprocess — calling app() in-process
raises SystemExit even on success (S1)."""

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
GOLDEN = Path(__file__).parent / "golden"
NATIVE = FIXTURES / "bundled_17.native.json"


def run(*args, stdin: str | None = None, cwd: Path | None = None):
    return subprocess.run(
        [sys.executable, "-m", "xword.cli", *args],
        input=stdin,
        capture_output=True,
        text=True,
        cwd=cwd,
    )


def test_bare_xword_prints_usage_to_stderr_nonzero():
    result = run()
    assert result.returncode != 0
    assert result.stdout == ""
    assert "Usage" in result.stderr or "usage" in result.stderr.lower()


def test_view_positional_file():
    result = run("view", str(NATIVE))
    assert result.returncode == 0
    assert result.stdout.startswith("┌")
    assert result.stderr == ""


def test_view_stdin_matches_positional():
    piped = run("view", stdin=NATIVE.read_text())
    positional = run("view", str(NATIVE))
    assert piped.returncode == 0
    assert piped.stdout == positional.stdout


def test_piped_stdout_has_no_ansi():
    result = run("view", str(NATIVE))
    assert "\x1b" not in result.stdout


def test_from_overrides_detection():
    # forcing the wrong format must fail loudly, proving --from wins
    result = run("view", "--from", "ipuz", str(NATIVE))
    assert result.returncode == 1
    assert result.stdout == ""
    assert "ipuz" in result.stderr


def test_bad_format_value_is_a_parse_error():
    result = run("view", "--from", "sudoku", str(NATIVE))
    assert result.returncode != 0
    assert result.stdout == ""


def test_unrecognised_input_suggests_from():
    result = run("view", stdin="hello world")
    assert result.returncode == 1
    assert result.stdout == ""
    assert "--from" in result.stderr


def test_missing_file_is_a_parse_error(tmp_path):
    result = run("view", str(tmp_path / "nope.json"))
    assert result.returncode != 0
    assert result.stdout == ""


def test_out_written_on_success(tmp_path):
    out = tmp_path / "view.txt"
    result = run("view", str(NATIVE), "--out", str(out))
    assert result.returncode == 0
    assert result.stdout == ""
    # --out uses the pinned-width capture, so it matches the golden exactly
    assert out.read_text() == (GOLDEN / "view_bundled_17_solved.txt").read_text()


def test_out_not_written_on_failure(tmp_path):
    out = tmp_path / "view.txt"
    result = run("view", "--out", str(out), stdin="garbage input")
    assert result.returncode == 1
    assert not out.exists()


def test_blank_flag(tmp_path):
    out = tmp_path / "blank.txt"
    assert run("view", "--blank", str(NATIVE), "--out", str(out)).returncode == 0
    assert out.read_text() == (GOLDEN / "view_bundled_17_blank.txt").read_text()


def test_per_verb_help():
    result = run("view", "--help")
    assert result.returncode == 0
    assert "--blank" in result.stdout


def test_render_svg_to_stdout():
    result = run("render", "--to", "svg", str(NATIVE))
    assert result.returncode == 0
    assert result.stdout.startswith("<svg")
    assert result.stderr == ""


def test_render_out_written_on_success(tmp_path):
    out = tmp_path / "puzzle.svg"
    result = run("render", "--to", "svg", str(NATIVE), "--out", str(out))
    assert result.returncode == 0
    assert result.stdout == ""
    assert out.read_text() == (GOLDEN / "render_bundled_17_solved.svg").read_text()


def test_render_html_to_stdout():
    result = run("render", "--to", "html", str(NATIVE))
    assert result.returncode == 0
    assert result.stdout.startswith("<!doctype html>")


def test_render_png_refuses_terminal():
    # png/pdf must not spew bytes at a real terminal (§6.3); the guard fires
    # before any raster work, so this holds even without the [raster] extra.
    import os
    import pty

    master, slave = pty.openpty()
    proc = subprocess.run(
        [sys.executable, "-m", "xword.cli", "render", "--to", "png", str(NATIVE)],
        stdin=subprocess.DEVNULL,
        stdout=slave,
        stderr=subprocess.PIPE,
        text=True,
    )
    os.close(slave)
    os.close(master)
    assert proc.returncode == 1
    assert "terminal" in proc.stderr.lower()


def test_render_png_out_writes_binary(tmp_path):
    pytest.importorskip("cairosvg")
    out = tmp_path / "puzzle.png"
    result = run("render", "--to", "png", str(NATIVE), "--out", str(out))
    assert result.returncode == 0
    assert result.stdout == ""
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
