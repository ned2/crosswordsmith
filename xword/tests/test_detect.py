"""Format-detection ladder (spec §7)."""

import pytest

from xword import XwordError
from xword.detect import detect


def test_native(fixtures):
    assert detect((fixtures / "bundled_17.native.json").read_text()) == "native"


def test_ipuz_checked_before_native(fixtures):
    assert detect((fixtures / "bundled_17.ipuz.json").read_text()) == "ipuz"
    assert detect((fixtures / "sample.ipuz.json").read_text()) == "ipuz"


def test_ipuz_signal_is_the_url_value_not_the_version_key():
    # a bare `version` key is NOT the ipuz signal
    with pytest.raises(XwordError):
        detect('{"version": "2.0", "cells": []}')
    assert detect('{"kind": ["http://ipuz.org/crossword#1"]}') == "ipuz"


def test_exolve_markers_indented(fixtures):
    assert detect((fixtures / "bundled_17.exolve").read_text()) == "exolve"
    assert detect("  \t exolve-begin\n") == "exolve"


def test_garbage_names_the_ladder():
    with pytest.raises(XwordError) as e:
        detect("hello world")
    message = str(e.value)
    for needle in ("ipuz", "native", "Exolve", "--from"):
        assert needle in message
