"""`render` — deterministic SVG/HTML goldens, the D6 id-lint, decorations, and
the raster (PNG/PDF) path (spec §6.3, §11; S3/S6).

Raster tests import cairosvg lazily via ``pytest.importorskip`` so the base
suite runs without the optional extra; run ``uv run --extra raster pytest`` to
exercise them.
"""

import builtins
import re
import sys

import pytest

from xword import XwordError
from xword.formats import parse_board
from xword.render import raster
from xword.render.geometry import board_geometry
from xword.render.html import html_document
from xword.render.svg import master_svg

# The one mandatory 4-digit constant is the SVG namespace URI; the D6 id-lint
# (no injected uuid/urn/year) must ignore it.
SVG_NS = "http://www.w3.org/2000/svg"
_NONDET = re.compile(r"uuid|urn:|\b(19|20)\d{2}\b", re.IGNORECASE)


def _geom(fixtures, name, *, blank=False):
    return board_geometry(parse_board((fixtures / name).read_text()), blank=blank)


# --- Goldens (byte-stable, deterministic) ------------------------------------


@pytest.mark.parametrize("blank", [False, True])
def test_svg_golden(fixtures, golden, blank):
    tag = "blank" if blank else "solved"
    out = master_svg(_geom(fixtures, "bundled_17.native.json", blank=blank))
    assert out == (golden / f"render_bundled_17_{tag}.svg").read_text()


@pytest.mark.parametrize("blank", [False, True])
def test_html_golden(fixtures, golden, blank):
    tag = "blank" if blank else "solved"
    out = html_document(_geom(fixtures, "bundled_17.native.json", blank=blank))
    assert out == (golden / f"render_bundled_17_{tag}.html").read_text()


def test_render_is_deterministic(fixtures):
    geom = _geom(fixtures, "bundled_17.native.json")
    assert master_svg(geom) == master_svg(geom)
    assert html_document(geom) == html_document(geom)


# --- D6: ids are pure coordinate functions; nothing injected -----------------


def test_ids_are_coordinate_derived(fixtures):
    geom = _geom(fixtures, "bundled_17.native.json")
    svg = master_svg(geom)
    html = html_document(geom)
    assert 'id="cell-r0-c0"' in svg
    assert 'id="word-1-across"' in svg
    assert 'id="word-1-across"' in html  # HTML clue <li> carries the same id


def test_no_injected_nondeterminism(fixtures):
    geom = _geom(fixtures, "bundled_17.native.json")
    for artefact in (master_svg(geom), html_document(geom)):
        scrubbed = artefact.replace(SVG_NS, "")  # the namespace URI is the only allowed year
        assert _NONDET.search(scrubbed) is None


# --- Structural fidelity: circles, rebus, bars -------------------------------


def test_circle_and_rebus_render(fixtures):
    # sample.ipuz has a circled cell and a "PH" rebus (S4 keepsake)
    svg = master_svg(_geom(fixtures, "sample.ipuz.json"))
    assert "<circle " in svg
    assert ">PH<" in svg  # the rebus glyph, drawn as one shrunk <text>


def test_bars_render(fixtures):
    # sample_decorated.exolve has an east bar (R|) and a south bar (E_) (S5 keepsake)
    svg = master_svg(_geom(fixtures, "sample_decorated.exolve"))
    assert svg.count("<line ") == 2


def test_blank_strips_letters_keeps_numbers(fixtures):
    solved = master_svg(_geom(fixtures, "bundled_17.native.json"))
    blank = master_svg(_geom(fixtures, "bundled_17.native.json", blank=True))
    assert ">O<" in solved and ">O<" not in blank  # OMEGA... letters gone
    assert 'fill="#0066cc">1<' in blank  # numbers survive the blank view


# --- Raster (PNG/PDF): dimensions + magic, not byte-golden (§11) --------------


def test_png_magic_and_dimensions(fixtures):
    pytest.importorskip("cairosvg")
    svg = master_svg(_geom(fixtures, "bundled_17.native.json"))
    png = raster.to_png(svg)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    declared = int(re.search(r'width="(\d+)"', svg).group(1))
    ihdr_width = int.from_bytes(png[16:20], "big")  # PNG IHDR width
    assert ihdr_width == declared


def test_pdf_magic(fixtures):
    pytest.importorskip("cairosvg")
    svg = master_svg(_geom(fixtures, "bundled_17.native.json"))
    assert raster.to_pdf(svg).startswith(b"%PDF")


# --- Raster: the two distinct failure modes (§8, S6) -------------------------


def _fake_import(monkeypatch, exc):
    real = builtins.__import__

    def fake(name, *args, **kwargs):
        if name == "cairosvg":
            raise exc
        return real(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "cairosvg", raising=False)
    monkeypatch.setattr(builtins, "__import__", fake)


def test_missing_raster_extra(monkeypatch):
    _fake_import(monkeypatch, ImportError("No module named 'cairosvg'"))
    with pytest.raises(XwordError) as excinfo:
        raster.to_png("<svg/>")
    assert "xword[raster]" in str(excinfo.value)


def test_missing_system_cairo(monkeypatch):
    _fake_import(monkeypatch, OSError("cannot load library 'libcairo.so.2'"))
    with pytest.raises(XwordError) as excinfo:
        raster.to_pdf("<svg/>")
    assert "libcairo2" in str(excinfo.value)
