"""Phase-2 `convert` (spec §6.2, §10/D7, §11).

Structural failures block and the error names the property + a capable
target; metadata drops (per-word `link`/`meta`, top-level `diagnostics`)
warn on stderr and `-q` silences them with stdout byte-identical; round-trips
are puzzle-lossless (the Q5 exolve title default is warned on the exolve hop
and, since P2, carried back through native rather than fail-stricking);
output structurally agrees with the engine's `export`
(fixtures `bundled_17.ipuz.json`/`bundled_17.exolve` are engine exports, and
`arrange_bundled_17_fixed.json` is the engine's arrange golden verbatim).
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from xword import XwordError
from xword.convert import convert_text
from xword.formats import parse_board

FIXTURES = Path(__file__).parent / "fixtures"
GOLDEN = Path(__file__).parent / "golden"
NATIVE = FIXTURES / "bundled_17.native.json"
IPUZ = FIXTURES / "bundled_17.ipuz.json"
EXOLVE = FIXTURES / "bundled_17.exolve"
ARRANGE = FIXTURES / "arrange_bundled_17_fixed.json"  # engine golden, has diagnostics


def run(*args, stdin: str | None = None):
    return subprocess.run(
        [sys.executable, "-m", "xword.cli", *args],
        input=stdin,
        capture_output=True,
        text=True,
    )


def mini_exolve(*rows: str, width: int = 3, height: int = 3) -> str:
    grid = "\n".join(f"    {row}" for row in rows)
    return (
        f"exolve-begin\n  exolve-width: {width}\n  exolve-height: {height}\n"
        f"  exolve-grid:\n{grid}\nexolve-end\n"
    )


def mini_ipuz(**overrides) -> str:
    obj = {
        "version": "http://ipuz.org/v2",
        "kind": ["http://ipuz.org/crossword#1"],
        "dimensions": {"width": 3, "height": 3},
        "puzzle": [[1, 0, 0], [0, "#", 0], [0, 0, 0]],
        "solution": [["C", "A", "B"], ["A", "#", "U"], ["B", "U", "S"]],
    }
    obj.update(overrides)
    return json.dumps(obj)


def native_without_links(src: str) -> dict:
    obj = json.loads(src)
    for word in obj["words"]:
        word["meta"].pop("link", None)
    return obj


class TestNoFeatureByteIdentity:
    """D6 acceptance gate for the serializer-gap work: `bundled_17` carries no
    rebus / prefill / shading, so its convert output must stay byte-identical
    as the three gaps land. The goldens were captured before any gap change;
    the new directives/fields must be emitted *only* when the feature is present.

    One intentional post-gap change: the exolve golden gained the Q5 default
    `exolve-title: Untitled` line (2026-07-07) — title-less boards targeting
    Exolve get the engine's default at the convert boundary (spec §14 Q5).
    """

    def test_no_feature_ipuz_byte_identical(self):
        out, _ = convert_text(NATIVE.read_text(), "ipuz")
        assert out == (GOLDEN / "convert_bundled_17_no_feature.ipuz.json").read_text()

    def test_no_feature_exolve_byte_identical(self):
        out, _ = convert_text(NATIVE.read_text(), "exolve")
        assert out == (GOLDEN / "convert_bundled_17_no_feature.exolve").read_text()


class TestStructuralFailures:
    """D7 strict side: each error names the property and a capable target."""

    def test_rectangular_blocks_native(self):
        src = mini_exolve("CAB", "OWL", height=2)
        with pytest.raises(XwordError, match="rectangular.*ipuz or exolve"):
            convert_text(src, "native")

    # NB no test_title_blocks_native: title/author are no longer structural
    # failures for native (P2 uplift). Their round-trip is pinned in
    # TestRoundTrips.test_titled_exolve_native and test_native_exolve_native.

    def test_circled_cell_blocks_native(self):
        src = mini_exolve("CAB@", "A.U", "BUS")
        with pytest.raises(XwordError, match="circled.*ipuz or exolve"):
            convert_text(src, "native")

    def test_barred_cell_blocks_native(self):
        src = mini_exolve("C|AB", "A.U", "BUS")
        with pytest.raises(XwordError, match="barred.*ipuz or exolve"):
            convert_text(src, "native")

    def test_prefilled_cell_blocks_native(self):
        src = mini_exolve("C!AB", "A.U", "BUS")
        with pytest.raises(XwordError, match="prefilled.*exolve"):
            convert_text(src, "native")

    def test_unfilled_cell_blocks_native(self):
        src = mini_exolve("0AB", "A.U", "BUS")
        with pytest.raises(XwordError, match="unfilled.*ipuz or exolve"):
            convert_text(src, "native")

    def test_foreign_styling_blocks_native_and_exolve(self):
        # `colortext` (text colour) is a real ipuz StyleSpec key with no Exolve
        # home — exolve-colour sets background only — so it stays fail-strict
        # even after Gap 3 taught the serializer the background `color` key.
        src = mini_ipuz(
            puzzle=[[1, 0, {"cell": 0, "style": {"colortext": "FF0000"}}], [0, "#", 0], [0, 0, 0]]
        )
        for target in ("native", "exolve"):
            with pytest.raises(XwordError, match="styling"):
                convert_text(src, target)

    def test_unencodable_enumeration_blocks_native(self):
        src = mini_ipuz(
            clues={"Across": [{"number": 1, "clue": "Taxi", "enumeration": "1.1.1"}], "Down": []}
        )
        with pytest.raises(XwordError, match="enumeration.*ipuz or exolve"):
            convert_text(src, "native")


class TestPrefilledToIpuz:
    """Gap 2 (closed): a prefilled/given cell serializes to ipuz as a puzzle
    `value` and survives a full serialize + parse-back round-trip."""

    def test_prefilled_serializes_value_and_round_trips(self):
        src = mini_exolve("C!AB", "A.U", "BUS")
        out, warnings = convert_text(src, "ipuz")
        assert warnings == []
        obj = json.loads(out)
        # the given cell carries its letter as a puzzle `value` (+ solution)
        assert obj["puzzle"][0][0] == {"cell": 1, "value": "C"}
        assert obj["solution"][0][0] == "C"
        # parse-back: the prefill flag + letter survive
        board = parse_board(out, "ipuz")
        given = board.grid[0][0]
        assert given.prefilled and given.letter == "C"
        others = [c for row in board.grid for c in row if c is not None and c is not given]
        assert not any(c.prefilled for c in others)  # only the given cell is prefilled

    def test_exolve_ipuz_exolve_preserves_prefill(self):
        src = mini_exolve("C!AB", "A.U", "BUS")
        via, _ = convert_text(src, "ipuz")
        back, _ = convert_text(via, "exolve")
        # the title-less source gains the Q5 default on the exolve hop — the
        # ONE sanctioned convert-time addition; everything else round-trips
        ours = parse_board(back, "exolve")
        assert ours.meta.pop("title") == "Untitled"
        assert ours == parse_board(src, "exolve")


class TestRebusToExolve:
    """Gap 1 (closed): a rebus (multi-char) cell serializes to Exolve via a
    whole-grid `exolve-option: rebus-cells` switch (space-separated tokens) and
    round-trips back through the rebus-mode tokeniser."""

    def test_rebus_serializes_and_parses_back(self, fixtures):
        src = (fixtures / "sample.ipuz.json").read_text()
        out, _ = convert_text(src, "exolve")
        assert "  exolve-option: rebus-cells" in out
        assert "    PH E F" in out  # grid switched to space-separated tokens
        board = parse_board(out, "exolve")
        assert board.grid[1][0].letter == "PH"  # rebus survives the parse-back
        assert board.grid[0][1].circle  # a decorator still rides its token (B@)

    def test_no_rebus_stays_compact(self):
        # a grid with no multi-char cell must NOT switch modes (byte-level D6)
        out, _ = convert_text(mini_exolve("CAB", "A.U", "BUS"), "exolve")
        assert "exolve-option" not in out
        assert "    CAB" in out  # compact, non-spaced row


class TestShadedToExolve:
    """Gap 3 (closed, serialize-only): a background-shade `color` style
    serializes to a deterministic exolve-colour directive. The reverse hop is
    puzzle-lossless-only — Exolve's colour model is coarser than the ipuz
    StyleSpec, so the exact StyleSpec is not reconstructed (no reverse reader)."""

    def test_color_style_serializes_exolve_colour(self, fixtures):
        src = (fixtures / "sample_shaded.ipuz.json").read_text()
        out, warnings = convert_text(src, "exolve")
        assert warnings == []  # a shade is structure, not metadata: never warns
        assert "  exolve-colour: pink r1c2" in out
        assert "exolve-option" not in out  # not rebus mode

    def test_multiple_shades_are_deterministic(self):
        src = mini_ipuz(
            puzzle=[
                [{"cell": 1, "style": {"color": "pink"}}, 0, {"cell": 0, "style": {"color": "aqua"}}],
                [0, "#", 0],
                [0, 0, 0],
            ]
        )
        out, _ = convert_text(src, "exolve")
        # one line per colour, colours sorted, cells sorted within a colour
        assert "  exolve-colour: aqua r1c3" in out
        assert "  exolve-colour: pink r1c1" in out
        assert out.index("aqua") < out.index("pink")  # deterministic colour order

    def test_shade_is_lossy_on_reverse_hop(self):
        # option (a): no reverse coordinate parser, so exolve→ipuz drops the
        # shade — the documented puzzle-lossless-only behaviour.
        src = mini_ipuz(
            puzzle=[[1, {"cell": 0, "style": {"color": "pink"}}, 0], [0, "#", 0], [0, 0, 0]]
        )
        exo, _ = convert_text(src, "exolve")
        assert "exolve-colour: pink" in exo
        back, _ = convert_text(exo, "ipuz")
        assert "color" not in back  # the StyleSpec is not reconstructed


class TestMetadataDrops:
    """D7 drop-and-warn side: native-only annotations, never structure."""

    def test_link_drop_warns_per_key(self):
        out, warnings = convert_text(NATIVE.read_text(), "ipuz")
        assert warnings == ["ipuz has no home for per-word meta key 'link'; dropped from 6 words"]
        assert "link" not in out

    def test_diagnostics_drop_warns(self):
        _, warnings = convert_text(ARRANGE.read_text(), "exolve")
        assert any("diagnostics" in w for w in warnings)
        assert any("'link'" in w for w in warnings)

    def test_native_target_drops_nothing(self):
        out, warnings = convert_text(ARRANGE.read_text(), "native")
        assert warnings == []
        got = json.loads(out)
        assert got == json.loads(ARRANGE.read_text())  # payload-lossless, diagnostics kept
        assert "diagnostics" in got


class TestRoundTrips:
    """§10: puzzle-lossless round-trips; the only loss is native's metadata."""

    def test_native_ipuz_native(self):
        src = NATIVE.read_text()
        via, _ = convert_text(src, "ipuz")
        back, warnings = convert_text(via, "native")
        assert warnings == []
        assert json.loads(back) == native_without_links(src)

    def test_native_exolve_native(self):
        # P2: native now carries title/author, so the way back SUCCEEDS. The
        # title-less source gains the Q5 default `Untitled` on the exolve hop
        # (Exet's Save requires a title); native carries that title back on the
        # return trip. It gains the ecosystem default and loses nothing
        # structural — NOT lossiness. (The exolve hop still drops per-word
        # `link` meta, exactly as ipuz does; the payload-lossless round-trip is
        # native→ipuz→native above.)
        src = NATIVE.read_text()
        via, _ = convert_text(src, "exolve")
        assert "  exolve-title: Untitled\n" in via
        back, warnings = convert_text(via, "native")
        assert warnings == []  # native holds everything the exolve hop produced
        back_obj = json.loads(back)
        assert back_obj["title"] == "Untitled"  # gained the default, carried back
        # structure survives the double hop: same grid letters + word answers
        assert parse_board(back, "native").grid == parse_board(src, "native").grid
        src_answers = {(w["number"], w["direction"]): w["answer"] for w in json.loads(src)["words"]}
        back_answers = {(w["number"], w["direction"]): w["answer"] for w in back_obj["words"]}
        assert back_answers == src_answers  # enumeration reconstruction preserved

    def test_titled_exolve_native(self):
        # A titled/authored exolve source (the engine's own export carries a
        # title) now converts to native carrying both anchors — the failure
        # this once raised is retired (P2). The engine exolve fixture is titled
        # "Untitled" and has no setter, so assert the title flows through.
        back, warnings = convert_text(EXOLVE.read_text(), "native")
        assert warnings == []
        assert json.loads(back)["title"] == "Untitled"

    def test_display_answers_reconstructed_from_enumeration(self):
        # the (5,5) hop through ipuz regains "OMEGA POINT", not "OMEGAPOINT"
        via, _ = convert_text(NATIVE.read_text(), "ipuz")
        back, _ = convert_text(via, "native")
        answers = {w["number"]: w["answer"] for w in json.loads(back)["words"] if w["direction"] == "across"}
        assert answers[1] == "OMEGA POINT"

    def test_enum_grid_mismatch_grid_wins(self):
        # a wrong ipuz enumeration is ignored; the grid letters are the answer
        src = mini_ipuz(
            clues={"Across": [{"number": 1, "clue": "Taxi", "enumeration": "(9)"}], "Down": []}
        )
        out, _ = convert_text(src, "native")
        across1 = next(w for w in json.loads(out)["words"] if w["number"] == 1 and w["direction"] == "across")
        assert across1["answer"] == "CAB"

    def test_decorations_survive_exolve_ipuz_exolve(self, fixtures):
        # bars/circles/unfilled/title all have ipuz homes (§10) — full fidelity
        src = (fixtures / "sample_decorated.exolve").read_text()
        via, warnings = convert_text(src, "ipuz")
        assert warnings == []
        assert '"barred"' in via  # synthesized StyleSpec, not a silent drop
        back, _ = convert_text(via, "exolve")
        assert parse_board(back, "exolve") == parse_board(src, "exolve")


class TestEngineCrossCheck:
    """§11: structural + anchor-field agreement with `crosswordsmith export`
    on the same layout — that structural agreement is the guarantee. After P1
    the engine retired its invented ipuz title and the constant `exolve-id`, so
    title/author now agree on both sides (ipuz: neither carries one for a
    title-less layout; exolve: both carry the Q5 `Untitled` default).

    Byte-level parity is BEST-EFFORT, not a claimed contract (spec §14). The
    Exolve header now leads with title/setter like the engine's, so in practice
    the two currently match byte-for-byte — but this suite does NOT pin that
    with a whole-file equality assertion (that would couple xword's tests to the
    engine's exact byte output and break on any engine formatting drift). ipuz
    differs only in JSON serializer formatting — SWI `json_write_dict` vs Python
    `json.dumps` — and no byte-compatibility claim is made either way. Hence the
    checks here are structural plus a few targeted, order-agnostic line checks.
    xword's own serializer determinism is pinned separately by its own goldens
    (TestNoFeatureByteIdentity)."""

    def _assert_structural_match(self, ours_text: str, engine_text: str, fmt: str):
        ours = parse_board(ours_text, fmt)
        engine = parse_board(engine_text, fmt)
        assert (ours.height, ours.width) == (engine.height, engine.width)
        assert ours.grid == engine.grid
        assert ours.words == engine.words  # numbering, clues, enumerations
        # Anchor-field parity now holds on BOTH sides (P1 retired the invented
        # ipuz title): ours and the engine agree on title AND author.
        assert ours.meta.get("title") == engine.meta.get("title")
        assert ours.meta.get("author") == engine.meta.get("author")
        if fmt == "exolve":  # title-less bundled -> both carry the Q5 default
            assert engine.meta.get("title") == "Untitled"
        else:  # ipuz: invent-nothing on both sides now (P1)
            assert "title" not in engine.meta
            assert "title" not in ours.meta

    def test_ipuz_matches_engine_export(self):
        ours, _ = convert_text(NATIVE.read_text(), "ipuz")
        self._assert_structural_match(ours, IPUZ.read_text(), "ipuz")

    def test_exolve_matches_engine_export(self):
        ours, _ = convert_text(NATIVE.read_text(), "exolve")
        self._assert_structural_match(ours, EXOLVE.read_text(), "exolve")
        # P1 retired the engine's constant exolve-id — neither side emits one now
        assert "exolve-id" not in ours
        assert "exolve-id" not in EXOLVE.read_text()
        # Header now leads with the title line as the engine does (best-effort
        # byte alignment). Order-agnostic presence check only — not a whole-file
        # equality contract (see class docstring; byte-parity is best-effort).
        assert "\n  exolve-title: Untitled\n" in ours
        assert "\n  exolve-title: Untitled\n" in EXOLVE.read_text()


class TestExolveTitleDefault:
    """Q5 (spec §14, resolved 2026-07-07): a title-less board targeting Exolve
    gains the engine's default `exolve-title: Untitled` at the convert
    boundary, with a warning — because Exet's Save crashes on a null title
    (exet-verification.md). ipuz stays invent-nothing; the library serializer
    stays faithful (the default is convert-layer only)."""

    def test_titleless_to_exolve_emits_default_and_warns(self):
        out, warnings = convert_text(mini_exolve("CAB", "A.U", "BUS"), "exolve")
        assert "\n  exolve-title: Untitled\n" in out  # byte-matches export.pl
        assert warnings == [
            "exolve requires a title (Exet's Save crashes without one); "
            "emitted default 'Untitled'"
        ]

    def test_titled_board_kept_verbatim_no_warning(self, fixtures):
        src = (fixtures / "sample_decorated.exolve").read_text()
        out, warnings = convert_text(src, "exolve")
        assert "  exolve-title: Decorated Sample" in out
        assert "Untitled" not in out
        assert not any("title" in w for w in warnings)

    def test_ipuz_stays_invent_nothing(self):
        out, warnings = convert_text(NATIVE.read_text(), "ipuz")
        assert "Untitled" not in out
        assert not any("title" in w for w in warnings)

    def test_library_serializer_stays_faithful(self):
        # invent-nothing holds at the format layer: serialize_board() adds no
        # title; only convert_text() injects the ecosystem default
        from xword.formats import serialize_board

        board = parse_board(mini_exolve("CAB", "A.U", "BUS"), "exolve")
        assert "exolve-title" not in serialize_board(board, "exolve")

    def test_cli_warns_on_stderr(self):
        result = run("convert", "--to", "exolve", str(NATIVE))
        assert result.returncode == 0
        assert "exolve requires a title" in result.stderr
        assert "  exolve-title: Untitled\n" in result.stdout


class TestCli:
    """§6.2/§8 stream contract: data on stdout, warnings on stderr, -q,
    --out only on success, detection with --from override."""

    def test_warning_on_stderr_quiet_silences_stdout_identical(self):
        loud = run("convert", "--to", "exolve", str(NATIVE))
        quiet = run("convert", "--to", "exolve", "-q", str(NATIVE))
        assert loud.returncode == 0 and quiet.returncode == 0
        assert "warning" in loud.stderr and "'link'" in loud.stderr
        assert quiet.stderr == ""
        assert loud.stdout == quiet.stdout
        assert loud.stdout.startswith("exolve-begin")

    def test_structural_failure_exits_nonzero_no_out_file(self, tmp_path):
        # A genuinely-structural failure: a rectangular source can't become
        # native (title no longer qualifies since P2 — native carries it now).
        src = tmp_path / "rect.exolve"
        src.write_text(mini_exolve("CAB", "OWL", height=2))
        out = tmp_path / "layout.json"
        result = run("convert", "--to", "native", "--out", str(out), str(src))
        assert result.returncode == 1
        assert result.stdout == ""
        assert "rectangular" in result.stderr and "ipuz or exolve" in result.stderr
        assert not out.exists()

    def test_stdin_detection_and_out(self, tmp_path):
        out = tmp_path / "puzzle.ipuz"
        result = run("convert", "--to", "ipuz", "--out", str(out), stdin=NATIVE.read_text())
        assert result.returncode == 0
        assert result.stdout == ""
        assert json.loads(out.read_text())["version"] == "http://ipuz.org/v2"

    def test_from_overrides_detection(self):
        result = run("convert", "--from", "ipuz", "--to", "exolve", str(NATIVE))
        assert result.returncode == 1
        assert "ipuz" in result.stderr

    def test_native_to_native_keeps_diagnostics(self):
        result = run("convert", "--to", "native", str(ARRANGE))
        assert result.returncode == 0
        assert result.stderr == ""
        expected = json.loads(ARRANGE.read_text())["diagnostics"]
        assert json.loads(result.stdout)["diagnostics"] == expected

    def test_diagnostics_drop_warns_on_stderr_only(self):
        result = run("convert", "--to", "ipuz", str(ARRANGE))
        assert result.returncode == 0
        assert "diagnostics" in result.stderr
        assert "diagnostics" not in result.stdout
