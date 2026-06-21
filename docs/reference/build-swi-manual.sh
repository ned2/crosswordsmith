#!/usr/bin/env bash
#
# build-swi-manual.sh - regenerate the local Markdown copy of the SWI-Prolog
# reference manual from the HTML docs shipped with the *installed* swipl, so the
# text always matches the version this repo runs against. Re-run after an
# SWI-Prolog upgrade. Requires `pandoc` and `swipl` on PATH.
#
# Output (committed, for agents to grep/read - see docs/reference/README.md):
#   docs/reference/swi-manual/*.md           core manual, one file per topic
#   docs/reference/swi-manual/packages/*.md  bundled package docs
#   docs/reference/swi-manual/INDEX.md       version header + file -> heading map
set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="$here/swi-manual"

command -v pandoc >/dev/null || { echo "error: pandoc not found" >&2; exit 1; }
command -v swipl  >/dev/null || { echo "error: swipl not found"  >&2; exit 1; }

# Source HTML dir and version straight from the installed system.
doc="$(swipl -q -g "absolute_file_name(swi(doc), D, [file_type(directory)]), writeln(D), halt" 2>/dev/null)"
ver="$(swipl --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
[ -d "$doc/Manual" ] || { echo "error: no Manual dir under $doc" >&2; exit 1; }

rm -rf "$out"
mkdir -p "$out/packages"

# HTML -> GitHub-flavoured Markdown. -native_divs/-native_spans (reader) and
# -raw_html (writer) strip the wrapper <div>/<span>/<a> HTML so only clean
# Markdown remains; the sed drops the nav image-link breadcrumb (the only .gif
# links); cat -s collapses the blank lines that leaves behind.
convert() {
    pandoc -f html-native_divs-native_spans -t gfm-raw_html --wrap=none "$1" 2>/dev/null \
        | sed -E '/\.gif\)/d' | cat -s > "$2"
}

count=0
for f in "$doc"/Manual/*.html; do
    convert "$f" "$out/$(basename "${f%.html}").md"
    count=$((count + 1))
done
for f in "$doc"/packages/*.html; do
    [ -e "$f" ] || continue
    convert "$f" "$out/packages/$(basename "${f%.html}").md"
    count=$((count + 1))
done

# INDEX.md: provenance header + filename -> first heading (its section title).
index="$out/INDEX.md"
heading() { grep -m1 -E '\S' "$1" 2>/dev/null | sed -E 's/^#+ *//; s/ +$//' | cut -c1-100; }
{
    echo "# SWI-Prolog reference manual - topic index"
    echo
    echo "Version **$ver** (matches the installed \`swipl\`). Generated from"
    echo "\`$doc/{Manual,packages}/*.html\` by \`docs/reference/build-swi-manual.sh\`."
    echo "Grep across this directory to find a predicate; open a topic file to read it."
    echo
    echo "## Core manual"
    for t in "$out"/*.md; do
        base="$(basename "$t")"
        [ "$base" = "INDEX.md" ] && continue
        echo "- [$base]($base) - $(heading "$t")"
    done
    echo
    echo "## Packages"
    for t in "$out"/packages/*.md; do
        base="$(basename "$t")"
        echo "- [packages/$base](packages/$base) - $(heading "$t")"
    done
} > "$index"

echo "converted $count files -> $out (SWI $ver)"
