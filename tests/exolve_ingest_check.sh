#!/usr/bin/env bash
#
# exolve_ingest_check.sh - automated partial check for AC-EXP-2: have the real
# Exolve engine (exolve-m.js, by Viresh Ratnakar - the same parser Exet is built
# on) INGEST a crosswordsmith Exolve export headlessly, then confirm it
# reconstructs the exact grid, entries, enumerations, and clue text of the source
# layout. This verifies the LOAD half of the Exet round-trip against the actual
# third-party engine; the full Exet UI Open->Save (the SAVE-back half) is still
# the manual step in docs/exet-verification.md (low-risk: Exet re-serialises its
# own parsed model).
#
# OPTIONAL / manual-tier - NOT part of `make test`. Needs: a Chrome/Chromium
# binary, curl + network (to fetch exolve-m.js, MIT), python3. Exits: 0 = PASS
# or SKIPPED (missing deps); 1 = a real ingestion mismatch.
#
# Usage: tests/exolve_ingest_check.sh [layout.json]   (default: the fixed golden)
set -u
cd "$(dirname "$0")/.."
LAYOUT="${1:-tests/golden/arrange_bundled_17_fixed.json}"
EXOLVE_URL="https://viresh-ratnakar.github.io/exolve-m.js"

CHROME="$(command -v google-chrome || command -v google-chrome-stable \
        || command -v chromium || command -v chromium-browser || true)"
[ -n "$CHROME" ] || { echo "SKIPPED: no Chrome/Chromium binary found"; exit 0; }
command -v python3 >/dev/null || { echo "SKIPPED: python3 not found"; exit 0; }

W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
curl -sSL -o "$W/exolve-m.js" "$EXOLVE_URL" 2>/dev/null
[ -s "$W/exolve-m.js" ] || { echo "SKIPPED: could not fetch exolve-m.js ($EXOLVE_URL)"; exit 0; }

./crosswordsmith export --to exolve "$LAYOUT" --out "$W/cs.exolve" 2>/dev/null \
    || { echo "FAIL: crosswordsmith export failed for $LAYOUT"; exit 1; }

B64="$(base64 -w0 "$W/cs.exolve")"
cat > "$W/h.html" <<HTML
<!doctype html><html><head><meta charset="utf-8"></head><body>
<div id="cs"></div><pre id="R"></pre>
<script src="exolve-m.js"></script>
<script>
try{
 const p=createExolve(atob("$B64"),"cs",false);
 const o={ok:true,w:p.gridWidth,h:p.gridHeight,grid:[],clues:[]};
 for(let r=0;r<p.gridHeight;r++){let s="";for(let c=0;c<p.gridWidth;c++){const x=p.grid[r][c];s+=(x&&x.isLight)?(x.solution||"?"):".";}o.grid.push(s);}
 for(const i in p.clues){const c=p.clues[i];if(!c||(c.dir!=='A'&&c.dir!=='D'))continue;o.clues.push({dir:c.dir,label:c.label,clue:c.clue,enumLen:c.enumLen,wordEndAfter:c.wordEndAfter,hyphenAfter:c.hyphenAfter,ncells:(c.cells||[]).length});}
 document.getElementById("R").textContent="RSTART"+JSON.stringify(o)+"REND";
}catch(e){document.getElementById("R").textContent="RSTART"+JSON.stringify({ok:false,error:String(e)})+"REND";}
</script></body></html>
HTML

"$CHROME" --headless=new --disable-gpu --no-sandbox --virtual-time-budget=8000 \
    --dump-dom "file://$W/h.html" > "$W/dump.html" 2>/dev/null
grep -oE 'RSTART.*REND' "$W/dump.html" | head -1 | sed 's/RSTART//; s/REND//' > "$W/result.json"
[ -s "$W/result.json" ] || { echo "FAIL: exolve engine produced no result (parse crash?)"; exit 1; }

RESULT_JSON="$W/result.json" LAYOUT_JSON="$LAYOUT" python3 - <<'PY'
import json, os, re
res=json.load(open(os.environ["RESULT_JSON"]))
src=json.load(open(os.environ["LAYOUT_JSON"]))
def enum_from_answer(a):
    segs=[]; seps=[]; n=0
    for ch in a:
        if ch in ' -': segs.append(n); seps.append(',' if ch==' ' else '-'); n=0
        else: n+=1
    segs.append(n)
    return "("+str(segs[0])+"".join(seps[k]+str(segs[k+1]) for k in range(len(seps)))+")"
def enum_from_exolve(L,we,hy):
    sep={**{i:',' for i in we}, **{i:'-' for i in hy}}
    segs=[]; sc=[]; st=0
    for i in range(L):
        if i in sep and i!=L-1: segs.append(i-st+1); sc.append(sep[i]); st=i+1
    segs.append(L-st)
    return "("+str(segs[0])+"".join(sc[k]+str(segs[k+1]) for k in range(len(sc)))+")"
def strip_enum(c): return re.sub(r'\s*\([0-9,\-]+\)\s*$','',c)
f=[]
if not res.get("ok"): f.append("exolve parse error: "+str(res.get("error")))
if not (res.get("w")==res.get("h")==src["gridLength"]): f.append("grid dimensions differ")
srcgrid=["".join('.' if c is None else c["letter"] for c in row) for row in src["grid"]]
if srcgrid!=res.get("grid"): f.append("grid letters/blocks differ")
exo={(c["dir"],str(c["label"])):c for c in res.get("clues",[])}
cnt={"across":0,"down":0}
for w in src["words"]:
    d='A' if w["direction"]=="across" else 'D'; cnt[w["direction"]]+=1
    c=exo.get((d,str(w["number"])))
    if not c: f.append(f'clue {d}{w["number"]} ({w["answer"]}) missing'); continue
    if enum_from_answer(w["answer"])!=enum_from_exolve(c["enumLen"],c["wordEndAfter"],c["hyphenAfter"]):
        f.append(f'{w["answer"]} enumeration differs')
    if c["ncells"]!=sum(1 for ch in w["answer"] if ch not in " -"): f.append(f'{w["answer"]} cell count differs')
    mc=(w.get("meta") or {}).get("clue","")
    if mc and strip_enum(c["clue"])!=mc: f.append(f'{w["answer"]} clue text differs')
exA=sum(1 for c in res.get("clues",[]) if c["dir"]=='A'); exD=sum(1 for c in res.get("clues",[]) if c["dir"]=='D')
if exA!=cnt["across"] or exD!=cnt["down"]: f.append("entry counts differ")
if f:
    print("EXOLVE INGEST CHECK: FAIL"); [print("  -",x) for x in f]; raise SystemExit(1)
print(f"EXOLVE INGEST CHECK: PASS - {res['w']}x{res['h']} grid, {exA} across + {exD} down, "
      f"grid+entries+enumerations+clue-text reconstructed exactly by the exolve engine.")
PY
