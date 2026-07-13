#!/usr/bin/env bash
#
# scripts/campaign-status.sh - read-only status of dispatched-agent worktrees.
#
# For each worktree under .claude/worktrees/: branch, HEAD, ahead/behind vs
# main, dirty-file count, most recent file write, and any live swipl/node
# process running inside it. This is the inspection an orchestrator needs to
# answer "is the agent still running?" WITHOUT spawning anything - it never
# mutates state, so it is always safe to run while agents are live.
#
#   scripts/campaign-status.sh

set -u
cd "$(dirname "$0")/.."

echo "main @ $(git rev-parse --short main) ($(git log -1 --format=%s main | cut -c1-60))"
echo

shopt -s nullglob
worktrees=(.claude/worktrees/*/)
if [[ ${#worktrees[@]} -eq 0 ]]; then
    echo "no worktrees under .claude/worktrees/"
    exit 0
fi

for wt in "${worktrees[@]}"; do
    name=$(basename "$wt")
    if ! git -C "$wt" rev-parse --git-dir >/dev/null 2>&1; then
        echo "$name: not a git worktree?"
        continue
    fi
    branch=$(git -C "$wt" branch --show-current 2>/dev/null)
    sha=$(git -C "$wt" rev-parse --short HEAD 2>/dev/null)
    counts=$(git -C "$wt" rev-list --left-right --count main...HEAD 2>/dev/null) || counts="? ?"
    behind=$(awk '{print $1}' <<<"$counts")
    ahead=$(awk '{print $2}' <<<"$counts")
    dirty=$(git -C "$wt" status --porcelain 2>/dev/null | wc -l | tr -d ' ')
    lastwrite=$(find "$wt" -name .git -prune -o -type f -printf '%TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null \
                | sort -r | head -1)
    echo "$name"
    echo "  branch: ${branch:-<detached>} @ ${sha:-?}  (+$ahead ahead / -$behind behind main)"
    echo "  dirty files: $dirty"
    echo "  last write:  ${lastwrite:-none}"
    pids=$(pgrep -af 'swipl|node' 2>/dev/null | grep -F "$(realpath "$wt")" || true)
    if [[ -n $pids ]]; then
        echo "  live processes:"
        sed 's/^/    /' <<<"$pids"
    else
        echo "  live processes: none"
    fi
    echo
done
