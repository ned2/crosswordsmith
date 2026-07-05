# Portability notes (not loaded as skill content)

This skill follows the Agent Skills open standard (agentskills.io) and is
readable natively by Codex CLI, Gemini CLI, and OpenCode (OpenCode reads
`.claude/skills/` in place; Codex/Gemini also search the `.agents/skills/`
interop alias). Full assessment: `docs/research/skill-portability-harnesses.md`.

## Model-class mapping per harness (update as families evolve)

| Class | Claude Code | Codex CLI | Gemini CLI |
|---|---|---|---|
| frontier (orchestrator) | Fable/Mythos | top GPT tier | Gemini Pro tier |
| strong (experiments, probes) | highest Opus | top or one-down codex tier | Gemini Pro tier |
| balanced (research) | highest Sonnet | mid tier | Gemini Flash tier |

OpenCode: per-agent `"model": "provider/model-id"` — can even mix vendors.

## Adapter deltas when running outside Claude Code

1. **Worktrees**: no harness else isolates sub-agents in worktrees natively.
   Fold an explicit `git worktree add` + cleanup step into the brief's
   base-check item and instruct the agent to work there.
2. **Stalled-agent recovery**: no resume-with-context elsewhere. Use the
   skill's stated fallback only: re-dispatch fresh from the last good base
   with a state recap; never partially trust an interrupted run.
3. **Nesting**: Gemini subagents cannot spawn subagents; keep the campaign
   one level deep (orchestrator → workers), which the skill already assumes.
4. **Triggering**: Gemini shows a user consent prompt on skill activation;
   OpenCode invokes via an explicit `skill` tool call. Don't rely on silent
   activation.
