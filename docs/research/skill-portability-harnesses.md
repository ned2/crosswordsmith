# Research: Agent Skills portability across coding-agent harnesses

Distilled from a research sweep (2026-07-05) assessing how reusable the
`.claude/skills/experiment-campaign/` skill is across OpenAI Codex CLI,
Gemini CLI, and OpenCode. Sources: official docs (developers.openai.com,
geminicli.com, opencode.ai, platform.claude.com, agentskills.io) fetched at
research time; dates from secondary coverage marked as approximate.

## The standard

- Anthropic shipped Agent Skills (SKILL.md folders, YAML `name`+`description`
  frontmatter, progressive disclosure, optional bundled files) Oct 2025, then
  released it as an OPEN STANDARD Dec 2025 — spec at agentskills.io,
  community governance via github.com/agentskills/agentskills.
- Adoption as of mid-2026: ~40 listed clients including Codex CLI/ChatGPT,
  Gemini CLI, VS Code/Copilot, OpenCode, Cursor, Amp, Goose, JetBrains
  Junie, Kiro, Factory, Roo Code.
- Interop path convention: `.agents/skills/` is the cross-tool alias
  (Codex and Gemini both search it). **OpenCode reads `.claude/skills/*`
  directly in place** — our skills work there with zero copying.

## Per-harness (as of 2026-07)

| | Codex CLI | Gemini CLI | OpenCode |
|---|---|---|---|
| Always-on context | AGENTS.md (hierarchical) | GEMINI.md | AGENTS.md |
| SKILL.md support | native (since Dec 2025), implicit trigger | native, but activation = `activate_skill` tool + USER CONSENT prompt | native; reads `.claude/skills/` too; invocation via explicit `skill` tool call |
| Subagents | native (~Mar 2026), parallel (max_threads default 6, max_depth 1) | native (~Apr 2026), parallel, NO nesting (subagents can't spawn subagents) | Task tool + `.opencode/agent/*.md` definitions |
| Per-dispatch model | yes (`model`, `model_reasoning_effort` per agent TOML) | yes (`model` per agent file, default inherit) | yes, most flexible: any `provider/model-id`, cross-vendor |
| Worktree isolation per agent | NOT in CLI (app-only; open issue openai/codex#23095) | absent (worktrees exist but not integrated with subagents) | core absent; third-party plugins only |
| Resume/steer a running or failed subagent | partial (steer/stop, `/agent` switch; no resume-with-context) | absent (single final-summary handoff only) | absent (session navigation is UI, not re-dispatch) |

## The verdict for our experiment-campaign skill

- **Format: fully portable.** All three read SKILL.md natively; our
  frontmatter uses only universally recognized fields; unknown fields are
  ignored where that matters (OpenCode).
- **Process content: fully portable.** The scientific backbone (phase-0
  substrate, probe-before-build, serial composition, brief checklist,
  acceptance ritual, mechanism-based closure) assumes only "some way to
  dispatch a separate context and get a report back" — all three have that.
- **Capability assumptions: degrade gracefully, two gaps everywhere.**
  (1) Worktree-per-agent isolation is native NOWHERE else — the adapter is
  to fold an explicit `git worktree add` step into the brief's base-check
  item (the skill already treats base setup as explicit instruction, so
  this is an addition, not a rewrite). (2) Resume-a-stalled-agent exists
  nowhere else — downgrade to "re-dispatch fresh from the last good base
  with a state recap", which the skill's Recovery guidance already names as
  the fallback.
- **Model classes: the vendor-neutral frontier/strong/balanced framing pays
  off** — every harness has per-agent model override, so porting is a
  mapping-table entry, not a rewrite (see the skill's notes.md).
- **Trigger friction differs**: silent auto-trigger (Claude Code, Codex) vs
  consent prompt (Gemini) vs agent-decided tool call (OpenCode). No content
  change; just don't rely on silent activation everywhere.

## Confidence

High (official docs): spec fields, search paths, subagent config fields,
Gemini consent flow, OpenCode `.claude/skills` reading. Medium (secondary):
exact subagent ship dates; adoption counts. Unresolved: whether Codex CLI
subagents will get native worktree isolation (open feature request).
