@AGENTS.md

# Claude Code notes

The shared agent guidance is imported above from [AGENTS.md](AGENTS.md). A few
reminders specific to working here with Claude Code:

- **Do not use Claude Code memories / any agent-private note store.** All lasting
  knowledge goes into version-controlled files (README.md, AGENTS.md, `docs/`).
  If you're tempted to "remember" something, write it to the appropriate tracked
  doc and commit it instead.
- **Long build/bench steps outrun the default Bash timeout.** The wasm build
  (`wasm/build/build-wasm.sh`) and any cmake/ninja step in `~/src/swipl-devel`
  run for minutes — launch them as background tasks, never under the default
  2-minute Bash timeout (see the Browser / WASM notes in [AGENTS.md](AGENTS.md)).
- **Never `Read` build/bench logs whole** — they can exceed the 256KB Read cap.
  Tail or grep them. To wait on a long step, use a background task or a Monitor
  until-loop, not `sleep N; tail`.
