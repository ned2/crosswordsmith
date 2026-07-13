# The experiment brief — every dispatch includes ALL eight parts

(Companion to the experiment-campaign skill. Read before every dispatch.)

1. **Base check first**: the exact commit the worktree must show, what to do
   if stale (merge from the main checkout), and a signature feature to
   verify exists in the code. Rationale: three silently-stale worktree
   incidents in one campaign; without this, results are corrupted invisibly.
   The ORCHESTRATOR provisions the branch/worktree at current mainline HEAD
   before dispatch (auto-isolation cuts from the session-start commit, not
   HEAD) and inlines the exact self-heal in the brief with the SHA filled
   in — `git checkout -B <branch> <current-mainline-sha>` — so a stale base
   costs the agent a self-heal, not a stop-and-report round-trip.

2. **Context pointers**: the ledger section and any research/plan docs the
   agent must read before starting.

3. **Hypothesis + mechanism**: what to change and WHY it should pay, in
   terms of the measured cost structure (cite the probe).

4. **Soundness reasoning — with a verify order**: include your worked
   argument for why the change preserves behavior, and explicitly instruct
   the agent to verify it in code rather than trust it ("verify this, don't
   take it on faith"). Orchestrator soundness sketches have been wrong
   before; the brief must license the agent to overturn you — the reference
   campaign's E-H10 was rescued exactly this way.

5. **Known traps**: forward every hazard learned in prior experiments that
   could bite this one, verbatim (e.g. "the cache must be backtrack-restored
   state, never nb_/assert"; "boundary cells must stay unbound in the letter
   grid").

6. **Pre-registered expectation + null-result permission**: "expected X; a
   null or negative result is legitimate — report it straight, do not force
   a win." Calibrates the agent against manufacturing wins and makes
   surprises legible (a −50% result against a pre-registered "low single
   digits" is a signal about YOUR model).

7. **Validation protocol**: tests green; equivalence oracle passes (never
   regenerate goldens — stop and report any diff); ratchet check against the
   CURRENT baseline with the expected per-workload reference numbers inlined
   in the brief; do NOT re-record the baseline; commit on a named branch.

8. **Deliverable contract**: a lettered list — per-workload results table
   with deltas, design decisions (including variants tried and dropped, with
   numbers), equivalence/test status, files + commit + branch, risks,
   verdict recommendation, and a draft ledger entry in house style.
   Decision-ready packets make adjudication cheap.

Also state the acceptance bar in the brief when it is strict ("zero
regressions across all rungs, or the avenue closes — either outcome is
decision-grade"), so the agent knows a clean negative is a success.

Watchdog note: long, silent measurement runs get killed by the ~600s output
watchdog. Instruct the agent to emit a per-rung/per-trial heartbeat (e.g.
`format(user_error, ...)`) or chunk the run so it never goes >10 minutes
without output — two Phase-0 stalls in one campaign came from healthy runs
that just went quiet.
