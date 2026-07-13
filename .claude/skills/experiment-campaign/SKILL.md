---
name: experiment-campaign
description: Run an experiment-driven, goal-oriented research campaign (e.g. performance optimization) as an orchestrator dispatching agents — for sustained campaigns that warrant multiple dispatched experiments and a regression ratchet, NOT one-off optimizations. Use when asked to systematically improve a measurable property of a code system through hypotheses, probes, and ratcheted experiments ("run an optimization campaign over X", "systematically push Y's performance").
---

# Experiment campaign: orchestrated, measured, ratcheted

A reusable process for improving a measurable property of a system through
dispatched experiments. Provenance: the crosswordsmith arrange campaign
(2026-07: −71..84% search cost, byte-identical output, ~10 experiments + 3
probes; worked example and post-mortem in
`docs/research/arrange-perf-campaign-2026-07.md` and `docs/experiments.md`).

Companion files (read at point of use, not up front):
- `brief-template.md` — the 8-part experiment-brief checklist (read before
  EVERY dispatch).
- `acceptance-checklist.md` — the accept/reject ritual (read at every
  adjudication).

## Phase 0 — the measurement substrate (adapt to the domain; build if absent)

1. **A trackable metric.** Prefer a deterministic, machine-independent proxy
   (instruction/inference/op counts, node expansions) where one exists —
   exact metrics make every later step cheaper. Where the goal is inherently
   noisy (latency, throughput), substitute a statistical protocol: fixed
   warmup, N trials, median + characterized spread, and treat sub-noise
   "wins" as null results.
2. **An equivalence check matched to the output.** Byte-identical goldens
   where output is deterministic; property-based or tolerance equivalence
   where it is not. This is what makes it safe to move fast with agents.
3. **A ratchet**: per-workload baseline + a check that fails on regression
   past a tolerance SET ABOVE CHARACTERIZED NOISE, re-records on accepted
   wins, and appends to an over-time history ledger.
4. **A workload ladder** spanning the scales the system actually runs at,
   along its MEASURED cost axes (let the first probe tell you what drives
   cost; don't guess the axes). Every rung must be a stable ratchet target,
   not a coin flip: confirm across seeds/instances. Exclude MARGINAL
   near-saturation cases (they flip on any change) but keep cases that pass
   ROBUSTLY across seeds — those are the envelope guards that protect wins.

## Roles, model classes, and dispatch

Four agent roles. The three DISPATCHED roles are parameters of this
workflow: each has a default model class below, and the user may override
any of them when invoking the skill. The orchestrator's model is NOT a
parameter — it is whatever model the skill was invoked on, and the process
assumes a frontier-class model there; if you are running on a lower tier,
flag the mismatch to the user before starting rather than proceeding
silently. The orchestrator may adapt a dispatched role's model class
mid-campaign when evidence warrants (e.g. a run-and-tabulate probe needing
no code changes can drop a tier), but must REPORT the change to the user.

| Role | Default model class |
|---|---|
| Orchestrator / chief scientist | frontier (e.g. Fable/Mythos class) — fixed at invocation |
| Experiment runners | strong (e.g. highest available Opus) |
| Probes | strong (e.g. highest available Opus) |
| Research | balanced (e.g. highest available Sonnet) |

The class logic, for when you adapt: assign by consequence and
detectability of error. Experiments and probes default strong because their
failure modes are silent — an experiment writes code with subtle soundness
traps, and a wrong probe number steers the whole campaign with nothing
downstream to catch it. Research defaults balanced because its claims are
cross-checked against the code before being acted on.

- **Orchestrator (you)**: do not run experiments inline — your context is
  for code understanding (done yourself, not delegated), hypothesis
  generation, brief writing, adjudication, and the ledger.
- **Experiment runners** (isolated worktree, one per experiment): implement
  + measure + verify + commit on a branch. Never merge themselves; never
  re-record the baseline.
- **Probes**: measurement-only dispatches (instrumentation, envelope
  mapping, premise checks). Highest-ROI class — a probe that kills a big
  refactor before it is built pays for the whole campaign.
- **Research**: literature/online sweeps returning ranked candidates with
  sources; distill into durable docs immediately.
- **Recovery**: an agent that stalls, exhausts budget, or returns an
  unusable/contradictory packet is resumed with "re-read your own edits,
  don't assume they landed" or re-dispatched from the last good base —
  never partially trusted. A confused packet usually means the brief
  under-specified; fix the brief, not just the agent.

Parallelism rule: research + probes may run alongside experiments, but never
run two candidate changes to the same code concurrently — serial composition
is what keeps attribution honest.

**One orchestrator, always.** Exactly one active orchestrator holds the
ledger; dispatched agents return DRAFT entries and never commit to the
mainline or touch the main checkout. Answer a status question ("is the agent
still running?") with READ-ONLY checks — the agent list/attach view,
`git -C <worktree> status`, `ps` — never by fork-resuming your own session
transcript into a background agent: an inherited orchestrator context plus
auto permissions does not stay a status check, it becomes a second
orchestrator writing to the mainline (the 2026-07 fill campaign lost ~1h of
reconciliation and a ledger correction exactly this way). In the
crosswordsmith repo, `scripts/campaign-status.sh` performs the full
read-only inspection in one command. If you discover a
twin, pause it via the agent view BEFORE any reconciliation; adopt its
verified results by appended correction, and never partially trust its
uncommitted state. To keep polling rare, emit a one-line status note after
every dispatch (who is running, on what branch, since when).

**Provision worktrees at current HEAD.** Auto-isolation cuts a worktree from
the commit the SESSION started on, not current mainline — so on any campaign
long enough for the mainline to advance, an auto-cut worktree starts silently
stale. Before every dispatch, the orchestrator creates the branch at current
mainline HEAD and hands the agent a ready worktree (see brief-template item 1
for the self-heal to inline in the brief).

## Campaign lifecycle

1. **Understand + plan playback**: read the target code fully yourself, then
   play the plan back to the user (goal, roles, stop criterion, how involved
   they want to be mid-campaign) before starting.
2. **Probe before building**: no experiment is dispatched until its premise
   has a measurement behind it — and a probe doesn't only falsify premises,
   it NOMINATES the next hypotheses by showing where cost actually lives.
   Beware profiler misattribution: inclusive/cumulative frames blame
   callers; confirm the hot line by self-time or a targeted micro-bench
   before optimizing it.
3. **Experiment loop** (serial): dispatch (per `brief-template.md`) →
   adjudicate (per `acceptance-checklist.md`) → ratchet → ledger → next.
   Follow surprises: an outsized win means your cost model was wrong — spawn
   the sequel hypothesis AND a probe re-validating any earlier conclusion
   the surprise casts doubt on.
4. **Close-out**: an envelope/impact probe answering the user's actual goal
   in product terms, then a campaign summary as a durable doc.
5. **New domain?** If the campaign moves to a structurally different system,
   write the plan as a doc and dispatch a RED-TEAM review before executing:
   instruct it to empirically refute specific load-bearing claims (file:line
   or measured numbers), assume the author is pattern-matching from the last
   campaign, and that "a review that just agrees is a failed review".

## Verification escalates with win size

Default: equivalence check + ratchet. For an implausibly large win, escalate
before accepting: stronger behavioral equivalence (e.g. full search-tree
solution counts across all modes, in the reference campaign), identity on
every workload, and the non-default code paths. A big win should raise
suspicion first. Also: verify tooling side effects — after any state-changing
tool run, check the state itself, not the success message (the reference
campaign found its own ratchet silently dropping new entries while logging
"recorded").

## Closing avenues and pausing

- Avenues close on **mechanism, not fatigue**: each closure cites a
  measurement (a hit rate, a diagnostic, a falsified premise).
- Scope verdicts precisely: "this SEAM is exhausted" is not "this CATEGORY
  is exhausted". If an earlier measurement still points somewhere, the
  search is not done — the reference campaign's biggest win came after two
  premature "mined out" declarations, from exactly where the first probe
  had been pointing.
- Corrections are APPENDED to the ledger, never rewritten into history.
- Weigh campaign cost: dispatching strong-model agents is expensive; stop
  for economic reasons when the remaining candidate wins are worth less
  than the tokens to find them, and SAY that's why.
- The pause is legitimate when every remaining avenue ends in either a
  measured dead-end or a named revisit trigger.

## Escalate policy, decide engineering

Decide reversible engineering questions yourself. Escalate to the user only
decisions that encode product values (e.g. relaxing the ratchet's gate to
admit an asymmetric win). When blocked on such a decision, PARK the work
shippable-on-a-branch with a written revisit trigger in the ledger — don't
force it either way.

## Durable artefacts (as you go, not at the end)

- An experiments ledger: one entry per accept/reject/probe, at decision
  time, rejections included WITH the measured mechanism (an unrecorded
  rejection will be re-litigated later).
- Research docs: distilled sweeps; campaign summary at close-out
  (interventions, intuitions, impact, lessons).
- Plans for new campaigns: written as docs, red-teamed before execution.
