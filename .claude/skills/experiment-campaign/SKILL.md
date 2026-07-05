---
name: experiment-campaign
description: Run an experiment-driven, goal-oriented research campaign (e.g. performance optimization) as an orchestrator dispatching agents. Use when asked to systematically improve a measurable property of a system through hypotheses, probes, and ratcheted experiments — "push the efficiency of X", "run an optimization campaign", "experiment-driven research on Y".
---

# Experiment campaign: orchestrated, measured, ratcheted

A reusable process for improving a measurable property of a system through
dispatched experiments, distilled from the crosswordsmith arrange campaign
(2026-07: −71..84% search cost, byte-identical output, ~10 experiments/probes;
see `docs/research/arrange-perf-campaign-2026-07.md` and `docs/experiments.md`
for the worked example and its post-mortem lessons).

## Prerequisites — build these FIRST if absent (they are phase 0, not optional)

1. **A deterministic, machine-independent metric** to hill-climb (e.g.
   inference counts, instruction counts, node expansions — not wall time).
2. **An identity oracle** proving behavior didn't change (golden outputs
   byte-compared; stronger oracles held in reserve — see Verification).
3. **A ratchet**: per-workload baseline file + check that fails on regression
   past a tolerance and re-records on accepted wins, plus an append-only
   history ledger. Workloads = a cost ladder spanning ~3 orders of magnitude,
   covering the axes the algorithm behaves differently on (size, density...).
   Exclude budget-saturating / near-cliff workloads — they flip on any change.

## Roles and dispatch

- **Orchestrator (you)**: chief scientist. Deep code understanding, hypothesis
  generation, brief writing, adjudication, ledger. Do not run experiments
  inline; your context is for judgment and synthesis.
- **Experiment agents** (strong coding model, isolated worktree, one per
  experiment): implement + measure + verify + commit on a branch. Never merge
  themselves; never re-record the baseline.
- **Research agents** (cheaper model): online/literature sweeps, returning
  ranked candidates with sources. Distill results into durable docs
  immediately — research that lives only in a conversation is lost.
- **Probe agents**: measurement-only dispatches (instrumentation, envelope
  mapping, premise checks). Cheapest and highest-ROI class — a probe that
  kills a big refactor before it is built pays for the whole campaign.

Parallelism rule: research + probes may run alongside experiments, but never
run two candidate changes to the same code concurrently — serial composition
is what keeps attribution honest.

## Campaign lifecycle

1. **Understand**: read the target code fully yourself. Play the plan back to
   the user before starting (goal, roles, stop criterion).
2. **Probe before building**: no experiment is dispatched until its premise
   has a measurement behind it (where does the cost actually go? does the
   assumed structure exist?). Re-probe when a surprise lands.
3. **Experiment loop** (serial): dispatch → adjudicate → accept/reject →
   ratchet → ledger → next. Follow surprises: an outsized win means your cost
   model was wrong — spawn the sequel hypothesis AND a probe re-validating
   any earlier conclusion the surprise casts doubt on.
4. **Close-out**: an envelope/impact probe answering the user's actual goal
   in product terms, then a campaign summary as a durable doc.
5. **New domain?** If the campaign moves to a structurally different system,
   write the plan as a doc and dispatch a RED-TEAM review before executing:
   instruct it to empirically refute specific load-bearing claims (file:line
   or measured numbers), assume the author is pattern-matching from the last
   campaign, and that "a review that just agrees is a failed review".

## The experiment brief (every dispatch includes ALL of these)

1. **Base check first**: the exact commit the worktree must show, what to do
   if stale (merge from the main checkout), and a signature feature to verify
   exists. (Three stale-base incidents in one campaign; silent corruption of
   results otherwise.)
2. **Context pointers**: the ledger section and any research docs to read.
3. **Hypothesis + mechanism**: what to change and WHY it should pay.
4. **Soundness reasoning — with a verify order**: include your worked
   argument, and explicitly instruct the agent to verify it in code rather
   than trust it. Orchestrator soundness sketches have been wrong before; the
   brief must license the agent to overturn you.
5. **Known traps**: forward every hazard learned in prior experiments that
   could bite this one, verbatim.
6. **Pre-registered expectation** + explicit permission for null results:
   "expected X; a null/negative result is legitimate — report it straight, do
   not force a win." Calibrates the agent and makes surprises legible.
7. **Validation protocol**: tests green, identity oracle byte-identical
   (never regenerate goldens; stop and report any diff), ratchet check
   against the CURRENT baseline with the expected reference numbers inlined,
   no baseline re-record, commit on a named branch.
8. **Deliverable contract**: a lettered list — per-workload results table,
   design decisions, identity status, files/commit/branch, risks, verdict
   recommendation, and a draft ledger entry in house style. Decision-ready
   packets make adjudication cheap.

## Acceptance ritual (fixed, serial)

merge branch → full test suite → ratchet re-record → ledger entry (adapted
from the agent's draft, with COMPOSED numbers measured against the current
baseline, not the agent's standalone numbers if bases diverged) → commit.
Rejections get ledger entries too, with the measured mechanism — a rejection
without a recorded reason will be re-litigated later. Never merge a rejected
branch for its docs; write the entry on the main branch and reference the
branch commit.

## Verification escalates with win size

Default: identity oracle + ratchet. For an implausibly large win, escalate
before accepting: full-behavior equivalence (e.g. complete search-tree
solution counts across all modes), byte-identity on every workload, and the
non-default code paths (seeded/alternative modes). A big win should raise
suspicion first, celebration second. Also: verify tooling side effects — after
any state-changing tool run, check the state itself, not the success message
(the campaign found its own ratchet silently dropping new entries while
logging "recorded").

## Closing avenues and pausing

- Avenues close on **mechanism, not fatigue**: each closure cites a
  measurement (a hit rate, a diagnostic, a falsified premise). "We tried
  things and they stopped working" is not a closure.
- Scope verdicts precisely: "this SEAM is exhausted" is not "this CATEGORY is
  exhausted". If an earlier measurement still points somewhere, the search is
  not done — the biggest win of the reference campaign came after two
  premature "mined out" declarations, from exactly where the first probe had
  been pointing all along.
- Corrections are APPENDED to the ledger, never rewritten into history.
- The pause is legitimate when every remaining avenue ends in either a
  measured dead-end or a named revisit trigger.

## Escalate policy, decide engineering

Decide reversible engineering questions yourself. Escalate to the user only
decisions that encode product values (e.g. relaxing the ratchet's gate to
admit an asymmetric win). When blocked on such a decision, PARK the work
shippable-on-a-branch with a written revisit trigger in the ledger — don't
force it either way.

## Durable artefacts (as you go, not at the end)

- `docs/experiments.md` (or equivalent): one entry per accept/reject/probe,
  at decision time.
- `docs/research/`: distilled research sweeps; campaign summary at close-out
  (interventions, intuitions, impact, lessons).
- Plans for new campaigns: `docs/plans/`, red-teamed before execution.
