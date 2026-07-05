# Acceptance ritual — fixed, serial, one experiment at a time

(Companion to the experiment-campaign skill. Read at every adjudication.)

## Accepting

1. **Verify on the branch first** (tests + equivalence), THEN merge — don't
   pollute the main branch to find out.
2. Full test suite on the merged result.
3. Ratchet re-record. Then **verify the baseline file actually changed as
   reported** — check the state, not the success message (a recorder bug
   once printed "recorded" while silently dropping entries).
4. Ledger entry, adapted from the agent's draft, using **COMPOSED numbers**
   measured against the current baseline — not the agent's standalone
   numbers if its base diverged from head. Composed and standalone can
   differ in either direction; composed is what ships.
5. Commit. Tool/infrastructure fixes discovered along the way go in
   separate commits from the experiment acceptance.

## Judging marginal results

Accept a marginal-but-monotonic win that composes cleanly and simplifies or
neutralizes the hot path — constant factors compound as the larger costs
around them shrink (the reference campaign's E-H5/E-H7 both grew in relative
value after later wins). Reject marginal results that add machinery or any
regression past tolerance.

## Rejecting

- Write the ledger entry ON THE MAIN BRANCH with the measured mechanism of
  failure, referencing the experiment branch commit. An unrecorded rejection
  will be re-litigated later.
- **Never merge a rejected branch to get its docs** — it contains
  ratchet-failing code. Leave the branch intact and referenced.
- If the rejection was caused by a wrong premise in YOUR brief, say so
  explicitly in the entry (whose error, what the correct fact is). The
  ledger's value is the honest trail.
- If a sound residue of the idea exists with a plausible fix for the
  blocker, one focused follow-up variant is justified — with a pre-declared
  strict bar ("zero regressions or the avenue closes").

## After either outcome

- Re-run the ratchet's reproduction check if determinism matters (repeat
  measurement should be exactly +0.00% on unchanged rungs).
- Ask: did this result contradict any standing verdict in the ledger? If
  yes, append a scoped correction now (seam vs category), not at close-out.
