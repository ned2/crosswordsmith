export const meta = {
  name: 'prolog-audit',
  description: 'Read-only SWI-Prolog 10.0.2 code audit (Map→Find→Verify→Synthesize) with manual-grounded findings',
  whenToUse: 'Run the audit defined in docs/prolog-audit-spec.md',
  phases: [
    { title: 'Map', detail: 'inventory predicates/imports/hotspots per Tier-1 file' },
    { title: 'Find', detail: 'one agent per (lens x file) + Tier-2/3 passes; manual-cited candidates' },
    { title: 'Triage', detail: 'merge/dedup candidates, assign stable ids' },
    { title: 'Verify', detail: 'adversarial perspective-diverse (correctness/idiom/perf-risk) per finding' },
    { title: 'Critic', detail: 'completeness critic: uncovered predicates/files/lenses' },
    { title: 'GapFind', detail: 'targeted second pass on critic gaps' },
    { title: 'Synthesize', detail: 'group by file, drop unverified, emit report markdown' },
  ],
}

// ----------------------------------------------------------------------------
// Shared context — the non-negotiables, injected verbatim into every agent.
// ----------------------------------------------------------------------------
const PRE = [
  'You are one agent in a READ-ONLY audit of the SWI-Prolog crossword generator in this repo.',
  'The authoritative brief is docs/prolog-audit-spec.md (already read by the orchestrator); these are its non-negotiables:',
  '',
  'GROUND TRUTH (mandatory): The version-matched manual is docs/reference/swi-manual/ (SWI-Prolog 10.0.2 — the EXACT version this repo runs).',
  '  - Every "use X instead" / "this is the idiom" / "this is a footgun" claim MUST cite a specific file under',
  '    docs/reference/swi-manual/ (e.g. lists.md, apply.md, pairs.md, aggregate.md, solutionsequences.md, ordsets.md,',
  '    assoc.md, error.md, format.md, control.md, strings.md, sort/predsort sections). GREP that tree to confirm the',
  '    predicate/section actually exists and says what you claim. swi-manual/INDEX.md is the topic map.',
  '  - Reason about SWI 10.0.2 semantics specifically: double_quotes defaults to `string`; sort/4 dedups only with @</@>',
  '    (@=</@>= keep duplicates and are stable); get_assoc/3 FAILS (not error) on missing key; findall/3 COPIES bindings.',
  '  - A claim the manual does not support is DROPPED or marked UNVERIFIED. NEVER assert idioms from model memory.',
  '',
  'READ-ONLY: Do NOT edit, write, or reformat any source file. You may run `swipl` for quick read-only reproduction',
  '  (e.g. `swipl -g "..." -t halt`) but must not modify any file. The only file the workflow writes is the final report.',
  '',
  'BEHAVIOUR FENCE: The 66 plunit tests (tests/crossword.plt), the golden file',
  '  (tests/golden/grid_17_topleft_across.txt), and the JSON input/output contract (docs/json-input-spec.md,',
  '  docs/json-output-spec.md) define observable behaviour. Any fix that COULD change observable output is behaviour=risk-HIGH.',
  '  Solver soundness/completeness must be preserved (strategies only reorder the same search tree).',
  '',
  'PERF FENCE: Performance is governed by docs/experiments.md + the benchmarks; the portable metric is INFERENCES.',
  '  A fix that trades performance for elegance is NOT auto-recommended — mark perf=needs-benchmark (or likely-worse).',
  '',
  'DO NOT RE-LITIGATE these settled decisions (cite, do not re-open). Raise a concern only as a one-shot CHALLENGE finding',
  '  (set challenge=true) backed by strong manual evidence:',
  '  - map_list_to_pairs (NOT findall) in select_inc/quality is INTENTIONAL: findall copies terms and breaks the',
  '    ==-based remove_x. (code comments + experiments E5).',
  '  - The quality engine (quality.pl) is DELIBERATELY cut-free (spec v1b.1).',
  '  - mrv_inc_deg / degree tie-break = REJECTED (experiments I6); value ordering = REJECTED (I2);',
  '    static length-order = REJECTED (E4).',
  '  - The I5 no_word_merge maximality rule is a real fix; the assoc-list grid is the chosen structure;',
  '    sort/4 with @>= is chosen deliberately where used (dup-keeping + stable).',
  '',
  'SCOPE: Tier 1 deep = crossword.pl, quality.pl. Tier 2 = tests/crossword.plt, tests/run_tests.pl,',
  '  benchmarks/run_benchmarks.pl, benchmarks/run_matrix.pl, benchmarks/start_sensitivity.pl, benchmarks/fixtures.pl.',
  '  Tier 3 light (bugs/footguns ONLY, not style) = fixtures/*.pl. EXCLUDED: docs/reference/swi-manual/**, mercury/, scratch/.',
].join('\n')

const LENSES = [
  {
    key: 'a',
    name: 'Most appropriate stdlib predicates / builtins',
    checklist: [
      '- Hand-rolled recursion a library predicate replaces (library(lists)/library(apply): maplist/foldl/include/exclude/partition; library(pairs); library(aggregate)).',
      '- Counting idiom: mrv_count/8 does findall(t, capped(...), Ts), length(Ts, Count). Is aggregate_all(count, Goal, N) or a counting fold a better fit, and does it compose with limit/2 capping? Weigh vs perf.',
      '- The recurring strip-spaces-and-measure idiom atom_chars(W,L), delete(L, \' \', L2), length(L2, N) appears in several predicates (quality.pl ~line 100,121; crossword.pl ~366,433,475,539). Is delete/3 right, is there a cleaner SWI idiom (exclude/split_string/atomic_list_concat), and should it be ONE shared helper?',
      '- Set/assoc/dict choices: library(ordsets) ops require sorted inputs — confirm ord_memberchk/ord_intersection inputs are genuine ordsets. assoc vs pairs vs dict appropriateness.',
      '- library(solution_sequences) (distinct/2, limit/2, order_by/2) where it would replace manual machinery (crossword.pl already imports it).',
      '- Randomisation: a custom shuffle where random_permutation/2 exists? Option handling: library(option) for memberchk(flag(X),Opts) patterns.',
    ].join('\n'),
  },
  {
    key: 'b',
    name: 'Idiomatic / elegant SWI-Prolog',
    checklist: [
      '- Determinism & control: appropriate ! vs ->/; vs once/1 vs \\+; green vs red cuts; predicates that should be deterministic leaving choicepoints (or vice-versa). NB quality.pl is intentionally cut-free.',
      '- Accumulator/tail-recursion vs foldl; maplist over manual recursion; library(yall) lambdas ([X]>>Goal) where they clarify a maplist/foldl goal.',
      '- Type checking: must_be/2, is_of_type/2 at boundaries (CLI/JSON input).',
      '- Error/message idiom: verify prolog:error_message//1 (crossword.pl ~297-311) is the correct CURRENT SWI extension hook (vs prolog:message//1) — cite the manual — and that throw(error(Formal,Context)) uses standard formal terms (type_error/3, domain_error/2, existence_error/2).',
      '- Module/script structure: crossword.pl is a non-module script loaded by tests; note (do not necessarily change) cost of user-module predicates / :- discontiguous smells. First-argument indexing on multi-clause predicates.',
      '- format/2,3 directive correctness; with_output_to/2 usage (crossword.pl ~218).',
    ].join('\n'),
  },
  {
    key: 'c',
    name: 'Gotchas / footguns (SWI 10.0.2)',
    checklist: [
      '- Term-copying: findall/3 copies bindings (already caused a fixed bug; map_list_to_pairs is the deliberate fix — do NOT flag it). Look for OTHER places a copy could bite (terms later compared with == or removed with an ==-based op).',
      '- sort/4 semantics: @>= keeps dups + stable; @> dedups. Confirm every sort/2|4, keysort/2, predsort/3 uses the variant its caller needs (dup-keeping where ties matter; stability).',
      '- double_quotes: default string in SWI — any "..." literal assumed to be codes/chars, or mixed with atom_chars/atom_codes?',
      '- Negation/instantiation: \\+/once over goals with unbound vars (floundering); is/2 and arithmetic comparison on possibly-unbound terms.',
      '- Resource/stream handling: --out file path (with_output_to, partial files on failure, closing/flushing), format(user_error,...) for reports.',
      '- Unbounded search / memory (project has an OOM history): any findall/aggregate_all or non-tail recursion that can blow stack on large inputs (cross-ref docs/experiments.md op-notes).',
      '- assoc: get_assoc/3 FAILS (not errors) on missing key — confirm callers rely on the right behaviour.',
      '- Harness: call_with_time_limit/2, statistics/2, nb_setval/global state and its interaction with backtracking.',
    ].join('\n'),
  },
]

const FINDING_ITEM = {
  type: 'object',
  required: ['lens', 'file', 'lines', 'severity', 'category', 'title', 'current', 'issue', 'proposal', 'manual_cite', 'behaviour', 'perf', 'confidence', 'steelman'],
  additionalProperties: false,
  properties: {
    lens: { type: 'string', enum: ['a', 'b', 'c'] },
    file: { type: 'string' },
    lines: { type: 'string', description: 'e.g. "433-444" or "100"' },
    severity: { type: 'string', enum: ['high', 'med', 'low', 'nit'] },
    category: { type: 'string' },
    title: { type: 'string' },
    current: { type: 'string', description: 'short snippet of current code' },
    issue: { type: 'string', description: 'what is wrong/suboptimal, concretely, on SWI 10.0.2' },
    proposal: { type: 'string', description: 'concrete minimal change; code if small' },
    manual_cite: { type: 'string', description: 'swi-manual/<file>.md (predicate/section) — MUST be a real file you grepped' },
    behaviour: { type: 'string', enum: ['preserving', 'risk-HIGH', 'uncertain'] },
    perf: { type: 'string', enum: ['none', 'likely-better', 'needs-benchmark', 'likely-worse'] },
    confidence: { type: 'string', enum: ['high', 'med', 'low'] },
    steelman: { type: 'string', description: 'why the current code might actually be fine' },
    challenge: { type: 'boolean', description: 'true only if this challenges a settled decision' },
  },
}

const FIND_SCHEMA = {
  type: 'object',
  required: ['candidates'],
  additionalProperties: false,
  properties: {
    notes_cleared: { type: 'array', items: { type: 'string' }, description: 'patterns you checked and found CORRECT/idiomatic (for the checked-and-cleared list), with a manual cite' },
    candidates: { type: 'array', items: FINDING_ITEM },
  },
}

const TRIAGE_ITEM = Object.assign({}, FINDING_ITEM, {
  required: ['id'].concat(FINDING_ITEM.required),
  properties: Object.assign({ id: { type: 'string' }, merged_from: { type: 'integer' } }, FINDING_ITEM.properties),
})
const TRIAGE_SCHEMA = {
  type: 'object',
  required: ['findings'],
  additionalProperties: false,
  properties: { findings: { type: 'array', items: TRIAGE_ITEM } },
}

const VERDICT_SCHEMA = {
  type: 'object',
  required: ['verdict', 'manual_supported', 'behaviour', 'perf', 'suggested_severity', 'confidence', 'reason'],
  additionalProperties: false,
  properties: {
    verdict: { type: 'string', enum: ['real', 'not_real', 'downgrade'] },
    manual_supported: { type: 'boolean', description: 'did you grep the cited manual file and confirm it supports the claim?' },
    manual_cite_check: { type: 'string', description: 'what the cited manual file actually says, or the correct cite if the original was wrong' },
    behaviour: { type: 'string', enum: ['preserving', 'risk-HIGH', 'uncertain'] },
    perf: { type: 'string', enum: ['none', 'likely-better', 'needs-benchmark', 'likely-worse'] },
    suggested_severity: { type: 'string', enum: ['high', 'med', 'low', 'nit', 'drop'] },
    confidence: { type: 'string', enum: ['high', 'med', 'low'] },
    reason: { type: 'string' },
    settled_decision_conflict: { type: 'boolean', description: 'true if this re-litigates a settled decision without strong manual evidence' },
  },
}

const CRITIC_SCHEMA = {
  type: 'object',
  required: ['coverage_ok', 'gaps'],
  additionalProperties: false,
  properties: {
    coverage_ok: { type: 'boolean' },
    summary: { type: 'string' },
    gaps: {
      type: 'array',
      items: {
        type: 'object',
        required: ['kind', 'target', 'why'],
        additionalProperties: false,
        properties: {
          kind: { type: 'string', enum: ['predicate', 'file', 'lens', 'pattern'] },
          target: { type: 'string' },
          file: { type: 'string' },
          lens: { type: 'string', enum: ['a', 'b', 'c'] },
          why: { type: 'string' },
          suggested_find: { type: 'string' },
        },
      },
    },
  },
}

// ----------------------------------------------------------------------------
// Prompt builders
// ----------------------------------------------------------------------------
function mapPrompt(file) {
  return [
    PRE, '',
    '=== TASK: Phase 0 MAP of ' + file + ' ===',
    'Read ' + file + ' in full. Produce a structured inventory to seed the finder fan-out:',
    '  - predicates: name/arity, line range, one-line role, and observed determinism (det/semidet/nondet) where clear.',
    '  - imports: every use_module / ensure_loaded.',
    '  - recurring_patterns: idioms repeated across predicates (e.g. the strip-spaces idiom, findall+sort scoring, keysort-by-pairs).',
    '  - hotspots: the highest-value audit targets (name, location, why) — at minimum check mrv_count, select_inc/select_mrv,',
    '    the strip-spaces idiom, the error_message//1 hook, --out / with_output_to handling, any aggregate/findall over large input.',
    'Do NOT propose fixes here — just map. Be precise about line numbers.',
  ].join('\n')
}

function finderPrompt(file, lens, mapSummary, tier) {
  return [
    PRE, '',
    '=== TASK: Phase 1 FIND — lens (' + lens.key + ') "' + lens.name + '" on ' + file + ' (' + tier + ') ===',
    'Audit ' + file + ' through THIS lens only. Read the file. For every candidate finding, GREP docs/reference/swi-manual/',
    'to back the claim and put the real file in manual_cite. Apply the behaviour & perf fences. Fill the full finding schema',
    '(current/issue/proposal/manual_cite/behaviour/perf/confidence/steelman). Always write a genuine steelman.',
    '',
    'Lens checklist (starting points, not exhaustive):',
    lens.checklist,
    '',
    'Map inventory for context:',
    mapSummary,
    '',
    'Rules: no model-memory idioms (drop or mark UNVERIFIED if the manual does not support it). Respect settled decisions',
    '(only challenge=true with strong manual evidence). Also populate notes_cleared with patterns you checked and found',
    'genuinely correct/idiomatic (with a manual cite) so the report can list them as checked-and-cleared.',
    'Quality over quantity, but be thorough — this lens/file is yours to cover completely.',
  ].join('\n')
}

function tier2Prompt(files, mapNote) {
  return [
    PRE, '',
    '=== TASK: Phase 1 FIND — Tier-2 review of ' + files.join(', ') + ' ===',
    'Review these test/benchmark harness files across all three lenses, with emphasis on (c) footguns and harness gotchas:',
    'call_with_time_limit/2, statistics/2 (inferences metric), nb_setval/b_setval global state vs backtracking, set_prolog_flag,',
    'plunit idioms (assertion/1, forall, setup/cleanup), with_output_to, stream handling, and any best-fit stdlib predicate',
    '(maplist/foldl/aggregate_all) that would replace hand-rolled harness loops. GREP docs/reference/swi-manual/ for every cite',
    '(plunit.md, statistics.md, builtin-statistics.md, time/thread libs, apply.md). Behaviour fence: tests must stay green;',
    'a benchmark change that alters the inferences metric is perf-relevant. Fill the full schema + notes_cleared.',
    mapNote,
  ].join('\n')
}

function tier3Prompt(files) {
  return [
    PRE, '',
    '=== TASK: Phase 1 FIND — Tier-3 LIGHT (bugs/footguns ONLY, not style) of fixture data files ===',
    'Files: ' + files.join(', ') + '. These are data/fixtures. Report ONLY genuine bugs or footguns (malformed terms, wrong',
    'arity, double_quotes surprises, data that could mislead the solver/benchmarks), NOT style. Most likely outcome is',
    '"clean" — if so, say so in notes_cleared with what you checked. GREP the manual for any cite. Fill the schema for any real bug.',
  ].join('\n')
}

function triagePrompt(candidatesJson) {
  return [
    PRE, '',
    '=== TASK: Phase 1.5 TRIAGE / DEDUP ===',
    'Below are raw candidate findings from all finders. Merge true duplicates (same file + overlapping lines + same root issue)',
    'into one record, keeping the strongest formulation and the best manual_cite; set merged_from to how many candidates fed it.',
    'Assign a stable id F001, F002, ... (zero-padded, in file then line order). Do NOT drop anything on quality grounds yet',
    '(verification is a separate next step) — only collapse genuine duplicates. Preserve challenge=true findings. Keep every',
    'distinct issue. Output the deduped findings array.',
    '',
    'RAW CANDIDATES (JSON):',
    candidatesJson,
  ].join('\n')
}

function verifyPrompt(finding, vlens) {
  return [
    PRE, '',
    '=== TASK: Phase 2 VERIFY (adversarial) — lens: ' + vlens.key + ' ===',
    'You are an independent SKEPTIC. Default to refuting unless the evidence is solid. Verify ONE candidate finding.',
    'Your verification lens: ' + vlens.focus,
    '',
    'Apply the bar: (1) issue is REAL on SWI 10.0.2 (manual-checked; reproduce with swipl if feasible);',
    '(2) the fix is correct AND behaviour-preserving (tests + golden + JSON contract stay intact);',
    '(3) it is genuinely more appropriate/idiomatic per the MANUAL, not merely different;',
    '(4) skeptic check: is the current code intentional (indexing, determinism, a documented/settled decision)?',
    'If it fails any check, vote not_real or downgrade. If the manual_cite does not actually support the claim, set',
    'manual_supported=false and give the correct cite (or none). If it re-litigates a settled decision without strong',
    'manual evidence, set settled_decision_conflict=true and verdict=not_real.',
    '',
    'CANDIDATE FINDING (JSON):',
    JSON.stringify(finding),
  ].join('\n')
}

function criticPrompt(mapSummary, coveredJson) {
  return [
    PRE, '',
    '=== TASK: COMPLETENESS CRITIC ===',
    'Given the Tier-1 map inventory and the set of findings produced so far, identify what the audit may have MISSED:',
    'predicates never examined, files/lenses under-covered, or recurring patterns not yet checked. Be concrete and point to',
    'specific predicates/lines. Only flag a gap if a real audit-worthy target was not addressed by any finding. For each gap',
    'give a suggested_find prompt focus. If coverage is solid, set coverage_ok=true with an empty/short gaps list.',
    '',
    'MAP INVENTORY:',
    mapSummary,
    '',
    'FINDINGS SO FAR (id/file/lines/lens/title):',
    coveredJson,
  ].join('\n')
}

function gapFindPrompt(gap, mapSummary) {
  return [
    PRE, '',
    '=== TASK: GAP-FIND (targeted second pass) ===',
    'The completeness critic flagged this uncovered target: ' + JSON.stringify(gap),
    'Audit specifically this target (' + (gap.file || '') + ' ' + (gap.target || '') + ', lens ' + (gap.lens || 'any') + ').',
    'Same rules: manual-cited candidates only, behaviour/perf fences, respect settled decisions, full schema + notes_cleared.',
    'If on inspection there is nothing to report, return an empty candidates array and say so in notes_cleared.',
    '',
    'Map inventory for context:',
    mapSummary,
  ].join('\n')
}

function synthPrompt(surviving, cleared, mapSummary, clearedNotes, criticOut, dateStr) {
  return [
    PRE, '',
    '=== TASK: Phase 3 SYNTHESIZE — write the report markdown ===',
    'Produce the COMPLETE contents of docs/prolog-audit-findings.md. Output ONLY the markdown (it becomes the file verbatim).',
    'It is ' + dateStr + '. This is a read-only audit of SWI-Prolog 10.0.2 code, generated by a multi-agent workflow.',
    '',
    'Structure (follow exactly):',
    '1. Title + one-paragraph preamble (read-only audit; manual-grounded; SWI 10.0.2; how to read severity/behaviour/perf).',
    '2. ## Executive summary:',
    '   - counts by severity (high/med/low/nit) of the FINAL surviving findings;',
    '   - **Top quick wins**: Med+ findings that are behaviour=preserving AND perf=none/likely-better (the safe recommends);',
    '   - **HIGH / behaviour-risk**: anything behaviour=risk-HIGH or severity=high, called out explicitly;',
    '   - **NEEDS-BENCHMARK**: every finding with perf=needs-benchmark or likely-worse (elegance-for-speed trades), listed',
    '     separately and NOT in the recommend list, each with an expected-impact note.',
    '3. ## Findings by file (group: crossword.pl, quality.pl, Tier-2, Tier-3). Render each surviving finding in the schema:',
    '   id | lens | file:lines | severity | category, then title, current, issue, proposal, manual_cite, behaviour, perf,',
    '   confidence, steelman. Use the verifier consensus (behaviour/perf/severity may have been adjusted; note downgrades).',
    '4. ## Checked and cleared — patterns audited and found correct/idiomatic (from cleared candidates + notes_cleared),',
    '   EXPLICITLY confirming the settled decisions still sound: map_list_to_pairs (term-copy fix), the cut-free quality engine,',
    '   sort/4 @>= dup-keeping choice, the rejected experiments I6/I2/E4, the I5 no_word_merge rule, the assoc-list grid.',
    '5. ## Coverage / residual gaps — from the completeness critic (what got light or no coverage, honestly).',
    'Drop UNVERIFIED claims. Keep manual cites concrete (swi-manual/<file>.md). Be concise but complete; no invented findings.',
    '',
    'SURVIVING FINDINGS (with verifier consensus) — JSON:',
    JSON.stringify(surviving),
    '',
    'CLEARED / DROPPED CANDIDATES (for the checked-and-cleared list; include why they were dropped) — JSON:',
    JSON.stringify(cleared),
    '',
    'CLEARED NOTES from finders (patterns confirmed correct) — JSON:',
    JSON.stringify(clearedNotes),
    '',
    'COMPLETENESS CRITIC OUTPUT — JSON:',
    JSON.stringify(criticOut),
    '',
    'TIER-1 MAP INVENTORY — JSON:',
    mapSummary,
  ].join('\n')
}

// ----------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------
function aggregate(finding, votes) {
  const v = votes.filter(Boolean)
  const real = v.filter(x => x.verdict === 'real').length
  const downg = v.filter(x => x.verdict === 'downgrade').length
  const notreal = v.filter(x => x.verdict === 'not_real').length
  const settledConflict = v.some(x => x.settled_decision_conflict) && !finding.challenge
  // survives if at least 2 of 3 keep it (real or downgrade) and not a settled re-litigation
  const survived = (real + downg) >= 2 && notreal < 2 && !settledConflict
  const behaviour = v.some(x => x.behaviour === 'risk-HIGH') ? 'risk-HIGH'
    : v.some(x => x.behaviour === 'uncertain') ? 'uncertain' : 'preserving'
  const perf = v.some(x => x.perf === 'likely-worse') ? 'likely-worse'
    : v.some(x => x.perf === 'needs-benchmark') ? 'needs-benchmark'
    : v.every(x => x.perf === 'likely-better') && v.length ? 'likely-better' : finding.perf || 'none'
  const order = ['high', 'med', 'low', 'nit']
  let sev = finding.severity
  if (downg > real) { const i = order.indexOf(sev); sev = order[Math.min(i + 1, 3)] }
  return Object.assign({}, finding, {
    severity: sev, behaviour, perf,
    _survived: survived, _votes: v, _real: real, _downg: downg, _notreal: notreal,
    _settledConflict: settledConflict,
  })
}

async function verifyAll(findings, VERIFY_LENSES) {
  const out = await parallel(findings.map(f => () =>
    parallel(VERIFY_LENSES.map(vl => () =>
      agent(verifyPrompt(f, vl), { label: 'verify:' + vl.key + ':' + f.id, phase: 'Verify', schema: VERDICT_SCHEMA, effort: 'high' })
    )).then(votes => aggregate(f, votes))
  ))
  return out.filter(Boolean)
}

// ----------------------------------------------------------------------------
// Orchestration
// ----------------------------------------------------------------------------
const DATE = '2026-06-22'

const VERIFY_LENSES = [
  { key: 'correctness', focus: 'CORRECTNESS. Is the issue REAL on SWI-Prolog 10.0.2? Reproduce with a quick read-only swipl query if feasible. Does the proposed fix actually work and preserve behaviour (66 plunit tests + golden + JSON contract)? Does it keep solver soundness/completeness?' },
  { key: 'idiom', focus: 'IDIOM + MANUAL GROUND-TRUTH. GREP docs/reference/swi-manual/ to confirm the cited file/predicate exists and actually says what the finding claims. Is the proposal genuinely MORE idiomatic per the manual, not merely different? Is the SWI-10.0.2-specific reasoning (double_quotes=string, sort/4 dedup rules, assoc fail-on-missing) correct?' },
  { key: 'perf-risk', focus: 'PERF-RISK + SETTLED-DECISION SKEPTIC. Could the fix change performance (inferences)? Consult docs/experiments.md and the benchmarks; if it trades elegance for speed or could regress, set perf=needs-benchmark or likely-worse. Also: does this re-litigate map_list_to_pairs / cut-free quality / I6 / I2 / E4 / I5 / sort@>= without strong manual evidence? If so set settled_decision_conflict=true.' },
]

// Phase 0 — Map
phase('Map')
log('Mapping Tier-1 files (crossword.pl, quality.pl)')
const maps = (await parallel([
  () => agent(mapPrompt('crossword.pl'), { label: 'map:crossword.pl', phase: 'Map', schema: { type: 'object', required: ['file', 'predicates', 'imports', 'hotspots'], additionalProperties: true, properties: { file: { type: 'string' }, summary: { type: 'string' }, predicates: { type: 'array', items: { type: 'object', additionalProperties: true } }, imports: { type: 'array', items: { type: 'string' } }, recurring_patterns: { type: 'array', items: { type: 'string' } }, hotspots: { type: 'array', items: { type: 'object', additionalProperties: true } } } }, effort: 'low' }),
  () => agent(mapPrompt('quality.pl'), { label: 'map:quality.pl', phase: 'Map', schema: { type: 'object', required: ['file', 'predicates', 'imports', 'hotspots'], additionalProperties: true, properties: { file: { type: 'string' }, summary: { type: 'string' }, predicates: { type: 'array', items: { type: 'object', additionalProperties: true } }, imports: { type: 'array', items: { type: 'string' } }, recurring_patterns: { type: 'array', items: { type: 'string' } }, hotspots: { type: 'array', items: { type: 'object', additionalProperties: true } } } }, effort: 'low' }),
])).filter(Boolean)
const mapSummary = JSON.stringify(maps)

// Phase 1 — Find (fan out): lens x Tier-1 file + Tier-2 + Tier-3
phase('Find')
log('Fanning out finders: ' + (LENSES.length * 2 + 3) + ' agents')
const tier1Finders = []
for (const file of ['crossword.pl', 'quality.pl']) {
  for (const lens of LENSES) {
    tier1Finders.push(() => agent(finderPrompt(file, lens, mapSummary, 'Tier 1 deep'),
      { label: 'find:' + lens.key + ':' + file, phase: 'Find', schema: FIND_SCHEMA }))
  }
}
const otherFinders = [
  () => agent(tier2Prompt(['tests/crossword.plt', 'tests/run_tests.pl'], 'Map note: tests load crossword.pl as a non-module script.'),
    { label: 'find:tier2:tests', phase: 'Find', schema: FIND_SCHEMA }),
  () => agent(tier2Prompt(['benchmarks/run_benchmarks.pl', 'benchmarks/run_matrix.pl', 'benchmarks/start_sensitivity.pl', 'benchmarks/fixtures.pl'], 'Map note: benchmarks measure inferences via statistics/2.'),
    { label: 'find:tier2:bench', phase: 'Find', schema: FIND_SCHEMA }),
  () => agent(tier3Prompt(['fixtures/bundled_17_clues.pl', 'fixtures/benchmark_08_words.pl', 'fixtures/benchmark_14_words.pl', 'fixtures/benchmark_16_dense_words.pl', 'fixtures/benchmark_20_words.pl', 'fixtures/benchmark_26_words.pl', 'fixtures/benchmark_70_mesh_words.pl', 'fixtures/quality_22_mesh.pl', 'fixtures/quality_61_mesh.pl', 'fixtures/toc_demo.pl']),
    { label: 'find:tier3:fixtures', phase: 'Find', schema: FIND_SCHEMA }),
]
const findResults = (await parallel(tier1Finders.concat(otherFinders))).filter(Boolean)
let rawCandidates = []
let clearedNotes = []
for (const r of findResults) {
  if (r && Array.isArray(r.candidates)) rawCandidates = rawCandidates.concat(r.candidates)
  if (r && Array.isArray(r.notes_cleared)) clearedNotes = clearedNotes.concat(r.notes_cleared)
}
log(rawCandidates.length + ' raw candidates, ' + clearedNotes.length + ' cleared-notes')

// Phase 1.5 — Triage / dedup
phase('Triage')
let findings = []
if (rawCandidates.length) {
  const triaged = await agent(triagePrompt(JSON.stringify(rawCandidates)), { label: 'triage', phase: 'Triage', schema: TRIAGE_SCHEMA })
  findings = (triaged && triaged.findings) || []
}
// ensure ids
findings = findings.map((f, i) => Object.assign({ id: f.id || ('F' + String(i + 1).padStart(3, '0')) }, f))
log(findings.length + ' deduped findings -> verify')

// Phase 2 — Verify (adversarial, perspective-diverse, 3 votes each)
phase('Verify')
let verified = await verifyAll(findings, VERIFY_LENSES)

// Phase — Completeness critic
phase('Critic')
const coveredJson = JSON.stringify(verified.map(f => ({ id: f.id, file: f.file, lines: f.lines, lens: f.lens, title: f.title, survived: f._survived })))
const critic = await agent(criticPrompt(mapSummary, coveredJson), { label: 'critic', phase: 'Critic', schema: CRITIC_SCHEMA, effort: 'high' })

// Phase — Gap-find (one bounded second pass on real gaps)
if (critic && critic.coverage_ok === false && Array.isArray(critic.gaps) && critic.gaps.length) {
  phase('GapFind')
  const gaps = critic.gaps.slice(0, 8)
  log('Critic flagged ' + critic.gaps.length + ' gaps; running ' + gaps.length + ' targeted finders')
  const gapResults = (await parallel(gaps.map((g, i) => () =>
    agent(gapFindPrompt(g, mapSummary), { label: 'gapfind:' + (i + 1), phase: 'GapFind', schema: FIND_SCHEMA })
  ))).filter(Boolean)
  let gapCands = []
  for (const r of gapResults) {
    if (r && Array.isArray(r.candidates)) gapCands = gapCands.concat(r.candidates)
    if (r && Array.isArray(r.notes_cleared)) clearedNotes = clearedNotes.concat(r.notes_cleared)
  }
  if (gapCands.length) {
    const t2 = await agent(triagePrompt(JSON.stringify(gapCands)), { label: 'triage:gap', phase: 'GapFind', schema: TRIAGE_SCHEMA })
    let gapFindings = ((t2 && t2.findings) || []).map((f, i) => Object.assign({ id: f.id || ('G' + String(i + 1).padStart(3, '0')) }, f))
    log(gapFindings.length + ' gap findings -> verify')
    const gapVerified = await verifyAll(gapFindings, VERIFY_LENSES)
    verified = verified.concat(gapVerified)
  }
}

// Phase 3 — Synthesize
phase('Synthesize')
const surviving = verified.filter(f => f._survived)
const cleared = verified.filter(f => !f._survived)
log('Synthesizing report: ' + surviving.length + ' surviving, ' + cleared.length + ' cleared/dropped')
const markdown = await agent(synthPrompt(surviving, cleared, mapSummary, clearedNotes, critic, DATE), { label: 'synthesize', phase: 'Synthesize', effort: 'high' })

// Stats for the orchestrator's summary
const sevCount = { high: 0, med: 0, low: 0, nit: 0 }
for (const f of surviving) if (sevCount[f.severity] !== undefined) sevCount[f.severity]++
const needsBench = surviving.filter(f => f.perf === 'needs-benchmark' || f.perf === 'likely-worse')
const highRisk = surviving.filter(f => f.severity === 'high' || f.behaviour === 'risk-HIGH')
const quickWins = surviving.filter(f => (f.severity === 'high' || f.severity === 'med') && f.behaviour === 'preserving' && (f.perf === 'none' || f.perf === 'likely-better'))

return {
  markdown,
  stats: {
    rawCandidates: rawCandidates.length,
    findings: findings.length,
    surviving: surviving.length,
    cleared: cleared.length,
    sevCount,
    needsBenchmark: needsBench.map(f => ({ id: f.id, file: f.file, lines: f.lines, title: f.title, perf: f.perf })),
    highOrRisk: highRisk.map(f => ({ id: f.id, file: f.file, lines: f.lines, title: f.title, severity: f.severity, behaviour: f.behaviour })),
    quickWins: quickWins.map(f => ({ id: f.id, file: f.file, lines: f.lines, title: f.title, severity: f.severity, manual_cite: f.manual_cite })),
    criticCoverageOk: critic ? critic.coverage_ok : null,
  },
}
