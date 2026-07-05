// crosswordsmith.d.ts - hand-written types for the client-side SDK
// (wasm-sdk-strategy §5.4, DEC-5). The wire shapes mirror the engine's
// envelope (browser.pl) and the arrange layout contract
// (docs/json-output-spec.md); tests/golden/arrange_*.json must structurally
// satisfy `Layout` — wasm/test/golden_type_check.ts is the drift lock.

// --- arrange input (docs/json-input-spec.md + DEC-6 params) ------------------

export interface ClueEntry {
  /** The answer to place. Letters as-placed (spaces are stripped). */
  answer: string;
  /** Optional passthrough metadata (by convention `clue` and `link`); the
   *  engine never inspects it and copies it verbatim into the output. */
  meta?: Record<string, unknown>;
}

export interface ArrangeInput {
  clues: ClueEntry[];
  /** Grid side N (an N x N grid). Positive integer; default 15 (the CLI's
   *  bare-invocation default). */
  size?: number;
  /** "fixed" (default): emit exactly N x N. "max": build on N x N, then crop
   *  to the tight enclosing square (the CLI's --max-size framing) — usually
   *  nicer for an HTML/CSS renderer than dead border rows. */
  mode?: "fixed" | "max";
  /** true: place a maximal subset, reporting drops in
   *  result.diagnostics.arrange.dropped. Default false (strict: all words or
   *  a failure envelope). Not combinable with seed. */
  bestEffort?: boolean;
  /** Perturb the strict search with RNG seed N (integer >= 0) for a
   *  reproducible pseudo-random layout; omit for the deterministic layout.
   *  Draw with randomSeed() for "regenerate". Not combinable with bestEffort.
   *  The engine owns its PRNG (portable splitmix64, not the VM RNG — which
   *  differs between the native GMP and wasm builds), so the same seed
   *  reproduces the same layout on the CLI, in the browser, and on any
   *  future build. Reproducibility is still ENGINE-VERSION-scoped (OQ-8):
   *  a future engine may change search heuristics, so record provenance
   *  (cw.version.engine) alongside any kept seed. */
  seed?: number;
}

// --- arrange result (docs/json-output-spec.md §6) ----------------------------

export interface GridCell {
  letter: string;
  /** Clue number when this is a start cell, else null. */
  number: number | null;
  /** Clue number of the across word covering this cell, else null. */
  across: number | null;
  /** Clue number of the down word covering this cell, else null. */
  down: number | null;
}

export interface PlacedWord {
  number: number;
  direction: "across" | "down";
  answer: string;
  /** [row, col] per letter, 0-based from the grid's top-left. */
  cells: [number, number][];
  /** The input entry's meta, verbatim ({} when none was given). */
  meta: Record<string, unknown>;
}

export interface ArrangeDiagnostics {
  /** true when no placed word reaches its checking target (the capped
   *  objective degenerated to plain total-crossings). */
  capInert: boolean;
  /** Input answers not placed (best-effort drops), in input order; [] when
   *  everything placed. */
  dropped: string[];
  /** The engine's objective value for this layout. */
  reward: number;
  /** Present ONLY when a seed perturbed the search: the seed that reproduces
   *  this exact layout. Absent = the deterministic search. */
  seed?: number;
}

export interface Layout {
  gridLength: number;
  /** gridLength rows x gridLength cells; null = empty cell. */
  grid: (GridCell | null)[][];
  words: PlacedWord[];
  /** Quality caveats (json-output-spec §6.4). Always present on arrange
   *  results; a consumer may ignore it wholesale and lose only provenance. */
  diagnostics: { arrange: ArrangeDiagnostics };
}

// --- lint (design-spec §8.1; strategy §6 phase 2) -----------------------------

export type LintProfile = "toc" | "blocked-uk" | "american" | "barred-ximenean";
export type LintSeverity = "PASS" | "WARN" | "FAIL";

export interface LintInput {
  /** A canonical layout — feed an arrange success `result` straight in. */
  layout: Layout;
  /** Required: the house-style rule profile (same domain as the CLI's
   *  --profile). */
  profile: LintProfile;
  /** Relax the symmetry rule to advisory (the CLI's --allow-asymmetry).
   *  Default false. */
  allowAsymmetry?: boolean;
}

export interface LintResult {
  rule: string;
  severity: LintSeverity;
  /** Present only on WARN/FAIL: what tripped, human-readable. */
  detail?: string;
}

export interface LintWordReport {
  number: number;
  direction: "across" | "down";
  answer: string;
  results: LintResult[];
}

/** The report IS the success result — a FAIL verdict is still a successful
 *  lint (the CLI's verdict-as-exit-code is a shell convention). */
export interface LintReport {
  profile: string;
  allowAsymmetry: boolean;
  /** Worst severity across all rule results. */
  verdict: LintSeverity;
  summary: { pass: number; warn: number; fail: number };
  grid: LintResult[];
  words: LintWordReport[];
}

// --- export (design-spec §8.2; strategy §6 phase 3) ---------------------------

export interface ExportInput {
  /** A canonical layout — feed an arrange success `result` straight in. */
  layout: Layout;
  /** ipuz: the result is the ipuz v2 JSON document itself (onward to
   *  .puz/.jpz/PDF via off-the-shelf tooling). exolve: the result is
   *  {format:"text", body} — Exolve plain text (round-trips to Exet). */
  to: "ipuz" | "exolve";
}

/** ipuz v2 crossword document. An external spec — typed openly on purpose;
 *  the engine emits at least these members. */
export interface IpuzDocument {
  version: string;
  kind: string[];
  [key: string]: unknown;
}

export interface ExolveText {
  format: "text";
  body: string;
}

// --- the response envelope (wasm-sdk-strategy §4) ----------------------------

export type FailureReason = "no_interlock" | "grid_too_small" | "unplaceable_words";

export interface FailureDetail {
  reason: FailureReason;
  /** Present only for unplaceable_words: answers with no possible crossing.
   *  The typical infeasible case (grid too tight for a full interlock) has
   *  no culprit list. */
  words?: string[];
}

export type CwErrorType =
  | "validation"          // bad params / clue schema (fix the request)
  | "budget_exceeded"     // search hit its inference budget; feasibility open
  | "resource_exhausted"  // the stack cap tripped (best-effort on mobile)
  | "unknown_verb"        // this engine (qlf) does not speak that verb
  | "unsupported_version" // envelope version mismatch (stale qlf or SDK)
  | "internal";           // engine bug or transport disaster — report it
// NOT "cancelled": a terminated worker posts nothing, so cancellation is the
// CancelledError rejection, never an envelope.

export interface CwError {
  type: CwErrorType;
  message: string;
  detail?: unknown;
}

export type Envelope<V extends string, R> =
  | { v: 1; id: string | null; verb: V; status: "success"; result: R }
  | { v: 1; id: string | null; verb: V; status: "failure"; detail: FailureDetail }
  | { v: 1; id: string | null; verb: V; status: "error"; error: CwError };

export type ArrangeEnvelope = Envelope<"arrange", Layout>;
export type LintEnvelope = Envelope<"lint", LintReport>;
export type ExportEnvelope<R = IpuzDocument | ExolveText> = Envelope<"export", R>;

export interface Capabilities {
  /** The verbs the LOADED engine answers for (engine-sourced, OQ-4). */
  verbs: string[];
  engine: { swipl: string };
}

// --- the facade --------------------------------------------------------------

export interface CreateOptions {
  /** Directory holding worker.js + swipl-web.{js,wasm,data} +
   *  crosswordsmith.qlf. Default: ../client/ relative to this module (the
   *  in-repo layout). */
  assetBaseUrl?: string | URL;
  /** Per-request logical stack cap in bytes (positive integer; default
   *  256_000_000). Keep it under the platform memory ceiling so a runaway
   *  search yields a resource_exhausted envelope, not a tab abort. */
  stackLimit?: number;
  /** Engine yield interval in inferences (not load-bearing; default 50_000). */
  heartbeat?: number;
}

export interface RequestOptions {
  /** Abort the request: a queued request is dequeued; an in-flight one
   *  hard-kills its worker (terminate + warm-spare promote). Either way the
   *  promise rejects with CancelledError. */
  signal?: AbortSignal;
}

/** @deprecated alias kept from the arrange-only slice; use RequestOptions. */
export type ArrangeOptions = RequestOptions;

export interface Crosswordsmith {
  /** Solve one arrange request. All three envelope statuses RESOLVE —
   *  discriminate on res.status; a legitimate "no layout" is a failure, not
   *  a rejection. Rejects only with CancelledError or a transport Error. */
  arrange(input: ArrangeInput, opts?: RequestOptions): Promise<ArrangeEnvelope>;
  /** Validate a layout against a house-style profile. Well-formed input
   *  always resolves to a success envelope carrying the report. */
  lint(input: LintInput, opts?: RequestOptions): Promise<LintEnvelope>;
  /** Transform a layout to ipuz or Exolve. The result type follows `to`. */
  export(input: ExportInput & { to: "ipuz" },
         opts?: RequestOptions): Promise<ExportEnvelope<IpuzDocument>>;
  export(input: ExportInput & { to: "exolve" },
         opts?: RequestOptions): Promise<ExportEnvelope<ExolveText>>;
  export(input: ExportInput, opts?: RequestOptions): Promise<ExportEnvelope>;
  /** Fresh engine round-trip (never a JS constant). */
  capabilities(): Promise<Capabilities>;
  readonly version: { sdk: string; engine: { swipl: string } };
  /** Terminate both workers and reject anything pending (CancelledError). */
  dispose(): void;
}

/** Spawn + warm the engine Worker (plus one standby spare) and resolve once
 *  the engine answered its first capabilities round-trip. */
export function createCrosswordsmith(opts?: CreateOptions): Promise<Crosswordsmith>;

/** An app-side seed for "regenerate": integer in [0, 1e9] (the CLI's
 *  --shuffle range — floats and >2^53 values are rejected by the engine). */
export function randomSeed(): number;

export class CancelledError extends Error {
  readonly name: "CancelledError";
}
