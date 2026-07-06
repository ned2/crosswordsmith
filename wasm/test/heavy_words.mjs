// heavy_words.mjs - the HEAVY 36-word set, single-sourced from the certified
// fixture. Parsed out of fixtures/ladder_15x15_36w.pl at import time (not a
// checked-in copy), so the probes CANNOT drift from the fixture the ratchet,
// gate #2 parity, and SEARCH_INF calibration all run on: ~38.3M inferences on
// a 15x15, satisfiable by construction. Previously this list was inline-
// triplicated in headless.mjs / error_probe.mjs / yield_probe.mjs (deployment
// plan §10.4) — import { HEAVY } from './heavy_words.mjs' instead.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const FIXTURE = join(dirname(fileURLToPath(import.meta.url)),
                     '../../fixtures/ladder_15x15_36w.pl');

// The fixture is a single clues/1 fact of bare-[Answer] entries: `['DFAD'],`.
const body = readFileSync(FIXTURE, 'utf8').match(/clues\(\[([\s\S]*?)\]\)\./);
if (!body) throw new Error(`heavy_words: no clues([...]) fact in ${FIXTURE}`);
export const HEAVY = [...body[1].matchAll(/\['([A-Z]+)'\]/g)].map(m => m[1]);

// Fail loudly on fixture-format drift — a silently short list would turn the
// heavy probes into toys (the searches would no longer exercise ~38M inf).
if (HEAVY.length !== 36)
  throw new Error(`heavy_words: expected 36 words from ${FIXTURE}, got ${HEAVY.length} — fixture format changed?`);
