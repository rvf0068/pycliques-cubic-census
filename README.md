# Cubic complement census

Research database of clique dynamics and flag-complex homotopy types of complements of finite simple cubic graphs. Companion to [pycliques-project](https://github.com/rvf0068/pycliques-project).

## Order 14 snapshot (26 November 2026)

| Scope | Number |
|---|---:|
| Cubic graphs | 540 |
| Connected cubic graphs | 509 |
| Disconnected cubic graphs | 31 |
| Recorded divergent complements | 502 |
| Recorded clique-null complements | 18 |
| Unresolved | 20 |

All 540 pairs of homotopy types of G and K(G) are recorded as exact wedges of spheres; 19 pairs agree and 521 differ. These are imported research results, not fresh package runs or claims of complete witness availability.

- [Behavior ledger](census/behavior_order14.csv)
- [Homotopy transitions](census/homotopy_order14.csv)
- [Unresolved IDs](census/unresolved_order14.txt): 94, 127, 182, 185, 206, 283, 287, 320, 344, 363, 368, 441, 442, 451, 462, 473, 497, 506, 508, 509
- [Evidence and numbering](docs/provenance.md)
- [Mathematical survey](docs/cubic_complements_survey.tex)
- [Org-mode reproduction workbook](docs/cubic_complements_reproduction.org)
- [Saved graph 308 cut certificate](certificates/308.json)

The current repository has [replayable divergence certificates](certificates/README.md) for **469 of 502** recorded divergent cases: 396 at depth zero, 72 at depth one, and graph 308 at depth two. Together with the 18 point reductions, this gives 487 behavior proofs replayable without search. The other 33 recorded divergent cases still lack a bundled witness.

Four historical divergent claims (150, 252, 270, 472) lack recovered maps even in the previous research ledger. Other positive witnesses were reported replayed in the research sessions, but not all of those witnesses are bundled yet. Do not confuse the ledger with a complete certificate archive.

## Quick checks

Run `uv run --with networkx python scripts/validate_census.py` and `uv run --with networkx python scripts/verify_coaffine.py`. See the workbook for the package-dependent searches. Search algorithms belong in the monorepo; this repository holds data, evidence and reproduction drivers. Timeouts are inconclusive.

Future additions: remaining archived certificates, independent verifiers, small-order per-ID ledgers, package-only rediscovery, and a browsable census.
