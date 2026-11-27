# Cubic complement census

Research database of clique dynamics and flag-complex homotopy types of complements of finite simple cubic graphs. Companion to [pycliques-project](https://github.com/rvf0068/pycliques-project).

## Order 14 census (dataset snapshot 26 November 2026; certificate update 27 November 2026)

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

The repository now contains independently replayable proofs for **all 502 divergent and all 18 clique-null cases**: all 520 determined behavior classifications. The divergence proofs comprise 396 stage-zero coaffine stable-image certificates, 72 stage-one certificates, 33 stage-two coaffine certificates, and graph 308's separate incidence-cut retraction. The 18 null cases are checked by their point-paring certificates. The remaining 20 graphs are mathematically unresolved in the current ledger.

The four verdicts that originally had only historical provenance (150, 252, 270, 472) now all have fresh witnesses. Graph 150 was recovered by the implicit incidence search. The new [deep native search record](research/stage2-discovery-deep-20261127.json) and [implicit incidence record](research/stage2-incidence-20261127.json) document the final eleven recoveries. The older source labels remain unchanged in the ledger for provenance.

## Quick checks

Run `uv run --with networkx python scripts/validate_census.py` and `uv run --with networkx python scripts/verify_coaffine.py`. See the workbook for the package-dependent searches. Search algorithms belong in the monorepo; this repository holds data, evidence and reproduction drivers. Timeouts are inconclusive.

Future additions: individual homotopy-certificate replay, complete small-order per-ID ledgers, and a browsable census. The order-14 behavior ledger is now fully certified.
