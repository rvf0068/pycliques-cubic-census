"""Independently replay the finite coaffine certificates.

The witness contains a map and invariant join factors, not a search log.
Reconstruct the host from the fixed cubic catalog; no SAT or package search.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import networkx as nx

from validate_census import ROOT, catalog, fold, pared_step, read_csv, relabel


def weak_morphism(g, mapping):
    assert set(mapping) == set(g)
    assert set(mapping.values()) <= set(g)
    assert all(
        mapping[u] == mapping[v] or g.has_edge(mapping[u], mapping[v])
        for u, v in g.edges
    )


def verify_expansive(g, certificate):
    """Check coaffinity, a stable-image retraction and expansive factors."""
    n = len(g)
    assert n == certificate["host_order"] and set(g) == set(range(n))
    images = certificate["map"]
    assert len(images) == n
    assert all(type(v) is int and 0 <= v < n for v in images)
    f = dict(enumerate(images))
    weak_morphism(g, f)
    assert all(u != f[u] and not g.has_edge(u, f[u]) for u in g)

    stable = set(g)
    steps = 0
    while (image := {f[v] for v in stable}) != stable:
        assert image < stable
        stable = image
        steps += 1
    assert stable and {f[v] for v in stable} == stable
    # The map restricts to a permutation (hence an automorphism) on its
    # finite stable image. Compute an idempotent power explicitly.
    permutation = nx.DiGraph((v, f[v]) for v in stable)
    order = math.lcm(*(len(c) for c in nx.strongly_connected_components(permutation)))
    exponent = max(1, (steps + order - 1) // order) * order
    r = {v: v for v in g}
    power = f.copy()
    k = exponent
    while k:
        if k & 1:
            r = {v: power[w] for v, w in r.items()}
        power = {v: power[w] for v, w in power.items()}
        k >>= 1
    weak_morphism(g, r)
    assert set(r.values()) == stable
    assert all(r[v] == v for v in stable)

    factors = [set(a) for a in certificate["factors"]]
    assert len(factors) >= 2
    assert all(a and a <= stable and {f[v] for v in a} == a for a in factors)
    for a, b in itertools.combinations(factors, 2):
        assert a.isdisjoint(b)
        assert all(g.has_edge(u, v) for u in a for v in b)
    assert len(factors) >= 3 or any(nx.is_connected(g.subgraph(a)) for a in factors)
    return {"stable_order": len(stable), "idempotent_exponent": exponent,
            "factor_orders": [len(a) for a in factors]}


def negative_controls():
    # O_3: a valid coaffine map with three invariant joined factors.
    g = nx.complete_multipartite_graph(2, 2, 2)
    proof = {"host_order": 6, "map": [1, 0, 3, 2, 5, 4],
             "factors": [[0, 1], [2, 3], [4, 5]]}
    verify_expansive(g, proof)
    # O_2: antipodal coaffinity alone does not meet the expansive rule.
    g2 = nx.complete_multipartite_graph(2, 2)
    bad_cases = [
        (g, {**proof, "map": list(range(6))}),
        (g, {**proof, "map": [1, 0, 3, 2, 5, 6]}),
        (g, {**proof, "factors": [[0, 1], [0, 1], [4, 5]]}),
        (g, {**proof, "factors": [[0], [2, 3], [4, 5]]}),
        (g2, {"host_order": 4, "map": [1, 0, 3, 2], "factors": [[0, 1], [2, 3]]}),
    ]
    for host, witness in bad_cases:
        try:
            verify_expansive(host, witness)
        except AssertionError:
            continue
        raise AssertionError("invalid control accepted")
    # A connected invariant factor may be used instead of three factors.
    g3 = nx.disjoint_union(nx.empty_graph(2), nx.cycle_graph(5))
    g3.add_edges_from((u, v) for u in range(2) for v in range(2, 7))
    verify_expansive(g3, {"host_order": 7, "map": [1, 0, 4, 5, 6, 2, 3],
                          "factors": [[0, 1], [2, 3, 4, 5, 6]]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="write per-ID replay report")
    args = parser.parse_args()
    negative_controls()
    graphs = catalog()
    ledger = read_csv("behavior_order14.csv")
    verified = {}
    for stage in (0, 1, 2):
        path = ROOT / f"certificates/coaffine_order14_stage{stage}.json"
        data = json.loads(path.read_text())
        assert data["schema"] == 1
        certs = data["certificates"]
        expected = {i for i, row in ledger.items()
                    if row["status"] == "DIVERGENT" and row["stage"] == f"{stage}.0"}
        if stage == 2:
            # This bundle contains only rediscovered cases, excluding the
            # independently verified cut certificate 308.
            expected.remove("308")
            report_path = ROOT / "research/stage2-discovery-20261126.json"
            report = json.loads(report_path.read_text())
            assert report["completed"] and set(report["attempts"]) == expected
            assert all(a["status"] in {"FOUND", "NOT_REDISCOVERED"}
                       for a in report["attempts"].values())
            found = {i for i, a in report["attempts"].items()
                     if a["status"] == "FOUND"}
            assert set(certs) == found
            assert report["found"] == len(found)
            assert report["not_rediscovered"] == len(expected - found)
            assert all(cert["method"] == report["attempts"][i]["method"]
                       for i, cert in certs.items())
        else:
            assert set(certs) == expected
        for i, cert in certs.items():
            assert i not in verified and cert["stage"] == stage
            g = relabel(fold(nx.complement(graphs[i])))
            for _ in range(stage):
                g = pared_step(g)
            verified[i] = {"stage": stage, **verify_expansive(g, cert)}
        print(f"Stage {stage}: {len(certs)} expansive coaffine certificates verified")
    assert len(verified) == 468 + len(found)
    # The separate validate_census.py checks the point reductions and 308.
    not_in_bundles = sorted(i for i, row in ledger.items()
                            if row["status"] == "DIVERGENT" and i not in verified)
    assert "308" in not_in_bundles
    missing = [i for i in not_in_bundles if i != "308"]
    historical_missing = [i for i in missing if
                          ledger[i]["source"] == "earlier_archived_certificate_not_replayed"]
    print(f"{len(verified)} divergence witnesses independently replayed")
    print(f"Recorded divergence verdicts without bundled witnesses: {len(missing)}")
    print("Historical-only remaining:", ", ".join(historical_missing))
    if args.report:
        args.report.write_text(json.dumps(
            {"verified_coaffine": verified,
             "divergent_not_in_coaffine_bundles": not_in_bundles,
             "unbundled_divergence": missing,
             "historical_only_remaining": historical_missing}, indent=2) + "\n")


if __name__ == "__main__":
    main()
