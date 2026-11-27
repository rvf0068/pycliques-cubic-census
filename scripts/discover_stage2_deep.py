"""Longer deterministic rediscovery of outstanding stage-two divergences."""
import argparse
import itertools
import json
import random

import networkx as nx
import pycliques as pc

from validate_census import catalog, fold, pared_step, relabel, read_csv
from verify_coaffine import verify_expansive


def factors(g, f):
    stable = set(g)
    while (image := {f[v] for v in stable}) != stable:
        stable = image
    orbits = list(nx.strongly_connected_components(
        nx.DiGraph((v, f[v]) for v in stable)))
    joined = nx.Graph()
    joined.add_nodes_from(range(len(orbits)))
    joined.add_edges_from(
        (i, j) for i, j in itertools.combinations(joined, 2)
        if all(g.has_edge(v, w) for v in orbits[i] for w in orbits[j])
    )
    for i, j in joined.edges:
        common = set(joined[i]) & set(joined[j])
        if common:
            return [orbits[i], orbits[j], orbits[min(common)]]
    for i in joined:
        B = set().union(*(orbits[j] for j in joined[i]))
        for C in nx.connected_components(g.subgraph(B)):
            if {f[v] for v in C} == C:
                return [orbits[i], C]
    return None


def witness(g, f, method):
    if f is None:
        return None
    F = factors(g, f)
    if F is None:
        return None
    c = dict(stage=2, host_order=len(g), map=[f[v] for v in range(len(g))],
             factors=[sorted(x) for x in F], method=method)
    verify_expansive(g, c)
    return c


def run_one(g, seconds, states, order, method, backend):
    if method == "coaffine":
        result = pc.find_coaffine_endomorphism(
            g, backend=backend, value_order=order,
            max_states=states, time_limit=seconds)
    else:
        options = dict(backend=backend, max_states=states, time_limit=seconds,
                       sat_max_clauses=6_000_000, sat_memory_limit_mb=3072)
        # The SAT backend rejects value_order='fiber'; it owns its branching policy.
        if backend != "sat":
            options["value_order"] = order
        result = pc.find_octahedral_retraction_coaffine(g, **options)
    return witness(g, result.endomorphism,
                   f"{backend}_{method}_{order}")


def search(g, seconds, sat_seconds, states, seeds):
    attempts = []
    for method in ("coaffine", "octahedral"):
        for order in ("fiber", "label"):
            name = f"native_{method}_{order}"
            try:
                c = run_one(g, seconds, states, order, method, "native")
            except (RuntimeError, ValueError) as exc:
                attempts.append({"method": name, "error": repr(exc)})
                continue
            attempts.append({"method": name, "found": c is not None})
            if c:
                return c, attempts
    for seed in seeds:
        labels = list(g)
        random.Random(seed).shuffle(labels)
        J = nx.Graph()
        J.add_nodes_from(labels)
        J.add_edges_from(g.edges)
        for order in ("fiber", "label"):
            name = f"native_reordered_{seed}_{order}"
            try:
                c = run_one(J, seconds, states, order, "octahedral", "native")
            except (RuntimeError, ValueError) as exc:
                attempts.append({"method": name, "error": repr(exc)})
                continue
            # The graph's labels are still the original 0,...,n-1 labels.
            attempts.append({"method": name, "found": c is not None})
            if c:
                return c, attempts
    try:
        c = run_one(g, sat_seconds, states * 2, "label", "octahedral", "sat")
    except (RuntimeError, ValueError) as exc:
        attempts.append({"method": "sat_octahedral", "error": repr(exc)})
        c = None
    attempts.append({"method": "sat_octahedral", "found": c is not None})
    return c, attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ids", nargs="+")
    parser.add_argument("--seconds", type=float, default=15)
    parser.add_argument("--sat-seconds", type=float, default=30)
    parser.add_argument("--states", type=int, default=1_000_000)
    parser.add_argument("--seeds", type=int, nargs="*", default=[0, 137, 271, 139])
    args = parser.parse_args()
    graphs, ledger = catalog(), read_csv("behavior_order14.csv")
    for i in args.ids:
        assert ledger[i]["status"] == "DIVERGENT" and ledger[i]["stage"] == "2.0"
        g = relabel(fold(nx.complement(graphs[i])))
        g = pared_step(pared_step(g))
        print("HOST", i, len(g), flush=True)
        try:
            c, attempts = search(g, args.seconds, args.sat_seconds,
                                 args.states, args.seeds)
        except (RuntimeError, ValueError) as exc:
            print("LIMIT", i, repr(exc), flush=True)
            continue
        if c:
            print("CERTIFICATE " + json.dumps(
                {"id": i, **c, "attempts": attempts}, separators=(",", ":")),
                flush=True)
        else:
            print("NOT_REDISCOVERED " + json.dumps(
                {"id": i, "attempts": attempts}, separators=(",", ":")),
                flush=True)


if __name__ == "__main__":
    main()
