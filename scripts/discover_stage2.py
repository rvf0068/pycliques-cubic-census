"""Optional bounded rediscovery with pinned pycliques; not part of fast CI."""
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
    joined.add_edges_from((i, j) for i, j in itertools.combinations(joined, 2)
                          if all(g.has_edge(v,w) for v in orbits[i] for w in orbits[j]))
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


def search(g, seconds):
    for order in ("fiber", "label"):
        r = pc.find_coaffine_endomorphism(g, backend="native", value_order=order,
                 max_states=200000, time_limit=seconds)
        c = witness(g, r.endomorphism, f"native_coaffine_{order}")
        if c:
            return c
    for order in ("fiber", "label"):
        r = pc.find_octahedral_retraction_coaffine(g, backend="native", value_order=order,
                 max_states=200000, time_limit=seconds)
        c = witness(g, r.endomorphism, f"native_octahedral_{order}")
        if c:
            return c
    for seed in (137,139):
        labels = list(g)
        random.Random(seed).shuffle(labels)
        J = nx.Graph()
        J.add_nodes_from(labels)
        J.add_edges_from(g.edges)
        r = pc.find_octahedral_retraction_coaffine(J, backend="native", value_order="fiber",
                 max_states=200000, time_limit=seconds)
        c = witness(g, r.endomorphism, f"native_reordered_{seed}")
        if c:
            return c
    r = pc.find_octahedral_retraction_coaffine(g, backend="sat",
             max_states=1000000, time_limit=2*seconds, sat_max_clauses=6000000,
             sat_memory_limit_mb=3072)
    return witness(g, r.endomorphism, "sat_octahedral")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ids", nargs="+")
    parser.add_argument("--seconds", type=float, default=5)
    args = parser.parse_args()
    graphs, ledger = catalog(), read_csv("behavior_order14.csv")
    for i in args.ids:
        assert ledger[i]["status"] == "DIVERGENT" and ledger[i]["stage"] == "2.0"
        g = relabel(fold(nx.complement(graphs[i])))
        g = pared_step(pared_step(g))
        print("HOST", i, len(g), flush=True)
        try:
            c = search(g, args.seconds)
        except (RuntimeError, ValueError) as exc:
            print("LIMIT", i, repr(exc), flush=True)
            continue
        if c:
            print("CERTIFICATE " + json.dumps({"id": i, **c}, separators=(",",":")), flush=True)
        else:
            print("NOT_REDISCOVERED", i, flush=True)


if __name__ == "__main__":
    main()
