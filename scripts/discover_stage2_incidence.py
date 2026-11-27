"""Incidence-cut rediscovery on implicit pared second clique iterates."""
import argparse
import json

import networkx as nx
import pycliques as pc

from validate_census import catalog, fold, pared_step, relabel, read_csv
from verify_coaffine import verify_expansive


def certificate(host, ans):
    if ans.status.name != "FOUND":
        return None
    graph = host.to_networkx(max_edges=8_000_000)
    if graph is None:
        return None
    assert pc.verify_induced_octahedron(graph, ans.antipodal_pairs)
    target_vertices = {v for pair in ans.antipodal_pairs for v in pair}
    target = graph.subgraph(target_vertices).copy()
    inclusion = {v: v for v in target}
    assert pc.verify_retraction_witness(
        graph, target, (ans.retraction, inclusion))
    assert pc.verify_coaffine_endomorphism(graph, ans.endomorphism)
    nodes = sorted(graph)
    relabel_map = {v: i for i, v in enumerate(nodes)}
    G = nx.relabel_nodes(graph, relabel_map, copy=True)
    f = {relabel_map[v]: relabel_map[w]
         for v, w in ans.endomorphism.items()}
    assert len(ans.antipodal_pairs) >= 3
    factors = [[relabel_map[v] for v in pair] for pair in ans.antipodal_pairs[:3]]
    proof = {"stage": 2, "host_order": len(G),
             "map": [f[v] for v in range(len(G))],
             "factors": factors, "method": "implicit_incidence_hybrid"}
    verify_expansive(G, proof)
    return proof


def build_host(cubic):
    p0 = relabel(fold(nx.complement(cubic)))
    p1 = pared_step(p0)
    raw = pc.ImplicitCliqueGraph.from_graph(p1, max_cliques=3000)
    if raw is None:
        return None, {"reason": "P2 clique count exceeded 3000"}
    pared = raw.pare()
    if not pc.verify_implicit_paring(raw, pared):
        raise AssertionError("implicit paring certificate failed")
    return pared.core, {"p0": len(p0), "p1": len(p1),
                        "raw_p2": len(raw), "pared_p2": len(pared.core),
                        "ground": len(pared.core.base_vertices)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ids", nargs="+")
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--steps", type=int, default=250000)
    parser.add_argument("--states", type=int, default=1000000)
    parser.add_argument("--seeds", type=int, nargs="*", default=[139, 137, 271])
    args = parser.parse_args()
    graphs, ledger = catalog(), read_csv("behavior_order14.csv")
    executable = pc.build_incidence_native()
    for identifier in args.ids:
        assert ledger[identifier]["status"] == "DIVERGENT"
        assert ledger[identifier]["stage"] == "2.0"
        host, info = build_host(graphs[identifier])
        print("HOST", identifier, json.dumps(info), flush=True)
        if host is None:
            print("LIMIT", identifier, json.dumps(info), flush=True)
            continue
        for seed in args.seeds:
            try:
                ans = pc.find_octahedral_retraction_incidence(
                    host, backend="hybrid", seed=seed,
                    max_steps=args.steps, max_states=args.states,
                    time_limit=args.seconds, executable=executable,
                    sat_max_clauses=6000000, sat_memory_limit_mb=3072)
                proof = certificate(host, ans)
            except (RuntimeError, ValueError) as exc:
                print("ERROR", identifier, seed, repr(exc), flush=True)
                continue
            if proof:
                print("CERTIFICATE " + json.dumps(
                    {"id": identifier, "seed": seed, **proof,
                     "dimension": len(ans.antipodal_pairs),
                     "states": ans.states_explored},
                    separators=(",", ":")), flush=True)
                break
            print("ATTEMPT", json.dumps(
                {"id": identifier, "seed": seed,
                 "status": ans.status.name,
                 "states": ans.states_explored,
                 "uncovered": ans.uncovered,
                 "partial_cuts": len(ans.partial_cuts)},
                separators=(",", ":")), flush=True)
        else:
            print("NOT_REDISCOVERED", identifier, flush=True)


if __name__ == "__main__":
    main()
