"""Incidence-cut rediscovery on implicit pared second clique iterates."""
import argparse
import json

import networkx as nx
import pycliques as pc

from validate_census import catalog, fold, pared_step, relabel, read_csv
from verify_coaffine import verify_expansive


def fold_active_ids(graph):
    """Return surviving original IDs under validate_census.fold's order."""
    closed = {v: set(graph[v]) | {v} for v in graph}
    remaining = set(graph)
    while remaining:
        move = next(
            ((v, w) for v in sorted(remaining) for w in sorted(graph[v])
             if closed[v] <= closed[w]),
            None,
        )
        if move is None:
            break
        v, _ = move
        for u in graph[v]:
            closed[u].remove(v)
        remaining.remove(v)
    return sorted(remaining)


def build_host(cubic):
    p0 = relabel(fold(nx.complement(cubic)))
    p1 = pared_step(p0)
    raw = pc.ImplicitCliqueGraph.from_graph(p1, max_cliques=3000)
    if raw is None:
        return None, None, None, {"reason": "P2 clique count exceeded 3000"}
    pared = raw.pare()
    if not pc.verify_implicit_paring(raw, pared):
        raise AssertionError("implicit paring certificate failed")

    # The implicit constructor sorts support masks numerically, whereas
    # validate_census.clique_masks sorts the corresponding vertex tuples.
    # Conjugate the implicit IDs to the canonical explicit P2 IDs by support.
    supports = clique_masks(p1, numeric=False)
    raw_graph = intersection_graph(supports)
    active_raw_ids = fold_active_ids(raw_graph)
    raw_id_by_mask = {mask: i for i, mask in enumerate(supports)}
    canonical_id = {raw_id: i for i, raw_id in enumerate(active_raw_ids)}
    host = pared.core
    to_canonical = {}
    for v in host:
        support_mask = sum(1 << u for u in host.clique(v))
        raw_id = raw_id_by_mask[support_mask]
        to_canonical[v] = canonical_id[raw_id]

    canonical_graph = pared_step(p1)
    if len(host) != len(canonical_graph):
        raise AssertionError("implicit and explicit pared P2 orders differ")
    if set(to_canonical.values()) != set(canonical_graph):
        raise AssertionError("support relabeling is not bijective")
    implicit_graph = host.to_networkx(max_edges=8_000_000)
    if implicit_graph is None:
        raise AssertionError("implicit host exceeds verification edge cap")
    for u in implicit_graph:
        for v in implicit_graph:
            if implicit_graph.has_edge(u, v) != canonical_graph.has_edge(
                to_canonical[u], to_canonical[v]
            ):
                raise AssertionError("support-based P2 relabeling is not an isomorphism")
    info = {"p0": len(p0), "p1": len(p1), "raw_p2": len(raw),
            "pared_p2": len(host), "ground": len(host.base_vertices)}
    return host, canonical_graph, to_canonical, info


def certificate(host, canonical_graph, to_canonical, ans):
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

    # Translate both the endomorphism and its invariant factors to the
    # explicit P2 labels used by the repository's solver-free verifier.
    f = {to_canonical[v]: to_canonical[w]
         for v, w in ans.endomorphism.items()}
    factors = [[to_canonical[v] for v in pair]
               for pair in ans.antipodal_pairs[:3]]
    proof = {"stage": 2, "host_order": len(canonical_graph),
             "map": [f[v] for v in range(len(canonical_graph))],
             "factors": factors, "method": "implicit_incidence_hybrid",
             "target_type": "octahedral",
             "octahedral_dimension": len(ans.antipodal_pairs)}
    verify_expansive(canonical_graph, proof)
    return proof

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
        host, canonical_graph, to_canonical, info = build_host(graphs[identifier])
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
                proof = certificate(host, canonical_graph, to_canonical, ans)
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
