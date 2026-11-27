"""Check imported ledgers and replay the available inexpensive witnesses.

This independent verifier uses NetworkX, not search engines or SAT.
It does not certify the unbundled divergence or homotopy witnesses.
"""
import csv
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]


def read_csv(name):
    with (ROOT / "census" / name).open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    result = {row["id"]: row for row in rows}
    assert len(rows) == len(result) == 540
    return result


def catalog():
    manifest = json.loads((ROOT / "data/manifest.json").read_text())
    data = {
        4: [nx.complete_graph(4)],
        6: [nx.complete_bipartite_graph(3, 3), nx.circular_ladder_graph(3)],
    }
    for entry in manifest["datasets"]:
        n = entry["order"]
        raw = (ROOT / f"data/cubic{n:02}.g6").read_bytes()
        digest = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
        assert digest == entry["sha"], n
        graphs = [nx.from_graph6_bytes(line) for line in raw.splitlines()]
        assert len(graphs) == entry["count"]
        assert all(
            len(g) == n and nx.is_connected(g)
            and set(dict(g.degree()).values()) == {3}
            for g in graphs
        )
        data[n] = graphs
    result = {str(i): g for i, g in enumerate(data[14], 1)}
    result.update({
        f"4+10:{i}": nx.disjoint_union(data[4][0], g)
        for i, g in enumerate(data[10], 1)
    })
    result.update({
        f"6+8:{i}:{j}": nx.disjoint_union(g, h)
        for i, g in enumerate(data[6], 1)
        for j, h in enumerate(data[8], 1)
    })
    result.update({
        f"4+4+6:{i}": nx.disjoint_union_all(data[4] * 2 + [g])
        for i, g in enumerate(data[6], 1)
    })
    assert len(result) == 540
    return result


def fold(graph, allowed=None):
    """Perform closed-neighborhood folds; optionally remove specified IDs."""
    g = graph.copy()
    closed = {v: set(g[v]) | {v} for v in g}
    remaining = set(g) if allowed is None else set(allowed)
    while remaining:
        move = next(
            ((v, w) for v in sorted(remaining) if v in g
             for w in sorted(g[v]) if closed[v] <= closed[w]),
            None,
        )
        if move is None:
            assert allowed is None, "specified fold sequence blocked"
            break
        v, _ = move
        for u in g[v]:
            closed[u].remove(v)
        del closed[v]
        g.remove_node(v)
        remaining.remove(v)
    return g


def relabel(g):
    return nx.relabel_nodes(g, {v: i for i, v in enumerate(sorted(g))})


def intersection_graph(supports):
    g = nx.Graph()
    g.add_nodes_from(range(len(supports)))
    g.add_edges_from(
        (i, j) for i, j in itertools.combinations(g, 2)
        if supports[i] & supports[j]
    )
    return g


def clique_masks(g, numeric=False):
    family = sorted(tuple(sorted(q)) for q in nx.find_cliques(g))
    masks = [sum(1 << v for v in q) for q in family]
    return sorted(masks) if numeric else masks


def pared_step(g):
    return relabel(fold(intersection_graph(clique_masks(g))))


def verify_308(cubic):
    c = json.loads((ROOT / "certificates/308.json").read_text())
    g = relabel(fold(nx.complement(cubic)))
    base = pared_step(g)
    assert len(base) == c["ground"] == 35
    supports = clique_masks(base, numeric=True)
    raw = intersection_graph(supports)
    assert len(raw) == 577
    active = set(c["active_ids"])
    host = fold(raw, set(raw) - active)
    assert set(host) == active
    assert len(host) == c["order"] == 571
    assert host.number_of_edges() == c["edges"] == 146127
    assert [supports[v] for v in c["active_ids"]] == c["supports"]
    cuts = c["cuts"]
    target = {v for a, b, _ in cuts for v in (a, b)}
    assert len(cuts) == c["dimension"] == 22 and len(target) == 44
    opposite = {u: v for a, b, _ in cuts for u, v in ((a, b), (b, a))}
    assert all(
        host.has_edge(a, b) == (opposite[a] != b)
        for a, b in itertools.combinations(target, 2)
    )
    for a, b, cut in cuts:
        assert supports[a] & ~cut == 0 and supports[b] & cut == 0
    assert all(
        any(supports[v] & ~cut == 0 or supports[v] & cut == 0
            for _, _, cut in cuts)
        for v in host
    )
    projection = {int(v): w for v, w in c["projection"].items()}
    assert set(projection) == active and set(projection.values()) == target
    assert all(projection[v] == v for v in target)
    assert all(
        projection[u] == projection[v]
        or host.has_edge(projection[u], projection[v])
        for u, v in host.edges
    )
    print("308: verified O_22 retraction")


def main():
    graphs = catalog()
    b = read_csv("behavior_order14.csv")
    h = read_csv("homotopy_order14.csv")
    assert graphs.keys() == b.keys() == h.keys()
    assert Counter(row["status"] for row in b.values()) == {
        "DIVERGENT": 502, "CONVERGENT": 18, "UNRESOLVED_IN_LEDGER": 20,
    }
    for i, row in b.items():
        for field in ("G", "KG"):
            profile = json.loads(row[field])
            assert isinstance(profile, dict)
            assert all(
                int(d) >= 0 and type(n) is int and n > 0
                for d, n in profile.items()
            )
            assert profile == json.loads(h[i][field])
        assert h[i]["same"] == str(json.loads(row["G"]) == json.loads(row["KG"]))
    assert sum(row["same"] == "True" for row in h.values()) == 19
    unresolved = {i for i, r in b.items() if r["status"] == "UNRESOLVED_IN_LEDGER"}
    assert unresolved == set((ROOT / "census/unresolved_order14.txt").read_text().split())
    assert {
        i for i, r in b.items()
        if r["source"] == "earlier_archived_certificate_not_replayed"
    } == {"150", "252", "270", "472"}
    print("540 data and ledger records checked")
    for i, row in b.items():
        if row["status"] != "CONVERGENT":
            continue
        g = relabel(fold(nx.complement(graphs[i])))
        for _ in range(int(float(row["stage"]))):
            g = pared_step(g)
        assert len(g) == 1, i
    print("18 point reductions independently replayed")
    verify_308(graphs["308"])


if __name__ == "__main__":
    main()
