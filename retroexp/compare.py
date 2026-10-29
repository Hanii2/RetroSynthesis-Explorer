"""Compare the library-blind engine with a textbook route."""
from . import engine
from .buildingblocks import is_available
from . import library as L


def textbook_leaves(tree, root=True):
    """Starting materials of a textbook route, cut where the SAME availability rule says 'purchasable'."""
    if not tree["step"] or (not root and is_available(tree["smiles"])):
        return [tree["smiles"]]
    return [x for c in tree["children"] for x in textbook_leaves(c, False)]


def compare(entry, route_index=0, depth=5):
    tree = L.route_tree(entry, route_index)
    first = engine.norm_set(tree["step"].reactants.split("."))
    leaves = engine.norm_set(textbook_leaves(tree))
    props = engine.propose(entry["smiles"])
    rank = next((i + 1 for i, p in enumerate(props) if engine.norm_set(p.precursors) == first), None)
    routes = engine.build_routes(entry["smiles"], depth=depth, n_routes=5)
    rr = next((i + 1 for i, n in enumerate(routes) if engine.norm_set(engine.starting_materials(n)) == leaves), None)
    return {"first_step_rank": rank, "route_rank": rr}


def run_all():
    rows = []
    for e in L.LIBRARY.values():
        for i, r in enumerate(e["routes"]):
            rows.append((e["name"], r["name"], compare(e, i)))
    n = len(rows)
    ok = lambda f: sum(1 for *_, c in rows if f(c))
    return {"n": n,
            "first_top1": ok(lambda c: c["first_step_rank"] == 1), "first_top3": ok(lambda c: c["first_step_rank"] and c["first_step_rank"] <= 3),
            "route_top1": ok(lambda c: c["route_rank"] == 1), "route_top5": ok(lambda c: c["route_rank"]), "rows": rows}
