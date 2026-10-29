"""How often does the (library-blind) engine find the textbook disconnection?  python tests/coverage.py"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from retroexp import engine, library as L
from retroexp.chem import canon
from retroexp.buildingblocks import is_available

rows = []
t0 = time.time()
for e in L.LIBRARY.values():
    for ri, r in enumerate(e["routes"]):
        tree = L.route_tree(e, ri)
        top_step = tree["step"]
        text_first = engine.norm_set(top_step.reactants.split("."))
        text_leaves = engine.norm_set(x for x in _leaves(tree)) if False else None
        def leaves(n, root=True):
            # textbook route cut where the SAME availability rule as the engine says "purchasable"
            return [n["smiles"]] if (not n["step"] or (not root and is_available(n["smiles"]))) else [x for c in n["children"] for x in leaves(c, False)]
        text_leaves = engine.norm_set(leaves(tree))
        props = engine.propose(e["smiles"])
        rank = next((i + 1 for i, p in enumerate(props) if engine.norm_set(p.precursors) == text_first), None)
        routes = engine.build_routes(e["smiles"], n_routes=5)
        rr = next((i + 1 for i, n in enumerate(routes) if engine.norm_set(engine.starting_materials(n)) == text_leaves), None)
        rows.append((e["name"], r["name"], rank, rr))
n = len(rows)
print(f"{'molecule':44s} {'route':26s} first-step rank   route rank(top5)")
for name, rn, rank, rr in rows:
    print(f"{name[:43]:44s} {rn[:25]:26s} {str(rank or '-'):>10s} {str(rr or '-'):>16s}")
f1 = sum(1 for *_, rk, _ in rows if rk == 1); f3 = sum(1 for *_, rk, _ in rows if rk and rk <= 3); f10 = sum(1 for *_, rk, _ in rows if rk)
r1 = sum(1 for *_, rr in rows if rr == 1); r5 = sum(1 for *_, rr in rows if rr)
print(f"\nroutes tested: {n}\nfirst disconnection = textbook:  top-1 {f1}/{n}   top-3 {f3}/{n}   anywhere in top 15 {f10}/{n}")
print(f"same starting materials as textbook route:  top-1 {r1}/{n}   top-5 {r5}/{n}      ({time.time()-t0:.0f}s)")
