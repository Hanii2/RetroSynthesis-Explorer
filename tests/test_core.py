"""Self-test:  python tests/test_core.py   (no internet needed)"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rdkit import Chem
from retroexp import compare, engine, library as L, quiz, safety, templates
from retroexp.buildingblocks import is_available
from retroexp.chem import canon

# 1) library consistency: parseable, every step links up, carbon never created from nothing
n_steps = 0
for e in L.LIBRARY.values():
    assert Chem.MolFromSmiles(e["smiles"]), e["name"]
    for r in e["routes"]:
        made = {canon(s.product) for s in r["steps"]}
        assert canon(e["smiles"]) in made, (e["name"], "target not produced")
        for s in r["steps"]:
            n_steps += 1
            reac = [Chem.MolFromSmiles(x) for x in s.reactants.split(".")]
            assert all(reac) and Chem.MolFromSmiles(s.product), (e["name"], s.reaction)
            c_in = sum(a.GetSymbol() == "C" for m in reac for a in m.GetAtoms())
            c_out = sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(s.product).GetAtoms())
            assert c_out <= c_in, (e["name"], s.reaction, "carbon balance")
            for x in s.reactants.split("."):
                assert is_available(x) or canon(x) in made, (e["name"], s.reaction, x, "reactant neither available nor made")
print(f"OK  library: {len(L.LIBRARY)} molecules, {n_steps} textbook steps, all consistent")

# 2) every template fires on at least one textbook target
fired = set()
for e in L.LIBRARY.values():
    for p in engine.propose(e["smiles"], 40):
        fired.add(p.template.id)
silent = [t.id for t in templates.TEMPLATES if t.id not in fired]
assert len(silent) <= 8, silent
print(f"OK  templates: {len(templates.TEMPLATES)} compile; {len(templates.TEMPLATES) - len(silent)} fire on library targets (silent: {silent})")

# 3) textbook examples that must come out right
def top(smi, k=1):
    return [(p.template.id, sorted(p.precursors)) for p in engine.propose(smi)[:k]]
assert top("CC(=O)Oc1ccccc1C(=O)O")[0][0] == "acetylate_aryl"
assert engine.norm_set(engine.propose("CC(C)(C)OC")[0].precursors) == engine.norm_set(["CC(C)(C)O", "CBr"])         # tertiary ALKOXIDE + methyl halide
assert not any("CC(C)(C)Br" in p.precursors for p in engine.propose("CC(C)(C)OC"))                                 # never t-BuBr + MeO-
assert engine.propose("c1ccc(-c2ccccc2)cc1")[0].template.id == "suzuki"
assert engine.propose("CC(=O)Nc1ccc([N+](=O)[O-])cc1")[0].template.id == "nitration"
print("OK  textbook cases: aspirin, MTBE (Williamson direction), biphenyl, 4-nitroacetanilide")

# 4) chemistry checks flag impossible disconnections
def has(smi, tid, rule, sev="fail"):
    return any(p.template.id == tid and any(f.rule == rule and f.severity == sev for f in p.findings) for p in engine.propose(smi, 40))
assert has("COC(=O)c1ccc([N+](=O)[O-])cc1", "nitration", "eas_regiochemistry")             # para-nitration of methyl benzoate
assert has("CC(C)Cc1ccccc1", "fc_alkylation", "friedel_crafts_alkylation")                    # isobutyl rearrangement
assert has("CC(=O)Oc1ccccc1C(=O)O", "grignard_co2", "organometallic_incompatible_groups")    # ester + Grignard
assert has("CNCCC(O)c1ccccc1", "grignard_carbonyl", "organometallic_acidic_protons")         # NH quenches RMgBr
print("OK  rules: regiochemistry, Friedel-Crafts rearrangement, Grignard incompatibilities")

# 5) route builder: complete, no failing step, scores sorted
r = engine.build_routes("CC(=O)Nc1ccc(O)cc1")
assert r and r[0].complete and not r[0].has_fail
print("OK  routes: paracetamol ->", " + ".join(sorted(engine.starting_materials(r[0]))))

# 6) agreement with the textbook (the engine never sees the library)
res = compare.run_all()
assert res["first_top3"] / res["n"] > 0.80 and res["route_top5"] / res["n"] > 0.70, res
print(f"OK  engine vs textbook: last step #1 in {res['first_top1']}/{res['n']}, top-3 {res['first_top3']}/{res['n']}; "
      f"same starting materials in top-5 routes {res['route_top5']}/{res['n']}")

# 7) safety filter
for smi in ["CNC(C)Cc1ccccc1", "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "CC(C)OP(C)(=O)F",
            "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"]:
    assert safety.check(smi)[0], smi
assert safety.check("x", "methamphetamine")[0]
assert not any(safety.check(e["smiles"], e["name"])[0] for e in L.LIBRARY.values())
print("OK  safety: harmful classes blocked, all 56 library molecules pass")

# 8) quiz
q = quiz.make_question(random.Random(1))
assert len(set(q["options"])) == 4 and q["options"][q["answer"]] in (q["step"].reaction, q["step"].reagents)
print("OK  quiz")
print("ALL TESTS PASSED")
