"""
Rule-based retrosynthesis engine.

  propose(target)  -> ranked one-step disconnections (Proposal)
  build_routes(target) -> complete routes down to available starting materials (Node trees)

It is a TEACHING tool built on textbook disconnection rules (see templates.py), not a trained AI model.
"""
from dataclasses import dataclass, field
from functools import lru_cache
from types import SimpleNamespace
from typing import List, Optional

from rdkit import Chem

from .buildingblocks import is_available, name_of
from .chem import canon
from .rules import PENALTY, Finding, check_step
from .templates import TEMPLATES, Template


@dataclass
class Proposal:
    template: Template
    target: str                 # the SMILES string that was used (atom indices of broken/changed refer to it)
    precursors: List[str]
    broken: List[int]
    changed: List[int]
    score: float = 0.0
    parts: dict = field(default_factory=dict)
    findings: list = field(default_factory=list)


@dataclass
class Node:
    smiles: str
    proposal: Optional[Proposal] = None
    children: list = field(default_factory=list)
    available: bool = False
    score: float = 0.0

    @property
    def n_steps(self):
        return (1 + sum(c.n_steps for c in self.children)) if self.proposal else 0

    def leaves(self):
        if not self.proposal:
            return [self]
        return [x for c in self.children for x in c.leaves()]

    @property
    def has_fail(self):
        own = bool(self.proposal and any(f.severity == "fail" for f in self.proposal.findings))
        return own or any(c.has_fail for c in self.children)

    @property
    def complete(self):
        return all(l.available for l in self.leaves())


def _changes(target_mol, prods):
    where = {}
    for pi, p in enumerate(prods):
        for a in p.GetAtoms():
            if a.HasProp("react_atom_idx"):
                where[int(a.GetProp("react_atom_idx"))] = (pi, a.GetIdx())
    broken, changed = [], []
    for b in target_mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in where and j in where:
            (pi, ai), (pj, aj) = where[i], where[j]
            if pi != pj:
                broken.append(b.GetIdx())
            else:
                pb = prods[pi].GetBondBetweenAtoms(ai, aj)
                if pb is None:
                    broken.append(b.GetIdx())
                elif pb.GetBondType() != b.GetBondType():
                    changed.append(b.GetIdx())
        elif (i in where) != (j in where):          # one end was removed / replaced
            broken.append(b.GetIdx())
    return broken, changed


EAS_TEMPLATES = {"nitration", "halogenation", "chlorosulfonation", "fc_acylation", "fc_acetylation_anhydride", "fc_alkylation"}
_META = [Chem.MolFromSmarts(x) for x in ("[N+](=O)[O-]", "[#6]=O", "C#N", "S(=O)=O", "C(F)(F)F", "[N+;!$([N+](=O)[O-])]")]


def _director(mol, ring_atom, sub_atom):
    """'op' (ortho/para director), 'meta' or None for the substituent at sub_atom."""
    sub = mol.GetAtomWithIdx(sub_atom)
    for q in _META:
        for match in mol.GetSubstructMatches(q):
            if match[0] == sub_atom:
                return "meta"
    if sub.GetSymbol() in ("O", "N", "F", "Cl", "Br", "I", "S"):
        return "op"
    if sub.GetSymbol() == "C" and (sub.GetIsAromatic() or sub.GetHybridization() == Chem.HybridizationType.SP3):
        return "op"
    return None


def eas_findings(mol, atom_idx):
    """Directing-effect check for an electrophilic aromatic substitution at ring atom `atom_idx` of the target."""
    ring = next((r for r in mol.GetRingInfo().AtomRings() if atom_idx in r and len(r) == 6
                 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)), None)
    if ring is None:
        return []
    pos = ring.index(atom_idx)
    out, supported = [], False
    names = {1: "ortho", 2: "meta", 3: "para"}
    for k, ra in enumerate(ring):
        if ra == atom_idx:
            continue
        for nb in mol.GetAtomWithIdx(ra).GetNeighbors():
            if nb.GetIdx() in ring or nb.GetAtomicNum() == 1:
                continue
            d = (k - pos) % 6
            rel = names[min(d, 6 - d)]
            kind = _director(mol, ra, nb.GetIdx())
            if kind == "meta" and rel != "meta":
                out.append(Finding(1, "fail", "eas_regiochemistry",
                                   f"the new group sits {rel} to a meta-directing group ({nb.GetSymbol()}-substituent): "
                                   f"electrophilic substitution would put it meta, so this disconnection does not work."))
            elif kind == "op" and rel == "meta":
                out.append(Finding(1, "warn", "eas_regiochemistry",
                                   f"the new group sits meta to an ortho/para-directing group ({nb.GetSymbol()}-substituent): disfavoured."))
            elif kind is not None:
                supported = True
    if supported and not out:
        out.append(Finding(1, "info", "eas_regiochemistry", "directing effects of the groups already on the ring agree with this position."))
    return out


@lru_cache(maxsize=4096)
def propose(target_smiles: str, limit: int = 15):
    mol = Chem.MolFromSmiles(target_smiles)
    if mol is None:
        return []
    mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol)
    t_heavy = mol.GetNumHeavyAtoms()
    t_can = canon(target_smiles)
    seen, out = set(), []
    for t in TEMPLATES:
        try:
            results = t.rxn.RunReactants((mol,))
        except Exception:
            continue
        for prods in results:
            try:
                for p in prods:
                    Chem.SanitizeMol(p)
                if t.check and not t.check(prods):
                    continue
                pre = [canon(Chem.MolToSmiles(p)) for p in prods]
            except Exception:
                continue
            if None in pre or t_can in pre or any(Chem.MolFromSmiles(x).GetNumHeavyAtoms() > t_heavy + 8 for x in pre):
                continue
            key = (t.id, frozenset(pre))
            if key in seen:
                continue
            seen.add(key)
            broken, changed = _changes(mol, prods)
            out.append(Proposal(t, target_smiles, sorted(set(pre)), broken, changed))
    for pr in out:
        avail = sum(is_available(x) for x in pr.precursors) / len(pr.precursors)
        biggest = max(Chem.MolFromSmiles(x).GetNumHeavyAtoms() for x in pr.precursors)
        simp = max(0.0, min(1.0, 1 - biggest / t_heavy))
        step = SimpleNamespace(product=target_smiles, reactants=pr.precursors)
        findings = check_step(1, step, False)
        if pr.template.id in EAS_TEMPLATES:
            for bi in pr.broken:
                b = mol.GetBondWithIdx(bi)
                for a in (b.GetBeginAtom(), b.GetEndAtom()):
                    if a.GetIsAromatic():
                        findings += eas_findings(mol, a.GetIdx())
        rules = max(0.0, 1.0 - sum(PENALTY[f.severity] for f in findings))
        findings = [f for f in findings if not (f.rule == "atom_bookkeeping" and
                                                (pr.template.kind in ("FGI", "Protecting") or pr.template.id in EAS_TEMPLATES))]
        pr.findings = findings
        pr.parts = {"reaction reliability": pr.template.weight, "precursors available": avail,
                    "simplification": simp, "chemistry check": rules}
        pr.score = 0.40 * pr.template.weight + 0.25 * avail + 0.15 * simp + 0.20 * rules
        if any(f.severity == "fail" for f in findings):
            pr.score *= 0.6                      # a chemically impossible disconnection must not rank high
    out.sort(key=lambda p: -p.score)
    return out[:limit]


def _route_score(node):
    if not node.proposal:
        return 1.0 if node.available else 0.0
    kids = [_route_score(k) for k in node.children]
    return 0.5 * node.proposal.score + 0.5 * sum(kids) / len(kids)


def _expand(smiles, depth, path, beam):
    c = canon(smiles)
    if is_available(c):
        return Node(c, available=True, score=1.0)
    if depth == 0 or c in path:
        return Node(c, available=False, score=0.0)
    best = None
    for pr in propose(c)[:beam]:
        kids = [_expand(x, depth - 1, path | {c}, beam) for x in pr.precursors]
        node = Node(c, pr, kids)
        node.score = _route_score(node) * (0.97 ** max(0, node.n_steps - 1))
        if best is None or (node.complete, not node.has_fail, node.score) > (best.complete, not best.has_fail, best.score):
            best = node
    return best or Node(c, available=False)


def build_routes(target_smiles, depth=5, beam=3, n_routes=5):
    """One complete route per top-level disconnection, best first."""
    c = canon(target_smiles)
    routes = []
    for pr in propose(c)[:max(n_routes * 2, 8)]:
        kids = [_expand(x, depth - 1, frozenset([c]), beam) for x in pr.precursors]
        node = Node(c, pr, kids)
        node.score = _route_score(node) * (0.97 ** max(0, node.n_steps - 1))
        routes.append(node)
    clean = [r for r in routes if not r.has_fail]
    routes = clean or routes                       # routes with a chemically impossible step are dropped (unless nothing else exists)
    routes.sort(key=lambda n: (-int(n.complete), -n.score))
    seen, uniq = set(), []
    for r in routes:
        key = frozenset(l.smiles for l in r.leaves())
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)
    return uniq[:n_routes]


def starting_materials(node):
    return {l.smiles for l in node.leaves()}


def normalize_for_compare(smiles):
    """Canonical SMILES in which every single aromatic/alkyl halogen (F, Cl, I) is written as Br, so that
    'Cl-CH2-oxirane' and 'Br-CH2-oxirane' count as the same type of building block when comparing routes."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return smiles
    m = Chem.RWMol(m)
    for a in m.GetAtoms():
        if a.GetSymbol() in ("F", "Cl", "I") and a.GetDegree() == 1:
            nb = a.GetNeighbors()[0]
            hal = [x for x in nb.GetNeighbors() if x.GetSymbol() in ("F", "Cl", "Br", "I")]
            if nb.GetSymbol() == "C" and len(hal) == 1 and not (nb.GetSymbol() == "C" and nb.GetIsAromatic() is False and False):
                a.SetAtomicNum(35)
    out = m.GetMol()
    Chem.SanitizeMol(out)
    return canon(Chem.MolToSmiles(out))


def norm_set(smiles_iter):
    return frozenset(normalize_for_compare(x) for x in smiles_iter)
