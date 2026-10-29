"""Colour schemes shared by the 2D drawing and the 3D viewer."""
import colorsys

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

SCHEMES = ["Vivid colours", "Functional groups", "Rainbow", "Classic (CPK)"]

VIVID = {"H": "#dfe6e9", "C": "#27ae60", "N": "#3498db", "O": "#e74c3c", "S": "#f1c40f", "P": "#e67e22",
         "F": "#1abc9c", "Cl": "#8e44ad", "Br": "#a0522d", "I": "#6c3483", "B": "#fd79a8", "Si": "#b2bec3"}
CLASSIC = {"H": "#ffffff", "C": "#909090", "N": "#3050f8", "O": "#ff0d0d", "S": "#ffff30", "P": "#ff8000",
           "F": "#90e050", "Cl": "#1ff01f", "Br": "#a62929", "I": "#940094"}
OTHER = "#b2bec3"
NEUTRAL, H_COLOR = "#cfd8dc", "#ecf0f1"

# (name, SMARTS, colour). Later entries override earlier ones, so specific groups win over general ones.
FUNCTIONAL_GROUPS = [
    ("Aromatic ring", "[a]", "#a29bfe"),
    ("Alkene / alkyne", "[CX3]=[CX3]", "#fdcb6e"),
    ("Alkyne", "[CX2]#[CX2]", "#fdcb6e"),
    ("Ether", "[OD2;$(O([#6])[#6]);!$(O-C=O)]", "#55efc4"),
    ("Alcohol / phenol (OH)", "[OX2H]", "#ff7675"),
    ("Amine", "[NX3;!$(NC=[O,S]);!$(NS(=O)=O);!$([N+])]", "#74b9ff"),
    ("Thiol / sulfide", "[SX2]", "#f1c40f"),
    ("Halogen", "[F,Cl,Br,I]", "#00b894"),
    ("Nitrile", "[CX2]#[NX1]", "#6c5ce7"),
    ("Nitro", "[N+](=O)[O-]", "#6c5ce7"),
    ("Carbonyl (C=O)", "[CX3]=[OX1]", "#e17055"),
    ("Ester", "[CX3](=[OX1])[OX2]", "#e84393"),
    ("Carboxylic acid", "[CX3](=[OX1])[OX2H1]", "#d63031"),
    ("Amide", "[CX3](=[OX1])[NX3]", "#0984e3"),
]
_PATTERNS = [(n, Chem.MolFromSmarts(s), c) for n, s, c in FUNCTIONAL_GROUPS]
assert all(p is not None for _, p, _ in _PATTERNS)


def _rgb(hex_):
    h = hex_.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _heavy_view(mol):
    """Heavy-atom copy (no H) and the map heavy_index -> index in the original molecule."""
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return Chem.RemoveHs(mol), idx


def group_assignment(mol):
    """Returns (atom_group {heavy idx -> (name, hex)}, bond_group {bond idx -> hex}, found [(name, hex, n)])."""
    heavy, _ = _heavy_view(mol) if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) else (mol, None)
    atom_group, bond_group, counts = {}, {}, {}
    for name, patt, hexc in _PATTERNS:
        for match in heavy.GetSubstructMatches(patt):
            for i in match:
                atom_group[i] = (name, hexc)
            for b in heavy.GetBonds():
                if b.GetBeginAtomIdx() in match and b.GetEndAtomIdx() in match:
                    bond_group[b.GetIdx()] = hexc
    for _, (name, hexc) in atom_group.items():
        counts[(name, hexc)] = counts.get((name, hexc), 0) + 1
    return atom_group, bond_group, [(n, c, k) for (n, c), k in counts.items()]


def atom_colors(mol, scheme):
    """Colour (hex) for EVERY atom of `mol` (hydrogens included) for the 3D viewer."""
    cols = {}
    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if scheme == "Functional groups":
        ag, _, _ = group_assignment(mol)
        for k, i in enumerate(heavy_idx):
            cols[i] = ag[k][1] if k in ag else NEUTRAL
    elif scheme == "Rainbow":
        n = max(1, len(heavy_idx) - 1)
        for k, i in enumerate(heavy_idx):
            r, g, b = colorsys.hsv_to_rgb(0.80 * k / n, 0.75, 0.95)
            cols[i] = "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
    else:
        table = CLASSIC if scheme == "Classic (CPK)" else VIVID
        for a in mol.GetAtoms():
            cols[a.GetIdx()] = table.get(a.GetSymbol(), OTHER)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            cols[a.GetIdx()] = H_COLOR if scheme in ("Functional groups", "Rainbow") else cols[a.GetIdx()]
    return cols


def legend(mol, scheme):
    """[(label, hex)] shown under the 2D picture."""
    if scheme == "Functional groups":
        _, _, found = group_assignment(mol)
        return [(f"{n} x{k}", c) for n, c, k in found] or [("no functional groups: only C-H skeleton", NEUTRAL)]
    if scheme == "Rainbow":
        return [("atom 1 (red) -> last atom (violet)", "#e74c3c")]
    table = CLASSIC if scheme == "Classic (CPK)" else VIVID
    syms = sorted({a.GetSymbol() for a in mol.GetAtoms()} | {"C", "H"}, key=lambda s: "CHNOSPFClBrI".find(s))
    return [(s, table.get(s, OTHER)) for s in syms if s in table]


def draw2d(smiles, scheme="Vivid colours", size=(460, 360), indices=False, stereo=True, legend="", hl_bonds=None, hl_color="#e63946", hl2_bonds=None):
    """PNG bytes of the 2D structure in the chosen colour scheme."""
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    o = d.drawOptions()
    o.addStereoAnnotation, o.addAtomIndices, o.bondLineWidth = stereo, indices, 2
    o.legendFontSize = 14
    hl_atoms, hl_bonds, ac, bc = [], [], {}, {}
    if scheme == "Classic (CPK)":
        pass                                             # RDKit default palette
    elif scheme == "Vivid colours":
        pal = {Chem.GetPeriodicTable().GetAtomicNumber(s): _rgb(c) for s, c in VIVID.items() if s != "H"}
        pal[6] = _rgb("#1e8449")                         # darker green: readable skeleton
        o.setAtomPalette(pal)
    elif scheme == "Functional groups":
        ag, bg, _ = group_assignment(mol)
        hl_atoms = list(ag)
        ac = {i: _rgb(c) for i, (_, c) in ag.items()}
        hl_bonds = list(bg)
        bc = {i: _rgb(c) for i, c in bg.items()}
        o.setAtomPalette({6: (0.1, 0.1, 0.1), 7: (0, 0.2, 0.8), 8: (0.8, 0, 0), 16: (0.6, 0.5, 0)})
    elif scheme == "Rainbow":
        n = max(1, mol.GetNumAtoms() - 1)
        for i in range(mol.GetNumAtoms()):
            r, g, b = colorsys.hsv_to_rgb(0.80 * i / n, 0.75, 0.95)
            ac[i] = (r, g, b)
        hl_atoms = list(ac)
        o.setAtomPalette({6: (0.1, 0.1, 0.1)})
    if hl_bonds or hl2_bonds:                              # disconnection pictures: broken bonds red, changed orange
        hl_atoms, ac, hl_bonds_all, bc = [], {}, [], {}
        for idxs, col in ((hl2_bonds or [], "#ff9f1c"), (hl_bonds or [], hl_color)):
            for bi in idxs:
                b = mol.GetBondWithIdx(bi)
                hl_bonds_all.append(bi)
                bc[bi] = _rgb(col)
                for ai in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
                    if ai not in ac:
                        hl_atoms.append(ai)
                    ac[ai] = _rgb(col)
        hl_bonds = hl_bonds_all
    d.DrawMolecule(mol, highlightAtoms=hl_atoms, highlightAtomColors=ac, highlightBonds=hl_bonds,
                   highlightBondColors=bc, legend=legend)
    d.FinishDrawing()
    return d.GetDrawingText()


def legend_html(items):
    chips = "".join(
        f'<span style="display:inline-block;background:{c};color:#111;padding:2px 9px;margin:2px 4px 2px 0;'
        f'border-radius:10px;font-size:0.85em">{label}</span>' for label, c in items)
    return chips
