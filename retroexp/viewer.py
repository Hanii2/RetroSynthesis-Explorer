"""3D viewer (py3Dmol) as an HTML string."""
import py3Dmol
from rdkit import Chem

from .colors import atom_colors

STYLES = {
    "Ball and stick": {"stick": {"radius": 0.14}, "sphere": {"scale": 0.28}},
    "Stick": {"stick": {"radius": 0.18}},
    "Space-filling": {"sphere": {"scale": 1.0}},
    "Wireframe": {"line": {}},
}
BACKGROUNDS = {"White": "#ffffff", "Light grey": "#eceff1", "Dark": "#1e272e"}


def viewer_html(molblock, style="Ball and stick", show_h=True, labels=False, width=620, height=460,
                spin=False, scheme="Vivid colours", background="White"):
    mol = Chem.MolFromMolBlock(molblock, removeHs=False)
    v = py3Dmol.view(width=width, height=height)
    v.addModel(molblock, "mol")
    groups = {}
    for idx, hexc in atom_colors(mol, scheme).items():
        groups.setdefault(hexc, []).append(idx)
    for hexc, idxs in groups.items():                    # one style per colour, atoms selected by index
        v.setStyle({"index": idxs}, {k: dict(val, color=hexc) for k, val in STYLES[style].items()})
    if not show_h:
        v.setStyle({"elem": "H"}, {})
    if labels:
        conf = mol.GetConformer()
        fc = "white" if background == "Dark" else "black"
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 1 and not show_h:
                continue
            p = conf.GetAtomPosition(a.GetIdx())
            v.addLabel(f"{a.GetSymbol()}{a.GetIdx() + 1}", {"position": {"x": p.x, "y": p.y, "z": p.z},
                                                             "fontSize": 11, "fontColor": fc, "backgroundOpacity": 0.0})
    v.setBackgroundColor(BACKGROUNDS.get(background, "#ffffff"))
    v.zoomTo()
    if spin:
        v.spin(True)
    return v._make_html()
