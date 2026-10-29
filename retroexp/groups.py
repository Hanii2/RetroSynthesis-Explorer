"""Functional-group SMARTS used by the rule checker and the precedent signature."""
from rdkit import Chem

FG_SMARTS = {
    "alcohol": "[CX4][OX2H]",
    "phenol": "c[OX2H]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6][CX3](=O)[OX2H0][#6;!$(C=O)]",
    "amide": "[CX3](=O)[NX3]",
    "amide_NH": "[NX3;H1,H2][CX3]=O",
    "amine_aliphatic": "[NX3;H2,H1;!$(NC=[O,S,N]);!$(NC=C);!$(Na);!$(NS(=O)=O)]",
    "amine_aromatic": "[NX3;H2,H1;!$(NC=[O,S,N]);!$(NS(=O)=O)]a",
    "ketone": "[#6][CX3](=O)[#6]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "aryl_ketone": "c[CX3](=O)[#6]",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "nitrile": "[CX2]#[NX1]",
    "alkene": "[CX3]=[CX3]",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "aryl_halide": "c[Cl,Br,I]",
    "acyl_halide": "[CX3](=O)[Cl,Br,I]",
    "anhydride": "[CX3](=O)[OX2][CX3](=O)",
    "aryl_alkyl": "c[CX4]",
    # protecting groups
    "benzyl_N": "[NX3][CH2]c",
    "benzyl_O": "[OX2][CH2]c",
    "boc": "[N,O]C(=O)OC(C)(C)C",
    "silyl": "[Si]",
}
PROTECTING = ["benzyl_N", "benzyl_O", "boc", "silyl"]

# extra patterns (not counted, only used by rules)
EXTRA = {
    "arene_CH": "[cH]",
    "ewg_on_arene": "c-[$([N+](=O)[O-]),$([NX3](=O)=O),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O)]",
    "activator_on_arene": "c-[$([OX2H]),$([OX2][#6]),$([NX3;!$(N=O);!$(NC=O);!$([N+])]),$([CX4])]",
    "isobutyl_halide": "[CH2]([Cl,Br,I])[CH]([CH3])[CH3]",
    "neopentyl_halide": "[CH2]([Cl,Br,I])C([CH3])([CH3])[CH3]",
    "long_primary_halide": "[CH2;!$(Cc)]([Cl,Br,I])[CH2][CX4]",
}

_P = {k: Chem.MolFromSmarts(v) for k, v in {**FG_SMARTS, **EXTRA}.items()}
for _k, _v in _P.items():
    assert _v is not None, f"bad SMARTS: {_k}"


def has(mol, name):
    return mol.HasSubstructMatch(_P[name])


def count_groups(mol):
    return {k: len(mol.GetSubstructMatches(_P[k])) for k in FG_SMARTS}
