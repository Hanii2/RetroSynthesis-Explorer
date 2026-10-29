"""3D structure + properties (RDKit ETKDG + MMFF)."""
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors


class InputError(Exception):
    pass


def make_3d(smiles, seed=7):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise InputError("That SMILES is not valid.")
    if m.GetNumHeavyAtoms() > 120:
        raise InputError("Molecule too large for the 3D viewer (> 120 heavy atoms).")
    mol = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    if AllChem.EmbedMolecule(mol, ps) < 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) < 0:
            raise InputError("Could not build a 3D structure for this molecule.")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=5000)
    return mol


def molblock(smiles):
    return Chem.MolToMolBlock(make_3d(smiles))


def properties(smiles):
    m = Chem.MolFromSmiles(smiles)
    return {"Formula": rdMolDescriptors.CalcMolFormula(m), "Molecular weight": round(Descriptors.MolWt(m), 2),
            "Heavy atoms": m.GetNumHeavyAtoms(), "logP (Crippen)": round(Crippen.MolLogP(m), 2),
            "H-bond donors": Lipinski.NumHDonors(m), "H-bond acceptors": Lipinski.NumHAcceptors(m),
            "Rings": rdMolDescriptors.CalcNumRings(m)}
