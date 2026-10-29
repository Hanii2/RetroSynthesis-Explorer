from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def to_mol(smi):
    return Chem.MolFromSmiles(smi) if smi else None


def canon(smi, stereo=False):
    m = to_mol(smi)
    if m is None:
        return None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m, isomericSmiles=stereo)
