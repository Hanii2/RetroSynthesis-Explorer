"""Compounds treated as 'available' (commercial / teaching-lab chemicals), with common names for display."""
from rdkit import Chem

from .chem import canon

NAMES = {
    # aromatics
    "benzene": "c1ccccc1", "toluene": "Cc1ccccc1", "phenol": "Oc1ccccc1", "aniline": "Nc1ccccc1",
    "nitrobenzene": "O=[N+]([O-])c1ccccc1", "benzaldehyde": "O=Cc1ccccc1", "benzoic acid": "OC(=O)c1ccccc1",
    "benzyl alcohol": "OCc1ccccc1", "benzyl bromide": "BrCc1ccccc1", "bromobenzene": "Brc1ccccc1",
    "chlorobenzene": "Clc1ccccc1", "acetophenone": "CC(=O)c1ccccc1", "benzoyl chloride": "O=C(Cl)c1ccccc1",
    "benzophenone": "O=C(c1ccccc1)c1ccccc1", "salicylic acid": "OC(=O)c1ccccc1O", "salicylaldehyde": "O=Cc1ccccc1O",
    "methyl salicylate": "COC(=O)c1ccccc1O", "4-aminophenol": "Nc1ccc(O)cc1", "4-nitrophenol": "O=[N+]([O-])c1ccc(O)cc1",
    "4-aminobenzoic acid (PABA)": "Nc1ccc(C(=O)O)cc1", "4-nitrobenzoic acid": "O=C(O)c1ccc([N+](=O)[O-])cc1",
    "4-nitrobenzoyl chloride": "O=C(Cl)c1ccc([N+](=O)[O-])cc1", "1-naphthol": "Oc1cccc2ccccc12",
    "2-nitrobenzaldehyde": "O=Cc1ccccc1[N+](=O)[O-]", "2,6-dimethylaniline": "Cc1cccc(C)c1N",
    "3-methylbenzoic acid": "Cc1cccc(C(=O)O)c1", "isobutylbenzene": "CC(C)Cc1ccccc1",
    "phenylboronic acid": "OB(O)c1ccccc1", "4-chlorobenzotrifluoride": "Clc1ccc(C(F)(F)F)cc1",
    "4-hydroxycoumarin": "O=C1C=C(O)c2ccccc2O1", "benzalacetone": "CC(=O)C=Cc1ccccc1",
    "acetanilide": "CC(=O)Nc1ccccc1", "anisole": "COc1ccccc1", "4-hydroxyacetophenone": "CC(=O)c1ccc(O)cc1",
    "1-(4-isobutylphenyl)ethanone": "CC(C)Cc1ccc(C(C)=O)cc1",
    # aliphatics, reagents, small building blocks
    "methanol": "CO", "ethanol": "CCO", "1-propanol": "CCCO", "isopropanol": "CC(C)O", "1-butanol": "CCCCO",
    "tert-butanol": "CC(C)(C)O", "isoamyl alcohol": "CC(C)CCO", "cyclohexanol": "OC1CCCCC1", "cyclohexanone": "O=C1CCCCC1",
    "acetone": "CC(C)=O", "2-butanone": "CCC(C)=O", "acetaldehyde": "CC=O", "formaldehyde": "C=O",
    "acetic acid": "CC(=O)O", "acetic anhydride": "CC(=O)OC(C)=O", "acetyl chloride": "CC(=O)Cl",
    "chloroacetyl chloride": "ClCC(=O)Cl", "ethyl acetate": "CCOC(C)=O", "methyl acetoacetate": "COC(=O)CC(C)=O",
    "malonic acid": "OC(=O)CC(=O)O", "diethylamine": "CCNCC", "dimethylamine": "CNC", "methylamine": "CN",
    "ammonia": "N", "isopropylamine": "CC(C)N", "tert-butylamine": "CC(C)(C)N",
    "N-benzyl-tert-butylamine": "CC(C)(C)NCc1ccccc1", "urea": "NC(N)=O", "2-(diethylamino)ethanol": "CCN(CC)CCO",
    "maleic anhydride": "O=C1C=CC(=O)O1", "1,3-butadiene": "C=CC=C", "cyclopentadiene": "C1=CC=CC1",
    "methyl vinyl ketone": "C=CC(C)=O", "epichlorohydrin": "ClCC1CO1", "epibromohydrin": "BrCC1CO1", "carbon dioxide": "O=C=O",
    "carbon monoxide": "[C-]#[O+]", "iodomethane": "CI", "iodoethane": "CCI", "bromoethane": "CCBr",
    "ethyl chloroacetate": "ClCC(=O)OCC", "2-methylcyclohexane-1,3-dione": "CC1C(=O)CCCC1=O",
    "1-bromobutane": "CCCCBr", "tert-butyl chloride": "CC(C)(C)Cl", "chloromethane": "CCl",
    "propene": "C=CC", "ethene": "C=C", "pyridine": "c1ccncc1", "diethyl ether": "CCOCC",
}
_CANON = {canon(s): n for n, s in NAMES.items()}
assert None not in _CANON


def name_of(smiles):
    return _CANON.get(canon(smiles))


def is_available(smiles):
    """Available = listed above, or a very small molecule (<= 7 heavy atoms), or a simple reagent-sized compound."""
    c = canon(smiles)
    if c in _CANON:
        return True
    m = Chem.MolFromSmiles(c)
    heavy = m.GetNumHeavyAtoms()
    if heavy <= 7:
        return True
    ri = m.GetRingInfo().NumRings()
    stereo = len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    reagent_like = (m.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[Cl,Br]")) or
                    m.HasSubstructMatch(Chem.MolFromSmarts("[SX4](=O)(=O)Cl")) or
                    m.HasSubstructMatch(Chem.MolFromSmarts("B(O)O")))
    return heavy <= 12 and ri <= 1 and stereo == 0 and (reagent_like or heavy <= 9)
