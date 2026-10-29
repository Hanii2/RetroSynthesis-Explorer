"""
Scope guard. This is a teaching tool for ordinary organic chemistry; it does not give routes to controlled drugs,
explosives or chemical-weapon agents. Best-effort filter (name list + structure patterns); not exhaustive.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

NAME_BLOCK = {
    "methamphetamine", "amphetamine", "mdma", "ecstasy", "heroin", "cocaine", "fentanyl", "lsd", "psilocybin",
    "psilocin", "dmt", "thc", "ketamine", "pcp", "ghb", "methadone", "oxycodone", "hydrocodone", "codeine",
    "morphine", "ephedrine", "pseudoephedrine", "phenylacetone", "p2p", "safrole", "mescaline", "cathinone",
    "mephedrone", "alprazolam", "diazepam", "tnt", "rdx", "hmx", "petn", "nitroglycerin", "tatp", "hmtd",
    "sarin", "soman", "tabun", "vx", "mustard gas", "sulfur mustard", "nitrogen mustard", "lewisite",
    "phosgene", "picric acid", "trinitrotoluene", "cyclonite",
}

_PATTERNS = {
    "amphetamine / cathinone-type stimulant": ["[CX4;H1]([CH3])([NX3;!$(NC=[O,S])])[CH2]c1ccccc1",
                                              "[CX4;H1]([CH3])([NX3;!$(NC=[O,S])])C(=O)c1ccccc1"],
    "fentanyl-type opioid": ["[#6]C(=O)N(c1ccccc1)C1CCN(CC1)"],
    "tropane alkaloid (cocaine-type)": ["C1CC2CCC(C1)N2"],
    "benzodiazepine": ["c1ccc2c(c1)C(c1ccccc1)=N[#6][#6](=O)N2"],
    "tryptamine hallucinogen (N,N-dialkyl)": ["[CH2]([CH2][NX3;H0;!$(NC=O)]([CH3,CH2])[CH3,CH2])c1c[nH]c2ccccc12"],
    "cannabinoid": ["c1cc2OC(C)(C)C3CCC(C)=CC3c2c(O)c1"],
    "organophosphorus nerve agent": ["[#15](=O)F", "[#15](F)", "[#15](=O)(C#N)", "[#15](=O)([#8])([#6])S[CH2][CH2][NX3]"],
    "mustard agent": ["[NX3](CC[Cl,Br,I])CC[Cl,Br,I]", "ClCCSCCCl", "[Cl,Br]CC[SX2]CC[Cl,Br]"],
    "arylcyclohexylamine (ketamine / PCP-type)": ["[NX3]C1(CCCCC1=O)c1ccccc1", "C1(CCCCC1)(c1ccccc1)N1CCCCC1"],
    "arsenical": ["[As]"],
    "phosgene / cyanogen halide": ["O=C(Cl)Cl", "[Cl,Br]C#N"],
}
_COMPILED = {k: [Chem.MolFromSmarts(s) for s in v] for k, v in _PATTERNS.items()}

# scaffolds derived from morphine / LSD (generic atoms + bonds so that substituents do not matter)
def _scaffold_query(smiles):
    from rdkit.Chem.Scaffolds import MurckoScaffold
    m = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
    p = Chem.AdjustQueryParameters()
    p.makeAtomsGeneric = True
    p.makeBondsGeneric = True
    p.adjustDegree = False
    return Chem.AdjustQueryProperties(m, p)


_SCAFFOLDS = {
    "opioid (morphinan)": _scaffold_query("CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"),
    "ergoline (LSD-type)": _scaffold_query("CCN(CC)C(=O)C1C=C2c3cccc4[nH]cc(c34)CC2N(C1)C"),
}


def check(smiles, text=""):
    """Returns (blocked: bool, reason: str)."""
    if text and text.strip().lower() in NAME_BLOCK:
        return True, "controlled substance, explosive or chemical-weapon agent"
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False, ""
    for label, patts in _COMPILED.items():
        if any(m.HasSubstructMatch(p) for p in patts):
            return True, label
    for label, q in _SCAFFOLDS.items():
        if m.HasSubstructMatch(q):
            return True, label
    nitro = len(m.GetSubstructMatches(Chem.MolFromSmarts("[N+](=O)[O-]")))
    nitrate = len(m.GetSubstructMatches(Chem.MolFromSmarts("[CX4]O[N+](=O)[O-]")))
    nitramine = len(m.GetSubstructMatches(Chem.MolFromSmarts("[NX3][N+](=O)[O-]")))
    if nitro >= 3 or nitrate >= 2 or nitramine >= 2:
        return True, "polynitro / nitrate-ester / nitramine (explosive-type) compound"
    if m.HasSubstructMatch(Chem.MolFromSmarts("[OX2][OX2]")):
        return True, "organic peroxide (explosive-type) compound"
    return False, ""


MESSAGE = ("This tool is for ordinary teaching molecules. It does not give synthesis routes for controlled drugs, "
           "explosives or chemical-weapon agents ({reason}). Try another compound.")


def name_blocked(text):
    return (text or "").strip().lower() in NAME_BLOCK
