"""Name or SMILES -> SMILES.  Order: library name -> common-chemical name -> valid SMILES -> PubChem -> NCI CIR."""
from urllib.parse import quote

import requests
from rdkit import Chem

from . import library
from .buildingblocks import NAMES

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/IsomericSMILES,CanonicalSMILES/JSON"
CIR = "https://cactus.nci.nih.gov/chemical/structure/{}/smiles"
KEYS = ("SMILES", "IsomericSMILES", "CanonicalSMILES", "ConnectivitySMILES")


class ResolveError(Exception):
    pass


def resolve(text, timeout=20):
    """Returns (smiles, source)."""
    text = (text or "").strip()
    if not text:
        raise ResolveError("Type a compound name (e.g. paracetamol) or a SMILES.")
    e = library.find(text)
    if e and e["name"].lower() == text.lower():
        return e["smiles"], "built-in library"
    for n, smi in NAMES.items():
        if n.lower() == text.lower():
            return smi, "common chemicals list"
    if " " not in text and Chem.MolFromSmiles(text) is not None:
        return text, "SMILES"
    try:
        r = requests.get(PUBCHEM.format(quote(text)), timeout=timeout)
        if r.ok:
            props = r.json()["PropertyTable"]["Properties"][0]
            for k in KEYS:
                if props.get(k):
                    return max(props[k].split("."), key=len), "PubChem"
    except Exception:
        pass
    try:
        r = requests.get(CIR.format(quote(text)), timeout=timeout)
        if r.ok and r.text.strip() and Chem.MolFromSmiles(r.text.split()[0]):
            return r.text.split()[0], "NCI CIR"
    except Exception:
        pass
    raise ResolveError(f"Could not find '{text}'. Check the English spelling, pick from the quick list, or paste a SMILES "
                       f"(internet is needed for names that are not in the built-in list).")
