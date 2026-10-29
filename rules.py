"""
Chemistry-rule validator  (the "textbook check").

AiZynthFinder only gives reactants -> product (no reagents), so every rule works on the
NET functional-group change of a step plus the groups that are "spectators" in the reactants.

To add your own rule (recommended for the report!):

    from retro.rules import register

    @register("my_rule")
    def my_rule(c):                      # c = StepContext (see below)
        if c.d["nitrile"] < 0 and c.R["alkene"] > 0:
            return [("warn", "nitrile reduction with an alkene present ...")]
        return []

Each rule returns a list of (severity, message); severity = "fail" | "warn" | "info".
"""
from dataclasses import dataclass

from rdkit import Chem

from .chem import to_mol
from .groups import FG_SMARTS, PROTECTING, count_groups, has

PENALTY = {"fail": 0.40, "warn": 0.12, "info": 0.0}
RULES = {}


def register(name):
    def deco(fn):
        RULES[name] = fn
        return fn
    return deco


@dataclass
class Finding:
    step: int
    severity: str
    rule: str
    message: str


class StepContext:
    def __init__(self, idx, reactants, product, target_has_stereo=False):
        self.idx, self.reactants, self.product = idx, reactants, product
        self.target_has_stereo = target_has_stereo
        self.rmols = [to_mol(r) for r in reactants]
        self.pmol = to_mol(product)
        self.valid = self.pmol is not None and all(m is not None for m in self.rmols)
        if not self.valid:
            return
        self.rcounts = [count_groups(m) for m in self.rmols]
        self.R = {k: sum(c[k] for c in self.rcounts) for k in FG_SMARTS}
        self.P = count_groups(self.pmol)
        self.d = {k: self.P[k] - self.R[k] for k in FG_SMARTS}

    def reactants_with(self, group):
        return [m for m, c in zip(self.rmols, self.rcounts) if c[group] > 0]


def _n_stereocenters(mol):
    return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))


# ------------------------------------------------------------------ rules
@register("reduction_chemoselectivity")
def reduction_chemoselectivity(c):
    d, R = c.d, c.R
    organometallic = d["alkyl_halide"] + d["aryl_halide"] < 0
    reduced = []
    if d["ketone"] + d["aldehyde"] < 0 and d["alcohol"] > 0 and not organometallic:
        reduced.append("ketone/aldehyde")
    if d["nitro"] < 0 and (d["amine_aromatic"] > 0 or d["amine_aliphatic"] > 0):
        reduced.append("nitro")
    if d["nitrile"] < 0 and d["amine_aliphatic"] > 0:
        reduced.append("nitrile")
    if d["ester"] < 0 and d["alcohol"] > 0 and d["carboxylic_acid"] <= 0 and not organometallic:
        reduced.append("ester")
    if d["carboxylic_acid"] < 0 and d["alcohol"] > 0 and not organometallic:
        reduced.append("acid")
    if not reduced:
        return []
    risk = {
        "ketone/aldehyde": ["ester", "nitrile", "nitro", "carboxylic_acid", "amide"],
        "nitro": ["alkene", "nitrile", "benzyl_N", "benzyl_O", "aryl_halide"],
        "nitrile": ["ester", "ketone", "aldehyde", "nitro", "amide"],
        "ester": ["ketone", "aldehyde", "nitrile", "nitro", "amide"],
        "acid": ["ketone", "aldehyde", "nitrile", "nitro", "amide", "ester"],
    }
    out = []
    for kind in reduced:
        survivors = [g for g in risk[kind] if R[g] > 0 and d[g] >= 0]
        if survivors:
            out.append(("warn", f"{kind} is reduced while {', '.join(survivors)} must survive: "
                                f"needs a chemoselective reagent (e.g. NaBH4 vs LiAlH4, Fe/HCl vs H2/Pd)."))
        else:
            out.append(("info", f"{kind} reduction; no sensitive group left untouched."))
    return out


@register("acylation_selectivity")
def acylation_selectivity(c):
    d, R, P = c.d, c.R, c.P
    donor_used = d["acyl_halide"] < 0 or d["anhydride"] < 0 or d["carboxylic_acid"] < 0
    if not donor_used or not (d["amide"] > 0 or d["ester"] > 0):
        return []
    out = []
    amine_in = R["amine_aliphatic"] + R["amine_aromatic"]
    amine_left = P["amine_aliphatic"] + P["amine_aromatic"]
    oh_in = R["alcohol"] + R["phenol"]
    oh_left = P["alcohol"] + P["phenol"]
    if d["ester"] > 0 and d["amide"] <= 0 and amine_in > 0 and amine_left > 0:
        out.append(("fail", "O-acylation drawn while a free amine is present: the amine is the better "
                            "nucleophile, so it must be protected (or the step is wrong)."))
    if d["amide"] > 0 and oh_in > 0 and oh_left > 0:
        out.append(("info", "N-acylation in presence of free OH: chemoselective (amine > OH), textbook-OK."))
    if d["ester"] > 0 and oh_in >= 2 and oh_left >= 1:
        out.append(("warn", "several OH groups, only one acylated: regioselectivity must be controlled."))
    return out


@register("friedel_crafts_acylation")
def friedel_crafts_acylation(c):
    if not (c.d["aryl_ketone"] > 0 and (c.d["acyl_halide"] < 0 or c.d["anhydride"] < 0)):
        return []
    out = []
    for m, cnt in zip(c.rmols, c.rcounts):
        if cnt["acyl_halide"] or cnt["anhydride"] or not has(m, "arene_CH"):
            continue
        ewg, act = has(m, "ewg_on_arene"), has(m, "activator_on_arene")
        if ewg and not act:
            out.append(("fail", "Friedel-Crafts acylation on a deactivated ring (electron-withdrawing group, "
                                "no activating group): this reaction does not work."))
        elif ewg and act:
            out.append(("info", "Friedel-Crafts on a ring with both EWG and activator: check regiochemistry."))
        if cnt["amine_aromatic"] or cnt["amine_aliphatic"]:
            out.append(("warn", "free amine complexes AlCl3 / gets acylated: protect the amine first."))
        if cnt["phenol"]:
            out.append(("warn", "free phenol: O-acylation / Fries rearrangement competes with C-acylation."))
    return out


@register("friedel_crafts_alkylation")
def friedel_crafts_alkylation(c):
    if not (c.d["aryl_alkyl"] > 0 and c.d["alkyl_halide"] < 0):
        return []
    out = []
    for m in c.reactants_with("alkyl_halide"):
        if has(m, "isobutyl_halide") or has(m, "neopentyl_halide"):
            out.append(("fail", "Friedel-Crafts alkylation with an isobutyl/neopentyl halide: carbocation "
                                "rearrangement gives the tert-alkyl product, not the drawn one."))
        elif has(m, "long_primary_halide"):
            out.append(("warn", "primary alkyl halide (>= 3 carbons) in Friedel-Crafts alkylation: "
                                "rearrangement / polyalkylation likely; acylation then reduction is safer."))
    return out


@register("organometallic_acidic_protons")
def organometallic_acidic_protons(c):
    d = c.d
    is_cc = (d["alkyl_halide"] + d["aryl_halide"] < 0 and d["alcohol"] > 0
             and d["ketone"] + d["aldehyde"] + d["ester"] < 0)
    if not is_cc:
        return []
    acidic = ["alcohol", "phenol", "carboxylic_acid", "amine_aliphatic", "amine_aromatic", "amide_NH"]
    out = []
    for cnt in c.rcounts:                       # ANY partner: R-MgBr is a strong base and is quenched by O-H / N-H
        bad = [g for g in acidic if cnt[g] > 0]
        if bad:
            who = "halide partner" if cnt["alkyl_halide"] + cnt["aryl_halide"] > 0 else "carbonyl partner"
            out.append(("fail", f"Grignard/organolithium-type C-C bond formation but the {who} carries acidic protons "
                                f"({', '.join(bad)}): the organometallic would be quenched; protect the group first."))
    return out


@register("stereochemistry")
def stereochemistry(c):
    n_p = _n_stereocenters(c.pmol)
    n_r = sum(_n_stereocenters(m) for m in c.rmols)
    if n_p > n_r:
        if c.target_has_stereo:
            return [("warn", f"{n_p - n_r} new stereocentre(s) and the target is stereodefined: "
                             f"the step needs chiral control (catalyst, auxiliary or resolution).")]
        return [("info", f"{n_p - n_r} new stereocentre(s): product is racemic unless chiral control is used.")]
    return []


@register("protecting_groups")
def protecting_groups(c):
    d, P, R = c.d, c.P, c.R
    other = sum(abs(v) for k, v in d.items() if k not in PROTECTING)
    pg = sum(abs(d[k]) for k in PROTECTING)
    out = []
    if pg and not other:
        out.append(("info", "protecting-group manipulation only (counts as an extra step)."))
    if d["benzyl_N"] + d["benzyl_O"] < 0:
        sens = [g for g in ["alkene", "nitro", "nitrile", "aryl_halide"] if P[g] > 0]
        if sens:
            out.append(("warn", f"benzyl removal usually uses H2/Pd, which would also affect {', '.join(sens)}."))
    return out


@register("atom_bookkeeping")
def atom_bookkeeping(c):
    def elems(m):
        return {a.GetSymbol() for a in m.GetAtoms()}
    new = elems(c.pmol) - set().union(*[elems(m) for m in c.rmols])
    if new:
        return [("warn", f"element(s) {', '.join(sorted(new))} appear in the product but not in the "
                         f"reactants: they must come from a reagent (check it).")]
    heavy_p = c.pmol.GetNumHeavyAtoms()
    heavy_r = sum(m.GetNumHeavyAtoms() for m in c.rmols)
    if heavy_p > heavy_r + 4:
        return [("warn", "product has far more heavy atoms than the reactants: a reagent contributes atoms.")]
    return []


# ------------------------------------------------------------------ driver
def check_step(idx, step, target_has_stereo=False):
    c = StepContext(idx, step.reactants, step.product, target_has_stereo)
    if not c.valid:
        return [Finding(idx, "fail", "valid_smiles", "a SMILES in this step could not be parsed")]
    out = []
    for name, fn in RULES.items():
        for sev, msg in fn(c):
            out.append(Finding(idx, sev, name, msg))
    return out


def check_route(route, target_has_stereo=False):
    """Returns (rule_score 0..1, findings list)."""
    findings = []
    for i, step in enumerate(route.steps, 1):
        findings += check_step(i, step, target_has_stereo)
    penalty = sum(PENALTY[f.severity] for f in findings)
    return max(0.0, 1.0 - penalty), findings


# ------------------------------------------------------------------ rules added for the retrosynthesis explorer
@register("organometallic_incompatible_groups")
def organometallic_incompatible_groups(c):
    """Grignard / organolithium partner (a halide) must not carry groups the organometallic would attack."""
    d = c.d
    if not (d["alkyl_halide"] + d["aryl_halide"] < 0 and (d["alcohol"] > 0 or d["carboxylic_acid"] > 0)):
        return []
    electrophiles = ["ester", "ketone", "aldehyde", "nitrile", "nitro", "amide", "acyl_halide", "anhydride"]
    out = []
    for cnt in c.rcounts:
        if cnt["alkyl_halide"] + cnt["aryl_halide"] > 0:
            bad = [g for g in electrophiles if cnt[g] > 0]
            if bad:
                out.append(("fail", f"The R-MgBr/R-Li formed from this halide would attack the {', '.join(bad)} in the same "
                                    f"molecule: choose a different disconnection or protect it."))
    return out


@register("amino_acid_chloride")
def amino_acid_chloride(c):
    out = []
    for cnt in c.rcounts:
        if cnt["acyl_halide"] and (cnt["amine_aliphatic"] or cnt["amine_aromatic"]):
            out.append(("warn", "an acid chloride with a free N-H amine in the same molecule would self-condense."))
    sulfonyl = Chem.MolFromSmarts("[SX4](=O)(=O)Cl")
    for m, cnt in zip(c.rmols, c.rcounts):
        if m.HasSubstructMatch(sulfonyl) and (cnt["amine_aliphatic"] or cnt["amine_aromatic"]):
            out.append(("fail", "a sulfonyl chloride with a free NH2/NHR in the same molecule polymerises: protect the amine "
                                "(e.g. as an acetamide) before chlorosulfonation."))
        if cnt["acyl_halide"] and m.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H0;!$(NC=O)]CC(=O)Cl")):
            out.append(("warn", "an alpha-amino acid chloride is unstable and awkward to handle; the alpha-halo acid chloride then "
                                "SN2 with the amine is the usual order."))
    return out
