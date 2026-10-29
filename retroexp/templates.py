"""
Retro-synthetic disconnection rules (reaction SMARTS written in the retro direction: target >> precursors).

Convention: precursors are the organic building blocks that END UP in the product (substrates, alkyl halides,
acid chlorides, amines, CO2 ...). Catalysts, acids/bases, oxidants and reductants are written in `reagents`.
`weight` = how typical/reliable the textbook reaction is (0..1). `kind` is C-C, C-X, FGI (functional group
interconversion), Ring or Protecting. Stereochemistry is NOT handled (it is removed before matching).
"""
from dataclasses import dataclass
from typing import Callable, Optional

from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class Template:
    id: str
    name: str
    kind: str
    retro: str
    reagents: str
    tip: str
    weight: float
    check: Optional[Callable] = None

    def __post_init__(self):
        self.rxn = AllChem.ReactionFromSmarts(self.retro)
        assert self.rxn is not None, self.id


def _two_identical(prods):          # for "ester + 2 R-MgBr": the two halides must be the same
    smi = [Chem.MolToSmiles(p) for p in prods[1:]]
    return len(set(smi)) == 1


NOT_ACYL = "!$(C=[O,S,N])"

TEMPLATES = [
    # ------------------------------------------------------------ C-X bond formation
    Template("fischer", "Fischer esterification", "C-X",
             "[#6:5][C:1](=[O:2])[O:3]-!@[CX4:4]>>[#6:5][C:1](=[O:2])O.[O:3][C:4]",
             "H2SO4 (cat.), heat; use excess alcohol or remove water (Dean-Stark)",
             "Disconnect the C(acyl)-O bond: carboxylic acid + alcohol. An equilibrium, so push it with excess alcohol.", 0.85),
    Template("acetylate_aryl", "O-Acetylation of a phenol", "C-X",
             "[CH3:5][C:1](=[O:2])[O:3]-!@[c:4]>>[CH3:5][C:1](=[O:2])OC(=O)C.[O:3][c:4]",
             "Acetic anhydride, H2SO4 or H3PO4 (cat.), 50-60 C",
             "Phenols are poor nucleophiles for Fischer esterification, so use the more reactive anhydride.", 0.85),
    Template("acylate_aryl", "Esterification of a phenol (acid chloride)", "C-X",
             "[#6:5][C:1](=[O:2])[O:3]-!@[c:4]>>[#6:5][C:1](=[O:2])Cl.[O:3][c:4]",
             "Pyridine or Et3N (base), 0 C to rt",
             "Acid chloride + phenol with a base to mop up HCl.", 0.65),
    Template("acylate_alkyl", "Ester from acid chloride + alcohol", "C-X",
             "[#6:5][C:1](=[O:2])[O:3]-!@[CX4:4]>>[#6:5][C:1](=[O:2])Cl.[O:3][C:4]",
             "Pyridine or Et3N, CH2Cl2, 0 C",
             "Faster and irreversible compared with Fischer, but needs the acid chloride.", 0.55),
    Template("acetylate_amine", "N-Acetylation of an amine", "C-X",
             "[CH3:5][C:1](=[O:2])-!@[NX3;!$(N-[N,O]);!$(N(-C=O)-C=O):3]>>[CH3:5][C:1](=[O:2])OC(=O)C.[N:3]",
             "Acetic anhydride, (NaOAc or pyridine), rt to 60 C",
             "Amide bond: disconnect C(=O)-N. The amine is far more nucleophilic than a phenol OH, so it is acylated selectively.", 0.9),
    Template("amide_chloride", "Amide from acid chloride + amine", "C-X",
             "[#6:5][C:1](=[O:2])-!@[NX3;!$(N-[N,O]);!$(N(-C=O)-C=O):3]>>[#6:5][C:1](=[O:2])Cl.[N:3]",
             "Et3N or aq. NaOH (Schotten-Baumann), 0 C",
             "Disconnect the amide C-N bond. Acid chloride + amine with a base to neutralise HCl.", 0.8),
    Template("sulfonamide", "Sulfonamide from sulfonyl chloride + amine", "C-X",
             "[#6:5][S:1](=[O:2])(=[O:4])-!@[NX3:3]>>[#6:5][S:1](=[O:2])(=[O:4])Cl.[N:3]",
             "Pyridine or excess amine, 0 C to rt", "Same logic as an amide, with S in place of C.", 0.8),
    Template("williamson", "Williamson ether synthesis", "C-X",
             "[#6;!$(C=[O,S,N]):1]-!@[OX2:2]-!@[CX4;H2,H3:3]>>[#6:1][O:2].Br[C:3]",
             "1) NaH or K2CO3 (make the alkoxide/phenoxide); 2) alkyl halide (SN2)",
             "Choose the PRIMARY (or methyl) halide: SN2 on a tertiary halide gives elimination instead.", 0.85),
    Template("amine_alkyl_activated", "N-Alkylation with an alpha-halo carbonyl", "C-X",
             "[NX3;!$(N-C=[O,S,N]);!$(N-S):1]-!@[CH2:2]-!@[C:3]=[O:4]>>[N:1].Cl[CH2:2][C:3]=[O:4]",
             "Excess amine or K2CO3, toluene or MeCN, heat",
             "alpha-Halo carbonyls are very good SN2 electrophiles.", 0.8),
    Template("amine_alkyl_benzylic", "N-Alkylation with a benzylic halide", "C-X",
             "[NX3;!$(N-C=[O,S,N]);!$(N-S):1]-!@[CH2:2]-!@[c:3]>>[N:1].Br[CH2:2][c:3]",
             "K2CO3, MeCN or DMF", "Benzylic halides are reactive SN2 electrophiles.", 0.65),
    Template("amine_alkyl", "N-Alkylation (SN2)", "C-X",
             "[NX3;!$(N-C=[O,S,N]);!$(N-S):1]-!@[CX4;H2,H3:2]>>[N:1].Br[C:2]",
             "Excess amine or K2CO3, heat",
             "Warning: over-alkylation is common (1 degree -> 2 -> 3 -> quaternary). Reductive amination is cleaner.", 0.45),
    Template("reductive_amination", "Reductive amination", "C-X",
             "[CX4;H1,H2;!$(C([#7,#8,#16])[#7,#8,#16]):1]-!@[NX3;!$(N-C=[O,S,N]);!$(N-S):2]>>[C:1]=O.[N:2]",
             "NaBH3CN or NaBH(OAc)3, mild acid (imine/iminium then reduction)",
             "Carbonyl + amine gives an imine that is reduced in situ; avoids over-alkylation.", 0.7),
    Template("epoxide_opening", "Epoxide ring opening by an amine", "C-X",
             "[OH:1][CH:2]([#6:3])[CH2:4]-!@[NX3;!$(N-C=[O,S]):5]>>[C:2]1([#6:3])[CH2:4][O:1]1.[N:5]",
             "Amine (excess), EtOH or iPrOH, heat",
             "The amine attacks the less hindered epoxide carbon (SN2), giving a beta-amino alcohol.", 0.8),
    Template("snar", "Nucleophilic aromatic substitution (SNAr)", "C-X",
             "[CX4:2][OX2:1][c:3]1[cH:4][cH:5][c:6](-[$(C(F)(F)F),$([N+](=O)[O-]),$(C#N):9])[cH:7][cH:8]1>>"
             "[CX4:2][OH:1].Cl[c:3]1[cH:4][cH:5][c:6](-[*:9])[cH:7][cH:8]1",
             "NaH (alkoxide), DMSO or DMF, heat",
             "Works because an electron-withdrawing group (CF3, NO2, CN) para to the leaving group stabilises the Meisenheimer intermediate.", 0.6),
    # ------------------------------------------------------------ C-C bond formation
    Template("grignard_carbonyl", "Grignard addition to a carbonyl", "C-C",
             "[OH:1][CX4;H0,H1;!$(C([OX2,N,S])[OX2,N,S]):2](-!@[#6:3])-[#6:4]>>[O:1]=[C:2]-[#6:4].Br-[#6:3]",
             "1) Mg, dry ether/THF (forms R-MgBr); 2) the aldehyde/ketone; 3) H3O+ work-up",
             "Break a C-C bond next to the C-OH: the carbonyl partner + an organomagnesium (from R-Br). The substrate must have NO acidic protons (OH, NH, COOH).", 0.8),
    Template("grignard_ester", "Double Grignard addition to an ester", "C-C",
             "[OH:1][CX4;H0:2](-!@[#6:3])(-!@[#6:4])-!@[#6:5]>>[#6:3][C:2](=[O:1])O[CH3].Br[#6:4].Br[#6:5]",
             "1) 2 equiv. R-MgBr, dry ether; 2) H3O+",
             "An ester reacts with TWO equivalents of Grignard (via a ketone), giving a tertiary alcohol with two identical R groups.", 0.8,
             check=_two_identical),
    Template("grignard_co2", "Carboxylation of a Grignard reagent", "C-C",
             "[c:3]-!@[C:1](=[O:2])[OH:4]>>Br[c:3].[O:2]=[C:1]=[O:4]",
             "1) Mg, dry ether; 2) CO2 (dry ice); 3) H3O+", "Aryl-MgBr + CO2 gives an aryl carboxylic acid (one extra carbon).", 0.75),
    Template("fc_acylation", "Friedel-Crafts acylation", "C-C",
             "[c:1]-!@[C:2](=[O:3])-[#6:4]>>[c:1].Cl[C:2](=[O:3])[#6:4]",
             "AlCl3 (stoichiometric), CH2Cl2 or CS2 (acid chloride or anhydride)",
             "Aryl ketone: disconnect Ar-C(=O). Fails on deactivated rings (NO2, C=O). The product ketone does not react again.", 0.8),
    Template("fc_acetylation_anhydride", "Friedel-Crafts acetylation (anhydride)", "C-C",
             "[c:1]-!@[C:2](=[O:3])-[CH3:4]>>[c:1].[CH3:4][C:2](=[O:3])OC(=O)C",
             "AlCl3 (2+ equiv.) or HF, acetic anhydride", "Acetic anhydride is cheaper and easier to handle than acetyl chloride.", 0.8),
    Template("fc_alkylation", "Friedel-Crafts alkylation", "C-C",
             "[c:1]-!@[CX4;!$(C-[!#6;!#1]):2]>>[c:1].Cl[C:2]",
             "AlCl3 or FeCl3",
             "Risky: carbocation rearrangement and poly-alkylation. Prefer acylation then reduction.", 0.3),
    Template("kolbe_schmitt", "Kolbe-Schmitt carboxylation", "C-C",
             "[c:1]([OH:5])[c:6]-!@[C:2](=[O:3])[OH:4]>>[c:1]([OH:5])[c:6].[O:3]=[C:2]=[O:4]",
             "1) NaOH (phenoxide); 2) CO2, 125 C, ~100 atm; 3) H+", "Phenoxide + CO2 puts the COOH ortho to the OH (salicylic acid).", 0.75),
    Template("suzuki", "Suzuki coupling", "C-C",
             "[c:1]-!@[c:2]>>[c:1]B(O)O.Br[c:2]",
             "Pd(PPh3)4 (cat.), K2CO3 or Na2CO3, toluene/EtOH/H2O, heat", "Biaryl bond = aryl boronic acid + aryl halide with a Pd catalyst.", 0.85),
    Template("wittig", "Wittig olefination", "C-C",
             "[#6:1][CH:2]=!@[CH:3][#6:4]>>[#6:1][CH:2]=O.Br[CH2:3][#6:4]",
             "1) PPh3 (makes the phosphonium salt from the bromide); 2) strong base (ylide); 3) the aldehyde",
             "The C=C bond is built exactly where the C=O was; no mixture of double-bond positions.", 0.75),
    Template("aldol_condensation", "Aldol (Claisen-Schmidt) condensation", "C-C",
             "[c:1][CH:2]=[CH:3]-!@[C:4](=[O:5])[#6:6]>>[c:1][CH:2]=O.[C:3][C:4](=[O:5])[#6:6]",
             "NaOH or KOH, EtOH/H2O, rt",
             "An enone = aromatic aldehyde (no alpha-H, cannot self-condense) + a ketone with alpha-H, then loss of water.", 0.85),
    Template("knoevenagel", "Knoevenagel-Doebner condensation", "C-C",
             "[c:1][CH:2]=[CH:3]-!@[C:4](=[O:5])[OH:6]>>[c:1][CH:2]=O.[OH:6][C:4](=[O:5])[CH2:3]C(=O)O",
             "Pyridine, piperidine (cat.), heat (decarboxylates)", "Aldehyde + malonic acid gives the cinnamic acid after decarboxylation.", 0.8),
    Template("aldol", "Aldol addition", "C-C",
             "[OH:1][CH:2]([#6:3])-!@[CH2:4][C:5]=[O:6]>>[O:1]=[CH:2][#6:3].[C:4][C:5]=[O:6]",
             "NaOH (cat.), cold", "beta-Hydroxy carbonyl: disconnect the Ca-Cb bond.", 0.7),
    Template("claisen", "Claisen condensation", "C-C",
             "[#6:1][C:2](=[O:3])[CH2:4]-!@[C:5](=[O:6])[O:7][CH2:8][CH3:9]>>[#6:1][C:2](=[O:3])OCC.[C:4][C:5](=[O:6])[O:7][CH2:8][CH3:9]",
             "NaOEt, EtOH, then H3O+", "beta-Keto ester = two esters; the alkoxide must match the ester (OEt with ethyl esters).", 0.7),
    Template("mannich", "Mannich reaction", "C-C",
             "[C:1](=[O:2])[CH2:3]-!@[CH2:4]-!@[NX3;!$(N-C=[O,S]):5]>>[C:1](=[O:2])[C:3].[C:4]=O.[N:5]",
             "Amine hydrochloride, formaldehyde (paraformaldehyde), EtOH, H+, heat",
             "beta-Amino ketone = ketone + formaldehyde + amine (three-component).", 0.65),
    Template("michael", "Michael addition", "C-C",
             "[#6;$([#6]-,:[#6]=O):1]-!@[CX4;H1,H2:2]-!@[CH2:3]-!@[C:4](=[O:5])[#6:6]>>[C:1].[C:2]=[C:3]-[C:4](=[O:5])[#6:6]",
             "Base (NaOEt, Et3N, pyrrolidine), rt", "Nucleophilic carbon (stabilised by C=O) adds to the beta-carbon of an enone (1,4-addition).", 0.65),
    Template("diels_alder", "Diels-Alder reaction", "Ring",
             "[CX4:1]1[CX3:2]=[CX3:3][CX4:4][CX4:5][CX4:6]1>>[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]",
             "Heat (or Lewis acid); diene + dienophile, [4+2]",
             "A cyclohexene is a diene + dienophile. The ring double bond marks the former diene centre.", 0.8),
    Template("hantzsch", "Hantzsch dihydropyridine synthesis", "Ring",
             "[#6:7][CH:1]1[C:2]([C:14](=[O:15])[O:16][#6:17])=[C:3]([CH3:8])[NH:4][C:5]([CH3:9])=[C:6]1[C:10](=[O:11])[O:12][#6:13]>>"
             "[#6:7][CH:1]=O.[CH3:8][C:3](=O)[CH2:2][C:14](=[O:15])[O:16][#6:17].[CH3:9][C:5](=O)[CH2:6][C:10](=[O:11])[O:12][#6:13].[NH3:4]",
             "EtOH or AcOH, reflux (one-pot, 4 components -> 3 here)", "Aldehyde + 2 beta-keto esters + ammonia builds the 1,4-dihydropyridine ring in one pot.", 0.8),
    Template("biltz", "Biltz hydantoin synthesis", "Ring",
             "[c:1][C:2]1([c:3])[NH:4][C:5](=[O:6])[NH:7][C:8]1=[O:9]>>[c:1][C:2](=O)[C:8](=[O:9])[c:3].[NH2:4][C:5](=[O:6])[NH2:7]",
             "NaOH or KOH, EtOH/H2O, reflux", "1,2-Diketone + urea; an aryl group migrates (benzilic-type rearrangement).", 0.6),
    Template("benzoin", "Benzoin condensation", "C-C",
             "[c:1][CH:2]([OH:3])[C:4](=[O:5])[c:6]>>[c:1][CH:2]=[O:3].[c:6][CH:4]=[O:5]",
             "NaCN or thiamine (cat.), EtOH/H2O", "Two aromatic aldehydes couple (umpolung); cyanide or thiamine makes the acyl anion equivalent.", 0.7),
    Template("carbonylation", "Pd-catalysed carbonylation (BHC ibuprofen step)", "C-C",
             "[c:1][CH:2]([CH3:3])[C:4](=[O:5])[OH:6]>>[c:1][CH:2]([CH3:3])[OH:6].[C-:4]#[O+:5]",
             "CO (high pressure), PdCl2(PPh3)2 (cat.), HCl", "Industrial step: benzylic alcohol + CO inserts the COOH group with 100% atom economy.", 0.6),
    # ------------------------------------------------------------ functional group interconversion (FGI)
    Template("nitro_reduction", "Reduction of a nitro group", "FGI",
             "[c:1][NH2:2]>>[c:1][N+:2](=O)[O-]",
             "H2, Pd/C  or  Sn/HCl  or  Fe/HCl", "An aromatic amine is made by reducing the nitro compound (nitration then reduction).", 0.85),
    Template("amine_debenzylation", "N-Benzyl protecting group (removed by hydrogenolysis)", "Protecting",
             "[NX3;H1,H2;!$(NC=[O,S]);!$(N-a):1]>>[N:1]Cc1ccccc1",
             "Deprotection: H2, Pd/C", "A benzyl group blocks over-alkylation of the amine and is removed in the last step.", 0.3),
    Template("amine_protect", "Protect an aromatic amine as an acetamide", "Protecting",
             "[c:1][NH2:2]>>[c:1][NH:2]C(C)=O",
             "Deprotection: aq. HCl or NaOH, reflux", "A free ArNH2 would be oxidised or acylated by harsh reagents (HNO3, ClSO3H); the acetamide tames it.", 0.4),
    Template("nitration", "Nitration of an arene", "FGI",
             "[c:1]-[N+:2](=[O:3])[O-:4]>>[c:1]",
             "HNO3, H2SO4, 0-50 C", "Electrophilic aromatic substitution; check the directing effects of groups already on the ring.", 0.8),
    Template("halogenation", "Halogenation of an arene", "FGI",
             "[c:1]-[Br,Cl:2]>>[c:1]",
             "Br2/FeBr3 or Cl2/FeCl3 (or Br2/AcOH for activated rings)", "Electrophilic aromatic substitution; activated rings (OH, NHAc) react without a catalyst.", 0.6),
    Template("chlorosulfonation", "Chlorosulfonation", "FGI",
             "[c:1]-[S:2](=[O:3])(=[O:4])Cl>>[c:1]",
             "ClSO3H (excess), 0 C to 60 C", "Electrophilic aromatic substitution that gives the sulfonyl chloride directly (para to an acetamido group).", 0.7),
    Template("reduce_ketone", "Reduction of a ketone", "FGI",
             "[CX4;H1:2]([OH:1])(-[#6:3])-[#6:4]>>[O:1]=[C:2](-[#6:3])-[#6:4]",
             "NaBH4, MeOH or EtOH (or H2, Raney Ni)", "A secondary alcohol is the reduction product of a ketone.", 0.65),
    Template("reduce_aldehyde", "Reduction of an aldehyde", "FGI",
             "[CH2:2]([OH:1])-[#6:3]>>[O:1]=[CH:2]-[#6:3]",
             "NaBH4, EtOH", "A primary alcohol is the reduction product of an aldehyde.", 0.55),
    Template("reduce_ester_arene", "Reduction of an aryl ester to a benzylic alcohol", "FGI",
             "[c:3][CH2:2][OH:1]>>[c:3][C:2](=[O:1])O[CH3]",
             "LiAlH4, dry THF; then H3O+", "LiAlH4 reduces esters (and ketones) to alcohols.", 0.45),
    Template("reduce_acid", "Reduction of a carboxylic acid", "FGI",
             "[CH2:2]([OH:1])-[#6:3]>>[O:1]=[C:2]([OH])-[#6:3]",
             "1) LiAlH4, dry THF; 2) H3O+", "LiAlH4 reduces acids and esters to primary alcohols (NaBH4 does not).", 0.4),
    Template("oxidise_alcohol_ketone", "Oxidation of a secondary alcohol", "FGI",
             "[#6:3][C:2](=[O:1])[#6:4]>>[#6:3][CH:2]([OH:1])[#6:4]",
             "PCC, CH2Cl2  or  Jones (CrO3, H2SO4, acetone)", "A ketone is the oxidation product of a secondary alcohol.", 0.55),
    Template("oxidise_alcohol_aldehyde", "Oxidation of a primary alcohol to an aldehyde", "FGI",
             "[#6:3][CH:2]=[O:1]>>[#6:3][CH2:2][OH:1]",
             "PCC, CH2Cl2 (anhydrous, stops at the aldehyde)", "Use PCC/Swern, not Jones or KMnO4, or you over-oxidise to the acid.", 0.45),
    Template("oxidise_side_chain", "Side-chain oxidation of an alkylbenzene", "FGI",
             "[c:3][C:2](=[O:1])[OH:4]>>[c:3][CH3:2]",
             "KMnO4, OH-, heat; then H3O+", "Any benzylic C-H carbon is oxidised all the way to COOH.", 0.5),
    Template("dehydration", "Dehydration of an alcohol (E1)", "FGI",
             "[CX3;!$(C[#6;a,$(C=O)]):1]=[CX3;!$(C[#6;a,$(C=O)]):2]>>[CX4:1]-[CX4:2]-[OH]",
             "H3PO4 or H2SO4, heat (distil the alkene out)", "Acid-catalysed loss of water; Zaitsev (more substituted) alkene is favoured.", 0.45),
    Template("halide_from_alcohol", "Alkyl halide from an alcohol", "FGI",
             "[CX4:1]-[Cl,Br,I:2]>>[CX4:1]-[OH]",
             "Tertiary: conc. HCl (SN1). Primary/secondary: SOCl2, PBr3 or HBr/H2SO4 (SN2)", "Convert the poor leaving group OH into a good one.", 0.6),
    Template("acid_chloride", "Acid chloride from a carboxylic acid", "FGI",
             "[#6:3][C:1](=[O:2])Cl>>[#6:3][C:1](=[O:2])[OH]",
             "SOCl2 or (COCl)2, DMF (cat.)", "Activate the acid before forming an amide or ester.", 0.5),
]
BY_ID = {t.id: t for t in TEMPLATES}
