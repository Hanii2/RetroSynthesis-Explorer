"""
Library of teaching molecules with textbook / literature routes (FORWARD direction, starting materials -> product).

The retrosynthesis ENGINE never sees this library: it is used for (a) the quick-pick list, (b) the "Textbook route"
tab and (c) to measure how often the engine finds the textbook disconnection.

A step is  S(reactants, product, reaction name, reagents/conditions).  `reactants` lists the organic building blocks
that end up in the product ("A.B"); catalysts, acids, bases, oxidants and reductants are written under reagents.
"""
from dataclasses import dataclass

from .chem import canon


@dataclass(frozen=True)
class Step:
    reactants: str
    product: str
    reaction: str
    reagents: str


def S(reactants, product, reaction, reagents):
    return Step(reactants, product, reaction, reagents)


# ------------------------------------------------------------------ reusable steps
NITRATE_BENZENE = S("c1ccccc1", "O=[N+]([O-])c1ccccc1", "Nitration", "HNO3, H2SO4, 50 °C")
REDUCE_NITROBENZENE = S("O=[N+]([O-])c1ccccc1", "Nc1ccccc1", "Nitro reduction", "Sn/HCl, then NaOH  (or H2, Pd/C)")
ACETANILIDE = S("Nc1ccccc1.CC(=O)OC(C)=O", "CC(=O)Nc1ccccc1", "N-Acetylation", "Ac2O, NaOAc, H2O, rt")
ANILINE_CHAIN = [NITRATE_BENZENE, REDUCE_NITROBENZENE]
ACETANILIDE_CHAIN = ANILINE_CHAIN + [ACETANILIDE]
PARACETAMOL_CHAIN = [
    S("Oc1ccccc1", "O=[N+]([O-])c1ccc(O)cc1", "Nitration (para isomer isolated)", "dilute HNO3, 20 °C; separate o-/p- isomers by steam distillation"),
    S("O=[N+]([O-])c1ccc(O)cc1", "Nc1ccc(O)cc1", "Nitro reduction", "H2, Pd/C  (or Sn/HCl)"),
    S("Nc1ccc(O)cc1.CC(=O)OC(C)=O", "CC(=O)Nc1ccc(O)cc1", "N-Acetylation", "Ac2O, warm water (the NH2 reacts, not the OH)"),
]
KOLBE = S("Oc1ccccc1.O=C=O", "OC(=O)c1ccccc1O", "Kolbe-Schmitt carboxylation", "1) NaOH; 2) CO2, 125 °C, ~100 atm; 3) H3O+")
BENZOIC_GRIGNARD = S("Brc1ccccc1.O=C=O", "OC(=O)c1ccccc1", "Grignard carboxylation", "1) Mg, dry ether; 2) CO2 (dry ice); 3) H3O+")
FC_BENZOPHENONE = S("c1ccccc1.O=C(Cl)c1ccccc1", "O=C(c1ccccc1)c1ccccc1", "Friedel-Crafts acylation", "AlCl3, CH2Cl2, 0 °C to rt")
FC_ACETOPHENONE = S("c1ccccc1.CC(=O)Cl", "CC(=O)c1ccccc1", "Friedel-Crafts acylation", "AlCl3, benzene (excess), 60 °C")
METHYL_BENZOATE = S("OC(=O)c1ccccc1.CO", "COC(=O)c1ccccc1", "Fischer esterification", "H2SO4 (cat.), MeOH (excess), reflux")
BENZOIN = S("O=Cc1ccccc1.O=Cc1ccccc1", "O=C(c1ccccc1)C(O)c1ccccc1", "Benzoin condensation", "NaCN (or thiamine, NaOH), EtOH/H2O")
BENZIL = S("O=C(c1ccccc1)C(O)c1ccccc1", "O=C(c1ccccc1)C(=O)c1ccccc1", "Oxidation of benzoin", "HNO3, or Cu(OAc)2/NH4NO3, AcOH, heat")
CLAISEN_SCHMIDT_BA = S("O=Cc1ccccc1.CC(C)=O", "CC(=O)C=Cc1ccccc1", "Claisen-Schmidt condensation", "NaOH, H2O/EtOH, 25 °C")
IBU_FC = S("CC(C)Cc1ccccc1.CC(=O)OC(C)=O", "CC(C)Cc1ccc(C(C)=O)cc1", "Friedel-Crafts acylation", "HF (cat.) or AlCl3, 0-80 °C")
ACETOPHENONE_MANNICH_PRE = FC_ACETOPHENONE

# ------------------------------------------------------------------ entries
# (name, level, category, smiles, [(route name, [steps])], note)
ENTRIES = [
    # ---------------- esters and ethers
    ("Aspirin", "Basic", "Esters and ethers", "CC(=O)Oc1ccccc1C(=O)O",
     [("From phenol", [KOLBE, S("OC(=O)c1ccccc1O.CC(=O)OC(C)=O", "CC(=O)Oc1ccccc1C(=O)O", "O-Acetylation of a phenol", "H2SO4 or H3PO4 (cat.), 50-60 °C")])],
     "Phenols esterify poorly under Fischer conditions, so the anhydride is used."),
    ("Salicylic acid", "Basic", "Aromatic substitution", "OC(=O)c1ccccc1O", [("Kolbe-Schmitt", [KOLBE])],
     "CO2 is a weak electrophile; the phenoxide is the nucleophile."),
    ("Methyl salicylate (wintergreen)", "Basic", "Esters and ethers", "COC(=O)c1ccccc1O",
     [("Fischer", [KOLBE, S("OC(=O)c1ccccc1O.CO", "COC(=O)c1ccccc1O", "Fischer esterification", "H2SO4 (cat.), MeOH (excess), reflux")])], ""),
    ("Isopentyl acetate (banana oil)", "Basic", "Esters and ethers", "CC(C)CCOC(C)=O",
     [("Fischer", [S("CC(=O)O.CC(C)CCO", "CC(C)CCOC(C)=O", "Fischer esterification", "H2SO4 (cat.), reflux")])], "Classic teaching-lab ester."),
    ("Benzyl acetate", "Basic", "Esters and ethers", "CC(=O)OCc1ccccc1",
     [("Fischer", [S("CC(=O)O.OCc1ccccc1", "CC(=O)OCc1ccccc1", "Fischer esterification", "H2SO4 (cat.), reflux")])], ""),
    ("Methyl benzoate", "Basic", "Esters and ethers", "COC(=O)c1ccccc1", [("Grignard + Fischer", [BENZOIC_GRIGNARD, METHYL_BENZOATE])], ""),
    ("Ethyl benzoate", "Basic", "Esters and ethers", "CCOC(=O)c1ccccc1",
     [("Grignard + Fischer", [BENZOIC_GRIGNARD, S("OC(=O)c1ccccc1.CCO", "CCOC(=O)c1ccccc1", "Fischer esterification", "H2SO4 (cat.), EtOH (excess), reflux")])], ""),
    ("Anisole", "Basic", "Esters and ethers", "COc1ccccc1",
     [("Williamson", [S("Oc1ccccc1.CI", "COc1ccccc1", "Williamson ether synthesis", "NaOH (phenoxide), then CH3I")])], ""),
    ("Phenetole (ethyl phenyl ether)", "Basic", "Esters and ethers", "CCOc1ccccc1",
     [("Williamson", [S("Oc1ccccc1.CCBr", "CCOc1ccccc1", "Williamson ether synthesis", "K2CO3, acetone, reflux")])], ""),
    ("tert-Butyl methyl ether (MTBE)", "Basic", "Esters and ethers", "CC(C)(C)OC",
     [("Williamson (right way round)", [S("CC(C)(C)O.CI", "CC(C)(C)OC", "Williamson ether synthesis", "NaH or K (make t-BuO-), then CH3I")])],
     "Exam favourite: use the TERTIARY alkoxide + METHYL halide. t-BuBr + MeO- would eliminate."),
    ("Phenacetin", "Intermediate", "Esters and ethers", "CCOc1ccc(NC(C)=O)cc1",
     [("From paracetamol", PARACETAMOL_CHAIN + [S("CC(=O)Nc1ccc(O)cc1.CCI", "CCOc1ccc(NC(C)=O)cc1", "Williamson ether synthesis", "K2CO3, acetone, reflux")])], ""),
    # ---------------- amides, amines, local anaesthetics
    ("Acetanilide", "Basic", "Amides and amines", "CC(=O)Nc1ccccc1", [("From benzene", ACETANILIDE_CHAIN)], "Protects the NH2 for later ring substitution."),
    ("Paracetamol", "Basic", "Amides and amines", "CC(=O)Nc1ccc(O)cc1", [("From phenol", PARACETAMOL_CHAIN)],
     "Chemoselective: the amine is acetylated, the phenol OH is not."),
    ("Benzocaine", "Basic", "Amides and amines", "CCOC(=O)c1ccc(N)cc1",
     [("Fischer", [S("Nc1ccc(C(=O)O)cc1.CCO", "CCOC(=O)c1ccc(N)cc1", "Fischer esterification", "EtOH (excess), conc. H2SO4 (the amine is protonated), reflux; neutralise")])], ""),
    ("Procaine", "Intermediate", "Amides and amines", "CCN(CC)CCOC(=O)c1ccc(N)cc1",
     [("Nitro route", [S("O=C(O)c1ccc([N+](=O)[O-])cc1", "O=C(Cl)c1ccc([N+](=O)[O-])cc1", "Acid chloride formation", "SOCl2, reflux"),
                       S("O=C(Cl)c1ccc([N+](=O)[O-])cc1.CCN(CC)CCO", "CCN(CC)CCOC(=O)c1ccc([N+](=O)[O-])cc1", "Ester formation (acid chloride + alcohol)", "Et3N, CH2Cl2, 0 °C"),
                       S("CCN(CC)CCOC(=O)c1ccc([N+](=O)[O-])cc1", "CCN(CC)CCOC(=O)c1ccc(N)cc1", "Nitro reduction", "H2, Pd/C")])],
     "The NO2 is carried through and reduced last, because a free NH2 would also react with the acid chloride."),
    ("Lidocaine", "Intermediate", "Amides and amines", "CCN(CC)CC(=O)Nc1c(C)cccc1C",
     [("Chloroacetamide route", [S("Cc1cccc(C)c1N.ClCC(=O)Cl", "Cc1cccc(C)c1NC(=O)CCl", "Amide formation", "AcOH, NaOAc, rt"),
                                 S("Cc1cccc(C)c1NC(=O)CCl.CCNCC", "CCN(CC)CC(=O)Nc1c(C)cccc1C", "N-Alkylation (SN2)", "toluene, reflux")])],
     "Amide first, then SN2: the alpha-halo amide is a good electrophile."),
    ("DEET", "Basic", "Amides and amines", "CCN(CC)C(=O)c1cccc(C)c1",
     [("Acid chloride route", [S("Cc1cccc(C(=O)O)c1", "Cc1cccc(C(=O)Cl)c1", "Acid chloride formation", "SOCl2, reflux"),
                               S("Cc1cccc(C(=O)Cl)c1.CCNCC", "CCN(CC)C(=O)c1cccc(C)c1", "Amide formation", "Et2NH (excess) or Et3N, 0 °C")])], ""),
    ("Sulfanilamide", "Intermediate", "Amides and amines", "Nc1ccc(S(N)(=O)=O)cc1",
     [("Acetanilide route", ACETANILIDE_CHAIN + [
         S("CC(=O)Nc1ccccc1", "CC(=O)Nc1ccc(S(=O)(=O)Cl)cc1", "Chlorosulfonation", "ClSO3H (excess), 0 → 60 °C"),
         S("CC(=O)Nc1ccc(S(=O)(=O)Cl)cc1.N", "CC(=O)Nc1ccc(S(N)(=O)=O)cc1", "Sulfonamide formation", "aq. NH3"),
         S("CC(=O)Nc1ccc(S(N)(=O)=O)cc1", "Nc1ccc(S(N)(=O)=O)cc1", "Amide hydrolysis (deprotection)", "aq. HCl, reflux; neutralise")])],
     "The acetamide protects the amine from ClSO3H and directs para."),
    ("Propranolol", "Intermediate", "Drugs (multi-step)", "CC(C)NCC(O)COc1cccc2ccccc12",
     [("Epichlorohydrin route", [S("Oc1cccc2ccccc12.ClCC1CO1", "C1(COc2cccc3ccccc23)CO1", "Ether formation (epichlorohydrin)", "NaOH, H2O/EtOH"),
                                 S("C1(COc2cccc3ccccc23)CO1.CC(C)N", "CC(C)NCC(O)COc1cccc2ccccc12", "Epoxide ring opening", "iPrNH2 (excess), EtOH, reflux")])], ""),
    ("Fluoxetine", "Advanced", "Drugs (multi-step)", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
     [("Mannich route", [FC_ACETOPHENONE,
                         S("CC(=O)c1ccccc1.C=O.CN", "CNCCC(=O)c1ccccc1", "Mannich reaction", "MeNH2·HCl, (CH2O)n, EtOH, HCl (cat.), reflux"),
                         S("CNCCC(=O)c1ccccc1", "CNCCC(O)c1ccccc1", "Ketone reduction", "NaBH4, MeOH"),
                         S("CNCCC(O)c1ccccc1.Clc1ccc(C(F)(F)F)cc1", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "SNAr (aryl ether formation)", "NaH, DMSO, 80-100 °C")])],
     "SNAr works because CF3 is electron-withdrawing and para to the leaving group."),
    ("Salbutamol", "Advanced", "Drugs (multi-step)", "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",
     [("Glaxo-type route from methyl salicylate", [KOLBE, S("OC(=O)c1ccccc1O.CO", "COC(=O)c1ccccc1O", "Fischer esterification", "H2SO4 (cat.), MeOH, reflux"),
         S("COC(=O)c1ccccc1O.ClCC(=O)Cl", "COC(=O)c1cc(C(=O)CCl)ccc1O", "Friedel-Crafts acylation", "AlCl3, CS2 or nitrobenzene (para to OH)"),
         S("COC(=O)c1cc(C(=O)CCl)ccc1O.CC(C)(C)NCc1ccccc1", "COC(=O)c1cc(C(=O)CN(Cc2ccccc2)C(C)(C)C)ccc1O", "N-Alkylation (SN2)", "excess amine or K2CO3, MeCN"),
         S("COC(=O)c1cc(C(=O)CN(Cc2ccccc2)C(C)(C)C)ccc1O", "CC(C)(C)N(Cc1ccccc1)CC(O)c1ccc(O)c(CO)c1", "Reduction of ketone and ester", "LiAlH4, THF; then H3O+"),
         S("CC(C)(C)N(Cc1ccccc1)CC(O)c1ccc(O)c(CO)c1", "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "Debenzylation (hydrogenolysis)", "H2, Pd/C")])],
     "The N-benzyl group blocks over-alkylation and is removed last; LiAlH4 reduces the ketone AND the ester."),
    ("Ibuprofen", "Intermediate", "Drugs (multi-step)", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
     [("BHC (green, 3 steps)", [IBU_FC, S("CC(C)Cc1ccc(C(C)=O)cc1", "CC(C)Cc1ccc(C(C)O)cc1", "Ketone reduction", "H2, Raney Ni"),
                                S("CC(C)Cc1ccc(C(C)O)cc1.[C-]#[O+]", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "Pd-catalysed carbonylation", "CO, PdCl2(PPh3)2, HCl, 130 °C")]),
      ("Boots (original, 6 steps)", [IBU_FC,
         S("CC(C)Cc1ccc(C(C)=O)cc1.ClCC(=O)OCC", "CCOC(=O)C1OC1(C)c1ccc(CC(C)C)cc1", "Darzens condensation", "NaOEt, EtOH"),
         S("CCOC(=O)C1OC1(C)c1ccc(CC(C)C)cc1", "CC(C=O)c1ccc(CC(C)C)cc1", "Hydrolysis and decarboxylation", "1) NaOEt/H2O; 2) H3O+, heat (-CO2)"),
         S("CC(C=O)c1ccc(CC(C)C)cc1.NO", "CC(C=NO)c1ccc(CC(C)C)cc1", "Oxime formation", "NH2OH·HCl, base"),
         S("CC(C=NO)c1ccc(CC(C)C)cc1", "CC(C#N)c1ccc(CC(C)C)cc1", "Dehydration of oxime to nitrile", "Ac2O, heat"),
         S("CC(C#N)c1ccc(CC(C)C)cc1", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "Nitrile hydrolysis", "H3O+, heat")])],
     "Same target, two strategies: compare number of steps and atom economy (a classic green-chemistry case)."),
    ("Phenytoin", "Advanced", "Drugs (multi-step)", "O=C1NC(c2ccccc2)(c2ccccc2)C(=O)N1",
     [("Benzoin -> benzil -> Biltz", [BENZOIN, BENZIL, S("O=C(c1ccccc1)C(=O)c1ccccc1.NC(N)=O", "O=C1NC(c2ccccc2)(c2ccccc2)C(=O)N1", "Biltz hydantoin synthesis", "KOH, EtOH/H2O, reflux; H3O+")])],
     "A phenyl group migrates (benzilic-acid-type rearrangement)."),
    ("Nifedipine", "Advanced", "Drugs (multi-step)", "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
     [("Hantzsch", [S("O=Cc1ccccc1[N+](=O)[O-].COC(=O)CC(C)=O.COC(=O)CC(C)=O.N", "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]", "Hantzsch dihydropyridine synthesis", "MeOH, reflux (protect from light)")])],
     "One pot: aldehyde + 2 equivalents of methyl acetoacetate + ammonia."),
    ("Warfarin", "Advanced", "Drugs (multi-step)", "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O",
     [("Michael addition", [CLAISEN_SCHMIDT_BA, S("O=C1C=C(O)c2ccccc2O1.CC(=O)C=Cc1ccccc1", "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O", "Michael addition", "pyridine or piperidine (cat.), H2O or MeOH, reflux")])],
     "Made as a racemate: the new stereocentre has no chiral control."),
    # ---------------- aromatic substitution
    ("Nitrobenzene", "Basic", "Aromatic substitution", "O=[N+]([O-])c1ccccc1", [("Nitration", [NITRATE_BENZENE])], ""),
    ("Aniline", "Basic", "Aromatic substitution", "Nc1ccccc1", [("Nitration then reduction", ANILINE_CHAIN)], ""),
    ("4-Nitroacetanilide", "Intermediate", "Aromatic substitution", "CC(=O)Nc1ccc([N+](=O)[O-])cc1",
     [("Nitrate acetanilide", ACETANILIDE_CHAIN + [S("CC(=O)Nc1ccccc1", "CC(=O)Nc1ccc([N+](=O)[O-])cc1", "Nitration (para)", "HNO3, H2SO4, 0-10 °C")])],
     "NHAc is an ortho/para director; direct nitration of aniline would oxidise it."),
    ("4-Bromoacetanilide", "Intermediate", "Aromatic substitution", "CC(=O)Nc1ccc(Br)cc1",
     [("Brominate acetanilide", ACETANILIDE_CHAIN + [S("CC(=O)Nc1ccccc1", "CC(=O)Nc1ccc(Br)cc1", "Bromination (para)", "Br2, AcOH")])], ""),
    ("Methyl 3-nitrobenzoate", "Intermediate", "Aromatic substitution", "COC(=O)c1cccc([N+](=O)[O-])c1",
     [("Nitrate methyl benzoate", [BENZOIC_GRIGNARD, METHYL_BENZOATE, S("COC(=O)c1ccccc1", "COC(=O)c1cccc([N+](=O)[O-])c1", "Nitration (meta)", "HNO3, H2SO4, below 15 °C")])],
     "The ester is meta-directing."),
    ("Acetophenone", "Basic", "Aromatic substitution", "CC(=O)c1ccccc1", [("Friedel-Crafts", [FC_ACETOPHENONE])], ""),
    ("Benzophenone", "Basic", "Aromatic substitution", "O=C(c1ccccc1)c1ccccc1", [("Friedel-Crafts", [FC_BENZOPHENONE])], ""),
    ("Biphenyl", "Intermediate", "C-C coupling", "c1ccc(-c2ccccc2)cc1",
     [("Suzuki", [S("OB(O)c1ccccc1.Brc1ccccc1", "c1ccc(-c2ccccc2)cc1", "Suzuki coupling", "Pd(PPh3)4, Na2CO3, toluene/EtOH/H2O, 80 °C")])], ""),
    ("Bisphenol A", "Intermediate", "Aromatic substitution", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
     [("Acid-catalysed condensation", [S("Oc1ccccc1.Oc1ccccc1.CC(C)=O", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", "Acid-catalysed condensation (para)", "HCl (cat.), phenol in excess, 50-70 °C")])],
     "The engine has no rule for this one: a good example of its limits."),
    # ---------------- Grignard and alcohols
    ("Benzoic acid", "Basic", "Grignard and alcohols", "OC(=O)c1ccccc1", [("Grignard + CO2", [BENZOIC_GRIGNARD])], ""),
    ("2-Phenyl-2-propanol", "Basic", "Grignard and alcohols", "CC(C)(O)c1ccccc1",
     [("PhMgBr + acetone", [S("Brc1ccccc1.CC(C)=O", "CC(C)(O)c1ccccc1", "Grignard addition", "1) Mg, dry ether; 2) acetone; 3) H3O+")])],
     "Three C-C disconnections are possible for a tertiary alcohol; the one giving the simplest ketone is neat."),
    ("Triphenylmethanol", "Intermediate", "Grignard and alcohols", "OC(c1ccccc1)(c1ccccc1)c1ccccc1",
     [("PhMgBr + benzophenone", [FC_BENZOPHENONE, S("Brc1ccccc1.O=C(c1ccccc1)c1ccccc1", "OC(c1ccccc1)(c1ccccc1)c1ccccc1", "Grignard addition", "1) Mg, dry ether; 2) benzophenone; 3) H3O+")])],
     "Alternative: methyl benzoate + 2 PhMgBr (symmetry: two identical groups)."),
    ("2-Methyl-2-butanol", "Basic", "Grignard and alcohols", "CCC(C)(C)O",
     [("EtMgBr + acetone", [S("CCBr.CC(C)=O", "CCC(C)(C)O", "Grignard addition", "1) Mg, dry ether; 2) acetone; 3) H3O+")])], ""),
    ("3-Methyl-3-pentanol", "Basic", "Grignard and alcohols", "CCC(C)(O)CC",
     [("EtMgBr + 2-butanone", [S("CCBr.CCC(C)=O", "CCC(C)(O)CC", "Grignard addition", "1) Mg, dry ether; 2) 2-butanone; 3) H3O+")])], ""),
    ("1-Phenylethanol", "Basic", "Grignard and alcohols", "CC(O)c1ccccc1",
     [("Reduce acetophenone", [FC_ACETOPHENONE, S("CC(=O)c1ccccc1", "CC(O)c1ccccc1", "Ketone reduction", "NaBH4, EtOH")])], ""),
    ("Diphenylmethanol (benzhydrol)", "Basic", "Grignard and alcohols", "OC(c1ccccc1)c1ccccc1",
     [("Reduce benzophenone", [FC_BENZOPHENONE, S("O=C(c1ccccc1)c1ccccc1", "OC(c1ccccc1)c1ccccc1", "Ketone reduction", "NaBH4, EtOH")])], ""),
    # ---------------- C=C formation and condensations
    ("trans-Chalcone", "Basic", "Condensations and alkenes", "O=C(C=Cc1ccccc1)c1ccccc1",
     [("Claisen-Schmidt", [S("O=Cc1ccccc1.CC(=O)c1ccccc1", "O=C(C=Cc1ccccc1)c1ccccc1", "Claisen-Schmidt condensation", "NaOH, EtOH/H2O, rt")])],
     "Benzaldehyde has no alpha-H, so it cannot self-condense."),
    ("Dibenzalacetone", "Intermediate", "Condensations and alkenes", "O=C(C=Cc1ccccc1)C=Cc1ccccc1",
     [("Double Claisen-Schmidt", [S("O=Cc1ccccc1.O=Cc1ccccc1.CC(C)=O", "O=C(C=Cc1ccccc1)C=Cc1ccccc1", "Double Claisen-Schmidt condensation", "NaOH, EtOH/H2O (2 PhCHO : 1 acetone)")])], ""),
    ("Cinnamic acid", "Basic", "Condensations and alkenes", "OC(=O)C=Cc1ccccc1",
     [("Knoevenagel-Doebner", [S("O=Cc1ccccc1.OC(=O)CC(=O)O", "OC(=O)C=Cc1ccccc1", "Knoevenagel-Doebner condensation", "pyridine, piperidine (cat.), reflux (-CO2)")])], ""),
    ("trans-Stilbene", "Intermediate", "Condensations and alkenes", "C(=Cc1ccccc1)c1ccccc1",
     [("Wittig", [S("O=Cc1ccccc1.BrCc1ccccc1", "C(=Cc1ccccc1)c1ccccc1", "Wittig olefination", "1) PPh3; 2) NaOH (ylide); 3) PhCHO")])], ""),
    ("Ethyl acetoacetate", "Intermediate", "Condensations and alkenes", "CCOC(=O)CC(C)=O",
     [("Claisen", [S("CCOC(C)=O.CCOC(C)=O", "CCOC(=O)CC(C)=O", "Claisen condensation", "NaOEt, EtOH; then H3O+")])], "Use the SAME alkoxide as the ester (OEt with ethyl esters)."),
    ("Benzoin", "Basic", "Condensations and alkenes", "O=C(c1ccccc1)C(O)c1ccccc1", [("Benzoin condensation", [BENZOIN])], ""),
    ("Benzil", "Intermediate", "Condensations and alkenes", "O=C(c1ccccc1)C(=O)c1ccccc1", [("Benzoin then oxidation", [BENZOIN, BENZIL])], ""),
    ("Coumarin", "Intermediate", "Condensations and alkenes", "O=C1C=Cc2ccccc2O1",
     [("Perkin", [S("O=Cc1ccccc1O.CC(=O)OC(C)=O", "O=C1C=Cc2ccccc2O1", "Perkin condensation", "NaOAc, 180 °C")])], ""),
    # ---------------- rings
    ("cis-4-Cyclohexene-1,2-dicarboxylic anhydride", "Basic", "Rings (Diels-Alder etc.)", "O=C1OC(=O)C2CC=CCC12",
     [("Diels-Alder", [S("C=CC=C.O=C1C=CC(=O)O1", "O=C1OC(=O)C2CC=CCC12", "Diels-Alder reaction", "heat (xylene, 110 °C); butadiene from 3-sulfolene")])], ""),
    ("Norbornene dicarboxylic anhydride", "Basic", "Rings (Diels-Alder etc.)", "O=C1OC(=O)C2C3C=CC(C3)C12",
     [("Diels-Alder (endo)", [S("C1=CC=CC1.O=C1C=CC(=O)O1", "O=C1OC(=O)C2C3C=CC(C3)C12", "Diels-Alder reaction (endo)", "rt; cyclopentadiene from cracking dicyclopentadiene")])], ""),
    ("Wieland-Miescher ketone", "Advanced", "Rings (Diels-Alder etc.)", "CC12CCC(=O)C=C1CCCC2=O",
     [("Robinson annulation", [S("CC1C(=O)CCCC1=O.C=CC(C)=O", "CC(=O)CCC1(C)C(=O)CCCC1=O", "Michael addition", "Et3N or KOH (cat.), MeOH"),
                               S("CC(=O)CCC1(C)C(=O)CCCC1=O", "CC12CCC(=O)C=C1CCCC2=O", "Intramolecular aldol condensation", "pyrrolidine or (S)-proline; then acid, heat")])],
     "Robinson annulation = Michael addition + intramolecular aldol condensation."),
    # ---------------- functional group interconversion
    ("Cyclohexene", "Basic", "Functional group changes", "C1=CCCCC1", [("Dehydration", [S("OC1CCCCC1", "C1=CCCCC1", "Dehydration of an alcohol (E1)", "85% H3PO4, heat; distil the alkene out")])], ""),
    ("Cyclohexanone", "Basic", "Functional group changes", "O=C1CCCCC1", [("Oxidation", [S("OC1CCCCC1", "O=C1CCCCC1", "Oxidation of a secondary alcohol", "Na2Cr2O7/H2SO4 (Jones) or PCC")])], ""),
    ("tert-Butyl chloride", "Basic", "Functional group changes", "CC(C)(C)Cl", [("SN1", [S("CC(C)(C)O", "CC(C)(C)Cl", "SN1 substitution with HCl", "conc. HCl, rt")])], ""),
    ("1-Bromobutane", "Basic", "Functional group changes", "CCCCBr", [("SN2", [S("CCCCO", "CCCCBr", "SN2 substitution", "NaBr, H2SO4, reflux")])], ""),
]


def _make():
    out = {}
    for name, level, cat, smi, routes, note in ENTRIES:
        out[name] = {"name": name, "level": level, "category": cat, "smiles": smi, "note": note,
                     "routes": [{"name": rn, "steps": list(steps)} for rn, steps in routes]}
    return out


LIBRARY = _make()
CATEGORIES = list(dict.fromkeys(e["category"] for e in LIBRARY.values()))
LEVELS = ["Basic", "Intermediate", "Advanced"]
BY_CANON = {canon(e["smiles"]): e for e in LIBRARY.values()}


def find(smiles_or_name):
    """Library entry for a name (any case) or a SMILES, else None."""
    for n, e in LIBRARY.items():
        if n.lower() == (smiles_or_name or "").strip().lower():
            return e
    c = canon(smiles_or_name)
    return BY_CANON.get(c) if c else None


def route_tree(entry, route_index=0):
    """Textbook route as a nested dict {smiles, step, children}, assembled from the step list."""
    steps = entry["routes"][route_index]["steps"]
    by_product = {canon(s.product): s for s in steps}

    def build(smi, seen):
        c = canon(smi)
        st = by_product.get(c)
        if st is None or c in seen:
            return {"smiles": c, "step": None, "children": []}
        kids = [build(r, seen | {c}) for r in dict.fromkeys(st.reactants.split("."))]
        return {"smiles": c, "step": st, "children": kids}

    return build(entry["smiles"], frozenset())
