# Retrosynthesis Explorer (BS Organic Chemistry project)

Kisi bhi molecule ka **naam ya SMILES** likho (jaise `Paracetamol`): app batati hai ke use kin building blocks se banaya ja sakta hai.
Disconnection approach (tumhare Organic Synthesis course wali) ke rules par chalti hai, **rangeen 2D** structures dikhati hai,
aur **56 molecules** pehle se added hain (textbook routes ke saath). Project poora **Python** mein hai, AiZynthFinder jaisa bhari model nahi.

## Install (ek baar)
1. `setup.bat` par double-click (Windows) ya `bash setup.sh` (Linux/WSL). Internet chahiye, chhoti download hai.
2. Phir `run_app.bat` par double-click. Browser mein app khul jati hai.
(Manual: `pip install -r requirements.txt` phir `streamlit run app.py`.)

## Kya kya hai
| Tab | Kya karta hai |
|---|---|
| **1. Disconnections** | Ek step peeche jane ke saare mumkin tareeqe, rank ke saath. Jo bond todna hai wo **laal** dikhta hai, precursors ki tasveer, reagents, "kyun" ka tip, aur chemistry checks (fail / warn / info) |
| **2. Full routes** | Poori route starting materials tak (forward direction mein step cards). Galat steps wali routes hata di jati hain |
| **3. Textbook route** | Library ke molecules ka textbook route + **engine ne ye khud dhoondha ya nahi** ka comparison |
| **4. Practice quiz** | Reaction ka naam ya reagents pehchano (textbook steps se), score ke saath |
| **5. 3D** | Target aur starting materials ko 3D mein ghumao |
| **About** | Limits, aur ek button jo engine ko poori library par check karta hai |

Sidebar: Level (Basic / Intermediate / Advanced), Category, Quick picks, Colour scheme (Vivid / **Functional groups** / Rainbow / Classic), route depth.

## Engine kaise kaam karta hai (viva ke liye)
- ~50 textbook disconnection rules (reaction SMARTS, retro direction): Fischer, Williamson, amide/sulfonamide, Grignard (3 qisam), Friedel-Crafts, aldol/Claisen-Schmidt,
  Claisen, Mannich, Michael, Wittig, Suzuki, Diels-Alder, Hantzsch, Biltz, Kolbe-Schmitt, epoxide opening, SNAr, nitration/halogenation, reductions/oxidations ...
- Har proposal ko score milta hai: **40% reaction reliability + 25% precursors available + 15% simplification + 20% chemistry checks**.
- Chemistry checks: chemoselectivity (reduction, N vs O acylation), Friedel-Crafts (deactivated ring, rearrangement), **directing effects** (ortho/meta/para),
  Grignard ke saath acidic protons / ester / nitrile, amino-acid chloride, protecting group ...
- Route builder tab tak dohrata hai jab tak sab starting materials "available" (common lab chemical ya chhota molecule) na ho jayen.
- Ye **AI model nahi**, rule-based teaching tool hai. Score probability nahi, ranking number hai. Stereochemistry ignore hoti hai.

## Test hua / nahi hua (honest)
- **Test hua** (`python tests/test_core.py`): library ke 56 molecules / 110 steps consistent; saare 50 templates compile aur 44 library targets par fire karte hain (baqi 6 apni misaal par alag se check kiye);
  textbook cases (MTBE mein tertiary alkoxide + methyl halide, biphenyl = Suzuki ...); chemistry checks; app ke saare tabs headless render; quiz aur self-check buttons.
- **Engine vs textbook (engine library nahi dekhta, comparison baad mein):** 57 routes mein textbook ka aakhri step engine ka **#1 37 mein**, **top-3 50 mein**;
  top-5 routes mein textbook jaise starting materials **47/57**. Misses jahan rule hi nahi (Robinson annulation, Perkin, Darzens/Boots, bisphenol A) — ye honest limits hain.
- **Test NAHI hua:** browser mein 3D viewer ki asli shakal (internet se 3Dmol.js load hota hai), PubChem se naam lookup (mere paas internet nahi tha; built-in list offline chalti hai), aur `setup.bat`/`run_app.bat`.

## Scope / safety
Sirf aam teaching molecules (60 heavy atoms tak). Controlled drugs, explosives aur chemical-weapon agents ke liye tool route nahi deta (best-effort filter, poori guarantee nahi).

## Report mein limitations likho
1. Rule-based engine hai, jo sirf shamil ~50 reactions jaanta hai; novel chemistry nahi.
2. Reagents generic hain, conditions textbook-level hain, yields/selectivity predict nahi hoti.
3. "Available" ka matlab common lab chemical ya chhota molecule; asli supplier check nahi.
4. Stereochemistry ignore hoti hai (racemic maan kar constitution dikhata hai).
5. Library ke routes literature/textbook ke hain; solid sources se verify karke citation do (Vardanyan & Hruby, *Synthesis of Essential Drugs*; Clayden; Warren *Organic Synthesis: The Disconnection Approach*).

## Files
`app.py` · `retroexp/templates.py` (disconnection rules) · `engine.py` (propose + route builder) · `rules.py` (chemistry checks) · `library.py` (56 molecules + routes) ·
`buildingblocks.py` (available chemicals + names) · `safety.py` · `compare.py` · `quiz.py` · `colors.py` / `viewer.py` / `build3d.py` / `resolve.py` · `tests/`
