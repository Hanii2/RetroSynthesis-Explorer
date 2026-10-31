"""Retrosynthesis Explorer for organic chemistry students.
Run: streamlit run app.py
"""
import random

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem

from retroexp import build3d, colors, compare, engine, library as L, quiz, safety, viewer
from retroexp.buildingblocks import name_of, is_available
from retroexp.chem import canon
from retroexp.resolve import ResolveError, resolve

st.set_page_config(
    page_title="Retrosynthesis Explorer",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ------------------------------------------------------------------ Pastel theme (readable text everywhere)
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@400;500;600;700;800&display=swap');

:root {
    --ink: #3b2a6b;          /* main dark text (high contrast on pastel) */
    --ink-soft: #6b5b95;     /* secondary text */
    --lav: #ede4ff;  --lav-2: #ddd0ff;
    --pink: #ffe1ee; --pink-2: #ffc9e0;
    --mint: #d9f7ec; --mint-2: #b8ecd8;
    --sky: #dcefff;  --sky-2: #bfe0ff;
    --peach: #ffe9d6;
    --line: #e4d8ff;
    --card: #ffffff;
}

/* ========== BACKGROUND ========== */
html, body, [data-testid="stAppViewContainer"], .stApp {
    background:
        radial-gradient(ellipse 70% 50% at 8% 0%, #ffe1ee 0%, transparent 55%),
        radial-gradient(ellipse 60% 45% at 95% 8%, #dcefff 0%, transparent 50%),
        radial-gradient(ellipse 55% 40% at 50% 100%, #d9f7ec 0%, transparent 50%),
        linear-gradient(180deg, #fdf7ff 0%, #f6f0ff 50%, #f2f8ff 100%) !important;
    color: var(--ink) !important;
    font-family: 'Poppins', system-ui, sans-serif;
}
[data-testid="stHeader"] { background: transparent !important; }
[data-testid="stToolbar"] { display: none; }
.stApp p, .stApp li, .stApp span, .stApp label, .stApp div[data-testid="stMarkdownContainer"] { color: var(--ink); }

/* ========== SIDEBAR ========== */
section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #fff1f8 0%, #f3ecff 45%, #e8f3ff 100%) !important;
    border-right: 1px solid var(--line);
    box-shadow: 2px 0 18px rgba(167, 139, 250, 0.10);
}
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] label p,
section[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p,
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {
    color: var(--ink) !important;
    font-weight: 600;
    font-size: 0.92rem;
}
section[data-testid="stSidebar"] h1, section[data-testid="stSidebar"] h2, section[data-testid="stSidebar"] h3 {
    background: linear-gradient(90deg, #8b6cf0, #e879a9);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
    font-weight: 700;
}

/* ========== TITLE ========== */
h1 {
    background: linear-gradient(90deg, #8b6cf0, #e879a9, #4fb6d9);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
    font-weight: 800 !important; letter-spacing: -0.6px; font-size: 2.35rem !important;
}
h2, h3 { color: #5b3fc4 !important; font-weight: 600 !important; }

/* ========== INPUTS: selectbox, text input (readable on every Streamlit version) ========== */
div[data-baseweb="select"] > div,
div[data-baseweb="select"] > div > div,
div[data-baseweb="input"],
div[data-baseweb="base-input"],
.stTextInput input, .stNumberInput input, .stTextArea textarea {
    background-color: #ffffff !important;
    color: var(--ink) !important;
    -webkit-text-fill-color: var(--ink) !important;
    border-color: #cdb8ff !important;
    border-radius: 12px !important;
}
div[data-baseweb="select"] * { color: var(--ink) !important; }
div[data-baseweb="select"] svg { fill: #8b6cf0 !important; color: #8b6cf0 !important; }
.stTextInput input::placeholder { color: #a99bcf !important; -webkit-text-fill-color: #a99bcf !important; }
div[data-baseweb="select"] > div:focus-within, .stTextInput input:focus {
    border-color: #a78bfa !important; box-shadow: 0 0 0 3px rgba(167, 139, 250, 0.25) !important;
}
/* dropdown list */
div[data-baseweb="popover"], div[data-baseweb="popover"] > div, ul[data-baseweb="menu"] {
    background: #ffffff !important; border-radius: 12px !important;
}
ul[data-baseweb="menu"] li, div[data-baseweb="popover"] li, li[role="option"] {
    background: #ffffff !important; color: var(--ink) !important;
}
ul[data-baseweb="menu"] li:hover, li[role="option"]:hover, li[aria-selected="true"] {
    background: var(--lav) !important;
}

/* ========== SLIDER ========== */
[data-testid="stSlider"] [data-testid="stTickBarMin"],
[data-testid="stSlider"] [data-testid="stTickBarMax"],
[data-testid="stSlider"] [data-testid="stThumbValue"] { color: var(--ink) !important; font-weight: 600; }
[data-testid="stSlider"] div[role="slider"] { background-color: #a78bfa !important; border-color: #fff !important; }

/* ========== TABS (pastel pills, text always visible) ========== */
div[data-testid="stTabs"] [data-baseweb="tab-list"] {
    gap: 0.45rem; background: rgba(255,255,255,0.75); border-radius: 16px;
    padding: 6px; border: 1px solid var(--line); flex-wrap: wrap;
}
button[data-baseweb="tab"] {
    background: var(--lav) !important; border-radius: 12px !important;
    padding: 8px 16px !important; height: auto !important; border: none !important;
    opacity: 1 !important;
}
button[data-baseweb="tab"] p, button[data-baseweb="tab"] div, button[data-baseweb="tab"] span {
    color: var(--ink) !important; font-weight: 500 !important; font-size: 0.95rem !important; opacity: 1 !important;
}
button[data-baseweb="tab"]:hover { background: var(--pink) !important; }
button[data-baseweb="tab"][aria-selected="true"] {
    background: linear-gradient(135deg, #cdbaff, #ffc3dc) !important;
    box-shadow: 0 3px 10px rgba(167, 139, 250, 0.35);
}
button[data-baseweb="tab"][aria-selected="true"] p,
button[data-baseweb="tab"][aria-selected="true"] div { color: #2d1b69 !important; font-weight: 700 !important; }
div[data-baseweb="tab-highlight"], div[data-baseweb="tab-border"] { display: none !important; }

/* ========== CARDS / EXPANDERS ========== */
div[data-testid="stExpander"] {
    background: var(--card) !important; border: 1px solid var(--line) !important;
    border-radius: 16px !important; margin-bottom: 1rem;
    box-shadow: 0 4px 20px rgba(167, 139, 250, 0.10);
}
div[data-testid="stExpander"] summary, div[data-testid="stExpander"] summary * { color: #4c2fb0 !important; font-weight: 600; }
div[data-testid="stExpander"]:hover { border-color: #c4b0ff !important; }

/* ========== METRICS / ALERTS ========== */
[data-testid="stMetric"] {
    background: linear-gradient(135deg, #fbf7ff, #f1f7ff); border: 1px solid var(--line);
    border-radius: 14px; padding: 0.8rem 1.1rem; box-shadow: 0 2px 12px rgba(167, 139, 250, 0.08);
}
[data-testid="stMetricValue"], [data-testid="stMetricLabel"] * { color: var(--ink) !important; }
.stAlert { border-radius: 12px !important; }
.stAlert p { color: var(--ink) !important; }

/* ========== BUTTONS (pastel gradient, dark text) ========== */
.stButton > button, .stDownloadButton > button {
    background: linear-gradient(135deg, #d9ccff, #ffd0e4) !important;
    border: 1px solid #cdb8ff !important; border-radius: 12px !important;
    font-weight: 600; padding: 0.5rem 1.4rem !important; transition: all 0.2s ease;
    box-shadow: 0 3px 12px rgba(167, 139, 250, 0.25);
}
.stButton > button p, .stDownloadButton > button p { color: #2d1b69 !important; font-weight: 600; }
.stButton > button:hover { transform: translateY(-2px); box-shadow: 0 8px 20px rgba(232, 121, 169, 0.30); }

/* ========== RADIO / CHECKBOX / CAPTIONS ========== */
[data-testid="stRadio"] label p, [data-testid="stCheckbox"] label p { color: var(--ink) !important; }
[data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] * { color: var(--ink-soft) !important; }

/* ========== DATAFRAMES ========== */
[data-testid="stDataFrame"] {
    border: 1px solid var(--line) !important; border-radius: 14px !important; overflow: hidden;
    box-shadow: 0 2px 12px rgba(167, 139, 250, 0.08); background: #ffffff !important;
}

/* ========== REFERENCE BADGES ========== */
.ref-badge { display: inline-block; padding: 6px 14px; border-radius: 999px; font-size: 0.85rem; font-weight: 600; margin: 6px 0; }
.ref-textbook { background: var(--mint); color: #0b6b4d !important; border: 1px solid var(--mint-2); }
.ref-engine   { background: var(--peach); color: #9a4a0a !important; border: 1px solid #ffd2a8; }
.pct-pill {
    display: inline-block; background: linear-gradient(135deg, #d9ccff, #ffd0e4); color: #2d1b69 !important;
    padding: 7px 16px; border-radius: 999px; font-weight: 700; font-size: 1.05rem;
    box-shadow: 0 4px 12px rgba(167, 139, 250, 0.28); border: 1px solid #cdb8ff;
}
.mol-card {
    background: var(--card); border: 1px solid var(--line); border-radius: 16px;
    padding: 1.1rem 1.4rem; margin-bottom: 1.2rem; box-shadow: 0 4px 18px rgba(167, 139, 250, 0.10);
}
.industrial-card {
    background: linear-gradient(135deg, #e9fbf3, #e6f4ff); border: 1px solid var(--mint-2);
    border-radius: 16px; padding: 1.2rem 1.5rem; margin-top: 1.5rem;
    box-shadow: 0 4px 16px rgba(94, 200, 160, 0.12);
}
.caption-soft { color: var(--ink-soft) !important; font-size: 0.95rem; line-height: 1.5; }
</style>
""", unsafe_allow_html=True)


# ------------------------------------------------------------------ helpers
def label(smiles):
    e = L.BY_CANON.get(canon(smiles))
    return name_of(smiles) or (e["name"] if e else None) or smiles


@st.cache_data(show_spinner=False)
def img(smiles, scheme, w=300, h=230, legend="", hl=(), hl2=()):
    return colors.draw2d(smiles, scheme, (w, h), legend=legend,
                         hl_bonds=list(hl) or None, hl2_bonds=list(hl2) or None)


@st.cache_data(show_spinner=False)
def c_resolve(q):
    return resolve(q)


@st.cache_data(show_spinner="Building routes ...")
def c_routes(smiles, depth, n):
    return engine.build_routes(smiles, depth=depth, n_routes=n)


@st.cache_data(show_spinner=False)
def c_block(smiles):
    return build3d.molblock(smiles)


@st.cache_data(show_spinner="Comparing with the textbook ...")
def c_compare(name, ri):
    return compare.compare(L.LIBRARY[name], ri)


def chips(findings):
    """Show only important chemistry warnings (fail / warn). Skip pure info."""
    for f in findings:
        if f.severity == "fail":
            st.error(f"**Problem:** {f.message}")
        elif f.severity == "warn":
            st.warning(f"**Note:** {f.message}")


def step_row(reactants, product, reaction, reagents, scheme, prefix=""):
    """One reaction step in the forward direction: reactants -> product."""
    cols = st.columns(len(reactants) * 2 + 1)
    for i, r in enumerate(reactants):
        cols[2 * i].image(img(r, scheme, 220, 170, legend=label(r) if label(r) != r else ""))
        if i < len(reactants) - 1:
            cols[2 * i + 1].markdown(
                "<div style='font-size:1.8em;text-align:center;padding-top:45px;color:#a855f7'>+</div>",
                unsafe_allow_html=True)
    cols[-1].image(img(product, scheme, 220, 170,
                       legend=label(product) if label(product) != product else ""))
    st.markdown(f"**{prefix}{reaction}**  \n*{reagents}*")


def forward_steps(node):
    """Post-order list of (reactants, product, proposal) from starting materials up to target."""
    out = []
    for c in node.children:
        out += forward_steps(c)
    if node.proposal:
        out.append((node.proposal.precursors, node.smiles, node.proposal))
    return out


def reliability_pct(node):
    """Simple percentage for display (based on the route score 0-1)."""
    return max(5, min(95, int(round(node.score * 100))))


# Industrial / commercial uses for common teaching molecules
INDUSTRIAL_USES = {
    "paracetamol": "Widely used as an analgesic and antipyretic (pain & fever relief). One of the most manufactured drugs worldwide. Found in tablets, syrups and combination medicines (e.g. with caffeine or codeine).",
    "acetaminophen": "Widely used as an analgesic and antipyretic (pain & fever relief). One of the most manufactured drugs worldwide. Found in tablets, syrups and combination medicines.",
    "aspirin": "Analgesic, antipyretic and anti-inflammatory drug. Also used in low doses as an antiplatelet agent to prevent heart attacks and strokes. Produced on a multi-thousand-tonne scale annually.",
    "ibuprofen": "Non-steroidal anti-inflammatory drug (NSAID) used for pain, inflammation and fever. Major commercial production via the BHC (Boots-Hoechst-Celanese) process. Available OTC worldwide.",
    "mtbe": "Fuel additive (oxygenate) that raises the octane number of gasoline and reduces emissions. Also used as a solvent in organic chemistry.",
    "biphenyl": "Used as a heat-transfer fluid, intermediate in the synthesis of dyes, polymers and pharmaceuticals, and as a fungistat for citrus fruits.",
    "phenol": "Large-scale industrial chemical used to make plastics (phenol-formaldehyde resins), bisphenol A, nylon precursors, detergents, and pharmaceuticals.",
    "benzene": "Fundamental petrochemical feedstock for styrene, phenol, cyclohexane, aniline and many other chemicals. Used in the production of plastics, resins, synthetic fibres and rubber.",
    "aniline": "Key intermediate for dyes, rubber chemicals, polyurethane precursors (MDI) and pharmaceuticals (e.g. paracetamol).",
    "salicylic acid": "Used in the manufacture of aspirin, and in cosmetics and dermatology (acne treatment, keratolytic agent).",
    "toluene": "Solvent and feedstock for benzene, xylenes, TNT and many fine chemicals. Used in paints, adhesives and as an octane booster.",
    "benzoic acid": "Food preservative (E210), precursor for plasticisers, benzoyl chloride and many pharmaceuticals.",
}


def get_industrial_uses(query: str, smiles: str = "") -> str:
    key = (query or "").strip().lower()
    if key in INDUSTRIAL_USES:
        return INDUSTRIAL_USES[key]
    # try partial match
    for name, text in INDUSTRIAL_USES.items():
        if name in key or key in name:
            return text
    # try library name
    entry = L.find(smiles) if smiles else None
    if entry:
        n = entry["name"].lower()
        if n in INDUSTRIAL_USES:
            return INDUSTRIAL_USES[n]
    return (
        "No specific industrial-use data is stored for this molecule in the teaching set. "
        "Many organic molecules are used as pharmaceutical intermediates, solvents, monomers or fine chemicals — "
        "check a reliable source (e.g. PubChem, Wikipedia, or a process chemistry reference) for commercial applications."
    )


def source_badge(is_textbook: bool, pct: int = None):
    """HTML badge + optional reliability percentage."""
    if is_textbook:
        html = (
            '<span class="ref-badge ref-textbook">📚 Textbook route '
            '(Clayden · Warren · standard undergrad examples)</span>'
        )
        if pct is not None:
            html += f' &nbsp; <span class="pct-pill">Confidence ~{pct}%</span>'
    else:
        html = (
            '<span class="ref-badge ref-engine">⚙️ Rule-based / engine-generated '
            '(no specific literature reference)</span>'
        )
        if pct is not None:
            html += f' &nbsp; <span class="pct-pill">Estimated reliability ~{pct}%</span>'
    return html


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.markdown("### Molecule selector")
    if "query" not in st.session_state:
        st.session_state.query = "Paracetamol"

    level = st.selectbox("Level", ["All"] + L.LEVELS)
    cats = ["All"] + [c for c in L.CATEGORIES
                      if any(e["category"] == c and (level == "All" or e["level"] == level)
                             for e in L.LIBRARY.values())]
    cat = st.selectbox("Category", cats)
    names = [n for n, e in L.LIBRARY.items()
             if (level == "All" or e["level"] == level) and (cat == "All" or e["category"] == cat)]

    def _picked():
        v = st.session_state.get("pick")
        if v and v != "—":
            st.session_state.query = v

    st.selectbox(f"Quick picks ({len(names)})", ["—"] + names, key="pick", on_change=_picked)
    st.text_input("Or type any name / SMILES", key="query")
    st.divider()
    scheme = st.selectbox("Colour scheme", colors.SCHEMES)
    depth = st.slider("Max steps in a route", 1, 6, 4)
    n_routes = st.slider("Number of routes to show", 1, 5, 3)
    st.markdown("---")
    st.caption("Teaching tool · rule-based · not an AI model")

# ------------------------------------------------------------------ target
st.title("Retrosynthesis Explorer")
st.markdown(
    '<p class="caption-soft">Find possible ways to make a molecule from simple starting materials. '
    'Built for organic chemistry students using classic textbook disconnection rules.</p>',
    unsafe_allow_html=True,
)

if safety.name_blocked(st.session_state.query):
    st.error(safety.MESSAGE.format(reason="controlled substance, explosive or chemical-weapon agent"))
    st.stop()
try:
    smiles, source = c_resolve(st.session_state.query)
except ResolveError as e:
    st.error(str(e))
    st.stop()
blocked, reason = safety.check(smiles, st.session_state.query)
if blocked:
    st.error(safety.MESSAGE.format(reason=reason))
    st.stop()
mol = Chem.MolFromSmiles(smiles)
if mol is None or mol.GetNumHeavyAtoms() > 60:
    st.error("Please choose a molecule with up to 60 heavy atoms.")
    st.stop()

entry = L.find(smiles)
st.markdown(
    f'<div class="mol-card"><strong>Target:</strong> {st.session_state.query} '
    f'&nbsp;·&nbsp; <code style="color:#7c3aed;background:#f5f3ff;padding:2px 8px;border-radius:6px">{smiles}</code></div>',
    unsafe_allow_html=True,
)

col1, col2 = st.columns([1, 1.4])
with col1:
    st.image(colors.draw2d(smiles, scheme, (380, 300)), use_container_width=True)
    st.markdown(colors.legend_html(colors.legend(mol, scheme)), unsafe_allow_html=True)
with col2:
    props = build3d.properties(smiles)
    rows = "".join(
        f'<tr style="background:{"#faf5ff" if i % 2 else "#ffffff"};">'
        f'<td style="padding:8px 14px;color:#5b21b6;font-weight:600;border-bottom:1px solid #ede9fe;">{k}</td>'
        f'<td style="padding:8px 14px;color:#1e1b4b;border-bottom:1px solid #ede9fe;">{v}</td></tr>'
        for i, (k, v) in enumerate(props.items())
    )
    st.markdown(
        f'<div style="background:#fff;border:1px solid #e9d5ff;border-radius:14px;overflow:hidden;'
        f'box-shadow:0 2px 12px rgba(139,92,246,0.06);">'
        f'<table style="width:100%;border-collapse:collapse;font-family:Poppins,sans-serif;font-size:0.92rem;">'
        f'<thead><tr style="background:#f5f3ff;">'
        f'<th style="padding:10px 14px;text-align:left;color:#5b21b6;font-weight:700;">Property</th>'
        f'<th style="padding:10px 14px;text-align:left;color:#5b21b6;font-weight:700;">Value</th>'
        f'</tr></thead><tbody>{rows}</tbody></table></div>',
        unsafe_allow_html=True,
    )
    if Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False):
        st.info("This molecule has stereocentres. The tool shows the constitution only (stereochemistry is ignored).")

tabs = st.tabs([
    "🔍 Possible last steps",
    "🛤️ Full routes",
    "📖 Textbook route",
    "📝 Practice quiz",
    "🌐 3D view",
    "ℹ️ About",
])

proposals = engine.propose(canon(smiles))

# ------------------------------------------------------------------ 1. one-step disconnections
with tabs[0]:
    st.markdown(
        "These are possible **last steps** that could give the target molecule.  \n"
        "The **red bond** is the one that is broken in the disconnection."
    )
    if not proposals:
        st.warning("No matching reaction rule found for this molecule.")
    else:
        top = proposals[:8]
        st.dataframe(
            pd.DataFrame([{
                "#": i + 1,
                "Reaction": p.template.name,
                "Precursors": " + ".join(label(x) for x in p.precursors),
                "Type": p.template.kind,
                "Score": f"{p.score:.2f}",
            } for i, p in enumerate(top)]),
            hide_index=True, use_container_width=True,
        )

        k = st.selectbox("Inspect this step", range(len(top)),
                         format_func=lambda i: f"{i + 1}. {top[i].template.name}")
        p = top[k]

        a, b, c = st.columns([2, 0.6, 2.4])
        a.image(img(p.target, scheme, 360, 260, "target", tuple(p.broken), tuple(p.changed)))
        b.markdown(
            "<div style='font-size:2.2em;text-align:center;padding-top:70px;color:#a855f7;font-weight:700'>⇒</div>",
            unsafe_allow_html=True,
        )
        with c:
            pcs = st.columns(max(1, len(p.precursors)))
            for col, x in zip(pcs, p.precursors):
                col.image(img(x, scheme, 200, 160, legend=label(x) if label(x) != x else ""))
                if is_available(x):
                    col.caption("✓ Available")
                else:
                    col.caption("Needs further steps")

        st.markdown(f"**Reaction:** {p.template.name}  \n**Typical conditions:** {p.template.reagents}")
        st.info(p.template.tip)
        chips(p.findings)

        # Reference note for one-step
        if entry:
            st.markdown(source_badge(True), unsafe_allow_html=True)
            st.caption(
                "This molecule is in the teaching library. The disconnections above are generated by the "
                "rule engine; the stored textbook route is shown in the ‘Textbook route’ tab."
            )
        else:
            st.markdown(source_badge(False, int(round(p.score * 100))), unsafe_allow_html=True)
            st.caption(
                "No literature reference is attached to this disconnection. "
                "It was produced by the rule-based engine (classic textbook reactions). "
                "The percentage is an internal ranking score, not a measured experimental yield."
            )

    # ---- Industrial / commercial uses (always show under this tab)
    st.markdown("---")
    uses = get_industrial_uses(st.session_state.query, smiles)
    st.markdown(
        f'<div class="industrial-card">'
        f'<h3 style="margin-top:0;color:#047857;">🏭 Industrial & Commercial Uses</h3>'
        f'<p style="margin-bottom:0;color:#134e4a;line-height:1.55;">{uses}</p>'
        f'</div>',
        unsafe_allow_html=True,
    )

# ------------------------------------------------------------------ 2. full routes
with tabs[1]:
    st.markdown(
        "Complete routes from simple starting materials to the target.  \n"
        "Steps are shown in the **forward** direction (starting materials → product)."
    )

    routes = c_routes(canon(smiles), depth, n_routes)
    if not routes:
        st.warning("No route found with the current settings.")
    else:
        for ri, node in enumerate(routes):
            sms = engine.starting_materials(node)
            complete = node.complete
            pct = reliability_pct(node)
            is_tb = bool(entry)

            title = f"Route {ri + 1}  ·  {node.n_steps} step(s)"
            if complete:
                title += "  ·  complete"
            else:
                title += "  ·  incomplete"

            with st.expander(title, expanded=(ri == 0)):
                st.markdown(source_badge(is_tb, pct), unsafe_allow_html=True)

                if is_tb:
                    st.success(
                        "Source: This route follows (or is close to) a textbook pathway stored in the library.  \n"
                        "Typical references: Clayden *Organic Chemistry*, "
                        "Warren *Organic Synthesis: The Disconnection Approach*, "
                        "or standard undergraduate synthesis examples. "
                        "These are **teaching routes**, not taken from a single patent or research paper."
                    )
                else:
                    st.warning(
                        "Source: **No authentic literature reference found** for this exact route.  \n"
                        "It was **generated by the rule-based engine** (classic textbook reactions only).  \n"
                        f"Estimated reliability: **~{pct}%** "
                        "(based on how common and reliable the reactions are, precursor availability, "
                        "and chemistry checks). This is **not** an experimental success probability."
                    )

                for si, (reac, prod, pr) in enumerate(forward_steps(node), 1):
                    step_row(reac, prod, pr.template.name, pr.template.reagents, scheme,
                             prefix=f"Step {si}: ")
                    chips([f for f in pr.findings if f.severity in ("fail", "warn")])
                    st.divider()

                st.markdown(
                    "**Starting materials:** " +
                    ", ".join(f"`{label(x)}`" for x in sorted(sms))
                )

# ------------------------------------------------------------------ 3. textbook route
with tabs[2]:
    if not entry:
        st.info(
            "This molecule is not in the built-in library of example compounds.  \n"
            "Choose one from the **Quick picks** list in the sidebar to see a stored textbook route "
            "with authentic teaching references."
        )
        st.markdown(source_badge(False), unsafe_allow_html=True)
        st.caption(
            "For molecules outside the library the tool only offers engine-generated routes "
            "(see the Full routes tab). No specific book or paper reference is claimed."
        )
    else:
        ri = 0
        if len(entry["routes"]) > 1:
            ri = st.radio(
                "Choose route", range(len(entry["routes"])),
                format_func=lambda i: entry["routes"][i]["name"], horizontal=True,
            )
        route = entry["routes"][ri]

        st.markdown(f"### {entry['name']} — {route['name']}")
        st.markdown(source_badge(True, 90), unsafe_allow_html=True)
        st.success(
            "Source: **Stored textbook route** from the project library.  \n"
            "These routes are based on standard undergraduate organic synthesis examples "
            "(Clayden *Organic Chemistry*, Warren’s *Disconnection Approach*, and common teaching routes). "
            "They are **not** taken from a single patent or research paper. "
            "For a formal report, verify and cite a solid source "
            "(e.g. Vardanyan & Hruby, *Synthesis of Essential Drugs*; Clayden; Warren)."
        )

        for si, s in enumerate(route["steps"], 1):
            step_row(s.reactants.split("."), s.product, s.reaction, s.reagents, scheme,
                     prefix=f"Step {si}: ")
            st.divider()

        # Simple comparison with engine
        cmp_ = c_compare(entry["name"], ri)
        st.subheader("Did the engine find a similar route?")
        fr = cmp_.get("first_step_rank")
        rr = cmp_.get("route_rank")
        c1, c2 = st.columns(2)
        c1.metric("Textbook last step ranked by engine", f"#{fr}" if fr else "Not found")
        c2.metric("Matching engine route", f"Route #{rr}" if rr else "Not in top 5")

# ------------------------------------------------------------------ 4. quiz
with tabs[3]:
    st.markdown("Test yourself with questions taken from the textbook routes.")
    qlevel = st.selectbox("Level", ["All"] + L.LEVELS, key="qlevel")

    if "quiz" not in st.session_state:
        st.session_state.quiz = {"seed": 1, "right": 0, "total": 0, "answered": False}
    qs = st.session_state.quiz

    q = quiz.make_question(random.Random(qs["seed"]), qlevel)
    s = q["step"]

    st.markdown(
        f"**Question:** {'Name this reaction' if q['kind'] == 'reaction' else 'Which reagents / conditions?'}  \n"
        f"*Target molecule of this route: {q['molecule']}*"
    )
    step_row(s.reactants.split("."), s.product, "?", "?", scheme)

    choice = st.radio(
        "Your answer", range(5), key=f"q{qs['seed']}",
        format_func=lambda i: "— choose —" if i == 0 else q["options"][i - 1],
    )

    col_a, col_b = st.columns(2)
    if col_a.button("Check answer") and choice > 0 and not qs["answered"]:
        qs["answered"] = True
        qs["total"] += 1
        correct = (choice - 1 == q["answer"])
        qs["right"] += int(correct)
        if correct:
            st.success(f"Correct!  {s.reaction} — {s.reagents}")
        else:
            st.error(f"Not quite.  Correct answer: {s.reaction} — {s.reagents}")

    if col_b.button("Next question"):
        qs["seed"] += 1
        qs["answered"] = False
        st.rerun()

    st.caption(f"Score: {qs['right']} / {qs['total']}")

# ------------------------------------------------------------------ 5. 3D
with tabs[4]:
    routes = c_routes(canon(smiles), depth, n_routes)  # ensure available
    mols = [smiles] + sorted({x for n in routes for x in engine.starting_materials(n)})
    pick = st.selectbox("Molecule to view", mols, format_func=lambda x: label(x))
    style = st.selectbox("Display style", list(viewer.STYLES))
    bg = st.selectbox("Background", list(viewer.BACKGROUNDS))

    try:
        html = viewer.viewer_html(c_block(pick), style, True, False, 640, 460, False, scheme, bg)
        components.html(html, height=480, scrolling=False)
    except Exception as e:
        st.warning(f"Could not generate 3D view: {e}")

# ------------------------------------------------------------------ about
with tabs[5]:
    st.markdown("""
### What this tool does
It suggests possible ways to make a molecule from simple starting materials using **classic textbook reactions**.

It uses about 50 well-known disconnection rules (Grignard, Friedel–Crafts, Williamson, amide formation, aldol, Wittig, Suzuki, Diels–Alder, etc.).

### References & reliability
- Routes for molecules **in the library** come from standard textbook examples
  (Clayden, Warren *Organic Synthesis: The Disconnection Approach*, common teaching routes).
  They are clearly marked with a green “Textbook route” badge.
- Routes for other molecules are **generated by the rule-based engine**.  
  They are marked with an amber badge and the message  
  **“No authentic literature reference found – rule-based / engine-generated”**.
- Every engine-generated route also shows an **estimated reliability percentage**
  (derived from reaction reliability, precursor availability, simplification and chemistry checks).
  This is a ranking score for teaching purposes, **not** an experimental yield or success probability.

### Important points
- This is a **teaching tool**, not a professional synthesis planner.
- It is **rule-based**, not a trained AI model.
- Stereochemistry is ignored.
- “Available” starting materials mean common laboratory chemicals or small molecules.
- Reaction conditions are typical textbook conditions; yields and selectivity are not predicted.

### Limits
- Only molecules up to ~60 heavy atoms.
- Does not suggest routes to controlled substances, explosives or chemical-weapon agents.
""")
    if st.button("Run self-check on the library"):
        r = compare.run_all()
        st.write(
            f"Checked **{r['n']}** textbook routes.  \n"
            f"Engine ranked the textbook last step as #1 in **{r['first_top1']}** cases, "
            f"and in the top 3 in **{r['first_top3']}** cases.  \n"
            f"An engine route reached the same starting materials in **{r['route_top5']}** of the top-5 results."
        )
