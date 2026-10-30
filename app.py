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
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ------------------------------------------------------------------ Modern colorful CSS (attractive multi-color theme)
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap');

/* ========== MAIN BACKGROUND ========== */
html, body, [data-testid="stAppViewContainer"] {
    background: 
        radial-gradient(ellipse 80% 50% at 20% -10%, rgba(99, 102, 241, 0.25), transparent),
        radial-gradient(ellipse 60% 40% at 80% 10%, rgba(236, 72, 153, 0.18), transparent),
        radial-gradient(ellipse 50% 30% at 50% 100%, rgba(20, 184, 166, 0.15), transparent),
        linear-gradient(165deg, #0f0c29 0%, #1a1a2e 40%, #16213e 70%, #0f3460 100%) !important;
    color: #f1f5f9 !important;
    font-family: 'Poppins', system-ui, sans-serif;
}
[data-testid="stHeader"] { background: transparent !important; }
[data-testid="stToolbar"] { display: none; }

/* subtle animated glow overlay */
[data-testid="stAppViewContainer"]::before {
    content: "";
    position: fixed;
    top: 0; left: 0; right: 0; bottom: 0;
    background: 
        radial-gradient(circle at 15% 20%, rgba(129, 140, 248, 0.08) 0%, transparent 40%),
        radial-gradient(circle at 85% 30%, rgba(244, 114, 182, 0.07) 0%, transparent 40%),
        radial-gradient(circle at 50% 80%, rgba(45, 212, 191, 0.06) 0%, transparent 40%);
    pointer-events: none;
    z-index: 0;
}

/* ========== SIDEBAR ========== */
section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #1e1b4b 0%, #312e81 50%, #1e3a5f 100%) !important;
    border-right: 1px solid rgba(165, 180, 252, 0.2);
    box-shadow: 4px 0 24px rgba(0,0,0,0.3);
}
section[data-testid="stSidebar"] * { color: #e0e7ff !important; }
section[data-testid="stSidebar"] .stSelectbox label,
section[data-testid="stSidebar"] .stSlider label,
section[data-testid="stSidebar"] .stTextInput label {
    color: #c7d2fe !important;
    font-weight: 600;
    font-size: 0.9rem;
}
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] h3 {
    background: linear-gradient(90deg, #a5b4fc, #f9a8d4);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
}

/* ========== TITLE ========== */
h1 {
    background: linear-gradient(90deg, #818cf8, #e879f9, #38bdf8, #2dd4bf);
    background-size: 300% 100%;
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    font-weight: 800 !important;
    letter-spacing: -0.8px;
    font-size: 2.4rem !important;
}
h2, h3 { 
    color: #c4b5fd !important; 
    font-weight: 600 !important;
}

/* ========== CARDS / EXPANDERS ========== */
div[data-testid="stExpander"] {
    background: linear-gradient(135deg, rgba(30, 27, 75, 0.7), rgba(49, 46, 129, 0.5)) !important;
    border: 1px solid rgba(165, 180, 252, 0.25) !important;
    border-radius: 16px !important;
    margin-bottom: 1rem;
    box-shadow: 0 8px 32px rgba(0,0,0,0.25);
    backdrop-filter: blur(8px);
}
div[data-testid="stExpander"] summary {
    color: #e0e7ff !important;
    font-weight: 600;
}
div[data-testid="stExpander"]:hover {
    border-color: rgba(192, 132, 252, 0.4) !important;
}

/* ========== METRICS ========== */
[data-testid="stMetric"] {
    background: linear-gradient(135deg, rgba(67, 56, 202, 0.4), rgba(126, 34, 206, 0.3));
    border: 1px solid rgba(167, 139, 250, 0.3);
    border-radius: 14px;
    padding: 0.8rem 1.1rem;
    box-shadow: 0 4px 16px rgba(0,0,0,0.2);
}
.stAlert { border-radius: 12px !important; }

/* ========== TABS ========== */
button[data-baseweb="tab"] {
    color: #94a3b8 !important;
    font-weight: 500;
    font-size: 0.95rem !important;
}
button[data-baseweb="tab"][aria-selected="true"] {
    color: #f0abfc !important;
    border-bottom-color: #e879f9 !important;
    font-weight: 600 !important;
}
div[data-testid="stTabs"] [data-baseweb="tab-list"] {
    gap: 0.5rem;
    background: rgba(30, 27, 75, 0.4);
    border-radius: 12px;
    padding: 4px;
    border: 1px solid rgba(165, 180, 252, 0.15);
}

/* ========== BUTTONS ========== */
.stButton > button {
    background: linear-gradient(135deg, #7c3aed, #db2777, #0891b2) !important;
    background-size: 200% 100% !important;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    font-weight: 600;
    padding: 0.5rem 1.4rem !important;
    transition: all 0.25s ease;
    box-shadow: 0 4px 15px rgba(124, 58, 237, 0.35);
}
.stButton > button:hover {
    background-position: 100% 0 !important;
    transform: translateY(-2px);
    box-shadow: 0 8px 25px rgba(219, 39, 119, 0.4);
}

/* ========== DATAFRAMES ========== */
[data-testid="stDataFrame"] {
    border: 1px solid rgba(165, 180, 252, 0.2);
    border-radius: 14px;
    overflow: hidden;
    box-shadow: 0 4px 20px rgba(0,0,0,0.2);
}

/* ========== REFERENCE BADGES ========== */
.ref-badge {
    display: inline-block;
    padding: 6px 14px;
    border-radius: 999px;
    font-size: 0.85rem;
    font-weight: 600;
    margin: 6px 0;
    letter-spacing: 0.2px;
}
.ref-textbook {
    background: linear-gradient(135deg, rgba(16, 185, 129, 0.2), rgba(20, 184, 166, 0.15));
    color: #6ee7b7;
    border: 1px solid #34d399;
    box-shadow: 0 0 12px rgba(52, 211, 153, 0.2);
}
.ref-engine {
    background: linear-gradient(135deg, rgba(245, 158, 11, 0.18), rgba(249, 115, 22, 0.12));
    color: #fcd34d;
    border: 1px solid #f59e0b;
    box-shadow: 0 0 12px rgba(245, 158, 11, 0.2);
}
.pct-pill {
    display: inline-block;
    background: linear-gradient(135deg, #7c3aed, #db2777);
    color: #fff;
    padding: 7px 16px;
    border-radius: 999px;
    font-weight: 700;
    font-size: 1.05rem;
    letter-spacing: 0.3px;
    box-shadow: 0 4px 14px rgba(124, 58, 237, 0.4);
    border: none;
}
.mol-card {
    background: linear-gradient(135deg, rgba(49, 46, 129, 0.5), rgba(30, 58, 138, 0.4));
    border: 1px solid rgba(165, 180, 252, 0.25);
    border-radius: 16px;
    padding: 1.1rem 1.4rem;
    margin-bottom: 1.2rem;
    box-shadow: 0 8px 32px rgba(0,0,0,0.25);
    backdrop-filter: blur(6px);
}
.caption-soft { 
    color: #a5b4fc !important; 
    font-size: 0.95rem; 
    line-height: 1.5;
}

/* nicer inputs */
.stTextInput input, .stSelectbox [data-baseweb="select"] {
    border-radius: 10px !important;
    border: 1px solid rgba(165, 180, 252, 0.3) !important;
    background: rgba(30, 27, 75, 0.5) !important;
}
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
                "<div style='font-size:1.8em;text-align:center;padding-top:45px;color:#e879f9'>+</div>",
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
    st.markdown("### 🧪 Molecule selector")
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
st.title("🧪 Retrosynthesis Explorer")
st.markdown(
    '<p class="caption-soft">Find possible ways to make a molecule from simple starting materials. '
    'Built for chemistry students using classic textbook disconnection rules.</p>',
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
    f'&nbsp;·&nbsp; <code style="color:#a7f3d0">{smiles}</code></div>',
    unsafe_allow_html=True,
)

col1, col2 = st.columns([1, 1.4])
with col1:
    st.image(colors.draw2d(smiles, scheme, (380, 300)), use_container_width=True)
    st.markdown(colors.legend_html(colors.legend(mol, scheme)), unsafe_allow_html=True)
with col2:
    props = build3d.properties(smiles)
    st.dataframe(pd.DataFrame(props.items(), columns=["Property", "Value"]).astype(str),
                 hide_index=True, use_container_width=True)
    if Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False):
        st.info("This molecule has stereocentres. The tool shows the constitution only (stereochemistry is ignored).")

tabs = st.tabs([
    "🔍 Disconnections",
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
            "<div style='font-size:2.2em;text-align:center;padding-top:70px;background:linear-gradient(90deg,#818cf8,#e879f9);-webkit-background-clip:text;-webkit-text-fill-color:transparent;font-weight:700'>⇒</div>",
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
