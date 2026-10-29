"""Retrosynthesis Explorer for organic chemistry students.   Run:  streamlit run app.py"""
import random

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem

from retroexp import build3d, colors, compare, engine, library as L, quiz, safety, viewer
from retroexp.buildingblocks import name_of
from retroexp.chem import canon
from retroexp.resolve import ResolveError, resolve

st.set_page_config(page_title="Retrosynthesis Explorer", page_icon="🧪", layout="wide")


# ------------------------------------------------------------------ helpers
def label(smiles):
    e = L.BY_CANON.get(canon(smiles))
    return name_of(smiles) or (e["name"] if e else None) or smiles


@st.cache_data(show_spinner=False)
def img(smiles, scheme, w=300, h=230, legend="", hl=(), hl2=()):
    return colors.draw2d(smiles, scheme, (w, h), legend=legend, hl_bonds=list(hl) or None, hl2_bonds=list(hl2) or None)


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
    for f in findings:
        {"fail": st.error, "warn": st.warning, "info": st.info}[f.severity](f"**{f.rule.replace('_', ' ')}**: {f.message}")


def step_row(reactants, product, reaction, reagents, scheme, prefix=""):
    """One reaction step in the forward direction: reactants -> product."""
    cols = st.columns(len(reactants) * 2 + 1)
    for i, r in enumerate(reactants):
        cols[2 * i].image(img(r, scheme, 230, 180, legend=label(r) if label(r) != r else ""))
        if i < len(reactants) - 1:
            cols[2 * i + 1].markdown("<div style='font-size:2em;text-align:center;padding-top:50px'>+</div>", unsafe_allow_html=True)
    cols[-1].image(img(product, scheme, 230, 180, legend=label(product) if label(product) != product else ""))
    st.markdown(f"**{prefix}{reaction}** — {reagents}")


def forward_steps(node):
    """Post-order list of (reactants, product, proposal) from the starting materials up to the target."""
    out = []
    for c in node.children:
        out += forward_steps(c)
    if node.proposal:
        out.append((node.proposal.precursors, node.smiles, node.proposal))
    return out


# ------------------------------------------------------------------ sidebar
st.sidebar.header("Choose a target molecule")
if "query" not in st.session_state:
    st.session_state.query = "Paracetamol"
level = st.sidebar.selectbox("Level", ["All"] + L.LEVELS)
cats = ["All"] + [c for c in L.CATEGORIES if any(e["category"] == c and (level == "All" or e["level"] == level) for e in L.LIBRARY.values())]
cat = st.sidebar.selectbox("Category", cats)
names = [n for n, e in L.LIBRARY.items() if (level == "All" or e["level"] == level) and (cat == "All" or e["category"] == cat)]


def _picked():
    v = st.session_state.get("pick")
    if v and v != "—":
        st.session_state.query = v


st.sidebar.selectbox(f"Quick picks ({len(names)} molecules)", ["—"] + names, key="pick", on_change=_picked)
st.sidebar.text_input("...or type ANY compound name / SMILES", key="query")
st.sidebar.divider()
scheme = st.sidebar.selectbox("Colour scheme", colors.SCHEMES,
                              help="Functional groups: OH, C=O, amine, ring ... get their own colours (legend under the picture).")
depth = st.sidebar.slider("Max route depth (steps)", 1, 6, 5)
n_routes = st.sidebar.slider("Routes to show", 1, 6, 3)

# ------------------------------------------------------------------ target
st.title("🧪 Retrosynthesis Explorer")
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
    st.error("Please choose a molecule with up to 60 heavy atoms (this is a teaching tool for ordinary organic targets).")
    st.stop()
entry = L.find(smiles)
st.caption(f"**{st.session_state.query}** → `{smiles}`  (source: {source})")
head_l, head_r = st.columns([1, 2])
head_l.image(colors.draw2d(smiles, scheme, (420, 320)), caption=f"Target ({scheme})")
head_l.markdown(colors.legend_html(colors.legend(mol, scheme)), unsafe_allow_html=True)
with head_r:
    st.dataframe(pd.DataFrame(build3d.properties(smiles).items(), columns=["Property", "Value"]).astype(str), hide_index=True)
    if entry:
        st.success(f"In the library: **{entry['level']}** · {entry['category']}. {entry['note']}")
    if Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False):
        st.info("This molecule has stereocentres. The engine ignores stereochemistry: it shows the constitution only.")

tabs = st.tabs(["1. Disconnections", "2. Full routes", "3. Textbook route", "4. Practice quiz", "5. 3D", "About"])
proposals = engine.propose(canon(smiles))

# ------------------------------------------------------------------ 1. one-step disconnections
with tabs[0]:
    st.write("Each row is ONE possible last step, written as a **disconnection** (retro direction). The red bond is the bond that is broken; "
             "orange = bond changed (functional-group change).")
    if not proposals:
        st.warning("No disconnection rule applies to this molecule. (The engine knows ~48 textbook reactions; see About.)")
    else:
        top = proposals[:10]
        st.dataframe(pd.DataFrame([{"#": i + 1, "Reaction": p.template.name, "Type": p.template.kind,
                                    "Precursors": " + ".join(label(x) for x in p.precursors), "Score": round(p.score, 2),
                                    "Checks": ", ".join(sorted({f.severity for f in p.findings if f.severity != 'info'})) or "ok"}
                                   for i, p in enumerate(top)]), hide_index=True)
        k = st.selectbox("Look at disconnection", range(len(top)), format_func=lambda i: f"{i + 1}. {top[i].template.name}")
        p = top[k]
        a, b, c = st.columns([2, 1, 3])
        a.image(img(p.target, scheme, 420, 300, "target", tuple(p.broken), tuple(p.changed)))
        b.markdown("<div style='font-size:2.5em;text-align:center;padding-top:80px'>⇒</div>", unsafe_allow_html=True)
        with c:
            pc = st.columns(len(p.precursors))
            for col, x in zip(pc, p.precursors):
                col.image(img(x, scheme, 240, 200, legend=label(x) if label(x) != x else ""))
                col.caption("✅ available" if engine.is_available(x) else "needs more steps")
        st.markdown(f"**Forward reaction:** {p.template.name}  \n**Reagents / conditions:** {p.template.reagents}")
        st.info(f"💡 {p.template.tip}")
        chips(p.findings)
        st.write("**Why this rank?** (score = 40% reaction reliability + 25% precursors available + 15% simplification + 20% chemistry checks)")
        st.bar_chart(pd.DataFrame({"part": list(p.parts), "value": list(p.parts.values())}).set_index("part"))

# ------------------------------------------------------------------ 2. full routes
with tabs[1]:
    st.write("Routes are built by applying the disconnections again and again until every starting material is available. "
             "Steps below are in the **forward** direction (starting materials → target).")
    routes = c_routes(canon(smiles), depth, n_routes)
    if not routes:
        st.warning("No route found.")
    for ri, node in enumerate(routes):
        sms = engine.starting_materials(node)
        with st.expander(f"Route {ri + 1}: {node.n_steps} step(s) · score {node.score:.2f} · "
                         f"{'complete' if node.complete else 'INCOMPLETE (some starting materials need more steps)'}", expanded=(ri == 0)):
            for si, (reac, prod, pr) in enumerate(forward_steps(node), 1):
                step_row(reac, prod, pr.template.name, pr.template.reagents, scheme, prefix=f"Step {si}: ")
                chips([f for f in pr.findings if f.severity != "info"])
                st.divider()
            st.write("**Starting materials:** " + ", ".join(f"{label(x)} {'✅' if engine.is_available(x) else '⚠️'}" for x in sorted(sms)))

# ------------------------------------------------------------------ 3. textbook route
with tabs[2]:
    if not entry:
        st.info("This molecule is not in the built-in library, so there is no stored textbook route. Pick one from **Quick picks** "
                "(56 molecules) to compare the engine with the textbook.")
    else:
        ri = 0
        if len(entry["routes"]) > 1:
            ri = st.radio("Route", range(len(entry["routes"])), format_func=lambda i: entry["routes"][i]["name"], horizontal=True)
        route = entry["routes"][ri]
        st.write(f"**{entry['name']}: {route['name']}**  ({len(route['steps'])} steps, forward direction)")
        for si, s in enumerate(route["steps"], 1):
            step_row(s.reactants.split("."), s.product, s.reaction, s.reagents, scheme, prefix=f"Step {si}: ")
            st.divider()
        cmp_ = c_compare(entry["name"], ri)
        st.subheader("Did the engine find this on its own?")
        st.caption("The engine never sees the library; this comparison is made afterwards.")
        c1, c2 = st.columns(2)
        fr = cmp_["first_step_rank"]
        c1.metric("Textbook LAST step among the engine's disconnections", f"#{fr}" if fr else "not found")
        rr = cmp_["route_rank"]
        c2.metric("Engine route with the same starting materials", f"route #{rr}" if rr else "not in top 5")

# ------------------------------------------------------------------ 4. quiz
with tabs[3]:
    st.write("Practice questions from the textbook steps in the library.")
    qlevel = st.selectbox("Level", ["All"] + L.LEVELS, key="qlevel")
    if "quiz" not in st.session_state:
        st.session_state.quiz = {"seed": 1, "right": 0, "total": 0, "answered": False}
    qs = st.session_state.quiz
    q = quiz.make_question(random.Random(qs["seed"]), qlevel)
    s = q["step"]
    st.write(f"**{'Name this reaction' if q['kind'] == 'reaction' else 'Which reagents / conditions?'}** (target of this route: {q['molecule']})")
    step_row(s.reactants.split("."), s.product, "?", "?", scheme)
    choice = st.radio("Your answer", range(5), key=f"q{qs['seed']}",
                      format_func=lambda i: "— choose an answer —" if i == 0 else q["options"][i - 1])
    if st.button("Check answer") and choice > 0 and not qs["answered"]:
        qs["answered"], qs["total"] = True, qs["total"] + 1
        qs["right"] += int(choice - 1 == q["answer"])
        (st.success if choice - 1 == q["answer"] else st.error)(
            ("Correct! " if choice - 1 == q["answer"] else "Not quite. ") + f"{s.reaction} — {s.reagents}")
    if st.button("Next question"):
        qs["seed"] += 1
        qs["answered"] = False
        st.rerun()
    st.caption(f"Score: {qs['right']} / {qs['total']}")

# ------------------------------------------------------------------ 5. 3D
with tabs[4]:
    best = routes[0] if routes else None
    mols = [smiles] + sorted({x for n in routes for x in engine.starting_materials(n)})
    pick = st.selectbox("Molecule", mols, format_func=lambda x: label(x))
    style = st.selectbox("Style", list(viewer.STYLES))
    bg = st.selectbox("Background", list(viewer.BACKGROUNDS))
    try:
        html = viewer.viewer_html(c_block(pick), style, True, False, 640, 460, False, scheme, bg)
        (st.iframe(html, height=470) if hasattr(st, "iframe") else components.html(html, height=470))
    except build3d.InputError as e:
        st.warning(str(e))

# ------------------------------------------------------------------ about
with tabs[5]:
    st.markdown("""
**What this is:** a *rule-based* retrosynthesis teaching tool. It applies ~48 textbook disconnections (Fischer, Williamson, amide formation,
Grignard, Friedel-Crafts, aldol, Wittig, Suzuki, Diels-Alder, Hantzsch ... plus functional-group changes) with RDKit reaction SMARTS, ranks them,
checks them with simple chemistry rules (chemoselectivity, directing effects, acidic protons ...), and repeats until it reaches available starting materials.
It is **not** a trained AI model and does not know every reaction.

**Limits:** stereochemistry is ignored; reagents are generic; "available" means a common lab chemical or a small molecule; scores are ranking numbers,
not probabilities; only molecules up to 60 heavy atoms; it will not give routes to controlled drugs, explosives or chemical-weapon agents.
""")
    if st.button("Run self-check against the 56-molecule library"):
        r = compare.run_all()
        st.write(f"Routes tested: **{r['n']}**.  Textbook last step is the engine's #1 disconnection in **{r['first_top1']}**, in the top 3 in **{r['first_top3']}**.  "
                 f"An engine route reaches the textbook's starting materials in **{r['route_top5']}** (top 5).")
