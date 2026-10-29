"""Retrosynthesis Explorer for organic chemistry students.
Run: streamlit run app.py
"""
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
                "<div style='font-size:1.8em;text-align:center;padding-top:45px'>+</div>",
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


# ------------------------------------------------------------------ sidebar
st.sidebar.header("Select molecule")
if "query" not in st.session_state:
    st.session_state.query = "Paracetamol"

level = st.sidebar.selectbox("Level", ["All"] + L.LEVELS)
cats = ["All"] + [c for c in L.CATEGORIES
                  if any(e["category"] == c and (level == "All" or e["level"] == level)
                         for e in L.LIBRARY.values())]
cat = st.sidebar.selectbox("Category", cats)
names = [n for n, e in L.LIBRARY.items()
         if (level == "All" or e["level"] == level) and (cat == "All" or e["category"] == cat)]


def _picked():
    v = st.session_state.get("pick")
    if v and v != "—":
        st.session_state.query = v


st.sidebar.selectbox(f"Quick picks ({len(names)})", ["—"] + names, key="pick", on_change=_picked)
st.sidebar.text_input("Or type any name / SMILES", key="query")
st.sidebar.divider()
scheme = st.sidebar.selectbox("Colour scheme", colors.SCHEMES)
depth = st.sidebar.slider("Max steps in a route", 1, 6, 4)
n_routes = st.sidebar.slider("Number of routes to show", 1, 5, 3)

# ------------------------------------------------------------------ target
st.title("🧪 Retrosynthesis Explorer")
st.caption("Find possible ways to Plan a backward Synthesis from the Target Molecule.")

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
st.markdown(f"**Target:** {st.session_state.query}  &nbsp;·&nbsp; `{smiles}`")

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

tabs = st.tabs(["Possible last steps", "Full routes", "Textbook route", "Practice quiz", "3D view", "About"])

proposals = engine.propose(canon(smiles))

# ------------------------------------------------------------------ 1. one-step disconnections
with tabs[0]:
    st.markdown("These are possible **last steps** that could give the target molecule.  \n"
                "The **red bond** is the one that is broken in the disconnection.")
    if not proposals:
        st.warning("No matching reaction rule found for this molecule.")
    else:
        top = proposals[:8]
        # Simple table – no technical score breakdown
        st.dataframe(
            pd.DataFrame([{
                "#": i + 1,
                "Reaction": p.template.name,
                "Precursors": " + ".join(label(x) for x in p.precursors),
                "Type": p.template.kind,
            } for i, p in enumerate(top)]),
            hide_index=True, use_container_width=True
        )

        k = st.selectbox("Inspect this step", range(len(top)),
                         format_func=lambda i: f"{i + 1}. {top[i].template.name}")
        p = top[k]

        a, b, c = st.columns([2, 0.6, 2.4])
        a.image(img(p.target, scheme, 360, 260, "target", tuple(p.broken), tuple(p.changed)))
        b.markdown("<div style='font-size:2.2em;text-align:center;padding-top:70px'>⇒</div>",
                   unsafe_allow_html=True)
        with c:
            pcs = st.columns(max(1, len(p.precursors)))
            for col, x in zip(pcs, p.precursors):
                col.image(img(x, scheme, 200, 160, legend=label(x) if label(x) != x else ""))
                if engine.is_available(x):
                    col.caption("✓ Available")
                else:
                    col.caption("Needs further steps")

        st.markdown(f"**Reaction:** {p.template.name}  \n**Typical conditions:** {p.template.reagents}")
        st.info(p.template.tip)
        chips(p.findings)

# ------------------------------------------------------------------ 2. full routes
with tabs[1]:
    st.markdown("Complete routes from simple starting materials to the target.  \n"
                "Steps are shown in the **forward** direction (starting materials → product).")

    routes = c_routes(canon(smiles), depth, n_routes)
    if not routes:
        st.warning("No route found with the current settings.")
    else:
        for ri, node in enumerate(routes):
            sms = engine.starting_materials(node)
            complete = node.complete
            pct = reliability_pct(node)

            title = f"Route {ri + 1}  ·  {node.n_steps} step(s)"
            if complete:
                title += f"  ·  complete"
            else:
                title += f"  ·  incomplete"

            with st.expander(title, expanded=(ri == 0)):
                # Reference / source note
                if entry:
                    st.success(
                        "Source: This route follows (or is close to) a textbook pathway stored in the library.  \n"
                        "Typical references for such routes: Clayden *Organic Chemistry*, "
                        "Warren *Organic Synthesis: The Disconnection Approach*, "
                        "or standard undergraduate synthesis examples."
                    )
                else:
                    st.warning(
                        f"Source: This route was **generated by the rule-based engine** "
                        f"(not taken from a specific book or paper).  \n"
                        f"Estimated reliability: **~{pct}%** "
                        f"(based on how common and reliable the reactions are)."
                    )

                for si, (reac, prod, pr) in enumerate(forward_steps(node), 1):
                    step_row(reac, prod, pr.template.name, pr.template.reagents, scheme,
                             prefix=f"Step {si}: ")
                    chips([f for f in pr.findings if f.severity in ("fail", "warn")])
                    st.divider()

                st.markdown("**Starting materials:** " +
                            ", ".join(f"`{label(x)}`" for x in sorted(sms)))

# ------------------------------------------------------------------ 3. textbook route
with tabs[2]:
    if not entry:
        st.info("This molecule is not in the built-in library of 56 example compounds.  \n"
                "Choose one from the **Quick picks** list in the sidebar to see a stored textbook route.")
    else:
        ri = 0
        if len(entry["routes"]) > 1:
            ri = st.radio("Choose route", range(len(entry["routes"])),
                          format_func=lambda i: entry["routes"][i]["name"], horizontal=True)
        route = entry["routes"][ri]

        st.markdown(f"### {entry['name']} — {route['name']}")
        st.success(
            "Source: Stored textbook route from the project library.  \n"
            "These routes are based on standard undergraduate organic synthesis examples "
            "(Clayden, Warren’s Disconnection Approach, and common teaching routes). "
            "They are **not** taken from a single patent or research paper."
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

    st.markdown(f"**Question:** {'Name this reaction' if q['kind'] == 'reaction' else 'Which reagents / conditions?'}  \n"
                f"*Target molecule of this route: {q['molecule']}*")
    step_row(s.reactants.split("."), s.product, "?", "?", scheme)

    choice = st.radio("Your answer", range(5), key=f"q{qs['seed']}",
                      format_func=lambda i: "— choose —" if i == 0 else q["options"][i - 1])

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
        # Always use components.html to avoid any Streamlit object leakage
        components.html(html, height=480, scrolling=False)
    except Exception as e:
        st.warning(f"Could not generate 3D view: {e}")

# ------------------------------------------------------------------ about
with tabs[5]:
    st.markdown("""
### What this tool does
It suggests possible ways to make a molecule from simple starting materials using **classic textbook reactions**.

It uses about 50 well-known disconnection rules (Grignard, Friedel–Crafts, Williamson, amide formation, aldol, Wittig, Suzuki, Diels–Alder, etc.).

### Important points
- This is a **teaching tool**, not a professional synthesis planner.
- It is **rule-based**, not a trained AI model.
- Stereochemistry is ignored.
- “Available” starting materials mean common laboratory chemicals or small molecules.
- Routes shown for molecules **in the library** come from standard textbook examples.
- Routes for other molecules are **generated by the engine**. In those cases the tool clearly states that the route is engine-generated and gives an estimated reliability percentage.

### Limits
- Only molecules up to ~60 heavy atoms.
- Does not suggest routes to controlled substances, explosives or chemical-weapon agents.
- Reaction conditions are typical textbook conditions; yields and selectivity are not predicted.
""")
    if st.button("Run self-check on the library"):
        r = compare.run_all()
        st.write(
            f"Checked **{r['n']}** textbook routes.  \n"
            f"Engine ranked the textbook last step as #1 in **{r['first_top1']}** cases, "
            f"and in the top 3 in **{r['first_top3']}** cases.  \n"
            f"An engine route reached the same starting materials in **{r['route_top5']}** of the top-5 results."
        )
