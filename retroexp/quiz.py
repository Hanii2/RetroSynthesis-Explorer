"""Practice questions made from the library's textbook steps (no engine involved, so every answer is unambiguous)."""
import random

from . import library as L


def all_steps(level=None):
    out = []
    for e in L.LIBRARY.values():
        if level and level != "All" and e["level"] != level:
            continue
        for r in e["routes"]:
            for s in r["steps"]:
                out.append((e["name"], s))
    return out


def make_question(rng, level="All", kind=None):
    steps = all_steps(level) or all_steps()
    pool = all_steps()
    name, step = rng.choice(steps)
    kind = kind or rng.choice(["reaction", "reagents"])
    key = "reaction" if kind == "reaction" else "reagents"
    correct = getattr(step, key)
    others = sorted({getattr(s, key) for _, s in pool if getattr(s, key) != correct})
    options = rng.sample(others, 3) + [correct]
    rng.shuffle(options)
    return {"kind": kind, "molecule": name, "step": step, "options": options, "answer": options.index(correct)}
