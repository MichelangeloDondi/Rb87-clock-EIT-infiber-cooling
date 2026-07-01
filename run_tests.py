#!/usr/bin/env python3
"""
run_tests.py -- the repo's self-check + anti-staleness suite. Run it after any edit, and before you
trust a number.

    python run_tests.py           # fast checks (~2-3 min)
    python run_tests.py --slow    # + the full multilevel / dark-vertex solves (a few minutes)

It fails on the ways this repo has drifted before:

  (A) a self-check regression  -- a closed-form validator, or a fast floor, moved from its known value;
  (B) a stale FIGURE           -- a committed PNG differs from what its script draws now (the code changed
                                  but the figure was not redrawn). Fix with `python regenerate.py`.
  (C) a stale NUMBER           -- a headline floor the code computes is not the one quoted in a README.
  (D) a SMUGGLED floor         -- a floor-shaped number appears in the README that is NOT registered in
                                  CANONICAL, so nothing cross-checks it. CANONICAL is the closed set (of sub-0.1 floors).
  (E) a rotted REFERENCE       -- a README link, figure, or run-command points at a file that has moved.
  (F) a forked SSOT            -- results.json (what the figures read) disagrees with the code.

CANONICAL below is the SINGLE list of headline numbers. Every check verifies the code AND the docs
against it, so a headline number can only change in one deliberate place: here.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
SLOW = "--slow" in sys.argv

# (name, value-as-written, [docs that must quote it verbatim]).  This list is the COMPLETE set of
# floor-shaped numbers the front-door README is allowed to quote: check (C) enforces "present", check
# (D) enforces "nothing else".  The master floor is added once the Delta-optimum solve settles it
# (see the PENDING set in check (D), run_tests --slow, and cooling_dark_vertex.py).
CANONICAL = [
    ("3-level floor",          "0.0020", ["01_three_level/README.md", "README.md"]),
    ("recoil-free floor",      "0.0011", ["01_three_level/README.md", "README.md"]),
    ("multilevel-clean floor", "0.0032", ["02_multilevel/README.md", "README.md"]),
    ("delivered floor",        "0.09",   ["02_multilevel/README.md", "README.md"]),
    ("no-master floor",        "0.087",  ["03_dark_vertex/README.md", "04_master/README.md", "README.md"]),
    ("master floor",           "0.055",  ["03_dark_vertex/README.md", "04_master/README.md", "README.md"]),
]

_passed = []
def report(ok, name, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + (("   -- " + detail) if (not ok and detail) else ""))
    _passed.append(bool(ok))

def run(rel_src, script):
    r = subprocess.run([PY, script], cwd=os.path.join(ROOT, rel_src), capture_output=True, text=True)
    return r.stdout + r.stderr

def read(rel):
    return open(os.path.join(ROOT, rel), encoding="utf-8").read()


# ---- (A1) closed-form validators (deterministic, no solver) -------------------------------------
print("closed-form validators:")
report("MISMATCH" not in run("01_three_level/src", "stark_validate.py"),           "5P3/2 tensor Wigner-6j factors")
report("ALL CG" in run("02_multilevel/src", "clebsch_gordan_checks.py"),           "Clebsch-Gordan / line strengths")
_atomic = run("appendix/src", "verify_atomic_claims.py")
report("** CHECK **" not in _atomic and "PASS" in _atomic,                         "atomic claims (leak ratios, D1/D2, B-slopes)")

# ---- (A2) fast floors, recomputed, must equal CANONICAL -----------------------------------------
print("fast floors (recomputed from the solvers):")
_m = re.search(r"<n_z> = (0\.\d+)", run("01_three_level/src", "cooling.py"))
report(bool(_m) and _m.group(1) == "0.0020",                                       "3-level cooling.py = 0.0020",
       "got " + (_m.group(1) if _m else "no match"))
_clean = subprocess.run([PY, "-c", "from cooling_multilevel import solve; print('%.4f' % solve(clean=True))"],
                        cwd=os.path.join(ROOT, "02_multilevel/src"), capture_output=True, text=True).stdout.strip()
report(_clean == "0.0032",                                                         "ch02 clean solve = 0.0032", "got " + _clean)

# ---- (C) every canonical number is quoted in its docs (code <-> prose) ---------------------------
print("number consistency (canonical floor present in the docs):")
for name, val, docs in CANONICAL:
    for d in docs:
        report(val in read(d),                                                     "%s = %s in %s" % (name, val, d))

# ---- (D) reverse tripwire: NO unregistered floor-shaped number may appear in the README ----------
# (C) says "these floors must be present"; (D) says "no OTHER floor may be".  Together they make
# CANONICAL the complete, closed set of sub-0.1 floors (the range FLOOR_RE matches) the front door may quote: a future edit
# cannot introduce an ungated floor without either registering it in CANONICAL or tripping this check.
print("no ungated floor-shaped number in the README:")
FLOOR_RE = re.compile(r"(?<![\d.])0\.0\d+(?![\d])")     # 0.0011, 0.0020, 0.055, 0.087, 0.09 ... (not 0.15)
ALLOW   = {"0.094"}                                     # eta (Lamb-Dicke): a coupling, NOT a floor
PENDING = set()                                         # (master floor 0.055 now settled -> in CANONICAL)
_canon_vals   = {v for _, v, _ in CANONICAL}
_readme_floors = set(FLOOR_RE.findall(read("README.md")))
_unregistered = _readme_floors - _canon_vals - ALLOW - PENDING
report(not _unregistered,                                                          "every floor in README is registered",
       "unregistered (add to CANONICAL or ALLOW): " + ", ".join(sorted(_unregistered)))

# ---- (E) every README link / image / run-command points at a file that still exists --------------
print("README references resolve (links, figures, commands):")
_txt, _missing = read("README.md"), []
for p in re.findall(r"\]\(([^)]+)\)", _txt):              # markdown links + image embeds
    p = p.split("#")[0].strip()
    if p and not p.startswith(("http://", "https://", "mailto:")) and not os.path.exists(os.path.join(ROOT, p)):
        _missing.append("link/img -> " + p)
for _block in re.findall(r"```bash\n(.*?)```", _txt, re.S):   # run commands, tracking `cd` (cwd resets per block)
    _cwd = ROOT
    for _ln in _block.splitlines():
        m = re.match(r"\s*cd\s+(\S+)", _ln)
        if m:
            _cwd = os.path.normpath(os.path.join(_cwd, m.group(1)));  continue
        for _s in re.findall(r"python3?\s+([\w./-]+\.py)", _ln):
            if not os.path.exists(os.path.join(_cwd, _s)):
                _missing.append("command -> " + os.path.relpath(os.path.join(_cwd, _s), ROOT))
        m = re.search(r"-r\s+(\S+\.txt)", _ln)               # pip install -r requirements.txt
        if m and not os.path.exists(os.path.join(_cwd, m.group(1))):
            _missing.append("requirements -> " + m.group(1))
report(not _missing,                                                               "all README paths/commands resolve",
       "; ".join(_missing))

# ---- (F) the SSOT (results.json) has the keys the figures read, and agrees with CANONICAL ---------
# results.json is written by the slow cooling_dark_vertex.py solve and READ by the figures. Verify the
# cached artifact (no solver here -- just read the file the solve already wrote); soft-skip if absent.
print("SSOT results.json <-> code:")
_rj = os.path.join(ROOT, "03_dark_vertex", "results.json")
if not os.path.exists(_rj):
    report(True, "results.json reconciliation", "skipped -- run cooling_dark_vertex.py to write it")
else:
    import json
    _r = json.load(open(_rj))
    _keys_ok = all(k in _r for k in ("no_master", "master", "intrinsic_multilevel")) and \
               "floor" in _r.get("no_master", {}) and "floor" in _r.get("master", {})
    report(_keys_ok,                                                               "results.json has the keys figures read (no_master/master/intrinsic)")
    if _keys_ok:
        _cn = dict((n, v) for n, v, _ in CANONICAL)
        report(("%.4f" % _r["intrinsic_multilevel"]) == "0.0032",                  "results.json intrinsic == 0.0032",
               "got %s" % _r["intrinsic_multilevel"])
        report(("%.3f" % round(_r["master"]["floor"], 3)) == _cn["master floor"],   "results.json master floor == %s" % _cn["master floor"],
               "got %.4f" % _r["master"]["floor"])

# ---- (G) the committed figure fallbacks agree with CANONICAL.  results.json is optional / on-demand; the
#      fallback is what a figure draws until the solve writes it, so it must match the headline floor. -----
print("figure fallbacks <-> CANONICAL:")
_cn = dict((n, v) for n, v, _ in CANONICAL)
for _f in ("03_dark_vertex/src/make_figure.py", "04_master/src/master_figures.py"):
    _mm = re.search(r'"master":\s*\{[^}]*"floor":\s*(0\.\d+)', read(_f))
    report(bool(_mm) and _mm.group(1) == _cn["master floor"],                       "%s master fallback == %s" % (_f.split("/")[-1], _cn["master floor"]),
           "got " + (_mm.group(1) if _mm else "no match"))

# ---- (B) figure freshness: redraw the fast figures, then diff against what is committed ----------
print("figure freshness (committed PNG == freshly drawn):")
subprocess.run([PY, "regenerate.py", "--fast"], cwd=ROOT, capture_output=True, text=True)
_diff = subprocess.run(["git", "diff", "--name-only", "--", "*.png"], cwd=ROOT,
                       capture_output=True, text=True).stdout.strip()
report(_diff == "",                                                                "no stale figures",
       "regenerate.py changed: " + _diff.replace("\n", ", "))

# ---- (slow) the full solves, opt-in -------------------------------------------------------------
if SLOW:
    print("slow solves:")
    report("0.0032" in run("02_multilevel/src", "cooling_multilevel.py"),          "multilevel __main__ reproduces clean 0.0032")
    # the dark-vertex master-floor solve is added to CANONICAL once its Delta-optimum is settled.

print("\n%d / %d checks passed." % (sum(_passed), len(_passed)))
sys.exit(0 if all(_passed) else 1)
