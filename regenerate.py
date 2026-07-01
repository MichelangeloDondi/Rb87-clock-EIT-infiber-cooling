#!/usr/bin/env python3
"""
regenerate.py -- redraw every figure in the repo from the current code + config, in one command.

Every figure here is a PNG produced by a script, not a hand-drawn image. So the way to keep the
figures honest is simple: after changing a physical number in any chapter's config.py (a detuning,
a Rabi frequency, the magnetic field), run this. Every figure that depends on that number is
redrawn from the code, so a committed figure can never silently disagree with the physics that made
it. `run_tests.py` then checks that nothing was left stale.

Usage:
    python regenerate.py                 # every chapter
    python regenerate.py 01 02           # only the listed chapters
    python regenerate.py --fast          # skip the slow qutip re-solves (redraw from a PRIOR full run's results.json)

Chapters 03 and 04 quote numbers that come from a qutip solve (minutes at N_fock=6). By default this
re-runs that solve so the figure matches it; --fast skips it and redraws from the cached results.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

# chapter -> (src dir, scripts to run IN ORDER). A solver listed before a figure recomputes the
# numbers the figure then draws (via the chapter's results file); figure-only scripts just redraw.
CHAPTERS = [
    ("01",       "01_three_level/src", ["plots.py"]),
    ("02",       "02_multilevel/src",  ["level_scheme.py"]),
    ("03",       "03_dark_vertex/src", ["cooling_dark_vertex.py", "make_figure.py"]),
    ("04",       "04_master/src",      ["master_figures.py"]),
    ("appendix", "appendix/src",       ["make_figure.py"]),
]
SLOW_SOLVERS = {"cooling_dark_vertex.py"}   # the ones --fast skips


def main(argv):
    fast = "--fast" in argv
    wanted = [a for a in argv if not a.startswith("--")]
    for name, srcdir, scripts in CHAPTERS:
        if wanted and name not in wanted:
            continue
        for script in scripts:
            if fast and script in SLOW_SOLVERS:
                print("[%s] %s  -- skipped (--fast)" % (name, script), flush=True)
                continue
            print("[%s] %s ..." % (name, script), flush=True)
            subprocess.run([sys.executable, script], cwd=os.path.join(ROOT, srcdir), check=True)
    print("\ndone -- every requested figure has been regenerated from the current code.")


if __name__ == "__main__":
    main(sys.argv[1:])
