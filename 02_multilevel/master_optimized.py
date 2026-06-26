"""
master_optimized.py -- run the OPTIMIZED master configuration.

The coherent-solver optimization (documented in ../03_master/) PROPOSED the master as a WEAK, FAR-DETUNED
F2->F'1 'tickle' with the EOM comb byproducts SUPPRESSED, claiming a floor ~0.023. This driver runs that same
operating point in THIS repo's engine (the detuned master is modeled as an incoherent F2->F'1 repumper, valid
at its low saturation) as an independent check.

  *** RESULT: this engine does NOT reproduce 0.023 -- it gives ~0.41 (comb off) / ~1.0 (comb on). ***

The ~18x gap is the F=1-recycling treatment: this engine uses the real off-resonant probe (F=1 dark ~0.36),
the coherent solver used a phenomenological F=1 repump (F=1 dark ~0.02). What both agree on: the master clears
|2,-2>, F=1 is the dominant residual, and the comb byproducts hurt. See ../03_master/README.md ("Optimized
master configuration") for the full triage. Operating point + master knobs are in config.py (master_*).

Run:  python master_optimized.py     (a few minutes; needs numpy, qutip, sympy)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as c
import cooling_multilevel as m

OP = dict(d2=c.master_d2, Delta=c.master_Delta, OmR=c.master_OmR, master=True)

def run(tag, **kw):
    R = m.solve(want=True, **dict(OP, **kw))
    pops = R['pops']
    f1 = sum(w for (F, mm), w in pops.items() if F == 1 and (F, mm) != (1, -1))
    print("  %-24s <n_z> = %.4f   P0 = %.3f   |2,-2> = %.4f   F=1 dark = %.3f"
          % (tag, R['nbar'], R['pn'][0], pops.get((2, -2), 0.0), f1), flush=True)
    return R['nbar']

if __name__ == "__main__":
    print("OPTIMIZED master config:  Delta=%.0f  OmR=%.2f  d2=%.2f  |  master: det=%.0f off F'1, Om=%.1f"
          % (c.master_Delta, c.master_OmR, c.master_d2, c.master_det, c.master_Om))
    run("clean (comb OFF)", rep_scale=0.0)    # the optimum: dedicated master, comb byproducts suppressed
    run("with comb byproducts", rep_scale=1.0)  # leftover EOM tones ON -- they hurt at this operating point
    print("DONE")
