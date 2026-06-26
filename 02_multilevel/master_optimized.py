"""
master_optimized.py -- run the OPTIMIZED master configuration.

The coherent-solver optimization (documented in ../upgrades/) found the master should be a WEAK, FAR-DETUNED
F2->F'1 'tickle' with the EOM comb byproducts SUPPRESSED. The optimum master is detuned (low saturation), so
it is modeled here as an incoherent F2->F'1 repumper in the same engine as the baseline (valid where the
coherent on-resonance master was not). Operating point + master knobs live in config.py (master_*).

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
