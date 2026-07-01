"""
03_dark_vertex -- fold the |F'1,0> leak into the solve, and compute the master floor.

The cooling pair |2,+1>, |1,-1> is two-photon resonant on EVERY m'=0 excited state it can reach --
|F'2,0> (the intended vertex) AND |F'1,0> (157 MHz below it). Two-photon resonance is fixed by the
two GROUND energies and the two laser frequencies; it does not care which excited state sits in
between. So the EIT dark state

        |D2> = (Op |2,+1> - Oc |1,-1>) / N            (dark on F'2 by construction)

is NOT dark on F'1: it keeps a residual coupling onto |F'1,0>, which decays 5/6 to F=1 and loads the
F=1 dark states. That leak -- not repumper inefficiency, not the master placement -- is the dominant,
master-proof floor term, and it is what sets the master floor at ~0.055 (computed here).

This module reuses 02_multilevel/src/cooling_multilevel.py verbatim. Its `with_e1` already carries the
F'1 spoiler edges COHERENTLY at the cooling-Lambda frequency (so the rotating frame still closes,
frame_conflict=0) -- which is why chapter 02's real-delivery headline was ~0.09 and not the clean-Lambda ~0.003. This module adds two things:

  (1) a DEDICATED master, a detuned F2->F'1 sigma+ repumper. Detuned (not parked on F'1) so the
      incoherent-rate model stays valid (det >> 3*Gamma, outside the engine's near-resonance guard) --
      which is also what the physics wants: the master should be the lightest touch that clears the
      one F=2 dark sublevel |2,-2>, not a sledgehammer on resonance.

  (2) a cancellation knob `fscale` that scales the PROBE's coherent |1,-1>->|F'1,0> edge, to test
      whether nulling the leak (constructive dark-state engineering) recovers the floor.

Everything else -- recoil, branching, tensor Stark, the Lamb-Dicke displacement -- is inherited
unchanged. All floors are single-atom, on-axis, radially-localized best cases (chapter 02's scope).

Run `python cooling_dark_vertex.py` to reproduce the chapter numbers (needs numpy, qutip, sympy).
"""

import os, sys, numpy as np, warnings
warnings.simplefilter("ignore")
import config as c                        # chapter-03 config (reuses chapter 02's numbers; edit config.py to explore)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "02_multilevel", "src"))
import cooling_multilevel as m            # the chapter-02 solver, untouched (it picks up the config above)
from sympy import S
from sympy.physics.wigner import wigner_6j, clebsch_gordan


# ----------------------------------------------------------------------------------
# 1.  The atomic ratios that set the leak.   <F||d||F'1> / <F||d||F'2>  (signed),
#     then the per-leg F'1/F'2 coupling ratio R = (reduced ratio)*(CG ratio).
# ----------------------------------------------------------------------------------
def clebsch(F, mq, q, Fp, mp):
    return float(clebsch_gordan(S(F), 1, S(Fp), S(mq), S(q), S(mp)))

def wigner6j(F, Fp):
    # {1/2 3/2 1 ; F' F 3/2}  -- the D2 hyperfine 6j (J=3/2 excited, J=1/2 ground, I=3/2)
    return float(wigner_6j(S(1)/2, S(3)/2, 1, S(Fp), S(F), S(3)/2))

def reduced_ratio(F):
    """<F||d||F'=1> / <F||d||F'=2>, SIGNED (the 6j sign is physical: the two legs interfere)."""
    num = np.sqrt((2*1 + 1)) * wigner6j(F, 1)
    den = np.sqrt((2*2 + 1)) * wigner6j(F, 2)
    return num / den

# control sigma- leg |2,+1>->F'(m'=0); probe sigma+ leg |1,-1>->F'(m'=0)
R_c = reduced_ratio(2) * (clebsch(2, 1, -1, 1, 0) / clebsch(2, 1, -1, 2, 0))   # ~ -0.346
R_p = reduced_ratio(1) * (clebsch(1, -1, 1, 1, 0) / clebsch(1, -1, 1, 2, 0))   # ~ +1.732
null_scale = R_c / R_p                                                     # probe scale that nulls the residual


# ----------------------------------------------------------------------------------
# 2.  Wrap the chapter-02 beam builders:  (a) scale the probe's F'1 edge by fscale,
#     (b) append a dedicated detuned master (F2->F'1, sigma+).
# ----------------------------------------------------------------------------------
_orig_beams = m.beams
_orig_repump = m.repump_beams
override_state = {"master": None, "fscale": 1.0}        # set by floor(); master = dict(det=, rabi=)

def _beams(Oc, Op, d2, with_e1, with_e3):
    bm = _orig_beams(Oc, Op, d2, with_e1, with_e3)
    if with_e1 and override_state["fscale"] != 1.0:
        for b in bm:
            if b["tag"] == "probe":
                b["edges"] = [(g, e, (cc * override_state["fscale"] if e == (1, 0) else cc))
                              for (g, e, cc) in b["edges"]]
    return bm

def _repump(Oc, Op, d2, repump_scale, shift=-1, tag_shift=None):
    bm = _orig_repump(Oc, Op, d2, repump_scale, shift, tag_shift)
    M = override_state["master"]
    if M is not None:
        # dedicated master: sigma+ on F=2->F'1 (the ladder also picks up the far-off F'2,F'3),
        # referenced to F2->F'1 so `det` is the master's red detuning from that line.
        bm.append(dict(edges=m.dipole_edges(2, +1, [1, 2, 3]), named=((2, 0), (1, 1)),
                       det=M["det"], Rabi=M["rabi"], kdir=+1, tag="master"))
    return bm

m.beams = _beams
m.repump_beams = _repump


# ----------------------------------------------------------------------------------
# 3.  Floor.  Minimize over a few two-photon detunings (the cooling sideband sits near
#     -nu_z; the multilevel a.c. shift pulls the exact spot).
# ----------------------------------------------------------------------------------
def floor(Delta=45.0, with_leak=True, comb=0.0, master=None, fscale=1.0,
          N_fock=6, d2s=(-0.16, -0.10, -0.04, 0.0), return_scan=False):
    """The cooling floor with the two-photon detuning delta2 AUTO-OPTIMISED (= min over d2s).

    Delta       : control detuning above |F'2,0> (in-trap), 2pi*MHz
    with_leak   : keep the F'1/F'3 spoiler edges (the honest model).  False = the old 'perfect dark'.
    comb        : amplitude of the leftover-comb repumpers (1.0 = single-EOM chain; 0 = suppressed)
    master      : None, or dict(det=<red detuning from F2->F'1>, rabi=<Rabi>) for the dedicated master
    fscale      : multiply the probe's |1,-1>->|F'1,0> edge (cancellation knob; 1.0 = physical)
    return_scan : also return the [(delta2, <n_z>), ...] scan (the delta2-vs-floor report curve)
    """
    c.Delta = float(Delta)
    override_state["master"] = master
    override_state["fscale"] = fscale
    scan = [(d, m.solve(d2=d, repump_scale=comb, with_e1=with_leak, with_e3=with_leak, N_fock=N_fock)) for d in d2s]
    nbar = min(n for _, n in scan)
    return (nbar, scan) if return_scan else nbar


# ----------------------------------------------------------------------------------
# 4.  The report: compute every floor (delta2 auto-optimised), find the master's Delta optimum,
#     and CACHE it to results.json -- the single source of truth that make_figure.py and chapter 04
#     read, so no floor is ever hard-coded. Re-run after changing a parameter to find the new optimum.
# ----------------------------------------------------------------------------------
MASTER = dict(det=-60.0, rabi=1.2)       # the dedicated master: light, detuned (det >> 3*Gamma)

def report(write=True):
    """Compute the chapter's floors, all with delta2 AUTO-OPTIMISED (min over the d2 scan), find the
    master's best Delta, and cache the lot to ../results.json. Returns the dict. Slow (~min) at N_fock=6."""
    import json
    master_by_delta = {D: floor(D, True, 0.0, MASTER) for D in (25, 30, 45)}     # master floor vs Delta
    m_delta = min(master_by_delta, key=master_by_delta.get)
    d2_grid = tuple(round(x, 3) for x in (-0.22, -0.18, -0.14, -0.10, -0.06, -0.02))
    _, d2_scan = floor(m_delta, True, 0.0, MASTER, d2s=d2_grid, return_scan=True)  # the delta2 servo curve
    cancel = {f: floor(45, True, 0.0, MASTER, fscale=f) for f in (1.0, 0.5, 0.0)}  # probe-F'1 cancellation
    res = {
        "R_c": round(R_c, 3), "R_p": round(R_p, 3), "null_scale": round(null_scale, 3),
        "no_master":     {"delta": 45,      "floor": round(floor(45, True, 1.0, None), 4)},
        "master":        {"delta": m_delta, "floor": round(master_by_delta[m_delta], 4)},
        "master_d45":    round(master_by_delta[45], 4),
        "no_leak_ideal": {"delta": 45,      "floor": round(floor(45, False, 0.0, MASTER), 4)},
        "intrinsic_multilevel": 0.0032,
        "master_delta_scan": {str(D): round(v, 4) for D, v in master_by_delta.items()},
        "master_d2_scan":    [[d, round(n, 4)] for d, n in d2_scan],
        "cancellation_scan": {("%.2f" % f): round(v, 4) for f, v in cancel.items()},
    }
    if write:
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results.json")
        json.dump(res, open(path, "w"), indent=2)
    return res


if __name__ == "__main__":
    res = report()   # computes + caches results.json
    print(__doc__.split("Run ")[0])
    print("=" * 82)
    print("1.  The leak is atomic.  Signed F'1/F'2 coupling ratios of the two dark-state legs:")
    print("      control |2,+1>:  R_c = % .3f ;  probe |1,-1>:  R_p = % .3f  (probe couples %.1fx stronger to F'1)"
          % (res["R_c"], res["R_p"], abs(res["R_p"])))
    print("      R_c != R_p  =>  the dark state cannot be dark on BOTH F'2 and F'1: that residual is the leak.")
    print("\n2.  The master floor (delta2 auto-optimised):")
    print("      minimal single-EOM chain (comb + leak), Delta=%s :  n_z = %.4f" % (res["no_master"]["delta"], res["no_master"]["floor"]))
    print("      + dedicated master, same Delta=45                :  n_z = %.4f" % res["master_d45"])
    print("      + dedicated master, leak-aware Delta=%-2s          :  n_z = %.4f   <-- the headline" % (res["master"]["delta"], res["master"]["floor"]))
    print("      no-leak ideal (F'1 leak cancelled)               :  n_z = %.4f   (master-proof gap = the leak)" % res["no_leak_ideal"]["floor"])
    print("\n2b. master floor vs Delta (delta2 auto-optimised at each -- the leak wants smaller Delta):")
    for D in sorted(res["master_delta_scan"], key=int):
        print("      Delta = %2s   n_z = %.4f" % (D, res["master_delta_scan"][D]))
    print("\n2c. floor vs the two-photon detuning delta2 at the optimal Delta=%s  (the servo curve):" % res["master"]["delta"])
    for d, n in res["master_d2_scan"]:
        print("      delta2 = %+.2f   n_z = %.4f" % (d, n))
    print("\n  All cached to results.json -- make_figure.py and chapter 04 read it, so no floor is hard-coded.")
