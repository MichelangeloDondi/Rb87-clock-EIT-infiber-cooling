"""
03_dark_vertex -- fold the |F'1,0> leak into the solve, and compute the master floor.

The cooling pair |2,+1>, |1,-1> is two-photon resonant on EVERY m'=0 excited state it can reach --
|F'2,0> (the intended vertex) AND |F'1,0> (157 MHz below it). Two-photon resonance is fixed by the
two GROUND energies and the two laser frequencies; it does not care which excited state sits in
between. So the EIT dark state

        |D2> = (Op |2,+1> - Oc |1,-1>) / N            (dark on F'2 by construction)

is NOT dark on F'1: it keeps a residual coupling onto |F'1,0>, which decays 5/6 to F=1 and loads the
F=1 dark states. That leak -- not repumper inefficiency, not the master placement -- is the dominant,
master-proof floor term, and it is what sets the master floor at ~0.055 (chapter 04's design only
called it a target; here it is computed).

This module reuses 02_multilevel/cooling_multilevel.py verbatim. Its `with_e1` already carries the
F'1 spoiler edges COHERENTLY at the cooling-Lambda frequency (so the rotating frame still closes,
frame_conflict=0) -- which is why chapter 02's headline was ~0.10 and not ~0.003. This module adds two things:

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
          N_fock=6, d2s=(-0.16, -0.10, -0.04, 0.0)):
    """
    Delta     : control detuning above |F'2,0> (in-trap), 2pi*MHz
    with_leak : keep the F'1/F'3 spoiler edges (the honest model).  False = the old 'perfect dark'.
    comb      : amplitude of the leftover-comb repumpers (1.0 = single-EOM chain; 0 = suppressed)
    master    : None, or dict(det=<red detuning from F2->F'1>, rabi=<Rabi>) for the dedicated master
    fscale    : multiply the probe's |1,-1>->|F'1,0> edge (cancellation knob; 1.0 = physical)
    """
    c.Delta = float(Delta)
    override_state["master"] = master
    override_state["fscale"] = fscale
    return min(m.solve(d2=d, repump_scale=comb, with_e1=with_leak, with_e3=with_leak, N_fock=N_fock)
               for d in d2s)


# ----------------------------------------------------------------------------------
# 4.  Reproduce the chapter.
# ----------------------------------------------------------------------------------
if __name__ == "__main__":
    MASTER = dict(det=-60.0, rabi=1.2)       # the dedicated master: light, detuned (det >> 3*Gamma)

    print(__doc__.split("Run ")[0])
    print("=" * 82)
    print("1.  The leak is atomic.  Signed F'1/F'2 coupling ratios of the two dark-state legs:")
    print("      control  |2,+1>:  R_c = % .3f" % R_c)
    print("      probe    |1,-1>:  R_p = % .3f   <-- probe couples %.2fx STRONGER to F'1 than to F'2"
          % (R_p, abs(R_p)))
    print("      R_c != R_p  =>  |D2> cannot be dark on BOTH F'2 and F'1.  The residual is the leak.")
    print("      (the probe scale that would null the residual:  f = R_c/R_p = % .3f)" % null_scale)

    print("\n" + "=" * 82)
    print("2.  The master floor  (master = a detuned dedicated F2->F'1 sigma+ repumper):\n")
    rows = [
        ("minimal single-EOM chain (comb only), leak ON ", floor(45, True,  1.0, None)),
        ("minimal single-EOM chain (comb only), leak OFF", floor(45, False, 1.0, None)),
        ("master config (comb suppressed),       leak ON ", floor(45, True,  0.0, MASTER)),
        ("master config (comb suppressed),       leak OFF", floor(45, False, 0.0, MASTER)),
        ("master config, leak ON, Delta = 30            ", floor(30, True,  0.0, MASTER)),
    ]
    for label, val in rows:
        print("      %-46s  n_z = %.4f" % (label, val))
    print("\n      => at fixed Delta=45 the master earns 0.088 -> 0.082 (clears |2,-2> + comb scatter); the")
    print("         LEAK (0.082 -> 0.029 with the leak off) is the ~0.05 that remains, and it is master-proof.")

    print("\n" + "=" * 82)
    print("2b. The game-changer test: best floor vs Delta (leak ON), no-master chain vs the master:\n")
    print("      Delta   no-master (comb)   + master")
    for D in (25, 30, 35, 45):
        print("      %4.0f       %.4f           %.4f" % (D, floor(D, True, 1.0, None), floor(D, True, 0.0, MASTER)))
    print("\n      => WITHOUT the master the comb chain is best at large Delta (~0.087, repump-limited: its")
    print("         comb tones sit near F'2 and can't follow the leak to small Delta). WITH the master the")
    print("         leak-favoured small Delta~25 is reachable (~0.055, leak-limited). So the master is a real")
    print("         but modest upgrade (0.087 -> 0.055), NOT a route to the 0.0032 ideal.")

    print("\n" + "=" * 82)
    print("3.  Can the leak be cancelled?  Scale the probe's F'1 edge (constructive dark engineering):\n")
    for f in (1.0, 0.5, 0.0, round(null_scale, 2), -0.5):
        tag = "   <- residual nulled" if abs(f - null_scale) < 0.02 else ""
        print("      probe-F'1 scale f = %+5.2f   n_z = %.4f%s" % (f, floor(45, True, 0.0, MASTER, fscale=f), tag))
    print("      no-leak ideal (with_e1 off):              n_z = %.4f" % floor(45, False, 0.0, MASTER))
    print("\n      Suppressing the dominant probe-F'1 term recovers ~0.08 -> ~0.04, toward the no-leak")
    print("      ideal.  But f<1 is not a knob: R_c, R_p are atomic constants, a resonant canceller adds")
    print("      a 3rd frequency to |F'1,0> (breaks the static frame), and a co-propagating tone at the")
    print("      probe frequency only rescales Op.  Genuine cancellation needs a time-dependent (Floquet)")
    print("      co-propagating tone -- a further frontier, outside this static-frame solver.")
