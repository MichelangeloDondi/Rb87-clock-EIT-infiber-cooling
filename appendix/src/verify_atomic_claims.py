"""
verify_atomic_claims.py -- recompute, from first principles, the atomic facts the appendix relies on.

Everything here is Clebsch-Gordan x signed hyperfine-6j (no solver, no fitting): the D1-vs-D2 leak
comparison, the selection rules behind the "no m=0 pi Lambda" and the |2,0> pi knob, and the magnetic
field slopes. Run it and read the PASS/-- column against the numbers quoted in the appendix notes.

  87Rb: I=3/2, ground 5S1/2 J=1/2 (g_F: F=1 -> -1/2, F=2 -> +1/2), muB = 1.399624 MHz/G.
  D2 = 5P3/2 (J'=3/2), F'=0..3, F'1-F'2 splitting 156.94 MHz, Gamma 6.07.
  D1 = 5P1/2 (J'=1/2), F'=1,2,  F'1-F'2 splitting 816.66 MHz (=2*A_{P1/2}), Gamma 5.746.

Run:  python verify_atomic_claims.py   (needs numpy, sympy)
"""
import numpy as np
from sympy import S
from sympy.physics.wigner import wigner_6j, clebsch_gordan

clebsch = lambda F, m, q, Fp, mp: float(clebsch_gordan(S(F), 1, S(Fp), S(m), S(q), S(mp)))   # <F m; 1 q | F' m'>
SPLIT = {"D2": 156.94, "D1": 816.66}      # F'1-F'2 hyperfine spacing (2pi MHz)
GAMMA = {"D2": 6.07,  "D1": 5.746}
JP    = {"D2": S(3)/2, "D1": S(1)/2}      # excited fine-structure J'
muB   = 1.399624                          # MHz/G

def line_amplitude(F, Fp, line):                     # signed hyperfine line amplitude (carries the 6j sign)
    return float(np.sqrt(2*Fp + 1) * wigner_6j(S(1)/2, JP[line], 1, S(Fp), S(F), S(3)/2))
def reduced_ratio(F, line):                # <F||d||F'1> / <F||d||F'2>, signed
    return line_amplitude(F, 1, line) / line_amplitude(F, 2, line)
def R_legs(line):                          # per-leg F'1/F'2 coupling ratio for the clock pair
    Rc = reduced_ratio(2, line) * (clebsch(2,  1, -1, 1, 0) / clebsch(2,  1, -1, 2, 0))   # control sigma- |2,+1>
    Rp = reduced_ratio(1, line) * (clebsch(1, -1,  1, 1, 0) / clebsch(1, -1,  1, 2, 0))   # probe   sigma+ |1,-1>
    return Rc, Rp

# check(label, computed, want): `want` = the target value quoted in the appendix .md tables / prose (the docs these claims back).
def check(label, got, want, tol=2e-3):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print("  %-52s %+11.5f   (want %+.4f)  %s" % (label, got, want, "PASS" if ok else "** CHECK **"))

def claim1_darkstate_leak_ratios():
    """1. The dark-state F'1 leak: per-leg coupling ratios R_c (control), R_p (probe)."""
    print("1. The dark-state F'1 leak: per-leg coupling ratios R_c (control), R_p (probe)")
    for line in ("D2", "D1"):
        Rc, Rp = R_legs(line)
        check("%s  R_c = control |2,+1> -> F'1/F'2" % line, Rc, -0.3464 if line == "D2" else -0.7746)
        check("%s  R_p = probe   |1,-1> -> F'1/F'2" % line, Rp, +1.7321 if line == "D2" else +0.7746)
        print("       |R_c - R_p| = %.4f   (leak amplitude; dark-on-F'2 is exactly 0)" % abs(Rc - Rp))
    Rc2, Rp2 = R_legs("D2"); Rc1, Rp1 = R_legs("D1")
    check("|R_c-R_p| ratio D1/D2 (leak NOT cancelled on D1)", abs(Rc1-Rp1)/abs(Rc2-Rp2), 0.745)


def claim2_detuning_not_amplitude():
    """2. Why D1's leak is small: it is the DETUNING, not the amplitude, not the branching."""
    print("2. Why D1's leak is small: it is the DETUNING, not the amplitude, not the branching")
    check("F'1-F'2 splitting ratio D1/D2", SPLIT["D1"]/SPLIT["D2"], 5.20)
    def leak_figure_of_merit(line, D):                          # leak figure of merit ~ Gamma (R_c-R_p)^2 / (Delta+split)^2
        Rc, Rp = R_legs(line); return GAMMA[line]*(Rc-Rp)**2 / (D + SPLIT[line])**2
    for D in (0.0, 25.0):
        check("leak suppression D2/D1 at Delta=%2.0f" % D, leak_figure_of_merit("D2", D)/leak_figure_of_merit("D1", D), 51.5 if D == 0 else 40.7)
    Rc2, Rp2 = R_legs("D2"); Rc1, Rp1 = R_legs("D1")
    check("  detuning factor (split_D1/split_D2)^2", (SPLIT["D1"]/SPLIT["D2"])**2, 27.08)
    check("  amplitude factor ((R_c-R_p)_D2/(R_c-R_p)_D1)^2",  (abs(Rc2-Rp2)/abs(Rc1-Rp1))**2, 1.800)


def claim3_no_pi_lambda_and_gem():
    """3. The pure m=0 pi pair has NO Lambda (so the scheme must use m=+-1)."""
    print("3. The pure m=0 pi pair has NO Lambda (so the scheme must use m=+-1)")
    reach = {g: [Fp for Fp in (0,1,2,3) if abs(Fp-g) <= 1 and abs(clebsch(g, 0, 0, Fp, 0)) > 1e-9] for g in (1, 2)}
    print("  |1,0> via pi reaches F' =", reach[1], "  (clebsch(|1,0>->F'1) = %.3f)" % clebsch(1,0,0,1,0))
    print("  |2,0> via pi reaches F' =", reach[2], "  (clebsch(|2,0>->F'2) = %.3f)" % clebsch(2,0,0,2,0))
    print("  common F' =", sorted(set(reach[1]) & set(reach[2])), " -> no shared excited state -> no dark state")
    print("  GEM: |2,0>->|F'1,0> pi = %+.4f (ALLOWED, pure F'1 knob);  |2,0>->|F'2,0> pi = %+.4f (FORBIDDEN)"
          % (clebsch(2,0,0,1,0), clebsch(2,0,0,2,0)))


def claim4_magnetic_slopes():
    """4. Magnetic field slopes (g_F*m*muB)."""
    print("4. Magnetic field slopes (g_F*m*muB)")
    gF_times_m = {(1,-1): -0.5*-1, (2,+1): 0.5*1, (2,-1): 0.5*-1, (2,-2): 0.5*-2, (2,0): 0.5*0}
    check("clock pair two-photon slope: gF_times_m(2,+1)-gF_times_m(1,-1)", (gF_times_m[(2,1)]-gF_times_m[(1,-1)])*muB, 0.0)
    check("stretched axial baseline |1,-1>+|2,-1> slope (MHz/G)", (gF_times_m[(1,-1)]-gF_times_m[(2,-1)])*muB, 1.40)
    check("stretched tilted-B (pi) |1,-1>+|2,-2| slope (MHz/G)", (gF_times_m[(1,-1)]-gF_times_m[(2,-2)])*muB, 2.10)
    check("|2,0> knob coherence drift vs clock dark (MHz/G)", abs(gF_times_m[(2,0)] - 0.5)*muB, 0.70)


if __name__ == "__main__":
    print(__doc__.split("Run ")[0])
    print("=" * 84)
    claim1_darkstate_leak_ratios()
    print("\n" + "=" * 84)
    claim2_detuning_not_amplitude()
    print("\n" + "=" * 84)
    claim3_no_pi_lambda_and_gem()
    print("\n" + "=" * 84)
    claim4_magnetic_slopes()
    print("\nDONE.")
