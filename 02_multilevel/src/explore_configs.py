"""
explore_configs.py -- which single-EOM / one-AOM configuration repumps best?

The repumpers are the leftover comb tones, so their detunings are FIXED by the EOM frequency f_mod
and the tag-AOM shift 2f_A (and which AOM order is used). The catch: ONE AOM moves the forward
sideband (repump1, on F=1) and the retro carrier (repump2, on F=2) in OPPOSITE directions --

    repump1 detuning (from F=1->F'2) = Delta + d2 - shift*2f_A
    repump2 detuning (from F=2->F'2) = Delta +        shift*2f_A      (shift = -1 down, +1 up)

-- so you cannot pull both onto the useful F'1 line at once. And the "other order" (up-shift)
configs push the tones toward F'0 (for repump1) and F'3 (for repump2), which decay BACK to the same
hyperfine (m-shuffle + recoil heating, no inter-F repump). This sweep computes the actual floor.

Run:  python explore_configs.py     (needs numpy, qutip, sympy; a few steady-state solves -> minutes)
"""
import sys
import numpy as np
import cooling_multilevel as m
import config as c

A_HFS = m.A_HFS
d2, N_fock = -0.10, 6

# (label, AOM order shift, 2f_A)   -- f_mod = A_HFS - shift*2f_A is the implied EOM frequency
CONFIGS = [
    ("A baseline (down, 200 MHz AOM, 2fA=400)", -1, 400),
    ("B up-shift (up,   200 MHz AOM, 2fA=400)", +1, 400),
    ("C up-shift (up,   150 MHz AOM, 2fA=300)", +1, 300),
    ("D near-res (down, 101 MHz AOM, 2fA=202)", -1, 202),
]


def out(s):
    print(s); sys.stdout.flush()


if __name__ == "__main__":
    out("settled multilevel solve, INCOHERENT off-resonant repumpers (frame-consistent at any power).")
    out(f"  d2={d2:+.2f}, repump_scale=1 (chain-natural power), N_fock={N_fock}")
    out(f"  clean 3-level Lambda floor = {m.solve(clean=True, N_fock=N_fock):.4f}\n")
    out("  rep1 = fwd sideband (F=1, sigma-): useful via F'1 (1/6->F2) or F'2 (1/2->F2); F'0 decays F1-only.")
    out("  rep2 = retro carrier (F=2, sigma+): useful via F'1 (5/6->F1) or F'2 (1/2); F'3 decays F2-only.\n")
    hdr = "%-42s f_mod   rep1: F'2/F'1/F'0    rep2: F'2/F'1/F'3    floor    dark  P(n=0)"
    out(hdr % "config")
    for (name, s, tfa) in CONFIGS:
        fmod = A_HFS - s * tfa
        r1 = c.Delta + d2 - s * tfa            # rep1 detuning from F=1->F'2
        r2 = c.Delta + s * tfa                 # rep2 detuning from F=2->F'2
        R = m.solve(d2=d2, repump_scale=1.0, shift=s, tag_shift=tfa, N_fock=N_fock, want=True)
        dark = sum(w for g, w in R['pops'].items() if g not in ((1, -1), (2, 1)))
        # 157/229/267 = the |F'2-F'1|/|F'2-F'0|/|F'3-F'2| hyperfine spacings (2pi MHz; cf. m.excited_hf_spacings)
        # a tone within ~3*Gamma of any line breaks the incoherent-rate model -> this row's floor is NOT trustworthy
        mind = min(abs(x) for x in (r1, r1 + 157, r1 + 229, r2, r2 + 157, r2 - 267))
        edge = "  (P_edge=%.0e)" % R['pn'][-1] if R['pn'][-1] > 1e-2 else ""
        flag = "  << rate model INVALID: a tone sits %.0f MHz off a line (treat coherently)" % mind if mind < 3 * m.GAMMA else ""
        out("%-42s %5.0f  %+5.0f/%+5.0f/%+5.0f    %+5.0f/%+5.0f/%+5.0f    %.4f   %.2f  %.3f%s%s"
            % (name, fmod, r1, r1 + 157, r1 + 229, r2, r2 + 157, r2 - 267,
               R['nbar'], dark, R['pn'][0], edge, flag))
    out("\n  (detunings in 2pi*MHz, signed, from each line. 'dark' = population outside |1,-1>,|2,+1>.")
    out("   F'0/F'3 columns are the intra-F decay-back lines -- a tone sitting on them heats, does not repump.)")
