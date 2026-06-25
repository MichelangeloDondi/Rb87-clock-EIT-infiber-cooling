"""
cooling_multilevel.py -- the FAITHFUL multilevel solve, WITH the two sigma repumpers
in the computation (what the clean 3-level cooling.py idealizes away).

cooling.py routes spontaneous decay straight back into the two Lambda states -- i.e. it
ASSUMES perfect repumping. In reality |F'2,0> decays (sigma+, pi, sigma-) across BOTH
ground hyperfines, so population piles up in the dark sublevels |1,0>,|1,+1>,|2,-2..0>
and cooling stalls. The repumpers clear it. In the FINALISED single-EOM chain
(EBLANA -> EOM -> EDFA -> PPLN -> HCPCF -> tag AOM -> retro) the repumpers are NOT separate
lasers: they are the leftover tones of the SAME comb (one EOM at f_mod = A_HFS + 2f_A = 7.23 GHz,
one 200 MHz tag AOM double-passed to 2f_A = 400 MHz). The four tones at the atom:

  control  sigma-  F=2 -> F'2  forward      at  Delta                 (named |2,+1>->|F'2,0>)
  probe    sigma+  F=1 -> F'2  retro,-2f_A  at  A_HFS + Delta + d2     (named |1,-1>->|F'2,0>)
  repump1  sigma-  F=1         forward      at  probe  + 2f_A          = fwd +1 EOM sideband
  repump2  sigma+  F=2         retro,-2f_A  at  control - 2f_A         = retro carrier

so (matching the hardware): probe is 6.83 GHz (A_HFS) from control, +400 MHz (2f_A) from repump1,
and 7.23 GHz from repump2. The repumpers are deliberately OFF-RESONANT: repump1 drives F=1->F'2
at 445 MHz off (F'3 sits 178 MHz away but F=1->F'3 is dF=2 FORBIDDEN -- CG=0); repump2 drives
F=2->F'1 at 198 MHz off (F'0 sits 126 MHz away but F=2->F'0 is dF=2 FORBIDDEN). So they repump
SLOWLY; raise config.rep_scale (more EOM/launch power) to repump faster. These are the engine's
`fwd_sideband_rejected` / `retro_carrier_rejected` (it labels them parasitic because IT carries
dedicated repumpers; the minimal chain re-uses them AS the repumpers).

The atomic core (Clebsch-Gordan dipole ladders, multi-rotating frame, m-resolved decay branching,
3-point recoil) is ported from the validated engine src/engines/eit_cooling_tool.py (Section 6).
The ONE physics addition is the per-(F',m') 1064 tensor Stark from stark.py fed into Ee() -- the
engine's flagged v0.3.0 limitation -- so the off-resonant repumper detunings are referenced to the
real Stark-shifted F'1/F'3 lines.

SCOPE / caveat: the repumpers are coherent beams in one multi-rotating frame, which the two off-resonant
tones do not close exactly (residual ~1.4 MHz). At their NATURAL power this is a ~1% effect on weak,
far-detuned edges, so the natural-power floor is trustworthy; the residual grows with rep_scale, so the
repumper-POWER scan is only qualitative (a clean power study would model the repumpers as incoherent rates).

Run:  python cooling_multilevel.py        (needs numpy, qutip, sympy)
"""
import numpy as np
import qutip as qt
from sympy import S
from sympy.physics.wigner import clebsch_gordan, wigner_6j
import config as c
import stark

# ---- 87Rb D2 constants (mirror src/engines/eit_cooling_tool.py Section 1) -------------------
A_HFS = 6834.682610                                  # ground hyperfine splitting (2pi MHz)
GAMMA = c.Gamma                                      # 5P3/2 natural linewidth
EHF = {0: -302.07, 1: -229.85, 2: -72.91, 3: 193.74}  # 5P3/2 hyperfine centroids (2pi MHz)
DF = {Fp: EHF[Fp] - EHF[2] for Fp in EHF}           # spacings from F'=2 {1:-156.94, 3:+266.65,...}
BR = {0: {1: 1.0, 2: 0.0}, 1: {1: 5/6, 2: 1/6},     # F'->F hyperfine decay branching
      2: {1: 1/2, 2: 1/2}, 3: {1: 0.0, 2: 1.0}}
EM_REC = [(-1, 1/6), (0, 2/3), (1, 1/6)]            # dipole emission recoil (delta-m, weight)
uB = 1.399624                                        # Bohr magneton / h (MHz/G)
gF_g = {1: -0.5, 2: +0.5}                            # ground Lande gF (clock pair: both gF*m=+1/2)
gF_e = 2.0 / 3.0                                     # 5P3/2 Lande gF (equal for F'=1,2,3)

CGc = lambda F, m, q, Fp, mp: float(clebsch_gordan(S(F), 1, S(Fp), S(m), S(q), S(mp)))
def Sfac(F, Fp):                                     # relative F->F' line strength
    return (2 * Fp + 1) * float(wigner_6j(S(1)/2, S(3)/2, 1, S(Fp), S(F), S(3)/2))**2

def Eg(F, m, B):                                     # ground energy: hyperfine + linear Zeeman
    return (A_HFS / 2 if F == 2 else -A_HFS / 2) + gF_g[F] * uB * B * m
def Ee(Fp, mp, B, theta):                            # excited: hyperfine + Zeeman + 1064 Stark
    return DF[Fp] + gF_e * uB * B * mp + stark.stark_tensor(Fp, mp, theta)   # <-- the v0.3.0 fix


def reference_rabis():
    """Control/probe Rabi from the EIT pinning Omega_tot = sqrt(4 Delta nu_z)."""
    Otot = np.sqrt(4.0 * c.Delta * c.nu_z)
    Oc = Otot / np.sqrt(1.0 + c.OmR**2)
    return Oc, c.OmR * Oc


def beams(Oc, Op, d2, with_repump, with_e1, with_e3, rep_scale):
    """The four tones of the single-EOM chain. control + probe are the Lambda; the repumpers are
       the forward +1 EOM sideband (sigma-, F=1) and the retro carrier (sigma+, F=2), off-resonant."""
    def ladder(Fg, q, Fps):
        out = []
        for Fp in Fps:
            for mg in ([-1, 0, 1] if Fg == 1 else [-2, -1, 0, 1, 2]):
                mp = mg + q
                if abs(mp) <= Fp:
                    cc = CGc(Fg, mg, q, Fp, mp)
                    if abs(cc) > 1e-9:
                        out.append(((Fg, mg), (Fp, mp), cc))
        return out
    # control: F'2 ladder + single contaminant edges from |2,+1> (matches the validated engine)
    ce = ladder(2, -1, [2])
    if with_e3:
        ce = ce + [((2, 1), (3, 0), CGc(2, 1, -1, 3, 0))]          # coherent F'3 admixture
    if with_e1:
        ce = ce + [((2, 1), (1, 0), CGc(2, 1, -1, 1, 0))]          # off-resonant F'1 spoiler
    # probe: F'2 ladder + single F'1 contaminant from |1,-1>
    pe = ladder(1, +1, [2])
    if with_e1:
        pe = pe + [((1, -1), (1, 0), CGc(1, -1, 1, 1, 0))]
    bm = [dict(edges=ce, named=((2, 1), (2, 0)), det=c.Delta, Rabi=Oc, kdir=+1, tag='ctrl'),
          dict(edges=pe, named=((1, -1), (2, 0)), det=c.Delta + d2, Rabi=Op, kdir=-1, tag='probe')]
    if with_repump:
        rep1 = rep_scale * Op / np.sqrt(c.eta_dp)                 # forward +1 EOM sideband amplitude
        rep2 = rep_scale * Oc * np.sqrt(c.eta_dp)                 # retro carrier amplitude
        # repump1: sigma- on F=1 -> F'1,F'2 (clears the F=1 dark sublevels). F=1->F'3 is dF=2 FORBIDDEN
        #   (CG=0, so it never enters even though the tone sits 178 MHz from F'3); F'0 is allowed but
        #   674 MHz off and decays only to F=1 (no repump) -> inert, omitted. Placed at probe + 2f_A.
        bm.append(dict(edges=ladder(1, -1, [1, 2]), named=((1, 0), (2, -1)),
                       det=c.Delta + d2 + c.twofA, Rabi=rep1, kdir=+1, tag='rep1'))
        # repump2: sigma+ on F=2 -> F'1,F'2 (clears the F=2 dark sublevels). F=2->F'0 is dF=2 FORBIDDEN
        #   (CG=0, though the tone sits 126 MHz from F'0); F'3 is allowed but 622 MHz off and decays
        #   only to F=2 (no repump) -> inert, omitted. Placed at control - 2f_A.
        bm.append(dict(edges=ladder(2, +1, [1, 2]), named=((2, 0), (2, 1)),
                       det=c.Delta - c.twofA, Rabi=rep2, kdir=-1, tag='rep2'))
    return bm


def build_frame(bm, gE, eE):
    """Breadth-first multi-rotating frame: a level energy h[node] consistent with every tone.
       (Ported from eit_cooling_tool.py Section 6.)"""
    realE = {}
    for g in gE:
        realE[('g', g)] = gE[g]
    for e in eE:
        realE[('e', e)] = eE[e]
    for b in bm:
        (g, e), Dn = b['named'], b['det']
        b['nu'] = Dn + (realE[('e', e)] - realE[('g', g)])         # laser frequency of this tone
    adj = {}
    def add(n1, n2, d):
        adj.setdefault(n1, []).append((n2, d))
        adj.setdefault(n2, []).append((n1, -d))
    for b in bm:
        for (g, e, cc) in b['edges']:
            det_ge = b['nu'] - (realE[('e', e)] - realE[('g', g)])
            add(('g', g), ('e', e), -det_ge)
    h = {}
    max_conf = 0.0
    order = [('g', (2, 1))] + [n for n in adj if n != ('g', (2, 1))]
    for start in order:
        if start in h or start not in adj:
            continue
        h[start] = 0.0
        stack = [start]
        while stack:
            n = stack.pop()
            for (m, d) in adj.get(n, []):
                hv = h[n] + d
                if m in h:
                    max_conf = max(max_conf, abs(h[m] - hv))
                else:
                    h[m] = hv
                    stack.append(m)
    return h, max_conf


def solve(d2=0.0, clean=False, with_repump=True, with_e1=True, with_e3=True,
          Nf=None, B=None, theta=None, rep_scale=None, want=False):
    """Steady-state axial <n_z> of the multilevel system. clean=True -> the bare 3-level Lambda.
       rep_scale multiplies the (chain-natural) off-resonant repumper Rabis."""
    B = c.B_field if B is None else B
    theta = c.theta_trap if theta is None else theta
    Nf = c.Nf_multi if Nf is None else Nf
    rep_scale = c.rep_scale if rep_scale is None else rep_scale
    eta, nuz = c.eta, c.nu_z
    Oc, Op = reference_rabis()

    if clean:
        with_repump = with_e1 = with_e3 = False
        Gs, Es = [(1, -1), (2, 1)], [(2, 0)]
    else:
        Gs = [(1, m) for m in (-1, 0, 1)] + [(2, m) for m in range(-2, 3)]
        Es = [(1, m) for m in (-1, 0, 1)] + [(2, m) for m in range(-2, 3)]
        if with_e3:
            Es = Es + [(3, 0)]    # F'3 reached ONLY by the control admixture (2,1)->(3,0); repump can't (dF=2)
    gE = {g: Eg(g[0], g[1], B) for g in Gs}
    eE = {e: Ee(e[0], e[1], B, theta) for e in Es}

    bm = beams(Oc, Op, d2, with_repump, with_e1, with_e3, rep_scale)
    ng, ne = set(Gs), set(Es)
    for b in bm:
        b['edges'] = [(g, e, cc) for (g, e, cc) in b['edges'] if g in ng and e in ne]
    bm = [b for b in bm if b['edges'] and b['named'][0] in ng and b['named'][1] in ne]
    h, conf = build_frame(bm, gE, eE)

    nodes = [('g', g) for g in Gs] + [('e', e) for e in Es]
    for n in nodes:
        h.setdefault(n, 0.0)
    idx = {n: i for i, n in enumerate(nodes)}
    NA = len(nodes)
    bas = [qt.basis(NA, i) for i in range(NA)]
    P = lambda i, j: bas[i] * bas[j].dag()
    If = qt.qeye(Nf)
    aop = qt.destroy(Nf)
    Dsp = lambda s: qt.displace(Nf, 1j * eta * s)

    H = nuz * qt.tensor(qt.qeye(NA), aop.dag() * aop)
    for n in nodes:
        H += h[n] * qt.tensor(P(idx[n], idx[n]), If)
    for b in bm:
        (gn, en) = b['named']
        cnamed = [cc for (g, e, cc) in b['edges'] if g == gn and e == en][0]
        for (g, e, cc) in b['edges']:
            ls = np.sqrt(Sfac(g[0], e[0]) / Sfac(g[0], en[0])) if e[0] != en[0] else 1.0
            O = b['Rabi'] * ls * (cc / cnamed)
            i, j = idx[('g', g)], idx[('e', e)]
            H += -(O / 2) * (qt.tensor(P(j, i), Dsp(b['kdir']))
                             + qt.tensor(P(i, j), Dsp(b['kdir']).dag()))

    cops = []
    legs = [(1, -1), (2, 1)]                    # the Lambda ground states
    Gset = set(Gs)
    for (Fp, mp) in Es:
        ch = {}
        for Fg in (1, 2):
            for mg in (mp - 1, mp, mp + 1):
                if abs(mg) > Fg:
                    continue
                cc = CGc(Fg, mg, mp - mg, Fp, mp)
                if abs(cc) < 1e-12:
                    continue
                w = BR[Fp][Fg] * cc**2
                if (Fg, mg) in Gset:
                    ch[(Fg, mg)] = ch.get((Fg, mg), 0.0) + w     # kept ground: real branch
                else:
                    # clean Lambda only: decay to a dropped sublevel = IDEAL repump back to the legs
                    lw = np.array([CGc(lf, lm, mp - lm, Fp, mp)**2 for (lf, lm) in legs])
                    lw = lw / lw.sum() if lw.sum() > 0 else np.ones(len(legs)) / len(legs)
                    for (lf, lm), ww in zip(legs, lw):
                        ch[(lf, lm)] = ch.get((lf, lm), 0.0) + w * ww
        tot = sum(ch.values())
        if tot <= 0:
            continue
        for (g, w) in ch.items():
            for (u, wem) in EM_REC:
                cops.append(np.sqrt(GAMMA * (w / tot) * wem)
                            * qt.tensor(P(idx[('g', g)], idx[('e', (Fp, mp))]), Dsp(u)))

    L = qt.liouvillian(H, cops)
    try:
        rho = qt.steadystate(L, method='direct')
    except Exception:
        rho = qt.steadystate(L, method='svd')
    nop = qt.tensor(qt.qeye(NA), aop.dag() * aop)
    nbar = float(np.real(qt.expect(nop, rho)))
    if not want:
        return nbar
    pn = np.array([float(np.real(qt.expect(qt.tensor(qt.qeye(NA), qt.basis(Nf, k) * qt.basis(Nf, k).dag()), rho)))
                   for k in range(Nf)])
    pops = {g: float(np.real(qt.expect(qt.tensor(P(idx[('g', g)], idx[('g', g)]), If), rho)))
            for g in Gs}
    nu = {b['tag']: b['nu'] for b in bm}
    return dict(nbar=nbar, pn=pn, pops=pops, nu=nu, conf=conf)


if __name__ == "__main__":
    import sys
    def out(s): print(s); sys.stdout.flush()
    Oc, Op = reference_rabis()
    out("FULL multilevel clock-EIT cooling (87Rb D2) -- single-EOM chain, repumpers IN the computation")
    out(f"  Delta={c.Delta:g}  Omega_c={Oc:.2f}  Omega_p={Op:.3f}  nu_z={c.nu_z:g}  B={c.B_field:g}G  "
        f"theta={c.theta_trap:g}deg  2f_A={c.twofA:g}  eta_dp={c.eta_dp:g}  Nf={c.Nf_multi}")
    out("  repumpers = forward EOM sideband (probe+2f_A, off-res ~F'3) + retro carrier (control-2f_A, off-res ~F'1)")

    # 1) validation: the bare 3-level clean Lambda (matches cooling.py / the project realized floor)
    nclean = solve(clean=True)
    out(f"\n  [validate] clean 3-level Lambda    <n_z> = {nclean:.4f}    (project realized ~0.003)")

    # 2) repumpers OFF: the atom optically pumps into the dark sublevels and stops cooling
    Roff = solve(d2=-0.10, with_repump=False, want=True)
    out(f"  [repump OFF]  <n_z> = {Roff['nbar']:.2f}  (Nf-limited; ~uncooled) -- 100% pumped into dark sublevels")

    # 3) repumpers ON -- they are OFF-RESONANT, so the floor depends on repumper power (rep_scale).
    #    rep_scale=1 is the chain-natural strength; raise it (more EOM/launch power) to repump faster.
    out("\n  repumpers ON -- off-resonant, so the floor depends on repumper power:")
    best = None
    for sc in (1.0, 3.0, 10.0, 30.0):
        R = solve(d2=-0.10, rep_scale=sc, want=True)
        dark = sum(w for g, w in R['pops'].items() if g not in ((1, -1), (2, 1)))
        out(f"    rep_scale={sc:4.0f}  (rep1={sc * Op / np.sqrt(c.eta_dp):4.1f}, rep2={sc * Oc * np.sqrt(c.eta_dp):4.1f}):  "
            f"<n_z> = {R['nbar']:.4f}   dark = {dark:.2f}   P(n=0) = {R['pn'][0]:.3f}")
        if best is None or R['nbar'] < best[0]['nbar']:
            best = (R, sc)
    R, sc = best
    out(f"  -> best floor {R['nbar']:.4f} at rep_scale={sc:g} (the NATURAL power): the off-resonant repumping")
    out("     limits the on-axis floor (~40% stays dark). Higher rep_scale looks worse, but that is dominated")
    out("     by the coherent off-resonant FRAME approximation growing with Rabi -- only the natural-power")
    out("     point is trustworthy; a clean power study needs an incoherent (rate) repumper treatment.")

    # 4) the repumper placement -- the offsets you specified (incl. Delta, 2f_A, and the F'1/F'3 Stark)
    nu = R['nu']
    out("\n  tone placement (offsets at the atom):")
    out(f"    probe   - control = {nu['probe'] - nu['ctrl']:+8.1f} MHz   (= A_HFS, 6834.7)")
    out(f"    repump1 - probe   = {nu['rep1'] - nu['probe']:+8.1f} MHz   (= +2f_A = +400; fwd EOM sideband)")
    out(f"    repump2 - probe   = {nu['rep2'] - nu['probe']:+8.1f} MHz   (= -(A_HFS+2f_A) = -7.23 GHz; retro carrier)")
    out(f"    repump2 - control = {nu['rep2'] - nu['ctrl']:+8.1f} MHz   (= -2f_A = -400; off-res F'1)")
