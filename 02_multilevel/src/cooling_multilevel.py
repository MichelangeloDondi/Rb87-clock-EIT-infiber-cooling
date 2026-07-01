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
SLOWLY; raise config.repump_scale (more EOM/launch power) to repump faster.

The atomic core is a standard multilevel Lindblad treatment of the 87Rb D2 manifold: Clebsch-Gordan
dipole ladders, a multi-rotating frame, m-resolved decay branching, and the 3-point photon recoil.
The per-(F',m') 1064 tensor Stark from stark.py is fed into excited_energy(), so the off-resonant repumper
detunings are referenced to the real Stark-shifted F'1/F'3 lines.

METHOD: the off-resonant repumpers are modeled as INCOHERENT scattering rates -- the virtual F'1/F'2
excited are adiabatically eliminated (R = Gamma (O/2)^2 / (d^2 + (Gamma/2)^2), with m-resolved decay and
the two-photon absorb+emit recoil). Unlike a coherent beam this needs no rotating-frame closure, so it is
frame-consistent at any power; control + probe stay coherent (their frame closes exactly, frame_conflict=0) and the
repumper excited never enter the Hilbert space (so the solve is smaller/faster).
  SCOPE: the rate above is the LOW-SATURATION limit -- it omits saturation (real scatter caps near Gamma/2)
  and the a.c.-Stark shift ~ Omega^2/d that grows with power. It is reliable only for repump_scale <~ 1
  (chain-natural power); the repump_scale sweep BELOW is therefore qualitative above ~natural power, and the
  high-power blow-up is the rate model breaking, not physics. Trust the repump_scale~1 point, not the trend.

Run:  python cooling_multilevel.py        (needs numpy, qutip, sympy)
"""
import numpy as np
import qutip as qt
from sympy import S
from sympy.physics.wigner import clebsch_gordan, wigner_6j
import config as c
import stark

# ---- 87Rb D2 constants ----------------------------------------------------------------------
A_HFS = 6834.682610                                  # ground hyperfine splitting (2pi MHz)
GAMMA = c.Gamma                                      # 5P3/2 natural linewidth
excited_hf_centroids = {0: -302.07, 1: -229.85, 2: -72.91, 3: 193.74}  # 5P3/2 hyperfine centroids (2pi MHz)
excited_hf_spacings = {Fp: excited_hf_centroids[Fp] - excited_hf_centroids[2] for Fp in excited_hf_centroids}           # spacings from F'=2 {1:-156.94, 3:+266.65,...}
decay_branching = {0: {1: 1.0, 2: 0.0}, 1: {1: 5/6, 2: 1/6},     # F'->F hyperfine decay branching
      2: {1: 1/2, 2: 1/2}, 3: {1: 0.0, 2: 1.0}}
emission_recoil = [(-1, 1/6), (0, 2/3), (1, 1/6)]            # dipole emission recoil (delta-m, weight)
muB = 1.399624                                        # Bohr magneton / frame_energy (MHz/G)
gF_ground = {1: -0.5, 2: +0.5}                            # ground Lande gF (clock pair: both gF*m=+1/2)
gF_excited = 2.0 / 3.0                                     # 5P3/2 Lande gF (equal for F'=1,2,3)

clebsch = lambda F, m, q, Fp, mp: float(clebsch_gordan(S(F), 1, S(Fp), S(m), S(q), S(mp)))
def line_strength(F, Fp):                                     # relative F->F' line strength (= line_amplitude**2)
    return (2 * Fp + 1) * float(wigner_6j(S(1)/2, S(3)/2, 1, S(Fp), S(F), S(3)/2))**2
def line_amplitude(F, Fp):                                     # SIGNED relative line amplitude (carries the 6j sign)
    return np.sqrt(2 * Fp + 1) * float(wigner_6j(S(1)/2, S(3)/2, 1, S(Fp), S(F), S(3)/2))

def ground_energy(F, m, B):                                     # ground energy: hyperfine + linear Zeeman
    return (A_HFS / 2 if F == 2 else -A_HFS / 2) + gF_ground[F] * muB * B * m
def excited_energy(Fp, mp, B, theta):                            # excited: hyperfine + Zeeman + 1064 Stark
    return excited_hf_spacings[Fp] + gF_excited * muB * B * mp + stark.stark_tensor(Fp, mp, theta)   # + per-(F',m') 1064 tensor Stark


def reference_rabis():
    """Control/probe Rabi from the EIT pinning Omega_tot = sqrt(4 Delta nu_z)."""
    Otot = np.sqrt(4.0 * c.Delta * c.nu_z)
    Oc = Otot / np.sqrt(1.0 + c.probe_control_ratio**2)
    return Oc, c.probe_control_ratio * Oc


def dipole_edges(Fg, q, Fps):
    """All CG-allowed dipole edges (Fg,mg)->(Fp,mp=mg+q) for the listed F'."""
    out = []
    for Fp in Fps:
        for mg in ([-1, 0, 1] if Fg == 1 else [-2, -1, 0, 1, 2]):
            mp = mg + q
            if abs(mp) <= Fp:
                cc = clebsch(Fg, mg, q, Fp, mp)
                if abs(cc) > 1e-9:
                    out.append(((Fg, mg), (Fp, mp), cc))
    return out


def beams(Oc, Op, d2, with_e1, with_e3):
    """The two COHERENT tones: control sigma- (F2->F'2) and probe sigma+ (F1->F'2), the Lambda.
       The off-resonant repumpers are NOT here -- they are handled as incoherent rates (repump_collapse_ops),
       so they never enter the rotating frame (which control+probe close exactly, frame_conflict=0)."""
    Dl = c.Delta
    ce = dipole_edges(2, -1, [2])
    if with_e3:
        ce = ce + [((2, 1), (3, 0), clebsch(2, 1, -1, 3, 0))]          # coherent F'3 admixture
    if with_e1:
        ce = ce + [((2, 1), (1, 0), clebsch(2, 1, -1, 1, 0))]          # F'1 spoiler (cooling-beam leak)
    pe = dipole_edges(1, +1, [2])
    if with_e1:
        pe = pe + [((1, -1), (1, 0), clebsch(1, -1, 1, 1, 0))]
    return [dict(edges=ce, named=((2, 1), (2, 0)), det=Dl, Rabi=Oc, kdir=+1, tag='ctrl'),
            dict(edges=pe, named=((1, -1), (2, 0)), det=Dl + d2, Rabi=Op, kdir=-1, tag='probe')]


def repump_beams(Oc, Op, d2, repump_scale, shift=-1, tag_shift=None):
    """The off-resonant repumpers (INCOHERENT). `shift` is the tag-AOM order:
         shift=-1 : retro DOWN-shifted by 2f_A  -> EOM at f_mod = A_HFS + 2f_A  (current config)
         shift=+1 : retro UP-shifted   by 2f_A  -> EOM at f_mod = A_HFS - 2f_A  (the 'other order')
    The forward sideband (rep1, F=1) and the retro carrier (rep2, F=2) move OPPOSITELY with `shift`:
         rep1 det (from F=1->F'2) = Delta + d2 - shift*2f_A
         rep2 det (from F=2->F'2) = Delta + shift*2f_A
    Ladders cover the nearby lines: rep1 -> F'0,F'1,F'2 (F=1->F'3 is dF=2 forbidden); rep2 -> F'1,F'2,F'3
    (F=2->F'0 is dF=2 forbidden). NB F'0 (from F=1) and F'3 (from F=2) decay back to the SAME hyperfine --
    intra-F m-shuffle + recoil heating, NO inter-F repump -- so a tone landing near them wastes power on heat."""
    Dl = c.Delta
    tag_shift = c.tag_shift if tag_shift is None else tag_shift
    rep1 = repump_scale * Op / np.sqrt(c.retro_efficiency)                      # forward +1 EOM sideband amplitude
    rep2 = repump_scale * Oc * np.sqrt(c.retro_efficiency)                      # retro carrier amplitude
    return [
        dict(edges=dipole_edges(1, -1, [0, 1, 2]), named=((1, 0), (2, -1)),
             det=Dl + d2 - shift * tag_shift, Rabi=rep1, kdir=+1, tag='rep1'),
        dict(edges=dipole_edges(2, +1, [1, 2, 3]), named=((2, 0), (2, 1)),
             det=Dl + shift * tag_shift, Rabi=rep2, kdir=-1, tag='rep2'),
    ]


def decay_channels(Fp, mp):
    """m-resolved spontaneous decay of a (virtual) |Fp,mp> to the ground sublevels: {(Fg,mg): weight}."""
    decay_targets = {}
    for Fg in (1, 2):
        for mg in (mp - 1, mp, mp + 1):
            if abs(mg) <= Fg:
                cc = clebsch(Fg, mg, mp - mg, Fp, mp)
                if abs(cc) > 1e-12:
                    decay_targets[(Fg, mg)] = decay_targets.get((Fg, mg), 0.0) + decay_branching[Fp][Fg] * cc**2
    return decay_targets


def repump_collapse_ops(repump_tones, ground_energies, idx, P, fock_identity, displacement, ground_set, B, theta):
    """Off-resonant repumpers as INCOHERENT scattering -- frame-free, but a LOW-SATURATION rate (valid repump_scale<~1).
       Each edge (gsrc)->(Fp,mp) at detuning d, edge-Rabi O scatters at R = Gamma (O/2)^2/(d^2+(Gamma/2)^2);
       the virtual excited |Fp,mp> is adiabatically eliminated and decays m-resolved, carrying the
       two-photon (absorb kdir + emit u) recoil. Also returns the off-resonant a.c.-Stark shift per ground."""
    cops, ac_stark_shift = [], {}
    shifted_excited_energy = lambda Fp, mp: excited_energy(Fp, mp, B, theta)                       # virtual excited energy (Stark-shifted)
    min_detuning = np.inf                                                  # closest any repumper tone comes to a line
    for b in repump_tones:
        (gn, en) = b['named']
        named_coupling = [cc for (g, e, cc) in b['edges'] if g == gn and e == en][0]
        nu = b['det'] + (shifted_excited_energy(*en) - ground_energies[gn])                        # tone frequency (centered-energy frame)
        b['nu'] = nu
        for (g, e, cc) in b['edges']:
            if g not in ground_set:
                continue
            # incoherent rate: sqrt(line_strength) magnitude is fine here -- R and ac_stark_shift below both go as O**2, so the
            # reduced-element SIGN is irrelevant (unlike the coherent H, which uses the signed line_amplitude ratio).
            line_scale = np.sqrt(line_strength(g[0], e[0]) / line_strength(g[0], en[0])) if e[0] != en[0] else 1.0
            O = b['Rabi'] * line_scale * (cc / named_coupling)                    # this edge's Rabi
            d = nu - (shifted_excited_energy(*e) - ground_energies[g])                             # this edge's detuning
            min_detuning = min(min_detuning, abs(d))
            R = GAMMA * (O / 2)**2 / (d**2 + (GAMMA / 2)**2)      # off-resonant scattering rate
            ac_stark_shift[g] = ac_stark_shift.get(g, 0.0) + (O / 2)**2 / d   # GROUND-only shift (eliminated excited's is dropped); physical as the leg DIFFERENTIAL at repump_scale<~1
            decay_branch = decay_channels(e[0], e[1])
            branch_total = sum(decay_branch.values())
            if branch_total <= 0:
                continue
            for (decay_ground, w) in decay_branch.items():
                if decay_ground not in ground_set:
                    continue
                for (u, emit_weight) in emission_recoil:
                    cops.append(np.sqrt(R * (w / branch_total) * emit_weight)
                                * qt.tensor(P(idx[('g', decay_ground)], idx[('g', g)]), displacement(u) * displacement(b['kdir'])))
    # Guard: the incoherent low-saturation rate (R, ac_stark_shift) is only valid far off resonance. If a config
    # (e.g. a small 2f_A in explore_configs) lands a comb tone within a few linewidths of a line, R approaches
    # the saturation limit and the ac_stark_shift ~ Omega^2/d diverges -- a near-resonant tone must be treated
    # coherently (in the Hilbert space), not as a rate. We WARN here (the solve still returns a number, but it
    # is NOT trustworthy for such a config); a caller that sweeps configs should flag or skip the warned rows.
    if min_detuning < 3.0 * c.Gamma:
        import warnings
        warnings.warn("a repumper tone sits %.0f MHz (< 3*Gamma) from a line: the incoherent rate model is "
                      "invalid here (treat it coherently); <n_z> is not trustworthy for this config." % min_detuning,
                      RuntimeWarning)
    return cops, ac_stark_shift


def build_rotating_frame(tones, ground_energies, excited_energies):
    """Breadth-first multi-rotating frame: a level energy frame_energy[node] consistent with every tone."""
    real_energies = {}
    for g in ground_energies:
        real_energies[('g', g)] = ground_energies[g]
    for e in excited_energies:
        real_energies[('e', e)] = excited_energies[e]
    for b in tones:
        (g, e), tone_detuning = b['named'], b['det']
        b['nu'] = tone_detuning + (real_energies[('e', e)] - real_energies[('g', g)])         # laser frequency of this tone
    adjacency = {}
    def add(n1, n2, d):
        adjacency.setdefault(n1, []).append((n2, d))
        adjacency.setdefault(n2, []).append((n1, -d))
    for b in tones:
        for (g, e, cc) in b['edges']:
            detuning_ge = b['nu'] - (real_energies[('e', e)] - real_energies[('g', g)])
            add(('g', g), ('e', e), -detuning_ge)
    frame_energy = {}
    max_frame_conflict = 0.0
    order = [('g', (2, 1))] + [n for n in adjacency if n != ('g', (2, 1))]
    for start in order:
        if start in frame_energy or start not in adjacency:
            continue
        frame_energy[start] = 0.0
        stack = [start]
        while stack:
            n = stack.pop()
            for (m, d) in adjacency.get(n, []):
                h_candidate = frame_energy[n] + d
                if m in frame_energy:
                    max_frame_conflict = max(max_frame_conflict, abs(frame_energy[m] - h_candidate))
                else:
                    frame_energy[m] = h_candidate
                    stack.append(m)
    return frame_energy, max_frame_conflict


def solve(d2=0.0, clean=False, with_repump=True, with_e1=True, with_e3=True,
          N_fock=None, B=None, theta=None, repump_scale=None, shift=-1, tag_shift=None, want=False):
    """Steady-state axial <n_z> of the multilevel system. clean=True -> the bare 3-level Lambda.
       repump_scale multiplies the (chain-natural) off-resonant repumper Rabis; shift/tag_shift pick the
       EOM/AOM configuration (see repump_beams)."""
    B = c.B_field if B is None else B
    theta = c.theta_trap if theta is None else theta
    N_fock = c.N_fock if N_fock is None else N_fock
    repump_scale = c.repump_scale if repump_scale is None else repump_scale
    eta, nu_z = c.eta, c.nu_z
    Oc, Op = reference_rabis()

    if clean:
        with_repump = with_e1 = with_e3 = False
        ground_states, excited_states = [(1, -1), (2, 1)], [(2, 0)]
    else:
        ground_states = [(1, m) for m in (-1, 0, 1)] + [(2, m) for m in range(-2, 3)]
        excited_states = [(2, m) for m in range(-2, 3)]              # F'2: the cooling Lambda's upper manifold
        if with_e1:
            excited_states = excited_states + [(1, 0)]                           # F'1 spoiler (cooling-beam off-res leak)
        if with_e3:
            excited_states = excited_states + [(3, 0)]                           # F'3 admixture (control)
        # the repumpers' F'1/F'2 are VIRTUAL (incoherent, eliminated) -- not in the Hilbert space
    ground_energies = {g: ground_energy(g[0], g[1], B) for g in ground_states}
    excited_energies = {e: excited_energy(e[0], e[1], B, theta) for e in excited_states}

    tones = beams(Oc, Op, d2, with_e1, with_e3)
    kept_grounds, kept_excited = set(ground_states), set(excited_states)
    for b in tones:
        b['edges'] = [(g, e, cc) for (g, e, cc) in b['edges'] if g in kept_grounds and e in kept_excited]
    tones = [b for b in tones if b['edges'] and b['named'][0] in kept_grounds and b['named'][1] in kept_excited]
    frame_energy, frame_conflict = build_rotating_frame(tones, ground_energies, excited_energies)
    # Guard: a stationary rotating frame exists only if every loop in the tone graph closes (frame_conflict=0).
    # If a future config (e.g. a microwave tone coupling the two grounds) closes a loop, frame_conflict != 0 and
    # the steady state below would be silently unphysical -- warn loudly rather than return a fake number.
    if frame_conflict > 1e-6:
        import warnings
        warnings.warn("rotating-frame conflict frame_conflict=%.3g (2pi MHz): no time-independent frame exists for this "
                      "tone set; the steady-state <n_z> is NOT trustworthy." % frame_conflict, RuntimeWarning)

    nodes = [('g', g) for g in ground_states] + [('e', e) for e in excited_states]
    for n in nodes:
        frame_energy.setdefault(n, 0.0)
    idx = {n: i for i, n in enumerate(nodes)}
    n_states = len(nodes)
    basis_kets = [qt.basis(n_states, i) for i in range(n_states)]
    P = lambda i, j: basis_kets[i] * basis_kets[j].dag()
    fock_identity = qt.qeye(N_fock)
    annihilation = qt.destroy(N_fock)
    displacement = lambda s: qt.displace(N_fock, 1j * eta * s)

    H = nu_z * qt.tensor(qt.qeye(n_states), annihilation.dag() * annihilation)
    for n in nodes:
        H += frame_energy[n] * qt.tensor(P(idx[n], idx[n]), fock_identity)
    for b in tones:
        (gn, en) = b['named']
        named_coupling = [cc for (g, e, cc) in b['edges'] if g == gn and e == en][0]
        for (g, e, cc) in b['edges']:
            # line_scale rescales a cross-manifold (e[0]!=en[0]) leak by the relative line amplitude. The SIGNED ratio
            # line_amplitude/line_amplitude carries the 6j sign of the reduced matrix element -- needed because this is a COHERENT
            # term whose relative phase matters for the with_e1/e3 admixture interference. The primary dark
            # state lives inside F'=2 and uses the signed CG ratio cc/named_coupling, which is exact regardless.
            line_scale = (line_amplitude(g[0], e[0]) / line_amplitude(g[0], en[0])) if e[0] != en[0] else 1.0
            O = b['Rabi'] * line_scale * (cc / named_coupling)
            i, j = idx[('g', g)], idx[('e', e)]
            H += -(O / 2) * (qt.tensor(P(j, i), displacement(b['kdir']))
                             + qt.tensor(P(i, j), displacement(b['kdir']).dag()))

    cops = []
    legs = [(1, -1), (2, 1)]                    # the Lambda ground states
    ground_set = set(ground_states)
    for (Fp, mp) in excited_states:
        decay_targets = {}
        for Fg in (1, 2):
            for mg in (mp - 1, mp, mp + 1):
                if abs(mg) > Fg:
                    continue
                cc = clebsch(Fg, mg, mp - mg, Fp, mp)
                if abs(cc) < 1e-12:
                    continue
                w = decay_branching[Fp][Fg] * cc**2
                if (Fg, mg) in ground_set:
                    decay_targets[(Fg, mg)] = decay_targets.get((Fg, mg), 0.0) + w     # kept ground: real branch
                else:
                    # clean Lambda only: decay to a dropped sublevel = IDEAL repump back to the legs
                    leg_weights = np.array([clebsch(lf, lm, mp - lm, Fp, mp)**2 for (lf, lm) in legs])
                    leg_weights = leg_weights / leg_weights.sum() if leg_weights.sum() > 0 else np.ones(len(legs)) / len(legs)
                    for (lf, lm), ww in zip(legs, leg_weights):
                        decay_targets[(lf, lm)] = decay_targets.get((lf, lm), 0.0) + w * ww
        branch_total = sum(decay_targets.values())
        if branch_total <= 0:
            continue
        for (g, w) in decay_targets.items():
            for (u, emit_weight) in emission_recoil:
                cops.append(np.sqrt(GAMMA * (w / branch_total) * emit_weight)
                            * qt.tensor(P(idx[('g', g)], idx[('e', (Fp, mp))]), displacement(u)))

    # the off-resonant repumpers: INCOHERENT scattering rates (no rotating-frame loop; low-saturation, valid repump_scale<~1)
    repump_freqs = {}
    if with_repump and not clean:
        repump_tones = repump_beams(Oc, Op, d2, repump_scale, shift, tag_shift)
        repump_ops, ac_stark_shift = repump_collapse_ops(repump_tones, ground_energies, idx, P, fock_identity, displacement, ground_set, B, theta)
        cops += repump_ops
        for g, sh in ac_stark_shift.items():
            H += sh * qt.tensor(P(idx[('g', g)], idx[('g', g)]), fock_identity)   # off-resonant a.c.-Stark shift
        repump_freqs = {b['tag']: b['nu'] for b in repump_tones}

    L = qt.liouvillian(H, cops)
    try:
        rho = qt.steadystate(L, method='direct')
    except Exception:
        rho = qt.steadystate(L, method='svd')
    number_op = qt.tensor(qt.qeye(n_states), annihilation.dag() * annihilation)
    nbar = float(np.real(qt.expect(number_op, rho)))
    if not want:
        return nbar
    pn = np.array([float(np.real(qt.expect(qt.tensor(qt.qeye(n_states), qt.basis(N_fock, k) * qt.basis(N_fock, k).dag()), rho)))
                   for k in range(N_fock)])
    pops = {g: float(np.real(qt.expect(qt.tensor(P(idx[('g', g)], idx[('g', g)]), fock_identity), rho)))
            for g in ground_states}
    nu = {**{b['tag']: b['nu'] for b in tones}, **repump_freqs}
    return dict(nbar=nbar, pn=pn, pops=pops, nu=nu, frame_conflict=frame_conflict)


if __name__ == "__main__":
    import sys
    def out(s): print(s); sys.stdout.flush()
    Oc, Op = reference_rabis()
    out("FULL multilevel clock-EIT cooling (87Rb D2) -- single-EOM chain, repumpers IN the computation")
    out(f"  Delta={c.Delta:g}  Omega_c={Oc:.2f}  Omega_p={Op:.3f}  nu_z={c.nu_z:g}  B={c.B_field:g}G  "
        f"theta={c.theta_trap:g}deg  2f_A={c.tag_shift:g}  retro_efficiency={c.retro_efficiency:g}  N_fock={c.N_fock}")
    out("  repumpers = forward EOM sideband + retro carrier, OFF-RESONANT, as INCOHERENT rates (frame-free)")

    # 1) validation: the bare 3-level clean Lambda (the 3-level floor, now with photon recoil)
    nclean = solve(clean=True)
    out(f"\n  [validate] clean 3-level Lambda    <n_z> = {nclean:.4f}    (recoil floor ~0.003)")

    # 2) repumpers OFF: the atom optically pumps into the dark sublevels and stops cooling
    Roff = solve(d2=-0.10, with_repump=False, want=True)
    out(f"  [repump OFF]  <n_z> = {Roff['nbar']:.2f}  (N_fock-limited; ~uncooled) -- 100% pumped into dark sublevels")

    # 2b) The operating point in d2. The full-manifold light shifts (repumper a.c.-Stark + the coherent F'1/F'3
    #     admixtures) move the dark resonance OFF the bare two-photon resonance d2=0. Scan d2 to see where the
    #     servo actually sits -- it is NOT at d2=0, and the floor is a steep function of it.
    out("\n  floor vs two-photon detuning d2 (repump_scale=1; the servo tracks the dark resonance, NOT d2=0):")
    d2scan = [(d, solve(d2=d, repump_scale=1.0)) for d in (0.0, -0.05, -0.10, -0.15, -0.20, -0.25)]
    for d, n in d2scan:
        out(f"    d2 = {d:+.2f}   <n_z> = {n:.4f}")
    D2_OP = min(d2scan, key=lambda t: t[1])[0]
    out(f"  -> operating point d2 = {D2_OP:+.2f} (dark-resonance minimum); at the bare d2=0 the floor is {d2scan[0][1]:.2f}.")

    # 3) repumpers ON -- INCOHERENT off-resonant rates. Low-saturation rate: TRUST repump_scale~1, NOT the
    #    high-power trend (the rate omits saturation + the a.c.-Stark shift, which break above ~natural power).
    out("\n  repumpers ON (incoherent off-resonant rates; low-saturation -> trust repump_scale~1 only) -- floor vs power:")
    headline = None
    for sc in (0.3, 1.0, 3.0, 10.0, 30.0):
        R = solve(d2=D2_OP, repump_scale=sc, want=True)
        dark = sum(w for g, w in R['pops'].items() if g not in ((1, -1), (2, 1)))
        out(f"    repump_scale={sc:5.1f}  (rep1={sc * Op / np.sqrt(c.retro_efficiency):5.1f}, rep2={sc * Oc * np.sqrt(c.retro_efficiency):5.1f}):  "
            f"<n_z> = {R['nbar']:.4f}   dark = {dark:.2f}   P(n=0) = {R['pn'][0]:.3f}")
        if sc == 1.0:
            headline = R
    R = headline
    out(f"  -> repump_scale=1 (chain-natural power) is the physical, TRUSTWORTHY point:  <n_z> = {R['nbar']:.3f}.")
    out("     Lower repump_scale under-repumps; the further fall at higher repump_scale is the low-saturation rate")
    out("     model breaking (it omits saturation + the a.c.-Stark shift), not physics -- see SCOPE.")

    # 4) the repumper placement -- the offsets you specified (incl. Delta, 2f_A, and the F'1/F'3 Stark)
    nu = R['nu']
    out("\n  tone placement (offsets at the atom):")
    out(f"    probe   - control = {nu['probe'] - nu['ctrl']:+8.1f} MHz   (= A_HFS, 6834.7)")
    out(f"    repump1 - probe   = {nu['rep1'] - nu['probe']:+8.1f} MHz   (= +2f_A = +400; fwd EOM sideband)")
    out(f"    repump2 - probe   = {nu['rep2'] - nu['probe']:+8.1f} MHz   (= -(A_HFS+2f_A) = -7.23 GHz; retro carrier)")
    out(f"    repump2 - control = {nu['rep2'] - nu['ctrl']:+8.1f} MHz   (= -2f_A = -400; off-res F'1)")
