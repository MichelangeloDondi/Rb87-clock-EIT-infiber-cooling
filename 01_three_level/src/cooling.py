"""
cooling.py -- the LAST-RESORT numerical confirmation of the analytic floor.

The physics and the algebra (README) get us to the EIT-cooling floor analytically:
  n_min  ~  (Gamma / 4 Delta)^2      [the off-resonant-scatter limit]
What has no clean closed form is the *exact* steady state of the driven 3-level Lambda
once it is dressed by the trap oscillator -- so we solve that one master equation
numerically and check it lands on the analytic estimate. That is the only thing this
script does, and it is the only code in Layer 1.

Model (README sections in [brackets]):
  - 3 atomic levels: 0=|g1>=|1,-1> (probe), 1=|g2>=|2,+1> (control), 2=|e>=|F'2,0>   [scheme]
  - (x) one harmonic oscillator (the axial motion), Fock cutoff N_fock                [scheme]
  - both legs detuned by Delta from |e>; two-photon detuning delta2 = (probe-control)  [op. point]
  - recoil on BOTH counter-propagating legs (probe kdir -1, control kdir +1) AND on the
    spontaneous emission (three-point dipole pattern), via the exact displacement operator [recoil]
  - spontaneous emission |e>-> |g1>,|g2> at Gamma/2 each                               [dissipation]

Run:  python cooling.py
Needs: numpy, qutip.
"""
import numpy as np
from qutip import basis, tensor, qeye, destroy, displace, steadystate, expect
import config as c


def nbar(delta2, want_p0=False):
    """Steady-state axial <n_z> for a two-photon detuning delta2 (probe - control).
       Recoil is on BOTH counter-propagating legs and on the spontaneous emission (see the header).
       want_p0=True also returns the true ground-band population P(n=0) from the steady state."""
    N = c.N_fock
    a = tensor(qeye(3), destroy(N))
    n = a.dag() * a
    g1, g2, e = basis(3, 0), basis(3, 1), basis(3, 2)

    P_g1 = tensor(g1 * g1.dag(), qeye(N))           # |g1><g1|
    P_g2 = tensor(g2 * g2.dag(), qeye(N))           # |g2><g2|

    # the two cooling legs counter-propagate along the trap axis -> opposite absorption recoil
    # (probe kdir -1, control kdir +1). Exact displacement, the same treatment as the multilevel solver.
    displacement = lambda s: tensor(qeye(3), displace(N, 1j * c.eta * s))
    probe   = tensor(e * g1.dag(), qeye(N)) * displacement(-1)     # |e><g1| * recoil  (probe leg)
    control = tensor(e * g2.dag(), qeye(N)) * displacement(+1)     # |e><g2| * recoil  (control leg)

    # rotating frame: ground states carry the laser detunings; |e> at 0.
    H = ((c.Delta + delta2) * P_g1 + c.Delta * P_g2 + c.nu_z * n
         + (c.Omega_c / 2.0) * (control + control.dag())
         + (c.Omega_p / 2.0) * (probe + probe.dag()))

    # spontaneous emission |e> -> |g1>,|g2>, each carrying the axial emission recoil: (u, weight), u = the k_z
    # projection of the emitted photon in units of eta (NOT delta-m). The weights give the isotropic axial variance
    # <u^2>=1/3 -- an approximation (true 1/5 pi, 2/5 sigma) that is harmless for the axial floor.
    emission_recoil = [(-1, 1 / 6), (0, 2 / 3), (1, 1 / 6)]
    c_ops = [np.sqrt(c.Gamma / 2.0 * w) * tensor(gk * e.dag(), qeye(N)) * displacement(u)
             for gk in (g1, g2) for (u, w) in emission_recoil]

    rho = steadystate(H, c_ops)
    nz = float(expect(n, rho))
    if not want_p0:
        return nz
    # true P(n=0): the motional ground-band diagonal of the steady state (NOT the thermal 1/(1+n) shortcut)
    p0 = float(np.real(rho.ptrace(1).diag()[0]))
    return nz, p0


if __name__ == "__main__":
    # delta2 is servoed in the experiment; here we scan it to FIND the floor (the servo point).
    scan = np.linspace(-0.5, 0.5, 41)
    results = [nbar(d, want_p0=True) for d in scan]   # (n_z, P(n=0)) at each servo point
    vals = [nz for nz, _ in results]
    i = int(np.argmin(vals))
    p0 = results[i][1]

    analytic = (c.Gamma / (4.0 * c.Delta)) ** 2
    print("single-atom, on-axis clock-EIT cooling -- Layer 1 (clean 3-level Lambda)")
    print(f"  Delta={c.Delta:g}  Omega_c={c.Omega_c:.3f}  Omega_p={c.Omega_p:.3f}  "
          f"nu_z={c.nu_z:g}  eta={c.eta:g}  (2pi*MHz)")
    print(f"  numeric floor  <n_z> = {vals[i]:.4f}  (full recoil: both legs + emission)  at delta2 = {scan[i]:+.3f}")
    print(f"  analytic floor (Gamma/4Delta)^2 = {analytic:.4f}  (recoil-free mechanism limit)")
    print(f"  ground-state population P(n=0) = {p0:.3f}  (from the steady-state motional diagonal)")
