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
  - the probe leg carries the photon recoil exp(i eta (a+a^dag)) ~ 1 + i eta (a+a^dag) [LD coupling]
  - spontaneous emission |e>-> |g1>,|g2> at Gamma/2 each                               [dissipation]

Run:  python cooling.py
Needs: numpy, qutip.
"""
import numpy as np
from qutip import basis, tensor, qeye, destroy, steadystate, expect
import config as c


def nbar(delta2):
    """Steady-state axial <n_z> for a two-photon detuning delta2 (probe - control)."""
    N = c.N_fock
    a = tensor(qeye(3), destroy(N))
    n = a.dag() * a
    g1, g2, e = basis(3, 0), basis(3, 1), basis(3, 2)

    P_g1 = tensor(g1 * g1.dag(), qeye(N))           # |g1><g1|
    P_g2 = tensor(g2 * g2.dag(), qeye(N))           # |g2><g2|
    seg1 = tensor(e * g1.dag(), qeye(N))            # |e><g1|  (probe leg)
    seg2 = tensor(e * g2.dag(), qeye(N))            # |e><g2|  (control leg)

    # probe leg dressed by the Lamb-Dicke recoil (first order in eta is plenty at eta=0.094)
    LD = tensor(qeye(3), qeye(N)) + 1j * c.eta * (a + a.dag())
    probe = seg1 * LD
  # control leg only diagonal in the harmonic oscillator space (FM, 30/6/2026) 

    # rotating frame: ground states carry the laser detunings; |e> at 0.
    H = ((c.Delta + delta2) * P_g1 + c.Delta * P_g2 + c.nu_z * n
         + (c.Omega_c / 2.0) * (seg2 + seg2.dag())
         + (c.Omega_p / 2.0) * (probe + probe.dag()))

    c_ops = [np.sqrt(c.Gamma / 2.0) * tensor(g1 * e.dag(), qeye(N)),   # |e>->|g1>
             np.sqrt(c.Gamma / 2.0) * tensor(g2 * e.dag(), qeye(N))]   # |e>->|g2>

    return float(expect(n, steadystate(H, c_ops)))


if __name__ == "__main__":
    # delta2 is servoed in the experiment; here we scan it to FIND the floor (the servo point).
    scan = np.linspace(-0.5, 0.5, 41)
    vals = [nbar(d) for d in scan]
    i = int(np.argmin(vals))

    analytic = (c.Gamma / (4.0 * c.Delta)) ** 2
    print("single-atom, on-axis clock-EIT cooling -- Layer 1 (clean 3-level Lambda)")
    print(f"  Delta={c.Delta:g}  Omega_c={c.Omega_c:.3f}  Omega_p={c.Omega_p:.3f}  "
          f"nu_z={c.nu_z:g}  eta={c.eta:g}  (2pi*MHz)")
    print(f"  numeric floor   <n_z> = {vals[i]:.4f}   at delta2 = {scan[i]:+.3f} (servo point)")
    print(f"  analytic floor  (Gamma/4Delta)^2 = {analytic:.4f}")
    print(f"  ground-state population P(n=0) ~ {1.0/(1.0+vals[i]):.3f}")
