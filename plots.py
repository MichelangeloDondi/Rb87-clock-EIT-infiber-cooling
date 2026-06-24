"""
plots.py -- the Layer 1 figures, all generated from the same physics as cooling.py / stark.py.
Run:  python plots.py     (writes the PNGs next to this file)

  eit_spectrum.png    EIT absorption vs probe detuning -- the cooling mechanism   (README section 3)
  cooling_curve.png   <n_z>(t) from a hot start, + the final Fock distribution    (README section 5)
  stark_manifold.png  the 1064 trap + the 5P_3/2 anti-trap shifts                 (README section 1)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qutip import basis, qeye, destroy, tensor, steadystate, expect, liouvillian
import config as c
import stark

BLUE, RED, GREY, GREEN = "#1565c0", "#c0392b", "#7f8c8d", "#2e7d32"


# ---------------------------------------------------------------- 1. EIT spectrum
def eit_spectrum():
    """Excited-state population (~absorption) vs probe two-photon detuning, electronic only."""
    g1, g2, e = basis(3, 0), basis(3, 1), basis(3, 2)
    Pe = e * e.dag()
    d2 = np.linspace(-3 * c.nu_z, 4 * c.nu_z, 600)
    pe = []
    for d in d2:
        H = ((c.Delta + d) * g1 * g1.dag() + c.Delta * g2 * g2.dag()
             + (c.Omega_c / 2) * (e * g2.dag() + g2 * e.dag())
             + (c.Omega_p / 2) * (e * g1.dag() + g1 * e.dag()))
        cops = [np.sqrt(c.Gamma / 2) * g1 * e.dag(), np.sqrt(c.Gamma / 2) * g2 * e.dag()]
        pe.append(expect(Pe, steadystate(H, cops)))
    pe = np.array(pe)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(d2, pe, color=BLUE, lw=2.2)
    for x, col, lab, ha in [(0.0, GREY, "carrier\n(dark — no heating)", "center"),
                            (+c.nu_z, RED, "cooling sideband\n(on the bright peak)", "left"),
                            (-c.nu_z, GREEN, "heating sideband\n(suppressed)", "right")]:
        ax.axvline(x, color=col, ls="--", lw=1.2, alpha=0.8)
        ax.annotate(lab, xy=(x, pe.max() * 0.92), fontsize=8.5, color=col, ha=ha)
    ax.set_xlabel(r"probe two-photon detuning  $\delta_2$  (2$\pi$·MHz)")
    ax.set_ylabel(r"excited population $P_e$  ($\propto$ absorption)")
    ax.set_title("EIT spectrum: dark carrier, bright cooling sideband, suppressed heating\n"
                 r"(bright peak placed at $+\nu_z$ by $\Omega_c^2/4\Delta=\nu_z$)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("eit_spectrum.png", dpi=150)
    print("wrote eit_spectrum.png")


# ---------------------------------------------------------------- 2. cooling curve + final Fock
def cooling_curve(n_init=5.0, N=12):
    """Steady state (fast) + the Liouvillian cooling rate -> the rate-equation cooling curve.
    A full time integration is needlessly slow here (cooling is slow, but the Delta=45 carrier
    dynamics are fast -- a stiff separation); the exponential below IS that solution."""
    a = tensor(qeye(3), destroy(N)); n = a.dag() * a
    g1, g2, e = basis(3, 0), basis(3, 1), basis(3, 2)
    LD = tensor(qeye(3), qeye(N)) + 1j * c.eta * (a + a.dag())
    probe = tensor(e * g1.dag(), qeye(N)) * LD
    H = (c.Delta * tensor(g1 * g1.dag(), qeye(N)) + c.Delta * tensor(g2 * g2.dag(), qeye(N))
         + c.nu_z * n
         + (c.Omega_c / 2) * (tensor(e * g2.dag(), qeye(N)) + tensor(g2 * e.dag(), qeye(N)))
         + (c.Omega_p / 2) * (probe + probe.dag()))
    cops = [np.sqrt(c.Gamma / 2) * tensor(g1 * e.dag(), qeye(N)),
            np.sqrt(c.Gamma / 2) * tensor(g2 * e.dag(), qeye(N))]

    rho_ss = steadystate(H, cops)
    n_ss = float(expect(n, rho_ss))
    pn = np.real(rho_ss.ptrace(1).diag())

    # cooling rate W = slowest non-zero decay mode of the Liouvillian (the relaxation bottleneck)
    decays = -liouvillian(H, cops).eigenenergies().real
    W = float(decays[decays > 1e-6].min())     # (2pi*MHz)
    tau_us = 0.159 / W                          # 0.159 us per 1/(2pi MHz)
    t_us = np.linspace(0, 10 * tau_us, 300)
    nt = n_ss + (n_init - n_ss) * np.exp(-t_us / tau_us)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.6, 4.2))
    axL.plot(t_us, nt, color=BLUE, lw=2.2)
    axL.axhline(n_ss, color=RED, ls="--", lw=1.2)
    axL.annotate(f"floor  $\\bar n_z$ ≈ {n_ss:.4f}", xy=(t_us[-1] * 0.4, n_ss),
                 xytext=(0, 8), textcoords="offset points", color=RED, fontsize=9.5)
    axL.set_xlabel(r"time ($\mu$s)"); axL.set_ylabel(r"$\langle n_z\rangle$")
    axL.set_title(f"cooling from a hot start ($\\bar n_0={n_init:g}$,  $\\tau\\approx${tau_us:.0f} $\\mu$s)")
    axL.set_yscale("log")

    axR.bar(range(8), pn[:8], color=BLUE, width=0.7)
    axR.annotate(f"$P(n{{=}}0)$ ≈ {pn[0]:.3f}", xy=(0, pn[0]), xytext=(1.5, pn[0] * 0.9),
                 color=RED, fontsize=10)
    axR.set_xlabel("Fock state $n$"); axR.set_ylabel("population")
    axR.set_title("steady-state motional distribution")
    for ax in (axL, axR):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("cooling_curve.png", dpi=150)
    print("wrote cooling_curve.png")


# ---------------------------------------------------------------- 3. Stark manifold
def stark_manifold():
    U0 = -stark.shift(c.alpha0_5S)                 # ground trap depth (down)
    sca = stark.shift(c.alpha0_5P32)               # 5P scalar (up)
    lo = stark.shift(c.alpha0_5P32 + c.alpha2_5P32)   # F'=3 stretched
    hi = stark.shift(c.alpha0_5P32 - c.alpha2_5P32)   # most anti-trapped

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    ax.axhline(0, color="k", lw=0.8, alpha=0.4)
    ax.text(0.5, 2, "free-atom levels (no light)", fontsize=8, color="k", alpha=0.5, ha="center")
    # ground
    ax.hlines(-U0, 0.1, 0.9, color=BLUE, lw=3)
    ax.text(0.5, -U0 - 3, f"5S$_{{1/2}}$ ground: trapped, $U_0$ = {U0:.0f} MHz", color=BLUE, ha="center", fontsize=9)
    # excited manifold
    ax.fill_between([1.1, 1.9], lo, hi, color=RED, alpha=0.12)
    ax.hlines(sca, 1.1, 1.9, color=RED, lw=3)
    ax.text(2.0, sca, f"|F'2,0⟩  pure scalar  +{sca:.0f} MHz  (EIT target)", color=RED, va="center", fontsize=9)
    ax.hlines(lo, 1.1, 1.9, color=GREY, lw=1.5, ls="--")
    ax.text(2.0, lo, f"F'=3 stretched  +{lo:.0f}", color=GREY, va="center", fontsize=8)
    ax.hlines(hi, 1.1, 1.9, color=GREY, lw=1.5, ls="--")
    ax.text(2.0, hi, f"most anti-trapped  +{hi:.0f}", color=GREY, va="center", fontsize=8)
    ax.text(1.5, hi + 6, "5P$_{3/2}$: anti-trapped\n(tensor splits F'=3;\nF'=2 tensor = 0)", color=RED, ha="center", fontsize=8)
    ax.set_xlim(0, 4.2); ax.set_ylim(-U0 - 10, hi + 16)
    ax.set_ylabel("light shift (2$\\pi$·MHz)")
    ax.set_title("1064 nm lattice (1 W + 1 W): ground trapped, 5P$_{3/2}$ expelled")
    ax.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("stark_manifold.png", dpi=150)
    print("wrote stark_manifold.png")


if __name__ == "__main__":
    eit_spectrum()
    cooling_curve()
    stark_manifold()
