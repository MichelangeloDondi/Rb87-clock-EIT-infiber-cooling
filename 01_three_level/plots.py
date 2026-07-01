"""
plots.py -- the 3-level figures, all generated from the same physics as cooling.py / stark.py.
Run:  python plots.py     (writes the PNGs next to this file)

  eit_spectrum.png    EIT absorption vs probe detuning -- the cooling mechanism   (README "How EIT cools")
  cooling_curve.png   <n_z>(t) from a hot start, + the final Fock distribution    (README "The number")
  stark_manifold.png  the 1064 trap + the 5P_3/2 anti-trap shifts                 (README "The trap")
"""
import os
import numpy as np
import scipy.linalg as sla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qutip import (basis, qeye, destroy, tensor, steadystate, expect, liouvillian,
                   operator_to_vector, thermal_dm)
import config as c
import stark

HERE = os.path.dirname(os.path.abspath(__file__))
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
    fig.savefig(os.path.join(HERE, "images", "eit_spectrum.png"), dpi=150)
    print("wrote eit_spectrum.png")


# ---------------------------------------------------------------- 2. cooling curve + final Fock
def cooling_curve(n_init=3.0, N=16):
    """EXACT <n_z>(t) from a hot start. The master equation is linear, so <n_z>(t) is an exact sum of
    decaying Liouvillian modes; we get it by eigendecomposition (no stiff Delta=45 time-stepping needed,
    and no single-exponential assumption). The relaxation is multi-modal, and the rate that matters is
    set by the modes <n_z> actually couples to -- NOT the slowest Liouvillian mode (which <n_z> barely
    touches). [A direct mesolve gives the same curve; it is just slow here because of the Delta=45 stiffness.]"""
    a = tensor(qeye(3), destroy(N)); n = a.dag() * a
    g1, g2, e = basis(3, 0), basis(3, 1), basis(3, 2)
    recoil_op = tensor(qeye(3), qeye(N)) + 1j * c.eta * (a + a.dag())
    probe = tensor(e * g1.dag(), qeye(N)) * recoil_op
    H = (c.Delta * tensor(g1 * g1.dag(), qeye(N)) + c.Delta * tensor(g2 * g2.dag(), qeye(N))
         + c.nu_z * n
         + (c.Omega_c / 2) * (tensor(e * g2.dag(), qeye(N)) + tensor(g2 * e.dag(), qeye(N)))
         + (c.Omega_p / 2) * (probe + probe.dag()))
    cops = [np.sqrt(c.Gamma / 2) * tensor(g1 * e.dag(), qeye(N)),
            np.sqrt(c.Gamma / 2) * tensor(g2 * e.dag(), qeye(N))]
    L = liouvillian(H, cops)
    rho_ss = steadystate(H, cops)
    n_ss = float(expect(n, rho_ss))
    pn = np.real(rho_ss.ptrace(1).diag())

    # exact modal reconstruction: rho(t) = sum_k exp(w_k t) |r_k><l_k|rho0>  =>  <n>(t) = sum_k amp_k exp(w_k t)
    w, Vr = sla.eig(L.full())
    rv = operator_to_vector(tensor(g1 * g1.dag(), thermal_dm(N, n_init))).full().ravel()
    nv = operator_to_vector(n).full().ravel()
    amp = (nv.conj() @ Vr) * (sla.solve(Vr, rv))
    t_us = np.linspace(0, 700, 400)
    nt = np.clip(np.real([np.sum(amp * np.exp(w * (t / 0.159))) for t in t_us]), 1e-6, None)
    tau_e = t_us[int(np.argmin(np.abs(nt - (n_ss + (nt[0] - n_ss) / np.e))))]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.6, 4.2))
    axL.plot(t_us, nt, color=BLUE, lw=2.2)
    axL.axhline(n_ss, color=RED, ls="--", lw=1.2)
    axL.annotate(f"floor  $\\bar n_z$ ≈ {n_ss:.4f}", xy=(t_us[-1] * 0.45, n_ss),
                 xytext=(0, 8), textcoords="offset points", color=RED, fontsize=9.5)
    axL.set_xlabel(r"time ($\mu$s)"); axL.set_ylabel(r"$\langle n_z\rangle$")
    axL.set_title(f"exact $\\langle n_z\\rangle(t)$ from $\\bar n_0\\approx{nt[0]:.1f}$ "
                  f"(1/e $\\approx${tau_e:.0f} $\\mu$s, multi-modal)")
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
    fig.savefig(os.path.join(HERE, "images", "cooling_curve.png"), dpi=150)
    print("wrote cooling_curve.png")


# ---------------------------------------------------------------- 3. Stark manifold
def stark_manifold():
    th = 90.0                                      # the REAL trap: transverse lattice (theta=90 deg)
    U0 = -stark.shift(c.alpha0_5S)                 # ground trap depth (down)
    sca = stark.stark_level(2, 0, th)              # F'2: pure scalar (tensor-null), any geometry
    f3hi = stark.stark_level(3, 3, th)             # F'3 stretched |3,+-3> -> HIGHEST at theta=90
    f3lo = stark.stark_level(3, 0, th)             # F'3 |3,0>            -> lowest at theta=90

    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.axhline(0, color="k", lw=0.8, alpha=0.4)
    ax.text(0.5, 2, "free-atom levels (no light)", fontsize=8, color="k", alpha=0.5, ha="center")
    # ground
    ax.hlines(-U0, 0.1, 0.9, color=BLUE, lw=3)
    ax.text(0.5, -U0 - 3, f"5S$_{{1/2}}$ ground: trapped, $U_0$ = {U0:.0f} MHz", color=BLUE, ha="center", fontsize=9)
    # excited manifold (F'3 band at the real theta=90 trap: stretched on top, |3,0> at the bottom)
    ax.fill_between([1.1, 1.9], f3lo, f3hi, color=RED, alpha=0.12)
    ax.hlines(sca, 1.1, 1.9, color=RED, lw=3)
    ax.text(2.0, sca, f"|F'2,0⟩  pure scalar  +{sca:.0f} MHz  (EIT target)", color=RED, va="center", fontsize=9)
    ax.hlines(f3hi, 1.1, 1.9, color=GREY, lw=1.5, ls="--")
    ax.text(2.0, f3hi, f"F'=3 stretched |3,±3⟩  +{f3hi:.0f}  (highest)", color=GREY, va="center", fontsize=8)
    ax.hlines(f3lo, 1.1, 1.9, color=GREY, lw=1.5, ls="--")
    ax.text(2.0, f3lo, f"F'=3 |3,0⟩  +{f3lo:.0f}  (lowest)", color=GREY, va="center", fontsize=8)
    ax.text(1.5, f3hi + 5, "5P$_{3/2}$: anti-trapped\n(tensor splits F'=3 — stretched UP\nat $\\theta{=}90°$; F'=2 tensor = 0)",
            color=RED, ha="center", fontsize=8)
    ax.set_xlim(0, 4.6); ax.set_ylim(-U0 - 10, f3hi + 18)
    ax.set_ylabel("light shift (2$\\pi$·MHz)")
    ax.set_title("1064 nm lattice (1 W + 1 W): ground trapped, 5P$_{3/2}$ expelled")
    ax.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "images", "stark_manifold.png"), dpi=150)
    print("wrote stark_manifold.png")


# ---------------------------------------------------------------- 4. the Lambda scheme
def lambda_scheme():
    """The clock-EIT Lambda: both legs to |F'2,0>, blue-detuned by Delta; the field-insensitive pair."""
    from matplotlib.patches import FancyArrowPatch
    ye, yv, yg = 3.0, 3.55, 0.0                       # excited, virtual (blue-detuned), ground heights
    fig, ax = plt.subplots(figsize=(7.8, 5.6))

    # excited |e> (solid) and the virtual level the blue-detuned beams point to (dashed)
    ax.plot([-0.6, 0.6], [ye, ye], color=RED, lw=3.4, solid_capstyle="round", zorder=4)
    ax.plot([-0.75, 0.75], [yv, yv], color="#9aa0aa", lw=1.3, ls=(0, (5, 3)), zorder=2)

    # ground states
    for xc, col in [(-2.0, GREEN), (2.0, BLUE)]:
        ax.plot([xc - 0.4, xc + 0.4], [yg, yg], color="#111", lw=3.0, solid_capstyle="round", zorder=3)
        ax.plot(xc, yg, "o", ms=11, color=col, mec="white", mew=1.2, zorder=5)

    # the two legs (to the virtual level, so the blue detuning is visible); probe weak (thin), control strong (thick)
    ax.add_patch(FancyArrowPatch((-2.0, yg), (-0.13, yv), arrowstyle="-|>", mutation_scale=18,
                 color=GREEN, lw=2.2, shrinkA=9, shrinkB=4, zorder=4))
    ax.add_patch(FancyArrowPatch((2.0, yg), (0.13, yv), arrowstyle="-|>", mutation_scale=18,
                 color=BLUE, lw=3.8, shrinkA=9, shrinkB=4, zorder=4))

    # detuning marker
    ax.annotate("", xy=(0.92, ye), xytext=(0.92, yv), arrowprops=dict(arrowstyle="<->", color="#444", lw=1.3))
    ax.text(1.02, (ye + yv) / 2, r"$\Delta=+45$" + "\n(both blue)", fontsize=9.3, va="center", ha="left", color="#444")

    # labels (excited-state text sits ABOVE the virtual level, clear of the converging beams)
    ax.text(0, yv + 0.34, r"$|e\rangle = |F'{=}2,\, m'{=}0\rangle$", ha="center", va="bottom", fontsize=12.5, color=RED)
    ax.text(0, yv + 0.12, "anti-trapped $+38$ MHz · pure scalar (tensor-null)", ha="center", va="bottom",
            fontsize=8.4, color=RED, alpha=0.85)
    ax.text(-1.62, 1.95, r"probe $\sigma^+$" + "\n(weak, $\Omega_p$)", color=GREEN, ha="center", va="center", fontsize=10.5)
    ax.text(1.64, 1.95, r"control $\sigma^-$" + "\n(strong, $\Omega_c$)", color=BLUE, ha="center", va="center", fontsize=10.5)
    ax.text(-2.0, yg - 0.26, r"$|g_1\rangle=|F{=}1,m{=}{-}1\rangle$", ha="center", va="top", fontsize=11, color="#111")
    ax.text(2.0, yg - 0.26, r"$|g_2\rangle=|F{=}2,m{=}{+}1\rangle$", ha="center", va="top", fontsize=11, color="#111")
    ax.text(-2.0, yg - 0.64, r"$g_F m_F=+\frac{1}{2}$", ha="center", va="top", fontsize=9.5, color=GREEN)
    ax.text(2.0, yg - 0.64, r"$g_F m_F=+\frac{1}{2}$", ha="center", va="top", fontsize=9.5, color=BLUE)

    # ground hyperfine spacing + the dark state / clock note
    ax.annotate("", xy=(-1.55, yg), xytext=(1.55, yg), arrowprops=dict(arrowstyle="<->", color="#bbb", lw=1.0))
    ax.text(0, yg + 0.12, "6.835 GHz (clock splitting)", ha="center", va="bottom", fontsize=8.2, color="#999")
    ax.text(0, -1.18, r"dark state  $|D\rangle \propto \Omega_c|g_1\rangle - \Omega_p|g_2\rangle$"
            "\nboth legs $g_F m_F=+\\frac{1}{2}$  $\\Rightarrow$  first-order $B$-insensitive  ($\\delta_2$ servoed to 0)",
            ha="center", va="top", fontsize=9.2, color="#333",
            bbox=dict(boxstyle="round,pad=0.5", fc="#f6f6f8", ec="#bbb", lw=1.0))

    ax.set_xlim(-3.4, 3.4); ax.set_ylim(-2.05, 4.6)
    ax.set_title("The clock-EIT $\\Lambda$: both legs to $|F'2,0\\rangle$, blue-detuned by $\\Delta$", fontsize=12.5, pad=8)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "images", "lambda_scheme.png"), dpi=150, bbox_inches="tight")
    print("wrote lambda_scheme.png")


if __name__ == "__main__":
    eit_spectrum()
    cooling_curve()
    stark_manifold()
    lambda_scheme()
