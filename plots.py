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


# ---------------------------------------------------------------- 4. the 24-level D2 scheme
def level_scheme():
    """The full 87Rb D2 manifold (8 ground + 16 excited) with the four tones the FULL solve uses:
    control + probe (the Lambda to |F'2,0>) and the two sigma repumpers on F->F'1."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    DF = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}       # 5P3/2 hyperfine, from |F'2>
    MP = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
    DIV, YC = 235.0, 5.6
    yE = lambda Fp: YC + DF[Fp] / DIV
    YF2, YF1, LW = 1.15, 0.0, 0.30
    ORANGE = "#e67e22"

    fig, ax = plt.subplots(figsize=(11.2, 8.4))
    # excited manifold (only |F'2,0>, the EIT target, emphasized; the rest are off-resonant)
    for Fp, ms in MP.items():
        for m in ms:
            eit = (Fp == 2 and m == 0)
            ax.plot([m - LW, m + LW], [yE(Fp), yE(Fp)], color=(RED if eit else "#9aa7b4"),
                    lw=3.0 if eit else 1.5, solid_capstyle="round", zorder=3)
        ax.text(3.6, yE(Fp), r"$F'=%d$" % Fp, va="center", ha="left", fontsize=12, color="#444")
    ax.annotate(r"$|F'2,0\rangle$  pure scalar  +38 MHz  (EIT target)",
                xy=(0, yE(2)), xytext=(1.15, yE(2) + 0.42), fontsize=10.5, color=RED,
                arrowprops=dict(arrowstyle="-", color="#999", lw=0.8))

    # ground manifold
    for (F, y, ms) in [(2, YF2, range(-2, 3)), (1, YF1, range(-1, 2))]:
        for m in ms:
            ax.plot([m - LW, m + LW], [y, y], color="#111", lw=1.6, solid_capstyle="round", zorder=3)
        ax.text(3.6, y, r"$F=%d$" % F, va="center", ha="left", fontsize=12, color="#444")
    ax.plot(-1, YF1, "o", ms=12, color=GREEN, zorder=5)
    ax.plot(+1, YF2, "o", ms=12, color=BLUE, zorder=5)
    ax.text(-1, YF1 - 0.26, r"$|1,-1\rangle$", ha="center", va="top", fontsize=10, color=GREEN)
    ax.text(+1, YF2 - 0.26, r"$|2,+1\rangle$", ha="center", va="top", fontsize=10, color=BLUE)
    ax.annotate("", xy=(-3.15, YF2), xytext=(-3.15, YF1), arrowprops=dict(arrowstyle="<->", color="#aaa"))
    ax.text(-3.3, (YF1 + YF2) / 2, "6.835 GHz", rotation=90, va="center", ha="right", fontsize=9, color="#888")
    ax.text(-3.0, YF2 + 0.5, "(gaps not to scale)", rotation=90, va="bottom", ha="right", fontsize=7, color="#bbb")

    # control + probe = the Lambda (solid); repumpers = off-resonant comb tones (dashed -> virtual level)
    def beam(x0, y0, x1, y1, col, lab, lx, ly, ha="center", dashed=False, term=False):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=15,
                     lw=2.2, color=col, zorder=4, shrinkA=8, shrinkB=6,
                     linestyle="--" if dashed else "-"))
        if term:
            ax.plot([x1 - 0.33, x1 + 0.33], [y1, y1], color=col, lw=1.5, ls=(0, (3, 2)), zorder=4)
        if lab:
            ax.text(lx, ly, lab, color=col, fontsize=10.0, ha=ha, va="center", fontweight="bold")
    beam(+1, YF2, 0, yE(2), BLUE, r"CONTROL $\sigma^-$" + "\nF2$\\to$F'2 (fwd)", 1.75, 3.3, "left")
    beam(-1, YF1, 0, yE(2), GREEN, r"PROBE $\sigma^+$" + "\nF1$\\to$F'2 (retro)", -1.8, 3.3, "right")
    # repump1 tone sits 178 MHz past F'3, but F=1->F'3 is dF=2 FORBIDDEN -> it drives F=1->F'2, 445 off
    yr1 = yE(3) + 178.0 / DIV
    beam(+1, YF1, 0.0, yr1, ORANGE, "", 0, 0, dashed=True, term=True)
    ax.text(1.45, yr1, "repump1 = fwd EOM sideband\n$\\sigma^-$, probe$+$400; drives F1$\\to$F'2 (445 off)\n"
            "(F1$\\not\\to$F'3: $\\Delta F{=}2$ forbidden)",
            color=ORANGE, fontsize=8.2, ha="left", va="center", fontweight="bold")
    # repump2 tone sits 126 MHz past F'0, but F=2->F'0 is dF=2 FORBIDDEN -> it drives F=2->F'1, 198 off
    yr2 = yE(1) - 198.0 / DIV
    beam(-2, YF2, -1.0, yr2, ORANGE, "", 0, 0, dashed=True, term=True)
    ax.text(-2.78, yr2, "repump2 = retro carrier\n$\\sigma^+$, control$-$400\ndrives F2$\\to$F'1 (198 off)\n"
            "(F2$\\not\\to$F'0: forbidden)",
            color=ORANGE, fontsize=8.2, ha="center", va="center", fontweight="bold")

    # the Delta marker
    ax.annotate("", xy=(0.4, yE(2)), xytext=(0.4, yE(2) + 0.45), arrowprops=dict(arrowstyle="<->", color="#444"))
    ax.text(0.5, yE(2) + 0.24, r"$\Delta\approx45$ (blue)", fontsize=9.5, va="center", ha="left", color="#444")

    # bottom annotation boxes
    b1 = ("Clock pair: both legs $g_F m_F=+\\frac{1}{2}$  $\\Rightarrow$  the dark state is\n"
          "first-order $B$-insensitive at any field.")
    ax.add_patch(FancyBboxPatch((-3.7, -1.85), 3.5, 0.95, boxstyle="round,pad=0.1",
                 lw=1.1, edgecolor=GREEN, facecolor="#eafaf0", zorder=6))
    ax.text(-1.95, -1.37, b1, fontsize=9, va="center", ha="center", color="#1b5e20", zorder=7)
    b2 = ("ONE EOM ($f_{\\rm mod}{=}6.83{+}0.40{=}7.23$ GHz) $+$ one 200 MHz tag AOM ($2f_A{=}400$).\n"
          "The repumpers are the leftover comb tones, deliberately OFF-resonant (slower;\n"
          "raise their power to speed up): probe$-$control$=$6.83 GHz, repump1$=$probe$+$400,\n"
          "repump2$=$probe$-$7.23 GHz.")
    ax.add_patch(FancyBboxPatch((0.05, -2.15), 4.9, 1.25, boxstyle="round,pad=0.1",
                 lw=1.1, edgecolor=ORANGE, facecolor="#fdf2e7", zorder=6))
    ax.text(2.5, -1.52, b2, fontsize=8.3, va="center", ha="center", color="#9c4d0a", zorder=7)

    ax.set_xlim(-4.0, 5.3); ax.set_ylim(-2.25, 8.4)
    ax.set_xlabel(r"$m_F$", fontsize=13); ax.set_xticks(range(-3, 4)); ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title(r"$^{87}$Rb D2 clock-EIT scheme (24 levels): control, probe, and the two $\sigma$ repumpers",
                 fontsize=13, pad=10)
    fig.tight_layout()
    fig.savefig("level_scheme.png", dpi=150, bbox_inches="tight")
    print("wrote level_scheme.png")


if __name__ == "__main__":
    eit_spectrum()
    cooling_curve()
    stark_manifold()
    level_scheme()
