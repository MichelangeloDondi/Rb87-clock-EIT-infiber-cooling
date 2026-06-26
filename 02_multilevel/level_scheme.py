"""
level_scheme.py -- the 24-level 87Rb D2 scheme with the four tones the multilevel solve uses:
control + probe (the Lambda to |F'2,0>) and the two off-resonant sigma repumpers (the leftover
comb tones).  Run:  python level_scheme.py     (writes level_scheme.png next to this file)

The hyperfine spacings are atomic data (Steck); the figure is a schematic of the tone placement,
not a Stark calculation (for the 1064 shifts see 01_three_level/stark.py and the README).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, RED, GREY, GREEN, ORANGE = "#1565c0", "#c0392b", "#7f8c8d", "#2e7d32", "#e67e22"


def level_scheme():
    """control + probe (the Lambda to |F'2,0>) and the two sigma repumpers on F->F'1/F'2."""
    DF = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}       # 5P3/2 hyperfine, from |F'2> (MHz)
    MP = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
    DIV, YC = 235.0, 5.6
    yE = lambda Fp: YC + DF[Fp] / DIV
    YF2, YF1, LW = 1.15, 0.0, 0.30

    fig, ax = plt.subplots(figsize=(11.2, 8.4))
    # excited manifold (|F'2,0>, the EIT target, emphasized; the rest are off-resonant)
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
    fig.savefig(os.path.join(HERE, "level_scheme.png"), dpi=150, bbox_inches="tight")
    print("wrote level_scheme.png")


if __name__ == "__main__":
    level_scheme()
