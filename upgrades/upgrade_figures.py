"""
upgrade_figures.py -- figures for the "possible upgrades" note.

  floor_ladder.png            the floor across delivery configs (the upgrade path), honestly sourced
  level_scheme_dedicated.png  the 24-level D2 scheme with dedicated sigma repumpers ON F'1 -- all
                              detunings + polarisations; the atomic scheme is the SAME for both upgrades
  bench_single_end.png        Upgrade A optical bench: single-end tagged retro + dedicated repumpers
  bench_dual_end.png          Upgrade B optical bench: dual-end double injection + dedicated repumpers

Pure matplotlib (no solves). Run:  python upgrade_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

BLUE, GREEN, RED, ORANGE, GREY = "#1565c0", "#2e7d32", "#c0392b", "#e67e22", "#7f8c8d"
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- 1. the floor ladder
def floor_ladder():
    # (label, floor, color, source)   source: where the number comes from -- be honest
    rows = [
        ("Minimal single-EOM chain\n(leftover-comb repumpers)",       0.103,  RED,    "computed here"),
        ("Upgrade A: single-end retro\n+ dedicated F'1 repumpers",     0.0072, ORANGE, "design target"),
        ("Upgrade B: dual-end (double injection)\n+ dedicated F'1 repumpers", 0.0048, GREEN, "design target"),
        ("Realistic total\n(Upgrade B + anti-trap squeezer)",          0.009,  BLUE,   "design target"),
    ]
    fig, ax = plt.subplots(figsize=(9.8, 5.0))
    y = np.arange(len(rows))[::-1]
    x0 = 0.0026
    for yi, (lab, val, col, src) in zip(y, rows):
        ax.hlines(yi, x0, val, color=col, lw=2.0, alpha=0.7)
        ax.plot(val, yi, "o", ms=13, color=col, zorder=5)
        ax.text(val * 1.13, yi, f"$\\bar n_z\\approx${val:.4f}   ({src})", va="center", ha="left",
                fontsize=9.5, color=col, fontweight="bold")
    # the target: the ideal mechanism floor (with recoil), computed in this repo
    ax.axvline(0.0032, color="#444", ls=(0, (4, 3)), lw=1.3)
    ax.text(0.0033, 1.5, "ideal mechanism floor  0.0032\n(with recoil, this repo)",
            color="#444", fontsize=8, va="center", ha="left", rotation=90)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xscale("log"); ax.set_xlim(x0, 0.55); ax.set_ylim(-0.6, 3.7)
    ax.set_xlabel(r"axial cooling floor  $\bar n_z$  (lower = colder)")
    ax.set_title("Delivery upgrades — recovering the floor with dedicated F'1 repumpers", pad=12)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "floor_ladder.png"), dpi=150, bbox_inches="tight")
    print("wrote floor_ladder.png")


# ---------------------------------------------------------------- 2. delivery architectures
def _fibre(ax, x0, x1, y):
    ax.add_patch(FancyBboxPatch((x0, y - 0.18), x1 - x0, 0.36,
                 boxstyle="round,pad=0.02,rounding_size=0.1", fc="#eef2f6", ec=GREY, lw=1.1))
    for xx in np.linspace(x0 + 0.3, x1 - 0.3, 12):
        ax.plot([xx, xx], [y - 0.13, y + 0.13], color=GREY, lw=0.5, alpha=0.5)
    ax.plot(0.5 * (x0 + x1), y, "o", ms=11, color="#222", zorder=5)            # the atom
    ax.text(0.5 * (x0 + x1), y - 0.32, "atom", ha="center", va="top", fontsize=8, color="#222")


def _arrow(ax, x0, y0, x1, y1, col, lab="", lx=None, ly=None, dashed=False, ha="center"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
                 lw=2.0, color=col, zorder=4, linestyle="--" if dashed else "-"))
    if lab:
        ax.text(lx, ly, lab, color=col, fontsize=8.2, ha=ha, va="center", fontweight="bold")


def _box(ax, x, y, w, h, text, ec=GREY, fc="#eef2f6", fs=8.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.07",
                 fc=fc, ec=ec, lw=1.5, zorder=4))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color="#1a1a1a", zorder=5, fontweight="bold")


def _mirror(ax, x, y0, y1):
    ax.plot([x, x], [y0, y1], color="#222", lw=4, zorder=4)
    for yy in np.linspace(y0 + 0.08, y1 - 0.08, 6):
        ax.plot([x, x + 0.28], [yy, yy - 0.17], color="#222", lw=0.8, zorder=4)


def _bench_axes(ax):
    ax.set_xlim(0, 14); ax.set_ylim(0, 8)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def bench_single_end():
    """Upgrade A bench: one 1560 EOM->SHG cooling chain into ONE fibre end, retro-tagged
    (lambda/4 + double-passed 200 MHz tag AOM + mirror), plus the dedicated 780 repumper source."""
    fig, ax = plt.subplots(figsize=(13.6, 7.2))
    _bench_axes(ax)
    # cooling-light source (top-left)
    _box(ax, 0.3, 6.05, 3.3, 1.05, "EBLANA 1560 nm  +  $\\phi$-EOM\n$f_{\\rm mod}$=6.835+0.40 = 7.235 GHz", ec=BLUE, fc="#e8f0fe")
    _box(ax, 3.95, 6.05, 2.55, 1.05, "EDFA  +  PPLN\nSHG 1560$\\to$780 nm", ec=BLUE, fc="#e8f0fe")
    _arrow(ax, 3.6, 6.575, 3.95, 6.575, BLUE)
    # dedicated repumper source (lower-left) -- the NEW hardware
    _box(ax, 0.3, 1.0, 2.55, 1.0, "780 master laser\n(F'1 reference)", ec=ORANGE, fc="#fdf2e7")
    _box(ax, 3.25, 1.18, 2.05, 0.66, "157 MHz AOM", ec=ORANGE, fc="#fdf2e7", fs=8.0)
    _arrow(ax, 2.85, 1.5, 3.25, 1.5, ORANGE)
    # combiner
    _box(ax, 6.75, 3.5, 1.1, 1.0, "dichroic\ncombiner", ec=GREY, fc="#f2f2f2", fs=7.7)
    _arrow(ax, 5.2, 6.05, 7.05, 4.5, BLUE, "$\\sigma^-$ carrier\n+ EOM comb", 4.6, 5.05, ha="right")
    _arrow(ax, 5.3, 1.55, 6.95, 3.55, ORANGE, "repump1 $\\sigma^-$ (157 AOM)\nrepump2 $\\sigma^+$ (master)", 5.05, 2.72, ha="right")
    # the fibre + atom
    _fibre(ax, 8.55, 11.35, 4.0)
    _arrow(ax, 7.85, 4.0, 8.55, 4.0, "#333", "ONE\nend", 8.2, 4.46)
    # retro assembly (right): lambda/4 + double-passed tag AOM + mirror
    _box(ax, 11.6, 3.66, 0.55, 0.68, "$\\lambda$/4", ec=RED, fc="#fdecea", fs=9)
    _box(ax, 12.3, 3.55, 1.05, 0.9, "tag AOM\n2$\\times$200=400", ec=RED, fc="#fdecea", fs=7.2)
    _mirror(ax, 13.6, 3.25, 4.75)
    _arrow(ax, 11.35, 4.0, 11.6, 4.0, "#333")
    _arrow(ax, 13.5, 3.5, 8.55, 3.5, RED, "retro $\\to$ $\\sigma^+$ probe (tagged +400) — also leaves rejected comb tones near F'2", 11.0, 3.18, dashed=True, ha="center")
    ax.text(7.0, 0.28, "ONE fibre end + retro mirror. The double-passed tag AOM (2$\\times$200 = 400 MHz) tags the retro as the $\\sigma^+$ probe;\n"
            "the rejected comb tones still sit near F'2 (small residual dark-state scatter)  $\\Rightarrow$  floor $\\approx$ 0.0072.",
            ha="center", va="bottom", fontsize=8.6, color="#333")
    ax.set_title("Upgrade A — single-end tagged retro  +  dedicated F'1 repumpers   ($\\bar n_z\\approx$ 0.0072)", fontsize=12.5, pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "bench_single_end.png"), dpi=150, bbox_inches="tight")
    print("wrote bench_single_end.png")


def bench_dual_end():
    """Upgrade B bench: the common 1560 seed split into two arms -- a clean control (SHG->780)
    into one fibre end, and a carrier-suppressed EOM comb (SHG->780) into the other. No retro, no tag."""
    fig, ax = plt.subplots(figsize=(13.6, 7.2))
    _bench_axes(ax)
    # common seed (top-center)
    _box(ax, 5.0, 6.75, 3.5, 0.8, "EBLANA 1560 nm  (common seed)", ec="#555", fc="#eef2f6", fs=8.6)
    # arm 1 -- clean control (left)
    _box(ax, 0.4, 4.7, 2.95, 0.95, "EDFA + PPLN\nSHG $\\to$ 780", ec=BLUE, fc="#e8f0fe")
    _arrow(ax, 5.3, 6.75, 1.9, 5.65, "#555", "control arm", 3.2, 6.55, ha="center")
    _arrow(ax, 1.9, 4.7, 4.65, 3.2, BLUE, "control  $\\sigma^-$\nfibre END 1", 1.45, 3.55, ha="left")
    # arm 2 -- carrier-suppressed comb (right)
    _box(ax, 9.85, 4.7, 3.55, 0.95, "$\\phi$-EOM (carrier-suppressed)\n+ EDFA + PPLN $\\to$ 780", ec=GREEN, fc="#eaf6ec", fs=8.0)
    _arrow(ax, 8.2, 6.75, 11.6, 5.65, "#555", "comb arm", 10.6, 6.55, ha="center")
    _arrow(ax, 11.6, 4.7, 8.85, 3.2, GREEN, "$\\sigma^+$ comb\nfibre END 2", 12.05, 3.55, ha="right")
    # the fibre + atom (both ends injected)
    _fibre(ax, 4.6, 8.9, 3.0)
    # dedicated repumpers (bottom)
    _box(ax, 4.3, 0.95, 3.4, 0.85, "780 master  +  157 MHz AOM\n(dedicated F'1 repumpers)", ec=ORANGE, fc="#fdf2e7", fs=8.0)
    _arrow(ax, 6.0, 1.8, 6.0, 2.82, ORANGE, "repump1 $\\sigma^-$,  repump2 $\\sigma^+$", 6.2, 2.28, ha="left")
    ax.text(7.0, 0.28, "BOTH fibre ends injected (double injection): control $\\sigma^-$ one end, the carrier-suppressed comb $\\sigma^+$ the other.\n"
            "NO retro mirror, NO tag AOM  $\\Rightarrow$  no rejected tones, full power each way  $\\Rightarrow$  floor $\\approx$ 0.0048 (preferred).",
            ha="center", va="bottom", fontsize=8.6, color="#333")
    ax.set_title("Upgrade B — dual-end double injection  +  dedicated F'1 repumpers   ($\\bar n_z\\approx$ 0.0048)", fontsize=12.5, pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "bench_dual_end.png"), dpi=150, bbox_inches="tight")
    print("wrote bench_dual_end.png")


# ---------------------------------------------------------------- 3. the dedicated-repumper level scheme
def level_scheme_dedicated():
    """24-level D2 scheme for BOTH upgrades: the control+probe Lambda to |F'2,0> PLUS the
    two DEDICATED repumpers, resonant ON F'1 (157 MHz below the cooling F'2). The atomic
    scheme is identical for single-end and dual-end -- only the delivery (the bench) differs."""
    DF = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}       # 5P3/2 hyperfine, from |F'2>
    MP = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
    DIV, YC = 235.0, 5.6
    yE = lambda Fp: YC + DF[Fp] / DIV
    YF2, YF1, LW = 1.15, 0.0, 0.30

    fig, ax = plt.subplots(figsize=(11.8, 8.8))
    # excited manifold: cooling target |F'2,0> (red) and the whole F'=1 row (repump target, orange)
    for Fp, ms in MP.items():
        for m in ms:
            cool = (Fp == 2 and m == 0)
            rep = (Fp == 1)
            col = RED if cool else (ORANGE if rep else "#9aa7b4")
            lw = 3.4 if (cool or rep) else 1.5
            ax.plot([m - LW, m + LW], [yE(Fp), yE(Fp)], color=col, lw=lw, solid_capstyle="round", zorder=3)
        lyoff = -0.16 if Fp == 0 else (0.06 if Fp == 1 else 0.0)
        ax.text(3.8, yE(Fp) + lyoff, r"$F'=%d$" % Fp, va="center", ha="left", fontsize=12, color="#444")
    ax.annotate(r"$|F'2,0\rangle$  $+38$ MHz — COOLING target", xy=(0.33, yE(2)),
                xytext=(1.45, yE(2) + 0.55), fontsize=10.5, color=RED, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color="#999", lw=0.8))
    ax.annotate(r"$F'\!=\!1$  $-157$ MHz — REPUMP target", xy=(1.05, yE(1)),
                xytext=(1.45, yE(1) - 0.58), fontsize=10.5, color="#b5651d", fontweight="bold",
                arrowprops=dict(arrowstyle="-", color="#d2a679", lw=0.8))

    # ground manifold
    for (F, y, ms) in [(2, YF2, range(-2, 3)), (1, YF1, range(-1, 2))]:
        for m in ms:
            ax.plot([m - LW, m + LW], [y, y], color="#111", lw=1.6, solid_capstyle="round", zorder=3)
        ax.text(3.8, y, r"$F=%d$" % F, va="center", ha="left", fontsize=12, color="#444")
    ax.plot(-1, YF1, "o", ms=12, color=GREEN, zorder=5)
    ax.plot(+1, YF2, "o", ms=12, color=BLUE, zorder=5)
    ax.text(-1, YF1 - 0.26, r"$|1,-1\rangle$", ha="center", va="top", fontsize=10, color=GREEN)
    ax.text(+1, YF2 - 0.26, r"$|2,+1\rangle$", ha="center", va="top", fontsize=10, color=BLUE)
    ax.annotate("", xy=(-3.15, YF2), xytext=(-3.15, YF1), arrowprops=dict(arrowstyle="<->", color="#aaa"))
    ax.text(-3.3, (YF1 + YF2) / 2, "6.835 GHz", rotation=90, va="center", ha="right", fontsize=9, color="#888")
    ax.text(-3.0, YF2 + 0.5, "(gaps not to scale)", rotation=90, va="bottom", ha="right", fontsize=7, color="#bbb")

    def beam(x0, y0, x1, y1, col, lab, lx, ly, ha="center"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=15,
                     lw=2.4, color=col, zorder=4, shrinkA=8, shrinkB=6))
        if lab:
            ax.text(lx, ly, lab, color=col, fontsize=9.6, ha=ha, va="center", fontweight="bold")
    # the cooling Lambda (sigma- control, sigma+ probe) -> |F'2,0>
    beam(+1, YF2, 0, yE(2), BLUE, r"CONTROL $\sigma^-$" + "\nF2$\\to$F'2", 1.55, 3.95, "left")
    beam(-1, YF1, 0, yE(2), GREEN, r"PROBE $\sigma^+$" + "\nF1$\\to$F'2", -1.6, 3.95, "right")
    # the dedicated repumpers (sigma- on F=1, sigma+ on F=2) -> ON F'1, stacked 157 MHz under the cooling target
    beam(+1, YF1, 0, yE(1), ORANGE, r"repump1 $\sigma^-$" + "\nF1$\\to$F'1", 1.55, 2.25, "left")
    beam(-1, YF2, 0, yE(1), ORANGE, r"repump2 $\sigma^+$" + "\nF2$\\to$F'1", -1.6, 2.25, "right")

    # detuning markers: Delta=+45 above F'2, and the 157 MHz F'1<->F'2 gap (the whole point)
    ax.annotate("", xy=(0.4, yE(2)), xytext=(0.4, yE(2) + 0.45), arrowprops=dict(arrowstyle="<->", color="#444"))
    ax.text(0.5, yE(2) + 0.24, r"$\Delta\approx+45$ (blue)", fontsize=9.5, va="center", ha="left", color="#444")
    ax.annotate("", xy=(2.7, yE(2)), xytext=(2.7, yE(1)), arrowprops=dict(arrowstyle="<->", color="#b5651d", lw=1.5))
    ax.text(2.82, (yE(1) + yE(2)) / 2, "157 MHz\n(F'1$\\leftrightarrow$F'2)", fontsize=8.6, va="center",
            ha="left", color="#b5651d", fontweight="bold")

    # boxes
    b1 = ("Clock pair: both legs $g_F m_F=+\\frac{1}{2}$  $\\Rightarrow$  the dark\nstate is first-order $B$-insensitive at any field.")
    ax.add_patch(FancyBboxPatch((-3.75, -2.0), 3.5, 0.95, boxstyle="round,pad=0.1",
                 lw=1.1, edgecolor=GREEN, facecolor="#eafaf0", zorder=6))
    ax.text(-2.0, -1.52, b1, fontsize=9, va="center", ha="center", color="#1b5e20", zorder=7)
    b2 = ("DEDICATED repumpers — near-resonant ON F'1 ($\\Delta_{\\rm rep}\\!\\approx\\!5$–$15$ MHz), 157 MHz\n"
          "below the cooling F'2: strong repump, but only weak F'2 dark-state scatter\n"
          "(ratio $\\approx$ 2700). F'1 decays to BOTH hyperfines (5/6$\\to$F1, 1/6$\\to$F2) $\\Rightarrow$ clears\n"
          "the dark sublevels. New 780 source needed (157 MHz is not a 7.23 GHz comb line).")
    ax.add_patch(FancyBboxPatch((-0.12, -2.32), 5.28, 1.5, boxstyle="round,pad=0.1",
                 lw=1.1, edgecolor=ORANGE, facecolor="#fdf2e7", zorder=6))
    ax.text(2.52, -1.57, b2, fontsize=7.9, va="center", ha="center", color="#9c4d0a", zorder=7)

    ax.set_xlim(-4.05, 5.5); ax.set_ylim(-2.4, 8.6)
    ax.set_xlabel(r"$m_F$", fontsize=13); ax.set_xticks(range(-3, 4)); ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title(r"Upgrade level scheme (both configs): the cooling $\Lambda$ + dedicated $\sigma$ repumpers ON F'1",
                 fontsize=12.5, pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "level_scheme_dedicated.png"), dpi=150, bbox_inches="tight")
    print("wrote level_scheme_dedicated.png")


if __name__ == "__main__":
    floor_ladder()
    level_scheme_dedicated()
    bench_single_end()
    bench_dual_end()
