"""
upgrade_figures.py -- figures for the "possible upgrades" note.

  floor_ladder.png         the floor across delivery configs (the upgrade path), honestly sourced
  delivery_comparison.png  single-end tagged retro  vs  dual-end double injection (the architectures)

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
        ("Minimal single-EOM chain\n(leftover-comb repumpers)",       0.103,  RED,    "this repo"),
        ("Upgrade A: single-end retro\n+ dedicated F'1 repumpers",     0.0072, ORANGE, "parent [V]"),
        ("Upgrade B: dual-end (double injection)\n+ dedicated F'1 repumpers", 0.0048, GREEN, "parent [V]"),
        ("Realistic total\n(Upgrade B + anti-trap squeezer)",          0.009,  BLUE,   "parent [V/I]"),
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


def delivery_comparison():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 4.8))

    # --- Panel A: single-end tagged retro ---
    axA.set_xlim(0, 10); axA.set_ylim(0, 6)
    _fibre(axA, 2.0, 8.0, 3.0)
    _arrow(axA, 0.4, 4.4, 3.0, 3.3, BLUE, "control + probe comb\n(one end, σ⁻)", 0.3, 5.3, ha="left")
    _arrow(axA, 0.4, 1.7, 3.0, 2.7, ORANGE, "dedicated F'1 repumpers\n(157 AOM + master)", 0.3, 0.9, ha="left")
    # retro mirror + tag AOM on the right
    axA.plot([8.7, 8.7], [2.2, 3.8], color="#222", lw=3)
    for yy in np.linspace(2.3, 3.7, 6):
        axA.plot([8.7, 9.0], [yy, yy - 0.18], color="#222", lw=0.7)
    axA.add_patch(FancyBboxPatch((7.7, 4.0), 1.6, 0.7, boxstyle="round,pad=0.05",
                  fc="#fdecea", ec=RED, lw=1.0))
    axA.text(8.5, 4.35, "tag AOM $2f_A$", ha="center", va="center", fontsize=8, color=RED)
    _arrow(axA, 8.6, 3.55, 5.2, 3.25, RED, "retro (λ/4 → σ⁺)\n+ rejected tones", 6.6, 4.6, dashed=True, ha="center")
    axA.text(5.0, 0.35, "ONE fibre end + retro. The tag AOM leaves rejected comb tones near F'2\n"
             "(small residual scatter) → floor ≈ 0.0072", ha="center", va="bottom", fontsize=8.5, color="#333")
    axA.set_title("Upgrade A — single-end tagged retro  ($\\bar n_z\\approx$0.0072)", fontsize=11)

    # --- Panel B: dual-end double injection ---
    axB.set_xlim(0, 10); axB.set_ylim(0, 6)
    _fibre(axB, 2.0, 8.0, 3.0)
    _arrow(axB, 0.4, 3.6, 2.0, 3.15, BLUE, "control  σ⁻\n(end 1)", 0.5, 4.5, ha="left")
    _arrow(axB, 0.4, 1.9, 2.0, 2.75, ORANGE, "repumpers", 0.5, 1.2, ha="left")
    _arrow(axB, 9.6, 3.6, 8.0, 3.15, GREEN, "EOM comb  σ⁺\n(carrier-suppressed,\nend 2)", 9.5, 4.7, ha="right")
    _arrow(axB, 9.6, 1.9, 8.0, 2.75, ORANGE, "repumpers", 9.5, 1.2, ha="right")
    axB.text(5.0, 0.35, "BOTH fibre ends injected (double injection); no retro, no tag AOM\n"
             "→ no rejected tones, full power each way → floor ≈ 0.0048", ha="center", va="bottom",
             fontsize=8.5, color="#333")
    axB.set_title("Upgrade B — dual-end double injection  ($\\bar n_z\\approx$0.0048)", fontsize=11)

    for ax in (axA, axB):
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("Both upgrades add DEDICATED F'1 repumpers (resonant on F'1, 157 MHz off the cooling F'2 — "
                 "strong repump, weak dark-state scatter)", fontsize=10.5, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "delivery_comparison.png"), dpi=150, bbox_inches="tight")
    print("wrote delivery_comparison.png")


if __name__ == "__main__":
    floor_ladder()
    delivery_comparison()
