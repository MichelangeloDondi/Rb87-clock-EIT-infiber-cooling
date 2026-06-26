"""
level_scheme.py -- the delivered tones on the 1064-shifted 87Rb 5P3/2 manifold, wired to stark.py.

Generates TWO figures (one script so the two stay consistent):
  level_scheme.png                          (here, 02_multilevel/) -- the BASELINE: 4 comb tones, no master
  ../upgrades/level_scheme_dedicated.png                            -- the MASTER upgrade: + the master F'1 repumper

Each beam carries TWO detunings, drawn as   in-trap (off):
  in-trap = detuning from the 1064-SHIFTED level (the real operating detuning, no parentheses);
  (off)   = detuning from the 1064-OFF (bare) transition, IN PARENTHESES.
Their difference is that transition's total 1064 shift = excited shift (UP) + ground shift U0 (DOWN, 22.7).
The grey reference lines show the shift: dotted = bare (no 1064), dashed = scalar-only, solid = full
(scalar+tensor). For F'2, F'0 and the ground the dashed sits on the solid -- those are PURE SCALAR.

All frequencies are 2*pi*MHz.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import stark
import config as c

HERE = os.path.dirname(os.path.abspath(__file__))
UPG  = os.path.join(os.path.dirname(HERE), "upgrades")
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

# ---- 1064 light shifts from stark.py (theta = the real transverse-lattice angle) ----
S = 0.62
def yc(f): return f * S                                    # true 2pi-MHz -> drawn y (compressed)
BARE = {0: -229.17, 1: -156.95, 2: 0.0, 3: 266.65}         # bare 87Rb 5P3/2 hyperfine vs F'2 (Steck)
TH = c.theta_trap
def shift(Fp, mp): return stark.stark_level(Fp, mp, TH)    # full excited 1064 shift = scalar + tensor
def lt(Fp, mp): return BARE[Fp] + shift(Fp, mp)            # in-trap (Stark-shifted) excited level
def yE(Fp, mp): return yc(lt(Fp, mp))
SCAL = shift(2, 0)                                         # +38 scalar (F'2 tensor-null)
U0   = -stark.shift(c.alpha0_5S)                           # +22.7 ground scalar DOWN-shift (both F)
Oc = np.sqrt(4.0 * c.Delta * c.nu_z); Op = c.OmR * Oc      # control / probe Rabis (~8.8 / 1.1)

em = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
gm = {1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2]}
DELTA = c.Delta
LV   = lt(2, 0) + DELTA                  # control/probe light (true MHz) = +83
L_R1 = LV + 400.0                        # repump1 = fwd +1 EOM sideband (probe + 2f_A)
L_R2 = LV - 400.0                        # repump2 = retro carrier      (control - 2f_A); tag DOWN-shifts
L_MF = lt(1, 1)                          # master fwd: ON the shifted |F'1,+1>  (the dedicated repump)
L_MR = L_MF - 400.0                      # master retro (down-shifted): benign off-resonant byproduct
yG2, yG1 = -430.0, -495.0               # ground rows (schematic; below the broken axis)
GB, GS = "#c4c4c4", "#9a9a9a"
C = dict(control="#1565c0", probe="#2e8b3d", rep1="#e8730c", rep2="#d32f2f", mast="#7b2fb5", four="#c2185b")

# two detunings per beam: (from the in-trap SHIFTED level | from the 1064-OFF BARE transition incl. ground U0)
def twonum(Llight, Fp, mp): return (Llight - lt(Fp, mp), Llight - BARE[Fp] + U0)
D = dict(control=twonum(LV, 2, 0), probe=twonum(LV, 2, 0), rep1=twonum(L_R1, 2, 0),
         rep2=twonum(L_R2, 1, 0), mfwd=twonum(L_MF, 1, 1), mret=twonum(L_MR, 1, 1))


def draw(with_master, outpath, title):
    fig, ax = plt.subplots(figsize=(18.0, 11.8))
    HW = 0.30
    def lvl(m, y, col="#9aa0aa", lw=2.0, z=2):
        ax.plot([m - HW, m + HW], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=z)
    def beam(p0, p1, col, lw=2.8, msc=16, dashed=False):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=msc, color=col, lw=lw,
                     zorder=5, shrinkA=6, shrinkB=5, linestyle=("--" if dashed else "-")))
    def tick(m, y, col):
        ax.plot([m - 0.33, m + 0.33], [y, y], ls=(0, (3, 2)), color=col, lw=1.7, zorder=4)
    def dl(x, y, key, col, ha="left", va="center"):
        it, off = D[key]
        ax.annotate("%+.0f (%+.0f)" % (it, off), xy=(x, y), color=col, fontsize=8.4,
                    fontweight="bold", ha=ha, va=va, zorder=9)

    # ---- grey reference lines: bare (dotted) + scalar-only (dashed) for EVERY level family ----
    for Fp, ms in em.items():
        x0, x1 = min(ms) - 0.42, max(ms) + 0.42
        ax.plot([x0, x1], [yc(BARE[Fp])] * 2, ls=(0, (1, 3)), color=GB, lw=1.0, zorder=1)        # bare
        ax.plot([x0, x1], [yc(BARE[Fp] + SCAL)] * 2, ls=(0, (6, 3)), color=GS, lw=0.9, zorder=1)  # scalar-only
    for F, ms in gm.items():                                  # GROUND: scalar pushes it DOWN by U0
        yg = yG2 if F == 2 else yG1
        x0, x1 = min(ms) - 0.42, max(ms) + 0.42
        ax.plot([x0, x1], [yg + yc(U0)] * 2, ls=(0, (1, 3)), color=GB, lw=1.0, zorder=1)          # bare (no 1064)
        ax.plot([x0, x1], [yg] * 2, ls=(0, (6, 3)), color=GS, lw=0.9, zorder=1)                   # scalar (= shifted)
    # mark the ground scalar shift once
    ax.annotate("", xy=(-2.0, yG2), xytext=(-2.0, yG2 + yc(U0)), arrowprops=dict(arrowstyle="<->", color="#777", lw=1.1))
    ax.text(-2.12, yG2 + yc(U0) / 2, "$-%.1f$\nscalar\n(both F)" % U0, color="#666", fontsize=7.8,
            va="center", ha="right", linespacing=1.1)

    # ---- solid (full scalar+tensor) levels ----
    for Fp, ms in em.items():
        for mp in ms:
            hot = (Fp == 2 and mp == 0)
            lvl(mp, yE(Fp, mp), "#c62828" if hot else "#9aa0aa", 4.2 if hot else 2.0, 4 if hot else 2)
    for F, ms in gm.items():
        yg = yG2 if F == 2 else yG1
        for m in ms:
            lvl(m, yg, "#2b2b2b", 2.6)
    ax.plot(-1, yG1, "o", color=C["probe"], ms=11, zorder=6, mec="white", mew=1.1)     # |1,-1>
    ax.plot(+1, yG2, "o", color=C["control"], ms=11, zorder=6, mec="white", mew=1.1)   # |2,+1>

    # ---- broken-axis marker (the 6.835 GHz gap to the ground is not to scale) ----
    ylow = yc(L_MR) if with_master else yc(L_R2)
    yb = (ylow + (yG2 + yc(U0))) / 2
    ax.text(-2.75, yb, "optical gap  +  6.835 GHz  —  not to scale", color="#aaa",
            fontsize=8.4, ha="center", va="center", style="italic")
    for yy in (yb + 11, yb - 11):
        ax.plot([-4.30, -4.10], [yy - 6, yy + 6], color="#999", lw=1.4, clip_on=False)

    # ---- the EIT Lambda (solid): control sigma- (dm=-1), probe sigma+ (dm=+1) ----
    ax.plot([-0.42, 0.42], [yc(LV), yc(LV)], ls=(0, (4, 3)), color="#777", lw=1.3, zorder=3)
    beam((-1, yG1), (-0.05, yc(LV)), C["probe"], 3.0)         # probe   |1,-1> -> |F'2,0>  dm=+1
    beam((+1, yG2), (0.05, yc(LV)), C["control"], 3.0)        # control |2,+1> -> |F'2,0>  dm=-1
    dl(-1.18, yc(40), "probe", C["probe"], ha="right")
    dl(1.20, yc(40), "control", C["control"], ha="left")

    # ---- the two EOM-comb repumpers (present in BOTH figures; off-res, dashed) ----
    tick(0.0, yc(L_R1), C["rep1"]); beam((1, yG1), (0.0, yc(L_R1) - 7), C["rep1"], 2.3, dashed=True)
    dl(0.0, yc(L_R1) + 16, "rep1", C["rep1"], ha="center", va="bottom")     # |1,+1>->|F'2,0> dm=-1
    tick(0.0, yc(L_R2), C["rep2"]); beam((-1, yG2), (0.0, yc(L_R2) + 7), C["rep2"], 2.3, dashed=True)
    dl(0.42, yc(L_R2) + 1, "rep2", C["rep2"], ha="left")                    # |2,-1>->|F'1,0> dm=+1

    if with_master:
        # master fwd: sigma+ ON F'1. Its KEY job is to clear |2,-2> -- the ONE F=2 sublevel the sigma-
        # control cannot reach (|2,-2>->|F'2,-3> is forbidden); with no repump 100% piles there.
        beam((-2, yG2), (-1.0, yE(1, -1)), C["mast"], 2.8)                  # |2,-2>->|F'1,-1> dm=+1 on-res
        dl(-2.0, (yG2 + yE(1, -1)) / 2, "mfwd", C["mast"], ha="right")
        ax.annotate("master clears $|2,-2\\rangle$ —\nthe one F=2 state the $\\sigma^-$\ncontrol can't reach",
                    xy=(-2.0, yG2 + 6), xytext=(-3.95, yG2 + 40), color=C["mast"], fontsize=7.6, ha="left",
                    va="center", arrowprops=dict(arrowstyle="-", color="#b9a0d0", lw=0.8))
        # |2,+2> is NOT a residual: the sigma- control clears it (|2,+2>->|F'2,+1>, near-resonant)
        ax.annotate("the $\\sigma^-$ control clears $|2,+2\\rangle$\n($\\to|F'2,+1\\rangle$) — not a residual",
                    xy=(2.0, yG2), xytext=(2.45, yG2 + 30), color="#666", fontsize=7.4, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.7))
        # master retro: sigma-, down-shifted 400 MHz below F'1 -> benign off-resonant byproduct
        tick(0.0, yc(L_MR), C["four"]); beam((1, yG2), (0.0, yc(L_MR) - 7), C["four"], 2.0, dashed=True)
        ax.text(0.0, yc(L_MR) - 14, "master retro: benign byproduct (400 off F'1)", color=C["four"],
                fontsize=7.4, ha="center", va="top")
        # the REAL residual: F=1 (|1,0>,|1,+1>), cleared only weakly by the off-resonant probe
        ax.annotate("real residual: F=1 ($|1,0\\rangle,|1,\\!+\\!1\\rangle$),\nonly weakly cleared by the probe",
                    xy=(0.5, yG1), xytext=(2.2, yG1 + 22), color="#b5651d", fontsize=7.6, ha="left",
                    va="center", arrowprops=dict(arrowstyle="-", color="#d2a679", lw=0.8))

    # ---- Delta bracket + |F'2,0> callout ----
    ax.annotate("", xy=(0.60, yE(2, 0)), xytext=(0.60, yc(LV)), arrowprops=dict(arrowstyle="<->", color="#333", lw=1.4))
    ax.text(0.70, (yE(2, 0) + yc(LV)) / 2, r"$\Delta=+45$", color="#111", fontsize=10.5, va="center", ha="left", fontweight="bold")
    ax.annotate(r"$|F'2,0\rangle$  $+38.0$ (pure scalar, tensor-null)", xy=(0.33, yE(2, 0)), xytext=(1.32, yE(2, 0) - 20),
                color="#c62828", fontsize=10, arrowprops=dict(arrowstyle="-", color="#c62828", lw=1.0))
    # F'1 / F'3 tensor fans (shown by solid pulling away from the scalar dashed)
    ax.annotate("F'1 fans: $|1,0\\rangle$ $+%.0f$ (up),\n$|1,\\pm1\\rangle$ $+%.0f$ (down)" % (shift(1, 0), shift(1, 1)),
                xy=(-1.0, yE(1, 0)), xytext=(-3.95, yE(1, 0) + 8), color="#666", fontsize=7.8, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.8))
    ax.annotate("F'3 fans: stretched $|3,\\pm3\\rangle$ $+%.0f$\n(highest), $|3,0\\rangle$ $+%.0f$ (lowest)"
                % (shift(3, 3), shift(3, 0)), xy=(3.0, yE(3, 3)), xytext=(4.25, yE(3, 3) + 8), color="#666",
                fontsize=7.8, ha="left", va="center", arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.8))

    # ---- level-family labels + ground kets ----
    ctr = {Fp: float(np.mean([yE(Fp, mp) for mp in em[Fp]])) for Fp in em}
    for Fp in em:
        ax.text(3.62, ctr[Fp], r"$F'=%d$" % Fp, va="center", ha="left", color="#666", fontsize=11)
    ax.text(3.62, yG2, r"$F=2$", va="center", ha="left", color="#222", fontsize=12)
    ax.text(3.62, yG1, r"$F=1$", va="center", ha="left", color="#222", fontsize=12)
    ax.text(-1, yG1 - 20, r"$|1,-1\rangle$", color=C["probe"], ha="center", va="top", fontsize=10.5)
    ax.text(1.12, yG2 - 20, r"$|2,+1\rangle$", color=C["control"], ha="center", va="top", fontsize=10.5)
    ax.annotate("", xy=(-3.45, yG2), xytext=(-3.45, yG1), arrowprops=dict(arrowstyle="<->", color="#999", lw=1.3))
    ax.text(-3.57, (yG1 + yG2) / 2, "6.835 GHz", rotation=90, va="center", ha="right", color="#777", fontsize=9)
    sb0 = yc(60)
    ax.plot([-3.95, -3.95], [sb0, sb0 + yc(100)], color="#444", lw=2.4)
    ax.text(-4.04, sb0 + yc(100) / 2, "100 MHz", rotation=90, va="center", ha="right", color="#444", fontsize=8.4)

    # ---- beam legend (with both detuning numbers) ----
    def lab(name, pol, tr, key, role):
        it, off = D[key]
        return r"%s $%s$ · %s · %+.0f (%+.0f)%s" % (name, pol, tr, it, off, role)
    leg = [
        Line2D([0], [0], color=C["control"], lw=3.0, label=lab("control", "\\sigma^-", "F2$\\to$F'2", "control", " · $\\Omega/2\\pi$=%.1f" % Oc)),
        Line2D([0], [0], color=C["probe"], lw=3.0, label=lab("probe", "\\sigma^+", "F1$\\to$F'2", "probe", " · $\\Omega/2\\pi$=%.1f" % Op)),
        Line2D([0], [0], color=C["rep1"], lw=2.4, ls="--", label=lab("repump1", "\\sigma^-", "F1$\\to$F'2", "rep1", " · fwd EOM sideband")),
        Line2D([0], [0], color=C["rep2"], lw=2.4, ls="--", label=lab("repump2", "\\sigma^+", "F2$\\to$F'1", "rep2", " · retro carrier")),
    ]
    if with_master:
        leg += [
            Line2D([0], [0], color=C["mast"], lw=2.7, label=lab("master fwd", "\\sigma^+", "F2$\\to$F'1", "mfwd", " · ON-res, cooler slave")),
            Line2D([0], [0], color=C["four"], lw=2.3, ls="--", label=lab("master retro", "\\sigma^-", "F2$\\to$F'1", "mret", " · benign byproduct")),
        ]
    ttl = "delivered tones   ·   number = $\\Delta$ from the in-trap level   (in parentheses: $\\Delta$ from the 1064-OFF transition)"
    bl = ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.560, 0.998), frameon=True, fontsize=8.6,
                   handlelength=2.4, borderpad=0.8, labelspacing=0.55, title=ttl, title_fontsize=9.0)
    ax.add_artist(bl)

    # ---- grey-line key ----
    ghost = [Line2D([0], [0], ls=(0, (1, 3)), color=GB, lw=1.4, label="bare — no 1064 shift"),
             Line2D([0], [0], ls=(0, (6, 3)), color=GS, lw=1.3, label="scalar only ($+38$ excited, $-23$ ground)"),
             Line2D([0], [0], ls="-", color="#9aa0aa", lw=2.2, label="full (scalar+tensor; $=$ scalar for F'2/F'0/ground)")]
    ax.legend(handles=ghost, loc="upper left", bbox_to_anchor=(0.002, 0.86), frameon=True, fontsize=8.2,
              handlelength=2.6, borderpad=0.5, labelspacing=0.45, title="grey reference lines", title_fontsize=8.6)

    # ---- Stark / delivery box ----
    common = ("$\\bf{single\\ EOM}$: $f_{mod}=A_{HFS}+2f_A=7.23$ GHz; $2f_A=400$ (tag $\\times2$, DOWN-shifts the retro).\n"
              "Each arrow's numbers $=$ detuning from the in-trap (Stark-SHIFTED) level $|$ ($\\Delta$ from the 1064-OFF\n"
              "   transition). Their difference is that transition's total 1064 shift $=$ excited (up) $+$ ground $U_0{=}{-}23$ (down):\n"
              "   $|F'2,0\\rangle$ $+38$ excited $+23$ ground $=+61$; tensor-null $\\Rightarrow$ same for any geometry.\n")
    if not with_master:
        txt = (common +
               "Repumpers are leftover comb tones, deliberately OFF-resonant: repump1 $=$ fwd $+1$ sideband\n"
               "   $\\to$ F1$\\to$F'2; repump2 $=$ retro carrier $\\to$ F2$\\to$F'1. Nearest lines F'3 / F'0 are\n"
               "   $\\Delta F=\\pm2$ dipole-FORBIDDEN $\\Rightarrow$ no coupling.")
    else:
        txt = (common +
               "$\\bf{master\\ upgrade}$: a dedicated F'1 repumper from a cooler-frequency slave. master fwd $\\sigma^+$\n"
               "   ON F'1 clears $|2,\\!-\\!2\\rangle$ — the ONE F=2 sublevel the $\\sigma^-$ control cannot reach\n"
               "   ($|2,\\!-\\!2\\rangle\\!\\to\\!|F'2,\\!-\\!3\\rangle$ forbidden); with no repump 100% piles there and cooling stops.\n"
               "$|2,\\!+\\!2\\rangle$ is NOT a residual — the $\\sigma^-$ control clears it ($\\to|F'2,\\!+\\!1\\rangle$, near-resonant).\n"
               "Folded into the tag arm, the master retro sits 400 MHz off F'1 $\\Rightarrow$ benign byproduct.\n"
               "$\\bf{The\\ real\\ floor\\ is\\ F{=}1\\!-\\!limited}$: $|1,0\\rangle,|1,\\!+\\!1\\rangle$ accumulate, cleared only weakly by\n"
               "   the off-resonant probe — the intrinsic cost of cooling the multilevel D2 line ($\\bar n_z\\!\\sim\\!0.1$).")
    ax.text(4.45, yc(-150), txt, fontsize=8.2, va="top", ha="left", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.7", fc="#f6f6f8", ec="#888", lw=1.0))

    ax.set_title(title, fontsize=13.2, pad=12)
    ax.set_xlabel("$m_F$", fontsize=12.5)
    ax.set_xticks(range(-3, 4)); ax.set_xticklabels(range(-3, 4))
    ax.set_xlim(-4.3, 10.0); ax.set_ylim(yG1 - 40, yc(L_R1) + 40)
    ax.set_yticks([]); ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(outpath, os.path.dirname(HERE)))


if __name__ == "__main__":
    print("1064 (theta=%.0f): scalar +%.1f, ground -%.1f ; F'1 +%.0f/+%.0f ; F'3 +%.0f..+%.0f"
          % (TH, SCAL, U0, shift(1, 1), shift(1, 0), shift(3, 0), shift(3, 3)))
    for k in ("control", "rep1", "rep2", "mfwd", "mret"):
        print("  %-8s in-trap %+7.1f   1064-OFF %+7.1f" % (k, D[k][0], D[k][1]))
    draw(False, os.path.join(HERE, "level_scheme.png"),
         r"$^{87}$Rb D2 clock-EIT — delivered tones (no master): the EIT $\Lambda$ + two comb repumpers")
    draw(True, os.path.join(UPG, "level_scheme_dedicated.png"),
         r"$^{87}$Rb D2 clock-EIT — the master upgrade: + a dedicated F'1 repumper (master fwd $\sigma^+$, retro benign)")
