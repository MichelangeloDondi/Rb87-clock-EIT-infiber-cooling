"""
level_scheme.py -- the delivered tones on the 1064-shifted 87Rb 5P3/2 manifold, wired to stark.py.

Generates TWO figures (one script so the two stay consistent):
  level_scheme.png                          (here, 02_multilevel/) -- the BASELINE: 4 comb tones, no master
  ../04_master/images/level_scheme_dedicated.png                            -- the MASTER upgrade: + the master F'1 repumper

Colour = comb line (same beam forward AND retro share a colour): carrier = blue, +1 sideband = green,
780 master = purple.  SOLID = forward pass, DASHED = backward (retro) pass.

Each beam's detuning label is the Stark decomposition  WW(-s-t-g=ZZ):
  WW = detuning from the BARE (1064-OFF) transition;
  s  = excited scalar 1064 shift (+38);   t = excited tensor 1064 shift (signed);
  g  = ground scalar 1064 shift (+23);    ZZ = the IN-TRAP detuning the atom sees.
Each Stark shift raises the F->F' transition, so each SUBTRACTS from the (blue) detuning: ZZ = WW - s - t - g.

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
UPG  = os.path.join(os.path.dirname(HERE), "04_master", "images")
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
Oc = np.sqrt(4.0 * c.Delta * c.nu_z); Op = c.probe_control_ratio * Oc      # control / probe Rabis (~8.8 / 1.1)

em = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
gm = {1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2]}
DELTA = c.Delta
LV   = lt(2, 0) + DELTA                  # control/probe light (true MHz) = +83
L_R1 = LV + 400.0                        # repump1 = fwd +1 EOM sideband (probe + 2f_A)
L_R2 = LV - 400.0                        # repump2 = retro carrier      (control - 2f_A); tag DOWN-shifts
L_MF = lt(1, 1)                          # master fwd: ON the shifted |F'1,+-1>  (the dedicated repump)
L_MR = L_MF - 400.0                      # master retro (down-shifted): benign off-resonant byproduct
yG2, yG1 = -430.0, -495.0               # ground rows (schematic; below the broken axis)
GB, GS = "#c4c4c4", "#9a9a9a"

# SAME colour per comb line (the same beam forward & retro); SOLID = forward, DASHED = backward (retro).
CARRIER, SIDEBAND, MASTER = "#1565c0", "#2e8b3d", "#7b2fb5"   # carrier / +1 sideband / 780 master
C   = dict(control=CARRIER, rep2=CARRIER, rep1=SIDEBAND, probe=SIDEBAND, mfwd=MASTER, mret=MASTER)
FWD = dict(control=True, rep1=True, mfwd=True, probe=False, rep2=False, mret=False)   # forward? else retro

# each beam's (light freq, target Fp, mp) for the Stark-decomposition label WW(-s-t-g=ZZ)
TGT = dict(control=(LV, 2, 0), probe=(LV, 2, 0), rep1=(L_R1, 2, 0),
           rep2=(L_R2, 1, 0), mfwd=(L_MF, 1, -1), mret=(L_MR, 1, 1))
def parts(key):
    """(WW bare-detuning, s exc-scalar, t exc-tensor, g gnd-scalar, ZZ in-trap detuning); ZZ = WW - s - t - g."""
    Llight, Fp, mp = TGT[key]
    s, t, g = SCAL, stark.stark_tensor(Fp, mp, TH), U0
    return (Llight - BARE[Fp] + U0, s, t, g, Llight - lt(Fp, mp))
def label(key):
    WW, s, t, g, ZZ = parts(key)
    tterm = -t if abs(t) > 0.5 else 0.0                    # tensor contribution (0.0, not -0.0, when tensor-null)
    return "%+.0f(%+.0f%+.0f%+.0f=%+.0f)" % (WW, -s, tterm, -g, ZZ)


def draw(with_master, outpath, title):
    fig, ax = plt.subplots(figsize=(18.0, 11.8))
    HW = 0.30
    def lvl(m, y, col="#9aa0aa", lw=2.0, z=2):
        ax.plot([m - HW, m + HW], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=z)
    def beam(p0, p1, key, lw=2.8, msc=16):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=msc, color=C[key], lw=lw,
                     zorder=5, shrinkA=6, shrinkB=5, linestyle=("-" if FWD[key] else "--")))
    def tick(m, y, col):
        ax.plot([m - 0.33, m + 0.33], [y, y], ls=(0, (3, 2)), color=col, lw=1.7, zorder=4)
    def dl(x, y, key, ha="left", va="center"):
        ax.annotate(label(key), xy=(x, y), color=C[key], fontsize=7.5, fontweight="bold", ha=ha, va=va, zorder=9)

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

    # ---- the EIT Lambda: control sigma- (carrier, FWD solid), probe sigma+ (sideband, RETRO dashed) ----
    ax.plot([-0.42, 0.42], [yc(LV), yc(LV)], ls=(0, (4, 3)), color="#777", lw=1.3, zorder=3)
    beam((+1, yG2), (0.05, yc(LV)), "control", 3.0)          # control |2,+1> -> |F'2,0>  dm=-1  (forward)
    beam((-1, yG1), (-0.05, yc(LV)), "probe", 2.6)           # probe   |1,-1> -> |F'2,0>  dm=+1  (retro)
    dl(1.20, yc(40), "control", ha="left")
    dl(-1.18, yc(40), "probe", ha="right")

    # ---- the two EOM-comb repumpers: rep1 (sideband, FWD solid), rep2 (carrier, RETRO dashed) ----
    tick(0.0, yc(L_R1), C["rep1"]); beam((1, yG1), (0.0, yc(L_R1) - 7), "rep1", 2.3)   # |1,+1>->F'2  (forward)
    dl(0.0, yc(L_R1) + 16, "rep1", ha="center", va="bottom")
    tick(0.0, yc(L_R2), C["rep2"]); beam((-1, yG2), (0.0, yc(L_R2) + 7), "rep2", 2.3)  # |2,-1>->F'1  (retro)
    dl(0.42, yc(L_R2) + 1, "rep2", ha="left")

    if with_master:
        # master fwd: sigma+ ON F'1, FORWARD (solid). Clears |2,-2> -- the one F=2 state the sigma- control misses.
        beam((-2, yG2), (-1.0, yE(1, -1)), "mfwd", 2.8)                     # |2,-2>->|F'1,-1>  (forward)
        dl(-2.0, (yG2 + yE(1, -1)) / 2, "mfwd", ha="right")
        ax.annotate("master clears $|2,-2\\rangle$ —\nthe one F=2 state the $\\sigma^-$\ncontrol can't reach",
                    xy=(-2.0, yG2 + 6), xytext=(-3.95, yG2 + 40), color=C["mfwd"], fontsize=7.6, ha="left",
                    va="center", arrowprops=dict(arrowstyle="-", color="#b9a0d0", lw=0.8))
        ax.annotate("the $\\sigma^-$ control clears $|2,+2\\rangle$\n($\\to|F'2,+1\\rangle$) — not a residual",
                    xy=(2.0, yG2), xytext=(2.45, yG2 + 30), color="#666", fontsize=7.4, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.7))
        # master retro: sigma-, BACKWARD (dashed), 400 MHz below F'1 -> benign byproduct
        tick(0.0, yc(L_MR), C["mret"]); beam((1, yG2), (0.0, yc(L_MR) - 7), "mret", 2.0)   # (retro)
        dl(1.20, yc(L_MR) + 1, "mret", ha="left")
        ax.text(0.0, yc(L_MR) - 14, "master retro: benign byproduct (400 off F'1)", color=C["mret"],
                fontsize=7.4, ha="center", va="top")
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

    # ---- beam legend (colour = comb line; solid/dashed = forward/retro; with the decomposition label) ----
    def leglab(name, pol, tr, key):
        return r"%s $%s$ · %s · %s" % (name, pol, tr, label(key))
    leg = [
        Line2D([0], [0], color=CARRIER, lw=3.0, ls="-",  label=leglab("control (fwd)", "\\sigma^-", "F2$\\to$F'2", "control")),
        Line2D([0], [0], color=CARRIER, lw=2.4, ls="--", label=leglab("repump2 (retro)", "\\sigma^+", "F2$\\to$F'1", "rep2")),
        Line2D([0], [0], color=SIDEBAND, lw=3.0, ls="-",  label=leglab("repump1 (fwd)", "\\sigma^-", "F1$\\to$F'2", "rep1")),
        Line2D([0], [0], color=SIDEBAND, lw=2.4, ls="--", label=leglab("probe (retro)", "\\sigma^+", "F1$\\to$F'2", "probe")),
    ]
    if with_master:
        leg += [
            Line2D([0], [0], color=MASTER, lw=2.7, ls="-",  label=leglab("master (fwd)", "\\sigma^+", "F2$\\to$F'1", "mfwd")),
            Line2D([0], [0], color=MASTER, lw=2.3, ls="--", label=leglab("master (retro)", "\\sigma^-", "F2$\\to$F'1", "mret")),
        ]
    ttl = ("beams: colour = comb line (carrier blue, +1 sideband green, master purple); SOLID = forward, DASHED = retro.\n"
           "detuning label  $WW(-s-t-g=ZZ)$:  WW = $\\Delta$ from the BARE (1064-OFF) transition; s,t = excited scalar/tensor\n"
           "shift; g = ground scalar shift; ZZ = the in-trap $\\Delta$.  Each shift raises the transition, so subtracts: ZZ = WW$-$s$-$t$-$g.")
    bl = ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.520, 0.998), frameon=True, fontsize=8.4,
                   handlelength=2.6, borderpad=0.8, labelspacing=0.55, title=ttl, title_fontsize=8.6)
    ax.add_artist(bl)

    # ---- grey-line key ----
    ghost = [Line2D([0], [0], ls=(0, (1, 3)), color=GB, lw=1.4, label="bare — no 1064 shift"),
             Line2D([0], [0], ls=(0, (6, 3)), color=GS, lw=1.3, label="scalar only ($+38$ excited, $-23$ ground)"),
             Line2D([0], [0], ls="-", color="#9aa0aa", lw=2.2, label="full (scalar+tensor; $=$ scalar for F'2/F'0/ground)")]
    ax.legend(handles=ghost, loc="upper left", bbox_to_anchor=(0.002, 0.86), frameon=True, fontsize=8.2,
              handlelength=2.6, borderpad=0.5, labelspacing=0.45, title="grey reference lines", title_fontsize=8.6)

    # ---- delivery / physics box ----
    common = ("$\\bf{single\\ EOM}$: $f_{mod}=A_{HFS}+2f_A=7.23$ GHz; $2f_A=400$ (tag $\\times2$, DOWN-shifts the retro).\n"
              "Each comb line appears twice: FORWARD (solid) and its $-2f_A$ RETRO (dashed), same colour. The carrier is\n"
              "   control (fwd) / repump2 (retro); the $+1$ sideband is repump1 (fwd) / probe (retro).\n")
    if not with_master:
        txt = (common +
               "Repumpers are leftover comb tones, deliberately OFF-resonant: repump1 $\\to$ F1$\\to$F'2; repump2 $\\to$ F2$\\to$F'1.\n"
               "Nearest lines are $\\Delta F=\\pm2$ dipole-FORBIDDEN per ground (F'3 from F=1, F'0 from F=2) $\\Rightarrow$ no coupling.")
    else:
        txt = (common +
               "$\\bf{master\\ upgrade}$: a dedicated F'1 repumper (cooler-frequency slave). master fwd $\\sigma^+$ ON F'1 clears\n"
               "   $|2,\\!-\\!2\\rangle$ — the ONE F=2 sublevel the $\\sigma^-$ control cannot reach ($\\to|F'2,\\!-\\!3\\rangle$ forbidden).\n"
               "$|2,\\!+\\!2\\rangle$ is NOT a residual — the control clears it ($\\to|F'2,\\!+\\!1\\rangle$). The retro is 400 off F'1 (benign).\n"
               "$\\bf{Floor\\ set\\ by\\ the\\ F'1\\ leak}$: the pair is 2-photon resonant on $|F'1,0\\rangle$ too, so the dark state\n"
               "   scatters there ($\\to$5/6 into F=1). The master can't fix it $\\Rightarrow \\bar n_z\\!\\approx\\!0.06$ (chapter 03).")
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
    print("1064 (theta=%.0f): scalar +%.1f, ground -%.1f" % (TH, SCAL, U0))
    for k in ("control", "probe", "rep1", "rep2", "mfwd", "mret"):
        print("  %-8s %s" % (k, label(k)))
    draw(False, os.path.join(HERE, "images", "level_scheme.png"),
         r"$^{87}$Rb D2 clock-EIT — delivered tones (no master): the EIT $\Lambda$ + two comb repumpers")
    draw(True, os.path.join(UPG, "level_scheme_dedicated.png"),
         r"$^{87}$Rb D2 clock-EIT — the master upgrade: + a dedicated F'1 repumper (master fwd $\sigma^+$, retro benign)")
