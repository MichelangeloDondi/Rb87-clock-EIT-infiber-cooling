"""
scheme_dry.py -- the "dry" clock-EIT level scheme: the six cooling-phase beams drawn on the
REAL 1064-light-shifted 5P_3/2 manifold (theta = 90 deg, the transverse lattice), with grey
reference lines marking the un-shifted positions.

All light shifts come from stark.py (which reads the polarizabilities in config.py) -- nothing
about the 1064 shift is hardcoded here, so the figure tracks the model. Run from the repo:
    python scheme_dry.py        # writes scheme_dry.png / scheme_dry.pdf

Three line styles per excited level:
    dotted grey = bare (no 1064)   dashed grey = scalar only (+scalar, no tensor)   solid = full
F'=2 and F'=0 get no dashed line: their tensor is null (6j{2 2 2;3/2..}=0), so scalar == full.

CAVEATS (see also the footnote on the figure):
  * The repumpers drawn are the DEDICATED near-resonant F'1 scheme (the planned upgrade):
    repump1 F1->F'1 (+15 above shifted F'1), repump2 F2->F'1 (+5). config.py's SOLVER still uses
    the off-resonant single-EOM comb (rep1 F1->F'2 445-off, rep2 F2->F'1 198-off) -- this figure
    reflects the corrected delivery, not the current solver model.
  * Repump detunings are quoted relative to the 1064-SHIFTED F'1 (what you tune to), not bare F'1.
  * OP helicity is a placeholder (drawn sigma-); confirm against the real state-prep sequence.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config as c          # noqa: E402
import stark                # noqa: E402

plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

# ---- 1064 nm light shifts FROM stark.py (theta = config.theta_trap = the real transverse lattice) ----
TH     = c.theta_trap
SCALAR = stark.shift(c.alpha0_5P32)        # 5P3/2 common scalar (~+38); was hardcoded 38.0
U0     = -stark.shift(c.alpha0_5S)         # ground trap depth (~22.7); ground pulled DOWN
BARE   = {0: -229.17, 1: -156.95, 2: 0.0, 3: 266.65}   # 5P3/2 hyperfine vs F'2 (fixed structure)

def yE(Fp, mp):  return BARE[Fp] + stark.stark_level(Fp, mp, TH)   # full = bare + scalar + tensor
L = lambda Fp, mp: stark.stark_level(Fp, mp, TH)                   # full shift (for the box)

em = {0: [0], 1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2], 3: [-3, -2, -1, 0, 1, 2, 3]}
gm = {1: [-1, 0, 1], 2: [-2, -1, 0, 1, 2]}
DELTA = c.Delta
yV = yE(2, 0) + DELTA
yG2, yG1 = -360.0, -420.0
REP1_DET, REP2_DET = 15.0, 5.0    # dedicated near-resonant F'1 repump detunings (the upgrade; see caveat)

C = dict(control="#1565c0", probe="#2e8b3d", rep1="#e8730c",
         rep2="#d32f2f", op="#7b2fb5", det="#00838f")
GB, GS = "#bcbcbc", "#8a8a8a"          # ghost colours: bare (light dotted), scalar-only (dark dashed)

fig, ax = plt.subplots(figsize=(16.8, 10.6))
HW = 0.30
def lvl(m, y, col="#9aa0aa", lw=2.0, z=2):
    ax.plot([m - HW, m + HW], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=z)

# ---- GHOST reference lines: where levels would sit WITHOUT the 1064 shifts ----
def gspan(ms): return min(ms) - 0.42, max(ms) + 0.42
for Fp, ms in em.items():
    x0, x1 = gspan(ms)
    ax.plot([x0, x1], [BARE[Fp]] * 2, ls=(0, (1, 2.5)), color=GB, lw=1.2, zorder=1)            # bare
    if abs(stark.stark_tensor(Fp, max(ms), TH)) > 1e-6:        # scalar-only line only where tensor != 0
        ax.plot([x0, x1], [BARE[Fp] + SCALAR] * 2, ls=(0, (6, 3)), color=GS, lw=1.1, zorder=1)
for F, ms in gm.items():                                       # ground: bare sits +U0 above the shifted
    x0, x1 = gspan(ms)
    ax.plot([x0, x1], [(yG2 if F == 2 else yG1) + U0] * 2, ls=(0, (1, 2.5)), color=GB, lw=1.2, zorder=1)

# ---- solid (full) levels ----
for Fp, ms in em.items():
    for mp in ms:
        hot = (Fp == 2 and mp == 0)
        lvl(mp, yE(Fp, mp), "#c62828" if hot else "#9aa0aa", 4.2 if hot else 2.0, 4 if hot else 2)
for F, ms in gm.items():
    yg = yG2 if F == 2 else yG1
    for m in ms:
        lvl(m, yg, "#2b2b2b", 2.6)
ax.plot(-1, yG1, "o", color=C["probe"], ms=11, zorder=6, mec="white", mew=1.1)
ax.plot(+1, yG2, "o", color=C["control"], ms=11, zorder=6, mec="white", mew=1.1)

# axis break on the LEFT EDGE only
ax.text(-2.6, -250, "optical gap  +  6.835 GHz  -  not to scale", color="#aaa", fontsize=8.6,
        ha="center", va="center", style="italic")
for yy in (-238, -262):
    ax.plot([-4.18, -4.00], [yy - 6, yy + 6], color="#999", lw=1.4, clip_on=False)

# ---- beams ----
def beam(p0, p1, col, lw=2.8, msc=17):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=msc, color=col, lw=lw,
                 zorder=5, shrinkA=7, shrinkB=6))
ax.plot([-0.42, 0.42], [yV, yV], ls=(0, (4, 3)), color="#777", lw=1.3, zorder=3)
beam((-1, yG1), (-0.05, yV), C["probe"], 3.0)                      # probe   sigma+  F1->F'2, +Delta
beam((+1, yG2), (0.05, yV), C["control"], 3.0)                     # control sigma-  F2->F'2, +Delta
beam((0, yG1), (-1, yE(1, -1) + REP1_DET), C["rep1"], 2.5)         # repump1 sigma-  F1->F'1, +15
beam((-1, yG2), (0, yE(1, 0) + REP2_DET), C["rep2"], 2.5)          # repump2 sigma+  F2->F'1, +5
beam((0, yG2), (-1, yE(2, -1)), C["op"], 2.0)                      # OP sigma (placeholder), on-res
beam((+2, yG2), (+3, yE(3, 3)), C["det"], 2.0)                     # detection sigma+ F2->F'3 cycling

# ---- detuning + shift brackets (values from stark.py / config.py) ----
ax.annotate("", xy=(0.60, yE(2, 0)), xytext=(0.60, yV), arrowprops=dict(arrowstyle="<->", color="#333", lw=1.4))
ax.text(0.70, (yE(2, 0) + yV) / 2, fr"$\Delta=+{DELTA:.0f}$", color="#111", fontsize=10.5,
        va="center", ha="left", fontweight="bold")
ax.annotate("", xy=(-1.95, 0.0), xytext=(-1.95, SCALAR), arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2))
ax.text(-2.05, SCALAR / 2, f"+{SCALAR:.0f}\nscalar", color="#444", fontsize=8.2, va="center", ha="right", linespacing=1.15)
ax.annotate(fr"$|F'2,0\rangle$  $+{SCALAR:.1f}$ (scalar, tensor-null)", xy=(0.33, yE(2, 0)),
            xytext=(1.05, yE(2, 0) - 24), color="#c62828", fontsize=10,
            arrowprops=dict(arrowstyle="-", color="#c62828", lw=1.0))
ax.annotate("", xy=(-1.62, yE(1, -1)), xytext=(-1.62, yE(1, 0)), arrowprops=dict(arrowstyle="<->", color="#666", lw=1.1))
ax.text(-1.72, (yE(1, -1) + yE(1, 0)) / 2, f"{L(1, 0) - L(1, 1):.0f} MHz\ntensor split\n(F'1)",
        color="#555", fontsize=7.6, va="center", ha="right", linespacing=1.15)

ctr = {Fp: float(np.mean([yE(Fp, mp) for mp in em[Fp]])) for Fp in em}
for Fp in em:
    ax.text(3.55, ctr[Fp], fr"$F'={Fp}$", va="center", ha="left", color="#666", fontsize=11)
ax.text(3.55, yG2, r"$F=2$", va="center", ha="left", color="#222", fontsize=12)
ax.text(3.55, yG1, r"$F=1$", va="center", ha="left", color="#222", fontsize=12)
ax.text(-1, yG1 - 24, r"$|1,-1\rangle$", color=C["probe"], ha="center", va="top", fontsize=10.5)
ax.text(1.1, yG2 - 24, r"$|2,+1\rangle$", color=C["control"], ha="center", va="top", fontsize=10.5)
ax.annotate("", xy=(-3.45, yG2), xytext=(-3.45, yG1), arrowprops=dict(arrowstyle="<->", color="#999", lw=1.3))
ax.text(-3.57, (yG1 + yG2) / 2, "6.835 GHz", rotation=90, va="center", ha="right", color="#777", fontsize=9)
ax.plot([-3.95, -3.95], [110, 160], color="#444", lw=2.4)
ax.text(-4.04, 135, "50 MHz", rotation=90, va="center", ha="right", color="#444", fontsize=8.4)

# ---- beam legend (top-right); Rabis from config.py ----
leg = [
 Line2D([0], [0], color=C["control"], lw=3.0, label=fr"control $\sigma^-$ . F2$\to$F'2 . $\Delta=+{DELTA:.0f}$ above $|F'2,0\rangle$ . $\Omega/2\pi={c.Omega_c:.1f}$ MHz"),
 Line2D([0], [0], color=C["probe"],   lw=3.0, label=fr"probe $\sigma^+$ . F1$\to$F'2 . $\Delta=+{DELTA:.0f}$ . $\Omega/2\pi={c.Omega_p:.1f}$ MHz ($\Omega_p/\Omega_c{{=}}{c.OmR:.2f}$)"),
 Line2D([0], [0], color=C["rep1"],    lw=2.6, label=fr"repump1 $\sigma^-$ . F1$\to$F'1 . $+{REP1_DET:.0f}$ above F'1 . $\Omega/2\pi\approx3$ MHz"),
 Line2D([0], [0], color=C["rep2"],    lw=2.6, label=fr"repump2 $\sigma^+$ . F2$\to$F'1 . $+{REP2_DET:.0f}$ above F'1 . $\Omega/2\pi\approx3$ MHz"),
 Line2D([0], [0], color=C["op"],      lw=2.2, label=r"OP $\sigma$ . F2$\to$F'2 . on-resonance . saturating (state prep)"),
 Line2D([0], [0], color=C["det"],     lw=2.2, label=r"detection $\sigma^+$ . F2$\to$F'3 . on-resonance . cycling (readout)"),
]
beam_leg = ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.665, 0.998), frameon=True, fontsize=9.4,
                     handlelength=2.0, borderpad=0.8, labelspacing=0.6,
                     title="beams  (EIT $\\Lambda$ + two $\\sigma$ repumpers + master prep/readout)", title_fontsize=10.0)
ax.add_artist(beam_leg)

# ---- ghost-line key (top-left) ----
ghost = [Line2D([0], [0], ls=(0, (1, 2.5)), color=GB, lw=1.6, label="bare - no 1064 shift"),
         Line2D([0], [0], ls=(0, (6, 3)), color=GS, lw=1.5, label=f"scalar only ($+{SCALAR:.0f}$, no tensor)"),
         Line2D([0], [0], ls="-", color="#9aa0aa", lw=2.4, label="full = scalar $+$ tensor (solid)")]
ax.legend(handles=ghost, loc="upper left", bbox_to_anchor=(0.002, 0.965), frameon=True, fontsize=8.8,
          handlelength=2.8, borderpad=0.6, labelspacing=0.5, title="grey reference lines", title_fontsize=9.0)

# ---- Stark box (values computed from stark.py, so they track config.py) ----
box = ("1064 nm light shifts  (atom at trap centre, theta=%.0f deg)\n"
       "ground 5S1/2:  -%.1f scalar - F=1 & F=2 together (6.835 GHz kept)\n"
       "excited 5P3/2:  +%.1f scalar (common)  +  tensor:\n"
       "   |F'2,m'> = +%.1f   (tensor null - clean target, any geometry)\n"
       "   |F'1,+/-1> = +%.1f ,  |F'1,0> = +%.1f   (%.0f MHz split)\n"
       "   |F'3,+/-3> = +%.1f  ->  |F'3,0> = +%.1f   (%.0f MHz spread)\n"
       "   |F'0,0> = +%.1f   (tensor null)\n"
       "=> F2->|F'2,0> sits +%.1f vs bare (+%.0f up, ground -%.1f down);\n"
       "    probe & control are +%.0f above that.\n"
       "=> split F'1 => repump sigma+- legs see unequal detunings\n"
       "    (the F'=1 inhomogeneity).   [all values 2*pi*MHz]") % (
       TH, U0, SCALAR, SCALAR, L(1, 1), L(1, 0), L(1, 0) - L(1, 1),
       L(3, 3), L(3, 0), L(3, 3) - L(3, 0), L(0, 0), SCALAR + U0, SCALAR, U0, DELTA)
ax.text(4.35, -120, box, fontsize=8.8, va="top", ha="left", linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.7", fc="#f6f6f8", ec="#888", lw=1.0))

# ---- caveat footnote (reconciliation items 1-3) ----
foot = ("Repumpers drawn = DEDICATED near-resonant F'1 scheme (planned upgrade); config.py's solver uses the "
        "off-resonant single-EOM comb (rep1 F1->F'2 445-off, rep2 F2->F'1 198-off).\n"
        "Repump detunings (+15, +5) are vs the 1064-SHIFTED F'1 (what you tune to), not bare F'1.   "
        "OP helicity drawn sigma- is a placeholder - confirm against the real state-prep sequence.")
ax.text(2.6, -462, foot, fontsize=7.8, va="bottom", ha="center", color="#777", style="italic", linespacing=1.4)

ax.set_title(r"$^{87}$Rb D2 clock-EIT cooling - the six beams on the 1064-shifted 5P$_{3/2}$ manifold (grey = un-shifted reference)",
             fontsize=13.3, pad=12)
ax.set_xlabel("$m_F$", fontsize=12.5)
ax.set_xticks(range(-3, 4)); ax.set_xticklabels(range(-3, 4))
ax.set_xlim(-4.3, 9.9); ax.set_ylim(-478, 360)
ax.set_yticks([]); ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "scheme_dry.png"), dpi=160, bbox_inches="tight")
fig.savefig(os.path.join(HERE, "scheme_dry.pdf"), bbox_inches="tight")
print("wrote scheme_dry.png / scheme_dry.pdf")
print(f"  SCALAR=+{SCALAR:.2f}  U0={U0:.2f}  |F'1,0>=+{L(1,0):.2f}  |F'1,1>=+{L(1,1):.2f}  "
      f"split={L(1,0)-L(1,1):.2f}  (theta={TH:.0f})")
