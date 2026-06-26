"""
scheme_delivery.py -- the delivered tones of the single-EOM tagged-retro chain, drawn on the
1064-shifted 87Rb 5P3/2 manifold.  Two figures:

  scheme_no_master    the four cooling tones: the EIT Lambda (probe + control) to |F'2,0>,
                      plus the two off-resonant comb-tone repumpers (rep1, rep2).
  scheme_with_master  the same, plus the 780 master folded onto F'2 (OP / F=2 repump) and its
                      retro-reflected copy -- the "4th repumper".

Everything is wired to stark.py / config.py: NO hand-typed light shifts.  The figure's point is
that |F'2,0> is a PURE-SCALAR target -- F'=2 is tensor-NULL, so every F'=2 m' sits at the same
+38, while F'=1 and F'=3 fan out under the tensor shift.  (stark.py confirms tensor(2,m')=0.)

Drawing conventions, made literal so they can be checked by eye:
  * sigma+ arrows step m -> m+1 ;  sigma- arrows step m -> m-1   (Delta m = +-1 in the drawing).
  * Detunings are COMPUTED (printed below) and labelled against the nearest dipole-ALLOWED line;
    dF=+-2 lines (F1->F'3, F2->F'0) are excluded -- they are closer but forbidden.
  * Every tone is drawn at its true level-frame height (to scale); only the 6.835 GHz ground gap
    is compressed (broken axis).  All frequencies are 2*pi*MHz.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import stark, config as c

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

# --- 1064 light shifts straight from stark.py (theta = real transverse-lattice angle) -------
TH     = c.theta_trap
SCALAR = stark.shift(c.alpha0_5P32)        # +38.0 common 5P3/2 scalar (every F',m')
U0     = -stark.shift(c.alpha0_5S)         # +22.7 ground down-shift (F-independent)
# bare 87Rb 5P3/2 hyperfine vs F'=2 (Steck) -- atomic data, not a light shift
BARE   = {0: -229.17, 1: -156.95, 2: 0.0, 3: 266.65}
def yE(Fp, mp): return BARE[Fp] + SCALAR + stark.stark_tensor(Fp, mp, TH)   # in-trap level
em = {0:[0], 1:[-1,0,1], 2:[-2,-1,0,1,2], 3:[-3,-2,-1,0,1,2,3]}
gm = {1:[-1,0,1], 2:[-2,-1,0,1,2]}

DELTA = c.Delta                            # +45, single-photon detuning blue of |F'2,0>
twofA = c.twofA                            # 400, tag-AOM double pass
yG2, yG1 = -470.0, -550.0                   # ground rows (schematic; 6.835 GHz is compressed)

# --- tone positions in the level frame (a tone "at P" is resonant with a level at height P) --
# probe & control sit DELTA above |F'2,0>; the comb repumpers are +-2fA from them.
P_LAMBDA = yE(2,0) + DELTA
P_rep1   = P_LAMBDA + twofA                 # forward +1 sideband  (probe + 2fA)
P_rep2   = P_LAMBDA - twofA                 # retro carrier        (control - 2fA)
P_mast   = yE(2,0)                          # master, resonant on F'2 (OP / F=2 repump)
P_4th    = P_mast + twofA                   # master retro         (master + 2fA)

# --- nearest dipole-ALLOWED line to a tone, for honest detuning labels ----------------------
ALLOWED = {1: (0,1,2), 2: (1,2,3)}          # F -> F' ; dF=+-2 (F1->F'3, F2->F'0) forbidden
def ctr(Fp): return float(np.mean([yE(Fp,mp) for mp in em[Fp]]))   # manifold centre
def nearest_allowed(P, Fg):
    Fp = min(ALLOWED[Fg], key=lambda f: abs(P - ctr(f)))
    return Fp, P - ctr(Fp)                  # (F', signed detuning to its centre)

C  = dict(control="#1565c0", probe="#2e8b3d", rep1="#e8730c",
          rep2="#d32f2f", mast="#7b2fb5", four="#c2185b")
GB, GS = "#c4c4c4", "#9a9a9a"               # ghost colours: bare (dotted), scalar-only (dashed)


def draw(with_master, fname, title):
    fig, ax = plt.subplots(figsize=(16.6, 12.8))
    HW = 0.30
    def lvl(m, y, col, lw, z=2):
        ax.plot([m-HW, m+HW], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=z)
    def beam(p0, p1, col, lw=2.8, dashed=False, msc=16):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=msc, color=col,
                     lw=lw, zorder=5, shrinkA=6, shrinkB=5, linestyle="--" if dashed else "-"))
    def tone(m, y, col, w=0.40):            # an off-resonant delivered tone (dashed marker)
        ax.plot([m-w, m+w], [y, y], ls=(0,(4,2)), color=col, lw=2.1, zorder=4)

    # ghost reference lines: bare (dotted) + scalar-only (dashed). For F'=0,2 scalar-only==full.
    for Fp, ms in em.items():
        x0, x1 = min(ms)-0.42, max(ms)+0.42
        ax.plot([x0,x1], [BARE[Fp]]*2, ls=(0,(1,3)), color=GB, lw=1.0, zorder=1)
        if Fp not in (0,2):
            ax.plot([x0,x1], [BARE[Fp]+SCALAR]*2, ls=(0,(6,3)), color=GS, lw=0.9, zorder=1)
    for F, ms in gm.items():
        x0, x1 = min(ms)-0.42, max(ms)+0.42
        ax.plot([x0,x1], [(yG2 if F==2 else yG1)+U0]*2, ls=(0,(1,3)), color=GB, lw=1.0, zorder=1)

    # solid in-trap levels (|F'2,0> highlighted) ----
    for Fp, ms in em.items():
        for mp in ms:
            hot = (Fp==2 and mp==0)
            lvl(mp, yE(Fp,mp), "#c62828" if hot else "#9aa0aa", 4.4 if hot else 2.2, 4 if hot else 2)
    for F, ms in gm.items():
        yg = yG2 if F==2 else yG1
        for m in ms: lvl(m, yg, "#2b2b2b", 2.6)
    ax.plot(-1, yG1, "o", color=C["probe"],   ms=11, zorder=6, mec="white", mew=1.1)
    ax.plot(+1, yG2, "o", color=C["control"], ms=11, zorder=6, mec="white", mew=1.1)

    # the EIT Lambda (resonant, solid; both legs land on the virtual level DELTA above |F'2,0>) -
    ax.plot([-0.45, 0.45], [P_LAMBDA]*2, ls=(0,(4,3)), color="#777", lw=1.3, zorder=3)
    beam((-1, yG1), (-0.05, P_LAMBDA), C["probe"],   3.0)   # probe   sigma+  |1,-1> -> |F'2,0>
    beam((+1, yG2), ( 0.05, P_LAMBDA), C["control"], 3.0)   # control sigma-  |2,+1> -> |F'2,0>

    # the two off-resonant comb repumpers (dashed) ----
    tone(0, P_rep1, C["rep1"]); beam((+1, yG1), (0, P_rep1-9), C["rep1"], 2.3, dashed=True)  # rep1 s- |1,+1>->F'2
    tone(0, P_rep2, C["rep2"]); beam((-1, yG2), (0, P_rep2+9), C["rep2"], 2.3, dashed=True)  # rep2 s+ |2,-1>->F'1
    F1, d1 = nearest_allowed(P_rep1, 1)
    F2, d2 = nearest_allowed(P_rep2, 2)
    ax.text(0, P_rep1+14, "rep1 $\\sigma^-$  %+.0f from F'%d" % (d1, F1), color=C["rep1"],
            fontsize=8.8, ha="center", va="bottom", fontweight="bold")
    ax.text(0.55, P_rep2, "rep2 $\\sigma^+$  %+.0f from F'%d" % (d2, F2), color=C["rep2"],
            fontsize=8.8, ha="left", va="center", fontweight="bold")

    if with_master:
        beam((0, yG2), (1.0, yE(2,1)), C["mast"], 2.6)                                   # master s+ |2,0>->|F'2,+1>
        tone(-1, P_4th, C["four"]); beam((0, yG2), (-1, P_4th-9), C["four"], 2.1, dashed=True)  # 4th s- |2,0>->F'3
        F4, d4 = nearest_allowed(P_4th, 2)
        ax.text(-1, P_4th+14, "4th $\\sigma^-$ (retro)  %+.0f from F'%d" % (d4, F4), color=C["four"],
                fontsize=8.8, ha="center", va="bottom", fontweight="bold")

    # detuning + shift brackets ----
    ax.annotate("", xy=(0.62, yE(2,0)), xytext=(0.62, P_LAMBDA), arrowprops=dict(arrowstyle="<->", color="#333", lw=1.4))
    ax.text(0.72, (yE(2,0)+P_LAMBDA)/2, r"$\Delta=+%.0f$" % DELTA, color="#111", fontsize=10.5,
            va="center", ha="left", fontweight="bold")
    ax.annotate("", xy=(-1.95, 0.0), xytext=(-1.95, SCALAR), arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2))
    ax.text(-2.05, SCALAR/2, "+%.0f\nscalar" % SCALAR, color="#444", fontsize=8.2, va="center", ha="right", linespacing=1.15)
    # the F'=2 tensor-null statement -- the message of the figure
    ax.annotate("$|F'2,0\\rangle$ at $+%.0f$ (pure scalar)\nF'=2 tensor-NULL: every m' degenerate"
                % SCALAR, xy=(2.05, yE(2,2)), xytext=(3.7, yE(2,0)+86), color="#c62828",
                fontsize=9.6, ha="left", va="center", arrowprops=dict(arrowstyle="-", color="#c62828", lw=1.0))
    ax.annotate("", xy=(-1.62, yE(1,-1)), xytext=(-1.62, yE(1,0)), arrowprops=dict(arrowstyle="<->", color="#666", lw=1.1))
    ax.text(-1.72, (yE(1,-1)+yE(1,0))/2, "%.0f MHz\nF'1 tensor\nsplit" % (yE(1,0)-yE(1,-1)),
            color="#555", fontsize=7.6, va="center", ha="right", linespacing=1.15)

    # row labels, ground splitting, scale bar ----
    for Fp in em: ax.text(3.55, ctr(Fp), fr"$F'={Fp}$", va="center", ha="left", color="#666", fontsize=11)
    ax.text(3.55, yG2, r"$F=2$", va="center", ha="left", color="#222", fontsize=12)
    ax.text(3.55, yG1, r"$F=1$", va="center", ha="left", color="#222", fontsize=12)
    ax.text(-1, yG1-26, r"$|1,-1\rangle$", color=C["probe"],   ha="center", va="top", fontsize=10.5)
    ax.text(1.12, yG2-26, r"$|2,+1\rangle$", color=C["control"], ha="center", va="top", fontsize=10.5)
    ax.annotate("", xy=(-3.45, yG2), xytext=(-3.45, yG1), arrowprops=dict(arrowstyle="<->", color="#999", lw=1.3))
    ax.text(-3.57, (yG1+yG2)/2, "6.835 GHz", rotation=90, va="center", ha="right", color="#777", fontsize=9)
    ax.plot([-3.95,-3.95], [110,160], color="#444", lw=2.4)
    ax.text(-4.04, 135, "50 MHz", rotation=90, va="center", ha="right", color="#444", fontsize=8.4)
    # broken-axis marks (the 6.835 GHz ground gap is the only compressed span; it sits
    # below the lowest off-resonant tone, rep2, and above the ground rows) ----
    ybrk = (P_rep2 + yG2) / 2
    ax.text(-2.55, ybrk, "6.835 GHz  —  not to scale", color="#aaa",
            fontsize=8.4, ha="center", va="center", style="italic")
    for dy in (-12, 12):
        ax.plot([-4.18,-4.00], [ybrk+dy-6, ybrk+dy+6], color="#999", lw=1.4, clip_on=False)

    # beam legend ----
    leg = [
     Line2D([0],[0], color=C["control"], lw=3.0, label=r"control $\sigma^-$ · F2$\to$F'2 · $\Delta=+45$ above $|F'2,0\rangle$ · $\Omega/2\pi=8.8$"),
     Line2D([0],[0], color=C["probe"],   lw=3.0, label=r"probe $\sigma^+$ · F1$\to$F'2 · $\Delta=+45$ · $\Omega/2\pi=1.1$ ($\Omega_p/\Omega_c{=}0.12$)"),
     Line2D([0],[0], color=C["rep1"], lw=2.3, ls="--", label=r"repump1 $\sigma^-$ · F1$\to$F'2 · fwd EOM sideband (probe$+2f_A$)"),
     Line2D([0],[0], color=C["rep2"], lw=2.3, ls="--", label=r"repump2 $\sigma^+$ · F2$\to$F'1 · retro carrier (control$-2f_A$)"),
    ]
    if with_master:
        leg += [
         Line2D([0],[0], color=C["mast"], lw=2.6, label=r"master $\sigma^+$ · F2$\to$F'2 · on-resonance · OP / F=2 repump"),
         Line2D([0],[0], color=C["four"], lw=2.1, ls="--", label=r"4th $\sigma^-$ · master retro (master$+2f_A$) · near F'3 (cycling)"),
        ]
    ttl = ("delivered tones — single seed + one EOM" if not with_master
           else "with the 780 master folded in (+ its retro = 4th repumper)")
    bl = ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.605, 0.998), frameon=True,
                   fontsize=9.2, handlelength=2.4, borderpad=0.8, labelspacing=0.6,
                   title=ttl, title_fontsize=9.8)
    ax.add_artist(bl)

    # ghost-line key ----
    ghost = [Line2D([0],[0], ls=(0,(1,3)), color=GB, lw=1.5, label="bare — no 1064 shift"),
             Line2D([0],[0], ls=(0,(6,3)), color=GS, lw=1.4, label="scalar only (+38, no tensor)"),
             Line2D([0],[0], ls="-", color="#9aa0aa", lw=2.2, label="full = scalar + tensor (solid)")]
    ax.legend(handles=ghost, loc="upper left", bbox_to_anchor=(0.002, 0.80), frameon=True,
              fontsize=8.4, handlelength=2.6, borderpad=0.5, labelspacing=0.45,
              title="grey reference lines", title_fontsize=8.7)

    # physics / delivery box ----
    extra = ("" if not with_master else
             "\n$\\bf{master\\ retro\\ (4th):}$ $\\sigma^-$ at master$+2f_A$.  Nearest ALLOWED line is\n"
             "F'3 (F2$\\to$F'3 cycling), not F'2 -- keep its power low.  It cannot reach the\n"
             "F=1 dark state (6.8 GHz away); the rep2-on-retro placement (master_aom.py)\n"
             "instead lands the retro usefully on F'1.")
    txt = ("$\\bf{single\\ EOM}$:  $f_{mod}=A_{HFS}+2f_A=6.835+0.400=7.235$ GHz;  $2f_A=400$ (tag $\\times2$).\n"
           "The repumpers are leftover comb tones, deliberately OFF-resonant:\n"
           "   rep1 = fwd $+1$ sideband (probe$+2f_A$) $\\to$ F1$\\to$F'2.\n"
           "   rep2 = retro carrier (control$-2f_A$) $\\to$ F2$\\to$F'1.\n"
           "Closer lines F'3 / F'0 are $\\Delta F=\\pm2$ dipole-FORBIDDEN $\\Rightarrow$ no coupling.\n"
           "$\\bf{1064\\ shifts}$ ($\\theta=90°$, from stark.py): ground $-%.1f$ (both F);\n"
           "   $|F'2,0\\rangle=+%.1f$ (6j-null, geometry-free); F'1 split $+%.0f/+%.0f$;\n"
           "   F'3 spread $+%.0f\\to+%.0f$.  [all $2\\pi\\cdot$MHz]"
           % (U0, SCALAR, yE(1,-1), yE(1,0), yE(3,3), yE(3,0)) + extra)
    ax.text(4.30, -150, txt, fontsize=8.5, va="top", ha="left", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.7", fc="#f6f6f8", ec="#888", lw=1.0))

    ax.set_title(title, fontsize=13.2, pad=12)
    ax.set_xlabel("$m_F$", fontsize=12.5)
    ax.set_xticks(range(-3,4)); ax.set_xticklabels(range(-3,4))
    ax.set_xlim(-4.3, 10.0); ax.set_ylim(yG1-55, P_rep1+55)
    ax.set_yticks([]); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, fname+".png"), dpi=158, bbox_inches="tight")
    fig.savefig(os.path.join(HERE, fname+".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(fname, "done")


if __name__ == "__main__":
    # echo the computed tone positions + nearest allowed lines (the detuning audit) ----
    print("1064 (stark.py, theta=%.0f):  scalar %+.2f,  ground -%.2f" % (TH, SCALAR, U0))
    print("|F'2,0> = %+.2f (tensor %.3f -> tensor-NULL, all F'=2 m' degenerate)"
          % (yE(2,0), stark.stark_tensor(2,0,TH)))
    for nm, P, Fg in [("probe/control", P_LAMBDA, 2), ("rep1", P_rep1, 1),
                      ("rep2", P_rep2, 2), ("master", P_mast, 2), ("4th (retro)", P_4th, 2)]:
        Fp, d = nearest_allowed(P, Fg)
        print("  %-14s P=%+7.1f   nearest allowed F'%d  (%+.0f from its centre)" % (nm, P, Fp, d))
    draw(False, "scheme_no_master",
         "$^{87}$Rb D2 clock-EIT — the four delivered tones (single seed + one EOM)")
    draw(True,  "scheme_with_master",
         "$^{87}$Rb D2 clock-EIT — master folded in + its retro (the 4th repumper)")
