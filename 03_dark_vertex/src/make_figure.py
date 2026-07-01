"""
make_figure.py -- the chapter-03 figure: (left) the floor ladder with the F'1 leak folded in,
(right) the leak-cancellation curve.

Pure matplotlib (no solve), like 04_master/master_figures.py. The plotted floors are the output
of cooling_dark_vertex.py (run it to reproduce them); they are written here as constants so the
figure regenerates without qutip. All values are ~10% Delta-/run-dependent -- read them as "a few
x 10^-2", not to the last digit.

Run:  python make_figure.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
IMAGES = os.path.join(os.path.dirname(HERE), "images")   # figures -> the chapter's images/ (this file lives in src/)
R_c, R_p = -0.35, 1.73          # the atomic F'1/F'2 coupling ratios (from cooling_dark_vertex.py)

# --- floors from cooling_dark_vertex.py (a few x 10^-2). Panel (a) compares at a FIXED Delta=45 so the
#     master's effect and the F'1 leak are not confounded with the Delta change; MASTER_OPT notes the
#     leak-aware Delta=30 optimum (the headline). ---
CHAIN      = 0.088  # minimal single-EOM chain (comb), leak ON, Delta=45   (chapter 02's ~0.10)
MASTER     = 0.082  # + dedicated master, comb off, leak ON, Delta=45      (same Delta as the chain)
MASTER_OPT = 0.058  # the master re-optimised to the leak-aware Delta=30   (the ~0.06 headline)
IDEAL      = 0.029  # the no-leak ideal (F'1 leak cancelled), Delta=45
RECOIL     = 0.0032 # intrinsic recoil limit (chapter 01)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.6, 5.0), gridspec_kw=dict(width_ratios=[1.32, 1]))

# ---------- panel (a): the floor ladder ----------
rungs = [
    ("minimal single-EOM chain\n(comb repumpers + F'1 leak), Δ=45",   CHAIN,  "#c0392b", "o", "computed"),
    ("+ dedicated master, Δ=45  (Δ≈30 optimum: %.3f)\n(clears |2,-2>, kills comb scatter)" % MASTER_OPT,
                                                                      MASTER, "#e67e22", "o", "computed"),
    ("if the F'1 leak were cancelled\n(no-leak ideal), Δ=45",          IDEAL,  "#2980b9", "o", "computed"),
    ("intrinsic recoil limit (ch. 01)",                               RECOIL, "#2e7d32", "X", "limit"),
]
y = np.arange(len(rungs))[::-1]
# the leak-limited band (at a FIXED Delta=45): the master is stuck above the no-leak ideal by the F'1 leak
axL.axvspan(IDEAL, MASTER, color="#e67e22", alpha=0.07, zorder=0)
axL.text(np.sqrt(IDEAL * MASTER), y.max() + 0.5,
         "the F'1 leak lives in this gap (Δ=45) —\nno repumper reaches below it",
         ha="center", va="bottom", fontsize=8.6, color="#a0610a", style="italic")
for yi, (label, val, col, mk, kind) in zip(y, rungs):
    axL.plot([RECOIL * 0.8, val], [yi, yi], color=col, lw=2.4, zorder=1, alpha=0.95)
    axL.scatter([val], [yi], s=190, color=col, zorder=3, edgecolor="white", linewidth=1.4, marker=mk)
    tag = "  (computed)" if kind == "computed" else "  (limit)"
    axL.text(val * 1.18, yi, r"$\bar n_z\!\approx\!%.3f$%s" % (val, tag), va="center", ha="left",
             fontsize=10.3, color=col, fontweight=("bold" if kind == "computed" else "normal"))
    axL.text(RECOIL * 0.74, yi, label, va="center", ha="right", fontsize=9.0)
axL.axvline(RECOIL, color="0.5", ls=":", lw=1.2)
axL.set_xscale("log"); axL.set_xlim(RECOIL * 0.62, 0.34); axL.set_ylim(-0.7, y.max() + 1.4)
axL.set_yticks([])
axL.set_xlabel(r"axial cooling floor  $\bar n_z$   (lower = colder)", fontsize=11)
axL.set_title("The master floor, with the F'1 leak folded in", fontsize=12.5)
for s in ("top", "right", "left"):
    axL.spines[s].set_visible(False)

# ---------- panel (b): the cancellation curve (swept at Delta=45) ----------
fs = np.array([-0.5, -0.20, 0.0, 0.5, 1.0])              # scale on the probe's |1,-1>->|F'1,0> edge
fl = np.array([0.050, 0.043, 0.043, 0.054, 0.082])       # floor at each scale (cooling_dark_vertex.py)
IDEAL = 0.029                                            # no-leak ideal (with_e1 off), Delta=45
axR.axhline(IDEAL, color="#2e7d32", ls="--", lw=1.6)
axR.text(0.96, IDEAL * 1.03, "no-leak ideal $\\approx %.3f$" % IDEAL, color="#2e7d32",
         fontsize=9.5, va="bottom", ha="right")
axR.plot(fs, fl, "-o", color="#2980b9", lw=2.0, ms=6, zorder=2)
axR.axvline(1.0, color="0.6", ls=":", lw=1.1)
axR.text(1.0, fl.max(), " physical\n (f = 1)", fontsize=8.6, color="0.4", va="top")
axR.axvline(R_c / R_p, color="#c0392b", ls=":", lw=1.1)
axR.annotate("residual nulled\n$f=R_c/R_p=%.2f$" % (R_c / R_p), xy=(R_c / R_p, 0.043),
             xytext=(R_c / R_p + 0.45, IDEAL + 0.022), fontsize=8.6, color="#c0392b",
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))
axR.scatter([1.0], [0.082], s=120, color="#2980b9", edgecolor="#c0392b", linewidth=2.0, zorder=4)
axR.set_xlabel(r"scale on the probe's $|1,-1\rangle\!\to\!|F'1,0\rangle$ coupling,  $f$", fontsize=10.5)
axR.set_ylabel(r"axial cooling floor  $\bar n_z$", fontsize=10.5)
axR.set_title("The leak is the floor: suppress the probe's\nF'1 coupling, recover the cold (swept at $\\Delta=45$)",
              fontsize=11.5)
axR.grid(alpha=0.25)
for s in ("top", "right"):
    axR.spines[s].set_visible(False)

fig.suptitle(r"Chapter 04 — the second dark vertex: $|D_2\rangle$ is not dark on $|F'1,0\rangle$  "
             r"($R_c=%.2f,\ R_p=%+.2f$)" % (R_c, R_p), fontsize=13.5, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(IMAGES, "dark_vertex_floor.png"), dpi=140, bbox_inches="tight")
print("wrote dark_vertex_floor.png")
