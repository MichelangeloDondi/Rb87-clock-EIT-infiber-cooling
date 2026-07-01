"""
make_figure.py -- the appendix figure: (left) the canceller verdict (static proxy works, real tone fails),
(right) the leak-aware Delta optimum (the coldest the D2 scheme delivers).

Pure matplotlib (no solve). The plotted values are the output of cancellation_floquet.py (panel A) and of
03_dark_vertex's floor() scanned over Delta (panel B); they are written here as constants so the figure
regenerates without qutip. Numbers are ~10% Delta-/run-dependent -- read them as "a few x 10^-2".

Run:  python make_figure.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
Om_p, R_c, R_p = 2.0, -0.346, 1.732
g_null = Om_p * (R_c - R_p)          # = -4.16, where the static proxy nulls the leak
BASE = 1.64e-4                       # no-canceller scatter (cancellation_floquet.py)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.4, 5.2))

# ---- panel A: the canceller verdict (from cancellation_floquet.py) ----
gx = np.linspace(-8.0, 1.5, 40)
# static proxy: a coherent null -> quadratic in (g - g_null), ~0 at g_null (scatter floor ~1e-16)
stat = np.clip(BASE * ((gx - g_null) / g_null) ** 2, 1e-16, None)
# resonant tone: adds its own off-resonant scatter, ~ base + B*g^2 (rises away from 0)
reso = BASE + 3.25e-5 * gx ** 2
axA.axhline(BASE, color="0.55", ls="--", lw=1.4)
axA.text(1.4, BASE * 1.18, "no canceller", color="0.4", fontsize=9, ha="right")
axA.semilogy(gx, reso, "-", color="#c0392b", lw=2.2, label="RESONANT tone (the real, physical one)")
axA.semilogy(gx, stat, "-", color="#2980b9", lw=2.2, label="static in-frame proxy (the unphysical knob)")
axA.axvline(g_null, color="#2980b9", ls=":", lw=1.1)
axA.annotate("proxy nulls the leak here\n(perfect dark state, scatter$\\to0$)", xy=(g_null, 2e-7),
             xytext=(g_null - 0.3, 3e-6), fontsize=8.4, color="#1a7a4a", ha="right",
             arrowprops=dict(arrowstyle="->", color="#1a7a4a", lw=1.0))
axA.annotate("the real tone only\nADDS scatter", xy=(-6.2, BASE + 3.25e-5 * 6.2 ** 2),
             xytext=(-7.7, 4e-3), fontsize=8.6, color="#c0392b",
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))
axA.set_ylim(1e-7, 6e-3)
axA.set_xlabel(r"canceller amplitude on $|1,-1\rangle\!\to\!|F'1,0\rangle$,  $g_X$", fontsize=10.5)
axA.set_ylabel("total steady-state scatter  (∝ leak)", fontsize=10.5)
axA.set_title("Can a tone cancel the leak? The Floquet verdict\n"
              "principle: yes.  physical tone: no — the $\\Delta{+}157$ beat defeats it", fontsize=11.5)
axA.legend(fontsize=8.6, loc="lower left")
axA.grid(alpha=0.25, which="both")

# ---- panel B: the leak-aware Delta optimum (03_dark_vertex floor scan) ----
Ds   = np.array([15, 20, 25, 30, 35, 45, 55])
fon  = np.array([0.0602, 0.0556, 0.0554, 0.0581, 0.0632, 0.0825, 0.1199])   # floor with the leak (honest)
foff = np.array([0.0394, 0.0327, 0.0301, 0.0292, 0.0290, 0.0298, 0.0312])   # no-leak floor (unreachable)
imin = int(np.argmin(fon))
axB.plot(Ds, fon, "-o", color="#e67e22", lw=2.2, ms=6, label="floor with the F'1 leak (honest)")
axB.plot(Ds, foff, "--o", color="#2e7d32", lw=1.8, ms=5, label="no-leak floor (unreachable)")
axB.fill_between(Ds, foff, fon, color="#e67e22", alpha=0.10)
axB.scatter([Ds[imin]], [fon[imin]], s=200, facecolor="none", edgecolor="#c0392b", linewidth=2.2, zorder=5)
axB.annotate("coldest the D2 scheme delivers\n$\\bar n_z\\approx%.3f$ at $\\Delta\\approx%d$" % (fon[imin], Ds[imin]),
             xy=(Ds[imin], fon[imin]), xytext=(Ds[imin] + 6, fon[imin] + 0.020), fontsize=9.2, color="#c0392b",
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.1))
axB.text(40, 0.107, "leak penalty grows\nwith $\\Delta$  (∝ $\\Delta/(\\Delta{+}157)^2$)",
         fontsize=8.6, color="#a0610a", style="italic", ha="center")
axB.set_xlabel(r"control detuning  $\Delta$  (2$\pi$·MHz)", fontsize=10.5)
axB.set_ylabel(r"axial cooling floor  $\bar n_z$", fontsize=10.5)
axB.set_title("The leak-aware operating point\nsmaller $\\Delta$ is colder once the leak dominates", fontsize=11.5)
axB.set_ylim(0, 0.13)
axB.legend(fontsize=8.8, loc="upper left")
axB.grid(alpha=0.25)

fig.suptitle("The leak frontier — cancellation fails; the D2 scheme bottoms out near 0.055", fontsize=13.5, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "images", "cancellation_verdict.png"), dpi=140, bbox_inches="tight")
print("wrote cancellation_verdict.png")
