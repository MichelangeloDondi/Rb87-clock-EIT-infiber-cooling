# 01 — the 3-level core

The hand-checkable heart of the scheme: a single ⁸⁷Rb atom as a clean 3-level Λ, cooled on its axial motion.
The physics is derived in the [top-level README](../README.md) (§1–§5); this folder is the code and figures.

**What it computes — the idealized EIT cooling floor:**

```
numeric floor  <n_z> = 0.0020  (full recoil: both legs + emission)   (servo point delta2 = 0)
analytic floor  (Gamma/4Delta)^2 = 0.0011
ground-state population P(n=0) ~ 0.999
```

The **0.0020** is the 3-level floor with the **full photon recoil** (both counter-propagating legs + spontaneous
emission) and perfect repumping; the recoil-free mechanism floor is (Γ/4Δ)² = 0.0011. The realistic intrinsic
limit, with the full ⁸⁷Rb decay branching, is 0.0032 — see [`../02_multilevel/`](../02_multilevel/).

**Files**

| file | what it does |
|---|---|
| `config.py` | every physical number (angular, 2π·MHz) |
| `stark.py` | trap depth + the 5P₃/₂ 1064 Stark shifts (closed form) |
| `stark_validate.py` | re-derives the Wigner-6j tensor factors and checks `stark.py` |
| `cooling.py` | the 3-level cooling floor (one `qutip` master equation, ~60 lines) |
| `plots.py` | the three figures (`eit_spectrum`, `cooling_curve`, `stark_manifold`) |

**Run** (seconds): `python src/stark.py` · `python src/stark_validate.py` · `python src/cooling.py` · `python src/plots.py`.
