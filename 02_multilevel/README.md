# 02 — the real ⁸⁷Rb manifold and the delivery

The baseline cooling, **without the master laser**: the full ⁸⁷Rb D2 manifold (8 ground sublevels + the 5P₃/₂
levels), real Clebsch–Gordan couplings and photon recoil, and the two repumpers of the single-EOM chain
(leftover comb tones, deliberately off-resonant). The physics is in the [top-level README](../README.md) (§6).

**What it computes — two floors:**

```
clean Lambda (manifold + recoil)              <n_z> = 0.0032   (the realistic intrinsic cooling limit)
real delivery (off-resonant comb repumpers)   <n_z> ~ 0.10     (~40% stuck in dark sublevels)
```

For this minimal chain the **repumping**, not the EIT mechanism, sets the floor. Going below ~0.1 needs
dedicated F′1 repumpers — see [`../04_master/`](../04_master/README.md).

**Files**

| file | what it does |
|---|---|
| `config.py` | every physical number (manifold + delivery: B, θ, 2f_A, …) |
| `stark.py` | the per-(F′,m′) 1064 Stark shifts (same closed form as 01; validated there) |
| `cooling_multilevel.py` | the multilevel Lindblad floor, repumpers in the computation (~1 min) |
| `clebsch_gordan_checks.py` | reconstructs known ⁸⁷Rb D2 facts from raw Clebsch–Gordan (checks the conventions) |
| `explore_configs.py` | the single-EOM / one-AOM configuration sweep (why the current placement wins) |
| `level_scheme.py` | the 24-level scheme figures on the 1064-shifted manifold — the baseline (no master) here, **and** the chapter-04 variant (written into `../04_master/`); colour = comb line, solid/dashed = forward/retro, each beam's label is the Stark decomposition `WW(−s−t−g=ZZ)` (bare detuning → in-trap) |

**Run:** `python src/level_scheme.py` (no solve) · `python src/cooling_multilevel.py` (~1 min) ·
`python src/clebsch_gordan_checks.py` · `python src/explore_configs.py`.

*Validity note:* the off-resonant repumpers are incoherent low-saturation rates — trust the natural-power point,
not the high-power trend (the rate omits saturation and the a.c.-Stark shift above ~natural power).
