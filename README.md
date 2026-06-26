# Clock-EIT cooling of a single ⁸⁷Rb atom in a 1064 nm fibre trap

EIT sideband cooling of a magnetic-field-insensitive ⁸⁷Rb **clock pair** on the **axial** motion of a 1064 nm
optical-lattice trap inside a hollow-core photonic-crystal fibre. The Λ legs are
|F=1,m=−1⟩ σ⁺ and |F=2,m=+1⟩ σ⁻, both to |F′=2,m′=0⟩; both have g_F·m_F = +½, so the dark state is
first-order field-insensitive (the "clock" property). **The question this repo answers: how low does the
axial motion cool?**

> **Status:** work in progress. The 3-level core is settled and hand-checkable; the multilevel layer is
> realistic but its repumper model has a stated validity limit (below). Numbers are single-atom and on-axis —
> the atom cloud is not modelled here.

## The numbers, honestly

| model | floor n̄_z | what it includes | where |
|---|---|---|---|
| 3-level Λ (idealized) | **0.0013** | recoil-free lower bound; perfect repumping | [`01_three_level/`](01_three_level/) |
| multilevel, clean Λ | **0.0032** | full ⁸⁷Rb manifold **+ photon recoil** — the realistic mechanism floor | [`02_multilevel/`](02_multilevel/) |
| multilevel, real delivery | **≈ 0.10** | **+ the real off-resonant repumping** (≈ 40 % stuck in dark sublevels) | [`02_multilevel/`](02_multilevel/) |

Quote **0.0032** as the mechanism floor and **≈ 0.10** as what the minimal single-EOM chain delivers. For this
chain the **repumping**, not the EIT mechanism, sets the floor — [`upgrades/`](upgrades/README.md) shows how
dedicated repumpers recover it. The 0.0013 is the idealized 3-level number: a lower bound, not a result.

All frequencies are angular, in 2π·MHz (a literal `6.07` means 2π·6.07 MHz). Every physical number lives in the
`config.py` of each folder.

---

## 1. The 1064 nm trap, and why the excited state is expelled

Two 1064 nm beams, **1 W each, counter-propagating**, make the lattice. The AC light shift of a level of
polarizability α in intensity I is

$$U = -\frac{\alpha\,I}{2\varepsilon_0 c}\,,$$

so α > 0 is pulled **down** (trapped) and α < 0 is pushed **up**. At a lattice antinode the two fields add, so
the intensity is 4× the single-beam peak 2P/πw₀²; for 1 W and w₀ = 19 µm that is I = 7.0×10⁹ W/m².

- **Ground 5S₁/₂**, α₀ = +687 a.u. → pulled down by **U₀ ≈ 22.7 MHz = 1.09 mK**. That is the trap depth; from
  it follow ν_z = 2π·430 kHz and η = 0.094.
- **Excited 5P₃/₂**, α₀ = −1149 a.u. (Chen–Raithel, PRA 92, 060501(R), 2015) — *negative* — so it is pushed
  **up** by 22.7 × (1149/687) ≈ **+38 MHz**. The excited state is **anti-trapped** at 1064 nm.

The tensor polarizability (α₂ = +563 a.u.) splits the excited manifold by m′ — with one clean exception the
whole scheme leans on:

> The cooling transition's upper state **|F′=2, m′=0⟩ is pure scalar**: the Wigner 6j {2 2 2; 3/2 3/2 3/2} = 0
> kills the entire F′=2 hyperfine tensor. So |F′2,0⟩ sits at **+38 MHz independent of polarization geometry** —
> a fixed, calculable shift, not a sublevel that wanders with the trap. (The F′=3 levels *do* split: +19 MHz
> for the stretched |3,±3⟩ up to +57 MHz.)

![1064 nm light shifts: ground trapped, the 5P₃/₂ manifold expelled](01_three_level/stark_manifold.png)

*The 1064 nm light shifts, every number from [`stark.py`](01_three_level/stark.py). The ground state is pulled
into a 23 MHz (1.1 mK) well; the whole 5P₃/₂ manifold is pushed up (anti-trapped). The EIT target |F′2,0⟩ sits
at the pure-scalar +38 MHz — fixed by the 6j-null, the same in any polarization geometry.*

Run [`python stark.py`](01_three_level/stark.py) for every number above; [`stark_validate.py`](01_three_level/stark_validate.py)
re-derives the Wigner-6j factors from scratch and checks them. "Δ = +45 MHz" is measured from the *in-trap*
|F′2,0⟩.

---

## 2. The Λ scheme

A Λ on the D2 line, both legs to **one** excited state:

![the clock-EIT Λ scheme](01_three_level/lambda_scheme.png)

*Both legs are blue-detuned by Δ = +45 MHz; the two-photon detuning δ₂ = (probe − control) is servoed to zero.
Both ground states have g_F·m_F = +½, so the dark state Ω_c|g₁⟩ − Ω_p|g₂⟩ is **first-order field-insensitive**
— the "clock" property, and the reason for this exact pair. From [`plots.py`](01_three_level/plots.py).*

---

## 3. How EIT cools

At two-photon resonance the atom falls into a **dark state** Ω_c|g1⟩ − Ω_p|g2⟩ that doesn't absorb. Scan the
probe and the absorption is **zero at δ₂ = 0** (the dark resonance) with a **narrow bright peak** displaced by
the control's AC-Stark shift Ω_c²/4Δ.

Add the motion: a trapped atom absorbs on sidebands at ±ν_z (red removes a phonon = cooling, blue adds one =
heating). EIT cools by lining the spectrum up so that

- the **carrier** sits on the dark resonance → no scattering, no carrier heating;
- the **cooling (red) sideband** sits on the bright peak → strong absorption;
- the **heating (blue) sideband** sits in the transparency window → suppressed.

Crucially this needs no resolved sideband (here ν_z/Γ ≈ 0.07): the *narrow EIT feature*, not the natural
linewidth, gives the selectivity. That is the whole point of EIT cooling.

![EIT spectrum: dark carrier, bright cooling sideband, suppressed heating](01_three_level/eit_spectrum.png)

*Absorption (excited population) vs the two-photon detuning δ₂. Zero at the carrier (the dark resonance), a
narrow bright peak parked on the cooling sideband at +ν_z, and the heating sideband at −ν_z left in the
transparency window. Computed by [`plots.py`](01_three_level/plots.py).*

---

## 4. The resonance condition, and the floor

The bright peak lands on the cooling sideband when its AC-Stark displacement equals the trap frequency:

$$\frac{\Omega_c^2}{4\Delta} = \nu_z \;\Rightarrow\; \Omega_c = \sqrt{4\,\Delta\,\nu_z} \approx 8.8 \text{ (2π·MHz)},$$

with the probe kept weak (Ω_p = 0.12 Ω_c). The motion then obeys a rate balance — cooling rate A₋, heating
rate A₊ — with steady state n̄_z = A₊ / (A₋ − A₊). With the cooling sideband on the bright peak, the leftover
heating is the natural-linewidth tail reaching back to the carrier, scaling as (Γ/4Δ)². So

$$\boxed{\;\bar n_{\min} \approx \left(\frac{\Gamma}{4\Delta}\right)^2 = \left(\frac{6.07}{180}\right)^2 \approx 0.0011\;}$$

— check it on a calculator. More detuning ⇒ lower floor, until photon recoil (~η² per scatter) takes over.
That one formula is the supervisable heart of the scheme.

---

## 5. The number ([`01_three_level/cooling.py`](01_three_level/cooling.py))

The exact steady state of the driven Λ dressed by the oscillator has no clean closed form, so **this is the one
place the 3-level core uses code** — a ~60-line `qutip` master equation (3 levels ⊗ oscillator), scanned over
the servo detuning δ₂:

```
numeric floor   <n_z> = 0.0013   at delta2 = +0.000 (servo point)
analytic floor  (Gamma/4Delta)^2 = 0.0011
ground-state population P(n=0) ~ 0.999
```

0.0013 sits just above the formula's 0.0011, and just below the full multilevel solver's clean-Λ **0.0032**
(§6) — the gap is the photon recoil this 3-level model leaves out.

![cooling curve from a hot start, and the final motional distribution](01_three_level/cooling_curve.png)

*Left: starting hot (n̄₀ ≈ 2.8), the axial motion cools to the floor in ~140 µs. Right: at steady state
essentially all the population is in the motional ground state, P(n=0) ≈ 0.999. From
[`plots.py`](01_three_level/plots.py).*

---

## 6. The real manifold and the delivery ([`02_multilevel/cooling_multilevel.py`](02_multilevel/cooling_multilevel.py))

The clean 3-level Λ idealises two things; this layer puts them back — the full ⁸⁷Rb D2 manifold (8 ground
sublevels + the 5P₃/₂ levels), real Clebsch–Gordan couplings and photon recoil. It is a standard multilevel
Lindblad solve; the CG / line-strength conventions are checked against the known D2 branching by
[`cg_validate.py`](02_multilevel/cg_validate.py), and the per-(F′,m′) 1064 Stark comes from the same
[`stark.py`](02_multilevel/stark.py) as §1. The tones are made and delivered by the finalised chain

```
EBLANA (1560) → EOM → EDFA → PPLN (SHG 780) → HCPCF (trap + delivery) → AOM (tag, ×2 pass) → retro
```

a **single seed and one EOM**: f_mod = A_HFS + 2f_A = 6.83 + 0.40 = 7.23 GHz, with a 200 MHz tag AOM
double-passed to 2f_A = 400 MHz. The tag **down-shifts** the retro (retro = forward − 2f_A).

![the 24-level scheme: control, probe, and the two off-resonant repumpers](02_multilevel/level_scheme.png)

*The full D2 manifold and the four tones. control σ⁻ and probe σ⁺ form the Λ to the pure-scalar |F′2,0⟩. The
two repumpers are not separate lasers — they are the **leftover comb tones**: the forward +1 EOM sideband
(σ⁻, F=1, at **probe + 400 MHz**) and the retro carrier (σ⁺, F=2, at **control − 400 MHz**), both deliberately
off-resonant. Their nearest **allowed** lines are F′2 (445 MHz, repump1) and F′1 (198 MHz, repump2); the closer
F′3 / F′0 are ΔF=±2 dipole-forbidden, so they don't couple. From [`level_scheme.py`](02_multilevel/level_scheme.py).*

**(i) Manifold + recoil.** With every m-sublevel, the full recoil, and the per-(F′,m′) 1064 Stark, the clean-Λ
floor is **0.0032** — just above the recoil-light 0.0013 of §5 (the difference *is* the recoil the 3-level
model dropped). The EIT mechanism and the (Γ/4Δ)² scaling are untouched.

**(ii) Repumping is essential — and it is the real cost.** Spontaneous decay from F′ spreads population across
both ground hyperfines into sublevels the Λ never addresses; with the repumpers off, the atom pumps **100 %
dark and cooling stops**. The comb-tone repumpers — modelled as **incoherent** off-resonant scattering (the
virtual F′ adiabatically eliminated, so no rotating-frame artifact) — do clear it, but only partly: at the
chain's **natural** power the on-axis floor settles at **≈ 0.10** (≈ 40 % of the population still in uncooled
dark sublevels). **For this minimal chain the repumping, not the EIT mechanism, is the limit.**
*(Scope: the rate Γ(Ω/2)²/(d²+(Γ/2)²) is the low-saturation limit — reliable only for repumper power ≲ natural.
Above that it omits saturation and the a.c.-Stark shift ∝ Ω²/d, so the high-power rise in the script's sweep is
the model breaking, not physics; trust only the natural-power point.)*

**(iii) Why the detunings are large — and why one EOM can't do better.** The repumper detunings are *fixed* by
f_mod and the tag shift 2f_A, and *one* AOM moves repump1 (F=1) and repump2 (F=2) in **opposite** directions,
so you cannot pull both onto a useful line. Worse, every leftover tone lives near the **cooling F′2 manifold**,
and a tone close to F′2 scatters the EIT dark state at a rate that **equals the cooling rate at δ ≈ 200 MHz off
F′2**. So the repumpers *must* sit ≳ 200 MHz off F′2 — the large detunings are that protection. A configuration
sweep ([`explore_configs.py`](02_multilevel/explore_configs.py)) confirms the current choice is the best of
them, and **caps the single-EOM chain near ~0.1**. The way below it is a *separate* manifold — dedicated
repumpers **on F′1** — laid out in [`upgrades/`](upgrades/README.md).

---

## Structure

Built up one stage at a time; each folder is self-contained (its own `config.py`) and readable on its own.

- **[`01_three_level/`](01_three_level/)** — the hand-checkable core: the trap + Stark shifts (`stark.py` +
  `stark_validate.py`), the 3-level cooling floor (`cooling.py`, ~60 lines), the figures (`plots.py`).
  **Start here.** Runs in seconds.
- **[`02_multilevel/`](02_multilevel/)** — the real ⁸⁷Rb manifold, photon recoil, and the two off-resonant
  comb repumpers of the single-EOM chain (`cooling_multilevel.py`), the configuration sweep
  (`explore_configs.py`), and the level scheme (`level_scheme.py`). **The baseline, without the master laser.**
- **[`upgrades/`](upgrades/README.md)** — *forward-looking*: how dedicated F′1 repumpers (and the master laser)
  recover the floor toward the mechanism limit. Not part of the baseline.

## How to run

```bash
pip install -r requirements.txt        # numpy, scipy, qutip, sympy, matplotlib

cd 01_three_level
python stark.py            # trap depth + the 5P3/2 Stark shifts (closed form)
python stark_validate.py   # re-derives the Wigner-6j factors and checks them
python cooling.py          # the 3-level floor  (< 10 s)
python plots.py            # the three figures (real solve, no drawn curves)

cd ../02_multilevel
python level_scheme.py        # the 24-level scheme figure (no solve)
python cooling_multilevel.py  # the realistic floor with recoil + repumping  (~1 min)
python explore_configs.py     # the single-EOM configuration sweep
```

There is no separate test runner: each script prints its own self-check, and the headline floor in §4–§5 is a
formula you can check by hand.

---

## What's beyond this repo

- **Dedicated repumpers / other delivery** ([`upgrades/`](upgrades/README.md)): on-resonant F′1 repumpers and
  dual-end delivery, with **design targets ~0.0048–0.0072** (these are *not* computed in this repo). §6 shows
  *why* you want them: off-resonant comb-tone repumping is the bottleneck.
- **The atom cloud:** atoms off-axis see a spread of ν_z and light shifts; the floor is then set by the radial
  temperature, removable with a flat-top trap. Deliberately out of this single-atom repo.

The §1–5 core is the supervisable heart: the EIT mechanism and the (Γ/4Δ)² floor. §6 is the honest price of the
real delivery. If any line of physics here doesn't follow, that's a writing bug — flag it.
