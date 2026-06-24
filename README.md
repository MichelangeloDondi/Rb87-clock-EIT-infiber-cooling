# Clock-EIT cooling — Layer 1: one atom at the trap centre

The physics core of the scheme, reduced to the case you can check by hand and with one short script:
**a single ⁸⁷Rb atom at the centre of the 1064 nm trap, cooled on its axial motion.** No cloud, no
off-axis spread, no excited-state contaminants, no retro parasitics — those are Layers 2–3 (last section).

**Result:** one atom on axis reaches **n̄_z ≈ 0.0013** — about **99.9 % axial ground state**. The hand
formula gives (Γ/4Δ)² ≈ 0.0011; `cooling.py` gives 0.0013; the full project solver gives 0.0014–0.0021.
Three independent routes, one floor.

Frequencies are angular, in 2π·MHz. All numbers live in [`config.py`](config.py).

---

## 1. The 1064 nm trap, and why the excited state is expelled

Two 1064 nm beams, **1 W each, counter-propagating**, make the lattice. The AC light shift of a level of
polarizability α in intensity I is

$$U = -\frac{\alpha\,I}{2\varepsilon_0 c}\,,$$

so α > 0 is pulled **down** (trapped) and α < 0 is pushed **up**. At a lattice antinode the two fields add,
so the intensity is 4× the single-beam peak 2P/πw₀²; for 1 W and w₀ = 19 µm that is I = 7.0×10⁹ W/m².

- **Ground 5S₁/₂**, α₀ = +687 a.u. → pulled down by **U₀ ≈ 22.8 MHz = 1.09 mK**. That is the trap depth;
  from it follow ν_z = 2π·430 kHz and η = 0.094.
- **Excited 5P₃/₂**, α₀ = −1149 a.u. (Chen–Raithel, PRA 92, 060501(R), 2015) — *negative* — so it is pushed
  **up** by 22.8 × (1149/687) ≈ **+38 MHz**. The excited state is **anti-trapped** at 1064 nm.

The tensor polarizability (α₂ = +563 a.u.) splits the excited manifold by m′ — with one clean exception
that the whole scheme leans on:

> The cooling transition's upper state **|F′=2, m′=0⟩ is pure scalar**: the Wigner 6j {2 2 2; 3/2 3/2 3/2}
> = 0 kills the entire F′=2 hyperfine tensor. So |F′2,0⟩ sits at **+38 MHz independent of polarization
> geometry** — a fixed, calculable shift, not a sublevel that wanders with the trap. (The F′=3 levels *do*
> split: +19 MHz for the stretched |3,±3⟩ up to +57 MHz.)

The whole thing is closed-form arithmetic — run [`python stark.py`](stark.py) to get every number above. The
practical consequences: "Δ = +45 MHz" is measured from the *in-trap* |F′2,0⟩, and the brief excursions onto
the anti-trapped excited state during cooling cost a little heating — the once-only "squeezer" term in Layer 2.

---

## 2. The Λ scheme

A Λ on the D2 line, both legs to **one** excited state:

```
                         |e> = |F'=2, m'=0>          (anti-trapped, +38 MHz; pure scalar)
                          /\
            probe σ⁺    /    \    control σ⁻
           (weak, Ω_p)/      \(strong, Ω_c)
                      /        \
   |g1> = |F=1, m=-1>          |g2> = |F=2, m=+1>
```

Both legs are blue-detuned by Δ = +45 MHz; the two-photon detuning δ₂ = (probe − control) is servoed to
zero. Both ground states have g_F·m_F = +½, so the dark state is **first-order field-insensitive** — the
"clock" property, and the reason for this exact pair. Here we take the pair as given and ask only: how low
does it cool?

---

## 3. How EIT cools

At two-photon resonance the atom falls into a **dark state** Ω_c|g1⟩ − Ω_p|g2⟩ that doesn't absorb. Scan
the probe and the absorption is **zero at δ₂ = 0** (the dark resonance) with a **narrow bright peak**
displaced by the control's AC-Stark shift Ω_c²/4Δ.

Add the motion: a trapped atom absorbs on sidebands at ±ν_z (red removes a phonon = cooling, blue adds one
= heating). EIT cools by lining the spectrum up so that

- the **carrier** sits on the dark resonance → no scattering, no carrier heating;
- the **cooling (red) sideband** sits on the bright peak → strong absorption;
- the **heating (blue) sideband** sits in the transparency window → suppressed.

Crucially this needs no resolved sideband (here ν_z/Γ ≈ 0.07): the *narrow EIT feature*, not the natural
linewidth, gives the selectivity. That is the whole point of EIT cooling.

---

## 4. The resonance condition, and the floor

The bright peak lands on the cooling sideband when its AC-Stark displacement equals the trap frequency:

$$\frac{\Omega_c^2}{4\Delta} = \nu_z \;\Rightarrow\; \Omega_c = \sqrt{4\,\Delta\,\nu_z} \approx 8.8 \text{ (2π·MHz)},$$

with the probe kept weak (Ω_p = 0.12 Ω_c). The motion then obeys a rate balance — cooling rate A₋, heating
rate A₊ — with steady state n̄_z = A₊ / (A₋ − A₊). With the cooling sideband on the bright peak, the leftover
heating is the natural-linewidth tail reaching back to the carrier, which scales as (Γ/4Δ)². So

$$\boxed{\;\bar n_{\min} \approx \left(\frac{\Gamma}{4\Delta}\right)^2 = \left(\frac{6.07}{180}\right)^2 \approx 0.0011\;}$$

— check it on a calculator. More detuning ⇒ lower floor, until photon recoil (~η² per scatter) takes over.
That one formula is the supervisable heart of Layer 1.

---

## 5. The number ([`cooling.py`](cooling.py))

The exact steady state of the driven Λ dressed by the oscillator has no clean closed form, so **this is the
one place we use code** — a 60-line `qutip` master equation (3 levels ⊗ oscillator), scanned over the servo
detuning δ₂:

```
numeric floor   <n_z> = 0.0013   at delta2 = +0.000 (servo point)
analytic floor  (Gamma/4Delta)^2 = 0.0011
ground-state population P(n=0) ~ 0.999
```

0.0013 sits just above the formula's 0.0011 (finite-η and recoil), and inside the full multilevel solver's
clean-Λ 0.0014–0.0021. Run: `pip install numpy qutip && python cooling.py` (< 10 s).

---

## 6. The chain, and what comes next

The tones are made and delivered by the finalised chain

```
EBLANA (1560) → EOM → EDFA → PPLN (SHG 780) → HCPCF (trap + delivery) → AOM (tag, ×2 pass) → retro
```

— the single-ended tagged-retro line. Layer 1 assumes it delivers an ideal Λ; the layers on from here:

- **Layer 2a** — the F′=1 contaminant the cooling light also drives: 0.0013 → **~0.0048**.
- **Layer 2b** — the tagged-retro rejected tones: → **~0.0072** (certified single-atom solve; + the
  anti-trap squeezer ≈ 0.003 once → all-in 0.008–0.010).
- **Layer 3** — the cloud: atoms off-axis see a spread of ν_z and light shifts; the floor is set by the
  radial temperature, removable with a flat-top trap. The hard, still-partly-open part — deliberately out
  of Layer 1.

Each re-attaches to this core; none changes the mechanism or the (Γ/4Δ)² floor. If any line of physics here
doesn't follow, that's a writing bug — flag it.
