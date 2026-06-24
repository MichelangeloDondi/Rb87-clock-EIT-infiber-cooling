# Clock-EIT cooling — Layer 1: one atom at the trap centre

This is the **physics core** of the clock-EIT cooling scheme, stripped to the one case you can hold
entirely in your head and check by hand: **a single ⁸⁷Rb atom sitting at the centre of the trap (on
axis), cooled on its axial motion.** No cloud, no off-axis averaging, no excited-state contaminants, no
retro-reflection parasitics — those are Layer 2 and Layer 3 (last section). The point of Layer 1 is that
*every number in it can be obtained two independent ways — a formula you can evaluate on paper, and a
~60-line script you can read top to bottom* — and the two agree.

**Result up front:** one atom on axis cools to **n̄_z ≈ 0.0013** (≈ **99.9 %** axial ground state). The
formula gives **(Γ/4Δ)² ≈ 0.0011**; the script gives **0.0013**. They agree — that agreement is the whole
deliverable.

---

## 1. What we compute, and the assumptions

We ask: *to what mean axial vibrational number n̄_z does clock-EIT cool one atom held at the trap centre?*

Assumptions, all of which Layer 2/3 will relax:
- One atom, **fixed at the trap centre** — so it sees the on-axis fields only (no radial inhomogeneity).
- We cool the **axial** mode: ν_z = 2π · 430 kHz (the stiff direction; the radial directions are ~80×
  softer and are a Layer-3 problem).
- The atom is a **clean 3-level Λ** — just the two cooling legs and the one excited state. The real F′=1
  neighbour and the delivery's stray tones are Layer 2.
- Lamb–Dicke regime: η = 0.094 ≪ 1 (the atom's motion is small compared with the optical wavelength).

All frequencies below are angular, quoted as ordinary frequencies in 2π·MHz (so "Γ = 6.07" means
2π · 6.07 MHz). Numbers live in [`config.py`](config.py).

---

## 2. The cooling scheme

A Λ on the ⁸⁷Rb D2 line — two ground states driven to **one** common excited state:

```
                         |e> = |F'=2, m'=0>
                          /\
            probe σ⁺    /    \    control σ⁻
           (weak, Ω_p)/      \(strong, Ω_c)
                      /        \
   |g1> = |F=1, m=-1>          |g2> = |F=2, m=+1>
   (Δ = +45 MHz blue of |e>;  two-photon detuning δ₂ = probe − control, servoed to 0)
```

Both ground legs have g_F·m_F = +½, so the dark superposition is **first-order magnetic-field
insensitive at any field** — that is the "clock" property, and the reason this exact pair was chosen over
a stretched pair. (The *why-this-pair* argument is its own topic; here we take the pair as given and ask
only about the cooling floor.)

---

## 3. Why EIT cools — the mechanism

This is the part worth reading slowly; everything quantitative follows from it.

**Electromagnetically induced transparency.** When the two lasers are two-photon resonant (δ₂ = 0) the
atom falls into a **dark state** |D⟩ ∝ Ω_c|g1⟩ − Ω_p|g2⟩ that the light cannot excite — the atom stops
scattering. As you scan the probe across the resonance, the excited-state population (the absorption) is
**zero exactly at δ₂ = 0** (the dark resonance) and has a **narrow bright peak** displaced from it by the
AC-Stark shift the strong control imposes on the bright state, δ_AC = Ω_c²/4Δ.

**Now switch on the motion.** A trapped atom absorbs not just at the carrier but on **motional sidebands**
at ±ν_z: the *red* sideband (−ν_z) removes one vibrational quantum (cooling), the *blue* sideband (+ν_z)
adds one (heating). EIT cools by **lining the spectrum up** so that:

- the **carrier** (no change in motion) sits on the **dark resonance** → it does not scatter → no carrier
  heating;
- the **cooling (red) sideband** sits on the **narrow bright peak** → strong absorption → strong cooling;
- the **heating (blue) sideband** sits in the **transparency window** → suppressed.

So absorption-while-cooling ≫ absorption-while-heating, and the atom ratchets down. EIT's advantage is
that it does **not** need a resolved sideband (ν_z ≪ Γ here — ν_z/Γ ≈ 0.07): the *narrow EIT feature*,
not the natural linewidth, provides the selectivity.

---

## 4. The cooling-resonance condition (one line of algebra)

The bright peak lands on the cooling sideband when its AC-Stark displacement equals the trap frequency:

$$\delta_{AC} = \frac{\Omega_c^2}{4\Delta} = \nu_z \quad\Longrightarrow\quad \Omega_c = \sqrt{4\,\Delta\,\nu_z}.$$

With Δ = 45 and ν_z = 0.430 this is Ω_c = √(4·45·0.430) ≈ **8.8** (2π·MHz). The probe is kept weak,
Ω_p = 0.12 Ω_c, so it perturbs the dark state but doesn't wash out the transparency. That is the entire
operating point: (Δ, Ω_c, Ω_p, δ₂) = (45, 8.8, 1.06, servoed to 0).

---

## 5. The floor (and the by-hand check)

In the Lamb–Dicke regime the motion obeys a simple rate equation: a cooling rate A₋ (red-sideband
absorption) and a heating rate A₊ (the residual blue-sideband absorption + off-resonant carrier scatter +
photon recoil). Detailed balance gives the steady state

$$\bar n_z = \frac{A_+}{A_- - A_+}.$$

With the cooling sideband on the bright peak, A₋ is large and A₊ is the small leftover. The leftover is set
by the natural-linewidth tail of the bright resonance reaching back to the carrier/blue positions, which
scales as (Γ/4Δ)². So the EIT cooling floor is

$$\boxed{\;\bar n_{\min} \;\approx\; \left(\frac{\Gamma}{4\Delta}\right)^2\;}
\qquad=\;\left(\frac{6.07}{4\cdot45}\right)^2 = \left(\frac{6.07}{180}\right)^2 \approx 0.0011.$$

**You can check that on a calculator.** It says: push Δ up (more detuning ⇒ less off-resonant scatter ⇒
lower floor), bounded eventually by photon recoil (each spontaneous emission kicks the atom by ~η²; small
here). That single formula is the supervisable heart of Layer 1.

---

## 6. The number — numerical confirmation ([`cooling.py`](cooling.py))

The *exact* steady state of the driven Λ once it is dressed by the oscillator has no clean closed form, so
**this is the one place we hand off to code.** [`cooling.py`](cooling.py) builds that one master equation
(3 levels ⊗ the oscillator), scans δ₂ to find the servo point, and reports the floor:

```
numeric floor   <n_z> = 0.0013   at delta2 = +0.000 (servo point)
analytic floor  (Gamma/4Delta)^2 = 0.0011
ground-state population P(n=0) ~ 0.999
```

The numeric 0.0013 sits just above the formula's 0.0011 (the difference is finite-η and recoil
corrections the formula drops) — and it matches the full multilevel solver's clean-Λ value of
0.0014–0.0021 in the parent project. Three independent routes, one floor.

Run it:
```bash
pip install numpy qutip      # the only dependencies
python cooling.py            # < 10 s
```

---

## 7. Where the light comes from (the finalised chain)

The two tones (and the repumps) are produced and delivered by the finalised optical chain

```
EBLANA (1560 nm seed) → EOM (sidebands) → EDFA (amplify) → PPLN (SHG → 780 nm)
        → HCPCF (the trap + delivery) → AOM (frequency tag, double-passed) → retro-reflection
```

i.e. the **single-ended tagged-retro** delivery. Layer 1 assumes this chain delivers an *ideal* Λ
(control + probe at the operating point); what the chain adds beyond the ideal Λ — the AOM-tagged
rejected return tones — is a Layer-2 heating term, below.

---

## 8. Scope — what Layer 1 is *not* (the road on from here)

Layer 1 is the clean on-axis floor, **0.0013**. The full experiment adds, in order:

- **Layer 2a — the F′=1 contaminant.** The cooling light also weakly drives the nearby |F′=1⟩ state, which
  leaks and heats. In the parent project this lifts the floor 0.0013 → **~0.0048**.
- **Layer 2b — the tagged-retro parasitics.** The `…→AOM→retro` rejected tones add off-resonant scatter,
  lifting it to **~0.0072** (single-ended tagged), the certified single-atom solve floor (+ a once-only
  anti-trap "squeezer" term ≈ 0.003 → the certified all-in 0.008–0.010).
- **Layer 3 — the cloud.** Atoms sit off-axis at a spread of radii and see a spread of trap frequencies
  and light shifts; the cloud floor is set by the radial temperature, and a flat-top trap removes that
  dependence. This is the hard, still-partly-open part — and it is deliberately **out of Layer 1**.

Each layer re-attaches to this core; none of it changes the mechanism or the (Γ/4Δ)² floor derived here.

---

*Layer 1 is meant to be read in an evening and checked with a calculator and one script. If any line of
the physics doesn't follow, that's a bug in the writing — flag it.*
