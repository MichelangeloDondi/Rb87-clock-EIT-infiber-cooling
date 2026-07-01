# Clock-EIT cooling of a single ⁸⁷Rb atom in a fibre trap — a numerical study

A **numerical modelling and design study** — theory plus simulation, with experimentally realistic parameters —
for a **planned** EIT sideband-cooling experiment: a single ⁸⁷Rb atom on the **axial** motion of a 1064 nm
optical-lattice trap inside a hollow-core photonic-crystal fibre. It fixes the operating point, the delivery, and
the achievable floor ahead of the build. The Λ legs are |F=1,m=−1⟩ σ⁺ and |F=2,m=+1⟩ σ⁻, both to
|F′=2,m′=0⟩. This is a **g_F·m_F-matched "clock" pair** (both legs have g_F·m_F = +½) — *not* the usual
m_F=0↔m_F=0 clock pair, but it serves the same purpose: a first-order magnetically insensitive two-photon
resonance (chapter 01). **The question: how low does this model predict the axial motion cools?**

> **Status:** work in progress. The 3-level core (ch. 01) is settled and hand-checkable; the multilevel layer
> (ch. 02) is realistic but its repumper model has a stated validity limit; the **F′1 leak** (ch. 03) sets the
> real floor (best no-master ≈ 0.087), and the optional **master** (ch. 04) pushes it to ≈ 0.055 (leak-limited).
> Numbers are single-atom and on-axis — the atom cloud is not modelled here.

## The numbers, honestly

| model | floor n̄_z | what it includes | where |
|---|---|---|---|
| 3-level Λ (idealized) | **0.0020** | full recoil (both legs + emission), perfect repumping; recoil-free mechanism floor (Γ/4Δ)² = 0.0011 | [`01_three_level/`](01_three_level/) |
| multilevel, clean Λ | **0.0032** | full ⁸⁷Rb manifold **+ photon recoil** — the realistic intrinsic cooling limit | [`02_multilevel/`](02_multilevel/) |
| multilevel, real delivery | **≈ 0.09** | **+ the real off-resonant repumping**, at the servoed δ₂-optimum (≈ 40 % stuck in dark sublevels) | [`02_multilevel/`](02_multilevel/) |

Quote **0.0032** as the intrinsic cooling limit and **≈ 0.09** as what the minimal single-EOM chain delivers. For this
chain the **repumping** sets the floor; [`03_dark_vertex/`](03_dark_vertex/README.md) exposes the deeper F′1-leak
limit underneath (best no-master floor ≈ 0.087), and [`04_master/`](04_master/README.md) how far an optional master
laser pushes it (≈ 0.055, leak-limited). The 0.0020 is the idealized 3-level number (recoil-free mechanism floor
0.0011): a lower bound, not a result.

All frequencies are angular, in 2π·MHz (a literal `6.07` means 2π·6.07 MHz). Every physical number lives in the
`src/config.py` of each chapter.

## Why these choices

The three questions a reader asks first:

- **Why EIT?** The trap is weak — ν_z/Γ ≈ 0.07, deep in the *unresolved*-sideband regime, so plain
  resolved-sideband cooling on the bare optical line is out. EIT sidesteps that: it builds its own narrow
  dark-resonance feature (its width set by the drive, not by Γ) and cools a broad range of n **continuously**,
  using the D2 beams themselves — no pulse train, no extra lasers. (Degenerate-Raman sideband cooling *also* works
  in this regime and is the nearest alternative, but it needs *separate* near-detuned Raman beams — see the
  [appendix](appendix/options.md).)
- **Why the D2 line (780 nm)?** A hardware choice: the delivery is a telecom chain seeded at 1560 nm and
  frequency-doubled (1560/2 = 780 nm = D2) — a mature, low-noise source — and on D2 the cooling target |F′2,0⟩ is
  **tensor-null** (chapter 01), a light-shift-stable upper state. The honest trade-off: the **D1** line (795 nm) has
  a ~41× *weaker* F′1 leak and would cool *colder* ([appendix](appendix/d1_comparison.md)), but 795 nm is off the
  telecom-doubling chain. D2 is the practical line; D1 is the physically cleaner one.
- **Why this g_F·m_F-matched pair, not the m_F=0↔m_F=0 clock pair? Geometry decides it.** With B along the fibre
  axis and the beams *also* along the axis, the light is **σ-only** (a transverse wave has no π component along B),
  so the Λ must use two m_F=±1 legs — |1,−1⟩ σ⁺ and |2,+1⟩ σ⁻; a pure m_F=0 Λ would need π light (a transverse B),
  which the geometry forbids. The **bonus**: both legs share the same linear Zeeman shift (g_F·m_F = +½), so it
  cancels in the two-photon resonance — first-order field-insensitive, *like* the m_F=0 clock states, with a small
  quadratic residual the δ₂ servo absorbs (chapter 01). Field-insensitivity is a bonus here, not the reason — the
  reason is the σ-only geometry.

## The chapters

Built up **one layer of complexity at a time** — each chapter adds a single piece of physics or hardware to the one
before, and is **self-contained**: its own `src/config.py`, its own physics write-up, runnable on its own. Read them
in order; each README derives the physics it needs, assuming only the chapters before it.

| # | chapter | what it adds | status |
|---|---|---|---|
| **01** | [`01_three_level/`](01_three_level/README.md) | the idealized 3-level Λ: the 1064 trap, the Stark shifts, the EIT mechanism, the (Γ/4Δ)² floor | **built** · n̄_z = 0.0020 |
| **02** | [`02_multilevel/`](02_multilevel/README.md) | the real ⁸⁷Rb D2 manifold + photon recoil + the single-EOM comb delivery (probe, control, retro) | **built** · 0.0032 clean / ≈ 0.09 real |
| **03** | [`03_dark_vertex/`](03_dark_vertex/README.md) | the second dark vertex — the cooling pair is two-photon resonant on F′1 too, so the dark state leaks — folds it in and finds the best floor **without** the master | **built** · ≈ 0.087 (no-master, repump-limited) |
| **04** | [`04_master/`](04_master/README.md) | the 780 master laser as a dedicated F′1 repumper — a possible upgrade; does it earn the extra hardware? | **built** · ≈ 0.055 (leak-limited) |
| 05 | *(planned)* | the anti-trapping heating from the expelled 5P₃/₂ excited state, computed explicitly (estimated < 2 % in the scope notes) | planned |
| 06 | *(planned)* | the atom cloud (frozen atoms): a spread of ν_z and light shifts off-axis | planned |
| 07 | *(planned)* | the full semiclassical Monte-Carlo simulation | planned |
| 08 | *(planned)* | beam depletion along the fibre (scattering, absorption, …) | planned |

Chapters 05–08 are the roadmap, not yet built.

**Beyond the chapters.** [`appendix/`](appendix/README.md) collects deep-dives on the chapter-03 **F′1 leak**:
the leak can't be cancelled by a co-propagating tone (Floquet test → D2 bottoms at ≈ 0.055); the **D1 line** has a
~41× weaker leak (detuning-driven, not branching — a future from-scratch chapter); the field-noise trade against a
leak-free but field-sensitive stretched scheme; and the full menu of alternatives (incl. a dRSC pivot) with verdicts.

## How to run

```bash
pip install -r requirements.txt        # numpy, scipy, qutip, sympy, matplotlib

cd 01_three_level
python src/stark.py            # trap depth + the 5P3/2 Stark shifts (closed form)
python src/stark_validate.py   # re-derives the Wigner-6j factors and checks them
python src/cooling.py          # the 3-level floor  (< 10 s)
python src/plots.py            # the three figures (real solve, no drawn curves)

cd ../02_multilevel
python src/level_scheme.py        # the 24-level scheme figure (no solve)
python src/cooling_multilevel.py  # the realistic floor with recoil + repumping  (~1 min)
python src/explore_configs.py     # the single-EOM configuration sweep

cd ../03_dark_vertex
python src/cooling_dark_vertex.py # the F′1 leak + the best no-master floor, scanned over Δ  (~minutes; qutip)
python src/make_figure.py         # the chapter-03 figure (no solve)

cd ../04_master
python src/master_figures.py     # the master-upgrade figures: floor ladder + benches (no solve)

cd ../appendix
python src/cancellation_floquet.py # the leak-cancellation Floquet test  (~minutes; qutip)
python src/make_figure.py          # the appendix figure (no solve)
```

Each chapter is `README.md` (the physics) + `src/` (the code) + `images/` (the figures). There is no separate test
runner: each script prints its own self-check, and the headline floor (chapter 01) is a formula you can check by hand.

---

## What this model does (and does not) include

So the numbers above are read with the right scope:

- **Single atom, on axis, axial motion only** — no atom cloud, no radial motion or radial–axial coupling.
  Off-axis the atom samples a weaker 1064 intensity, and since Ω_c² ∝ I but ν_z ∝ √I the EIT condition
  Ω_c²/4Δ = ν_z drifts, walking the bright peak off the sideband — so **every number here is a radially-localized
  best case**; a radially-hot atom cools worse (the radial layer is deliberately out of this repo).
- **The fibre enters as waist and delivery, not as light–matter physics.** The hollow-core fibre sets the mode
  waist (19 µm) and carries the single-EOM delivery chain (comb, retro, tag); the cooling physics computed here is
  that of the 1064 nm lattice and is otherwise identical in free space. Fibre-specific effects — the guided mode's
  longitudinal (E_z) field and its vector shift, surface / Casimir–Polder forces, and propagation loss / beam
  depletion — are out of scope (beam depletion is chapter 08).
- **Anti-trapped excited state — estimated negligible (sudden approximation).** The 5P₃/₂ is anti-trapped (+38),
  but its 26 ns lifetime is ≪ the 2.3 µs trap period, so the atom is *frozen* during the excited excursion; a
  sudden-approximation estimate (setting the excited curvature to ±ν_z) moves the floor by < 2 %. The
  shared-potential approximation the solver uses is therefore safe. (Chapter 05 is slated to compute this heating
  explicitly rather than estimate it.)
- **Linear polarization assumed** — the vector (circular) light shift is dropped. Along the quantization axis it
  acts like a fictitious B-field, which the g_F·m_F-matched pair cancels just like a real one (chapter 01); the
  residual is the transverse/spatially-varying part — again a radial effect.
- **Off-resonant tones treated incoherently** (chapter 02) — they are in fact phase-locked to the Λ; the incoherent
  rate is valid because they sit 100s of MHz off (interference suppressed), while the near-resonant master repumper
  of chapter 03 is treated coherently.
- **Recoil is on every channel, but only its axial projection.** In *both* solvers the two counter-propagating
  Λ legs (absorption) and the spontaneous emission carry the exact recoil (via the displacement operator, Lamb–Dicke
  η = 0.094); the multilevel adds it to the repumper scattering too. All of it acts on the single **axial**
  oscillator. The atom's **radial** motion is not modelled, so the *radial* component of the emission recoil —
  the 3D heating of the transverse motion — is out of scope here; the parent study treats that radial walk
  separately, and in this repo it is the planned chapter 06 (the off-axis cloud). The on-axis axial floor quoted
  here is therefore a radially-cold best case (see the first bullet).
- **Perfect two-photon servo** — the servo is taken to hold δ₂ exactly on the dark resonance (δ₂ = 0 in the
  3-level, ≈ −0.15 in the multilevel; chapter 02) with zero two-photon linewidth. Because the floor is steep in δ₂
  (chapter 02), this is the load-bearing idealisation: in practice the Raman-coherence linewidth — laser phase
  noise, servo jitter — sets the delivered floor.
- **No technical noise** — no laser intensity/phase noise, no magnetic-field noise, no trap-frequency jitter.
- **Repumper model** (chapter 02) is a low-saturation *incoherent* rate — trustworthy only near natural power (see
  its scope note).
- **Detuning reference** — the per-(F′,m′) 1064 tensor Stark on F′1/F′3 *is* included (via `stark.py`); sub-MHz
  residuals in the bare-F′ hyperfine reference are not.
- **The quoted digits are computed values** at the `config.py` parameters, not a claim of physical precision to
  that many figures — read 0.0020 as ≈ 2×10⁻³.

**Sensitivity.** The floor scales as (Γ/4Δ)², so a ±10 % drift in Δ shifts it by ≈ ∓20 %. The cooling itself
relies on the AC-Stark condition Ω_c²/4Δ = ν_z holding, so the bright peak stays on the cooling sideband; a
few-% drift in Ω_c or Δ is tolerable, and a full Δ/Ω_c sensitivity sweep is a natural next check (not yet done).

Chapter 01 is the supervisable heart: the EIT mechanism and the (Γ/4Δ)² floor. Chapters 02–03 are the honest price
of the real delivery and the F′1 leak; the chapter map above (04–08) is what is *not* modelled yet. If any line of
physics here doesn't follow, that's a writing bug — flag it.
