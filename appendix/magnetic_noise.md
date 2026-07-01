# Magnetic noise: the clock-vs-stretched crossover

The F′1 leak is a property of the **clock pair**: |1,−1⟩ + |2,+1⟩ → |F′2,0⟩ lives at m′=0, where the spoiler
|F′1,0⟩ sits 157 MHz away. One way to dodge it is geometric — cool on a *different* Λ that has no second dark
vertex — but the clean candidate is **field-sensitive**, so the choice becomes a magnetic-noise trade-off.

## 1. The two D2 schemes

- **Clock pair** (this repo): |1,−1⟩(σ⁻) + |2,+1⟩(σ⁻… σ⁺) → |F′2,0⟩. Both legs have g_F·m_F = +½, so the
  two-photon resonance is **first-order field-insensitive** (slope 0, exact). It pays the F′1 leak — floor
  ≈ 0.055–0.086 (chapters 04–05 / [cancellation.md](cancellation.md)).
- **Stretched** ([colleague's scheme](https://github.com/valepioli/EIT_cooling_simulation), build
  `EIT_cooling_Rb_F1_axial`): |1,−1⟩(σ⁻) + |2,−2⟩(π) → |F′2,−2⟩, the dark state at the edge m′=−2. F′=1 has no
  m′=−2 sublevel (max |m′|=1), so neither leg can reach it — **no second dark vertex, no leak**. But the legs
  have g_F·m_F = +½ and −1, a differential of **+3/2 → 2.10 MHz/G**: strongly **field-sensitive** (run at B = 4 G;
  the static offset is compensated, field *noise* jitters the dark resonance).

| | clock pair | stretched |
|---|---|---|
| F′1 leak | yes (157 MHz) → ≈ 0.055–0.086 | none (edge m′=−2) → ≈ 0.041 |
| field sensitivity | 0 (immune) | 2.10 MHz/G |
| limited by | the leak | magnetic noise |

(The 0 and 2.10 MHz/G slopes are verified in [`verify_atomic_claims.py`](verify_atomic_claims.py); the ≈ 0.086 /
≈ 0.041 anchors are the chapter-02 leak-on / leak-off floors.)

## 2. The crossover: σ_B ≈ 40 mG  *(model-dependent)*

The clock floor is flat in field noise (field-insensitive); the stretched floor is the leak-free EIT floor with
its two-photon detuning jittered by σ_δ = 2.10 MHz/G · σ_B, averaged over a quasi-static Gaussian. The two cross
where the stretched scheme's noise penalty overtakes the clock pair's leak penalty:

```
  clock pair (flat)                    n_z ≈ 0.086
  stretched, on resonance              n_z ≈ 0.041
  stretched, sigma_B = 40 mG           n_z ≈ 0.090   <- crosses the clock pair
```

**σ_B\* ≈ 39 mG (1σ)** ≈ 0.2 ν_z of two-photon jitter. Below it the stretched scheme is colder; above it the
field-insensitive clock pair wins. In a typical shielded lab (σ_B ~ few–tens of mG), that puts the *cooling*
advantage with the stretched scheme. **This crossover is a noise-model result, not a selection rule** — it
depends on the quasi-static-Gaussian averaging and on the stretched scheme realizing its full leak-free floor —
so read 39 mG as an **upper bound**; the honest range is ~25–40 mG.

## 3. Caveats that decide the call

- **Cooling-floor comparison only.** The clock pair's field-insensitivity is for the *interrogation*; if the
  science needs field-insensitive states, that is a separate and possibly decisive reason for it, independent of
  this crossover (you can also cool stretched and transfer to clock states for readout).
- **The stretched floor is a best case.** It is the leak-free EIT floor used as a proxy. Interior / swapped dark
  legs can suffer repump-clearability failures (unclearable |2,+2⟩-type leaks → heating) that a clean Λ avoids;
  if the stretched scheme inherits any of that, its floor rises and the crossover moves **lower** (clock pair
  favoured at less noise). A faithful solve of the actual stretched repump topology would pin it.
- **Quasi-static noise assumed** (shot-to-shot). Genuinely fast noise broadens rather than jitters the dark
  resonance — same order, different model.

## 4. The D1 field-bonus is robust *(pending the D1 floor chapter)*

The D1 line ([d1_comparison.md](d1_comparison.md)) has essentially no leak, and a separate observation is that its
accidental-CPT field-bonus (a small extra channel near B ≈ 0.10 G) is **not** razor-thin: averaging over field
noise, n̄ stays below ~0.005 for σ_B ≲ 30 mG and floors out at the robust no-bonus value — never worse. (The main
D1 cooling is field-insensitive at any B; only the bonus channel is B-sensitive.) The specific n̄ numbers here
need the from-scratch D1 solve and are part of the **deferred D1 floor chapter**, not computed in this repo.

## 5. What it says to do

- **If 795 nm is available, D1 removes the dilemma** — field-insensitive *and* leak-free, so neither the
  stretched scheme's field sensitivity nor the clock pair's leak is needed.
- **If stuck on D2,** the stretched scheme is the colder *cooler* at realistic noise (σ_B ≲ 40 mG), if its repump
  topology is clean. The clock pair's field-insensitivity wins on the cooling floor only above ~40 mG — so its D2
  case rests on the *interrogation* needing field-insensitive states, not on the cooling number.
