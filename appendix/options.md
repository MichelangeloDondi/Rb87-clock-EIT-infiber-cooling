# Options beyond the F′1 leak — the full menu, with verdicts

Cancellation fails ([cancellation.md](cancellation.md)) and D1 is the clean exit
([d1_comparison.md](d1_comparison.md)). For completeness, here are the other routes that were considered on D2,
each with a verdict. The selection-rule claims are reproduced in
[`verify_atomic_claims.py`](verify_atomic_claims.py); floor numbers that need the full solver are flagged.

## W1 — a pure m=0 (π) pair: no Λ exists

Could the leak be sidestepped by cooling on the m=0 pair under π light? **No — there is no Λ.** Under π (q=0),
|1,0⟩ reaches only F′={0,2} and |2,0⟩ only F′={1,3} (the 0→0, ΔF=0 rule makes CG(|1,0⟩→F′1) and CG(|2,0⟩→F′2)
*exactly* zero). The two legs share **no common excited state**, so no dark superposition can form; rotating B
does not help (it is an atomic-frame selection rule). This is precisely *why* the scheme uses m=±1.

The one useful by-product — call it the **GEM**: |2,0⟩→|F′1,0⟩ (π) is allowed (CG = −0.632) while |2,0⟩→|F′2,0⟩
(π) is forbidden (0). So a π beam on |2,0⟩ is a **pure F′1 knob** that does not perturb the F′2 cooling vertex —
the seed of W2.

## W2 — tripod / auxiliary-ground hybrid: promising but unproven

A *field-insensitive* tripod is impossible: only |1,−1⟩ and |2,+1⟩ have g_F·m_F = +½, so there is no third
matched leg. The interesting case is the **hybrid** with W1: keep the clock pair and add the |2,0⟩ π knob to
supply, from a *separate ground state*, the static |F′1,0⟩ coupling that no co-propagating tone can
(cancellation.md §3). |2,0⟩ is itself field-insensitive (g_F·m_F = 0); its cancellation coherence drifts only at
0.70 MHz/G (a third of the stretched scheme's), and crucially the **cooling stays field-insensitive**, so it would
degrade gracefully to the 0.055 clock floor rather than diverge. In principle this evades the cancellation no-go
(it is a new ground state with a forbidden F′2 coupling — a clean knob the §3 tones could not be).

**But a quick solver bolt-on did not cancel:** adding the |2,0⟩ π beam (Rabi ±3..10, detuning across/off the F′1
resonance) moved the floor < 1 % (≈ 0.047 either way), most likely because the |2,0⟩-dark coherence is not
established by the dissipative dynamics (the comb's retro carrier also drives |2,0⟩ and competes). **Verdict: the
most promising of these four and the only one worth a dedicated chapter, but UNPROVEN** — it needs the repump
rebuilt so |2,0⟩ is reserved for the knob, the coherence prepared, and detuning/Rabi/phase scanned; it may yet hit
a cancellation.md-style obstruction. *(The < 1 % figure is a solver result, not reproduced here.)*

## W3 — Raman pulses to clear F′=1: does not help

Coherently transferring the F=1 population back does **not** lower the floor. The leak's dominant cost is recoil
**heating** — the dark state *scattering* on |F′1,0⟩ — which a coherent transfer cannot undo; the population part
is already handled (better, with entropy removal) by the incoherent repumpers.

## W4 — high field (~200 G): wrong direction

The m′=0 leak is itself **field-insensitive** (m′=0 does not Zeeman-shift), so it does not move at low field;
reaching it would need Paschen-Back fields (hundreds of G), by which point the clock pair's first-order
insensitivity is swamped by the quadratic Zeeman — noise-sensitive. (F′=0 schemes, the other high-field idea, can
only be built from F=1 grounds, which are field-sensitive, and carry a *closer* 72 MHz F′1 leak.) Self-defeating.

## The dRSC pivot — sidesteps the whole class

The strongest alternative leaves EIT entirely: degenerate-Raman sideband cooling has **no dark superposition on a
resonance** — a coherent far-detuned Raman |m,n⟩→|m−1,n−1⟩ (no resonant excited coupling) plus a separate
optical-pumping step. Notes:

- It does not dodge anything D1-EIT does not also dodge; D1-EIT vs dRSC is a hardware call, not a physics one.
- The 1064 nm lattice is too far-detuned to drive the Raman — it needs **separate near-detuned Raman beams** plus
  a polarizer injected down the fibre.
- **Geometry is the crux:** B‖axis gives σ± only (no π), so it cannot drive the Δm=1 ladder. Either tilt B
  (standard) or use the Δm=2 σ⁺/σ⁻ ladder (drivable with k‖B, at the cost of even/odd-m sublattice bookkeeping; in
  F=1 the m=0 state must be fed back by the pump).
- It needs a specific *stable* field — ≈ 0.61 G (Δm=1) or ≈ 0.31 G (Δm=2) — with degeneracy held within the Raman
  linewidth if the field is stable to ~60 mG. Not field-insensitive, but the cooling field need not be the
  interrogation field.
- Floor is η²-limited (η² = 0.0085), typically n̄ ~ 0.02–0.05 *(a solver/estimate band, not computed here)* —
  competitive with D1-EIT.

## Decision summary

The F′1 leak is "solved" on D2 only at ≈ 0.055 (the field-insensitive clock pair) or ≈ 0.041 (the colleague's
field-sensitive stretched scheme, below ~40 mG noise — [magnetic_noise.md](magnetic_noise.md)). Two clean exits:

- **(a) D1-EIT** — the leak is ~41× off-resonance, floor expected ~0.02 (robust), reuses the EIT apparatus, needs
  795 nm. The recommended path if the wavelength is available.
- **(b) dRSC** — sidesteps the leak class entirely, n̄ ~ 0.02–0.05, but needs Raman + polarizer optics and the
  tilted-B / Δm=2 fibre geometry.

The **|2,0⟩ hybrid** (W2) is the speculative D2 fix that might preserve field-insensitivity without 795 nm — worth
a chapter, but unproven.
