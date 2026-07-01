# A different transition: why the D1 line has a far smaller F′1 leak

The F′1 leak ([cancellation.md](cancellation.md)) is set by the atomic ratio R_c − R_p and the F′1–F′2 detuning.
Both differ between transitions, so the obvious lever is the **D1 line** (5P₁/₂). This note computes the leak on
both lines from Clebsch–Gordan × signed hyperfine-6j ([`verify_atomic_claims.py`](verify_atomic_claims.py)
reproduces every number) and pins down *why* D1 wins. The headline correction: it is the **detuning**, not the
branching.

## The leak still exists on D1 — it is not cancelled

On D1 the two dark-state legs are exactly equal and opposite, R_c = −0.775, R_p = +0.775 (= ∓√15/5, from the
J′=1/2 structure), versus D2's lopsided R_c = −0.346, R_p = +1.732. So the leak amplitude |R_c − R_p| drops from
**2.08 (D2) to 1.55 (D1)** — only ~26 % smaller. **It is not a selection-rule null:** both legs remain
dipole-allowed, R_c − R_p is firmly nonzero, and the dark state is still bright on the D1 spectator |F′1,0⟩. Any
account that says D1 "avoids" the leak by a cancellation is wrong.

## The real reason is the detuning (5.2× further off resonance)

The leak is the dark state's *off-resonant* scatter onto the spectator, ~ Γ·(R_c − R_p)²/(Δ + split)². What
changes decisively on D1 is **split**: the F′1–F′2 hyperfine spacing is **816.7 MHz on D1 vs 156.9 MHz on D2**
(5.2×), because 5P₁/₂ is far more spread out (816.7 = 2·A_{P½}) than 5P₃/₂. Squared, that is a ~27× suppression;
the amplitude buys a further 1.8×. Net the D1 leak scatter is **~41× weaker at Δ = 25 (≈51× at Δ = 0)**:

| line | R_c | R_p | \|R_c−R_p\| | F′1–F′2 split | leak scatter (rel.) |
|---|---|---|---|---|---|
| D2 (5P₃/₂) | −0.346 | +1.732 | 2.08 | 156.9 MHz | 1 |
| D1 (5P₁/₂) | −0.775 | +0.775 | 1.55 | 816.7 MHz | ~1/41 (Δ=25) |

with the suppression decomposing as **27× (detuning) × 1.8× (amplitude)**: detuning carries >90 % of it. The
amplitude and the (true) inverted branching are secondary.

## Inverted branching — a secondary bonus, not the mechanism

A leak event only costs you if the F′1 spoiler decays into the *uncooled* hyperfine. On D2 it dumps **5/6 into
F=1** (the floor — needs aggressive repumping); on D1 the J′=1/2 6j inverts this to **1/6 into F=1, 5/6 back to
F=2** (a cooling manifold, lightly repumped). That is a genuine bonus, but it is a *second-order* benefit on top
of the ~41× scatter suppression — not the reason D1 wins. There is also no F′3 on D1 to add a second spoiler, and
J′=1/2 carries no tensor light shift, so the magic |F′2,0⟩ vertex D2 needed becomes automatic.

## What this does *not* yet establish — the floor

The atomic comparison above is exact. The *floor* it implies (the leak suppressed ~41× would drop below the
recoil/repump floor, leaving D1 **repump-limited near n̄_z ~ 0.02**, no longer leak-limited) needs a from-scratch
5P₁/₂ port of the chapter-02 solver — different excited manifold (F′=1,2 only), D1 branching and g-factors, the
817 MHz splitting, no tensor, 795 nm recoil. **That floor is not computed in this repo** (it is the deferred D1
chapter; the result would only *improve* on the D2 floor, since it removes a present-but-weaker leak). The price,
either way, is hardware: D1 is **795 nm**, off the experiment's 780.

**Verify:** `python verify_atomic_claims.py` reproduces R_c/R_p on both lines, the 156.9/816.7 MHz splittings,
|R_c−R_p|, and the 41×/51× leak suppression with its 27×-detuning × 1.8×-amplitude decomposition.
