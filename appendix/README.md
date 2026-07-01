# Appendix — the leak frontier

This appendix sits **outside the numbered 01–08 chapters**. Those build the model up in realism; this is a set of
focused deep-dives on the question chapter 03 raised — the **F′1 leak** that caps the D2 master floor at ≈ 0.06 —
and on what, if anything, beats it.

| note | question | verdict |
|---|---|---|
| [`cancellation.md`](cancellation.md) | Can a co-propagating tone cancel the leak? | **No** (Floquet test). The D2 scheme bottoms at **n̄_z ≈ 0.055** (Δ ≈ 25). |
| [`d1_comparison.md`](d1_comparison.md) | Is there a transition with a smaller leak? | **D1**: leak ~41× weaker — and it is the **detuning** (817 vs 157 MHz), not branching. Floor = the deferred D1 chapter. |
| [`magnetic_noise.md`](magnetic_noise.md) | Is the field-insensitive clock pair worth its leak? | A field-noise trade vs a leak-free but field-sensitive scheme; crossover ~40 mG (model-dependent). |
| [`options.md`](options.md) | What else was considered? | W1–W4 + a dRSC pivot, each with a verdict; the clean exits are D1-EIT and dRSC. |

The atomic facts these rest on (R_c/R_p on both lines, the splittings, the no-Λ and GEM selection rules, the
field slopes) are recomputed from first principles in [`verify_atomic_claims.py`](verify_atomic_claims.py) —
`python verify_atomic_claims.py` prints them against the quoted values.

**Bottom line.** On D2 the leak is irreducible: the clock pair delivers ≈ 0.055 (field-insensitive) and the only
colder D2 option is the colleague's field-sensitive stretched scheme (≈ 0.041 below ~40 mG noise). The clean way
past the leak is a different transition — **D1**, where the spoiler is 5× further off resonance — or **dRSC**,
which sidesteps the dark-resonance class altogether. The floor numbers for those are future work (the D1 chapter
is deferred until its from-scratch solver is verified).
