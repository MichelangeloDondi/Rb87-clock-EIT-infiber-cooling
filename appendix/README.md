# Appendix — the leak frontier

This appendix sits **outside the numbered chapters** (01–04 are built; 05–08 are deferred). Those build the model up in realism; this is a set of
focused deep-dives on the question chapter 03 raised — the **F′1 leak** that caps the D2 master floor at ≈ 0.055 —
and on what, if anything, beats it.

| note | question | verdict |
|---|---|---|
| [`cancellation.md`](cancellation.md) | Can a co-propagating tone cancel the leak? | **No** (Floquet test). The D2 scheme bottoms at **n̄_z ≈ 0.055** (Δ ≈ 25). |
| [`d1_comparison.md`](d1_comparison.md) | Is there a transition with a smaller leak? | **D1**: leak ~41× weaker — and it is the **detuning** (817 vs 157 MHz), not branching. Floor = the deferred D1 chapter. |
| [`magnetic_noise.md`](magnetic_noise.md) | Is the field-insensitive clock pair worth its leak? | A field-noise trade vs a leak-free but field-sensitive scheme; crossover ~40 mG (model-dependent). |
| [`options.md`](options.md) | What else was considered? | W1–W4 + a dRSC pivot, each with a verdict; the clean exits are D1-EIT and dRSC. |

The atomic facts these rest on (R_c/R_p on both lines, the splittings, the no-Λ and GEM selection rules, the
field slopes) are recomputed from first principles in [`verify_atomic_claims.py`](src/verify_atomic_claims.py) —
`python src/verify_atomic_claims.py` prints them against the quoted values.

**Bottom line.** On D2 the leak is irreducible: the clock pair delivers ≈ 0.055 (field-insensitive) and the only
colder D2 option is a field-sensitive stretched scheme (a different Lambda) (≈ 0.041 below ~40 mG noise). The clean way
past the leak is a different transition — **D1**, where the spoiler is 5× further off resonance — or **dRSC**,
which sidesteps the dark-resonance class altogether. The floor numbers for those are future work (the D1 chapter
is deferred until its from-scratch solver is verified).

## References

Full entries in [`../references/`](../references/README.md).

- **Chow *et al.* (2023)**, [arXiv:2312.06438](https://arxiv.org/abs/2312.06438) — F′-manifold off-resonant scattering as a heating floor (the D2 leak; motivates the D1 comparison here).
- **Steck**, *Rubidium 87 D Line Data* — the D1/D2 line data behind the "spoiler 5× further off resonance" argument.
- **Chicireanu *et al.* (2011)**, PRL **106**, 063002 — the ⁸⁷Rb field-insensitive clock pair, the alternative to the field-sensitive stretched scheme weighed here.
