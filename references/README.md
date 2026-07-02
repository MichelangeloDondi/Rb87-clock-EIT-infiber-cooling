# References

The sources cited across this repository, grouped by what they support — the works the built
chapters (01–04) and the appendix use. The machine-readable form is [`references.bib`](references.bib); each
chapter's `README.md` ends with a short **References** list of the works it cites.

Two kinds of citation appear in the repo:

- **Formulas and mechanism** are cited in the chapter prose at the equation.
- **Borrowed numeric constants** (the ⁸⁷Rb atomic data, the polarizabilities) are cited both in
  the prose and in the `config.py` / `stark.py` comment next to the value, so a line-by-line
  reader sees the provenance without leaving the code.

Design choices that are *this study's own* (the operating point Δ = 45, the Lamb–Dicke η = 0.094,
the delivery detunings) carry **no** citation — they are ours, not borrowed.

---

## Method — EIT / sideband cooling theory

- **Morigi, Eschner & Keitel (2000)** — *Ground-state laser cooling using electromagnetically
  induced transparency*, Phys. Rev. Lett. **85**, 4458. [doi](https://doi.org/10.1103/PhysRevLett.85.4458)
  — the EIT-cooling mechanism (dark-state/Fano interference) and the `(Γ/4Δ)²` floor. **The bedrock cite** (chapters 01–02).
- **Eschner, Morigi, Schmidt-Kaler & Blatt (2003)** — *Laser cooling of trapped ions*, J. Opt. Soc.
  Am. B **20**, 1003. [doi](https://doi.org/10.1364/JOSAB.20.001003) — places EIT cooling among the methods (chapter 01/02).
- **Wineland & Itano (1979)** — *Laser cooling of atoms*, Phys. Rev. A **20**, 1521.
  [doi](https://doi.org/10.1103/PhysRevA.20.1521) — the recoil/sideband framework and the recoil limit (chapters 01–02).

## Atomic data & angular-momentum algebra — the borrowed constants

- **Steck** — *Rubidium 87 D Line Data*, [steck.us/alkalidata](http://steck.us/alkalidata) (rev. 2.3.3).
  — the source of `A_HFS = 6.835 GHz`, `Γ = 2π·6.07 MHz`, the 5P₃/₂ hyperfine centroids and the 26 ns lifetime (config.py, cooling_multilevel.py).
- **Chen & Raithel (2015)** — Phys. Rev. A **92**, 060501(R). [doi](https://doi.org/10.1103/PhysRevA.92.060501)
  — the 1064 nm polarizabilities: 5S₁/₂ `α₀ = +687`, 5P₃/₂ `α₀ = −1149` (the +38 anti-trap), tensor `α₂ = +563` a.u. (stark.py). *⚠ confirm full author list/title.*
- **Grimm, Weidemüller & Ovchinnikov (2000)** — *Optical dipole traps for neutral atoms*, Adv. At.
  Mol. Opt. Phys. **42**, 95. [doi](https://doi.org/10.1016/S1049-250X(08)60186-X) — the dipole/AC-Stark trap-depth relation `U = −αI/(2ε₀c)` (chapter 01, stark.py).
- **Edmonds (1957)** — *Angular Momentum in Quantum Mechanics*, Princeton Univ. Press. — the Wigner-6j
  / Clebsch–Gordan conventions behind the dipole line strengths and the tensor-Stark 6j-null (clebsch_gordan_checks.py, stark_validate.py).

## Prior art — EIT / dark-state cooling of neutral ⁸⁷Rb

- **Huang, Chai & Lan (2021)** — Phys. Rev. A **103**, 013305. [doi](https://doi.org/10.1103/PhysRevA.103.013305) ·
  [arXiv](https://arxiv.org/abs/2101.04309) — dark-state (EIT) sideband cooling of an ⁸⁵Rb ensemble to sub-recoil. The closest prior art for the cooling claim.
- **Xin *et al.* (2025)** — *Fast quantum gas formation via EIT cooling*, Nat. Phys. **21**, 63.
  [doi](https://doi.org/10.1038/s41567-024-02677-9) — lattice EIT cooling of neutral ⁸⁷Rb to degeneracy.
- **Chow *et al.* (2023)** — *Fano resonance in excitation spectroscopy and cooling of an optically
  trapped single atom*, [arXiv:2312.06438](https://arxiv.org/abs/2312.06438) — single-neutral-⁸⁷Rb EIT cooling; identifies F′-manifold off-resonant scattering as a heating floor (chapter 02, appendix).

## Platform & delivery — HCPCF and the single-EOM chain

- **Bajcsy *et al.* (2011)** — Phys. Rev. A **83**, 063830. [doi](https://doi.org/10.1103/PhysRevA.83.063830) ·
  [arXiv](https://arxiv.org/abs/1104.5220) — the seminal laser-cooled atoms inside a single-mode hollow-core PCF. Platform origin.
- **Xin, Leong, Chen, Wang & Lan (2018)** — *An atom interferometer inside a hollow-core PCF*, Sci.
  Adv. **4**, eaat9989. [doi](https://doi.org/10.1126/sciadv.1701723) · [arXiv](https://arxiv.org/abs/1705.08062)
  — the fiber-EOM counter-propagating orthogonal-polarization delivery this chain inherits.
- **Marchesini, Dondi *et al.* (2024)** — *All-fiber, near-infrared laser system at 780 nm for atom cooling*, Opt.
  Continuum **3**, 1868. — the all-fiber 780 nm system for the **MOT + molasses** atom-preparation stage. (The EOM-based two-photon delivery of ch02 is a later addition, not in this paper.)
- **Agnew, Lowit & Arnold (2024)** — *Simple tunable phase-locked lasers for quantum technologies*,
  [arXiv:2404.16806](https://arxiv.org/abs/2404.16806) — single fiber-EOM common-mode generation, demonstrated at 6.83 GHz (= the ⁸⁷Rb clock splitting), no OPLL (chapter 02 delivery).
- **Robust single-sideband Raman generation by FBG filtration (2022)** — Opt. Express,
  [arXiv:2205.04875](https://arxiv.org/abs/2205.04875) — the 1560→780 doubled common-mode two-photon source (chapter 02 delivery).
- **Lin *et al.* (2013)** — *Sympathetic EIT laser cooling of motional modes in an ion chain*, Phys.
  Rev. Lett. **110**, 153002. [doi](https://doi.org/10.1103/PhysRevLett.110.153002) — origin of sympathetic-EIT; the spectator-mode mechanism.

## In-fibre cooling & the RSC comparison

- **Leong *et al.* (2020)** — *Large array of Schrödinger-cat states facilitated by an optical
  waveguide*, Nat. Commun. **11**, 5295 — in-fibre resolved-sideband cooling of ⁸⁵Rb to the axial ground state (the RSC contrast to this EIT scheme).
- **Wang *et al.* (2022)** — Phys. Rev. Research **4**, L022058. [doi](https://doi.org/10.1103/PhysRevResearch.4.L022058) ·
  [arXiv](https://arxiv.org/abs/2112.10088) — in-fibre Λ-gray-molasses / sub-Doppler cooling.
- **Kaufman, Lester & Regal (2012)** — *Cooling a single atom in an optical tweezer to its quantum
  ground state*, Phys. Rev. X **2**, 041014. [doi](https://doi.org/10.1103/PhysRevX.2.041014) · [arXiv](https://arxiv.org/abs/1209.3028) — the canonical single-atom RSC baseline.

## Clock pair, DLS, and deferred-chapter context

- **Chicireanu *et al.* (2011)** — *Differential light-shift cancellation in a magnetic-field-
  insensitive transition of ⁸⁷Rb*, Phys. Rev. Lett. **106**, 063002. [doi](https://doi.org/10.1103/PhysRevLett.106.063002) ·
  [arXiv](https://arxiv.org/abs/1010.1520) — the closest paper to this exact system: an ⁸⁷Rb field-insensitive clock pair with DLS cancellation (chapter 01).
- **Scharnhorst *et al.* (2018)** — Phys. Rev. A **98**, 023424. [doi](https://doi.org/10.1103/PhysRevA.98.023424) ·
  [arXiv](https://arxiv.org/abs/1711.00732) — a hot spectator (radial) mode degrades the cooled (axial) mode's n̄ (the radial-spoils-axial premise of the deferred chapter 06).
- **Roßnagel *et al.* (2015)** — *Fast thermometry for trapped ions using dark resonances*, New J.
  Phys. **17**, 045004. [doi](https://doi.org/10.1088/1367-2630/17/4/045004) · [arXiv](https://arxiv.org/abs/1412.5014) — dark-resonance thermometry (the deferred chapter 07).
