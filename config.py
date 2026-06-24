"""
config.py -- physical parameters for the single-atom, on-axis clock-EIT cooling floor.

Layer 1: ONE 87Rb atom, fixed at the trap centre (on axis), cooled on the axial mode
by the clean 3-level clock-EIT Lambda. No contaminants, no retro parasitics, no cloud
average -- those are Layer 2 / Layer 3. This file is just the numbers; the physics is
in README.md, and the only script that reads this is cooling.py.

CONVENTION: every frequency is ANGULAR, written as an ordinary frequency in 2*pi*MHz
(a literal 6.07 means 2*pi * 6.07 MHz). Same convention as the main repo.
"""
import numpy as np

# --- the atom: 87Rb D2, clock-EIT Lambda  |1,-1> (probe, s+) , |2,+1> (control, s-) -> |F'2,0> ---
Gamma = 6.07        # 5P_3/2 natural linewidth        (2pi*MHz)
nu_z  = 0.430       # axial trap frequency            (2pi*MHz)  = 2pi * 430 kHz
eta   = 0.094       # axial Lamb-Dicke parameter (one 780 nm photon)

# --- operating point (the finalised single-ended tagged-retro point; see README) ---
Delta = 45.0        # single-photon detuning, blue of |F'2,0>   (2pi*MHz)
OmR   = 0.12        # probe/control Rabi ratio  Omega_p / Omega_c

# EIT-cooling resonance condition (derived in README): the control AC-Stark-shifts the
# bright state by Omega_c^2/(4 Delta); set that equal to nu_z so the bright (absorbing)
# resonance lands on the COOLING sideband.  ->  Omega_c = sqrt(4 Delta nu_z).
Omega_c = np.sqrt(4.0 * Delta * nu_z)   # ~ 8.8  (2pi*MHz)
Omega_p = OmR * Omega_c                  # weak probe

# --- numerics (used ONLY by cooling.py, the last-resort confirmation) ---
N_fock = 12         # motional Fock cutoff (floor ~1e-3, so 12 states is ample)
