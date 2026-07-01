"""
config.py -- every physical number for the multilevel solve, in one place.
Read by stark.py (Stark shifts), cooling_multilevel.py (the floor), level_scheme.py (figure),
explore_configs.py (the configuration sweep).

Frequencies are ANGULAR, in 2*pi*MHz (a "6.07" means 2*pi*6.07 MHz).
Polarizabilities are in atomic units. The physics is in the top-level README.
"""
import numpy as np

# --- the 1064 nm lattice: two counter-propagating beams, 1 W each ---
P_beam = 1.0          # power per beam (W)
w0     = 19e-6        # 1/e^2 intensity waist (m) -- the kagome HCPCF mode at 1064

# --- dynamic polarizabilities at 1064 nm (a.u.), Chen-Raithel PRA 92, 060501(R) (2015) ---
alpha0_5S   =  687.3  # ground 5S_1/2 scalar   (>0 -> trapped)
alpha0_5P32 = -1149.0 # excited 5P_3/2 scalar  (<0 -> anti-trapped)
alpha2_5P32 =  563.0  # excited 5P_3/2 tensor

# --- the atom / cooling transition (87Rb D2; clock-EIT Lambda |1,-1>,|2,+1> -> |F'2,0>) ---
Gamma = 6.07          # 5P_3/2 natural linewidth (2pi*MHz)
nu_z  = 0.430         # axial trap frequency (2pi*MHz) = 2pi*430 kHz  (measured)
eta   = 0.094         # axial Lamb-Dicke parameter

# --- operating point ---
Delta = 45.0          # single-photon detuning, blue of |F'2,0>  (2pi*MHz)
OmR   = 0.12          # probe/control Rabi ratio

# --- the FULL multilevel solve (cooling_multilevel.py): real 87Rb D2 manifold + repumpers ---
#   8 ground sublevels (F=1: m=-1..1; F=2: m=-2..2) + 5P3/2 F'=1,2(,3); the same Lambda,
#   PLUS the two sigma repumpers that clear the dark ground states the Lambda alone leaks into.
B_field      = 1.0    # magnetic field (G); cooling field, axial. Floor is B-insensitive (clock pair).
theta_trap   = 90.0   # trap linear-pol angle to the axial B (deg). 90 = transverse lattice (real);
                      #   0 = pol||B (reproduces the polarizability-authority F' Stark shifts).
# Single-EOM tagged-retro delivery: ONE seed -> EOM (f_mod = A_HFS + 2f_A = 6.83 + 0.40 = 7.23 GHz)
# -> ... -> tag AOM (2f_A, DOWN-shifts the retro) -> retro. The repumpers are NOT separate lasers --
# they are the leftover tones of the SAME comb, deliberately left OFF-RESONANCE:
#   repump1 = forward +1 EOM sideband (sigma-, F=1) at probe + 2f_A
#             -> drives F=1->F'2, 445 MHz off. (F'3 sits 178 MHz away, but F=1->F'3 is dF=2 FORBIDDEN.)
#   repump2 = retro carrier           (sigma+, F=2) at control - 2f_A
#             -> drives F=2->F'1, 198 MHz off. (F'0 sits 126 MHz away, but F=2->F'0 is dF=2 FORBIDDEN.)
twofA        = 400.0  # double-passed tag-AOM total shift 2*f_A (2pi MHz) = 200 MHz AOM, x2
eta_dp       = 0.30   # retro double-pass power efficiency (sets the repumper field amplitudes)
rep_scale    = 1.0    # multiply the chain-natural repumper Rabis -- raise it to repump FASTER
                      #   (off-resonant repumping is slow; more power/depth speeds it up)
Nf_multi     = 5      # Fock cutoff for the multilevel solve (kept small; the cold floor is low-n). The cooled
                      #   floor is Nf-converged (it changes by <~4% up to Nf=13); explore_configs / chapter 04 use Nf=6.

# The master laser (chapter 03) is a DESIGN, not a number computed here: it sits on the F'1 resonance, where the
# off-resonant (incoherent) repumper model below breaks down. Pinning its floor needs the F'1 vertex treated
# coherently -- that is chapter 04. The master scheme figure is drawn by level_scheme.py; see ../03_master/.
