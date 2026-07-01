"""
config.py -- every physical number for the 3-level core, in one place.
Three scripts read it:  stark.py (trap + Stark shifts), cooling.py (the floor), plots.py (figures).

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
Omega_c = np.sqrt(4.0 * Delta * nu_z)   # EIT-cooling condition (README) -> ~8.8
Omega_p = OmR * Omega_c
# NB the 3-level pins Om_c; the multilevel (02_multilevel) pins the TOTAL drive Om_tot^2 = Om_c^2 + Om_p^2 -- the
# physical convention (the a.c.-Stark shift goes as Om_tot^2). With Om_p = 0.12 Om_c the two differ by 0.7%, negligible for the floor.

N_fock = 12           # motional Fock cutoff (cooling.py)
